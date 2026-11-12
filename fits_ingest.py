#!/usr/bin/env python3
"""Safely import visibility FITS data into a CASA Measurement Set.

This module deliberately performs a small, dependency-free FITS structure
inspection before loading a CASA task.  A normal image FITS file is therefore
never handed speculatively to a visibility importer.  Classic random-groups
UVFITS and FITS-IDI are supported; ambiguous inputs fail closed unless their
format is selected explicitly.

The importer writes to a sibling staging directory, validates the resulting
Measurement Set with ``casatools.table``, atomically publishes it, writes an
atomic provenance report, and creates the completion marker last.

A GMRT LTA file (``.lta``, ``.ltb``, ``.lta_1`` ...) is first converted to
UVFITS by ``lta_convert.py`` with NCRA's listscan and gvfits, in
``<outdir>/ingest/lta/``, and that FITS is then imported with ``importgmrt``
like any uGMRT FITS. The report's ``source_fits`` is then the LTA file
itself, so a resumed run is bound to the archive file, and
``lta_conversion`` records the conversion. The converted FITS is removed once
the Measurement Set is validated, unless it is asked to be kept.
"""

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent

# Run as `python3 fits_ingest.py` (a modular CASA), the pipeline folder would
# come first on sys.path, where its own inspect.py stands in for the standard
# module. The other stage scripts drop it the same way.
if sys.path and sys.path[0] and Path(sys.path[0]).resolve() == SCRIPT_DIR:
    sys.path.pop(0)


def _load_casa_compat():
    """Load the sibling compatibility module without relying on sys.path."""

    module_path = SCRIPT_DIR / "casa_compat.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_fits_ingest_casa_compat",
        str(module_path),
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load CASA compatibility module: {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


import_casa_task = _load_casa_compat().import_casa_task


def _load_lta_convert():
    """Load the LTA converter only when an LTA file is to be imported."""

    module_path = SCRIPT_DIR / "lta_convert.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_fits_ingest_lta_convert",
        str(module_path),
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load the LTA converter: {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# 1: GMRT LTA input, converted with listscan and gvfits before importgmrt.
PIPELINE_INTERFACE_VERSION = 1

FITS_BLOCK_SIZE = 2880
FITS_CARD_SIZE = 80
CANONICAL_FORMATS = {"auto", "uvfits", "fitsidi", "lta"}
TELESCOPE_CHOICES = {"auto", "ugmrt", "meerkat", "lofar"}
REQUIRED_MS_SUBTABLES = (
    "ANTENNA",
    "DATA_DESCRIPTION",
    "FIELD",
    "OBSERVATION",
    "POLARIZATION",
    "SPECTRAL_WINDOW",
)
REQUIRED_MAIN_COLUMNS = {
    "ANTENNA1",
    "ANTENNA2",
    "DATA_DESC_ID",
    "FIELD_ID",
    "TIME",
}
VISIBILITY_COLUMNS = {"DATA", "FLOAT_DATA", "CORRECTED_DATA"}


class FITSDetectionError(ValueError):
    """Raised when an input is not decisively a supported visibility FITS."""


class MeasurementSetValidationError(RuntimeError):
    """Raised when a CASA importer did not materialize a usable MS."""


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def normalize_fits_format(value):
    """Return the canonical CLI spelling for a FITS visibility format."""

    normalized = str(value or "auto").strip().lower().replace("_", "-")
    if normalized == "fits-idi":
        normalized = "fitsidi"

    if normalized not in CANONICAL_FORMATS:
        choices = ", ".join(sorted(CANONICAL_FORMATS))
        raise ValueError(
            f"Unsupported FITS format '{value}'. Choose one of: {choices}."
        )

    return normalized


def is_lta_path(path):
    """True for a GMRT LTA name: .lta, .ltb or a numbered part like .lta_1."""

    return bool(re.search(r"\.lt[ab](_\d+)?$", Path(path).name, re.IGNORECASE))


def normalize_telescope(value):
    normalized = str(value or "auto").strip().lower()
    if normalized == "gmrt":
        normalized = "ugmrt"

    if normalized not in TELESCOPE_CHOICES:
        choices = ", ".join(sorted(TELESCOPE_CHOICES))
        raise ValueError(
            f"Unsupported telescope '{value}'. Choose one of: {choices}."
        )

    return normalized


def _unquoted_value(text):
    """Remove a FITS comment while respecting quoted string values."""

    quoted = False
    index = 0
    while index < len(text):
        character = text[index]
        if character == "'":
            if quoted and index + 1 < len(text) and text[index + 1] == "'":
                index += 2
                continue
            quoted = not quoted
        elif character == "/" and not quoted:
            return text[:index].strip()
        index += 1
    return text.strip()


def _parse_fits_value(card):
    if len(card) < 10 or card[8] != "=":
        return None

    raw = _unquoted_value(card[10:])
    if not raw:
        return None

    if raw.startswith("'"):
        value = []
        index = 1
        while index < len(raw):
            character = raw[index]
            if character == "'":
                if index + 1 < len(raw) and raw[index + 1] == "'":
                    value.append("'")
                    index += 2
                    continue
                break
            value.append(character)
            index += 1
        return "".join(value).rstrip()

    if raw == "T":
        return True
    if raw == "F":
        return False

    numeric = raw.replace("D", "E").replace("d", "e")
    try:
        if any(token in numeric for token in (".", "E", "e")):
            return float(numeric)
        return int(numeric)
    except ValueError:
        return raw.strip()


def _read_header(handle, offset, file_size):
    if offset < 0 or offset + FITS_BLOCK_SIZE > file_size:
        raise FITSDetectionError(
            "The FITS file ends before the next complete header block."
        )

    handle.seek(offset)
    cards = []
    blocks_read = 0
    found_end = False

    # A 2.8 MiB header is already far beyond ordinary visibility FITS usage.
    # The limit also prevents a damaged file from causing an unbounded scan.
    while blocks_read < 1024:
        block = handle.read(FITS_BLOCK_SIZE)
        if len(block) != FITS_BLOCK_SIZE:
            raise FITSDetectionError("The FITS header is truncated.")
        blocks_read += 1

        for start in range(0, FITS_BLOCK_SIZE, FITS_CARD_SIZE):
            raw_card = block[start:start + FITS_CARD_SIZE]
            try:
                card = raw_card.decode("ascii")
            except UnicodeDecodeError as exc:
                raise FITSDetectionError(
                    "The FITS header contains non-ASCII card data."
                ) from exc
            cards.append(card)
            if card[:8].strip() == "END":
                found_end = True
                break

        if found_end:
            break

    if not found_end:
        raise FITSDetectionError(
            "No END card was found within the FITS header safety limit."
        )

    header = {}
    for card in cards:
        keyword = card[:8].strip().upper()
        if not keyword or keyword in {"COMMENT", "HISTORY", "END"}:
            continue
        if card[8:10] != "= ":
            continue
        # Preserve the first occurrence; duplicate structural cards are not
        # allowed to silently change how the file is classified or skipped.
        header.setdefault(keyword, _parse_fits_value(card))

    return header, blocks_read * FITS_BLOCK_SIZE


def _positive_integer(header, name, default=0):
    try:
        value = int(header.get(name, default))
    except (TypeError, ValueError) as exc:
        raise FITSDetectionError(
            f"FITS header keyword {name} is not an integer."
        ) from exc
    if value < 0:
        raise FITSDetectionError(
            f"FITS header keyword {name} cannot be negative."
        )
    return value


def _hdu_data_size(header):
    try:
        bitpix = abs(int(header.get("BITPIX", 0)))
    except (TypeError, ValueError) as exc:
        raise FITSDetectionError("FITS BITPIX is not an integer.") from exc

    if bitpix == 0 or bitpix % 8:
        raise FITSDetectionError("FITS BITPIX must be a non-zero multiple of 8.")

    naxis = _positive_integer(header, "NAXIS")
    axes = [_positive_integer(header, f"NAXIS{axis}") for axis in range(1, naxis + 1)]
    pcount = _positive_integer(header, "PCOUNT")
    gcount = _positive_integer(header, "GCOUNT", 1)
    if gcount == 0:
        return 0

    if bool(header.get("GROUPS")) and axes and axes[0] == 0:
        elements = math.prod(axes[1:]) if len(axes) > 1 else 0
    else:
        elements = math.prod(axes) if axes else 0

    return (bitpix // 8) * gcount * (pcount + elements)


def _padded_size(size):
    if size == 0:
        return 0
    return ((size + FITS_BLOCK_SIZE - 1) // FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE


def classify_telescope_names(values):
    """Map observatory name strings onto supported telescope identities."""

    detected = set()
    for value in values or ():
        upper = str(value).upper()
        if "MEERKAT" in upper:
            detected.add("meerkat")
        elif "LOFAR" in upper:
            detected.add("lofar")
        elif "GMRT" in upper:
            detected.add("ugmrt")

    return sorted(detected)


def _header_telescope(headers):
    values = []
    for header in headers:
        for keyword in ("TELESCOP", "ARRNAM"):
            value = str(header.get(keyword, "")).strip()
            if value:
                values.append(value)

    return values, classify_telescope_names(values)


def inspect_fits_structure(path, max_hdus=256):
    """Read enough FITS metadata to classify visibility containers.

    No array or visibility payload is read.  The returned dictionary is JSON
    serializable and suitable for inclusion in the ingest provenance report.
    """

    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Input FITS file not found: {path}")

    file_size = path.stat().st_size
    if file_size < FITS_BLOCK_SIZE:
        raise FITSDetectionError("The input is too small to contain a FITS HDU.")

    headers = []
    offset = 0
    with open(path, "rb") as handle:
        while offset < file_size and len(headers) < max_hdus:
            header, header_size = _read_header(handle, offset, file_size)
            if not headers and header.get("SIMPLE") is not True:
                raise FITSDetectionError(
                    "The input does not begin with a valid SIMPLE FITS header."
                )
            if headers and not str(header.get("XTENSION", "")).strip():
                raise FITSDetectionError(
                    "A non-primary FITS HDU is missing XTENSION."
                )

            data_size = _hdu_data_size(header)
            next_offset = offset + header_size + _padded_size(data_size)
            if next_offset > file_size:
                raise FITSDetectionError(
                    "A FITS HDU declares data beyond the end of the file."
                )

            headers.append(header)
            if next_offset == offset:
                raise FITSDetectionError("FITS HDU scanning made no progress.")
            offset = next_offset

        if offset < file_size and len(headers) >= max_hdus:
            raise FITSDetectionError(
                f"The FITS file exceeds the safety limit of {max_hdus} HDUs."
            )
        if offset != file_size:
            raise FITSDetectionError(
                "Bytes remain after the final complete FITS HDU; the file may "
                "be truncated or malformed."
            )

    primary = headers[0]
    naxis = _positive_integer(primary, "NAXIS")
    primary_axes = [
        _positive_integer(primary, f"NAXIS{axis}")
        for axis in range(1, naxis + 1)
    ]
    random_groups_uvfits = bool(
        primary.get("GROUPS") is True
        and naxis >= 2
        and primary_axes
        and primary_axes[0] == 0
        and _positive_integer(primary, "PCOUNT") > 0
        and _positive_integer(primary, "GCOUNT", 1) > 0
        and math.prod(primary_axes[1:]) > 0
    )

    extension_names = []
    fitsidi_uv_data = False
    for header in headers[1:]:
        extension_name = str(header.get("EXTNAME", "")).strip().upper()
        if extension_name:
            extension_names.append(extension_name)
        if (
            str(header.get("XTENSION", "")).strip().upper() == "BINTABLE"
            and extension_name == "UV_DATA"
            and _positive_integer(header, "NAXIS1") > 0
            and _positive_integer(header, "NAXIS2") > 0
        ):
            fitsidi_uv_data = True

    primary_axes_meta = []
    for axis in range(1, naxis + 1):
        primary_axes_meta.append({
            "axis": axis,
            "ctype": str(primary.get(f"CTYPE{axis}", "")).strip(),
            "naxis": primary_axes[axis - 1],
            "crval": primary.get(f"CRVAL{axis}"),
            "cdelt": primary.get(f"CDELT{axis}"),
            "crpix": primary.get(f"CRPIX{axis}"),
        })

    telescope_values, detected_telescopes = _header_telescope(headers)
    primary_image = bool(
        not random_groups_uvfits
        and naxis > 0
        and primary_axes
        and math.prod(primary_axes) > 0
    )

    return {
        "hdu_count": len(headers),
        "extension_names": extension_names,
        "random_groups_uvfits": random_groups_uvfits,
        "fitsidi_uv_data": fitsidi_uv_data,
        "primary_image": primary_image,
        "telescope_header_values": telescope_values,
        "detected_telescopes": detected_telescopes,
        "primary_gcount": _positive_integer(primary, "GCOUNT", 1),
        "primary_pcount": _positive_integer(primary, "PCOUNT"),
        "primary_axes": primary_axes_meta,
    }


def _axis_by_type(summary, token):
    """Return the primary-HDU axis whose CTYPE starts with *token*."""

    for axis in summary.get("primary_axes", []) or ():
        if str(axis.get("ctype", "")).upper().startswith(token):
            return axis
    return None


def compare_fits_and_ms(summary, validation):
    """Cross-check the imported Measurement Set against the source FITS.

    The validation contract asks for the imported row count, correlation and
    spectral structure to be compared with the FITS metadata before an ingest
    is accepted.  Nothing performed that comparison, so a partial or
    restructured import looked identical to a complete one.

    Every check reports ``consistent``, ``inconsistent`` or
    ``not_comparable``; the last is used wherever the FITS convention cannot
    determine the answer on its own, rather than guessing.
    """

    checks = {}

    if not summary.get("random_groups_uvfits"):
        return {
            "status": "not_comparable",
            "reason": (
                "structural cross-checks are defined here for random-groups "
                "UVFITS only"
            ),
            "checks": checks,
        }

    groups = int(summary.get("primary_gcount") or 0)
    if_axis = _axis_by_type(summary, "IF")
    stokes_axis = _axis_by_type(summary, "STOKES")
    freq_axis = _axis_by_type(summary, "FREQ")

    spectral = validation.get("spectral", {}) or {}
    polarization = validation.get("polarization", {}) or {}
    main_rows = int(validation.get("main_rows") or 0)
    spw_count = int(spectral.get("spectral_window_count") or 0)

    # One MS row per (random group, IF).
    if groups > 0 and spw_count > 0:
        expected_rows = groups * (
            int(if_axis["naxis"]) if if_axis else spw_count
        )
        checks["main_rows"] = {
            "status": (
                "consistent" if main_rows == expected_rows else "inconsistent"
            ),
            "fits_groups": groups,
            "fits_if_axis": int(if_axis["naxis"]) if if_axis else None,
            "measurement_set_spectral_windows": spw_count,
            "expected_rows": expected_rows,
            "measurement_set_rows": main_rows,
        }
    else:
        checks["main_rows"] = {
            "status": "not_comparable",
            "reason": "no usable GCOUNT or spectral-window count",
        }

    correlations = polarization.get("correlation_count")
    if stokes_axis and correlations:
        checks["correlations"] = {
            "status": (
                "consistent"
                if int(stokes_axis["naxis"]) == int(correlations)
                else "inconsistent"
            ),
            "fits_stokes_axis": int(stokes_axis["naxis"]),
            "measurement_set_correlations": int(correlations),
        }
    else:
        checks["correlations"] = {
            "status": "not_comparable",
            "reason": "no STOKES axis or correlation count",
        }

    # The primary FREQ axis describes a single IF; with several IFs the true
    # coverage lives in the FQ table, so do not pretend to compare it.
    channels = spectral.get("channel_counts") or []
    if (
        freq_axis and
        spw_count == 1 and
        (if_axis is None or int(if_axis["naxis"]) == 1) and
        len(channels) == 1
    ):
        checks["channels"] = {
            "status": (
                "consistent"
                if int(freq_axis["naxis"]) == int(channels[0])
                else "inconsistent"
            ),
            "fits_freq_axis": int(freq_axis["naxis"]),
            "measurement_set_channels": int(channels[0]),
        }
    else:
        checks["channels"] = {
            "status": "not_comparable",
            "reason": (
                "the primary FREQ axis describes one IF; multi-IF coverage "
                "is defined by the FQ table"
            ),
        }

    inconsistent = sorted(
        name
        for name, item in checks.items()
        if item.get("status") == "inconsistent"
    )

    return {
        "status": "inconsistent" if inconsistent else (
            "consistent"
            if any(
                item.get("status") == "consistent"
                for item in checks.values()
            )
            else "not_comparable"
        ),
        "inconsistent_checks": inconsistent,
        "checks": checks,
    }


def resolve_visibility_format(summary, requested="auto"):
    requested = normalize_fits_format(requested)
    signatures = []
    if summary.get("random_groups_uvfits"):
        signatures.append("uvfits")
    if summary.get("fitsidi_uv_data"):
        signatures.append("fitsidi")

    if requested != "auto":
        if requested not in signatures:
            observed = ", ".join(signatures) if signatures else "none"
            raise FITSDetectionError(
                f"--fits-format {requested} does not match the FITS "
                f"visibility structure (detected signatures: {observed})."
            )
        return requested

    if len(signatures) == 1:
        return signatures[0]
    if len(signatures) > 1:
        raise FITSDetectionError(
            "The file contains both random-groups UVFITS and FITS-IDI "
            "UV_DATA signatures. Select --fits-format explicitly after "
            "verifying the producer."
        )
    if summary.get("primary_image"):
        raise FITSDetectionError(
            "The input is an ordinary image FITS, not a visibility FITS. "
            "Provide random-groups UVFITS, FITS-IDI, or a Measurement Set."
        )
    raise FITSDetectionError(
        "The file has no supported visibility signature. Expected a "
        "random-groups UVFITS primary HDU or a FITS-IDI UV_DATA table."
    )


def detect_visibility_fits(path, requested="auto"):
    """Return ``uvfits`` or ``fitsidi`` after strict content detection."""

    return resolve_visibility_format(inspect_fits_structure(path), requested)


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def build_input_identity(path):
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "sha256": sha256_file(path),
    }


def atomic_write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}.{time.time_ns()}")
    try:
        with open(temporary, "w") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_write_checkpoint(path):
    path = Path(path)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}.{time.time_ns()}")
    try:
        with open(temporary, "w") as handle:
            handle.write(utc_now() + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def read_json(path, default=None):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (FileNotFoundError, OSError, ValueError, TypeError):
        return default


def _table_summary(path, table_factory):
    tool = table_factory()
    opened = False
    try:
        tool.open(str(path), nomodify=True)
        opened = True
        return int(tool.nrows()), sorted(str(name) for name in tool.colnames())
    finally:
        if opened:
            try:
                tool.close()
            except Exception:
                pass
        done = getattr(tool, "done", None)
        if callable(done):
            try:
                done()
            except Exception:
                pass


def _iter_scalars(value):
    """Yield scalars from a CASA/NumPy array or nested container."""

    flat = getattr(value, "flat", None)
    if flat is not None:
        try:
            yield from flat
            return
        except TypeError:
            pass

    if isinstance(value, dict):
        for item in value.values():
            yield from _iter_scalars(item)
        return

    if isinstance(value, (list, tuple)):
        for item in value:
            yield from _iter_scalars(item)
        return

    yield value


def _table_column(path, table_factory, column, default=None):
    """Read one column, tolerating variable-shaped CASA columns."""

    tool = table_factory()
    opened = False
    try:
        tool.open(str(path), nomodify=True)
        opened = True
        if column not in {str(name) for name in tool.colnames()}:
            return default
        try:
            return tool.getcol(column)
        except Exception:
            getvarcol = getattr(tool, "getvarcol", None)
            if callable(getvarcol):
                return list(getvarcol(column).values())
            return default
    except Exception:
        return default
    finally:
        if opened:
            try:
                tool.close()
            except Exception:
                pass
        done = getattr(tool, "done", None)
        if callable(done):
            try:
                done()
            except Exception:
                pass


def validate_measurement_set(path, table_factory=None):
    """Validate a materialized MS and return a JSON-safe structural summary."""

    path = Path(path)
    if path.is_symlink():
        raise MeasurementSetValidationError(
            f"Refusing a symlink as the imported Measurement Set: {path}"
        )
    if not path.is_dir() or not (path / "table.dat").is_file():
        raise MeasurementSetValidationError(
            f"Imported output is not a CASA table (table.dat missing): {path}"
        )

    for name in REQUIRED_MS_SUBTABLES:
        subtable = path / name
        if not subtable.is_dir() or not (subtable / "table.dat").is_file():
            raise MeasurementSetValidationError(
                f"Imported Measurement Set is missing required subtable {name}."
            )

    if table_factory is None:
        try:
            from casatools import table as table_factory
        except Exception as exc:
            raise MeasurementSetValidationError(
                "casatools.table is required to validate an imported "
                "Measurement Set. Run FITS ingestion inside CASA 6."
            ) from exc

    try:
        main_rows, main_columns = _table_summary(path, table_factory)
    except Exception as exc:
        raise MeasurementSetValidationError(
            f"CASA could not open the imported Measurement Set MAIN table: {exc}"
        ) from exc

    if main_rows <= 0:
        raise MeasurementSetValidationError(
            "The imported Measurement Set contains no visibility rows."
        )
    missing_columns = sorted(REQUIRED_MAIN_COLUMNS.difference(main_columns))
    if missing_columns:
        raise MeasurementSetValidationError(
            "The imported Measurement Set MAIN table is missing columns: "
            + ", ".join(missing_columns)
        )
    data_columns = sorted(VISIBILITY_COLUMNS.intersection(main_columns))
    if not data_columns:
        raise MeasurementSetValidationError(
            "The imported Measurement Set has no visibility data column."
        )

    subtables = {}
    for name in REQUIRED_MS_SUBTABLES:
        try:
            rows, columns = _table_summary(path / name, table_factory)
        except Exception as exc:
            raise MeasurementSetValidationError(
                f"CASA could not open imported MS subtable {name}: {exc}"
            ) from exc
        if rows <= 0:
            raise MeasurementSetValidationError(
                f"Imported MS subtable {name} contains no rows."
            )
        subtables[name] = {"rows": rows, "columns": columns}

    # A structurally complete MS can still carry an unusable spectral or
    # polarization axis. Catching that here fails closed at ingest with an
    # actionable message, rather than three stages later as an obscure band
    # detection or correlation-selection error.
    channel_counts = []
    for value in _iter_scalars(
        _table_column(path / "SPECTRAL_WINDOW", table_factory, "NUM_CHAN", [])
    ):
        try:
            channel_counts.append(int(value))
        except (TypeError, ValueError):
            continue

    frequencies = []
    for value in _iter_scalars(
        _table_column(path / "SPECTRAL_WINDOW", table_factory, "CHAN_FREQ", [])
    ):
        try:
            frequencies.append(float(value))
        except (TypeError, ValueError):
            continue

    if not channel_counts or any(count < 1 for count in channel_counts):
        raise MeasurementSetValidationError(
            "The imported Measurement Set has a spectral window with no "
            f"channels: NUM_CHAN={channel_counts}"
        )
    if not frequencies:
        raise MeasurementSetValidationError(
            "The imported Measurement Set has no channel frequencies."
        )
    if not all(math.isfinite(value) and value > 0.0 for value in frequencies):
        raise MeasurementSetValidationError(
            "The imported Measurement Set contains non-finite or "
            "non-positive channel frequencies."
        )

    correlation_counts = []
    for value in _iter_scalars(
        _table_column(path / "POLARIZATION", table_factory, "NUM_CORR", [])
    ):
        try:
            correlation_counts.append(int(value))
        except (TypeError, ValueError):
            continue

    correlation_types = []
    for value in _iter_scalars(
        _table_column(path / "POLARIZATION", table_factory, "CORR_TYPE", [])
    ):
        try:
            correlation_types.append(int(value))
        except (TypeError, ValueError):
            continue

    if not correlation_counts or any(
        count < 1 for count in correlation_counts
    ):
        raise MeasurementSetValidationError(
            "The imported Measurement Set has a polarization setup with no "
            f"correlations: NUM_CORR={correlation_counts}"
        )

    telescopes = sorted({
        str(value).strip()
        for value in _iter_scalars(
            _table_column(
                path / "OBSERVATION",
                table_factory,
                "TELESCOPE_NAME",
                [],
            )
        )
        if str(value).strip()
    })

    return {
        "main_rows": main_rows,
        "main_columns": main_columns,
        "data_columns": data_columns,
        "subtables": subtables,
        "spectral": {
            "spectral_window_count": len(channel_counts),
            "channel_counts": channel_counts,
            "frequency_min_hz": min(frequencies),
            "frequency_max_hz": max(frequencies),
        },
        "polarization": {
            "setup_count": len(correlation_counts),
            "correlation_count": max(correlation_counts),
            "correlation_types": sorted(set(correlation_types)),
        },
        "observation": {
            "telescope_names": telescopes,
        },
    }


def _output_path(path):
    raw = Path(path).expanduser()
    if not raw.is_absolute():
        raw = Path.cwd() / raw
    # Resolve only the parent. Resolving an existing final symlink would make
    # a force operation act on its target instead of rejecting the symlink.
    return raw.parent.resolve() / raw.name


def _archive_name(path):
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    candidate = path.with_name(f"{path.name}.replaced.{timestamp}.{os.getpid()}")
    if candidate.exists():
        raise RuntimeError(f"Replacement archive path already exists: {candidate}")
    return candidate


class FITSIngestor:
    """Restart-safe CASA visibility-FITS importer."""

    def __init__(
        self,
        input_fits,
        output_ms,
        outdir="Products",
        fits_format="auto",
        telescope="auto",
        gmrt_flagfiles=None,
        force=False,
        task_loader=None,
        ms_validator=None,
        listscan="listscan",
        gvfits="gvfits",
        keep_converted_fits=False,
        lta_converter=None,
    ):
        self.input_fits = Path(input_fits).expanduser().resolve()
        self.output_ms = _output_path(output_ms)
        self.outdir = _output_path(outdir)
        self.requested_format = normalize_fits_format(fits_format)
        self.requested_telescope = normalize_telescope(telescope)
        self.gmrt_flagfiles = [
            Path(item).expanduser().resolve()
            for item in (gmrt_flagfiles or [])
        ]
        self.force = bool(force)
        self.task_loader = task_loader or import_casa_task
        self.ms_validator = ms_validator or validate_measurement_set
        self.ingest_dir = self.outdir / "ingest"
        self.report_path = self.ingest_dir / "ingest_report.json"
        self.checkpoint_path = self.ingest_dir / "ingest.done"
        # An LTA file is chosen by --fits-format lta or, with auto, by its
        # name: its content is GMRT's own format, which only listscan reads.
        self.input_is_lta = self.requested_format == "lta" or (
            self.requested_format == "auto" and is_lta_path(self.input_fits)
        )
        self.listscan = listscan
        self.gvfits = gvfits
        self.keep_converted_fits = bool(keep_converted_fits)
        self.lta_converter = lta_converter

        try:
            self.ingest_dir.relative_to(self.output_ms)
        except ValueError:
            pass
        else:
            raise ValueError(
                "The output Measurement Set cannot contain the ingest "
                "report/checkpoint directory. Choose a destination below "
                "the ingest directory (the default) or elsewhere."
            )

    def _resolve_telescope(self, summary):
        detected = list(summary.get("detected_telescopes", []))
        if len(detected) > 1:
            raise FITSDetectionError(
                "FITS headers identify multiple supported telescopes: "
                + ", ".join(detected)
            )
        detected_name = detected[0] if detected else None
        if (
            self.requested_telescope != "auto"
            and detected_name
            and detected_name != self.requested_telescope
        ):
            raise FITSDetectionError(
                f"--telescope {self.requested_telescope} conflicts with FITS "
                f"header telescope {detected_name}."
            )
        if self.requested_telescope != "auto":
            return self.requested_telescope
        return detected_name or "auto"

    def _select_task(self, fits_format, telescope):
        if fits_format == "fitsidi":
            return "importfitsidi"
        if fits_format == "uvfits" and telescope == "ugmrt":
            # CASA provides this task specifically for GMRT/uGMRT UVFITS and
            # optionally applies observatory flag files. Do not silently fall
            # back to the generic importer for explicitly identified GMRT data.
            return "importgmrt"
        return "importuvfits"

    def _task_kwargs(self, task_name, staging, source=None):
        source = self.input_fits if source is None else source
        if task_name == "importfitsidi":
            return {"fitsidifile": str(source), "vis": str(staging)}
        values = {"fitsfile": str(source), "vis": str(staging)}
        if task_name == "importgmrt" and self.gmrt_flagfiles:
            values["flagfile"] = [str(path) for path in self.gmrt_flagfiles]
        return values

    def _matching_report(
        self,
        report,
        source_identity,
        flag_identities,
        fits_format,
        task_name,
        telescope,
    ):
        return bool(
            isinstance(report, dict)
            and report.get("schema_version") == 1
            and report.get("stage") == "ingest"
            and report.get("status") == "completed"
            and report.get("source_fits") == source_identity
            and report.get("gmrt_flagfiles", []) == flag_identities
            and report.get("requested_fits_format") == self.requested_format
            and report.get("fits_format") == fits_format
            and report.get("casa_task") == task_name
            and report.get("requested_telescope", "auto")
            == self.requested_telescope
            and report.get("telescope") == telescope
            and report.get("measurement_set") == str(self.output_ms)
            and report.get("source_format", "fits")
            == ("lta" if self.input_is_lta else "fits")
        )

    def converted_fits_path(self):
        """Where the FITS converted from an LTA file is written."""

        converter = _load_lta_convert()
        return self.ingest_dir / "lta" / f"{converter.lta_stem(self.input_fits)}.FITS"

    def _lta_telescope(self):
        """LTA files come only from the GMRT; refuse any other telescope."""

        if self.requested_telescope not in ("auto", "ugmrt"):
            raise FITSDetectionError(
                f"--telescope {self.requested_telescope} conflicts with an "
                "LTA input: LTA is the GMRT's own format."
            )
        return "ugmrt"

    def _convert_lta(self, source_identity):
        """Convert the LTA input to UVFITS; return the conversion record."""

        if self.lta_converter is not None:
            convert = self.lta_converter
        else:
            convert = _load_lta_convert().convert
        record = convert(
            self.input_fits,
            self.converted_fits_path(),
            listscan=self.listscan,
            gvfits=self.gvfits,
            lta_identity=source_identity,
        )
        converted = Path(record["fits"]["path"])
        summary = inspect_fits_structure(converted)
        try:
            resolve_visibility_format(summary, "uvfits")
        except FITSDetectionError as exc:
            raise FITSDetectionError(
                f"gvfits did not write random-groups UVFITS ({converted}): {exc}"
            ) from exc
        self._resolve_telescope(summary)
        return record, converted, summary

    def _remove_converted_fits(self, record):
        """Remove the converted FITS and its record; keep the text logs."""

        path = Path(record["fits"]["path"])
        expected = self.converted_fits_path()
        if path != expected:
            raise RuntimeError(f"Refusing to remove unexpected file: {path}")
        converter = _load_lta_convert()
        for item in (path.with_name(path.name + converter.RECORD_SUFFIX), path):
            if item.is_file() and not item.is_symlink():
                item.unlink()

    @staticmethod
    def _remove_owned_staging(path, output_name):
        path = Path(path)
        prefix = f".{output_name}.partial.{os.getpid()}."
        if path.parent == path.parent.resolve() and path.name.startswith(prefix):
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
            elif path.exists() or path.is_symlink():
                path.unlink()
            return
        raise RuntimeError(f"Refusing to clean unexpected staging path: {path}")

    def _run_locked(self):
        kind = "LTA" if self.input_is_lta else "FITS"
        if not self.input_fits.is_file():
            raise FileNotFoundError(f"Input {kind} file not found: {self.input_fits}")
        if self.output_ms.is_symlink():
            raise RuntimeError(
                f"Refusing a symlink as the output Measurement Set: {self.output_ms}"
            )
        if self.output_ms == self.input_fits:
            raise RuntimeError(f"Input {kind} and output MS paths must be different.")

        initial_stat = self.input_fits.stat()
        if self.input_is_lta:
            # Nothing in an LTA file can be read before listscan has run, but
            # what it becomes is fixed: GMRT UVFITS, read by importgmrt. The
            # FITS is inspected once it exists.
            summary = None
            fits_format = "uvfits"
            telescope = self._lta_telescope()
        else:
            summary = inspect_fits_structure(self.input_fits)
            fits_format = resolve_visibility_format(summary, self.requested_format)
            telescope = self._resolve_telescope(summary)
        task_name = self._select_task(fits_format, telescope)
        source_identity = build_input_identity(self.input_fits)
        if (
            source_identity["size_bytes"] != initial_stat.st_size
            or source_identity["mtime_ns"] != initial_stat.st_mtime_ns
        ):
            raise RuntimeError(f"Input {kind} changed while it was being inspected.")

        flag_identities = [build_input_identity(path) for path in self.gmrt_flagfiles]
        if flag_identities and not (
            fits_format == "uvfits" and telescope == "ugmrt"
        ):
            raise ValueError(
                "GMRT flag files require UVFITS with --telescope ugmrt "
                "(or matching GMRT FITS telescope metadata)."
            )

        report = read_json(self.report_path, {})
        report_matches = self._matching_report(
            report,
            source_identity,
            flag_identities,
            fits_format,
            task_name,
            telescope,
        )

        if self.output_ms.exists() and not self.force:
            if report_matches and self.checkpoint_path.is_file():
                self.ms_validator(self.output_ms)
                print(f"Reusing validated imported Measurement Set: {self.output_ms}")
                return self.output_ms
            raise RuntimeError(
                "Imported Measurement Set output already exists without a "
                "matching completed ingest provenance record:\n"
                f"{self.output_ms}\n"
                "Choose another --imported-ms path or use --force after "
                "checking the existing data."
            )

        self.output_ms.parent.mkdir(parents=True, exist_ok=True)
        self.ingest_dir.mkdir(parents=True, exist_ok=True)
        if not self.output_ms.exists() and self.checkpoint_path.exists():
            # A completion marker cannot remain authoritative after its
            # materialized product has disappeared. Clearing only the marker
            # makes any failed re-import visibly incomplete while preserving
            # the old report for diagnosis.
            self.checkpoint_path.unlink()
        staging = self.output_ms.with_name(
            f".{self.output_ms.name}.partial.{os.getpid()}.{time.time_ns()}"
        )
        if staging.exists() or staging.is_symlink():
            raise RuntimeError(f"Ingest staging path already exists: {staging}")

        # The CASA task is loaded before an LTA conversion, which can take
        # hours, so that a missing task stops the run first.
        task = self.task_loader(task_name, required=True)
        archive = None
        conversion = None
        try:
            source = self.input_fits
            if self.input_is_lta:
                conversion, source, summary = self._convert_lta(source_identity)
            task(**self._task_kwargs(task_name, staging, source))
            self.ms_validator(staging)

            current_stat = self.input_fits.stat()
            if (
                current_stat.st_size != source_identity["size_bytes"]
                or current_stat.st_mtime_ns != source_identity["mtime_ns"]
            ):
                raise RuntimeError(f"Input {kind} changed while CASA was importing it.")
            for identity in flag_identities:
                flag_stat = Path(identity["path"]).stat()
                if (
                    flag_stat.st_size != identity["size_bytes"]
                    or flag_stat.st_mtime_ns != identity["mtime_ns"]
                ):
                    raise RuntimeError(
                        "A GMRT flag file changed while CASA was importing it."
                    )

            # The old checkpoint remains valid if a forced replacement fails
            # before this point. Remove it only when the new MS is ready to be
            # published, then publish report and a fresh marker in that order.
            if self.checkpoint_path.exists():
                self.checkpoint_path.unlink()

            if self.output_ms.exists():
                archive = _archive_name(self.output_ms)
                os.replace(self.output_ms, archive)
            try:
                os.replace(staging, self.output_ms)
            except Exception:
                if archive is not None and not self.output_ms.exists():
                    os.replace(archive, self.output_ms)
                    archive = None
                raise

            validation = self.ms_validator(self.output_ms)

            # The FITS headers state what should have been imported; the MS
            # states what was. Compare them before the ingest is accepted.
            structure = compare_fits_and_ms(summary, validation)
            if structure.get("status") == "inconsistent":
                print(
                    "WARNING: the imported Measurement Set does not match the "
                    "source FITS structure ("
                    + ", ".join(structure.get("inconsistent_checks", []))
                    + "). Compare ingest_report.json structure_comparison "
                    "against the archive product before using this data."
                )

            imported_names = (
                validation.get("observation", {}) or {}
            ).get("telescope_names", [])
            imported_telescopes = classify_telescope_names(imported_names)
            if (
                self.requested_telescope != "auto" and
                imported_telescopes and
                self.requested_telescope not in imported_telescopes
            ):
                raise MeasurementSetValidationError(
                    f"--telescope {self.requested_telescope} conflicts with "
                    "the telescope CASA wrote into the imported Measurement "
                    "Set: " + ", ".join(imported_names)
                )

            completed = {
                "schema_version": 1,
                "stage": "ingest",
                "status": "completed",
                "created_utc": utc_now(),
                # The positional input: for LTA input, the LTA file itself.
                "source_fits": source_identity,
                "source_format": "lta" if self.input_is_lta else "fits",
                "gmrt_flagfiles": flag_identities,
                "requested_fits_format": self.requested_format,
                "fits_format": fits_format,
                "requested_telescope": self.requested_telescope,
                "telescope": telescope,
                "casa_task": task_name,
                "measurement_set": str(self.output_ms),
                "fits_summary": summary,
                "validation": validation,
                "structure_comparison": structure,
                "imported_telescope_names": imported_names,
                "imported_telescopes": imported_telescopes,
            }
            if archive is not None:
                completed["replaced_measurement_set_archive"] = str(archive)
            if conversion is not None:
                completed["lta_conversion"] = conversion
                completed["converted_fits_kept"] = self.keep_converted_fits
            atomic_write_json(self.report_path, completed)
            atomic_write_checkpoint(self.checkpoint_path)
            if conversion is not None and not self.keep_converted_fits:
                # The Measurement Set now holds everything the FITS held. A
                # FITS left by a failed import is kept, so a retry reuses it.
                try:
                    self._remove_converted_fits(conversion)
                except OSError as exc:
                    print(f"WARNING: could not remove the converted FITS: {exc}")
            return self.output_ms
        finally:
            if staging.exists() or staging.is_symlink():
                self._remove_owned_staging(staging, self.output_ms.name)

    def run(self):
        """Acquire an output-scoped lock and execute one ingest transaction."""

        self.output_ms.parent.mkdir(parents=True, exist_ok=True)
        lock_path = self.output_ms.with_name(
            f".{self.output_ms.name}.ingest.lock"
        )
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        try:
            descriptor = os.open(lock_path, flags, 0o600)
        except FileExistsError as exc:
            raise RuntimeError(
                "Another FITS ingest appears to own the output lock:\n"
                f"{lock_path}\n"
                "Wait for it to finish. If no ingest process is running, "
                "inspect and remove this stale lock manually."
            ) from exc

        lock_identity = os.fstat(descriptor)
        try:
            payload = (
                f"pid={os.getpid()}\n"
                f"created_utc={utc_now()}\n"
                f"output_ms={self.output_ms}\n"
            ).encode("utf-8")
            os.write(descriptor, payload)
            os.fsync(descriptor)
            return self._run_locked()
        finally:
            os.close(descriptor)
            try:
                current_identity = os.stat(lock_path, follow_symlinks=False)
            except FileNotFoundError:
                current_identity = None
            if (
                current_identity is not None
                and current_identity.st_dev == lock_identity.st_dev
                and current_identity.st_ino == lock_identity.st_ino
            ):
                lock_path.unlink()


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Import random-groups UVFITS, FITS-IDI or GMRT LTA visibility "
            "data into a validated CASA Measurement Set."
        )
    )
    parser.add_argument(
        "input_fits",
        help="Input visibility FITS file, or a GMRT LTA file (.lta, .ltb, .lta_1 ...).",
    )
    parser.add_argument("output_ms", help="Output Measurement Set directory.")
    parser.add_argument("--outdir", default="Products", help="Pipeline output root.")
    parser.add_argument(
        "--fits-format",
        default="auto",
        choices=("auto", "uvfits", "fitsidi", "fits-idi", "lta"),
        help=(
            "Visibility FITS format; auto uses content-based detection, and "
            "a name ending in .lta, .ltb or .lta_N for GMRT LTA files."
        ),
    )
    parser.add_argument(
        "--telescope",
        default="auto",
        choices=tuple(sorted(TELESCOPE_CHOICES)),
        help="Telescope hint; uGMRT UVFITS uses CASA importgmrt.",
    )
    parser.add_argument(
        "--gmrt-flagfile",
        action="append",
        default=[],
        help="Optional GMRT online flag file (repeat for multiple files).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-import and preserve any existing output under a timestamped name.",
    )
    parser.add_argument(
        "--listscan",
        default="listscan",
        help="NCRA's listscan program (or a command running it), for LTA input.",
    )
    parser.add_argument(
        "--gvfits",
        default="gvfits",
        help="NCRA's gvfits program (or a command running it), for LTA input.",
    )
    parser.add_argument(
        "--keep-lta-fits",
        action="store_true",
        help="Keep the FITS converted from an LTA file after the import.",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        result = FITSIngestor(
            input_fits=args.input_fits,
            output_ms=args.output_ms,
            outdir=args.outdir,
            fits_format=args.fits_format,
            telescope=args.telescope,
            gmrt_flagfiles=args.gmrt_flagfile,
            force=args.force,
            listscan=args.listscan,
            gvfits=args.gvfits,
            keep_converted_fits=args.keep_lta_fits,
        ).run()
        print(f"Imported Measurement Set: {result}")
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
