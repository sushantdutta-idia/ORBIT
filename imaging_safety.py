#!/usr/bin/env python3
"""CASA-free helpers for final-imaging resume and cache safety."""

import hashlib
import math
import re
from pathlib import Path


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of failing partway through a
# stage. Increase it whenever pipeline.py starts to rely on something new here.
PIPELINE_INTERFACE_VERSION = 2


PBCOR_DEPENDENT_CHECKPOINTS = (
    "pbcor",
    "common_beam",
    "metrics",
    "qa_plots",
    "fits",
    "imaging",
)


# killMS writes its solution archive below the input MS even when the output
# solution name is supplied explicitly.  DD calibration uses this reserved,
# content-addressed namespace and records the solution file separately.  It is
# not part of the CASA table and must not make an otherwise immutable imaging
# MS appear to have changed after DD calibration.
_PIPELINE_KILLMS_ARTIFACT = re.compile(
    r"^killMS\.radio_continuum_dd_[0-9a-f]{64}\..+$"
)


def ugmrt_primary_beam_response(radial_argument, coefficients):
    """Evaluate the NCRA 2023 eighth-order uGMRT power-beam fit.

    ``radial_argument`` is angular separation in arcminutes multiplied by
    observing frequency in GHz.  It may be a scalar or a NumPy array.
    """

    values = list(coefficients or ())
    if len(values) != 4:
        raise ValueError("uGMRT primary-beam fit requires four coefficients.")
    a, b, c, d = (float(value) for value in values)
    x = radial_argument
    return (
        1.0 +
        (a / 1.0e3) * x ** 2 +
        (b / 1.0e7) * x ** 4 +
        (c / 1.0e10) * x ** 6 +
        (d / 1.0e13) * x ** 8
    )


def first_outward_beam_cutoff(coefficients, cutoff, step=0.001, limit=200.0):
    """Find the first main-lobe radius where the PB drops to *cutoff*.

    The fitted eighth-order polynomial is not globally monotonic and can turn
    upward well outside its measured main lobe.  Returning the first outward
    crossing prevents disconnected far-field polynomial lobes from being
    treated as usable primary-beam response.
    """

    threshold = float(cutoff)
    increment = float(step)
    maximum = float(limit)
    if not math.isfinite(threshold) or not 0.0 < threshold < 1.0:
        raise ValueError("Primary-beam cutoff must lie strictly between 0 and 1.")
    if not math.isfinite(increment) or increment <= 0.0:
        raise ValueError("Primary-beam cutoff search step must be positive.")

    previous_x = 0.0
    previous_value = float(ugmrt_primary_beam_response(0.0, coefficients))
    count = int(math.ceil(maximum / increment))
    for index in range(1, count + 1):
        current_x = index * increment
        current_value = float(
            ugmrt_primary_beam_response(current_x, coefficients)
        )
        if previous_value > threshold and current_value <= threshold:
            low, high = previous_x, current_x
            for _ in range(64):
                middle = (low + high) / 2.0
                value = float(
                    ugmrt_primary_beam_response(middle, coefficients)
                )
                if value > threshold:
                    low = middle
                else:
                    high = middle
            return high
        previous_x = current_x
        previous_value = current_value

    raise ValueError(
        "uGMRT primary-beam polynomial has no outward cutoff crossing "
        f"below radial argument {maximum}"
    )


def angular_separation_rad(ra1, dec1, ra2, dec2):
    """Great-circle separation between two directions, in radians."""

    delta_ra = ra2 - ra1
    numerator = math.hypot(
        math.cos(dec2) * math.sin(delta_ra),
        math.cos(dec1) * math.sin(dec2) -
        math.sin(dec1) * math.cos(dec2) * math.cos(delta_ra),
    )
    denominator = (
        math.sin(dec1) * math.sin(dec2) +
        math.cos(dec1) * math.cos(dec2) * math.cos(delta_ra)
    )
    return math.atan2(numerator, denominator)


def nearest_pointing_offset_arcmin(image_direction_rad, pointing_directions_rad):
    """Smallest separation between an image centre and any antenna pointing.

    A primary-beam model is centred on the image reference pixel.  That is
    only correct while the image reference direction is the direction the
    antennas were actually pointing at.  A re-phased Measurement Set breaks
    that assumption silently and mis-centres the beam, which biases every
    off-axis flux density.
    """

    try:
        ra, dec = (float(value) for value in image_direction_rad)
    except (TypeError, ValueError) as exc:
        raise ValueError("Image direction must be two numeric radians.") from exc

    offsets = []
    for pointing in pointing_directions_rad or ():
        try:
            pointing_ra, pointing_dec = (float(value) for value in pointing)
        except (TypeError, ValueError):
            continue
        if not (math.isfinite(pointing_ra) and math.isfinite(pointing_dec)):
            continue
        offsets.append(
            math.degrees(
                angular_separation_rad(ra, dec, pointing_ra, pointing_dec)
            ) * 60.0
        )

    if not offsets:
        raise ValueError("No usable pointing direction was supplied.")

    return min(offsets)


def normalize_frequency_hz(value, unit, converter):
    """Convert one native spectral coordinate to a finite positive Hz value.

    ``converter`` is deliberately injected so this module remains importable
    outside CASA. In production it is ``casatools.quanta.convert``.
    """

    numeric = float(value)
    native_unit = str(unit or "").strip()
    if not math.isfinite(numeric):
        raise ValueError("Spectral reference value is not finite.")
    if not native_unit:
        raise ValueError("Spectral coordinate has no native frequency unit.")

    converted = converter(
        {"value": numeric, "unit": native_unit},
        "Hz",
    )
    if not isinstance(converted, dict) or "value" not in converted:
        raise ValueError("Frequency converter did not return a quantity value.")

    frequency_hz = float(converted["value"])
    if not math.isfinite(frequency_hz) or frequency_hz <= 0.0:
        raise ValueError(
            "Spectral reference frequency must be finite and greater than zero."
        )
    return frequency_hz


def image_geometry_mismatch(source, candidate):
    """Return why two CASA image grids are not scientifically conformant.

    CASA may normalize a coordinate-system record when ``image.fromarray``
    writes it (record ordering, redundant WCS fields, or numeric rendering can
    change).  A byte hash of ``coordsys.torecord()`` is therefore a provenance
    identity, not a valid geometry-equivalence test.  This comparison uses the
    live grid quantities that determine pixel alignment while keeping strict,
    sub-pixel tolerances.
    """

    if not isinstance(source, dict) or not isinstance(candidate, dict):
        return "one or both geometry descriptors are not objects"

    for key in ("shape", "pixel_axes"):
        try:
            left = [int(value) for value in source.get(key, [])]
            right = [int(value) for value in candidate.get(key, [])]
        except (TypeError, ValueError):
            return f"{key} is not an integer sequence"
        if left != right or not left:
            return f"{key} differs: {left!r} != {right!r}"

    def numeric_sequence(key, expected_length=2):
        try:
            left = [float(value) for value in source.get(key, [])]
            right = [float(value) for value in candidate.get(key, [])]
        except (TypeError, ValueError):
            return None, None, f"{key} is not numeric"
        if len(left) != expected_length or len(right) != expected_length:
            return left, right, f"{key} does not contain {expected_length} values"
        if not all(math.isfinite(value) for value in left + right):
            return left, right, f"{key} contains a non-finite value"
        return left, right, None

    increments_a, increments_b, reason = numeric_sequence(
        "direction_increment_rad"
    )
    if reason:
        return reason
    for left, right in zip(increments_a, increments_b):
        if not math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-15):
            return (
                "direction increments differ: "
                f"{increments_a!r} != {increments_b!r}"
            )

    pixels_a, pixels_b, reason = numeric_sequence(
        "direction_reference_pixel"
    )
    if reason:
        return reason
    for left, right in zip(pixels_a, pixels_b):
        if not math.isclose(left, right, rel_tol=0.0, abs_tol=1e-7):
            return (
                "direction reference pixels differ: "
                f"{pixels_a!r} != {pixels_b!r}"
            )

    values_a, values_b, reason = numeric_sequence(
        "direction_reference_value_rad"
    )
    if reason:
        return reason
    longitude_delta = math.remainder(
        values_a[0] - values_b[0],
        2.0 * math.pi,
    )
    if abs(longitude_delta) > 1e-10 or not math.isclose(
        values_a[1],
        values_b[1],
        rel_tol=0.0,
        abs_tol=1e-10,
    ):
        return (
            "direction reference values differ: "
            f"{values_a!r} != {values_b!r}"
        )

    try:
        frequency_a = float(source.get("frequency_hz"))
        frequency_b = float(candidate.get("frequency_hz"))
    except (TypeError, ValueError):
        return "spectral reference frequency is not numeric"
    if not (
        math.isfinite(frequency_a) and
        math.isfinite(frequency_b) and
        frequency_a > 0.0 and
        frequency_b > 0.0
    ):
        return "spectral reference frequency is invalid"
    if not math.isclose(
        frequency_a,
        frequency_b,
        rel_tol=1e-12,
        abs_tol=1e-3,
    ):
        return (
            "spectral reference frequencies differ: "
            f"{frequency_a!r} != {frequency_b!r} Hz"
        )

    return None


def image_geometries_match(source, candidate):
    """Return whether two live CASA image descriptors share one pixel grid."""

    return image_geometry_mismatch(source, candidate) is None


def image_root_partial_patterns(imagename):
    """Return force-cleanup globs scoped to one final-image basename."""

    name = Path(imagename).name
    return (
        f".{name}*.partial.*",
        f".{name}*.tmp.*",
    )


def path_stat_fingerprint(path):
    """Bind provenance to a materialized CASA table without hashing pixels.

    CASA images can be many gigabytes.  A stable inventory of relative file
    names, sizes, and nanosecond modification times catches replacement or
    mutation while ignoring transient table-lock files.
    """

    root = Path(path).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"CASA table directory is missing: {root}")

    digest = hashlib.sha256()
    file_count = 0
    for item in sorted(root.rglob("*"), key=lambda value: str(value)):
        if not item.is_file() or item.name.endswith(".lock"):
            continue
        relative = item.relative_to(root)
        if (
            len(relative.parts) == 1 and
            _PIPELINE_KILLMS_ARTIFACT.fullmatch(item.name)
        ):
            continue
        stat = item.stat()
        digest.update(str(relative).encode("utf-8"))
        digest.update(str(stat.st_size).encode("ascii"))
        digest.update(str(stat.st_mtime_ns).encode("ascii"))
        file_count += 1

    if file_count == 0:
        raise ValueError(f"CASA table contains no durable files: {root}")
    return digest.hexdigest()


def casa_table_content_identity(path, chunk_size=8 * 1024 * 1024):
    """Identify a CASA table directory by its file contents.

    ``path_stat_fingerprint`` binds file sizes and modification times. A
    product fingerprinted straight after it is written can later read back
    with different modification times although nothing wrote to it. That
    happened to a uGMRT primary-beam image on a network filesystem between
    the moment imaging recorded it and the orchestrator's check. Contents do
    not change that way, so image products are bound by content. Measurement
    Sets are far too large to hash on every check and keep their stat
    fingerprint.
    """

    root = Path(path).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"CASA table directory is missing: {root}")

    digest = hashlib.sha256()
    files = []
    for item in sorted(root.rglob("*"), key=lambda value: str(value)):
        if not item.is_file() or item.name.endswith(".lock"):
            continue
        relative = item.relative_to(root)
        if (
            len(relative.parts) == 1 and
            _PIPELINE_KILLMS_ARTIFACT.fullmatch(item.name)
        ):
            continue
        file_digest = hashlib.sha256()
        with item.open("rb") as handle:
            for chunk in iter(lambda: handle.read(chunk_size), b""):
                file_digest.update(chunk)
        size = int(item.stat().st_size)
        files.append({
            "path": str(relative),
            "size_bytes": size,
            "sha256": file_digest.hexdigest(),
        })
        digest.update(
            f"{relative}\0{size}\0{file_digest.hexdigest()}\0".encode("utf-8")
        )

    if not files:
        raise ValueError(f"CASA table contains no durable files: {root}")
    return {"schema_version": 1, "sha256": digest.hexdigest(), "files": files}


# A CASA table's metadata files, as opposed to the files holding its data.
# CASA rewrites table.dat whenever any keyword changes -- units, image info, a
# mask entry, the coordinate record -- and the history subtable records what
# has opened the image. Neither says anything about the pixels. A primary beam
# whose table.dat was rewritten while its pixels stayed identical is still the
# beam that was applied. Coordinate metadata is covered separately: the
# geometry is recorded in provenance, and imaging re-reads it live and
# compares it with the restored image before reusing a beam.
_METADATA_MEMBERS = frozenset({"table.dat", "table.info", "table.lock"})
_HISTORY_SUBTABLE = "logtable"


def _is_data_member(relative):
    parts = str(relative).split("/")
    if parts[0] == _HISTORY_SUBTABLE:
        return False
    return parts[-1] not in _METADATA_MEMBERS


def data_identity(identity):
    """Reduce a content identity to the files that hold the data."""

    digest = hashlib.sha256()
    members = {}
    for item in (identity or {}).get("files", []):
        if not isinstance(item, dict) or not _is_data_member(item.get("path")):
            continue
        members[item["path"]] = item
        digest.update(
            f"{item['path']}\0{item.get('size_bytes')}\0"
            f"{item.get('sha256')}\0".encode("utf-8")
        )
    return digest.hexdigest(), members


def describe_content_change(recorded, live, data_only=False):
    """Name the files whose contents differ between two content identities."""

    if data_only:
        before = data_identity(recorded)[1]
        after = data_identity(live)[1]
    else:
        before = {
            item.get("path"): item
            for item in (recorded or {}).get("files", [])
            if isinstance(item, dict)
        }
        after = {item["path"]: item for item in (live or {}).get("files", [])}

    parts = []
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    changed = sorted(
        name
        for name in set(before) & set(after)
        if before[name].get("sha256") != after[name].get("sha256")
    )
    if added:
        parts.append("added " + ", ".join(added))
    if removed:
        parts.append("removed " + ", ".join(removed))
    if changed:
        parts.append("rewritten " + ", ".join(changed))
    return "; ".join(parts) or "the overall content hash differs"


def artifact_matches_provenance(provenance, prefix, path):
    """Whether an image still matches the identity recorded for it.

    Provenance written with content identities (``<prefix>_content_identity``)
    is compared over the files that hold the data. Older provenance carries
    only ``<prefix>_stat_sha256`` and is compared by that. Returns
    ``(True, "")`` or ``(False, reason)``; a missing or unreadable image
    raises.
    """

    provenance = provenance if isinstance(provenance, dict) else {}
    recorded = provenance.get(f"{prefix}_content_identity")
    if isinstance(recorded, dict) and recorded.get("files"):
        live = casa_table_content_identity(path)
        recorded_data, recorded_members = data_identity(recorded)
        if not recorded_members:
            # Nothing was recorded that holds data, so there is nothing to
            # compare but the whole table. Fail closed on any difference.
            if live.get("sha256") == recorded.get("sha256"):
                return True, ""
            return False, (
                "its contents changed: " +
                describe_content_change(recorded, live)
            )
        if data_identity(live)[0] == recorded_data:
            return True, ""
        return False, (
            "its pixels changed: " +
            describe_content_change(recorded, live, data_only=True)
        )

    if provenance.get(f"{prefix}_stat_sha256") == path_stat_fingerprint(path):
        return True, ""
    return False, (
        "its file sizes or modification times differ from provenance recorded "
        "before content identities were kept"
    )


def measurement_set_instance_identity(ms_path):
    """Bind an MS directory instance without hashing mutable MAIN columns.

    FLAG, MODEL_DATA, and CORRECTED_DATA legitimately change throughout a
    reduction.  The MS root inode plus required immutable schema subtables
    still detects a same-path table swap while remaining stable across those
    intended visibility writes.
    """

    path = Path(ms_path).expanduser().resolve()
    if not (path.is_dir() and (path / "table.dat").is_file()):
        raise ValueError(f"Not a CASA Measurement Set table: {path}")

    stable_subtables = (
        "ANTENNA",
        "DATA_DESCRIPTION",
        "OBSERVATION",
        "POLARIZATION",
        "SPECTRAL_WINDOW",
    )
    root_stat = path.stat()
    digest = hashlib.sha256()
    digest.update(
        repr((int(root_stat.st_dev), int(root_stat.st_ino))).encode("utf-8")
    )
    members = 0

    for subtable_name in stable_subtables:
        subtable = path / subtable_name
        if not (subtable.is_dir() and (subtable / "table.dat").is_file()):
            raise ValueError(
                "Measurement Set is missing required subtable: "
                f"{subtable}"
            )
        subtable_stat = subtable.stat()
        digest.update(
            repr(
                (
                    subtable_name,
                    int(subtable_stat.st_dev),
                    int(subtable_stat.st_ino),
                )
            ).encode("utf-8")
        )
        members += 1

        for filename in ("table.dat", "table.info"):
            member = subtable / filename
            if not member.is_file():
                if filename == "table.info":
                    continue
                raise ValueError(
                    "Measurement Set subtable is incomplete: "
                    f"{member}"
                )
            stat = member.stat()
            digest.update(
                repr(
                    (
                        subtable_name,
                        filename,
                        int(stat.st_dev),
                        int(stat.st_ino),
                        int(stat.st_size),
                    )
                ).encode("utf-8")
            )
            if stat.st_size <= 4 * 1024 * 1024:
                with member.open("rb") as handle:
                    for chunk in iter(
                        lambda: handle.read(1024 * 1024), b""
                    ):
                        digest.update(chunk)
            members += 1

    return {
        "schema_version": 1,
        "root_device": int(root_stat.st_dev),
        "root_inode": int(root_stat.st_ino),
        "bound_members": members,
        "sha256": digest.hexdigest(),
    }


def file_content_identity(path, chunk_size=8 * 1024 * 1024):
    """Return a SHA-256 identity for one published regular-file product."""

    product = Path(path).expanduser().resolve()
    if not product.is_file():
        raise ValueError(f"Published file is missing: {product}")

    digest = hashlib.sha256()
    with product.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)

    stat = product.stat()
    return {
        "path": str(product),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "sha256": digest.hexdigest(),
    }


def validate_file_content_identities(paths, recorded):
    """Validate an exact published-file inventory against recorded hashes."""

    expected_paths = [
        str(Path(path).expanduser().resolve())
        for path in (paths or [])
    ]
    if len(expected_paths) != len(set(expected_paths)):
        return False, "published FITS inventory contains duplicate paths"
    if not isinstance(recorded, dict):
        return False, "published FITS content identities are missing"
    if set(recorded) != set(expected_paths):
        return False, "published FITS content-identity inventory is mismatched"

    for product in expected_paths:
        expected = recorded.get(product)
        if not isinstance(expected, dict):
            return False, f"published FITS identity is missing: {product}"
        try:
            actual = file_content_identity(product)
        except (OSError, ValueError) as exc:
            return False, str(exc)
        # Modification times are recorded for reference only. They can read
        # back differently on a network filesystem with no write in between;
        # the path, size and SHA-256 are what bind the product.
        if any(
            expected.get(key) != actual[key]
            for key in ("path", "size_bytes", "sha256")
        ):
            return False, f"published FITS content changed: {product}"

    return True, "published FITS content identities match"
