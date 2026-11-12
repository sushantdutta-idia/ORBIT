#!/usr/bin/env python3
"""Convert a GMRT LTA file to random-groups UVFITS with listscan and gvfits.

The GMRT Online Archive delivers raw visibilities in GMRT's own LTA format.
NCRA's two programs convert them: ``listscan FILE.lta`` writes ``FILE.log``,
a text list of the file's scans that names ``TEST.FITS`` as the output, and
``gvfits FILE.log`` writes that FITS file into its working directory. CASA's
``importgmrt`` then reads the FITS; ``fits_ingest.py`` does that part.

Both programs run in a new, empty staging folder on a link to the LTA file,
so nothing is written beside the archive file, a failed conversion leaves no
half-written FITS, and two conversions cannot meet in one folder. The log is
used exactly as listscan wrote it: every scan, antenna and channel goes into
the FITS, because the pipeline flags data and chooses scans itself. The FITS
gets its final name only after gvfits has finished, together with a record
of the LTA file it came from, so a later run reuses it instead of converting
again.

The programs are not part of CASA. NCRA distributes them, and the
CAPTURE-CASA6 repository (https://github.com/ruta-k/CAPTURE-CASA6) carries
Linux builds of both. Give their paths with ``--listscan`` and ``--gvfits``,
or put them on PATH. A value with spaces is a whole command, for example one
that runs the program in a container; give absolute paths in it, because
the programs run inside the staging folder.

Only the standard library is used, so this runs under CASA's Python or any
other Python 3.
"""

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


# 1: first release. 2: libraries beside the programs (or in lib/ there) are
# found first; a missing library is explained.
PIPELINE_INTERFACE_VERSION = 2

# ".lta" holds GWB data or the upper sideband of GSB data, ".ltb" the GSB
# lower sideband; a long recording may be split into ".lta_1", ".lta_2" ...
LTA_SUFFIX = re.compile(r"\.lt[ab](_\d+)?$", re.IGNORECASE)
TOOL_SOURCE = "https://github.com/ruta-k/CAPTURE-CASA6"
RECORD_SUFFIX = ".conversion.json"
FITS_NAME = re.compile(r"[^\s'\"=,;]+\.fits\b", re.IGNORECASE)
# The FITS holds the visibilities with a weight for each, so it is at least
# as large as the LTA file; an MS of similar size follows it.
SPACE_WARNING_FACTOR = 3.0


class LTAConversionError(RuntimeError):
    """Raised when listscan or gvfits did not produce a usable FITS file."""


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def is_lta_path(path):
    """True for a name ending in .lta, .ltb or a numbered part such as .lta_1."""

    return bool(LTA_SUFFIX.search(Path(path).name))


def lta_stem(path):
    """A plain name for the products of one LTA file.

    ``X.lta`` gives ``X``; the other suffixes stay in the name, so ``X.ltb``
    gives ``X_ltb`` and ``X.lta_1`` gives ``X_lta_1`` and no two parts of one
    observation share a product name.
    """

    name = Path(path).name
    match = LTA_SUFFIX.search(name)
    base = name[:match.start()] if match else name
    suffix = match.group(0).lower() if match else ""
    if suffix and suffix != ".lta":
        base = f"{base}_{suffix[1:]}"
    clean = "".join(
        character if character.isalnum() or character in "-_" else "_"
        for character in base
    ).strip("_")
    return clean or "observation"


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def file_identity(path, checksum=True):
    path = Path(path).expanduser().resolve()
    stat = path.stat()
    identity = {
        "path": str(path),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }
    if checksum:
        identity["sha256"] = sha256_file(path)
    return identity


def read_json(path, default=None):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError, TypeError):
        return default


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


def tool_command(value, name):
    """The command that runs one of the two programs, as a list.

    The program is resolved to an absolute path now, because it will be
    started from inside the staging folder.
    """

    if isinstance(value, (list, tuple)):
        parts = [str(part) for part in value]
    else:
        parts = shlex.split(str(value or name))
    if not parts:
        raise ValueError(f"No command was given for {name}.")
    program = parts[0]
    found = shutil.which(program)
    if found is None:
        candidate = Path(program).expanduser()
        if candidate.is_file():
            raise FileNotFoundError(
                f"{name} is not executable: {candidate}. Run "
                f"'chmod +x {candidate}' and try again."
            )
        raise FileNotFoundError(
            f"{name} was not found ({program}). LTA files are converted "
            f"with NCRA's listscan and gvfits; Linux builds of both are in "
            f"{TOOL_SOURCE}. Download them, make them executable "
            f"(chmod +x), and give their paths with --listscan and --gvfits "
            f"(gmrt_listscan and gmrt_gvfits in config.yaml)."
        )
    parts[0] = str(Path(found).resolve())
    return parts


def tool_identity(command):
    identity = {"command": list(command)}
    try:
        identity["sha256"] = sha256_file(command[0])
    except OSError:
        pass
    return identity


def tool_library_dirs(command):
    """Folders searched first for the programs' shared libraries.

    The NCRA builds are old: listscan needs libgfortran.so.3, which current
    systems and CASA containers no longer have. A copy put beside the
    program, or in a lib/ folder beside it, is found there.
    """

    folders = []
    for part in command:
        path = Path(part)
        if not path.is_absolute() or not path.is_file():
            continue
        for folder in (path.parent, path.parent / "lib"):
            if folder.is_dir() and str(folder) not in folders:
                folders.append(str(folder))
    return folders


def tool_environment(command, environment=None):
    environment = dict(os.environ if environment is None else environment)
    folders = tool_library_dirs(command)
    if folders:
        existing = environment.get("LD_LIBRARY_PATH", "")
        environment["LD_LIBRARY_PATH"] = ":".join(
            folders + [entry for entry in existing.split(":") if entry]
        )
    return environment


def run_tool(command, cwd, transcript):
    """Run one program with no terminal input; return its exit status."""

    with open(transcript, "wb") as output:
        completed = subprocess.run(
            command,
            cwd=str(cwd),
            stdin=subprocess.DEVNULL,
            stdout=output,
            stderr=subprocess.STDOUT,
            env=tool_environment(command),
        )
    return completed.returncode


MISSING_LIBRARY = re.compile(r"error while loading shared libraries: (\S+?):")


def failure_message(name, status, transcript, command):
    """What a failed program said, and what to do about a missing library."""

    text = tail(transcript)
    message = f"{name} exited with status {status}. Its last lines:\n{text}"
    missing = MISSING_LIBRARY.search(text)
    if missing:
        folder = Path(command[0]).parent / "lib"
        message += (
            f"\n{name} needs {missing.group(1)}, which this system does not "
            f"have (inside a CASA container, that container). Put a copy of "
            f"it in {folder}/ (or beside {name}); the converters look there "
            "first. libgfortran.so.3 is in conda-forge's libgfortran 3.0.0 "
            "package: https://conda.anaconda.org/conda-forge/linux-64/"
            "libgfortran-3.0.0-1.tar.bz2"
        )
    return message


def tail(path, lines=15):
    try:
        text = Path(path).read_text(errors="replace")
    except OSError:
        return ""
    return "\n".join(text.splitlines()[-lines:])


def fits_names_in_log(text):
    """The FITS file names a listscan log mentions, in order, once each."""

    names = []
    for match in FITS_NAME.finditer(text):
        if match.group(0) not in names:
            names.append(match.group(0))
    return names


def find_log(staging, link, before):
    expected = staging / f"{link.name.split('.')[0]}.log"
    if expected.is_file():
        return expected
    created = sorted(
        path for path in staging.iterdir()
        if path.name not in before and path.suffix.lower() == ".log"
    )
    if len(created) == 1:
        return created[0]
    return None


def find_fits(staging, named, before):
    for name in named:
        candidate = staging / name
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    created = sorted(
        path for path in staging.iterdir()
        if path.name not in before
        and path.suffix.lower() == ".fits"
        and path.is_file()
        and path.stat().st_size > 0
    )
    if len(created) == 1:
        return created[0]
    if len(created) > 1:
        raise LTAConversionError(
            "gvfits wrote more than one FITS file ("
            + ", ".join(path.name for path in created)
            + "); the log should name one output."
        )
    return None


def reusable_record(record_path, output_fits, lta_identity):
    """The record of an earlier conversion of this very LTA file, if valid."""

    record = read_json(record_path)
    if not isinstance(record, dict) or record.get("schema_version") != 1:
        return None
    if record.get("lta") != lta_identity:
        return None
    fits = record.get("fits") or {}
    try:
        stat = Path(output_fits).stat()
    except OSError:
        return None
    if (
        fits.get("path") != str(output_fits)
        or fits.get("size_bytes") != stat.st_size
        or fits.get("mtime_ns") != stat.st_mtime_ns
    ):
        return None
    return record


def remove_staging(path, output_name):
    path = Path(path)
    prefix = f".{output_name}.partial.{os.getpid()}."
    if not path.name.startswith(prefix) or path.is_symlink():
        raise RuntimeError(f"Refusing to clean unexpected staging path: {path}")
    # rmtree removes the link to the LTA file, never the file itself.
    shutil.rmtree(path)


def convert(
    lta,
    output_fits,
    listscan="listscan",
    gvfits="gvfits",
    reuse=True,
    lta_identity=None,
    echo=print,
):
    """Convert ``lta`` to ``output_fits``; return the conversion record.

    ``lta_identity`` may carry the LTA file's identity with its SHA-256 when
    the caller has already computed it (hashing a large file takes minutes).
    """

    lta = Path(lta).expanduser().resolve()
    raw_output = Path(output_fits).expanduser()
    if not raw_output.is_absolute():
        raw_output = Path.cwd() / raw_output
    output_fits = raw_output.parent.resolve() / raw_output.name
    if not lta.is_file():
        raise FileNotFoundError(f"LTA file not found: {lta}")
    if output_fits == lta:
        raise ValueError("The FITS output must not be the LTA file itself.")

    listscan_command = tool_command(listscan, "listscan")
    gvfits_command = tool_command(gvfits, "gvfits")
    identity = lta_identity or file_identity(lta)
    record_path = output_fits.with_name(output_fits.name + RECORD_SUFFIX)

    if reuse:
        record = reusable_record(record_path, output_fits, identity)
        if record is not None:
            echo(f"Reusing the FITS converted earlier from this LTA file: {output_fits}")
            return dict(record, reused=True)

    output_fits.parent.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(output_fits.parent).free
    size = identity["size_bytes"]
    if free < size:
        raise LTAConversionError(
            f"Only {free / 1e9:.1f} GB is free in {output_fits.parent}; the "
            f"FITS from this {size / 1e9:.1f} GB LTA file will not fit."
        )
    if free < SPACE_WARNING_FACTOR * size:
        echo(
            f"WARNING: {free / 1e9:.1f} GB free in {output_fits.parent}. The "
            f"FITS and then the Measurement Set made from this "
            f"{size / 1e9:.1f} GB LTA file may need about "
            f"{SPACE_WARNING_FACTOR * size / 1e9:.0f} GB together."
        )

    staging = output_fits.with_name(
        f".{output_fits.name}.partial.{os.getpid()}.{time.time_ns()}"
    )
    staging.mkdir()
    try:
        match = LTA_SUFFIX.search(lta.name)
        suffix = match.group(0) if match else ".lta"
        link = staging / f"{lta_stem(lta)}{suffix}"
        os.symlink(lta, link)
        before = set(os.listdir(staging))

        echo(f"listscan {lta.name}")
        started = time.monotonic()
        status = run_tool(
            listscan_command + [link.name], staging, staging / "listscan.out"
        )
        listscan_seconds = time.monotonic() - started
        if status != 0:
            raise LTAConversionError(failure_message(
                "listscan", status, staging / "listscan.out", listscan_command
            ))
        log = find_log(staging, link, before)
        if log is None:
            beside = lta.with_name(f"{lta.name.split('.')[0]}.log")
            note = (
                f" A log beside the LTA file exists ({beside}); this listscan "
                "may write there instead of its working folder."
                if beside.exists() else ""
            )
            raise LTAConversionError(
                "listscan finished but wrote no .log file." + note
                + " Its last lines:\n" + tail(staging / "listscan.out")
            )
        log_text = log.read_text(errors="replace")
        named = fits_names_in_log(log_text)
        outside = [name for name in named if "/" in name or "\\" in name]
        if outside:
            raise LTAConversionError(
                "listscan's log names a FITS output outside the working "
                f"folder ({', '.join(outside)}); gvfits was not run."
            )

        echo(f"gvfits {log.name}  (this takes a while for a large file)")
        before = set(os.listdir(staging))
        started = time.monotonic()
        status = run_tool(
            gvfits_command + [log.name], staging, staging / "gvfits.out"
        )
        gvfits_seconds = time.monotonic() - started
        if status != 0:
            raise LTAConversionError(failure_message(
                "gvfits", status, staging / "gvfits.out", gvfits_command
            ))
        fits = find_fits(staging, named, before)
        if fits is None:
            raise LTAConversionError(
                "gvfits finished but wrote no FITS file"
                + (f" (the log names {', '.join(named)})" if named else "")
                + ". Its last lines:\n" + tail(staging / "gvfits.out")
            )

        now = lta.stat()
        if (
            now.st_size != identity["size_bytes"]
            or now.st_mtime_ns != identity["mtime_ns"]
        ):
            raise LTAConversionError(
                "The LTA file changed while it was being converted."
            )

        if record_path.exists():
            record_path.unlink()
        os.replace(fits, output_fits)
        kept = {}
        for source, ending in (
            (log, ".listscan.log"),
            (staging / "listscan.out", ".listscan.out"),
            (staging / "gvfits.out", ".gvfits.out"),
            # gvfits's own record of what it did ("See gvfits.log").
            (staging / "gvfits.log", ".gvfits.log"),
        ):
            if not source.is_file():
                continue
            destination = output_fits.with_name(output_fits.name + ending)
            os.replace(source, destination)
            kept[ending[1:]] = str(destination)

        record = {
            "schema_version": 1,
            "created_utc": utc_now(),
            "lta": identity,
            "fits": file_identity(output_fits, checksum=False),
            "fits_name_in_log": fits.name,
            "listscan": tool_identity(listscan_command),
            "gvfits": tool_identity(gvfits_command),
            "kept_files": kept,
            "seconds": {
                "listscan": round(listscan_seconds, 1),
                "gvfits": round(gvfits_seconds, 1),
            },
        }
        atomic_write_json(record_path, record)
        echo(f"Converted {lta.name} -> {output_fits}")
        return dict(record, reused=False)
    finally:
        if staging.exists():
            remove_staging(staging, output_fits.name)


def load_fits_ingest():
    path = Path(__file__).resolve().with_name("fits_ingest.py")
    spec = importlib.util.spec_from_file_location("lta_convert_fits_ingest", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Convert a GMRT LTA file to random-groups UVFITS with NCRA's "
            "listscan and gvfits."
        )
    )
    parser.add_argument("lta", help="Input LTA file (.lta, .ltb, .lta_1 ...).")
    parser.add_argument("output_fits", help="Output FITS file.")
    parser.add_argument("--listscan", default="listscan", help="listscan program or command.")
    parser.add_argument("--gvfits", default="gvfits", help="gvfits program or command.")
    parser.add_argument(
        "--no-reuse",
        action="store_true",
        help="Convert again even if this LTA file was converted to this FITS before.",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        record = convert(
            args.lta,
            args.output_fits,
            listscan=args.listscan,
            gvfits=args.gvfits,
            reuse=not args.no_reuse,
        )
        fits_ingest = load_fits_ingest()
        summary = fits_ingest.inspect_fits_structure(record["fits"]["path"])
        fits_ingest.resolve_visibility_format(summary, "uvfits")
        print(
            f"FITS {record['fits']['path']} is random-groups UVFITS "
            f"({', '.join(summary.get('telescope_header_values', [])) or 'telescope not named'})."
        )
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
