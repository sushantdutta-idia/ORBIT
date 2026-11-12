#!/usr/bin/env python3
# ============================================================
# pipeline.py
#
# Orchestrator for uGMRT, MeerKAT, and LOFAR continuum reduction.
#
# The individual stage scripts remain the authoritative implementations. This file only builds and runs them in the correct order.
#
# Standard order:
#   fits_ingest.py (visibility FITS or GMRT LTA only) -> inspect.py -> flagging.py
#   -> calibration.py -> selfcal.py
#   -> imaging.py -> source_detection.py
#
# DD calibration is intentionally opt-in:
#   -> ddcal.py
#
# Author: Sushant Dutta
# ============================================================

import os
import sys

SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

# Avoid shadowing Python's standard-library inspect module with
# the local inspect.py when argparse/dataclasses are imported.
if (
    sys.path and
    os.path.abspath(sys.path[0]) == SCRIPT_DIR
):

    sys.path.pop(0)

import argparse
import contextlib
import fcntl
import hashlib
import importlib.util
import json
import math
import platform
import re
import signal
import shlex
import shutil
import socket
import subprocess
import threading
import time
from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# CONFIG LOADING
# ============================================================

def load_config_module():

    config_path = Path(SCRIPT_DIR) / "config.py"

    if not config_path.exists():

        return None

    spec = importlib.util.spec_from_file_location(
        "ugmrt_pipeline_config",
        str(config_path),
    )

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


def load_telescope_profiles_module():

    module_path = Path(SCRIPT_DIR) / "telescope_profiles.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_telescope_profiles",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_imaging_safety_module():

    module_path = Path(SCRIPT_DIR) / "imaging_safety.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_imaging_safety",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_memory_diagnostics_module():

    module_path = Path(SCRIPT_DIR) / "memory_diagnostics.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_memory_diagnostics",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_fits_ingest_module():

    module_path = Path(SCRIPT_DIR) / "fits_ingest.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_fits_ingest",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONFIG = load_config_module()
TELESCOPE_PROFILE_MODULE = load_telescope_profiles_module()
IMAGING_SAFETY = load_imaging_safety_module()
MEMORY_DIAGNOSTICS = load_memory_diagnostics_module()

# Stages that run inside CASA, and so read casacore's resource settings.
CASA_STAGES = frozenset({
    "inspection",
    "flagging",
    "calibration",
    "selfcal",
    "imaging",
})

# Stages producing one set of products per science target. The stages before
# them describe the observation as a whole and run once for all targets.
PER_TARGET_STAGES = ("selfcal", "imaging", "source_detection", "ddcal")

# Reports describing the observation rather than one target. Under a target
# scope every other report lives in that target's own product directory.
SHARED_REPORTS = frozenset({
    "ingest",
    "inspection",
    "flagging_pre",
    "flagging_post",
    "bad_antennas",
    "recommended_refant",
    "calibration",
    "split_targets",
})


class StageInterrupted(Exception):
    """A stage was stopped from outside: cancellation, time limit, session end.

    Distinct from a stage failure, and deliberately not a RuntimeError: the
    per-target loop continues after a target fails, but an interruption
    applies to the whole run, so it must not be caught there.
    """


# Conventional exit code for SIGINT, which the stage scripts use when they
# are interrupted.
STAGE_INTERRUPTED_EXIT_CODE = 130


def target_directory_name(target):
    """A filesystem-safe directory name for a science target."""

    safe = "".join(
        character if character.isalnum() or character in "+-._" else "_"
        for character in str(target).strip()
    )
    return safe.strip("._") or "target"

# Files whose current behaviour this orchestrator relies on, with the minimum
# PIPELINE_INTERFACE_VERSION each must declare. Files are copied to the
# cluster one at a time, and a stale copy otherwise shows up only as a failure
# from old code partway through a stage.
# Calibration settings that decide what happens when a solve cannot use an
# antenna at all. They govern failure handling, not the science: an antenna
# they exclude is recorded as such, and every solve it took part in is redone
# without it. So a change to them alone must not make a calibration attempt,
# or its finished products, look like a different science recipe.
CALIBRATION_FAILURE_POLICY_DEFAULTS = {
    "auto_quarantine_failed_delay_antennas": True,
    "max_auto_quarantine_delay_antennas": 1,
    "max_auto_quarantine_delay_fraction": 0.10,
    "min_antennas_after_delay_quarantine": 8,
}

SOURCE_INTERFACE_VERSIONS = {
    # 3: flags written by the provisional (round-1) solutions are released
    # before the round-2 solve instead of being kept for good.
    # 4: calibrated targets averaged to the band's widest imaging channel.
    # 5: isolated failed solutions bridged in the tables targets receive.
    "calibration.py": 5,
    # 2: the PSF-aware auto-multithresh threshold and zero-mask test that
    # selfcal.py and imaging.py call.
    "casa_compat.py": 2,
    # 2: MeerKAT's WSClean settings. Without them MeerKAT would be imaged as
    # one flat-spectrum channel, so a stale copy must stop the run.
    # 3: self-calibration coverage limits, whole-image deconvolution and the
    # uGMRT Band 4 field. A stale copy would self-calibrate nothing again.
    # 4: Band 4's target channel width (target_max_channel_width_khz) and
    # the self-calibration major-cycle bound (clean_max_major_cycles).
    # 5: uGMRT target-transfer bridging, and the config.yaml switches for
    # amplitude self-calibration and self-calibration scales.
    "config.py": 5,
    "ddcal.py": 1,
    "flagging.py": 1,
    # 7: the masked final-image clean is bounded by major cycles, as
    # self-calibration's is.
    "imaging.py": 7,
    "imaging_safety.py": 2,
    # 3: the gain-outlier summary is no longer lost to numpy name arrays.
    "qa_products.py": 3,
    # 8: the masked self-calibration clean is bounded by major cycles
    # (clean_max_major_cycles) and is not continued past that bound.
    # 9: the report's final_metrics are measured again.
    # 10: amplitude self-calibration is tested on an image before it is
    # kept, and residual flagging subtracts an accepted loop's model.
    # 11: the amplitude solve's threshold, gain type and reference antennas.
    "selfcal.py": 12,
    "telescope_profiles.py": 1,
}

# Needed only when the run can use WSClean, and only the second of them when
# it can be MeerKAT: a uGMRT run, imaged with CASA, never loads either.
WSCLEAN_SOURCE_INTERFACE_VERSIONS = {"wsclean_imaging.py": 1}
# 1: GMRT LTA input, converted with listscan and gvfits before importgmrt.
LTA_SOURCE_INTERFACE_VERSIONS = {"fits_ingest.py": 1, "lta_convert.py": 2}
# The GMRT's own raw format: .lta (GWB, or the GSB upper sideband), .ltb (the
# GSB lower sideband), and numbered parts of a split recording (.lta_1 ...).
LTA_NAME_PATTERN = re.compile(r"\.lt[ab](_\d+)?$", re.IGNORECASE)
MEERKAT_WSCLEAN_SOURCE_INTERFACE_VERSIONS = {"meerkat_beam.py": 1}


def outdated_source_files(root, required=None):
    """List pipeline files older than this orchestrator expects."""

    problems = []
    for name, needed in (required or SOURCE_INTERFACE_VERSIONS).items():
        path = Path(root) / name
        try:
            lines = path.read_text().splitlines()
        except OSError:
            problems.append(f"{name} is missing")
            continue
        declared = None
        for line in lines:
            key, _, value = line.partition("=")
            if key.strip() == "PIPELINE_INTERFACE_VERSION" and not line[:1].isspace():
                try:
                    declared = int(value.strip())
                except ValueError:
                    declared = None
                break
        if declared is None or declared < needed:
            problems.append(
                f"{name} is from an earlier release "
                f"(interface {declared if declared is not None else 'unmarked'}, "
                f"needs {needed})"
            )
    return problems
FITS_INGEST = load_fits_ingest_module()

if CONFIG is not None:

    PIPELINE_STAGES = CONFIG.stage_order(
        include_ddcal=True
    )
    DEFAULT_OUTDIR = CONFIG.DEFAULT_OUTDIR
    DEFAULT_NCPU = CONFIG.DEFAULT_NCPU
    EXECUTION_PROFILE_NAMES = getattr(
        CONFIG,
        "EXECUTION_PROFILE_NAMES",
        ("science",),
    )

else:

    PIPELINE_STAGES = [
        "ingest",
        "inspection",
        "flagging",
        "calibration",
        "selfcal",
        "imaging",
        "source_detection",
        "ddcal",
    ]
    DEFAULT_OUTDIR = "Products"
    DEFAULT_NCPU = 32
    EXECUTION_PROFILE_NAMES = ("science",)


RUN_CONFIG = (
    getattr(CONFIG, "PIPELINE_RUN_CONFIG", {})
    if CONFIG is not None
    else {}
)

IMAGING_TASK_CONFIG = (
    getattr(CONFIG, "FINAL_IMAGING_TASK_CONFIG", {})
    if CONFIG is not None
    else {}
)

SOURCE_DETECTION_TASK_CONFIG = (
    getattr(CONFIG, "SOURCE_DETECTION_CONFIG", {})
    if CONFIG is not None
    else {}
)

DD_TASK_CONFIG = (
    getattr(CONFIG, "DD_TASK_CONFIG", {})
    if CONFIG is not None
    else {}
)


def get_execution_profile(name):
    """Return the selected profile settings without changing config.py."""

    normalized = str(name or "science").strip().lower()

    if CONFIG is not None and hasattr(CONFIG, "get_execution_profile"):
        return CONFIG.get_execution_profile(normalized)

    if normalized == "science":
        return {}

    raise ValueError(
        "Execution profiles require config.py; use the 'science' profile "
        "when running pipeline.py without it."
    )


STAGE_TO_SCRIPT = {
    "ingest": "fits_ingest.py",
    "inspection": "inspect.py",
    "flagging": "flagging.py",
    "calibration": "calibration.py",
    "selfcal": "selfcal.py",
    "imaging": "imaging.py",
    "source_detection": "source_detection.py",
    "ddcal": "ddcal.py",
}

# ============================================================
# HELPERS
# ============================================================

def command_to_string(command):

    return " ".join(
        shlex.quote(str(part))
        for part in command
    )


def read_json(path, default=None):

    path = Path(path)

    if not path.exists():

        return default

    with open(path) as f:

        return json.load(f)


def utc_now():

    return datetime.now(timezone.utc).isoformat()


def atomic_write_json(path, data):
    """Replace a JSON state file only after the new document is complete."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(
        f".{path.name}.tmp.{os.getpid()}"
    )

    with open(temporary, "w") as handle:

        json.dump(data, handle, indent=2, default=str)
        handle.flush()
        os.fsync(handle.fileno())

    os.replace(temporary, path)


def parse_stage_list(value):

    if not value:

        return None

    if isinstance(value, (list, tuple)):

        return [
            str(item).strip()
            for item in value
            if str(item).strip()
        ]

    return [
        item.strip()
        for item in value.split(",")
        if item.strip()
    ]


class TrackingArgumentParser(argparse.ArgumentParser):
    """Record which option destinations were explicitly supplied."""

    def parse_known_args(self, args=None, namespace=None):

        tokens = list(sys.argv[1:] if args is None else args)
        explicit_destinations = set()
        for token in tokens:
            option = str(token).split("=", 1)[0]
            action = self._option_string_actions.get(option)
            if action is not None:
                explicit_destinations.add(action.dest)

        parsed, extras = super().parse_known_args(tokens, namespace)
        explicit_destinations.update(
            getattr(parsed, "_explicit_cli_dests", ())
        )
        parsed._explicit_cli_dests = frozenset(explicit_destinations)
        return parsed, extras


def csv_value(value):

    if value is None:

        return None

    if isinstance(value, (list, tuple)):

        values = [
            str(item).strip()
            for item in value
            if str(item).strip()
        ]

        return ",".join(values) or None

    value = str(value).strip()

    return value or None


def mapping_value(value):

    if value is None:

        return None

    if isinstance(value, dict):

        pairs = []
        for target, calibrator in value.items():
            target_text = str(target).strip()
            calibrator_text = str(calibrator).strip()
            if target_text and calibrator_text:
                pairs.append(f"{target_text}:{calibrator_text}")

        return ",".join(pairs) or None

    if isinstance(value, (list, tuple)):

        pairs = []
        for item in value:
            if isinstance(item, dict):
                for target, calibrator in item.items():
                    target_text = str(target).strip()
                    calibrator_text = str(calibrator).strip()
                    if target_text and calibrator_text:
                        pairs.append(
                            f"{target_text}:{calibrator_text}"
                        )
            else:
                item_text = str(item).strip()
                if item_text:
                    pairs.append(item_text)

        return ",".join(pairs) or None

    value = str(value).strip()

    return value or None


def sha256_file_or_directory(path):
    """Return the durable content identity of an external CASA model.

    The calibration stage uses the same tree-hash scheme in its provenance.
    Keeping it in the orchestration fingerprint prevents a same-path model
    replacement from reusing old calibration products before CASA is invoked.
    CASA table locks are intentionally excluded because they are transient
    reader state, not model content.
    """

    path = Path(path).expanduser().resolve()
    if not path.exists():
        return None

    digest = hashlib.sha256()

    def update_file(filename):
        with open(filename, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)

    if path.is_file():
        digest.update(b"external-model-file\0")
        update_file(path)
        return digest.hexdigest()

    if not path.is_dir():
        return None

    digest.update(b"external-model-directory-v1\0")
    for entry in sorted(path.rglob("*"), key=lambda item: item.as_posix()):
        if entry.name == ".DS_Store" or entry.name.endswith(".lock"):
            continue
        relative = entry.relative_to(path).as_posix().encode(
            "utf-8", errors="surrogateescape"
        )
        if entry.is_symlink():
            digest.update(b"L\0" + relative + b"\0")
            digest.update(os.readlink(entry).encode(
                "utf-8", errors="surrogateescape"
            ))
        elif entry.is_dir():
            digest.update(b"D\0" + relative + b"\0")
        elif entry.is_file():
            digest.update(b"F\0" + relative + b"\0")
            update_file(entry)

    return digest.hexdigest()


def validate_stage(stage):

    if stage not in PIPELINE_STAGES:

        raise ValueError(
            f"Unknown stage '{stage}'. Valid stages: "
            f"{', '.join(PIPELINE_STAGES)}"
        )


def select_stage_window(stages,
                        from_stage=None,
                        to_stage=None):

    selected = list(stages)

    if from_stage:

        validate_stage(from_stage)
        selected = selected[
            selected.index(from_stage):
        ]

    if to_stage:

        validate_stage(to_stage)
        selected = selected[
            : selected.index(to_stage) + 1
        ]

    return selected


# ============================================================
# WHERE A RUN HAPPENS: THIS NODE OR A SLURM JOB
# ============================================================

# A run holds <outdir>/.pipeline.lock and refreshes its time every
# LOCK_HEARTBEAT_SECONDS. flock alone may not be seen from another node of a
# shared filesystem, so a run started on the Jupyter node and a SLURM job on a
# compute node also read the record: a lock refreshed within
# LOCK_STALE_SECONDS, or naming a SLURM job squeue still lists, is in use.
LOCK_HEARTBEAT_SECONDS = 60
LOCK_STALE_SECONDS = 300
BACKGROUND_LOG_NAME = "pipeline_run.log"
# Lower priority for a run on a shared interactive node, so that the Jupyter
# server stays responsive while CASA works (as `nice -n 19` did on ilifu).
BACKGROUND_NICENESS = 19
MODULE_LOAD_FAILED = 97
MODULE_ENVIRONMENT_SKIP = {"_", "SHLVL", "PWD", "OLDPWD"}


def module_environment(module, runner=subprocess.run, timeout=180):
    """The environment a login shell has after ``module load MODULE``.

    On clusters such as ilifu CASA comes from an environment module, which
    only a shell can load. The pipeline asks one for its environment once
    and runs every CASA stage in it, so neither the terminal nor the SLURM
    job has to load the module first.
    """

    script = (
        f"module load {shlex.quote(str(module))} >/dev/null 2>&1 "
        f"|| exit {MODULE_LOAD_FAILED}; env -0"
    )
    try:
        completed = runner(
            ["bash", "-lc", script],
            capture_output=True,
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"bash is not available to load the module {module}."
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"Loading the module {module} took more than {timeout} s."
        ) from exc
    if completed.returncode == MODULE_LOAD_FAILED:
        raise RuntimeError(
            f"'module load {module}' failed here. 'module avail casa' lists "
            "the CASA modules; or load CASA yourself and leave casa_module "
            "empty."
        )
    if completed.returncode != 0:
        detail = (completed.stderr or b"").decode(errors="replace").strip()
        raise RuntimeError(
            f"Loading the module {module} failed (status "
            f"{completed.returncode}): {detail[-500:]}"
        )
    environment = {}
    for entry in (completed.stdout or b"").split(b"\0"):
        key, separator, value = entry.partition(b"=")
        if key and separator:
            environment[key.decode(errors="replace")] = value.decode(
                errors="replace"
            )
    if not environment.get("PATH"):
        raise RuntimeError(
            f"Loading the module {module} gave no environment back."
        )
    return environment


def apply_module_environment(captured, environment=None):
    """Make this process's environment the module's; return what changed.

    A login shell may rebuild PATH from the system defaults, so entries of
    the current PATH that the module's PATH lacks are kept after its own.
    """

    environment = os.environ if environment is None else environment
    changed = []
    for key, value in captured.items():
        if key in MODULE_ENVIRONMENT_SKIP or key.startswith("BASH_FUNC_"):
            continue
        if key == "PATH":
            entries = [entry for entry in value.split(":") if entry]
            for entry in environment.get("PATH", "").split(":"):
                if entry and entry not in entries:
                    entries.append(entry)
            value = ":".join(entries)
        if environment.get(key) != value:
            environment[key] = value
            changed.append(key)
    return changed


CONTAINER_RUN = re.compile(
    r"\b(singularity|apptainer)\s+exec\s+(?:-\S+\s+)*(\S+\.(?:sif|simg|img))\b"
)


def module_casa_container(module, runner=subprocess.run, timeout=180):
    """The CASA container a module runs, as a command for its Python; or None.

    Some sites (ilifu among them) ship CASA as a modular container whose
    module defines no `casa` program, only shell functions such as
    `casapy () { singularity exec IMAGE.sif xvfb-run -a ipython "$@"; }`.
    A shell function cannot be started from this program, but the image it
    names can: its Python runs the CASA stages, as `casa -c` would.
    """

    script = (
        f"module load {shlex.quote(str(module))} >/dev/null 2>&1 "
        f"|| exit {MODULE_LOAD_FAILED}; type casapy casa 2>/dev/null; "
        "alias 2>/dev/null; true"
    )
    try:
        completed = runner(
            ["bash", "-lc", script],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if completed.returncode != 0:
        return None
    match = CONTAINER_RUN.search(completed.stdout or "")
    if not match:
        return None
    return f"{match.group(1)} exec {match.group(2)} python3"


def casa_python_version(command, runner=subprocess.run, timeout=600):
    """CASA's version from the Python that runs the stages; error if none."""

    probe = (
        "import casatasks, casatools; "
        "print('CASA', casatools.version_string())"
    )
    try:
        completed = runner(
            [*command, "-c", probe],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"Cannot start {command_to_string(command)}: {exc}"
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"{command_to_string(command)} did not start CASA within "
            f"{timeout} s."
        ) from exc
    lines = [
        line for line in (completed.stdout or "").splitlines()
        if line.startswith("CASA ")
    ]
    if completed.returncode != 0 or not lines:
        detail = "\n".join(
            ((completed.stderr or "") + (completed.stdout or "")).splitlines()[-12:]
        )
        raise RuntimeError(
            f"{command_to_string(command)} cannot import casatasks and "
            f"casatools (status {completed.returncode}):\n{detail}"
        )
    return lines[-1][len("CASA "):].strip()


def slurm_job_active(job_id, runner=subprocess.run):
    """True if squeue lists the job, False if not, None if it cannot tell."""

    try:
        completed = runner(
            ["squeue", "-h", "-j", str(job_id), "-o", "%T"],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if completed.returncode != 0:
        # squeue reports an invalid job id once the job has left the queue.
        return False if "Invalid job id" in (completed.stderr or "") else None
    return bool(completed.stdout.strip())


def process_alive(pid):
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except (OSError, TypeError, ValueError):
        return False
    return True


def live_lock_owner(text, lock_path, host=None, now=None, job_active=None):
    """Who holds a run lock, from its record; None if no one does.

    On the same machine, whether the process lives decides. A run on another
    machine is seen only through its record: it is gone once the SLURM job
    it ran in has left the queue, and otherwise it holds the lock while it
    keeps refreshing the lock's time. A listed SLURM job alone proves
    nothing: on ilifu a Jupyter session is a SLURM job too, and it outlives
    a run started in it.
    """

    try:
        record = json.loads(text) if str(text or "").strip() else None
    except ValueError:
        return None
    if not isinstance(record, dict):
        return None
    owner_host = record.get("host")
    if not owner_host:
        # A record from before hosts were written: flock decides.
        return None
    host = host or socket.gethostname()
    pid = record.get("pid")
    if owner_host == host:
        try:
            other = int(pid)
        except (TypeError, ValueError):
            return None
        if other != os.getpid() and process_alive(other):
            return f"process {other} on {host}"
        return None
    job = record.get("slurm_job_id")
    if job and (job_active or slurm_job_active)(job) is False:
        return None
    try:
        age = (now if now is not None else time.time()) - (
            Path(lock_path).stat().st_mtime
        )
    except OSError:
        return None
    if age < LOCK_STALE_SECONDS:
        return (
            f"process {pid} on {owner_host}"
            + (f" (SLURM job {job})" if job else "")
            + f", which refreshed its lock {age:.0f} s ago"
        )
    return None


def current_run_owner(outdir, job_active=None):
    """Who is running the pipeline in ``outdir`` now, or None."""

    lock_path = Path(outdir) / ".pipeline.lock"
    if not lock_path.is_file():
        return None
    try:
        with open(lock_path, "a+") as handle:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                handle.seek(0)
                record = handle.read().strip()
                return f"a process on this machine ({record or 'no record'})"
            handle.seek(0)
            text = handle.read()
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    except OSError:
        return None
    return live_lock_owner(text, lock_path, job_active=job_active)


def memory_mb(value):
    """A SLURM memory size ('120GB', '120G', '237568', '500M') in MB."""

    text = str(value or "").strip().upper().rstrip("+")
    if text.endswith("B"):
        text = text[:-1]
    match = re.fullmatch(r"(\d+(?:\.\d+)?)([KMGT]?)", text)
    if not match:
        return None
    factor = {"K": 1 / 1024, "M": 1, "G": 1024, "T": 1024 * 1024}
    return float(match.group(1)) * factor[match.group(2) or "M"]


def slurm_request_problem(slurm, runner=subprocess.run):
    """Why no node of the partition can take this job, or None.

    SLURM's own refusal ("Requested node configuration is not available")
    does not say which request is too large; this does, before submitting.
    """

    try:
        tasks = int(slurm.get("ntasks") or 1)
        cpus = int(slurm.get("cpus_per_task") or 1)
        nodes = max(1, int(slurm.get("nodes") or 1))
    except (TypeError, ValueError):
        return None
    need_cpus = math.ceil(tasks * cpus / nodes)
    need_memory = memory_mb(slurm.get("mem"))
    partition = slurm.get("partition")
    command = ["sinfo", "-h", "-o", "%c %m"]
    if partition:
        command.extend(["-p", str(partition)])
    try:
        completed = runner(command, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if completed.returncode != 0:
        return None
    sizes = []
    for line in completed.stdout.splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            node_cpus = int(parts[0].rstrip("+"))
        except ValueError:
            continue
        node_memory = memory_mb(parts[1])
        if node_memory is not None:
            sizes.append((node_cpus, node_memory))
    if not sizes:
        return None
    if any(
        node_cpus >= need_cpus and (need_memory is None or node_memory >= need_memory)
        for node_cpus, node_memory in sizes
    ):
        return None
    most_cpus = max(node_cpus for node_cpus, _ in sizes)
    most_memory = max(node_memory for _, node_memory in sizes)
    return (
        f"The job asks for {tasks} task(s) x {cpus} CPUs = {need_cpus} CPUs"
        f" and {slurm.get('mem')} per node, but the nodes of partition "
        f"{partition or '(default)'} have at most {most_cpus} CPUs and "
        f"{most_memory / 1024:.0f} GB. Lower slurm_ntasks, "
        "slurm_cpus_per_task or slurm_mem in config.yaml. CASA stages other "
        "than final imaging run serially: slurm_ntasks: 1 with 16 CPUs is "
        "enough for them."
    )


def strip_launch_options(argv):
    """The command line without --background, --submit or --write-sbatch."""

    kept = []
    skip_value = False
    for index, part in enumerate(argv):
        if skip_value:
            skip_value = False
            continue
        if part in ("--background", "--submit"):
            continue
        if part == "--write-sbatch":
            following = argv[index + 1] if index + 1 < len(argv) else ""
            skip_value = bool(following) and not following.startswith("-")
            continue
        if part.startswith("--write-sbatch="):
            continue
        kept.append(part)
    return kept


def start_in_background(args, argv, popen=subprocess.Popen, wait_seconds=20,
                        job_active=None):
    """Start this run on this node, detached from the terminal.

    The run continues when the terminal or the browser tab is closed. It is
    started at low priority, and its output is appended to
    <outdir>/pipeline_run.log. A restart of the Jupyter server itself still
    stops it; a SLURM job (--submit) does not have that weakness.
    """

    if args.dry_run or args.list_stages or args.resume_status:
        raise RuntimeError(
            "--background starts a run; use --dry-run, --list-stages or "
            "--resume-status without it."
        )
    outdir = Path(args.outdir).expanduser().resolve()
    owner = current_run_owner(outdir, job_active=job_active)
    if owner:
        raise RuntimeError(
            f"A run is already using {outdir}: {owner}. Wait for it to "
            "finish, or stop it first."
        )
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / BACKGROUND_LOG_NAME
    command = [
        sys.executable,
        "-u",
        str(Path(__file__).resolve()),
        *strip_launch_options(argv),
    ]
    host = socket.gethostname()
    with open(log_path, "ab") as log:
        log.write(
            (
                f"\n===== {utc_now()} started on {host}: "
                f"{command_to_string(command)}\n"
            ).encode()
        )
        log.flush()
        environment = dict(os.environ, PYTHONUNBUFFERED="1")
        process = popen(
            command,
            stdin=subprocess.DEVNULL,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            env=environment,
            preexec_fn=lambda: os.nice(BACKGROUND_NICENESS),
        )

    # Most mistakes (no CASA, a missing input, another run) stop a run in
    # its first seconds: wait for them here rather than leave them in a log.
    deadline = time.monotonic() + wait_seconds
    while time.monotonic() < deadline:
        status = process.poll()
        if status is not None:
            print(
                f"The run stopped at once (status {status}). The end of "
                f"{log_path}:"
            )
            try:
                lines = log_path.read_text(errors="replace").splitlines()
            except OSError:
                lines = []
            print("\n".join(lines[-20:]))
            return status or 1
        time.sleep(1)

    print(f"Started on {host} as process {process.pid}, at low priority.")
    session_job = os.environ.get("SLURM_JOB_ID")
    if session_job:
        print(
            f"This session is itself SLURM job {session_job} (on ilifu the "
            "Jupyter server is one); the run lasts as long as that job."
        )
    print(f"Log      : {log_path}")
    print(f"Follow   : tail -f {shlex.quote(str(log_path))}")
    print(f"Stages   : {shlex.quote(str(outdir / 'pipeline_state.json'))}")
    print(f"Stop     : kill {process.pid}")
    print(
        "It keeps running if this terminal or the browser tab is closed, "
        "but not if the Jupyter server restarts; starting the same command "
        "again then continues where it stopped."
    )
    return 0


def submit_to_slurm(args, runner=subprocess.run, which=shutil.which,
                    job_active=None):
    """Write this run's SLURM job, check it fits a node, and submit it."""

    sbatch = which("sbatch")
    if sbatch is None:
        raise RuntimeError(
            "SLURM's sbatch is not available on this machine. Run the "
            "pipeline here with --background instead."
        )
    if CONFIG is None:
        raise RuntimeError("--submit needs config.py for the job settings.")
    outdir = Path(args.outdir).expanduser().resolve()
    owner = current_run_owner(outdir, job_active=job_active)
    if owner:
        raise RuntimeError(
            f"A run is already using {outdir}: {owner}. Wait for it to "
            "finish, or stop it first."
        )
    slurm = CONFIG.ILIFU_SLURM_CONFIG
    problem = slurm_request_problem(slurm, runner=runner)
    if problem:
        raise RuntimeError(problem)
    if int(slurm.get("ntasks") or 1) > 1 and args.mpi == "off":
        print(
            f"Note: the job asks for {slurm.get('ntasks')} tasks, but with "
            "--mpi off only one is used; slurm_ntasks: 1 would start sooner."
        )

    path = write_sbatch_from_config(args)
    completed = runner(
        [sbatch, "--parsable", str(path)],
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        print((completed.stderr or completed.stdout or "").strip())
        print(
            "sbatch refused the job. 'sinfo -s' lists the partitions; "
            "change the slurm_ settings in config.yaml and submit again."
        )
        return completed.returncode or 1
    job = completed.stdout.strip().split(";")[0]
    log_dir = Path(slurm.get("log_dir", "slurm_logs"))
    name = str(slurm.get("job_name") or "pipeline")
    output = str(log_dir / slurm.get("output", "%x-%j.out"))
    error = str(log_dir / slurm.get("error", "%x-%j.err"))
    output = output.replace("%j", job).replace("%x", name)
    error = error.replace("%j", job).replace("%x", name)
    print(f"Submitted SLURM job {job} ({path}).")
    print("Queue    : squeue -u $USER")
    print(f"Log      : tail -f {output}   (errors: {error})")
    print(f"Stages   : {outdir / 'pipeline_state.json'}")
    print(f"Stop     : scancel {job}")
    return 0


def write_sbatch_from_config(args):

    if CONFIG is None or not hasattr(CONFIG, "write_sbatch"):

        raise RuntimeError(
            "config.py does not provide write_sbatch()."
        )

    path = (
        None
        if args.write_sbatch == "__CONFIG_DEFAULT__"
        else args.write_sbatch
    )

    embedded = []

    if args.ms:
        embedded.append(args.ms)

    embedded.extend(["--outdir", args.outdir])

    for option, value in [
        ("--fits-format", args.fits_format),
        ("--imported-ms", args.imported_ms),
        ("--telescope", args.telescope),
        ("--receiver-band", args.receiver_band),
        ("--data-state", args.data_state),
        ("--imager", args.imager),
        ("--wsclean", args.wsclean),
        ("--wsclean-data-column", args.wsclean_data_column),
        ("--lofar-beam-mode", args.lofar_beam_mode),
        ("--target", args.target),
        ("--targets", ",".join(args.targets) if args.targets else None),
        ("--target-selection", args.target_selection),
        (
            "--auto-target-max-flagged-fraction",
            args.auto_target_max_flagged_fraction,
        ),
        ("--flux-calibrator", args.flux_calibrator),
        ("--bandpass-calibrator", args.bandpass_calibrator),
        ("--flux-standard", args.flux_standard),
        ("--flux-model", args.flux_model),
        ("--flux-model-type", args.flux_model_type),
        (
            "--secondary-flux-calibrator",
            args.secondary_flux_calibrator,
        ),
        ("--phase-calibrator", args.phase_calibrator),
        ("--target-phase-map", mapping_value(args.target_phase_map)),
        ("--target-fields", csv_value(args.target_fields)),
        ("--reference-antenna", args.reference_antenna),
        ("--stages", csv_value(args.stages)),
        ("--from-stage", args.from_stage),
        ("--to-stage", args.to_stage),
        ("--working-ms", args.working_ms),
        ("--selfcal-ms", args.selfcal_ms),
        ("--imaging-ms", args.imaging_ms),
        ("--dd-ms", args.dd_ms),
        ("--casa", args.casa),
        ("--casa-module", getattr(args, "casa_module", None)),
        ("--casa-python", getattr(args, "casa_python", None)),
        ("--python", args.python),
        ("--ncpu", args.ncpu),
        ("--mpicasa", args.mpicasa),
        ("--progress-interval", args.progress_interval),
        ("--solver-command-template", args.solver_command_template),
        ("--apply-command-template", args.apply_command_template),
    ]:
        if value not in (None, ""):
            embedded.extend([option, str(value)])

    for flagfile in getattr(args, "gmrt_flagfile", []) or []:
        embedded.extend(["--gmrt-flagfile", str(flagfile)])

    if args.ms and (
        args.fits_format == "lta" or LTA_NAME_PATTERN.search(Path(args.ms).name)
    ):
        embedded.extend(["--listscan", args.listscan, "--gvfits", args.gvfits])
        if args.keep_lta_fits:
            embedded.append("--keep-lta-fits")

    embedded.append(
        "--copy-input-ms"
        if args.copy_input_ms
        else "--no-copy-input-ms"
    )
    embedded.append(
        "--include-ddcal"
        if args.include_ddcal
        else "--no-ddcal"
    )
    embedded.append("--resume" if args.resume else "--no-resume")
    if args.continue_on_error:
        embedded.append("--continue-on-error")
    if args.subbands is not None:
        embedded.extend(["--subbands", str(args.subbands)])
    if args.no_subbands is True:
        embedded.append("--no-subbands")
    elif args.no_subbands is False:
        embedded.append("--make-subbands")
    embedded.extend(
        [
            "--dd-mode",
            args.dd_mode,
            "--dd-engine",
            args.dd_engine,
            (
                "--prepare-ddcal"
                if args.dd_prepare_only
                else "--execute-ddcal"
            ),
            (
                "--skip-dd-post-imaging"
                if args.skip_dd_post_imaging
                else "--run-dd-post-imaging"
            ),
            "--profile",
            args.profile,
            "--mpi",
            args.mpi,
        ]
    )

    if args.threads is not None:
        embedded.extend(["--threads", str(args.threads)])
    if args.mpi_processes is not None:
        embedded.extend(
            ["--mpi-processes", str(args.mpi_processes)]
        )

    embedded.extend(args.sbatch_extra_arg)

    outfile = CONFIG.write_sbatch(
        path=path,
        pipeline_script=Path(SCRIPT_DIR) / "pipeline.py",
        dry_run=args.sbatch_dry_run,
        extra_args=embedded,
    )

    print(f"Slurm sbatch file written: {outfile}")
    return outfile


# ============================================================
# PIPELINE RUNNER
# ============================================================

class PipelineRunner:

    def __init__(self, args):

        self.args = args
        self._input_explicit = bool(getattr(self.args, "ms", None))
        # A repeated invocation of the exact same command must be a resume
        # when this output directory already has a run manifest.  This is
        # intentionally independent of the YAML default: operators often
        # launch CASA with only ``-c pipeline.py`` and expect a killed job to
        # continue without having to add ``--resume`` on the retry.  An
        # explicit ``--no-resume`` remains authoritative because the tracking
        # parser records both resume flags under the ``resume`` destination.
        explicit_resume = "resume" in set(
            getattr(self.args, "_explicit_cli_dests", ())
        )
        configured_outdir = Path(
            getattr(self.args, "outdir", DEFAULT_OUTDIR)
        ).expanduser().resolve()
        if (
            not explicit_resume and
            (configured_outdir / "pipeline_run_manifest.json").is_file()
        ):
            self.args.resume = True
        configured_input = RUN_CONFIG.get("ms") or None
        self._input_from_config = bool(
            not getattr(self.args, "ms", None) and configured_input
        )
        if self._input_from_config:
            self.args.ms = configured_input
        self.args.target_selection = str(
            getattr(self.args, "target_selection", "auto") or "auto"
        ).strip().lower()
        if self.args.target_selection not in {"all", "auto", "explicit"}:
            raise ValueError(
                "--target-selection must be 'auto' or 'explicit'."
            )
        self.args.auto_target_max_flagged_fraction = float(
            getattr(
                self.args,
                "auto_target_max_flagged_fraction",
                0.80,
            )
        )
        if not (
            math.isfinite(self.args.auto_target_max_flagged_fraction) and
            0.0 <= self.args.auto_target_max_flagged_fraction <= 1.0
        ):
            raise ValueError(
                "--auto-target-max-flagged-fraction must be between "
                "0 and 1."
            )
        self._target_selection_origin = (
            "explicit" if getattr(self.args, "target", None) else None
        )
        requested_targets = getattr(self.args, "targets", None)
        if isinstance(requested_targets, str):
            requested_targets = [
                item.strip() for item in requested_targets.split(",")
            ]
        self.args.targets = [
            str(item).strip()
            for item in (requested_targets or [])
            if str(item).strip()
        ]
        self.profile_settings = get_execution_profile(args.profile)
        self.root = Path(SCRIPT_DIR)
        self.outdir = Path(args.outdir).expanduser().resolve()
        self.args.outdir = str(self.outdir)
        # Per-target stages write below this root, which is the output
        # directory itself until a target scope is entered.
        self.stage_outdir = self.outdir
        self.target_scope_name = None

        self.input_data = None
        self.input_kind = None
        self.imported_ms = None
        self.args.fits_format = self.normalize_fits_format(
            getattr(self.args, "fits_format", "auto")
        )
        self.args.gmrt_flagfile = [
            str(Path(item).expanduser().resolve())
            for item in (getattr(self.args, "gmrt_flagfile", []) or [])
        ]

        if self.args.ms:
            self.input_data = Path(
                self.args.ms
            ).expanduser().resolve()
            self.input_kind = self.classify_input_path(
                self.input_data,
                self.args.fits_format,
            )

            if self.input_kind == "fits":
                imported_override = getattr(
                    self.args,
                    "imported_ms",
                    None,
                )
                self.imported_ms = (
                    self.resolve_output_path(imported_override)
                    if imported_override
                    else self.default_imported_ms_path(self.input_data)
                )
                self.args.imported_ms = str(self.imported_ms)
                self.args.ms = str(self.imported_ms)
            else:
                self.args.ms = str(self.input_data)
                if getattr(self.args, "imported_ms", None):
                    raise ValueError(
                        "--imported-ms is valid only when the positional "
                        "input is a visibility FITS file."
                    )

        self.commands = []
        self.source_ms = None
        self.working_ms = None
        self.resume_state = None
        self.run_lock_handle = None
        self._available_cpu_count = None
        self._thread_limit = None
        self._configuration_digest = None
        self._stage_fingerprints = {}
        self._stage_recipe_fingerprints = {}
        # Recorded recipe fingerprints proven equal to the computed ones; see
        # recorded_recipe_matches_after_provenance_refreshes.
        self._stage_recipe_pins = {}
        self._mpi_launcher = None
        self._mpi_launcher_checked = False
        self._casa_module_loaded = False
        self._lock_heartbeat = None
        self._recorded_source_ms_identity = None
        self._prepared_working_from_manifest = False
        self._calibration_split_repair_required = False

    # --------------------------------------------------------
    # Input classification
    # --------------------------------------------------------

    @staticmethod
    def normalize_fits_format(value):

        normalized = str(value or "auto").strip().lower().replace("_", "-")
        aliases = {
            "auto": "auto",
            "uvfits": "uvfits",
            "fitsidi": "fitsidi",
            "fits-idi": "fitsidi",
            "lta": "lta",
        }
        if normalized not in aliases:
            raise ValueError(
                "--fits-format must be auto, uvfits, fitsidi, or lta."
            )
        return aliases[normalized]

    # --------------------------------------------------------

    @staticmethod
    def is_lta_path(path):
        """True for a GMRT LTA file name (.lta, .ltb, .lta_1 ...)."""

        return bool(LTA_NAME_PATTERN.search(Path(path).name))

    # --------------------------------------------------------

    def input_is_lta(self):
        """The positional input is a GMRT LTA file to convert and import."""

        return (
            self.input_kind == "fits" and
            self.input_data is not None and
            (
                self.args.fits_format == "lta" or (
                    self.args.fits_format == "auto" and
                    self.is_lta_path(self.input_data)
                )
            )
        )

    # --------------------------------------------------------

    @staticmethod
    def classify_input_path(path, fits_format="auto"):
        """Classify a positional input without opening CASA tables."""

        path = Path(path)
        if fits_format != "auto":
            return "fits"
        if LTA_NAME_PATTERN.search(path.name):
            # GMRT LTA goes through the ingest stage like visibility FITS.
            return "fits"
        if path.name.lower().endswith(".ms"):
            return "ms"
        if path.is_dir():
            return "ms"
        if path.is_file():
            return "fits"
        if path.name.lower().endswith(
            (
                ".fits",
                ".fit",
                ".fts",
                ".uvfits",
                ".fitsidi",
                ".fits-idi",
                ".idi",
            )
        ):
            return "fits"
        return "ms"

    # --------------------------------------------------------

    @staticmethod
    def resolve_output_path(path):
        """Resolve an output parent without following the final component."""

        raw = Path(path).expanduser()
        if not raw.is_absolute():
            raw = Path.cwd() / raw
        return raw.parent.resolve() / raw.name

    # --------------------------------------------------------

    def default_imported_ms_path(self, input_path):

        name = Path(input_path).name
        match = LTA_NAME_PATTERN.search(name)
        if match:
            # X.lta -> X.ms; X.ltb -> X_ltb.ms; X.lta_1 -> X_lta_1.ms, the
            # names lta_convert.py gives the converted FITS.
            suffix = match.group(0).lower()
            name = name[:match.start()] + (
                "" if suffix == ".lta" else "_" + suffix[1:]
            )
        lowered = name.lower()
        for suffix in [
            ".fits-idi",
            ".fitsidi",
            ".uvfits",
            ".fits",
            ".fit",
            ".fts",
            ".idi",
        ]:
            if lowered.endswith(suffix):
                name = name[:-len(suffix)]
                break
        clean = "".join(
            character if character.isalnum() or character in "-_" else "_"
            for character in name
        ).strip("_") or "visibility_input"
        return (self.outdir / "ingest" / f"{clean}.ms").resolve()

    # --------------------------------------------------------
    # Telescope/backend context
    # --------------------------------------------------------

    def resolved_telescope(self):

        if self.args.telescope not in (None, "", "auto"):
            return TELESCOPE_PROFILE_MODULE.normalize_telescope(
                self.args.telescope
            )

        inspection = read_json(
            self.report_path("inspection"),
            {},
        )
        return TELESCOPE_PROFILE_MODULE.telescope_name_from_report(
            inspection
        )

    # --------------------------------------------------------

    def resolved_receiver_band(self):

        if self.args.receiver_band:
            telescope = self.resolved_telescope()
            if telescope:
                return TELESCOPE_PROFILE_MODULE.normalize_band(
                    telescope,
                    self.args.receiver_band,
                )
            return str(self.args.receiver_band)

        inspection = read_json(
            self.report_path("inspection"),
            {},
        )
        band = inspection.get("band")
        if isinstance(band, dict):
            band = band.get("name") or band.get("band")
        return str(band).upper() if band else None

    # --------------------------------------------------------

    def wsclean_data_column(self, telescope):
        """The data column WSClean reads.

        LOFAR images a calibrated MS it was handed, and keeps the configured
        column (CORRECTED_DATA unless set otherwise). MeerKAT is imaged from
        this pipeline's own self-calibrated or split MS, whose calibrated
        visibilities are in DATA, so its column is chosen the way the CASA
        backend chooses it. Setting wsclean_data_column to anything but its
        default selects that column explicitly for either telescope.
        """

        configured = self.args.wsclean_data_column
        if telescope == "meerkat" and str(configured).upper() == (
            "CORRECTED_DATA"
        ):
            return "auto"
        return configured

    def imaging_backend(self):

        if self.args.imager != "auto":
            return self.args.imager

        telescope = self.resolved_telescope()
        if telescope:
            return TELESCOPE_PROFILE_MODULE.default_imager(telescope)

        # Before a first inspection the original CASA backend is the only
        # backward-compatible choice.  A detected LOFAR or MeerKAT report
        # switches this to WSClean before the imaging stage is built.
        return "casa"

    # --------------------------------------------------------
    # Process resources
    # --------------------------------------------------------

    @staticmethod
    def _positive_environment_integer(name):

        value = os.environ.get(name)

        if value is None:
            return None

        # SLURM_CPUS_ON_NODE can contain values such as ``64(x2)``.
        value = str(value).strip().split("(", 1)[0]

        try:
            value = int(value)
        except (TypeError, ValueError):
            return None

        return value if value > 0 else None

    # --------------------------------------------------------

    def available_cpu_count(self):
        """Return CPUs usable by this process, respecting Slurm/cgroups."""

        if self._available_cpu_count is not None:
            return self._available_cpu_count

        limits = []

        try:
            limits.append(len(os.sched_getaffinity(0)))
        except (AttributeError, OSError):
            pass

        if os.environ.get("SLURM_JOB_ID"):
            slurm_totals = []

            for variable in [
                "SLURM_CPUS_ON_NODE",
                "SLURM_JOB_CPUS_PER_NODE",
            ]:
                value = self._positive_environment_integer(variable)

                if value is not None:
                    slurm_totals.append(value)

            tasks = self._positive_environment_integer("SLURM_NTASKS")
            cpus_per_task = self._positive_environment_integer(
                "SLURM_CPUS_PER_TASK"
            )

            if tasks and cpus_per_task:
                slurm_totals.append(tasks * cpus_per_task)
            elif cpus_per_task:
                slurm_totals.append(cpus_per_task)

            if slurm_totals:
                limits.append(max(slurm_totals))

        detected = os.cpu_count()

        if detected:
            limits.append(detected)

        self._available_cpu_count = max(1, min(limits or [1]))
        return self._available_cpu_count

    # --------------------------------------------------------

    def thread_limit(self):
        """Return the numerical thread count for serial CASA/Python work."""

        if self._thread_limit is not None:
            return self._thread_limit

        requested = self.args.threads

        if requested is None:
            requested = self.profile_settings.get("threads")

        if requested is None:
            requested = self.available_cpu_count()

        try:
            requested = int(requested)
        except (TypeError, ValueError) as exc:
            raise ValueError("--threads must be a positive integer.") from exc

        if requested < 1:
            raise ValueError("--threads must be a positive integer.")

        self._thread_limit = min(
            requested,
            self.available_cpu_count(),
        )
        return self._thread_limit

    # --------------------------------------------------------

    def stage_thread_limit(self, stage):
        """Apply memory-safe execution caps without changing science."""

        limit = self.thread_limit()
        if stage != "flagging":
            return limit

        requested = getattr(self.args, "flagging_threads", None)
        if requested in (None, ""):
            return limit
        try:
            requested = int(requested)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "--flagging-threads must be a positive integer."
            ) from exc
        if requested < 1:
            raise ValueError(
                "--flagging-threads must be a positive integer."
            )
        return min(limit, requested)

    # --------------------------------------------------------

    def mpi_process_count(self):
        """Return the requested CASA MPI process count, including rank 0."""

        requested = self.args.mpi_processes

        if requested is None:
            # Continuum tclean parallelizes its major cycle. Using every CPU
            # for MPI can waste memory and time in gather operations; 16 is a
            # portable high-throughput default and remains user-overridable.
            requested = min(16, self.available_cpu_count())

        try:
            requested = int(requested)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "--mpi-processes must be an integer greater than one."
            ) from exc

        if requested < 2:
            raise ValueError(
                "--mpi-processes must be an integer greater than one."
            )

        count = min(
            requested,
            self.mpi_rank_capacity(),
            self.available_cpu_count(),
        )

        if count < 2:
            raise RuntimeError(
                "CASA MPI needs at least two allocated Slurm tasks/ranks. "
                "Request --ntasks=2 or more, or use --mpi off."
            )

        return count

    # --------------------------------------------------------

    def mpi_rank_capacity(self):
        """Return MPI ranks exposed by Slurm, or local CPUs when forced."""

        if os.environ.get("SLURM_JOB_ID"):
            return (
                self._positive_environment_integer("SLURM_NTASKS") or
                1
            )

        return self.available_cpu_count()

    # --------------------------------------------------------

    def mpi_threads_per_process(self):
        """Use Slurm's hybrid MPI/OpenMP layout without oversubscription."""

        if not os.environ.get("SLURM_JOB_ID"):
            return 1

        requested = (
            self._positive_environment_integer("SLURM_CPUS_PER_TASK") or
            1
        )
        maximum = max(
            1,
            self.available_cpu_count() // self.mpi_process_count(),
        )
        return min(requested, maximum)

    # --------------------------------------------------------

    @staticmethod
    def executable_path(command):

        if not command:
            return None

        expanded = Path(str(command)).expanduser()

        if expanded.is_file() and os.access(expanded, os.X_OK):
            return str(expanded.resolve())

        return shutil.which(str(command))

    # --------------------------------------------------------

    def load_casa_module(self):
        """Load the configured CASA module into this run's environment.

        Every stage the run starts inherits it, so CASA is found whether the
        run was started from a terminal or by a SLURM job. A dry run only
        warns when the module cannot be loaded.
        """

        module = getattr(self.args, "casa_module", None)
        if not module or self._casa_module_loaded:
            return
        try:
            captured = module_environment(module)
        except RuntimeError as exc:
            if self.args.dry_run:
                print(f"WARNING: {exc}")
                return
            raise
        apply_module_environment(captured)
        self._casa_module_loaded = True
        # mpicasa is looked for beside the CASA the module provides.
        self._mpi_launcher_checked = False
        if self.casa_python_command():
            print(
                f"Loaded module {module}; CASA stages run with "
                f"{self.args.casa_python}"
            )
            return
        casa = shutil.which(self.args.casa)
        if casa:
            print(f"Loaded module {module}: casa is {casa}")
            return
        container = module_casa_container(module)
        if container:
            self.args.casa_python = container
            print(
                f"Loaded module {module}: it has no `casa` program but runs "
                f"CASA from a container, so the CASA stages run with "
                f"{container}"
            )
            return
        print(f"Loaded module {module}: casa is still not found")

    # --------------------------------------------------------

    def mpi_launcher(self):
        """Find the mpicasa paired with CASA where possible."""

        if self._mpi_launcher_checked:
            return self._mpi_launcher

        self._mpi_launcher_checked = True
        requested = self.args.mpicasa
        casa_path = self.executable_path(self.args.casa)

        if (
            requested in (None, "", "mpicasa") and
            casa_path is not None
        ):
            sibling = Path(casa_path).with_name("mpicasa")

            if sibling.is_file() and os.access(sibling, os.X_OK):
                self._mpi_launcher = str(sibling)
                return self._mpi_launcher

        self._mpi_launcher = self.executable_path(requested or "mpicasa")
        return self._mpi_launcher

    # --------------------------------------------------------

    def mpi_enabled(self):
        """Resolve auto/on/off CASA MPI mode without affecting macOS."""

        mode = self.args.mpi

        if mode == "off":
            return False

        if self.casa_python_command():
            # A modular CASA has no mpicasa; --mpi on is refused in preflight.
            return False

        launcher = self.mpi_launcher()

        if mode == "on":
            if platform.system() == "Darwin":
                raise RuntimeError(
                    "CASA mpicasa is not supported on macOS. Use --mpi off "
                    "(or the default auto mode) for local Mac runs."
                )
            return launcher is not None

        # CASA's mpicasa launcher is not supported on macOS. On Linux, auto
        # mode is intentionally conservative outside a scheduler allocation.
        return bool(
            platform.system() != "Darwin" and
            os.environ.get("SLURM_JOB_ID") and
            launcher is not None and
            self.mpi_rank_capacity() >= 2 and
            self.available_cpu_count() >= 2
        )

    # --------------------------------------------------------

    def stage_uses_mpi(self, stage):

        return bool(
            stage == "imaging" and
            self.imaging_backend() == "casa" and
            self.mpi_enabled()
        )

    # --------------------------------------------------------

    def casa_memory_budget_mb(self, snapshot=None):
        """The memory figure this run hands CASA, in megabytes.

        The casarc writes it and the pre-flight check plans from it, so both
        describe the same run.
        """

        fraction = RUN_CONFIG.get("casa_memory_fraction", 0.5)
        try:
            fraction = float(fraction)
        except (TypeError, ValueError):
            fraction = 0.0
        if not 0.0 < fraction <= 1.0:
            return None

        if snapshot is None:
            snapshot = self.stage_memory_snapshot()
        limit = MEMORY_DIAGNOSTICS.effective_limit_bytes(snapshot)
        if not limit:
            return None

        return max(512, int(limit * fraction / (1024 ** 2)))

    def casa_resource_rc(self, stage, threads, snapshot=None):
        """Tell CASA how much memory this job may actually use.

        CASA sizes some of its work from the memory it believes the machine
        has. casacore reads a job's cgroup limit only from a v1 hierarchy and
        otherwise reports the node's free memory, so under Slurm it can plan
        against far more than the job is allowed. w-projection is the sharp
        case: it computes its convolution functions in chunks of up to 100
        planes sized from that figure, which for an 8192-pixel image is one
        allocation of about 19 GiB against a node's 176 GiB of free memory,
        and about 3 GiB when it knows the job has 28 GiB.

        casacore reads these keys from the files named by CASARCFILES, in
        order, first match winning, so the run's file goes first and the
        user's own files keep working behind it.
        """

        megabytes = self.casa_memory_budget_mb(snapshot)
        if not megabytes:
            return None

        path = self.outdir / "casa" / f"{stage}_resources.rc"
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(
                "# Written by pipeline.py for this run.\n"
                f"system.resources.memory: {megabytes}\n"
                f"system.resources.cores: {max(1, int(threads))}\n"
            )
        except OSError as exc:
            print(f"WARNING: could not write {path}: {exc}", flush=True)
            return None

        limit = MEMORY_DIAGNOSTICS.effective_limit_bytes(
            snapshot if snapshot is not None else self.stage_memory_snapshot()
        )
        print(
            f"CASA memory budget: {megabytes / 1024:.1f} GiB of the "
            f"{MEMORY_DIAGNOSTICS.format_bytes(limit)} this job may use.",
            flush=True,
        )
        return path

    @staticmethod
    def casa_resource_file_list(rc_path, environ):
        """CASARCFILES with the run's file first, then casacore's defaults."""

        existing = environ.get("CASARCFILES")
        if existing:
            return f"{rc_path}:{existing}"

        files = [str(rc_path)]
        home = environ.get("HOME")
        if home:
            files.extend([
                f"{home}/.casarc",
                f"{home}/.casa/rc",
                f"{home}/.aipsrc",
            ])

        # CASAPATH is "root arch site host"; casacore reads .aipsrc/aipsrc
        # from each of those in this order.
        parts = (environ.get("CASAPATH") or "").split()
        if parts:
            root = parts[0]
            arch = f"{root}/{parts[1]}" if len(parts) > 1 else None
            site = f"{arch}/{parts[2]}" if arch and len(parts) > 2 else None
            host = f"{site}/{parts[3]}" if site and len(parts) > 3 else None
            files.append(f"{root}/.aipsrc")
            for directory in (host, site, arch):
                if directory:
                    files.append(f"{directory}/aipsrc")

        return ":".join(files)

    def stage_environment(self, stage, memory_snapshot=None):
        """Configure numerical runtimes without CPU oversubscription."""

        environment = os.environ.copy()
        limit = (
            1
            if stage == "ddcal"
            else
            self.mpi_threads_per_process()
            if self.stage_uses_mpi(stage)
            else self.stage_thread_limit(stage)
        )

        value = str(limit)

        for variable in [
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
        ]:
            environment[variable] = value

        environment["OMP_DYNAMIC"] = "FALSE"
        environment["PYTHONUNBUFFERED"] = "1"

        if stage in CASA_STAGES:
            rc_path = self.casa_resource_rc(stage, limit, memory_snapshot)
            if rc_path is not None:
                environment["CASARCFILES"] = self.casa_resource_file_list(
                    rc_path,
                    environment,
                )

        scratch = os.environ.get("SLURM_TMPDIR")

        if scratch and Path(scratch).is_dir() and os.access(scratch, os.W_OK):
            environment["TMPDIR"] = scratch

        return environment

    # --------------------------------------------------------
    # Product/report helpers
    # --------------------------------------------------------

    def report_path(self, name):

        root = self.outdir if name in SHARED_REPORTS else self.stage_outdir

        if CONFIG is not None:

            return root / CONFIG.REPORT_FILES[name]

        fallback = {
            "ingest": "ingest/ingest_report.json",
            "inspection": "inspection/inspection_report.json",
            "flagging_pre": "flagging/preflag_summary.json",
            "flagging_post": "flagging/postflag_summary.json",
            "calibration": "calibration/calibration_report.json",
            "split_targets": "calibration/split_targets.json",
            "selfcal": "selfcal/reports/selfcal_report.json",
            "imaging": "imaging/reports/imaging_report.json",
            "source_detection": (
                "source_detection/reports/source_detection_report.json"
            ),
            "ddcal": "ddcal/reports/ddcal_report.json",
        }

        return root / fallback[name]

    # --------------------------------------------------------

    def stage_script(self, stage):

        if stage == "imaging" and self.imaging_backend() == "wsclean":
            return self.root / "wsclean_imaging.py"

        return self.root / STAGE_TO_SCRIPT[stage]

    # --------------------------------------------------------

    def casa_python_command(self):
        """The command that starts CASA's Python, as a list; None for `casa`.

        A modular CASA (a Python with casatasks and casatools, often in a
        container) runs each stage as `<python> script args`; the full CASA
        program runs it as `casa --nologger --nogui -c script args`.
        """

        value = getattr(self.args, "casa_python", None)
        if not value:
            return None
        return shlex.split(str(value))

    # --------------------------------------------------------

    def casa_command(self, script, *stage_args, mpi=False):

        python = self.casa_python_command()
        if python:
            if mpi:
                raise RuntimeError(
                    "CASA MPI needs mpicasa, which a modular CASA (casa_python) "
                    "does not have. Use --mpi off."
                )
            return [
                *python,
                str(script),
                *[str(arg) for arg in stage_args],
            ]

        casa = self.args.casa
        command = []

        if mpi:
            launcher = self.mpi_launcher()

            if launcher is None:
                raise RuntimeError(
                    "CASA MPI was selected but mpicasa was not found. "
                    "Load the matching CASA installation, pass --mpicasa, "
                    "or use --mpi off."
                )

            casa = self.executable_path(casa) or casa
            command.extend([
                launcher,
                "-n",
                str(self.mpi_process_count()),
            ])

        command.extend([
            casa,
            "--nologger",
            "--nogui",
            "-c",
            str(script),
            *[
                str(arg)
                for arg in stage_args
            ],
        ])

        return command

    # --------------------------------------------------------

    def python_command(self, script, *stage_args):

        return [
            self.args.python,
            str(script),
            *[
                str(arg)
                for arg in stage_args
            ],
        ]

    # --------------------------------------------------------
    # Input resolution
    # --------------------------------------------------------

    @staticmethod
    def raw_ms_stages():

        return {
            "inspection",
            "flagging",
            "calibration",
        }

    # --------------------------------------------------------

    def stage_requires_primary_ms(self, stage):
        """Return whether *stage* consumes the positional input MS.

        Later stages can be driven by an explicit target MS or by a durable
        upstream report.  A FITS import must not be injected merely because
        such a stage happens to operate on some MS.
        """

        if stage in self.raw_ms_stages():
            return True

        if stage == "selfcal":
            return bool(
                self.args.data_state == "calibrated"
                and not self.args.selfcal_ms
            )

        if stage == "imaging":
            if self.args.imaging_ms:
                return False

            selfcal_report = read_json(
                self.report_path("selfcal"),
                {},
            )
            final_ms = selfcal_report.get("final_ms")
            if (
                final_ms
                and Path(final_ms).exists()
                and (
                    not self.args.target
                    or selfcal_report.get("target") == self.args.target
                )
            ):
                return False

            return self.args.data_state == "calibrated"

        return False

    def stages_require_primary_ms(self, stages):

        return any(
            self.stage_requires_primary_ms(stage)
            for stage in stages
        )

    # --------------------------------------------------------

    def discard_unused_configured_primary_input(self, stages):
        """Do not turn a config default into an explicit downstream input."""

        if (
            not self._input_from_config or
            self.stages_require_primary_ms(stages)
        ):
            return

        configured = self.input_data
        self.args.ms = None
        self.input_data = None
        self.input_kind = None
        self.imported_ms = None
        self.args.imported_ms = None
        self.source_ms = None
        self.working_ms = None
        self._stage_fingerprints.clear()
        self._stage_recipe_fingerprints.clear()
        print(
            "Ignoring configured primary input for this report-driven "
            "downstream-only run: "
            f"{configured}"
        )

    # --------------------------------------------------------

    def resume_option_was_explicit(self, destination):

        return destination in set(
            getattr(self.args, "_explicit_cli_dests", ())
        )

    # --------------------------------------------------------

    def restore_manifest_science_context(self, manifest):
        """Restore omitted science choices without hiding CLI overrides."""

        request = manifest.get("instrument_request", {})
        if isinstance(request, dict):
            for key in [
                "telescope",
                "receiver_band",
                "data_state",
                "imager",
                "lofar_beam_mode",
            ]:
                if (
                    key in request and
                    not self.resume_option_was_explicit(key)
                ):
                    setattr(self.args, key, request.get(key))

        profile = manifest.get("execution_profile")
        if profile and not self.resume_option_was_explicit("profile"):
            self.args.profile = profile
            self.profile_settings = get_execution_profile(profile)

        fields = manifest.get("field_selection", {})
        if isinstance(fields, dict):
            field_mapping = {
                "flux_calibrator": "flux_calibrator",
                "bandpass_calibrator": "bandpass_calibrator",
                "secondary_flux_calibrator": "secondary_flux_calibrator",
                "phase_calibrator": "phase_calibrator",
                "target_phase_map": "target_phase_map",
                "target_fields": "target_fields",
                "reference_antenna": "reference_antenna",
            }
            for manifest_key, argument_name in field_mapping.items():
                if (
                    manifest_key in fields and
                    not self.resume_option_was_explicit(argument_name)
                ):
                    setattr(
                        self.args,
                        argument_name,
                        fields.get(manifest_key),
                    )

        downstream = manifest.get("downstream_selection", {})
        if isinstance(downstream, dict):
            if (
                "target" in downstream and
                not self.resume_option_was_explicit("target")
            ):
                self.args.target = downstream.get("target")
            if (
                "target_selection" in downstream and
                not self.resume_option_was_explicit("target_selection")
            ):
                self.args.target_selection = downstream.get(
                    "target_selection"
                )
            if (
                "auto_target_max_flagged_fraction" in downstream and
                not self.resume_option_was_explicit(
                    "auto_target_max_flagged_fraction"
                )
            ):
                self.args.auto_target_max_flagged_fraction = float(
                    downstream.get("auto_target_max_flagged_fraction")
                )
            if downstream.get("mode") in {"explicit", "automatic"}:
                self._target_selection_origin = downstream.get("mode")

        calibrator_model = manifest.get("calibrator_model", {})
        if isinstance(calibrator_model, dict):
            for key in [
                "flux_standard",
                "flux_model",
                "flux_model_type",
            ]:
                if (
                    key in calibrator_model and
                    not self.resume_option_was_explicit(key)
                ):
                    setattr(self.args, key, calibrator_model.get(key))

        if not self.resume_option_was_explicit("gmrt_flagfile"):
            self.args.gmrt_flagfile = list(
                manifest.get("gmrt_flagfiles", []) or []
            )

    # --------------------------------------------------------

    def rehydrate_resume_run_context(self, stages):
        """Restore the effective MS for an input-less stage-scoped resume."""

        if (
            not self.args.resume or
            self._input_explicit
        ):
            return

        manifest_path = self.run_manifest_path()
        if not manifest_path.is_file():
            return
        try:
            manifest = read_json(manifest_path, {})
        except (OSError, ValueError, TypeError) as exc:
            raise RuntimeError(
                "Cannot restore the resume context because the run manifest "
                f"is unreadable: {exc}"
            ) from exc
        if not isinstance(manifest, dict):
            raise RuntimeError(
                "Cannot restore the resume context from a non-object run "
                "manifest."
            )

        self.restore_manifest_science_context(manifest)
        effective_stages = self.selected_stages()
        if not self.stages_require_primary_ms(effective_stages):
            self._stage_fingerprints.clear()
            self._stage_recipe_fingerprints.clear()
            return

        input_kind = str(manifest.get("input_kind") or "ms").lower()
        source_value = manifest.get("input_data") or manifest.get("source_ms")
        if not source_value or input_kind not in {"ms", "fits"}:
            raise RuntimeError(
                "The run manifest cannot reconstruct the original input "
                "kind/path needed by the requested resume stage."
            )
        source = Path(source_value).expanduser().resolve()

        if input_kind == "ms":
            working_value = manifest.get("working_ms")
            working = (
                Path(working_value).expanduser().resolve()
                if working_value else None
            )
            effective = (
                working
                if working is not None and self.casa_table_exists(working)
                else source
            )
            if not self.casa_table_exists(effective):
                raise RuntimeError(
                    "Neither the manifest-managed working MS nor the "
                    "original source MS is available for this resume: "
                    f"{effective}"
                )
            self.input_data = source
            self.input_kind = "ms"
            self.source_ms = Path(
                manifest.get("source_ms") or source
            ).expanduser().resolve()
            self.working_ms = working
            self.args.ms = str(effective)
            self._recorded_source_ms_identity = manifest.get(
                "source_ms_identity"
            )
            self._prepared_working_from_manifest = bool(
                working is not None and effective == working
            )
        else:
            imported_value = manifest.get("imported_ms")
            if not imported_value:
                raise RuntimeError(
                    "The FITS run manifest does not record its imported MS."
                )
            imported = Path(imported_value).expanduser().resolve()
            if "ingest" not in stages and not self.casa_table_exists(imported):
                raise RuntimeError(
                    "The manifest-imported Measurement Set is unavailable: "
                    f"{imported}"
                )
            self.input_data = source
            self.input_kind = "fits"
            self.imported_ms = imported
            self.args.imported_ms = str(imported)
            self.args.ms = str(imported)
            self._prepared_working_from_manifest = (
                "ingest" not in effective_stages
            )

        self._stage_fingerprints.clear()
        self._stage_recipe_fingerprints.clear()
        print(
            "Restored resume context from the run manifest; effective "
            f"Measurement Set: {self.args.ms}"
        )

    # --------------------------------------------------------

    @staticmethod
    def mutable_ms_stages():

        return {
            "flagging",
            "calibration",
        }

    # --------------------------------------------------------

    @staticmethod
    def is_within(path, parent):

        try:
            path.relative_to(parent)
        except ValueError:
            return False

        return True

    # --------------------------------------------------------

    def working_ms_path(self, source_ms):

        if self.args.working_ms:
            return Path(self.args.working_ms).expanduser().resolve()

        return (
            self.outdir.resolve() /
            "working_ms" /
            source_ms.name
        )

    # --------------------------------------------------------

    def working_ms_manifest_path(self, working_ms):

        return (
            working_ms.parent /
            f"{working_ms.name}.provenance.json"
        )

    # --------------------------------------------------------

    @staticmethod
    def measurement_set_identity(ms_path):
        """Return a bounded identity for an immutable source MS.

        Visibility columns can be terabytes, so hashing every byte is not a
        practical restart check.  This identity binds the directory inode and
        every table member's relative path, inode, size and modification time,
        while content-hashing the small CASA table metadata files.  Replacing
        or editing an MS at the same pathname therefore invalidates a managed
        working copy without reading the full visibility payload.
        """

        path = Path(ms_path).expanduser().resolve()
        if not (path.is_dir() and (path / "table.dat").is_file()):
            raise ValueError(f"Not a CASA Measurement Set table: {path}")

        digest = hashlib.sha256()
        root_stat = path.stat()
        entry_count = 0

        for directory, dirnames, filenames in os.walk(path):
            # CASA creates transient table-lock entries while a table is
            # opened.  They are not science data and must not make a stable
            # source identity change merely because another reader ran.
            dirnames[:] = sorted(
                name for name in dirnames
                if not name.lower().endswith(".lock")
            )
            filenames = sorted(
                name for name in filenames
                if not name.lower().endswith(".lock")
            )
            directory_path = Path(directory)

            for name in dirnames + filenames:
                member = directory_path / name
                member_stat = member.stat(follow_symlinks=False)
                relative = member.relative_to(path).as_posix()
                record = (
                    relative,
                    int(member_stat.st_mode),
                    int(member_stat.st_dev),
                    int(member_stat.st_ino),
                    int(member_stat.st_size),
                    int(member_stat.st_mtime_ns),
                )
                digest.update(repr(record).encode("utf-8"))
                entry_count += 1

                if (
                    not member.is_symlink() and
                    member.is_file() and
                    name in {"table.dat", "table.info"} and
                    member_stat.st_size <= 4 * 1024 * 1024
                ):
                    with member.open("rb") as handle:
                        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                            digest.update(chunk)

        return {
            "schema_version": 2,
            "root_device": int(root_stat.st_dev),
            "root_inode": int(root_stat.st_ino),
            "entry_count": entry_count,
            "stat_sha256": digest.hexdigest(),
        }

    # --------------------------------------------------------

    @staticmethod
    def measurement_set_instance_identity(ms_path):
        return IMAGING_SAFETY.measurement_set_instance_identity(ms_path)

    # --------------------------------------------------------

    def source_ms_identity_for_provenance(self):
        """Return the live or manifest-anchored immutable source identity."""

        if self._recorded_source_ms_identity is not None:
            return self._recorded_source_ms_identity
        if self.input_kind != "ms" or self.input_data is None:
            return None
        try:
            return self.measurement_set_identity(self.input_data)
        except (OSError, ValueError):
            return None

    # --------------------------------------------------------

    @staticmethod
    def _remove_copy_staging(path):
        """Remove only a pipeline-owned, uniquely named copy staging path."""

        path = Path(path)

        if f".copying.{os.getpid()}" not in path.name:
            raise RuntimeError(
                f"Refusing to clean an unexpected copy path: {path}"
            )

        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()

    # --------------------------------------------------------

    def copy_measurement_set(self, source_ms, working_ms):
        """Copy an MS with native/reflink tools and an atomic destination."""

        staging = working_ms.with_name(
            f".{working_ms.name}.copying.{os.getpid()}"
        )

        if staging.exists():
            raise RuntimeError(
                f"Measurement Set copy staging path already exists: {staging}"
            )

        methods = []

        if sys.platform == "darwin" and Path("/bin/cp").is_file():
            # APFS clonefile copy-on-write; if unsupported, ordinary cp is
            # attempted next. Both preserve the source MS.
            methods.extend([
                (
                    "macOS cp clonefile",
                    ["/bin/cp", "-cRp", str(source_ms), str(staging)],
                ),
                (
                    "macOS cp",
                    ["/bin/cp", "-Rp", str(source_ms), str(staging)],
                ),
            ])

        elif shutil.which("cp"):
            methods.extend([
                (
                    "GNU cp reflink=auto",
                    [
                        shutil.which("cp"),
                        "-a",
                        "--reflink=auto",
                        "--",
                        str(source_ms),
                        str(staging),
                    ],
                ),
                (
                    "native cp",
                    [
                        shutil.which("cp"),
                        "-a",
                        "--",
                        str(source_ms),
                        str(staging),
                    ],
                ),
            ])

        errors = []

        for method, command in methods:

            result = subprocess.run(
                command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )

            if result.returncode == 0 and self.casa_table_exists(staging):
                os.replace(staging, working_ms)
                return method

            if result.returncode == 0:
                errors.append(
                    f"{method}: copied directory is not a CASA table"
                )
            else:
                errors.append(
                    f"{method}: "
                    f"{result.stderr.strip() or result.returncode}"
                )

            if staging.exists():
                self._remove_copy_staging(staging)

        try:
            shutil.copytree(
                source_ms,
                staging,
                copy_function=shutil.copy2,
            )

            if not self.casa_table_exists(staging):
                raise RuntimeError(
                    "copied directory is not a CASA table"
                )

            os.replace(staging, working_ms)
            return "shutil.copytree"

        except Exception as exc:

            if staging.exists():
                self._remove_copy_staging(staging)

            detail = "; ".join(errors)
            raise RuntimeError(
                "Could not create the working Measurement Set copy. "
                f"Native-copy attempts: {detail or 'none'}"
            ) from exc

    # --------------------------------------------------------

    def adopt_completed_calibration_working_copy(
        self,
        source_ms,
        working_ms,
        source_identity,
        existing_provenance,
    ):
        """Repair legacy working-copy provenance using independent evidence.

        This migration is intentionally limited to an MS whose completed
        calibration report fingerprints its *live* post-calibration contents.
        Path agreement alone is not enough to adopt a mutable CASA table.
        """

        if not self.args.resume or not self.casa_table_exists(working_ms):
            return False, "resume is disabled or the working MS is invalid"

        if existing_provenance and (
            existing_provenance.get("source_ms") not in (None, str(source_ms))
            or existing_provenance.get("working_ms") not in (
                None,
                str(working_ms),
            )
        ):
            return False, "existing provenance names a different source or MS"

        run_manifest_path = self.run_manifest_path()
        try:
            run_manifest = read_json(run_manifest_path, {})
        except (OSError, ValueError, TypeError) as exc:
            return False, f"run manifest is unreadable: {exc}"
        if not isinstance(run_manifest, dict):
            return False, "run manifest is not a JSON object"

        manifest_source = run_manifest.get("source_ms")
        manifest_working = run_manifest.get("working_ms")
        try:
            manifest_source = str(
                Path(manifest_source).expanduser().resolve()
            )
            manifest_working = str(
                Path(manifest_working).expanduser().resolve()
            )
        except (TypeError, ValueError):
            return False, "run manifest has no valid source/working-MS paths"
        if (
            manifest_source != str(source_ms) or
            manifest_working != str(working_ms) or
            run_manifest.get("source_ms_identity") != source_identity or
            str(run_manifest.get("input_kind") or "ms").lower() != "ms"
        ):
            return False, (
                "run manifest does not bind the live immutable source to "
                "this working MS"
            )

        calibration_marker = (
            self.outdir / "calibration" / "calibration.done"
        )
        calibration_report_path = self.report_path("calibration")
        try:
            calibration_report = read_json(calibration_report_path, {})
        except (OSError, ValueError, TypeError) as exc:
            return False, f"calibration report is unreadable: {exc}"
        if (
            not calibration_marker.is_file() or
            not isinstance(calibration_report, dict)
        ):
            return False, "completed calibration evidence is unavailable"

        reported_input = (
            calibration_report.get("measurement_set_used") or
            calibration_report.get("original_measurement_set")
        )
        try:
            reported_input = str(
                Path(reported_input).expanduser().resolve()
            )
        except (TypeError, ValueError):
            return False, "calibration report has no valid working-MS path"
        if reported_input != str(working_ms):
            return False, "calibration report names a different working MS"

        try:
            live_working_identity = IMAGING_SAFETY.path_stat_fingerprint(
                working_ms
            )
        except (OSError, TypeError, ValueError) as exc:
            return False, f"working MS cannot be fingerprinted: {exc}"
        if (
            calibration_report.get("measurement_set_stat_sha256") !=
            live_working_identity
        ):
            return False, (
                "live working-MS contents do not match the completed "
                "calibration report"
            )

        try:
            inspection_report = read_json(
                self.report_path("inspection"),
                {},
            )
        except (OSError, ValueError, TypeError) as exc:
            return False, f"inspection report is unreadable: {exc}"
        inspected_ms = (
            inspection_report.get("measurement_set")
            if isinstance(inspection_report, dict)
            else None
        )
        try:
            inspected_ms = str(Path(inspected_ms).expanduser().resolve())
        except (TypeError, ValueError):
            return False, "inspection report has no valid working-MS path"
        if inspected_ms != str(working_ms):
            return False, "inspection report names a different working MS"

        provenance_path = self.working_ms_manifest_path(working_ms)
        previous_sha256 = None
        previous_archive = None
        if provenance_path.is_file():
            try:
                previous_sha256 = hashlib.sha256(
                    provenance_path.read_bytes()
                ).hexdigest()
            except OSError:
                previous_sha256 = None
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            archive_path = provenance_path.with_name(
                f"{provenance_path.name}.legacy.{stamp}"
            )
            counter = 1
            while archive_path.exists():
                archive_path = provenance_path.with_name(
                    f"{provenance_path.name}.legacy.{stamp}.{counter}"
                )
                counter += 1
            shutil.copy2(provenance_path, archive_path)
            previous_archive = str(archive_path)
        repaired = {
            "schema_version": 2,
            "source_ms": str(source_ms),
            "working_ms": str(working_ms),
            "source_identity": source_identity,
            "created_utc": (
                existing_provenance.get("created_utc")
                if isinstance(existing_provenance, dict)
                else None
            ) or run_manifest.get("created_utc") or utc_now(),
            "copy_method": (
                existing_provenance.get("copy_method")
                if isinstance(existing_provenance, dict)
                else None
            ) or "legacy pipeline-managed copy",
            "adopted_utc": utc_now(),
            "adoption": {
                "method": "completed-calibration-evidence-v1",
                "run_manifest": str(run_manifest_path),
                "calibration_report": str(calibration_report_path),
                "calibration_measurement_set_stat_sha256": (
                    live_working_identity
                ),
                "previous_provenance_sha256": previous_sha256,
                "previous_provenance_archive": previous_archive,
            },
        }
        atomic_write_json(provenance_path, repaired)
        print(
            "Repaired legacy working-MS provenance using the immutable "
            "source manifest and matching completed-calibration fingerprint:\n"
            f"{provenance_path}"
        )
        return True, "completed calibration evidence matches"

    # --------------------------------------------------------

    def prepare_working_copy(self, stages):
        """Copy the raw MS and use that copy for mutable CASA stages."""

        if self.input_kind == "fits":
            if self.args.working_ms:
                raise RuntimeError(
                    "--working-ms applies to an MS input. For visibility "
                    "FITS, choose the conversion destination with "
                    "--imported-ms instead."
                )
            if (
                self.args.copy_input_ms and
                self.mutable_ms_stages().intersection(stages)
            ):
                print(
                    "Visibility FITS or LTA input is read-only; the newly imported "
                    "MS is already the pipeline working copy, so no second "
                    "MS copy is needed."
                )
            return

        if not self.args.copy_input_ms:
            return

        if not self.mutable_ms_stages().intersection(stages):
            return

        if self._prepared_working_from_manifest:
            source_ms = self.source_ms
            working_ms = self.working_ms
            if (
                source_ms is None or
                working_ms is None or
                not self.casa_table_exists(working_ms)
            ):
                raise RuntimeError(
                    "The run manifest does not identify a usable managed "
                    "working Measurement Set for this resume."
                )
            manifest_path = self.working_ms_manifest_path(working_ms)
            try:
                provenance = read_json(manifest_path, {})
            except (OSError, ValueError, TypeError):
                provenance = {}
            source_identity = self._recorded_source_ms_identity
            if not isinstance(source_identity, dict):
                raise RuntimeError(
                    "The run manifest has no immutable source-MS identity; "
                    "this managed working copy cannot be resumed safely."
                )
            if (
                isinstance(provenance, dict) and
                provenance.get("source_ms") == str(source_ms) and
                provenance.get("working_ms") == str(working_ms) and
                provenance.get("source_identity") == source_identity
            ):
                print(
                    "Reusing manifest-managed working Measurement Set:\n"
                    f"{working_ms}"
                )
                return
            adopted, reason = self.adopt_completed_calibration_working_copy(
                source_ms,
                working_ms,
                source_identity,
                provenance if isinstance(provenance, dict) else {},
            )
            if adopted:
                print(
                    "Reusing validated manifest-managed working "
                    f"Measurement Set:\n{working_ms}"
                )
                return
            raise RuntimeError(
                "The manifest-managed working Measurement Set cannot be "
                f"safely resumed because {reason}."
            )

        source_ms = Path(self.args.ms).expanduser().resolve()
        working_ms = self.working_ms_path(source_ms)

        if (
            source_ms == working_ms or
            self.is_within(working_ms, source_ms) or
            self.is_within(source_ms, working_ms)
        ):
            raise RuntimeError(
                "The working MS must be separate from the input MS:\n"
                f"input:   {source_ms}\n"
                f"working: {working_ms}"
            )

        self.source_ms = source_ms
        self.working_ms = working_ms
        self.args.ms = str(working_ms)

        if self.args.dry_run:
            return

        manifest_path = self.working_ms_manifest_path(working_ms)
        source_identity = self.measurement_set_identity(source_ms)

        if working_ms.exists():
            try:
                manifest = read_json(manifest_path, {})
            except (OSError, ValueError, TypeError):
                manifest = {}

            if (
                manifest.get("source_ms") == str(source_ms) and
                manifest.get("working_ms") == str(working_ms) and
                manifest.get("source_identity") == source_identity and
                self.casa_table_exists(working_ms)
            ):
                print(
                    "Reusing existing working Measurement Set:\n"
                    f"{working_ms}"
                )
                return

            adopted, adoption_reason = (
                self.adopt_completed_calibration_working_copy(
                    source_ms,
                    working_ms,
                    source_identity,
                    manifest if isinstance(manifest, dict) else {},
                )
            )
            if adopted:
                print(
                    "Reusing validated existing working Measurement Set:\n"
                    f"{working_ms}"
                )
                return

            raise RuntimeError(
                "Working Measurement Set already exists without matching "
                "input provenance:\n"
                f"{working_ms}\n"
                f"Safe adoption was rejected because {adoption_reason}.\n"
                "Choose --working-ms with a new path or remove the "
                "incomplete/stale working copy after checking it."
            )

        working_ms.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        print(
            "Creating working Measurement Set copy:\n"
            f"  source:  {source_ms}\n"
            f"  working: {working_ms}"
        )

        copy_method = self.copy_measurement_set(
            source_ms,
            working_ms,
        )

        after_copy_identity = self.measurement_set_identity(source_ms)
        if after_copy_identity != source_identity:
            raise RuntimeError(
                "The source Measurement Set changed while its managed "
                "working copy was being created. The unprovenanced copy was "
                f"preserved for diagnosis at {working_ms}; do not use it."
            )

        if not self.casa_table_exists(working_ms):
            raise RuntimeError(
                "The copied Measurement Set is not a readable CASA table: "
                f"{working_ms}"
            )

        print(f"Measurement Set copy method: {copy_method}")

        manifest = {
            "schema_version": 2,
            "source_ms": str(source_ms),
            "working_ms": str(working_ms),
            "source_identity": source_identity,
            "created_utc": datetime.now(
                timezone.utc,
            ).isoformat(),
            "copy_method": copy_method,
        }

        atomic_write_json(manifest_path, manifest)

    # --------------------------------------------------------

    def target_products_root(self, target):
        """Where one target's per-target stages keep their products.

        A run made before per-target directories existed keeps its products
        where they are, so that its completed self-calibration and imaging are
        still validated and reused. Every other target gets its own directory.
        """

        if not target:
            return self.outdir

        for name in ("selfcal", "imaging"):
            report = read_json(
                self.outdir / (
                    CONFIG.REPORT_FILES[name]
                    if CONFIG is not None else
                    {
                        "selfcal": "selfcal/reports/selfcal_report.json",
                        "imaging": "imaging/reports/imaging_report.json",
                    }[name]
                ),
                {},
            )
            if (
                isinstance(report, dict) and
                str(report.get("target") or "").strip() == str(target).strip()
            ):
                return self.outdir

        return self.outdir / "targets" / target_directory_name(target)

    def science_target_scopes(self, stages):
        """Return the (target, products root) pairs the per-target stages run for."""

        if not set(PER_TARGET_STAGES).intersection(stages):
            return []

        # An explicitly supplied downstream Measurement Set, or an already
        # calibrated input, names one dataset rather than a target list.
        if (
            self.args.data_state == "calibrated" or
            self.args.selfcal_ms or
            self.args.imaging_ms
        ):
            return [(self.args.target, self.outdir)]

        available = self.available_science_targets()
        if not available:
            return [(self.args.target, self.outdir)]

        requested = [
            item for item in (getattr(self.args, "targets", None) or [])
        ]
        if self.args.target and self.args.target not in requested:
            requested = [self.args.target] + requested

        if requested:
            unknown = [item for item in requested if item not in available]
            if unknown:
                raise RuntimeError(
                    "Requested target(s) "
                    + ", ".join(repr(item) for item in unknown)
                    + " are not among the inspected/calibrated science "
                    "targets. Available targets: " + ", ".join(available)
                )
            chosen = [item for item in available if item in requested]
        elif self.args.target_selection == "all":
            chosen = list(available)
        elif len(available) == 1:
            chosen = list(available)
        elif self.args.target_selection == "auto":
            chosen = [self.automatically_select_science_target(available)]
        else:
            raise RuntimeError(
                "Multiple science targets are available for downstream "
                "self-calibration/imaging: " + ", ".join(available) + ". "
                "Target selection is configured as explicit. Rerun with "
                "--target TARGET_NAME, --targets A,B, or "
                "--target-selection all."
            )

        return [(target, self.target_products_root(target)) for target in chosen]

    @contextlib.contextmanager
    def target_scope(self, target, root):
        """Run the per-target stages for one target, under its own root."""

        previous = (
            self.args.target,
            self.stage_outdir,
            self.target_scope_name,
        )
        self.args.target = target
        self.stage_outdir = Path(root)
        self.target_scope_name = target
        self._stage_fingerprints.clear()
        self._stage_recipe_fingerprints.clear()
        try:
            yield
        finally:
            (
                self.args.target,
                self.stage_outdir,
                self.target_scope_name,
            ) = previous
            self._stage_fingerprints.clear()
            self._stage_recipe_fingerprints.clear()

    def resolve_calibrated_target_ms(self):

        if self.args.selfcal_ms:

            return self.args.selfcal_ms

        if self.args.data_state == "calibrated" and self.args.ms:

            return self.args.ms

        split_report = read_json(
            self.report_path("split_targets"),
            {},
        )

        if split_report:

            if self.args.target:
                if self.args.target in split_report:
                    return split_report[self.args.target]
                raise RuntimeError(
                    f"Requested target {self.args.target!r} is not present "
                    "in the calibrated target products. Available targets: "
                    + ", ".join(map(str, split_report))
                )

            if len(split_report) != 1:

                raise RuntimeError(
                    "Multiple calibrated target MSs are available. "
                    "Available targets: "
                    + ", ".join(map(str, split_report))
                    + ". "
                    "Pass --target TARGET_NAME (or set target in "
                    "config.yaml), or pass --selfcal-ms explicitly."
                )

            only_target = next(iter(split_report))

            return split_report[only_target]

        if self.args.dry_run:

            return (
                "<target MS from "
                f"{self.report_path('split_targets')}>"
            )

        raise RuntimeError(
            "Cannot resolve self-calibration input MS. "
            "Run calibration first or pass --selfcal-ms."
        )

    # --------------------------------------------------------

    def available_science_targets(self):
        """Return ordered target candidates from durable products/reports."""

        split_report = read_json(
            self.report_path("split_targets"),
            {},
        )
        if isinstance(split_report, dict) and split_report:
            return list(split_report)

        inspection = read_json(
            self.report_path("inspection"),
            {},
        )
        fields = inspection.get("fields", {})
        targets = fields.get("targets", []) if isinstance(fields, dict) else []
        if not isinstance(targets, list):
            return []
        return list(dict.fromkeys(
            str(target).strip()
            for target in targets
            if str(target).strip()
        ))

    # --------------------------------------------------------

    def target_selection_report_path(self):

        return self.outdir / "target_selection.json"

    # --------------------------------------------------------

    def persist_downstream_target_selection(self, report):
        """Commit an automatic target decision to the run provenance."""

        if self.args.dry_run:
            return

        atomic_write_json(self.target_selection_report_path(), report)

        manifest_path = self.run_manifest_path()
        if not manifest_path.is_file():
            return

        manifest = read_json(manifest_path, {})
        if not isinstance(manifest, dict):
            raise RuntimeError(
                f"Run manifest is not a JSON object: {manifest_path}"
            )
        manifest["downstream_selection"] = {
            "target": self.args.target,
            "mode": "automatic",
            "target_selection": self.args.target_selection,
            "auto_target_max_flagged_fraction": (
                self.args.auto_target_max_flagged_fraction
            ),
            "policy": report["policy"],
            "max_flagged_fraction": report["max_flagged_fraction"],
            "report": str(self.target_selection_report_path()),
            "source_statistics_sha256": report[
                "source_statistics_sha256"
            ],
        }
        manifest["updated_utc"] = utc_now()
        atomic_write_json(manifest_path, manifest)

    # --------------------------------------------------------

    def automatically_select_science_target(
        self,
        targets,
        *,
        defer_until_flagging=False,
    ):
        """Choose one target from validated post-flag field statistics.

        The score is deliberately independent of field ordering and source
        names.  A target must have finite, internally consistent CASA summary
        counts and remain below the configured flagged-fraction ceiling.
        """

        statistics_path = self.outdir / "flagging" / "field_statistics.json"
        if not statistics_path.is_file():
            if defer_until_flagging:
                return None
            raise RuntimeError(
                "Cannot safely choose among multiple science targets because "
                "post-flagging field statistics are missing: "
                f"{statistics_path}. Run/validate the flagging stage, or pass "
                "--target TARGET_NAME explicitly."
            )

        try:
            raw_bytes = statistics_path.read_bytes()
            statistics = json.loads(raw_bytes)
        except (OSError, ValueError, TypeError) as exc:
            raise RuntimeError(
                "Cannot safely choose a science target because the field "
                f"statistics are unreadable: {statistics_path}: {exc}"
            ) from exc

        by_field = (
            statistics.get("by_field", {})
            if isinstance(statistics, dict)
            else {}
        )
        if not isinstance(by_field, dict):
            raise RuntimeError(
                "Cannot safely choose a science target because "
                f"{statistics_path} has no object-valued 'by_field' map."
            )

        maximum = self.args.auto_target_max_flagged_fraction
        ranking = []
        invalid = []

        for target in targets:
            summary = by_field.get(target)
            if not isinstance(summary, dict):
                invalid.append(f"{target}: missing summary")
                continue

            total = summary.get("total")
            flagged = summary.get("flagged")
            if (
                isinstance(total, bool) or
                isinstance(flagged, bool)
            ):
                invalid.append(f"{target}: non-numeric counts")
                continue
            try:
                total = float(total)
                flagged = float(flagged)
            except (TypeError, ValueError):
                invalid.append(f"{target}: non-numeric counts")
                continue
            if not (math.isfinite(total) and math.isfinite(flagged)):
                invalid.append(f"{target}: non-finite counts")
                continue
            if total <= 0 or flagged < 0 or flagged > total:
                invalid.append(
                    f"{target}: inconsistent total={total:g}, "
                    f"flagged={flagged:g}"
                )
                continue

            usable = total - flagged
            flagged_fraction = flagged / total
            eligible = usable > 0 and flagged_fraction <= maximum
            reason = (
                "eligible"
                if eligible
                else (
                    "no usable visibilities"
                    if usable <= 0
                    else "flagged fraction exceeds safety ceiling"
                )
            )
            ranking.append({
                "target": target,
                "total_visibilities": total,
                "flagged_visibilities": flagged,
                "usable_visibilities": usable,
                "flagged_fraction": flagged_fraction,
                "eligible": eligible,
                "reason": reason,
            })

        if invalid:
            raise RuntimeError(
                "Cannot safely choose a science target because field-level "
                "flagging metrics are incomplete or invalid: "
                + "; ".join(invalid)
                + ". Pass --target TARGET_NAME explicitly only after "
                "reviewing the data."
            )

        ranking.sort(key=lambda item: (
            not item["eligible"],
            -item["usable_visibilities"],
            item["flagged_fraction"],
            item["target"],
        ))
        eligible = [item for item in ranking if item["eligible"]]
        if not eligible:
            fractions = ", ".join(
                f"{item['target']}={item['flagged_fraction']:.1%}"
                for item in ranking
            )
            raise RuntimeError(
                "No science target passes the automatic-selection safety "
                f"ceiling ({maximum:.1%} flagged). Observed: {fractions}. "
                "Inspect the flagging products and pass --target explicitly "
                "only if the data are scientifically usable."
            )

        selected = eligible[0]["target"]
        self.args.target = selected
        self._target_selection_origin = "automatic"
        self._stage_fingerprints.clear()
        self._stage_recipe_fingerprints.clear()

        report = {
            "schema_version": 1,
            "selected_utc": utc_now(),
            "selected_target": selected,
            "policy": "max_usable_postflag_visibilities",
            "max_flagged_fraction": maximum,
            "source_statistics": str(statistics_path),
            "source_statistics_sha256": hashlib.sha256(
                raw_bytes
            ).hexdigest(),
            "tie_breaking": [
                "largest usable visibility count",
                "lowest flagged fraction",
                "lexicographically smallest target name",
            ],
            "ranking": ranking,
        }
        self.persist_downstream_target_selection(report)

        winner = eligible[0]
        print(
            "Automatically selected science target "
            f"{selected!r}: {winner['usable_visibilities']:g} usable "
            f"visibilities, {winner['flagged_fraction']:.1%} flagged. "
            "The ranked decision is recorded in "
            f"{self.target_selection_report_path()}. Use --target NAME to "
            "override this choice explicitly."
        )
        return selected

    # --------------------------------------------------------

    def validate_downstream_target_selection(
        self,
        stages,
        *,
        defer_until_flagging=False,
    ):
        """Resolve one safe target before target-specific downstream work."""

        downstream = {
            "selfcal", "imaging", "source_detection", "ddcal"
        }.intersection(stages)
        if not downstream or self.args.data_state == "calibrated":
            return

        targets = self.available_science_targets()
        if not targets:
            return

        for requested in (
            [self.args.target]
            if self.args.target and
            self._target_selection_origin != "automatic" else []
        ) + list(getattr(self.args, "targets", None) or []):
            if requested not in targets:
                raise RuntimeError(
                    f"Requested target {requested!r} is not one of "
                    "the inspected/calibrated science targets. Available "
                    "targets: " + ", ".join(targets)
                )
        if (
            self.args.target and
            self._target_selection_origin != "automatic"
        ) or getattr(self.args, "targets", None):
            return

        # Every target is imaged unless one was named or another mode is set.
        if self.args.target_selection == "all":
            return

        # An explicitly supplied downstream MS is already an unambiguous
        # selection. Otherwise the orchestrator must never guess which target
        # receives expensive self-calibration and publication imaging.
        explicit_input = (
            ("selfcal" in downstream and self.args.selfcal_ms) or
            ("selfcal" not in downstream and self.args.imaging_ms)
        )
        if explicit_input or len(targets) == 1:
            return

        if self.args.target_selection == "auto":
            return self.automatically_select_science_target(
                targets,
                defer_until_flagging=defer_until_flagging,
            )

        raise RuntimeError(
            "Multiple science targets are available for downstream "
            "self-calibration/imaging: " + ", ".join(targets) + ". "
            "Target selection is configured as explicit. Rerun the same "
            "output directory with --target TARGET_NAME, or opt into the "
            "validated selector with --target-selection auto. "
            "Completed inspection, flagging, and calibration products will "
            "be validated and reused."
        )

    # --------------------------------------------------------

    def resolve_imaging_ms(self):

        if self.args.imaging_ms:

            return self.args.imaging_ms

        selfcal_report = read_json(
            self.report_path("selfcal"),
            {},
        )

        final_ms = selfcal_report.get("final_ms")
        selfcal_target = selfcal_report.get("target")

        if (
            final_ms and
            Path(final_ms).exists() and
            (
                not self.args.target or
                selfcal_target == self.args.target
            )
        ):

            return final_ms

        if self.args.data_state == "calibrated" and self.args.ms:

            return self.args.ms

        split_report = read_json(
            self.report_path("split_targets"),
            {},
        )
        if split_report:
            if self.args.target:
                if self.args.target in split_report:
                    return split_report[self.args.target]
                raise RuntimeError(
                    f"Requested imaging target {self.args.target!r} is not "
                    "present in split_targets.json. Available targets: "
                    + ", ".join(map(str, split_report))
                )
            if len(split_report) == 1:
                return next(iter(split_report.values()))
            raise RuntimeError(
                "Multiple calibrated target MSs are available for imaging. "
                "Pass --target TARGET_NAME (or set target in config.yaml), "
                "or pass --imaging-ms explicitly."
            )

        if final_ms and Path(final_ms).exists() and self.args.target:
            raise RuntimeError(
                f"The available self-calibrated MS belongs to target "
                f"{selfcal_target!r}, not requested target "
                f"{self.args.target!r}. Use the matching target or pass a "
                "validated --imaging-ms from a separate target run."
            )

        return None

    # --------------------------------------------------------

    def resolve_dd_ms(self):

        if self.args.dd_ms:

            return self.args.dd_ms

        imaging_report = read_json(
            self.report_path("imaging"),
            {},
        )

        ms = imaging_report.get("measurement_set")

        if ms and Path(ms).exists():

            if (
                self.args.target and
                imaging_report.get("target") != self.args.target
            ):
                raise RuntimeError(
                    f"The imaging report belongs to target "
                    f"{imaging_report.get('target')!r}, not requested DD "
                    f"target {self.args.target!r}. Use a target-specific "
                    "output directory."
                )

            return ms

        selfcal_report = read_json(
            self.report_path("selfcal"),
            {},
        )

        final_ms = selfcal_report.get("final_ms")

        if final_ms and Path(final_ms).exists():

            if (
                self.args.target and
                selfcal_report.get("target") != self.args.target
            ):
                raise RuntimeError(
                    f"The self-calibration report belongs to target "
                    f"{selfcal_report.get('target')!r}, not requested DD "
                    f"target {self.args.target!r}."
                )

            return final_ms

        return None

    # --------------------------------------------------------
    # Stage command builders
    # --------------------------------------------------------

    def command_for_stage(self, stage):

        script = self.stage_script(stage)

        if not script.exists():

            raise FileNotFoundError(
                f"Stage script not found: {script}"
            )

        if stage == "ingest":
            if self.input_kind != "fits" or self.input_data is None:
                raise RuntimeError(
                    "The ingest stage requires a UVFITS, FITS-IDI or GMRT "
                    "LTA file as the positional input."
                )
            stage_args = [
                str(self.input_data),
                str(self.imported_ms),
                "--outdir",
                self.args.outdir,
                "--fits-format",
                self.args.fits_format,
                "--telescope",
                self.args.telescope,
            ]
            for flagfile in getattr(self.args, "gmrt_flagfile", []) or []:
                stage_args.extend(["--gmrt-flagfile", str(flagfile)])
            if self.input_is_lta():
                # Where the converters live is not science; the resume
                # comparison of commands leaves these options out.
                stage_args.extend([
                    "--listscan", self.args.listscan,
                    "--gvfits", self.args.gvfits,
                ])
                if self.args.keep_lta_fits:
                    stage_args.append("--keep-lta-fits")
            return self.casa_command(script, *stage_args)

        if stage in [
            "inspection",
            "flagging",
            "calibration",
        ]:

            if not self.args.ms:

                raise RuntimeError(
                    f"Stage '{stage}' requires an input MS."
                )

            stage_args = [
                self.args.ms,
                "--outdir",
                self.args.outdir,
            ]

            if stage == "inspection":
                stage_args.extend(
                    [
                        "--telescope",
                        self.args.telescope,
                        "--data-state",
                        self.args.data_state,
                    ]
                )
                if self.args.receiver_band:
                    stage_args.extend(
                        ["--receiver-band", self.args.receiver_band]
                    )

            if stage == "flagging":
                stage_args.extend(
                    [
                        "--profile",
                        self.args.profile,
                        "--provenance-fingerprint",
                        self.stage_fingerprint("flagging"),
                    ]
                )

            if stage == "calibration":

                if self._calibration_split_repair_required:
                    stage_args.append("--repair-split-targets-only")

                calibration_options = [
                    ("--flux-calibrator", self.args.flux_calibrator),
                    (
                        "--bandpass-calibrator",
                        self.args.bandpass_calibrator,
                    ),
                    (
                        "--secondary-flux-calibrator",
                        self.args.secondary_flux_calibrator,
                    ),
                    ("--phase-calibrator", self.args.phase_calibrator),
                    (
                        "--target-phase-map",
                        mapping_value(self.args.target_phase_map),
                    ),
                    ("--targets", csv_value(self.args.target_fields)),
                    ("--refant", self.args.reference_antenna),
                    ("--telescope", self.args.telescope),
                    ("--flux-standard", self.args.flux_standard),
                    ("--flux-model", self.args.flux_model),
                    ("--flux-model-type", self.args.flux_model_type),
                    (
                        "--provenance-fingerprint",
                        self.stage_fingerprint("calibration"),
                    ),
                    (
                        "--recipe-fingerprint",
                        self.stage_recipe_fingerprint("calibration"),
                    ),
                ]

                for option, value in calibration_options:

                    if value:

                        stage_args.extend([option, value])

            return self.casa_command(script, *stage_args)

        if stage == "selfcal":
            stage_args = [
                self.resolve_calibrated_target_ms(),
                "--outdir",
                str(self.stage_outdir),
                "--upstream-outdir",
                str(self.outdir),
                "--profile",
                self.args.profile,
                "--telescope",
                self.resolved_telescope() or self.args.telescope,
            ]
            receiver_band = (
                self.resolved_receiver_band() or self.args.receiver_band
            )
            if receiver_band:
                stage_args.extend(["--receiver-band", receiver_band])
            if self.args.reference_antenna:
                stage_args.extend(
                    ["--reference-antenna", self.args.reference_antenna]
                )
            if self.args.target:
                stage_args.extend(["--target", self.args.target])
            stage_args.extend([
                "--provenance-fingerprint",
                self.stage_fingerprint("selfcal"),
            ])
            return self.casa_command(script, *stage_args)

        if stage == "imaging":

            if self.imaging_backend() == "wsclean":
                imaging_ms = self.resolve_imaging_ms()
                if not imaging_ms:
                    raise RuntimeError(
                        "WSClean imaging requires an explicit/calibrated MS. "
                        "Pass --imaging-ms, or use --data-state calibrated "
                        "with the target MS."
                    )
                wsclean_telescope = (
                    self.resolved_telescope() or self.args.telescope
                )
                # MeerKAT's primary beam is the measured katbeam model,
                # applied to the result; there is no beam state to choose.
                # Everything else that reaches this backend is LOFAR, whose
                # EveryBeam correction needs one stated explicitly.
                meerkat_beam = wsclean_telescope == "meerkat"
                if not meerkat_beam and not self.args.lofar_beam_mode:
                    raise RuntimeError(
                        "LOFAR imaging requires an explicit beam state. Pass "
                        "--lofar-beam-mode full to trust MS metadata, "
                        "differential only for a known phase-centre-corrected "
                        "MS that needs forcing, or none for quick-look output."
                    )

                stage_args = [
                    imaging_ms,
                    "--outdir",
                    str(self.stage_outdir),
                    "--upstream-outdir",
                    str(self.outdir),
                    "--profile",
                    self.args.profile,
                    "--telescope",
                    wsclean_telescope,
                    "--wsclean",
                    self.args.wsclean,
                    "--threads",
                    self.thread_limit(),
                    *(
                        [] if meerkat_beam else
                        ["--beam-mode", self.args.lofar_beam_mode]
                    ),
                    "--datacolumn",
                    self.wsclean_data_column(wsclean_telescope),
                ]
                if self.args.receiver_band:
                    stage_args.extend(
                        ["--receiver-band", self.args.receiver_band]
                    )
                if self.args.target:
                    stage_args.extend(["--target", self.args.target])
                stage_args.extend([
                    "--provenance-fingerprint",
                    self.stage_fingerprint("imaging"),
                ])

                return self.python_command(script, *stage_args)

            stage_args = []
            imaging_ms = self.resolve_imaging_ms()

            if imaging_ms:

                stage_args.append(imaging_ms)

            stage_args.extend(
                [
                    "--outdir",
                    str(self.stage_outdir),
                    "--upstream-outdir",
                    str(self.outdir),
                    "--profile",
                    self.args.profile,
                    "--telescope",
                    self.resolved_telescope() or self.args.telescope,
                ]
            )

            receiver_band = (
                self.resolved_receiver_band() or self.args.receiver_band
            )
            if receiver_band:
                stage_args.extend(["--receiver-band", receiver_band])

            if self.args.target:

                stage_args.extend(
                    [
                        "--target",
                        self.args.target,
                    ]
                )

            stage_args.extend([
                "--provenance-fingerprint",
                self.stage_fingerprint("imaging"),
            ])

            if self.args.subbands is not None:

                stage_args.extend(
                    [
                        "--subbands",
                        self.args.subbands,
                    ]
                )

            # Sub-band imaging is opt-in: --make-subbands on the command line,
            # or make_subbands: true in config.yaml for a uGMRT run. Anything
            # else, including a bare --subbands N, leaves it off.
            if self.args.no_subbands is True:

                stage_args.append("--no-subbands")

            elif self.args.no_subbands is False:

                stage_args.append("--make-subbands")

            elif (
                self.args.profile != "laptop" and
                self.resolved_telescope() in {None, "ugmrt"} and
                bool(IMAGING_TASK_CONFIG.get("make_subbands", False))
            ):

                stage_args.append("--make-subbands")

            else:

                if self.args.subbands is not None:
                    print(
                        f"Note: --subbands {self.args.subbands} has no "
                        "effect because sub-band imaging is off. Add "
                        "--make-subbands to enable it."
                    )

                stage_args.append("--no-subbands")

            use_mpi = self.stage_uses_mpi(stage)

            if use_mpi:
                stage_args.append("--parallel")

            return self.casa_command(
                script,
                *stage_args,
                mpi=use_mpi,
            )

        if stage == "source_detection":

            imaging_report = read_json(
                self.report_path("imaging"),
                {},
            )
            if (
                self.args.target and
                imaging_report and
                imaging_report.get("target") != self.args.target
            ):
                raise RuntimeError(
                    f"The available imaging report belongs to target "
                    f"{imaging_report.get('target')!r}, not requested source "
                    f"detection target {self.args.target!r}. Run imaging for "
                    "the requested target in a separate output directory."
                )

            stage_args = [
                "--outdir",
                str(self.stage_outdir),
            ]

            if self.args.target:

                stage_args.extend(
                    [
                        "--target",
                        self.args.target,
                    ]
                )

            if not bool(
                SOURCE_DETECTION_TASK_CONFIG.get("enabled", True)
            ):

                stage_args.append("--disable")

            return self.python_command(
                script,
                *stage_args,
            )

        if stage == "ddcal":

            stage_args = []
            dd_ms = self.resolve_dd_ms()

            if dd_ms:

                stage_args.append(dd_ms)

            stage_args.extend(
                [
                    "--outdir",
                    str(self.stage_outdir),
                    "--engine",
                    self.args.dd_engine,
                    "--mode",
                    self.args.dd_mode,
                    "--ncpu",
                    self.args.ncpu,
                    "--provenance-fingerprint",
                    self.stage_fingerprint("ddcal"),
                ]
            )

            if self.args.target:

                stage_args.extend(
                    [
                        "--target",
                        self.args.target,
                    ]
                )

            if self.args.dd_prepare_only:

                stage_args.append("--prepare-only")

            else:

                stage_args.append("--execute")

            if self.args.skip_dd_post_imaging:

                stage_args.append("--skip-imaging")

            else:

                stage_args.append("--post-imaging")

            if self.args.solver_command_template:

                stage_args.extend(
                    [
                        "--solver-command-template",
                        self.args.solver_command_template,
                    ]
                )

            if self.args.apply_command_template:

                stage_args.extend(
                    [
                        "--apply-command-template",
                        self.args.apply_command_template,
                    ]
                )

            return self.python_command(
                script,
                *stage_args,
            )

        raise RuntimeError(
            f"No command builder for stage: {stage}"
        )

    # --------------------------------------------------------
    # Stage selection
    # --------------------------------------------------------

    def selected_stages(self):

        if self.args.stages:

            stages = parse_stage_list(
                self.args.stages
            )

            for stage in stages:

                validate_stage(stage)

            # Run an explicit selection in pipeline order regardless of how it
            # was typed. Dependency checking only blocks a stage whose
            # upstream stage failed *in this run*, so a list such as
            # "imaging,calibration" would previously have executed imaging
            # first and either consumed the previous run's products or failed
            # with an unrelated input-resolution error.
            ordered = sorted(
                dict.fromkeys(stages),
                key=PIPELINE_STAGES.index,
            )

            if ordered != list(dict.fromkeys(stages)):
                print(
                    "Reordered the requested stages into pipeline order: "
                    + ", ".join(ordered)
                )

            stages = ordered

        else:

            if self.args.data_state == "calibrated":
                # A calibrated target MS has already crossed the telescope-
                # specific calibration boundary.  Do not mutate it with the
                # raw-data flag/calibration or generic CASA selfcal stages.
                stages = [
                    "ingest",
                    "inspection",
                    "imaging",
                    "source_detection",
                ]
            else:
                stages = [
                    stage
                    for stage in PIPELINE_STAGES
                    if stage != "ddcal"
                ]

            if self.args.include_ddcal:

                stages.append("ddcal")

        stages = select_stage_window(
            stages,
            from_stage=self.args.from_stage,
            to_stage=self.args.to_stage,
        )

        if (
            "ddcal" in stages and
            not self.args.include_ddcal and
            not self.args.stages
        ):

            stages.remove("ddcal")

        if self.input_kind == "fits":
            needs_import = self.stages_require_primary_ms(stages)
            if needs_import:
                # Ingest is a strict prerequisite regardless of the order in
                # which a custom comma-separated stage list mentioned it.
                stages = [stage for stage in stages if stage != "ingest"]
                stages.insert(0, "ingest")
        elif "ingest" in stages:
            explicitly_requested = bool(
                self.args.stages or
                self.args.from_stage == "ingest" or
                self.args.to_stage == "ingest"
            )
            if explicitly_requested:
                raise ValueError(
                    "The ingest stage accepts visibility FITS or GMRT LTA "
                    "input only; a Measurement Set starts at inspection."
                )
            stages.remove("ingest")

        return stages

    # --------------------------------------------------------
    # Execution
    # --------------------------------------------------------

    def acquire_run_lock(self):
        """Prevent two pipeline jobs from mutating one output tree."""

        if self.args.dry_run:

            return

        lock_path = self.outdir / ".pipeline.lock"
        handle = open(lock_path, "a+")

        try:

            fcntl.flock(
                handle.fileno(),
                fcntl.LOCK_EX | fcntl.LOCK_NB,
            )

        except BlockingIOError as exc:

            handle.seek(0)
            owner = handle.read().strip() or "unknown owner"
            handle.close()
            raise RuntimeError(
                "Another pipeline process is already using this output "
                f"directory:\n{self.outdir}\nLock owner: {owner}"
            ) from exc

        # flock may not reach across nodes of a shared filesystem: a run on
        # the Jupyter node and a SLURM job elsewhere also check the record.
        handle.seek(0)
        owner = live_lock_owner(handle.read(), lock_path)
        if owner:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            handle.close()
            raise RuntimeError(
                "Another pipeline run is already using this output "
                f"directory:\n{self.outdir}\nLock owner: {owner}"
            )

        handle.seek(0)
        handle.truncate()
        json.dump(
            {
                "pid": os.getpid(),
                "host": socket.gethostname(),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                "started_utc": utc_now(),
            },
            handle,
        )
        handle.flush()
        self.run_lock_handle = handle
        self.start_lock_heartbeat(lock_path)

    # --------------------------------------------------------

    def start_lock_heartbeat(self, lock_path):
        """Refresh the lock's time while this run lives (see live_lock_owner)."""

        stop = threading.Event()

        def beat():
            while not stop.wait(LOCK_HEARTBEAT_SECONDS):
                try:
                    os.utime(lock_path)
                except OSError:
                    pass

        threading.Thread(
            target=beat,
            name="pipeline-lock-heartbeat",
            daemon=True,
        ).start()
        self._lock_heartbeat = stop

    # --------------------------------------------------------

    def release_run_lock(self):

        if self.run_lock_handle is None:

            return

        heartbeat = getattr(self, "_lock_heartbeat", None)
        if heartbeat is not None:
            heartbeat.set()
            self._lock_heartbeat = None
        # An empty record: no one holds the directory any more.
        try:
            self.run_lock_handle.seek(0)
            self.run_lock_handle.truncate()
            self.run_lock_handle.flush()
        except OSError:
            pass
        fcntl.flock(
            self.run_lock_handle.fileno(),
            fcntl.LOCK_UN,
        )
        self.run_lock_handle.close()
        self.run_lock_handle = None

    # --------------------------------------------------------

    def state_path(self):

        return self.outdir / "pipeline_state.json"

    # --------------------------------------------------------

    def load_resume_state(self):

        if self.resume_state is not None:

            return self.resume_state

        path = self.state_path()

        if path.exists():

            try:

                state = read_json(path, {})

            except (OSError, ValueError, TypeError) as exc:

                print(
                    "WARNING: pipeline_state.json is unreadable; stage "
                    f"products will be inspected directly: {exc}"
                )
                state = {}

        else:

            state = {}

        if not isinstance(state, dict):

            state = {}

        state.setdefault("schema_version", 1)
        state.setdefault("created_utc", utc_now())
        state.setdefault("stages", {})

        if not isinstance(state["stages"], dict):

            state["stages"] = {}

        self.resume_state = state

        return state

    # --------------------------------------------------------

    def save_resume_state(self):

        if self.args.dry_run:

            return

        state = self.load_resume_state()
        state["updated_utc"] = utc_now()
        state["configuration_sha256"] = self.configuration_digest()
        atomic_write_json(self.state_path(), state)

    # --------------------------------------------------------

    def stage_fingerprint(self, stage):

        if stage in self._stage_fingerprints:
            return self._stage_fingerprints[stage]

        digest = hashlib.sha256()

        fingerprint_paths = [
            self.stage_script(stage),
            self.root / "config.py",
            self.root / "config.yaml",
            self.root / "telescope_profiles.py",
            self.root / "casa_compat.py",
            self.root / "qa_products.py",
        ]
        if stage == "flagging":
            fingerprint_paths.append(self.root / "flagging_selectors.py")
        if stage == "ingest" and self.input_is_lta():
            fingerprint_paths.append(self.root / "lta_convert.py")
        if stage == "imaging":
            fingerprint_paths.append(self.root / "imaging_safety.py")
            if (
                self.imaging_backend() == "wsclean" and
                self.resolved_telescope() == "meerkat"
            ):
                # The beam correction sets MeerKAT's published flux scale.
                fingerprint_paths.append(self.root / "meerkat_beam.py")

        for path in fingerprint_paths:

            if path.exists():

                digest.update(path.name.encode("utf-8"))
                digest.update(path.read_bytes())

        # Telescope/band are represented below in canonical resolved form.
        # Raw aliases (or auto versus an equivalent verified explicit value)
        # must not invalidate otherwise identical products.
        common_argument_names = {"data_state"}
        stage_argument_names = {
            "ingest": {
                "fits_format",
                "imported_ms",
                "telescope",
                "gmrt_flagfile",
            },
            "inspection": {"ms"},
            "flagging": {"ms", "copy_input_ms", "working_ms", "profile"},
            "calibration": {
                "ms",
                "flux_calibrator",
                "bandpass_calibrator",
                "secondary_flux_calibrator",
                "phase_calibrator",
                "target_phase_map",
                "target_fields",
                "reference_antenna",
                "flux_standard",
                "flux_model",
                "flux_model_type",
            },
            "selfcal": {
                "selfcal_ms",
                "reference_antenna",
                "target",
                "profile",
            },
            "imaging": {
                "imaging_ms",
                "target",
                "wsclean",
                "wsclean_data_column",
                "lofar_beam_mode",
                "subbands",
                "no_subbands",
                "profile",
            },
            "source_detection": {"target"},
            "ddcal": {
                "dd_ms",
                "target",
                "dd_mode",
                "dd_engine",
                "dd_prepare_only",
                "skip_dd_post_imaging",
                "solver_command_template",
                "apply_command_template",
            },
        }
        selected_argument_names = (
            common_argument_names | stage_argument_names.get(stage, set())
        )

        stage_input = None
        if stage == "ingest":
            stage_input = str(self.input_data) if self.input_data else None
        elif stage in {"inspection", "flagging", "calibration"}:
            stage_input = self.args.ms
        elif stage == "selfcal":
            stage_input = self.resolve_calibrated_target_ms()
        elif stage == "imaging":
            stage_input = self.resolve_imaging_ms()
        elif stage == "ddcal":
            stage_input = self.resolve_dd_ms()

        stage_input_stat_sha256 = None
        if stage in {"imaging", "ddcal"} and stage_input:
            try:
                stage_input_stat_sha256 = (
                    IMAGING_SAFETY.path_stat_fingerprint(stage_input)
                )
            except (OSError, ValueError):
                stage_input_stat_sha256 = None

        stage_input_instance_identity = None
        if stage != "ingest" and stage_input:
            try:
                stage_input_instance_identity = (
                    self.measurement_set_instance_identity(stage_input)
                )
            except (OSError, ValueError):
                stage_input_instance_identity = None

        input_identity = None
        if stage == "ingest" and self.input_data is not None:
            try:
                stat = self.input_data.stat()
                input_identity = {
                    "size_bytes": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                }
            except OSError:
                input_identity = None

        source_ms_identity = self.source_ms_identity_for_provenance()

        flux_model_sha256 = None
        if stage == "calibration" and self.args.flux_model:
            flux_model_sha256 = sha256_file_or_directory(
                self.args.flux_model
            )

        dd_source_catalog_content_identity = None
        if stage == "ddcal":
            source_report = read_json(
                self.report_path("source_detection"),
                {},
            )
            source_catalog = self.resolve_product_path(
                source_report.get("products", {}).get(
                    "source_catalog_csv"
                )
            )
            if source_catalog is not None and source_catalog.is_file():
                try:
                    dd_source_catalog_content_identity = (
                        IMAGING_SAFETY.file_content_identity(source_catalog)
                    )
                except (OSError, ValueError):
                    dd_source_catalog_content_identity = None

        # Ingestion happens before inspection can discover the telescope,
        # receiver band, or imaging backend.  Binding its fingerprint to
        # those downstream discoveries would make a valid import appear
        # stale on the very next --resume run.  The ingest-specific science
        # arguments above still bind every option that can affect the
        # conversion itself.
        telescope_context = (
            None if stage == "ingest" else self.resolved_telescope()
        )
        receiver_band_context = (
            None if stage == "ingest" else self.resolved_receiver_band()
        )
        imaging_backend_context = (
            None if stage == "ingest" else self.imaging_backend()
        )

        digest.update(
            json.dumps(
                {
                    "stage": stage,
                    "telescope": telescope_context,
                    "receiver_band": receiver_band_context,
                    "data_state": self.args.data_state,
                    "imaging_backend": imaging_backend_context,
                    "stage_input": stage_input,
                    "stage_input_stat_sha256": stage_input_stat_sha256,
                    "stage_input_instance_identity": (
                        stage_input_instance_identity
                    ),
                    "input_identity": input_identity,
                    "source_ms_identity": source_ms_identity,
                    "flux_model_sha256": flux_model_sha256,
                    "dd_source_catalog_content_identity": (
                        dd_source_catalog_content_identity
                    ),
                    "science_arguments": {
                        key: value
                        for key, value in vars(self.args).items()
                        if key in selected_argument_names
                    },
                },
                sort_keys=True,
                default=str,
            ).encode("utf-8")
        )

        upstream_reports = {
            "inspection": ["ingest"],
            "flagging": ["inspection"],
            "calibration": ["inspection", "flagging_post"],
            "selfcal": ["inspection", "calibration", "split_targets"],
            "imaging": ["inspection", "calibration", "selfcal"],
            "source_detection": ["imaging"],
            "ddcal": ["imaging", "source_detection"],
        }
        for report_name in upstream_reports.get(stage, []):
            report_path = self.report_path(report_name)
            if report_path.is_file():
                digest.update(report_name.encode("utf-8"))
                digest.update(report_path.read_bytes())

        self._stage_fingerprints[stage] = digest.hexdigest()
        return self._stage_fingerprints[stage]

    # --------------------------------------------------------

    def stage_recipe_fingerprint(self, stage):
        """Return the science-recipe fingerprint used for this stage.

        Normally the computed value. When a recorded fingerprint has been
        proven to describe the same recipe (see
        ``recorded_recipe_matches_after_provenance_refreshes``) the recorded
        value is used instead, so internal checkpoints bound to it stay valid.
        The pin applies only while the computed value is unchanged.
        """

        value = self.computed_stage_recipe_fingerprint(stage)
        pin = self.__dict__.setdefault("_stage_recipe_pins", {}).get(stage)
        if pin and pin.get("computed") == value:
            return pin["recorded"]
        return value

    # --------------------------------------------------------

    def imaging_recipe_band_configuration(self, receiver_band):
        """The band's imaging settings as they enter the science recipe.

        ``wsclean_channels_out`` configures only WSClean. A CASA image does
        not depend on it, so it stays out of a CASA recipe: otherwise a
        MeerKAT output directory made with CASA before MeerKAT gained WSClean
        settings could not be resumed with ``--imager casa``.
        """

        if CONFIG is None:
            return None
        band = getattr(CONFIG, "FINAL_IMAGING_CONFIG", {}).get(receiver_band)
        if band is None or self.imaging_backend() != "casa":
            return band
        return {
            key: value for key, value in band.items()
            if key != "wsclean_channels_out"
        }

    def imaging_recipe_task_configuration(self, telescope):
        """The imaging task settings that belong in the science recipe.

        A setting that cannot affect a run's science must not change that
        run's recipe, or a release that adds it would make every existing
        output directory refuse to resume. The WSClean memory cap is a
        resource limit and never science. The MeerKAT MFS beam mode shapes
        only MeerKAT's WSClean products.
        """

        excluded = {"wsclean_memory_fraction"}
        if not (
            telescope == "meerkat" and
            self.imaging_backend() == "wsclean"
        ):
            excluded.add("meerkat_mfs_beam_mode")
        return {
            key: value for key, value in IMAGING_TASK_CONFIG.items()
            if key not in excluded
        }

    def computed_stage_recipe_fingerprint(self, stage, report_bytes=None,
                                          calibration_task_overrides=None):
        """Fingerprint science inputs/settings separately from source code.

        Completed products remain bound to ``stage_fingerprint``, including
        every implementation file.  This recipe fingerprint has a narrower
        purpose: a failed calibration may continue from individually
        validated task checkpoints after a code-only bug fix, while any
        changed MS, calibrator, model, band, upstream report, or effective
        calibration setting still fails closed.

        ``report_bytes`` substitutes upstream report content for a calibration
        fingerprint; such values are never cached.
        """

        if (
            report_bytes is not None or calibration_task_overrides
        ) and stage != "calibration":
            raise ValueError(
                "report substitution is supported only for calibration"
            )

        substituted = (
            report_bytes is not None or bool(calibration_task_overrides)
        )
        if (
            not substituted and
            stage in self._stage_recipe_fingerprints
        ):
            return self._stage_recipe_fingerprints[stage]

        if stage == "flagging":
            stage_input = self.args.ms
            stage_input_instance_identity = None
            if stage_input:
                try:
                    stage_input_instance_identity = (
                        self.measurement_set_instance_identity(stage_input)
                    )
                except (OSError, TypeError, ValueError):
                    pass

            telescope = self.resolved_telescope()
            receiver_band = self.resolved_receiver_band()
            payload = {
                "schema_version": 2,
                "stage": "flagging",
                "telescope": telescope,
                "receiver_band": receiver_band,
                "data_state": self.args.data_state,
                "stage_input": stage_input,
                "stage_input_instance_identity": (
                    stage_input_instance_identity
                ),
                "science_arguments": {
                    key: value
                    for key, value in vars(self.args).items()
                    if key in {
                        "profile",
                        "target_fields",
                        "flux_calibrator",
                        "phase_calibrator",
                    }
                },
                "flagging_configuration": (
                    getattr(CONFIG, "FLAGGING_CONFIG", {})
                    if CONFIG is not None else {}
                ),
                "telescope_flagging_configuration": (
                    getattr(CONFIG, "TELESCOPE_FLAGGING_CONFIG", {}).get(
                        telescope,
                        {},
                    )
                    if CONFIG is not None and telescope else {}
                ),
                "profile_flagging_configuration": (
                    self.profile_settings.get("flagging", {})
                ),
            }
            digest = hashlib.sha256()
            digest.update(json.dumps(
                payload,
                sort_keys=True,
                default=str,
            ).encode("utf-8"))
            inspection_path = self.report_path("inspection")
            if inspection_path.is_file():
                digest.update(b"inspection")
                digest.update(inspection_path.read_bytes())
            value = digest.hexdigest()
            self._stage_recipe_fingerprints[stage] = value
            return value

        if stage == "selfcal":
            # Self-calibration needs a recipe of its own. The implementation
            # fingerprint used for the remaining stages changes with every
            # code update, so a killed self-calibration could otherwise never
            # continue from its accepted loops once a fix was deployed.
            stage_input = self.resolve_calibrated_target_ms()
            stage_input_instance_identity = None
            if stage_input:
                try:
                    stage_input_instance_identity = (
                        self.measurement_set_instance_identity(stage_input)
                    )
                except (OSError, TypeError, ValueError):
                    pass

            telescope = self.resolved_telescope()
            receiver_band = self.resolved_receiver_band()
            band_profile = None
            try:
                if telescope and receiver_band:
                    band_profile = TELESCOPE_PROFILE_MODULE.get_band_profile(
                        telescope,
                        receiver_band,
                    )
            except (TypeError, ValueError):
                band_profile = None

            payload = {
                "schema_version": 2,
                "stage": "selfcal",
                "telescope": telescope,
                "receiver_band": receiver_band,
                "data_state": self.args.data_state,
                "stage_input": stage_input,
                "stage_input_instance_identity": (
                    stage_input_instance_identity
                ),
                "science_arguments": {
                    key: value
                    for key, value in vars(self.args).items()
                    if key in {
                        "selfcal_ms",
                        "reference_antenna",
                        "target",
                        "profile",
                    }
                },
                "effective_configuration": {
                    "band": (
                        getattr(CONFIG, "SELFCAL_CONFIG", {}).get(
                            receiver_band
                        )
                        if CONFIG is not None else None
                    ),
                    "imaging_band": (
                        getattr(CONFIG, "SELFCAL_IMAGING_CONFIG", {}).get(
                            receiver_band
                        )
                        if CONFIG is not None else None
                    ),
                    "task": (
                        getattr(CONFIG, "SELFCAL_TASK_CONFIG", {})
                        if CONFIG is not None else {}
                    ),
                    "flat_imaging_overrides": (
                        getattr(CONFIG, "FLAT_IMAGING_OVERRIDES", [])
                        if CONFIG is not None else []
                    ),
                    "execution_profile": self.profile_settings,
                },
                "band_profile": band_profile,
            }
            digest = hashlib.sha256()
            digest.update(json.dumps(
                payload,
                sort_keys=True,
                default=str,
            ).encode("utf-8"))
            for report_name in ["inspection", "calibration", "split_targets"]:
                content = (
                    self.calibration_science_content()
                    if report_name == "calibration" else
                    self.recipe_report_content(report_name)
                )
                if content is not None:
                    digest.update(report_name.encode("utf-8"))
                    digest.update(content)
            value = digest.hexdigest()
            self._stage_recipe_fingerprints[stage] = value
            return value

        if stage not in {"calibration", "imaging"}:
            value = self.stage_fingerprint(stage)
            self._stage_recipe_fingerprints[stage] = value
            return value

        if stage == "imaging":
            stage_input = self.resolve_imaging_ms()
            stage_input_stat_sha256 = None
            stage_input_instance_identity = None
            if stage_input:
                try:
                    stage_input_stat_sha256 = (
                        IMAGING_SAFETY.path_stat_fingerprint(stage_input)
                    )
                    stage_input_instance_identity = (
                        self.measurement_set_instance_identity(stage_input)
                    )
                except (OSError, TypeError, ValueError):
                    pass

            telescope = self.resolved_telescope()
            receiver_band = self.resolved_receiver_band()
            argument_names = {
                "data_state",
                "imaging_ms",
                "target",
                "wsclean",
                "wsclean_data_column",
                "lofar_beam_mode",
                "subbands",
                "no_subbands",
                "profile",
            }
            payload = {
                "schema_version": 1,
                "stage": stage,
                "telescope": telescope,
                "receiver_band": receiver_band,
                "imaging_backend": self.imaging_backend(),
                "stage_input": stage_input,
                "stage_input_stat_sha256": stage_input_stat_sha256,
                "stage_input_instance_identity": (
                    stage_input_instance_identity
                ),
                "science_arguments": {
                    key: value
                    for key, value in vars(self.args).items()
                    if key in argument_names
                },
                "effective_configuration": {
                    "band": self.imaging_recipe_band_configuration(
                        receiver_band
                    ),
                    "task": self.imaging_recipe_task_configuration(
                        telescope
                    ),
                    "telescope_task": (
                        getattr(
                            CONFIG,
                            "TELESCOPE_IMAGING_TASK_CONFIG",
                            {},
                        ).get(telescope)
                        if CONFIG is not None else None
                    ),
                    "execution_profile": self.profile_settings,
                },
            }
            katbeam_wsclean = (
                telescope == "meerkat" and
                self.imaging_backend() == "wsclean"
            )
            if katbeam_wsclean:
                # The katbeam correction weights the MFS beam by WSClean's
                # imaging weights and evaluates each sub-band at the centre
                # of what survived the flagging stage's frequency masks
                # (finding 165). Both set the published flux scale off-axis,
                # so the scheme, and the masks the flagging stage applied,
                # belong to the recipe.
                payload["meerkat_beam_correction"] = {
                    "mfs_beam_weights": "wsclean_image_weight",
                    "subband_frequency": "unflagged_centre",
                }
            digest = hashlib.sha256()
            digest.update(json.dumps(
                payload,
                sort_keys=True,
                default=str,
            ).encode("utf-8"))
            for report_name in ["inspection", "calibration", "selfcal"]:
                report_path = self.report_path(report_name)
                if report_path.is_file():
                    digest.update(report_name.encode("utf-8"))
                    digest.update(report_path.read_bytes())
            mask_report = (
                self.outdir / "flagging" / "known_frequency_mask_report.json"
            )
            if katbeam_wsclean and mask_report.is_file():
                digest.update(b"known_frequency_masks")
                digest.update(mask_report.read_bytes())
            value = digest.hexdigest()
            self._stage_recipe_fingerprints[stage] = value
            return value

        stage_input = self.args.ms
        stage_input_instance_identity = None
        if stage_input:
            try:
                stage_input_instance_identity = (
                    self.measurement_set_instance_identity(stage_input)
                )
            except (OSError, ValueError):
                pass

        source_ms_identity = self.source_ms_identity_for_provenance()

        flux_model_sha256 = None
        if self.args.flux_model:
            flux_model_sha256 = sha256_file_or_directory(
                self.args.flux_model
            )

        telescope = self.resolved_telescope()
        receiver_band = self.resolved_receiver_band()
        telescope_profile = None
        band_profile = None
        try:
            if telescope:
                telescope_profile = (
                    TELESCOPE_PROFILE_MODULE.get_telescope_profile(telescope)
                )
            if telescope and receiver_band:
                band_profile = TELESCOPE_PROFILE_MODULE.get_band_profile(
                    telescope,
                    receiver_band,
                )
        except (TypeError, ValueError):
            telescope_profile = band_profile = None

        argument_names = {
            "data_state",
            "ms",
            "flux_calibrator",
            "bandpass_calibrator",
            "secondary_flux_calibrator",
            "phase_calibrator",
            "target_phase_map",
            "target_fields",
            "reference_antenna",
            "flux_standard",
            "flux_model",
            "flux_model_type",
        }
        effective_configuration = {
            "band": (
                getattr(CONFIG, "CAL_CONFIG", {}).get(receiver_band)
                if CONFIG is not None else None
            ),
            "task": (
                dict(
                    getattr(CONFIG, "CALIBRATION_TASK_CONFIG", {}),
                    **{
                        # Only settings the configuration has: a key it
                        # lacks was absent from the recorded payload too.
                        key: value
                        for key, value in (
                            calibration_task_overrides or {}
                        ).items()
                        if key in getattr(CONFIG, "CALIBRATION_TASK_CONFIG", {})
                    },
                )
                if CONFIG is not None else {}
            ),
            "telescope_task": (
                getattr(CONFIG, "TELESCOPE_CALIBRATION_TASK_CONFIG", {}).get(
                    telescope
                )
                if CONFIG is not None else None
            ),
            "execution_profile": self.profile_settings,
        }

        # Only calibration reaches this payload: flagging, self-calibration
        # and imaging have their own branches above, and the other stages use
        # the implementation fingerprint. The empty key is kept so recorded
        # calibration recipe fingerprints stay valid.
        effective_configuration["stage_specific"] = {}

        payload = {
            "schema_version": 2,
            "stage": stage,
            "telescope": telescope,
            "receiver_band": receiver_band,
            "stage_input": stage_input,
            "stage_input_instance_identity": stage_input_instance_identity,
            "source_ms_identity": source_ms_identity,
            "flux_model_sha256": flux_model_sha256,
            "science_arguments": {
                key: value
                for key, value in vars(self.args).items()
                if key in argument_names
            },
            "effective_configuration": effective_configuration,
            "telescope_profile": telescope_profile,
            "band_profile": band_profile,
        }

        digest = hashlib.sha256()
        digest.update(json.dumps(
            payload,
            sort_keys=True,
            default=str,
        ).encode("utf-8"))
        for report_name in ["inspection", "flagging_post"]:
            if report_bytes is not None and report_name in report_bytes:
                content = report_bytes[report_name]
            else:
                report_path = self.report_path(report_name)
                content = (
                    report_path.read_bytes()
                    if report_path.is_file() else None
                )
            if content is not None:
                digest.update(report_name.encode("utf-8"))
                digest.update(content)

        value = digest.hexdigest()
        if not substituted:
            self._stage_recipe_fingerprints[stage] = value
        return value

    # --------------------------------------------------------

    # Report keys that record which implementation produced or last refreshed
    # a report, rather than anything about the science recipe.
    RECIPE_VOLATILE_REPORT_KEYS = frozenset({
        "provenance_fingerprint",
        "provenance_migration_history",
    })

    def recipe_report_content(self, report_name):
        """Upstream report content for a recipe, without implementation keys."""

        path = self.report_path(report_name)
        if not path.is_file():
            return None
        raw = path.read_bytes()
        try:
            value = json.loads(raw)
        except ValueError:
            return raw
        if isinstance(value, dict):
            value = {
                key: item
                for key, item in value.items()
                if key not in self.RECIPE_VOLATILE_REPORT_KEYS
            }
        return json.dumps(value, sort_keys=True, default=str).encode("utf-8")

    # --------------------------------------------------------

    # What a calibration run determines about the data self-calibration starts
    # from. Everything else in the report is specific to one run and changes
    # whenever the stage is re-entered, even with identical solutions: the
    # implementation and recipe fingerprints, the working MS's file-stat hash
    # (sizes and modification times), checkpoint records, the QA plot path and
    # the flux model's migration history.
    CALIBRATION_SCIENCE_REPORT_KEYS = (
        "original_measurement_set",
        "measurement_set_used",
        "telescope",
        "band",
        "flux_calibrator",
        "flux_calibrators",
        "secondary_flux_calibrator",
        "bandpass_calibrator",
        "phase_calibrator",
        "phase_calibrators",
        "target_phase_map",
        "targets",
        "reference_antenna",
        "bad_antennas",
        "delay_quarantine",
        "hanning_smoothing",
        "calibration_parameters",
        "calibration_tables",
        "calibration_table_qa",
        "split_target_identities",
    )

    def calibration_science_content(self):
        """The calibration report reduced to what defines its solutions."""

        path = self.report_path("calibration")
        if not path.is_file():
            return None
        raw = path.read_bytes()
        try:
            report = json.loads(raw)
        except ValueError:
            return raw
        if not isinstance(report, dict):
            return json.dumps(report, sort_keys=True, default=str).encode(
                "utf-8"
            )

        science = {
            key: report[key]
            for key in self.CALIBRATION_SCIENCE_REPORT_KEYS
            if key in report
        }
        identity = (report.get("flux_model_provenance") or {}).get("identity")
        if isinstance(identity, dict):
            science["flux_model_identity"] = {
                key: value
                for key, value in identity.items()
                if not key.startswith("pipeline_")
            }
        return json.dumps(science, sort_keys=True, default=str).encode("utf-8")

    # --------------------------------------------------------

    def supersede_failed_attempt(
        self,
        stage,
        previous,
        current_fingerprint,
        current_recipe_fingerprint,
        reason,
    ):
        """Restart a failed per-target stage whose recipe has changed.

        Calibration cannot do this: it writes flags and corrected data into
        the working MS, so an attempt under another recipe has already
        changed its own input. Self-calibration and imaging keep every
        product under their own directory. Self-calibration alters its input
        MS only by rewriting MODEL_DATA, which the initial image of a fresh
        run rewrites in full before any gain solve, and imaging does not
        write to its input at all. Setting the old directory aside therefore
        removes all of the old recipe's products while keeping flagging and
        calibration, where the alternative was a new output directory that
        repeats both.
        """

        stage_dir = self.stage_outdir / stage
        archived = None
        if stage_dir.exists():
            index = 1
            while True:
                archived = stage_dir.with_name(f"{stage}.superseded.{index}")
                if not archived.exists():
                    break
                index += 1
            if not self.args.dry_run:
                stage_dir.rename(archived)

        previous.setdefault("superseded_attempts", []).append({
            "superseded_utc": utc_now(),
            "archived_directory": str(archived) if archived else None,
            "previous_attempt_stage_fingerprint": previous.get(
                "attempt_stage_fingerprint"
            ),
            "previous_attempt_recipe_fingerprint": previous.get(
                "attempt_recipe_fingerprint"
            ),
            "current_recipe_fingerprint": current_recipe_fingerprint,
            "reason": reason,
        })
        previous["attempt_stage_fingerprint"] = current_fingerprint
        previous["attempt_recipe_fingerprint"] = current_recipe_fingerprint
        self.save_resume_state()

        if archived and not self.args.dry_run:
            where = f" Its products were moved to {archived}"
        elif archived:
            where = f" A real run would move its products to {archived}"
        else:
            where = " It left no products"
        print(
            f"Stage '{stage}' was previously attempted under a different "
            f"science recipe ({reason}) and did not complete.{where}, and "
            "it will start again. Flagging and calibration are kept.",
            flush=True,
        )

    # --------------------------------------------------------

    def post_flagging_report_before_provenance_refreshes(self):
        """Rebuild the post-flag report as it was before each refresh.

        When flagging is re-entered after it has completed, flagging.py only
        rebinds the report to the current implementation: it records the
        replaced fingerprint in ``provenance_migration_history`` and rewrites
        the file with the same ``json.dump(indent=4)`` call. Undoing those
        refreshes, newest first, therefore reproduces the earlier bytes
        exactly. Returns one candidate per refresh.
        """

        path = self.report_path("flagging_post")
        try:
            report = json.loads(path.read_bytes())
        except (OSError, ValueError):
            return []
        if not isinstance(report, dict):
            return []

        candidates = []
        current = report
        while True:
            history = current.get("provenance_migration_history")
            if not isinstance(history, list) or not history:
                break
            entry = history[-1]
            if not isinstance(entry, dict) or "previous_fingerprint" not in entry:
                break

            previous = dict(current)
            if len(history) > 1:
                previous["provenance_migration_history"] = history[:-1]
            else:
                del previous["provenance_migration_history"]
            previous["provenance_fingerprint"] = entry["previous_fingerprint"]
            candidates.append(
                json.dumps(previous, indent=4, default=str).encode("utf-8")
            )
            if entry["previous_fingerprint"] is None:
                # A report written before provenance was recorded had no key.
                without_key = dict(previous)
                del without_key["provenance_fingerprint"]
                candidates.append(
                    json.dumps(
                        without_key,
                        indent=4,
                        default=str,
                    ).encode("utf-8")
                )
            current = previous

        return candidates

    @staticmethod
    def calibration_failure_policy():
        """The antenna-quarantine settings currently in effect."""

        task = (
            getattr(CONFIG, "CALIBRATION_TASK_CONFIG", {})
            if CONFIG is not None else {}
        )
        return {
            key: task.get(key, default)
            for key, default in CALIBRATION_FAILURE_POLICY_DEFAULTS.items()
        }

    def recorded_recipe_matches_apart_from_failure_policy(
        self,
        stage,
        recorded,
        previous=None,
    ):
        """Prove a calibration recipe that differs only in quarantine policy.

        The recorded SHA-256 is reproduced from the current payload with the
        antenna-quarantine settings put back to the values the attempt ran
        under (recorded since this check exists) or, for attempts recorded
        before it, the values the pipeline has always shipped with. The
        post-flag report's provenance refreshes may be undone at the same
        time. Any other difference still fails. On success the recorded value
        is pinned, as for a provenance refresh, so the calibration
        checkpoints bound to it stay valid.
        """

        if stage != "calibration" or not recorded:
            return False

        current = self.calibration_failure_policy()
        candidates = []
        for key in ("attempt_failure_policy", "product_failure_policy"):
            value = (previous or {}).get(key)
            if isinstance(value, dict):
                candidates.append(dict(current, **value))
        candidates.append(dict(CALIBRATION_FAILURE_POLICY_DEFAULTS))

        reports = [None] + [
            {"flagging_post": candidate}
            for candidate in self.post_flagging_report_before_provenance_refreshes()
        ]
        tried = []
        for policy in candidates:
            if policy == current or policy in tried:
                continue
            tried.append(policy)
            for report_bytes in reports:
                reconstructed = self.computed_stage_recipe_fingerprint(
                    stage,
                    report_bytes=report_bytes,
                    calibration_task_overrides=policy,
                )
                if reconstructed == recorded:
                    self.__dict__.setdefault("_stage_recipe_pins", {})[
                        stage
                    ] = {
                        "computed": self.computed_stage_recipe_fingerprint(
                            stage
                        ),
                        "recorded": recorded,
                        "failure_policy": policy,
                    }
                    return True
        return False

    def recipe_failure_policy(self, stage):
        """The quarantine policy the stage's recipe fingerprint was made with.

        Normally the current settings; while a recorded fingerprint is pinned
        through a policy-only difference, the policy that reproduced it, so
        the state always pairs a fingerprint with the settings behind it.
        """

        pin = self.__dict__.get("_stage_recipe_pins", {}).get(stage) or {}
        if (
            pin.get("failure_policy") and
            pin.get("computed") == self.computed_stage_recipe_fingerprint(stage)
        ):
            return dict(pin["failure_policy"])
        return self.calibration_failure_policy()

    def recorded_recipe_matches_after_provenance_refreshes(
        self,
        stage,
        recorded,
    ):
        """Prove a recorded calibration recipe across flagging refreshes.

        The calibration recipe hashes the post-flag report byte for byte, and
        a re-entered flagging stage rewrites only that report's implementation
        provenance. Such a code-only refresh must not look like a science
        change. The recorded SHA-256 is reproduced from the current payload,
        inspection report and post-flag report with the refreshes undone; any
        other difference still fails. On success the recorded value is pinned
        so the calibration checkpoints bound to it remain valid.
        """

        if stage != "calibration" or not recorded:
            return False

        for candidate in self.post_flagging_report_before_provenance_refreshes():
            reconstructed = self.computed_stage_recipe_fingerprint(
                stage,
                report_bytes={"flagging_post": candidate},
            )
            if reconstructed == recorded:
                self.__dict__.setdefault("_stage_recipe_pins", {})[stage] = {
                    "computed": self.computed_stage_recipe_fingerprint(stage),
                    "recorded": recorded,
                }
                return True

        return False

    # --------------------------------------------------------

    def report_provenance_is_current(self, stage, reported):
        """Whether a report's implementation fingerprint is current.

        Adopting a completed stage after a code update records the new
        fingerprint in the pipeline state but leaves the stage's own report
        untouched. The report is still current when the adoption history
        leads from its fingerprint to the current one; otherwise every
        adopted stage would be judged stale on the following run and
        re-entered, rewriting reports that downstream recipes depend on.
        """

        current = self.stage_fingerprint(stage)
        if not reported:
            return False
        if reported == current:
            return True

        entry = self.load_resume_state().get("stages", {}).get(
            self.stage_state_key(stage),
            {},
        )
        recorded = (
            entry.get("product_stage_fingerprint") or
            entry.get("stage_fingerprint")
        )
        if recorded != current:
            return False

        successors = {}
        for record in entry.get("implementation_adoption_history") or []:
            if not isinstance(record, dict):
                continue
            before = record.get("previous_stage_fingerprint")
            after = record.get("current_stage_fingerprint")
            if before and after:
                successors.setdefault(before, set()).add(after)

        seen = set()
        pending = [reported]
        while pending:
            value = pending.pop()
            if value == current:
                return True
            if value in seen:
                continue
            seen.add(value)
            pending.extend(successors.get(value, ()))
        return False

    # --------------------------------------------------------

    @staticmethod
    def normalized_checkpoint_resume_command(command):
        """Remove implementation provenance tokens from a stage command."""

        ignored = {
            "--provenance-fingerprint",
            "--recipe-fingerprint",
            # The LTA converters' location and whether their FITS is kept.
            "--listscan",
            "--gvfits",
        }
        ignored_flags = {"--keep-lta-fits"}
        normalized = []
        parts = [str(part) for part in (command or [])]
        index = 0
        while index < len(parts):
            part = parts[index]
            if part in ignored_flags:
                index += 1
                continue
            if part in ignored:
                index += 2
                continue
            if any(part.startswith(option + "=") for option in ignored):
                index += 1
                continue
            normalized.append(part)
            index += 1
        return normalized

    # --------------------------------------------------------

    def legacy_calibration_attempt_is_compatible(self, previous):
        """One-time migration for state written before recipe fingerprints."""

        old_command = self.normalized_checkpoint_resume_command(
            previous.get("command")
        )
        new_command = self.normalized_checkpoint_resume_command(
            self.command_for_stage("calibration")
        )
        if not old_command or old_command != new_command:
            return False, "the recorded calibration command changed"

        # An external model cannot be proven unchanged from legacy state,
        # which did not store a separate attempt-recipe fingerprint.
        if self.args.flux_model:
            return False, (
                "legacy state cannot prove an external flux model unchanged"
            )

        started_text = previous.get("started_utc")
        try:
            started = datetime.fromisoformat(
                str(started_text).replace("Z", "+00:00")
            )
            if started.tzinfo is None:
                started = started.replace(tzinfo=timezone.utc)
            started_timestamp = started.timestamp()
        except (TypeError, ValueError, OverflowError):
            return False, "the prior calibration start time is unavailable"

        protected_paths = [
            self.root / "config.py",
            self.root / "config.yaml",
            self.root / "telescope_profiles.py",
            self.report_path("inspection"),
            self.report_path("flagging_post"),
        ]
        for path in protected_paths:
            if path.exists() and path.stat().st_mtime > started_timestamp + 1.0:
                return False, f"science input/configuration changed: {path}"

        return True, (
            "legacy command, input manifest, upstream reports, and science "
            "configuration are unchanged"
        )

    # --------------------------------------------------------

    def legacy_flagging_attempt_is_compatible(self, previous):
        """Prove that a pre-recipe flagging attempt is safe to resume.

        Older pipeline states recorded only the stage implementation
        fingerprint.  Do not infer a recipe from that hash: require the
        normalized command, the durable flagging-parameter snapshot, the
        current MS/report identity, and unchanged science configuration to
        agree.  This permits a code-only upgrade to continue a killed
        flagging stage while still rejecting an edited recipe or input.
        """

        old_command = self.normalized_checkpoint_resume_command(
            previous.get("command")
        )
        try:
            new_command = self.normalized_checkpoint_resume_command(
                self.command_for_stage("flagging")
            )
        except Exception as exc:
            return False, f"the current flagging command cannot be resolved: {exc}"
        if not old_command or old_command != new_command:
            return False, "the recorded flagging command changed"

        parameters_path = self.outdir / "flagging" / "flagging_parameters.json"
        parameters = read_json(parameters_path, None)
        if not isinstance(parameters, dict):
            return False, "the legacy flagging-parameter record is unavailable"

        expected_telescope = self.resolved_telescope()
        recorded_telescope = parameters.get("telescope")
        if (
            expected_telescope and recorded_telescope and
            str(recorded_telescope).strip().lower() != expected_telescope
        ):
            return False, "the recorded flagging telescope changed"

        expected_band = self.resolved_receiver_band()
        recorded_band = parameters.get("band")
        if (
            expected_band and recorded_band and
            str(recorded_band).strip().upper() != str(expected_band).upper()
        ):
            return False, "the recorded flagging receiver band changed"

        if (
            parameters.get("execution_profile") and
            parameters.get("execution_profile") != self.args.profile
        ):
            return False, "the recorded flagging execution profile changed"

        post_report = read_json(
            self.report_path("flagging_post"),
        )
        if isinstance(post_report, dict):
            reported_ms = post_report.get("measurement_set")
            if reported_ms and self.args.ms:
                if Path(reported_ms).expanduser().resolve() != Path(
                    self.args.ms
                ).expanduser().resolve():
                    return False, "the recorded flagging Measurement Set changed"
            if (
                expected_telescope and post_report.get("telescope") and
                str(post_report.get("telescope")).strip().lower() !=
                expected_telescope
            ):
                return False, "the completed flagging telescope changed"
            if (
                expected_band and post_report.get("band") and
                str(post_report.get("band")).strip().upper() !=
                str(expected_band).upper()
            ):
                return False, "the completed flagging receiver band changed"

        started_text = previous.get("started_utc")
        try:
            started = datetime.fromisoformat(
                str(started_text).replace("Z", "+00:00")
            )
            if started.tzinfo is None:
                started = started.replace(tzinfo=timezone.utc)
            started_timestamp = started.timestamp()
        except (TypeError, ValueError, OverflowError):
            return False, "the prior flagging start time is unavailable"

        protected_paths = [
            self.root / "config.py",
            self.root / "config.yaml",
            self.root / "telescope_profiles.py",
            self.report_path("inspection"),
        ]
        for path in protected_paths:
            if path.exists() and path.stat().st_mtime > started_timestamp + 1.0:
                return False, f"science input/configuration changed: {path}"

        return True, (
            "legacy flagging command, parameter snapshot, input context, "
            "and science configuration are unchanged"
        )

    # --------------------------------------------------------

    def legacy_imaging_attempt_is_compatible(self, previous):
        """Prove a failed pre-recipe imaging attempt is safe to continue."""

        old_command = self.normalized_checkpoint_resume_command(
            previous.get("command")
        )
        new_command = self.normalized_checkpoint_resume_command(
            self.command_for_stage("imaging")
        )
        if not old_command or old_command != new_command:
            return False, "the recorded imaging command changed"

        stage_input = self.resolve_imaging_ms()
        try:
            live_input_identity = IMAGING_SAFETY.path_stat_fingerprint(
                stage_input
            )
        except (OSError, TypeError, ValueError) as exc:
            return False, f"the imaging MS cannot be fingerprinted: {exc}"

        selfcal_report = read_json(self.report_path("selfcal"), {})
        reported_final_ms = self.resolve_product_path(
            selfcal_report.get("final_ms")
            if isinstance(selfcal_report, dict)
            else None
        )
        if (
            reported_final_ms is None or
            str(reported_final_ms.resolve()) !=
            str(Path(stage_input).expanduser().resolve()) or
            selfcal_report.get("final_ms_stat_sha256") != live_input_identity
        ):
            return False, (
                "the current imaging MS is not bound to the unchanged "
                "self-calibration report"
            )

        started_text = previous.get("started_utc")
        try:
            started = datetime.fromisoformat(
                str(started_text).replace("Z", "+00:00")
            )
            if started.tzinfo is None:
                started = started.replace(tzinfo=timezone.utc)
            started_timestamp = started.timestamp()
        except (TypeError, ValueError, OverflowError):
            return False, "the prior imaging start time is unavailable"

        protected_paths = [
            self.root / "config.py",
            self.root / "config.yaml",
            self.root / "telescope_profiles.py",
            self.report_path("inspection"),
            self.report_path("calibration"),
            self.report_path("selfcal"),
        ]
        for path in protected_paths:
            if path.exists() and path.stat().st_mtime > started_timestamp + 1.0:
                return False, f"science input/configuration changed: {path}"

        return True, (
            "the recorded command, final self-calibrated MS, upstream "
            "reports, and imaging configuration are unchanged"
        )

    # --------------------------------------------------------

    def stage_state_key(self, stage):
        """Resume-state key for a stage, qualified by its target scope.

        Products for one target keep the unqualified key so that a run made
        before per-target directories existed is still recognised.
        """

        if not self.target_scope_name or self.stage_outdir == self.outdir:
            return stage
        return f"{stage}@{self.target_scope_name}"

    def update_stage_state(self, stage, status, **values):

        state = self.load_resume_state()
        entry = state["stages"].setdefault(self.stage_state_key(stage), {})

        if status == "running":

            entry["attempts"] = int(entry.get("attempts", 0)) + 1
            entry["started_utc"] = utc_now()
            entry.pop("finished_utc", None)
            entry.pop("validation_error", None)

        if (
            status in {"completed", "failed", "stale", "interrupted"} and
            not (
                status == "completed" and
                entry.get("status") == "completed"
            )
        ):

            entry["finished_utc"] = utc_now()

        entry["status"] = status
        current_fingerprint = self.stage_fingerprint(stage)
        current_recipe_fingerprint = self.stage_recipe_fingerprint(stage)
        if status == "running":
            entry["attempt_stage_fingerprint"] = current_fingerprint
            entry["attempt_recipe_fingerprint"] = (
                current_recipe_fingerprint
            )
            if stage == "calibration":
                entry["attempt_failure_policy"] = (
                    self.recipe_failure_policy(stage)
                )
        elif status == "completed":
            if not values.get("inferred", False):
                entry["product_stage_fingerprint"] = current_fingerprint
                entry["product_recipe_fingerprint"] = (
                    current_recipe_fingerprint
                )
                if stage == "calibration":
                    entry["product_failure_policy"] = (
                        self.recipe_failure_policy(stage)
                    )
            else:
                entry.setdefault(
                    "product_stage_fingerprint",
                    current_fingerprint,
                )
        entry["validation_stage_fingerprint"] = current_fingerprint
        entry["validation_recipe_fingerprint"] = (
            current_recipe_fingerprint
        )
        # Retain the old field for readers of schema-v1 state, but never
        # overwrite it when merely validating an existing product.
        if "product_stage_fingerprint" in entry:
            entry.setdefault(
                "stage_fingerprint",
                entry["product_stage_fingerprint"],
            )
        entry.update(values)
        self.save_resume_state()

    # --------------------------------------------------------

    def resolve_product_path(self, value):

        if not value:

            return None

        path = Path(str(value)).expanduser()

        if path.is_absolute():

            return path

        candidates = [
            Path.cwd() / path,
            self.root / path,
            self.stage_outdir / path,
            self.outdir / path,
            self.outdir.parent / path,
        ]

        for candidate in candidates:

            if candidate.exists():

                return candidate.resolve()

        return candidates[0]

    # --------------------------------------------------------

    @staticmethod
    def casa_table_exists(path):

        return bool(
            path is not None and
            path.is_dir() and
            (path / "table.dat").is_file()
        )

    # --------------------------------------------------------

    def read_stage_report(self, name):

        path = self.report_path(name)

        if not path.is_file():

            return None, f"missing report: {path}"

        try:

            report = read_json(path, None)

        except (OSError, ValueError, TypeError) as exc:

            return None, f"invalid report {path}: {exc}"

        if not isinstance(report, dict):

            return None, f"report is not a JSON object: {path}"

        return report, ""

    # --------------------------------------------------------

    def validate_meerkat_wsclean_report(self, products, parameters,
                                        band_profile):
        """Check a MeerKAT WSClean report's primary-beam provenance.

        The counterpart of the CASA backend's checks, on FITS products: the
        beam model and version, the reference frequency, that only Taylor-0
        was corrected, and that the restored image, the beam and the
        corrected image still are the files the correction was made from.
        """

        status = products.get("pbcor_status") or {}
        if not products.get("pbcor_image"):
            return False, "MeerKAT imaging report has no PBcor image product"
        if products.get("pbcor_method") != "katbeam_wsclean_mfs":
            return False, (
                "MeerKAT WSClean imaging report does not record the "
                "katbeam PB-correction method"
            )
        if status.get("katbeam_model") != band_profile.get("katbeam_model"):
            return False, (
                "MeerKAT imaging report has missing or mismatched katbeam "
                "model provenance"
            )
        if str(status.get("katbeam_version", "")).strip().lower() in {
            "", "unknown", "none",
        }:
            return False, (
                "MeerKAT katbeam provenance has no trustworthy package "
                "version"
            )
        try:
            beam_frequency = float(status.get("katbeam_frequency_mhz"))
        except (TypeError, ValueError):
            beam_frequency = 0.0
        if beam_frequency <= 0:
            return False, (
                "MeerKAT imaging report has no valid katbeam reference "
                "frequency"
            )
        if status.get("spectral_terms_corrected") is not False:
            return False, (
                "MeerKAT PB report must explicitly record that only the "
                "flux scale, not spectral terms, was corrected"
            )
        mode = status.get("mfs_beam_mode")
        if (
            mode not in {"band_average", "reference_frequency"} or
            parameters.get("mfs_beam_mode") != mode
        ):
            return False, (
                "MeerKAT imaging report has a missing or inconsistent "
                "MFS beam mode"
            )

        restored_path = self.resolve_product_path(
            products.get("restored_image")
        )
        pbcor_path = self.resolve_product_path(products.get("pbcor_image"))
        pb_path = self.resolve_product_path(status.get("katbeam_pb_image"))
        provenance_path = self.resolve_product_path(
            status.get("pbcor_provenance")
        )
        for label, path in (
            ("restored image", restored_path),
            ("PB-corrected image", pbcor_path),
            ("katbeam primary-beam image", pb_path),
            ("PB-correction provenance", provenance_path),
        ):
            if path is None or not path.is_file():
                return False, f"MeerKAT {label} is missing"

        # Schema 2 weights the MFS beam by WSClean's imaging weights and puts
        # each sub-band where its unflagged data lie (finding 165). Schema 1
        # products are 12% low 35 arcmin out on MeerKAT L-band, so they are
        # not adopted.
        weighting = status.get("mfs_weighting")
        if weighting not in {
            "wsclean_image_weight", "covered_fraction",
            "uniform_band_samples",
        }:
            return False, (
                "MeerKAT imaging report does not record how its MFS beam "
                "was weighted; it predates finding 165 and must be "
                "corrected again"
            )
        provenance = read_json(provenance_path, {})
        expected = {
            "schema_version": 2,
            "provider": "katbeam",
            "method": "katbeam_wsclean_mfs",
            "model": band_profile.get("katbeam_model"),
            "katbeam_version": status.get("katbeam_version"),
            "mfs_beam_mode": mode,
            "mfs_weighting": weighting,
            "pbmin": parameters.get("pbmin"),
            "frequency_mhz": beam_frequency,
            "source_image": str(restored_path.resolve()),
            "pb_image": str(pb_path.resolve()),
            "pbcor_image": str(pbcor_path.resolve()),
        }
        for key, value in expected.items():
            actual = provenance.get(key)
            if key.endswith("_image") and actual:
                actual = str(Path(actual).expanduser().resolve())
            if actual != value:
                return False, (
                    "MeerKAT PB-correction provenance is missing or "
                    f"mismatched for {key}: {provenance_path}"
                )

        for label, artifact, key in (
            ("restored image", restored_path, "source"),
            ("katbeam primary-beam image", pb_path, "pb"),
            ("PB-corrected image", pbcor_path, "pbcor"),
        ):
            recorded = provenance.get(key)
            try:
                live = IMAGING_SAFETY.file_content_identity(artifact)
            except (OSError, ValueError) as exc:
                return False, f"Could not fingerprint MeerKAT {label}: {exc}"
            if not isinstance(recorded, dict) or any(
                recorded.get(field) != live[field]
                for field in ("size_bytes", "sha256")
            ):
                return False, (
                    f"MeerKAT {label} no longer matches its PB-correction "
                    "provenance"
                )
        return True, "MeerKAT WSClean PB-correction provenance is valid"

    def validate_report_instrument(self, report, stage):

        expected_telescope = self.resolved_telescope()
        reported_telescope = (
            TELESCOPE_PROFILE_MODULE.telescope_name_from_report(report)
        )

        if expected_telescope:
            if reported_telescope is None:
                if expected_telescope != "ugmrt":
                    return False, (
                        f"{stage} report has no telescope provenance and "
                        f"cannot be reused for {expected_telescope}"
                    )
            elif reported_telescope != expected_telescope:
                return False, (
                    f"{stage} report telescope {reported_telescope} does "
                    f"not match {expected_telescope}"
                )

        expected_band = self.resolved_receiver_band()
        reported_band = report.get("band")
        if isinstance(reported_band, dict):
            reported_band = (
                reported_band.get("name") or reported_band.get("band")
            )
        if expected_band:
            if not reported_band:
                if expected_telescope != "ugmrt" or self.args.receiver_band:
                    return False, (
                        f"{stage} report has no receiver-band provenance "
                        f"and cannot be reused for {expected_band}"
                    )
            elif str(expected_band).upper() != str(reported_band).upper():
                return False, (
                    f"{stage} report band {reported_band} does not match "
                    f"{expected_band}"
                )

        return True, "instrument context is compatible"

    # --------------------------------------------------------

    @staticmethod
    def postcal_chunk_is_scan_local(parameters):
        """Accept whole-scan chunks, or a recorded memory-bounded chunk.

        The corrected-data detectors must never combine data across scans.
        When a whole scan does not fit in memory the calibration stage bounds
        flagdata's time chunk instead; with combinescans=False such a chunk
        still never crosses a scan boundary. It is accepted only when the
        report's own chunk plan records that deliberate bounded decision.
        """

        ntime = parameters.get("ntime")
        if ntime == "scan":
            return True

        text = str(ntime or "")
        if not text.endswith("s"):
            return False
        try:
            seconds = float(text[:-1])
        except ValueError:
            return False
        if not seconds > 0:
            return False

        plan = parameters.get("flagdata_chunk_plan") or {}
        return (
            isinstance(plan, dict) and
            plan.get("status") == "bounded" and
            plan.get("ntime") == text
        )

    def validate_stage_completion(
        self,
        stage,
        allow_implementation_mismatch=False,
    ):
        """Validate durable products before a stage is skipped on resume.

        ``allow_implementation_mismatch`` relaxes only the embedded pipeline
        implementation fingerprint. All live input, telescope/band, target,
        report, table, image, and content-identity checks still run. It is
        used solely while migrating a completed stage across a code update
        whose recorded command or recipe proves the science request unchanged.
        """

        if stage == "ingest":

            marker = self.outdir / "ingest" / "ingest.done"
            report, reason = self.read_stage_report("ingest")
            if not marker.is_file():
                return False, f"missing completion marker: {marker}"
            if report is None:
                return False, reason
            if report.get("stage") != "ingest" or report.get("status") != "completed":
                return False, "ingest report is not marked completed"

            source = report.get("source_fits", {})
            if self.input_data is None:
                return False, "no visibility FITS input is bound to this run"
            expected_source = str(self.input_data.resolve())
            reported_source = source.get("path")
            if reported_source:
                reported_source = str(
                    Path(reported_source).expanduser().resolve()
                )
            if reported_source != expected_source:
                return False, (
                    "ingest report belongs to a different visibility FITS "
                    f"file: {reported_source}"
                )

            try:
                stat = self.input_data.stat()
            except OSError as exc:
                return False, f"visibility FITS input is unavailable: {exc}"
            if (
                source.get("size_bytes") != stat.st_size or
                source.get("mtime_ns") != stat.st_mtime_ns
            ):
                return False, (
                    "visibility FITS size or modification time changed since "
                    "ingestion"
                )

            reported_sha256 = source.get("sha256")
            if not isinstance(reported_sha256, str) or len(reported_sha256) != 64:
                return False, "ingest report has no valid source SHA-256"
            try:
                live_sha256 = FITS_INGEST.sha256_file(self.input_data)
                after_hash_stat = self.input_data.stat()
            except OSError as exc:
                return False, f"could not fingerprint visibility FITS: {exc}"
            if (
                after_hash_stat.st_size != stat.st_size or
                after_hash_stat.st_mtime_ns != stat.st_mtime_ns
            ):
                return False, "visibility FITS changed while it was fingerprinted"
            if live_sha256.lower() != reported_sha256.lower():
                return False, (
                    "visibility FITS SHA-256 does not match the imported "
                    "Measurement Set provenance"
                )

            requested = self.normalize_fits_format(
                report.get("requested_fits_format", "auto")
            )
            if requested != self.args.fits_format:
                return False, (
                    f"ingest FITS format request {requested!r} does not match "
                    f"{self.args.fits_format!r}"
                )
            resolved_format = self.normalize_fits_format(
                report.get("fits_format")
            )
            if resolved_format == "auto":
                return False, "ingest report did not resolve the FITS dialect"

            requested_telescope = str(
                report.get("requested_telescope", "")
            ).strip().lower()
            if requested_telescope == "gmrt":
                requested_telescope = "ugmrt"
            current_telescope = str(
                self.args.telescope or "auto"
            ).strip().lower()
            if requested_telescope != current_telescope:
                return False, (
                    "ingest telescope request does not match this run: "
                    f"{requested_telescope!r} != {current_telescope!r}"
                )

            imported_telescope = str(
                report.get("telescope", "")
            ).strip().lower()
            if imported_telescope == "gmrt":
                imported_telescope = "ugmrt"
            if (
                current_telescope != "auto" and
                imported_telescope != current_telescope
            ):
                return False, (
                    f"ingest report telescope {imported_telescope!r} does "
                    f"not match {current_telescope!r}"
                )

            expected_task = (
                "importfitsidi"
                if resolved_format == "fitsidi"
                else "importgmrt"
                if imported_telescope == "ugmrt"
                else "importuvfits"
            )
            if report.get("casa_task") != expected_task:
                return False, (
                    f"ingest report CASA task {report.get('casa_task')!r} "
                    f"does not match expected {expected_task!r}"
                )
            try:
                expected_flagfiles = [
                    FITS_INGEST.build_input_identity(Path(path))
                    for path in self.args.gmrt_flagfile
                ]
            except (OSError, ValueError) as exc:
                return False, f"could not fingerprint a GMRT flag file: {exc}"
            if report.get("gmrt_flagfiles", []) != expected_flagfiles:
                return False, (
                    "ingest report GMRT flag-file provenance does not "
                    "match this run"
                )

            measurement_set = self.resolve_product_path(
                report.get("measurement_set")
            )
            expected_ms = Path(self.args.ms).expanduser().resolve()
            if measurement_set is None or measurement_set.resolve() != expected_ms:
                return False, (
                    "ingest report Measurement Set does not match the planned "
                    f"import destination: {expected_ms}"
                )
            if not self.casa_table_exists(measurement_set):
                return False, (
                    "ingested Measurement Set is missing or invalid: "
                    f"{measurement_set}"
                )

            return True, (
                f"{resolved_format} provenance and imported Measurement Set "
                "are valid"
            )

        if stage == "inspection":

            report, reason = self.read_stage_report("inspection")

            if report is None:

                return False, reason

            if not report.get("measurement_set"):

                return False, "inspection report has no measurement_set"

            if self.args.ms:
                reported_ms = Path(
                    report["measurement_set"]
                ).expanduser().resolve()
                expected_ms = Path(self.args.ms).expanduser().resolve()
                if reported_ms != expected_ms:
                    return False, (
                        "inspection report belongs to a different Measurement "
                        f"Set: {reported_ms}"
                    )

            compatible, reason = self.validate_report_instrument(
                report,
                "inspection",
            )
            if not compatible:
                return False, reason

            for filename in [
                "fields.json",
                "antennas.json",
                "qa_plot_status.json",
            ]:

                product = self.outdir / "inspection" / filename

                if not product.is_file():

                    return False, f"missing inspection product: {product}"

                try:

                    read_json(product, None)

                except (OSError, ValueError, TypeError) as exc:

                    return False, (
                        f"invalid inspection product {product}: {exc}"
                    )

            return True, "inspection report is valid"

        if stage == "flagging":

            donefile = self.outdir / "flagging" / "flagging.done"
            report, reason = self.read_stage_report("flagging_post")

            if not donefile.is_file():

                return False, f"missing completion marker: {donefile}"

            if report is None:

                return False, reason

            compatible, context_reason = self.validate_report_instrument(
                report,
                "flagging",
            )
            if not compatible:
                return False, context_reason

            reported_ms = report.get("measurement_set")
            if not reported_ms or (
                Path(reported_ms).expanduser().resolve() !=
                Path(self.args.ms).expanduser().resolve()
            ):
                return False, (
                    "flagging report belongs to a different Measurement Set"
                )

            if (
                not allow_implementation_mismatch and
                not self.report_provenance_is_current(
                    "flagging",
                    report.get("provenance_fingerprint"),
                )
            ):
                return False, (
                    "flagging report has no matching input/recipe "
                    "provenance fingerprint"
                )

            return True, "flagging marker and post-flag report are valid"

        if stage == "calibration":

            checkpoint = self.outdir / "calibration" / "calibration.done"
            report, reason = self.read_stage_report("calibration")

            if not checkpoint.is_file():

                return False, f"missing completion marker: {checkpoint}"

            if report is None:

                return False, reason

            compatible, context_reason = self.validate_report_instrument(
                report,
                "calibration",
            )
            if not compatible:
                return False, context_reason

            reported_input = (
                report.get("original_measurement_set") or
                report.get("measurement_set")
            )
            if not reported_input or (
                Path(reported_input).expanduser().resolve() !=
                Path(self.args.ms).expanduser().resolve()
            ):
                return False, (
                    "calibration report belongs to a different input "
                    "Measurement Set"
                )

            if (
                not allow_implementation_mismatch and
                not self.report_provenance_is_current(
                    "calibration",
                    report.get("provenance_fingerprint"),
                )
            ):
                return False, (
                    "calibration report has no matching input/recipe "
                    "provenance fingerprint"
                )

            calibrated_working_ms = report.get("measurement_set_used")
            if not calibrated_working_ms:
                return False, (
                    "calibration report does not identify the Measurement "
                    "Set on which solutions were applied"
                )
            try:
                current_calibration_ms_identity = (
                    IMAGING_SAFETY.path_stat_fingerprint(
                        calibrated_working_ms
                    )
                )
            except (OSError, ValueError) as exc:
                return False, (
                    "calibration input Measurement Set cannot be "
                    f"fingerprinted: {exc}"
                )
            if (
                report.get("measurement_set_stat_sha256") !=
                current_calibration_ms_identity
            ):
                return False, (
                    "calibration report is not bound to the current "
                    "post-calibration Measurement Set contents"
                )

            calibration_telescope = (
                self.resolved_telescope() or
                TELESCOPE_PROFILE_MODULE.telescope_name_from_report(report)
            )
            if calibration_telescope in {"ugmrt", "meerkat"}:
                telescope_label = {
                    "ugmrt": "uGMRT",
                    "meerkat": "MeerKAT",
                }[calibration_telescope]
                postcal_checkpoint = (
                    f"{calibration_telescope}_postcal_flagging"
                )
                round2_checkpoint = (
                    f"{calibration_telescope}_round2_complete"
                )
                parameters = report.get("calibration_parameters", {})
                checkpoints = report.get("checkpoints", {})
                try:
                    solution_rounds = int(
                        parameters.get("solution_rounds", 0)
                    )
                except (TypeError, ValueError):
                    solution_rounds = 0

                if solution_rounds != 2:
                    return False, (
                        f"{telescope_label} calibration report does not "
                        "record the "
                        "required two solution rounds"
                    )
                if not parameters.get("postcal_flagging"):
                    return False, (
                        f"{telescope_label} calibration report does not "
                        "record guarded "
                        "corrected-data flagging"
                    )
                if parameters.get("prebp_gaintype") != "G":
                    return False, (
                        f"{telescope_label} calibration report does not "
                        "record the required receptor-dependent "
                        "pre-bandpass G solve"
                    )
                for name in [
                    postcal_checkpoint,
                    round2_checkpoint,
                ]:
                    if not checkpoints.get(name):
                        return False, (
                            f"{telescope_label} calibration checkpoint is "
                            "missing: "
                            f"{name}"
                        )
                    durable_marker = (
                        self.outdir / "calibration" / f"{name}.done"
                    )
                    if not durable_marker.is_file():
                        return False, (
                            f"MeerKAT durable checkpoint is missing: {name}"
                        )

                postcal_report = read_json(
                    self.outdir /
                    "calibration" /
                    f"{postcal_checkpoint}_report.json",
                    {},
                )
                if postcal_report.get("status") != "completed":
                    return False, (
                        f"{telescope_label} corrected-data flagging report "
                        "is missing "
                        "or incomplete"
                    )
                if not str(postcal_report.get("correlation", "")).startswith(
                    "ABS_"
                ):
                    return False, (
                        f"{telescope_label} corrected-data flagging report "
                        "does not record an absolute parallel-hand selector"
                    )
                postcal_parameters = postcal_report.get("parameters", {})
                if (
                    postcal_parameters.get("datacolumn") != "CORRECTED" or
                    not self.postcal_chunk_is_scan_local(postcal_parameters) or
                    postcal_parameters.get("combinescans") is not False or
                    postcal_parameters.get("extendflags") is not False
                ):
                    return False, (
                        f"{telescope_label} corrected-data flagging report "
                        "does not preserve the scan-local, non-extending "
                        "detector configuration"
                    )

                extend = postcal_parameters.get("extend")
                if isinstance(extend, dict) and extend.get("applied"):
                    # Growing flags around what the detectors found is a
                    # separate pass, and it must stay inside one scan too.
                    if extend.get("combinescans") is not False:
                        return False, (
                            f"{telescope_label} corrected-data flag "
                            "extension is not scan-local"
                        )

                round2_report = read_json(
                    self.outdir /
                    "calibration" /
                    f"{calibration_telescope}_round2_report.json",
                    {},
                )
                if round2_report.get("status") != "completed":
                    return False, (
                        f"{telescope_label} round-2 calibration report is "
                        "missing or incomplete"
                    )

                for directory in [
                    self.outdir / "calibration",
                    self.outdir / "calibration" / "round1",
                ]:
                    for filename in [
                        "delay.cal",
                        "prebp.cal",
                        "bandpass.cal",
                        "phase.cal",
                        "phase_scan.cal",
                        "amp.cal",
                        "flux.cal",
                    ]:
                        table_path = directory / filename
                        if not self.casa_table_exists(table_path):
                            return False, (
                                f"{telescope_label} round-1/round-2 "
                                "calibration table "
                                f"is missing: {table_path}"
                            )

            split_report, split_reason = self.read_stage_report(
                "split_targets"
            )

            if split_report is None:

                return False, split_reason

            if not split_report:

                return False, "split_targets.json contains no target MS"

            split_identity_report = report.get("split_target_identities", {})
            recorded_split_targets = (
                split_identity_report.get("targets", {})
                if isinstance(split_identity_report, dict)
                else {}
            )
            if (
                not isinstance(recorded_split_targets, dict) or
                set(recorded_split_targets) != set(split_report)
            ):
                return False, (
                    "calibration report has no exact split-target "
                    "Measurement Set identity inventory"
                )

            for target, product in split_report.items():

                path = self.resolve_product_path(product)

                if not self.casa_table_exists(path):

                    return False, (
                        f"calibrated target MS is missing for {target}: "
                        f"{product}"
                    )

                recorded = recorded_split_targets.get(target, {})
                try:
                    current_instance = (
                        self.measurement_set_instance_identity(path)
                    )
                except (OSError, ValueError) as exc:
                    return False, (
                        f"calibrated target MS identity is unreadable for "
                        f"{target}: {exc}"
                    )
                if (
                    not isinstance(recorded, dict) or
                    Path(recorded.get("path", "")).expanduser().resolve() !=
                    path.resolve() or
                    recorded.get("instance_identity") != current_instance
                ):
                    return False, (
                        "calibrated target MS was replaced or is not bound "
                        f"to the calibration report for {target}: {path}"
                    )

            return True, "calibration report and target MS products are valid"

        if stage == "selfcal":

            checkpoint = (
                self.stage_outdir / "selfcal" / "checkpoints" / "selfcal.done"
            )
            report, reason = self.read_stage_report("selfcal")

            if not checkpoint.is_file():

                return False, f"missing completion marker: {checkpoint}"

            if report is None:

                return False, reason

            compatible, context_reason = self.validate_report_instrument(
                report,
                "selfcal",
            )
            if not compatible:
                return False, context_reason

            if (
                self.args.target and
                report.get("target") != self.args.target
            ):
                return False, (
                    f"selfcal report target {report.get('target')!r} does "
                    f"not match {self.args.target!r}"
                )

            try:
                expected_input = Path(
                    self.resolve_calibrated_target_ms()
                ).expanduser().resolve()
            except Exception as exc:
                return False, (
                    "could not resolve the current self-calibration input "
                    f"MS: {exc}"
                )
            reported_input = self.resolve_product_path(
                report.get("input_measurement_set") or
                report.get("measurement_set")
            )
            if (
                reported_input is None or
                reported_input.expanduser().resolve() != expected_input
            ):
                return False, (
                    "selfcal report input MS does not match the current "
                    f"calibrated target MS: {reported_input} != "
                    f"{expected_input}"
                )

            if (
                not allow_implementation_mismatch and
                not self.report_provenance_is_current(
                    "selfcal",
                    report.get("provenance_fingerprint"),
                )
            ):
                return False, (
                    "selfcal report provenance fingerprint does not match "
                    "the current input, recipe, and upstream products"
                )

            for key in ["final_ms", "final_image"]:

                product = report.get(key)
                path = self.resolve_product_path(product)

                if not self.casa_table_exists(path):

                    return False, f"selfcal {key} is missing: {product}"

            try:
                input_identity = IMAGING_SAFETY.path_stat_fingerprint(
                    expected_input
                )
                final_identity = IMAGING_SAFETY.path_stat_fingerprint(
                    self.resolve_product_path(report.get("final_ms"))
                )
            except (OSError, TypeError, ValueError) as exc:
                return False, f"could not fingerprint selfcal MS product: {exc}"
            if report.get("input_measurement_set_stat_sha256") != input_identity:
                return False, "selfcal input MS no longer matches its provenance"
            if report.get("final_ms_stat_sha256") != final_identity:
                return False, "selfcal final MS no longer matches its provenance"

            return True, "selfcal report, MS, and final image are valid"

        if stage == "imaging":

            checkpoint = (
                self.stage_outdir / "imaging" / "checkpoints" / "imaging.done"
            )
            report, reason = self.read_stage_report("imaging")

            if not checkpoint.is_file():

                return False, f"missing completion marker: {checkpoint}"

            if report is None:

                return False, reason

            compatible, context_reason = self.validate_report_instrument(
                report,
                "imaging",
            )
            if not compatible:
                return False, context_reason

            if (
                not allow_implementation_mismatch and
                not self.report_provenance_is_current(
                    "imaging",
                    report.get("provenance_fingerprint"),
                )
            ):
                return False, (
                    "imaging report provenance fingerprint does not match "
                    "the current input, recipe, and upstream products"
                )

            if (
                self.args.target and
                report.get("target") != self.args.target
            ):
                return False, (
                    f"imaging report target {report.get('target')!r} does "
                    f"not match {self.args.target!r}"
                )

            products = report.get("products", {})
            parameters = report.get("imaging_parameters", {})
            backend = report.get("imaging_backend", "casa")
            imaging_telescope = (
                self.resolved_telescope() or
                TELESCOPE_PROFILE_MODULE.telescope_name_from_report(report)
            )
            reported_band = report.get("band")
            if isinstance(reported_band, dict):
                reported_band = (
                    reported_band.get("name") or
                    reported_band.get("band")
                )

            if self.args.imager != "auto":
                expected_backend = self.args.imager
            elif imaging_telescope:
                expected_backend = TELESCOPE_PROFILE_MODULE.default_imager(
                    imaging_telescope
                )
            else:
                expected_backend = "casa"
            if backend != expected_backend:
                hint = ""
                if (
                    imaging_telescope == "meerkat" and
                    backend == "casa" and
                    expected_backend == "wsclean"
                ):
                    hint = (
                        ". MeerKAT now defaults to WSClean; pass --imager "
                        "casa to continue an output directory imaged with "
                        "CASA"
                    )
                return False, (
                    f"imaging backend {backend!r} does not match the "
                    f"required {expected_backend!r} backend for "
                    f"{imaging_telescope}{hint}"
                )

            try:
                band_profile = TELESCOPE_PROFILE_MODULE.get_band_profile(
                    imaging_telescope,
                    reported_band,
                )
            except (TypeError, ValueError):
                return False, (
                    "imaging report has no valid telescope-scoped receiver "
                    f"band: {reported_band!r}"
                )

            expected_provider = band_profile.get("primary_beam_provider")
            reported_provider = parameters.get("primary_beam_provider")
            beam_mode = parameters.get("beam_mode")

            if imaging_telescope == "lofar":
                requested_beam = self.args.lofar_beam_mode
                if requested_beam and beam_mode != requested_beam:
                    return False, (
                        f"LOFAR imaging beam mode {beam_mode!r} does not "
                        f"match requested mode {requested_beam!r}"
                    )
                if beam_mode == "none":
                    expected_provider = "none"
                    pb_required = False
                elif beam_mode in {"full", "differential"}:
                    expected_provider = "wsclean_everybeam"
                    pb_required = True
                else:
                    return False, (
                        "LOFAR imaging report has no explicit full, "
                        "differential, or none beam mode"
                    )
            elif imaging_telescope == "meerkat":
                if expected_provider != "katbeam":
                    return False, (
                        f"No validated final MeerKAT primary-beam provider "
                        f"is available for {reported_band}"
                    )
                pb_required = True
            else:
                pb_required = True

            if reported_provider != expected_provider:
                return False, (
                    f"imaging report PB provider {reported_provider!r} does "
                    f"not match required provider {expected_provider!r}"
                )

            if imaging_telescope == "meerkat" and backend == "wsclean":
                meerkat_valid, meerkat_reason = (
                    self.validate_meerkat_wsclean_report(
                        products,
                        parameters,
                        band_profile,
                    )
                )
                if not meerkat_valid:
                    return False, meerkat_reason

            elif imaging_telescope == "meerkat":
                pb_status = products.get("pbcor_status", {})
                if not products.get("pbcor_image"):
                    return False, (
                        "MeerKAT imaging report has no PBcor image product"
                    )
                if products.get("pbcor_method") != "katbeam_impbcor_tt0":
                    return False, (
                        "MeerKAT imaging report does not record the katbeam "
                        "Taylor-0 PB-correction method"
                    )
                if pb_status.get("katbeam_model") != band_profile.get(
                    "katbeam_model"
                ):
                    return False, (
                        "MeerKAT imaging report has missing or mismatched "
                        "katbeam model provenance"
                    )
                try:
                    beam_frequency = float(
                        pb_status.get("katbeam_frequency_mhz")
                    )
                except (TypeError, ValueError):
                    beam_frequency = 0.0
                if beam_frequency <= 0:
                    return False, (
                        "MeerKAT imaging report has no valid katbeam "
                        "reference frequency"
                    )
                if pb_status.get("spectral_terms_corrected") is not False:
                    return False, (
                        "MeerKAT PB report must explicitly record that only "
                        "Taylor-0, not spectral terms, was corrected"
                    )

                katbeam_pb = self.resolve_product_path(
                    pb_status.get("katbeam_pb_image")
                )
                if not self.casa_table_exists(katbeam_pb):
                    return False, (
                        "MeerKAT katbeam primary-beam artifact is missing or "
                        f"invalid: {pb_status.get('katbeam_pb_image')}"
                    )
                provenance_path = Path(
                    str(katbeam_pb) + ".provenance.json"
                )
                provenance = read_json(provenance_path, {})
                if provenance.get("schema_version") != 2:
                    return False, (
                        "MeerKAT katbeam provenance is missing the current "
                        f"schema: {provenance_path}"
                    )
                if str(provenance.get("katbeam_version", "")).strip().lower() in {
                    "",
                    "unknown",
                    "none",
                }:
                    return False, (
                        "MeerKAT katbeam provenance has no trustworthy "
                        "package version"
                    )
                restored_path = self.resolve_product_path(
                    products.get("restored_image")
                )
                provenance_expected = {
                    "provider": "katbeam",
                    "model": band_profile.get("katbeam_model"),
                    "source_image": (
                        str(restored_path.resolve())
                        if restored_path is not None
                        else None
                    ),
                    "imsize": parameters.get("imsize"),
                    "cell": parameters.get("cell"),
                    "pb_image": str(katbeam_pb.resolve()),
                }
                for key, value in provenance_expected.items():
                    actual = provenance.get(key)
                    if key.endswith("_image") and actual:
                        actual = str(Path(actual).expanduser().resolve())
                    if actual != value:
                        return False, (
                            "MeerKAT katbeam provenance is missing or "
                            f"mismatched for {key}: {provenance_path}"
                        )
                try:
                    provenance_frequency = float(
                        provenance.get("frequency_mhz")
                    )
                except (TypeError, ValueError):
                    provenance_frequency = 0.0
                if provenance_frequency != beam_frequency:
                    return False, (
                        "MeerKAT katbeam provenance reference frequency "
                        "does not match the imaging report"
                    )

                source_geometry = provenance.get("source_geometry")
                pb_geometry = provenance.get("pb_geometry")
                if (
                    not isinstance(source_geometry, dict) or
                    not isinstance(pb_geometry, dict) or
                    not IMAGING_SAFETY.image_geometries_match(
                        source_geometry,
                        pb_geometry,
                    )
                ):
                    return False, (
                        "MeerKAT katbeam provenance has missing or mismatched "
                        "source/PB geometry"
                    )
                if not str(provenance.get("katbeam_version") or "").strip():
                    return False, (
                        "MeerKAT katbeam provenance has no package version"
                    )

                for label, artifact, prefix in [
                    ("restored image", restored_path, "source"),
                    ("katbeam PB", katbeam_pb, "pb"),
                ]:
                    try:
                        matches, detail = (
                            IMAGING_SAFETY.artifact_matches_provenance(
                                provenance,
                                prefix,
                                artifact,
                            )
                        )
                    except (OSError, TypeError, ValueError) as exc:
                        return False, (
                            f"Could not fingerprint MeerKAT {label}: {exc}"
                        )
                    if not matches:
                        return False, (
                            f"MeerKAT {label} no longer matches its katbeam "
                            f"provenance: {detail}"
                        )

                pbcor_path = self.resolve_product_path(
                    products.get("pbcor_image")
                )
                corrected_provenance_path = self.resolve_product_path(
                    pb_status.get("pbcor_provenance")
                )
                if (
                    corrected_provenance_path is None or
                    not corrected_provenance_path.is_file()
                ):
                    return False, (
                        "MeerKAT PB-corrected image provenance is missing: "
                        f"{pb_status.get('pbcor_provenance')}"
                    )
                corrected_provenance = read_json(
                    corrected_provenance_path,
                    {},
                )
                corrected_expected = {
                    "schema_version": 1,
                    "provider": "katbeam",
                    "method": "katbeam_impbcor_tt0",
                    "model": band_profile.get("katbeam_model"),
                    "source_image": str(restored_path.resolve()),
                    "pb_image": str(katbeam_pb.resolve()),
                    "pbcor_image": str(pbcor_path.resolve()),
                    "pbmin": parameters.get("pbmin"),
                }
                for key, value in corrected_expected.items():
                    actual = corrected_provenance.get(key)
                    if key.endswith("_image") and actual:
                        actual = str(Path(actual).expanduser().resolve())
                    if actual != value:
                        return False, (
                            "MeerKAT PB-corrected image provenance is missing "
                            f"or mismatched for {key}: "
                            f"{corrected_provenance_path}"
                        )

                for label, artifact, prefix in [
                    ("restored image", restored_path, "source"),
                    ("katbeam PB", katbeam_pb, "pb"),
                    ("PB-corrected image", pbcor_path, "pbcor"),
                ]:
                    try:
                        matches, detail = (
                            IMAGING_SAFETY.artifact_matches_provenance(
                                corrected_provenance,
                                prefix,
                                artifact,
                            )
                        )
                    except (OSError, TypeError, ValueError) as exc:
                        return False, (
                            f"Could not fingerprint MeerKAT {label}: {exc}"
                        )
                    if not matches:
                        return False, (
                            f"MeerKAT {label} no longer matches the "
                            f"PB-correction provenance: {detail}"
                        )

            if imaging_telescope == "ugmrt":
                pb_status = products.get("pbcor_status", {})
                expected_coefficients = band_profile.get(
                    "primary_beam_coefficients"
                )
                expected_range = band_profile.get(
                    "primary_beam_valid_frequency_mhz"
                )
                if products.get("pbcor_method") != (
                    "ugmrt_polynomial_impbcor_tt0"
                ):
                    return False, (
                        "uGMRT imaging report does not record measured "
                        "polynomial Taylor-0 PB correction"
                    )
                if (
                    pb_status.get("ugmrt_pb_model") !=
                    "NCRA-uGMRT-PB-2023-11-29" or
                    pb_status.get("ugmrt_pb_coefficients") !=
                    expected_coefficients or
                    pb_status.get("ugmrt_pb_valid_frequency_mhz") !=
                    expected_range
                ):
                    return False, (
                        "uGMRT primary-beam model provenance is missing or "
                        "does not match the telescope profile"
                    )
                if pb_status.get("spectral_terms_corrected") is not False:
                    return False, (
                        "uGMRT PB report must explicitly record that only "
                        "Taylor-0, not spectral terms, was corrected"
                    )

                restored_path = self.resolve_product_path(
                    products.get("restored_image")
                )
                pb_path = self.resolve_product_path(
                    pb_status.get("ugmrt_pb_image")
                )
                pbcor_path = self.resolve_product_path(
                    products.get("pbcor_image")
                )
                if not self.casa_table_exists(restored_path):
                    return False, "uGMRT restored image is missing"
                if not self.casa_table_exists(pbcor_path):
                    return False, "uGMRT PB-corrected image is missing"
                if not self.casa_table_exists(pb_path):
                    return False, (
                        "uGMRT measured primary-beam artifact is missing: "
                        f"{pb_status.get('ugmrt_pb_image')}"
                    )

                beam_provenance_path = Path(
                    str(pb_path) + ".provenance.json"
                )
                beam_provenance = read_json(beam_provenance_path, {})
                expected_beam_values = {
                    "schema_version": 1,
                    "provider": "ugmrt_polynomial",
                    "model": "NCRA-uGMRT-PB-2023-11-29",
                    "band": reported_band,
                    "coefficients": expected_coefficients,
                    "valid_frequency_mhz": expected_range,
                    "source_image": str(restored_path.resolve()),
                    "pb_image": str(pb_path.resolve()),
                    "imsize": parameters.get("imsize"),
                    "cell": parameters.get("cell"),
                }
                for key, value in expected_beam_values.items():
                    actual = beam_provenance.get(key)
                    if key.endswith("_image") and actual:
                        actual = str(Path(actual).expanduser().resolve())
                    if actual != value:
                        return False, (
                            "uGMRT PB provenance is missing or mismatched "
                            f"for {key}: {beam_provenance_path}"
                        )
                source_geometry = beam_provenance.get("source_geometry")
                pb_geometry = beam_provenance.get("pb_geometry")
                if (
                    not isinstance(source_geometry, dict) or
                    not isinstance(pb_geometry, dict) or
                    not IMAGING_SAFETY.image_geometries_match(
                        source_geometry,
                        pb_geometry,
                    )
                ):
                    return False, (
                        "uGMRT source/PB geometry provenance is missing or "
                        "mismatched"
                    )

                try:
                    beam_frequency = float(
                        beam_provenance.get("frequency_mhz")
                    )
                    frequency_valid = (
                        float(expected_range[0]) <= beam_frequency <=
                        float(expected_range[1])
                    )
                except (TypeError, ValueError, IndexError):
                    frequency_valid = False
                if not frequency_valid:
                    return False, (
                        "uGMRT PB reference frequency is outside the "
                        "measured model range"
                    )

                corrected_provenance_path = self.resolve_product_path(
                    pb_status.get("pbcor_provenance")
                )
                if (
                    corrected_provenance_path is None or
                    not corrected_provenance_path.is_file()
                ):
                    return False, (
                        "uGMRT PB-corrected image provenance is missing"
                    )
                corrected_provenance = read_json(
                    corrected_provenance_path,
                    {},
                )
                expected_corrected_values = {
                    "schema_version": 1,
                    "provider": "ugmrt_polynomial",
                    "method": "ugmrt_polynomial_impbcor_tt0",
                    "model": "NCRA-uGMRT-PB-2023-11-29",
                    "band": reported_band,
                    "source_image": str(restored_path.resolve()),
                    "pb_image": str(pb_path.resolve()),
                    "pbcor_image": str(pbcor_path.resolve()),
                    "pbmin": parameters.get("pbmin"),
                }
                for key, value in expected_corrected_values.items():
                    actual = corrected_provenance.get(key)
                    if key.endswith("_image") and actual:
                        actual = str(Path(actual).expanduser().resolve())
                    if actual != value:
                        return False, (
                            "uGMRT PB-corrected image provenance is missing "
                            f"or mismatched for {key}: "
                            f"{corrected_provenance_path}"
                        )

                for label, artifact, prefix in [
                    ("restored image", restored_path, "source"),
                    ("primary beam", pb_path, "pb"),
                ]:
                    try:
                        matches, detail = (
                            IMAGING_SAFETY.artifact_matches_provenance(
                                beam_provenance,
                                prefix,
                                artifact,
                            )
                        )
                        if matches:
                            corrected_matches, corrected_detail = (
                                IMAGING_SAFETY.artifact_matches_provenance(
                                    corrected_provenance,
                                    prefix,
                                    artifact,
                                )
                            )
                    except (OSError, TypeError, ValueError) as exc:
                        return False, f"Could not fingerprint uGMRT {label}: {exc}"
                    if not matches:
                        return False, (
                            f"uGMRT {label} no longer matches PB provenance: "
                            f"{detail}"
                        )
                    if not corrected_matches:
                        return False, (
                            f"uGMRT {label} no longer matches PB-correction "
                            f"provenance: {corrected_detail}"
                        )
                try:
                    matches, detail = (
                        IMAGING_SAFETY.artifact_matches_provenance(
                            corrected_provenance,
                            "pbcor",
                            pbcor_path,
                        )
                    )
                except (OSError, TypeError, ValueError) as exc:
                    return False, (
                        f"Could not fingerprint uGMRT PB-corrected image: {exc}"
                    )
                if not matches:
                    return False, (
                        "uGMRT PB-corrected image no longer matches its "
                        f"provenance: {detail}"
                    )

            expected_ms = self.resolve_imaging_ms()
            reported_ms = report.get("measurement_set")
            if expected_ms and reported_ms:
                if (
                    Path(expected_ms).expanduser().resolve() !=
                    Path(reported_ms).expanduser().resolve()
                ):
                    return False, (
                        "imaging report Measurement Set does not match the "
                        "currently resolved target input"
                    )
                try:
                    current_ms_identity = IMAGING_SAFETY.path_stat_fingerprint(
                        expected_ms
                    )
                except (OSError, ValueError) as exc:
                    return False, f"could not fingerprint imaging input MS: {exc}"
                if (
                    report.get("measurement_set_stat_sha256") !=
                    current_ms_identity
                ):
                    return False, (
                        "imaging input MS no longer matches report provenance"
                    )

            keys = ["restored_image"]
            if pb_required or products.get("pbcor_image"):
                keys.append("pbcor_image")

            for key in keys:

                product = products.get(key)
                path = self.resolve_product_path(product)

                exists = (
                    path is not None and path.is_file()
                    if backend == "wsclean"
                    else self.casa_table_exists(path)
                )

                if not exists:

                    return False, f"imaging {key} is missing: {product}"

            fits_products = products.get("fits", [])

            if not fits_products:

                return False, "imaging report contains no FITS products"

            fits_valid, fits_reason = (
                IMAGING_SAFETY.validate_file_content_identities(
                    fits_products,
                    products.get("fits_content_identities"),
                )
            )
            if not fits_valid:
                return False, fits_reason

            for product in fits_products:

                path = self.resolve_product_path(product)

                if path is None or not path.is_file():

                    return False, f"imaging FITS product is missing: {product}"

            return True, "imaging report and primary image products are valid"

        if stage == "source_detection":

            checkpoint = (
                self.stage_outdir /
                "source_detection" /
                "checkpoints" /
                "source_detection.done"
            )
            report, reason = self.read_stage_report("source_detection")

            if not checkpoint.is_file():

                return False, f"missing completion marker: {checkpoint}"

            if report is None:

                return False, reason

            compatible, context_reason = self.validate_report_instrument(
                report,
                "source detection",
            )
            if not compatible:
                return False, context_reason

            imaging_report, imaging_reason = self.read_stage_report(
                "imaging"
            )
            if imaging_report is None:
                return False, imaging_reason

            expected_target = (
                self.args.target or imaging_report.get("target")
            )
            if report.get("target") != expected_target:
                return False, (
                    f"source-detection report target "
                    f"{report.get('target')!r} does not match current "
                    f"imaging target {expected_target!r}"
                )
            if report.get("imaging_backend", "casa") != imaging_report.get(
                "imaging_backend",
                "casa",
            ):
                return False, (
                    "source-detection report imaging backend does not match "
                    "the current imaging report"
                )

            status = report.get("status")

            if status in {"disabled", "skipped"}:

                if not report.get("reason"):

                    return False, (
                        "source-detection report records a skip without a "
                        "reason"
                    )

                return True, f"source detection {status}: {report['reason']}"

            if status != "success":

                return False, (
                    "source-detection report status is not successful: "
                    f"{status}"
                )

            source_input = self.resolve_product_path(
                report.get("input_image")
            )
            if source_input is None or not source_input.is_file():
                return False, (
                    "source-detection input image is missing: "
                    f"{report.get('input_image')}"
                )

            preference = report.get("parameters", {}).get(
                "image_preference",
                "pbcor",
            )
            role_key = (
                "pbcor_image" if preference == "pbcor"
                else "restored_image"
            )
            imaging_products = imaging_report.get("products", {})
            role_product = imaging_products.get(role_key)
            if not role_product:
                return False, (
                    f"current imaging report has no {role_key} for the "
                    "recorded source-detection preference"
                )

            expected_images = []
            direct_role = self.resolve_product_path(role_product)
            if direct_role is not None and direct_role.is_file():
                expected_images.append(direct_role.resolve())

            expected_fits = (
                self.stage_outdir /
                "imaging" /
                "fits" /
                f"{Path(str(role_product)).name}.fits"
            )
            if expected_fits.is_file():
                expected_images.append(expected_fits.resolve())

            for fits_product in imaging_products.get("fits", []) or []:
                fits_path = self.resolve_product_path(fits_product)
                if (
                    fits_path is not None and
                    fits_path.is_file() and
                    fits_path.name == expected_fits.name
                ):
                    expected_images.append(fits_path.resolve())

            if not expected_images:
                return False, (
                    f"no current FITS image resolves the imaging role "
                    f"{role_key}"
                )
            if source_input.resolve() not in set(expected_images):
                return False, (
                    "source-detection catalogue was made from a different "
                    f"image than the current {role_key}"
                )

            if preference == "pbcor":
                source_pb_method = report.get("primary_beam", {}).get(
                    "method"
                )
                current_pb_method = imaging_products.get("pbcor_method")
                if source_pb_method != current_pb_method:
                    return False, (
                        "source-detection PB method does not match the "
                        "current imaging report"
                    )

            products = report.get("products", {})

            product_paths = [
                value
                for value in products.values()
                if value
            ]
            product_identities = report.get("product_content_identities")
            identities_valid, identities_reason = (
                IMAGING_SAFETY.validate_file_content_identities(
                    product_paths,
                    product_identities,
                )
            )
            if not identities_valid:
                return False, (
                    "source-detection published-product provenance is "
                    f"invalid: {identities_reason}"
                )

            for key in [
                "source_catalog_csv",
                "source_catalog_fits",
                "source_regions",
            ]:

                product = products.get(key)
                path = self.resolve_product_path(product)

                if path is None or not path.is_file():

                    return False, (
                        f"source-detection product is missing ({key}): "
                        f"{product}"
                    )

            return True, "source-detection report and catalogues are valid"

        if stage == "ddcal":

            checkpoint = (
                self.stage_outdir / "ddcal" / "checkpoints" / "ddcal.done"
            )
            report, reason = self.read_stage_report("ddcal")

            if not checkpoint.is_file():

                return False, f"missing completion marker: {checkpoint}"

            if report is None:

                return False, reason

            compatible, context_reason = self.validate_report_instrument(
                report,
                "DD calibration",
            )
            if not compatible:
                return False, context_reason
            if report.get("schema_version") != 2:
                return False, "DD report is missing the current provenance schema"
            if (
                not allow_implementation_mismatch and
                not self.report_provenance_is_current(
                    "ddcal",
                    report.get("provenance_fingerprint"),
                )
            ):
                return False, (
                    "DD report provenance fingerprint does not match the "
                    "current MS, catalogue, recipe, and upstream products"
                )
            if (
                self.args.target and
                report.get("target") != self.args.target
            ):
                return False, (
                    f"DD report target {report.get('target')!r} does not "
                    f"match {self.args.target!r}"
                )

            expected_ms = self.resolve_dd_ms()
            reported_ms = self.resolve_product_path(
                report.get("measurement_set")
            )
            if (
                not expected_ms or
                reported_ms is None or
                Path(expected_ms).expanduser().resolve() !=
                reported_ms.expanduser().resolve()
            ):
                return False, (
                    "DD report Measurement Set does not match the current "
                    "DD input"
                )
            try:
                live_ms_identity = IMAGING_SAFETY.path_stat_fingerprint(
                    expected_ms
                )
            except (OSError, TypeError, ValueError) as exc:
                return False, f"could not fingerprint DD input MS: {exc}"
            if report.get("measurement_set_stat_sha256") != live_ms_identity:
                return False, "DD input MS no longer matches report provenance"

            status = report.get("status")
            if report.get("requested_mode") != self.args.dd_mode:
                return False, (
                    "DD mode changed since the recorded run: "
                    f"{report.get('requested_mode')!r} -> "
                    f"{self.args.dd_mode!r}"
                )

            if report.get("requested_engine") != self.args.dd_engine:
                return False, (
                    "DD engine changed since the recorded run: "
                    f"{report.get('requested_engine')!r} -> "
                    f"{self.args.dd_engine!r}"
                )

            products = report.get("products", {})
            source_catalog = self.resolve_product_path(
                products.get("source_catalog_csv")
            )
            if source_catalog is not None and source_catalog.is_file():
                try:
                    live_catalog_identity = (
                        IMAGING_SAFETY.file_content_identity(source_catalog)
                    )
                except (OSError, ValueError) as exc:
                    return False, f"could not fingerprint DD catalogue: {exc}"
                if (
                    report.get("source_catalog_content_identity") !=
                    live_catalog_identity
                ):
                    return False, (
                        "DD source catalogue no longer matches report "
                        "provenance"
                    )

            published_files = report.get("published_files")
            if not isinstance(published_files, list) or not published_files:
                return False, "DD report has no published-file inventory"
            identities_valid, identities_reason = (
                IMAGING_SAFETY.validate_file_content_identities(
                    published_files,
                    report.get("product_content_identities"),
                )
            )
            if not identities_valid:
                return False, (
                    "DD published-product provenance is invalid: "
                    f"{identities_reason}"
                )

            solver_config = self.resolve_product_path(
                products.get("solver_config")
            )
            solver_configuration = (
                read_json(solver_config, {})
                if solver_config is not None and solver_config.is_file()
                else {}
            )
            solver_signature = report.get("solver_signature")
            if (
                not isinstance(solver_signature, str) or
                len(solver_signature) != 64 or
                solver_configuration.get("solver_signature") !=
                solver_signature or
                solver_configuration.get("measurement_set_stat_sha256") !=
                live_ms_identity or
                solver_configuration.get(
                    "source_catalog_content_identity"
                ) != report.get("source_catalog_content_identity")
            ):
                return False, (
                    "DD solver configuration does not match the report's "
                    "MS/catalogue/signature provenance"
                )

            dd_checkpoint_provenance = self.resolve_product_path(
                products.get("ddcal_checkpoint_provenance")
            )
            dd_checkpoint = (
                read_json(dd_checkpoint_provenance, {})
                if (
                    dd_checkpoint_provenance is not None and
                    dd_checkpoint_provenance.is_file()
                )
                else {}
            )
            if (
                dd_checkpoint.get("checkpoint") != "ddcal" or
                dd_checkpoint.get("solver_signature") != solver_signature or
                dd_checkpoint.get("measurement_set_stat_sha256") !=
                live_ms_identity or
                dd_checkpoint.get("source_catalog_content_identity") !=
                report.get("source_catalog_content_identity")
            ):
                return False, (
                    "DD completion checkpoint has missing or mismatched "
                    "input provenance"
                )

            if int(report.get("solver_direction_count", 0) or 0) < 1:
                return False, (
                    "DD report contains no solver calibration directions"
                )

            for key in [
                "directions_csv",
                "directions_region",
                "solver_directions_region",
                "solver_config",
                "command_script",
            ]:
                path = self.resolve_product_path(products.get(key))
                if path is None or not path.is_file():
                    return False, (
                        f"DD preparation product is missing ({key}): "
                        f"{products.get(key)}"
                    )

            if self.args.dd_prepare_only:
                if status not in {"prepared", "success", "skipped"}:
                    return False, (
                        f"DD preparation report status is not complete: "
                        f"{status!r}"
                    )
            else:
                if status != "success":
                    return False, (
                        "DD execution was requested but the report status is "
                        f"{status!r}, not 'success'"
                    )

                solver_checkpoint = (
                    self.stage_outdir /
                    "ddcal" /
                    "checkpoints" /
                    "solver.done"
                )
                if not solver_checkpoint.is_file():
                    return False, (
                        f"missing solver completion marker: "
                        f"{solver_checkpoint}"
                    )
                solver_provenance_path = self.resolve_product_path(
                    products.get("solver_checkpoint_provenance")
                )
                solver_provenance = (
                    read_json(solver_provenance_path, {})
                    if (
                        solver_provenance_path is not None and
                        solver_provenance_path.is_file()
                    )
                    else {}
                )
                if (
                    solver_provenance.get("checkpoint") != "solver" or
                    solver_provenance.get("solver_signature") !=
                    solver_signature or
                    solver_provenance.get(
                        "measurement_set_stat_sha256"
                    ) != live_ms_identity or
                    solver_provenance.get(
                        "source_catalog_content_identity"
                    ) != report.get("source_catalog_content_identity")
                ):
                    return False, (
                        "DD solver checkpoint has missing or mismatched "
                        "input provenance"
                    )

                restored = products.get("dd_restored_fits", [])
                restored_paths = [
                    self.resolve_product_path(item)
                    for item in restored
                ] if isinstance(restored, list) else []
                if not (
                    restored_paths and
                    all(
                        item is not None and item.is_file()
                        for item in restored_paths
                    )
                ):
                    return False, (
                        "DD report has no valid restored FITS products"
                    )

                if not self.args.skip_dd_post_imaging:
                    pbcor = products.get("dd_pbcor_fits", [])
                    pbcor_paths = [
                        self.resolve_product_path(item)
                        for item in pbcor
                    ] if isinstance(pbcor, list) else []
                    if not (
                        pbcor_paths and
                        all(
                            item is not None and item.is_file()
                            for item in pbcor_paths
                        )
                    ):
                        return False, (
                            "DD report has no valid PB-corrected FITS product"
                        )

                    metrics = self.resolve_product_path(
                        products.get("dd_image_metrics")
                    )
                    if metrics is None or not metrics.is_file():
                        return False, (
                            "DD PB-corrected image metrics report is missing"
                        )

                    post_provenance_path = self.resolve_product_path(
                        products.get("post_dd_checkpoint_provenance")
                    )
                    post_provenance = (
                        read_json(post_provenance_path, {})
                        if (
                            post_provenance_path is not None and
                            post_provenance_path.is_file()
                        )
                        else {}
                    )
                    if (
                        post_provenance.get("checkpoint") !=
                        "post_dd_imaging" or
                        post_provenance.get("solver_signature") !=
                        solver_signature
                    ):
                        return False, (
                            "DD post-imaging checkpoint has missing or "
                            "mismatched solver provenance"
                        )

            return True, (
                f"DD-cal checkpoint and {status} report are valid"
            )

        return False, f"no completion validator for stage {stage}"

    # --------------------------------------------------------

    def completed_stage_implementation_is_compatible(
        self,
        stage,
        previous,
        current_recipe_fingerprint,
    ):
        """Prove that a completed stage changed implementation, not recipe."""

        recorded_recipe = previous.get("product_recipe_fingerprint")
        if (
            recorded_recipe and
            recorded_recipe == current_recipe_fingerprint
        ):
            return True, "the recorded science-recipe fingerprint matches"

        if (
            recorded_recipe and
            self.recorded_recipe_matches_after_provenance_refreshes(
                stage,
                recorded_recipe,
            )
        ):
            return True, (
                "the recorded science-recipe fingerprint matches once the "
                "post-flag report's implementation-provenance refreshes are "
                "set aside"
            )

        if (
            recorded_recipe and
            self.recorded_recipe_matches_apart_from_failure_policy(
                stage,
                recorded_recipe,
                previous,
            )
        ):
            return True, (
                "the recorded science-recipe fingerprint matches once the "
                "antenna-quarantine limits are set aside"
            )

        # A stage with a science recipe of its own recorded it apart from
        # its implementation fingerprint. When that record no longer
        # matches, the science changed, whatever the command line says: on
        # the Band 5 run of 2026-10-25 a flux model and a secondary flux
        # calibrator set in config.yaml changed the calibration recipe and
        # left its command as it was, and the old calibration was adopted.
        # The command is compared only for records without such a recipe.
        product_fingerprint = (
            previous.get("product_stage_fingerprint") or
            previous.get("stage_fingerprint")
        )
        if (
            recorded_recipe and
            stage in {"flagging", "calibration", "selfcal", "imaging"} and
            recorded_recipe != product_fingerprint
        ):
            return False, "the recorded science-recipe fingerprint changed"

        old_command = self.normalized_checkpoint_resume_command(
            previous.get("command")
        )
        if not old_command:
            return False, "the completed stage has no recorded command"
        try:
            new_command = self.normalized_checkpoint_resume_command(
                self.command_for_stage(stage)
            )
        except Exception as exc:
            return False, f"the current stage command cannot be resolved: {exc}"
        if old_command != new_command:
            return False, "the recorded science command changed"

        return True, (
            "the normalized science command matches and only embedded "
            "implementation provenance changed"
        )

    # --------------------------------------------------------

    def adopt_completed_stage_after_implementation_change(
        self,
        stage,
        previous,
        current_fingerprint,
        current_recipe_fingerprint,
        compatibility_reason,
        validation_reason,
    ):
        """Atomically migrate a fully validated completed-stage record."""

        history = previous.setdefault("implementation_adoption_history", [])
        history.append({
            "adopted_utc": utc_now(),
            "previous_stage_fingerprint": (
                previous.get("product_stage_fingerprint") or
                previous.get("stage_fingerprint")
            ),
            "current_stage_fingerprint": current_fingerprint,
            "current_recipe_fingerprint": current_recipe_fingerprint,
            "reason": compatibility_reason,
            "validation": validation_reason,
        })
        previous["status"] = "completed"
        previous["product_stage_fingerprint"] = current_fingerprint
        previous["product_recipe_fingerprint"] = (
            current_recipe_fingerprint
        )
        previous["stage_fingerprint"] = current_fingerprint
        previous["validation_stage_fingerprint"] = current_fingerprint
        previous["validation_recipe_fingerprint"] = (
            current_recipe_fingerprint
        )
        previous["validation"] = validation_reason
        previous["validated_utc"] = utc_now()
        previous.pop("validation_error", None)
        self.save_resume_state()

    # --------------------------------------------------------

    def should_skip_stage(self, stage):

        if not self.args.resume:

            return False

        state = self.load_resume_state()
        previous = state["stages"].get(self.stage_state_key(stage), {})
        current_fingerprint = self.stage_fingerprint(stage)
        product_fingerprint = (
            previous.get("product_stage_fingerprint") or
            previous.get("stage_fingerprint")
        )
        attempt_fingerprint = previous.get("attempt_stage_fingerprint")
        current_recipe_fingerprint = self.stage_recipe_fingerprint(stage)
        attempt_recipe_fingerprint = previous.get(
            "attempt_recipe_fingerprint"
        )

        if (
            product_fingerprint and
            product_fingerprint != current_fingerprint
        ):
            compatible, compatibility_reason = (
                self.completed_stage_implementation_is_compatible(
                    stage,
                    previous,
                    current_recipe_fingerprint,
                )
            )
            # A proven-equivalent recorded recipe may now be pinned.
            current_recipe_fingerprint = self.stage_recipe_fingerprint(stage)
            valid = False
            validation_reason = ""
            if compatible:
                valid, validation_reason = self.validate_stage_completion(
                    stage,
                    allow_implementation_mismatch=True,
                )
            if compatible and valid:
                if stage == "inspection":
                    self.refresh_run_manifest_instrument_context()
                self.adopt_completed_stage_after_implementation_change(
                    stage,
                    previous,
                    current_fingerprint,
                    current_recipe_fingerprint,
                    compatibility_reason,
                    validation_reason,
                )
                print(
                    f"Stage '{stage}' was completed by an older pipeline "
                    "implementation; its unchanged science command and all "
                    "durable products were validated and adopted."
                )
                return True

            split_target_repair = (
                stage == "calibration" and
                any(
                    marker in validation_reason.lower()
                    for marker in [
                        "split_target",
                        "split-target",
                        "split targets",
                        "split target",
                        "calibrated target ms",
                    ]
                )
            )
            if compatible and (stage == "inspection" or split_target_repair):
                self.update_stage_state(
                    stage,
                    "stale",
                    validation_error=validation_reason,
                )
                if split_target_repair:
                    self._calibration_split_repair_required = True
                    print(
                        "Calibration solutions and checkpoints are valid, "
                        "but one or more split-target MS products are stale "
                        f"({validation_reason}). The calibration stage will "
                        "preserve each unbound target MS, rebuild only the "
                        "split-target products, refresh its report, and "
                        "continue."
                    )
                else:
                    print(
                        "The old inspection products are incomplete or "
                        f"stale ({validation_reason}). Inspection is "
                        "read-only, so the stage will overwrite its reports "
                        "and continue."
                    )
                return False

            mismatch_reason = (
                validation_reason
                if compatible and validation_reason else
                compatibility_reason
            )
            if stage == "imaging" and self.imaging_backend() == "wsclean":
                # WSClean imaging has no checkpoint that could carry an old
                # recipe into a new one: its only one binds WSClean's outputs
                # to the command, Measurement Set and data column that made
                # them, and its report to the stage fingerprint. Rerunning in
                # place redoes what the change affects, which for a beam
                # correction fix (finding 165) is the correction alone; a new
                # --outdir would repeat flagging and calibration for nothing.
                self.update_stage_state(
                    stage,
                    "stale",
                    validation_error=mismatch_reason,
                )
                print(
                    f"Stage '{stage}' changed recipe since it completed "
                    f"({mismatch_reason}). WSClean imaging keeps no "
                    "checkpoint across recipes, so it reruns in place: "
                    "WSClean runs again only if its command, data or "
                    "outputs changed."
                )
                return False
            raise RuntimeError(
                f"Refusing to rerun stage '{stage}' over its old internal "
                "checkpoints because the completed products cannot be "
                "proven to have the same science recipe and valid live "
                f"artifacts ({mismatch_reason}). Use a new --outdir for a "
                "genuine input/recipe change; code-only upgrades are now "
                "adopted automatically when validation succeeds."
            )

        if (
            stage != "inspection" and
            previous.get("status") in {"running", "failed", "interrupted"} and
            attempt_fingerprint and
            attempt_fingerprint != current_fingerprint
        ):
            compatible = False
            compatibility_reason = ""
            if stage == "calibration" and attempt_recipe_fingerprint:
                compatible = (
                    attempt_recipe_fingerprint == current_recipe_fingerprint
                )
                compatibility_reason = (
                    "the stored science-recipe fingerprint matches"
                )
                if (
                    not compatible and
                    self.recorded_recipe_matches_after_provenance_refreshes(
                        stage,
                        attempt_recipe_fingerprint,
                    )
                ):
                    compatible = True
                    current_recipe_fingerprint = (
                        self.stage_recipe_fingerprint(stage)
                    )
                    compatibility_reason = (
                        "the stored science-recipe fingerprint matches once "
                        "the post-flag report's implementation-provenance "
                        "refreshes are set aside"
                    )
                elif (
                    not compatible and
                    self.recorded_recipe_matches_apart_from_failure_policy(
                        stage,
                        attempt_recipe_fingerprint,
                        previous,
                    )
                ):
                    compatible = True
                    current_recipe_fingerprint = (
                        self.stage_recipe_fingerprint(stage)
                    )
                    compatibility_reason = (
                        "the stored science-recipe fingerprint matches once "
                        "the antenna-quarantine limits, which govern failure "
                        "handling rather than the science, are set aside"
                    )
                elif not compatible:
                    compatibility_reason = (
                        "the stored science-recipe fingerprint changed"
                    )
            elif stage == "calibration":
                compatible, compatibility_reason = (
                    self.legacy_calibration_attempt_is_compatible(previous)
                )
            elif stage == "flagging":
                if (
                    attempt_recipe_fingerprint and
                    attempt_recipe_fingerprint == current_recipe_fingerprint
                ):
                    compatible = True
                    compatibility_reason = (
                        "the stored science-recipe fingerprint matches"
                    )
                else:
                    # Pre-fix flagging states used the implementation hash as
                    # their recipe hash. Re-prove those attempts from the
                    # durable parameter snapshot rather than treating a
                    # code-only upgrade as a science change.
                    compatible, compatibility_reason = (
                        self.legacy_flagging_attempt_is_compatible(previous)
                    )
                    if compatible:
                        compatibility_reason = (
                            "the durable legacy flagging recipe matches"
                        )
            elif stage == "selfcal":
                # The recipe binds the calibrated target MS instance, target,
                # reference antenna, profile, self-calibration configuration
                # and the calibration and split-target reports. Attempts
                # recorded before self-calibration had a recipe of its own
                # stored the implementation fingerprint here and cannot match.
                compatible = bool(
                    attempt_recipe_fingerprint and
                    attempt_recipe_fingerprint == current_recipe_fingerprint
                )
                compatibility_reason = (
                    "the stored science-recipe fingerprint matches"
                    if compatible else
                    "the stored science-recipe fingerprint changed"
                    if attempt_recipe_fingerprint and
                    attempt_recipe_fingerprint != attempt_fingerprint else
                    "the attempt predates a separate self-calibration recipe"
                )
            elif stage == "imaging":
                if (
                    attempt_recipe_fingerprint and
                    attempt_recipe_fingerprint == current_recipe_fingerprint
                ):
                    # The recipe covers the imaging MS instance, target, band
                    # and task configuration, and the upstream reports, so it
                    # settles this on its own.
                    compatible = True
                    compatibility_reason = (
                        "the stored science-recipe fingerprint matches"
                    )
                else:
                    # Imaging states written before a separate recipe
                    # fingerprint stored the implementation hash in both
                    # fields. A strict command/MS/upstream/config audit
                    # permits a code-only PB fix to resume without reusing
                    # products from a changed target or imaging recipe.
                    compatible, compatibility_reason = (
                        self.legacy_imaging_attempt_is_compatible(previous)
                    )

            if not compatible and stage in {"selfcal", "imaging"}:
                # Nothing from the old attempt is reused; see the method.
                self.supersede_failed_attempt(
                    stage,
                    previous,
                    current_fingerprint,
                    current_recipe_fingerprint,
                    compatibility_reason,
                )

            elif not compatible:
                raise RuntimeError(
                    f"Stage '{stage}' was previously attempted with a "
                    "different implementation and its science recipe cannot "
                    "be proven unchanged ("
                    f"{compatibility_reason or 'no matching recipe record'}"
                    "). Use a new --outdir rather than mixing recipes."
                )

            else:
                history = previous.setdefault(
                    "implementation_resume_history",
                    [],
                )
                history.append({
                    "resumed_utc": utc_now(),
                    "previous_stage_fingerprint": attempt_fingerprint,
                    "current_stage_fingerprint": current_fingerprint,
                    "recipe_fingerprint": current_recipe_fingerprint,
                    "reason": compatibility_reason,
                })
                previous["attempt_stage_fingerprint"] = current_fingerprint
                previous["attempt_recipe_fingerprint"] = (
                    current_recipe_fingerprint
                )
                self.save_resume_state()
                print(
                    f"{stage.capitalize()} implementation changed, but its "
                    "Measurement Set, effective science recipe, and upstream "
                    "reports are unchanged. Resuming from individually "
                    "validated internal checkpoints."
                )

        valid, reason = self.validate_stage_completion(stage)

        if valid:

            if stage == "inspection":
                self.refresh_run_manifest_instrument_context()

            self.update_stage_state(
                stage,
                "completed",
                inferred=previous.get("status") != "completed",
                validation=reason,
                validated_utc=utc_now(),
            )

            print(
                f"Stage '{stage}' already completed; validated products "
                "will be reused."
            )

            return True

        if (
            stage == "calibration" and
            any(
                marker in reason.lower()
                for marker in [
                    "split_target",
                    "split-target",
                    "split targets",
                    "split target",
                    "calibrated target ms",
                ]
            )
        ):
            self._calibration_split_repair_required = True

        if previous.get("status") == "completed":

            if (
                stage == "calibration" and
                any(
                    marker in reason.lower()
                    for marker in [
                        "split_target",
                        "split-target",
                        "split targets",
                        "split target",
                        "calibrated target ms",
                    ]
                )
            ):
                self._calibration_split_repair_required = True

            self.update_stage_state(
                stage,
                "stale",
                validation_error=reason,
            )

            print(
                f"Stage '{stage}' completion record is stale: {reason}. "
                "The stage will resume from its internal checkpoints."
            )

        elif previous.get("status") in {"running", "failed", "interrupted"}:

            print(
                f"Resuming stage '{stage}' after previous status "
                f"'{previous.get('status')}'."
            )

        return False

    # --------------------------------------------------------

    def print_resume_status(self, stages):

        state = self.load_resume_state()

        print("\nResume status")
        print("=" * 70)

        for stage in stages:

            entry = state["stages"].get(self.stage_state_key(stage), {})
            valid, reason = self.validate_stage_completion(stage)
            recorded = entry.get("status", "not-recorded")
            reusable = "yes" if valid else "no"

            print(
                f"{stage:12s} recorded={recorded:12s} "
                f"reusable={reusable}  {reason}"
            )

        print("=" * 70)

    # --------------------------------------------------------

    @staticmethod
    def stage_memory_snapshot():
        """Read memory accounting without ever failing the stage."""

        try:
            return MEMORY_DIAGNOSTICS.snapshot()
        except Exception as exc:
            return {"error": str(exc)}

    def casa_log_tail(self, started_unix, lines=200):
        """The end of the CASA log a stage was writing when it died.

        A crash inside CASA leaves its account in that log and nowhere else:
        the task that was running, and the memory figures CASA planned from.
        The orchestrator cannot see any of it, so it collects it here rather
        than asking for it after the fact.
        """

        candidates = []
        for directory in (Path.cwd(), self.outdir, self.stage_outdir):
            try:
                candidates.extend(Path(directory).glob("casa*.log"))
            except OSError:
                continue

        recent = []
        for path in set(candidates):
            try:
                stat = path.stat()
            except OSError:
                continue
            if stat.st_mtime >= started_unix - 60.0:
                recent.append((stat.st_mtime, path))
        if not recent:
            return None

        _, newest = max(recent)
        try:
            with newest.open(errors="replace") as handle:
                tail = handle.readlines()[-int(lines):]
        except OSError:
            return None

        return {"path": str(newest), "lines": [line.rstrip() for line in tail]}

    def write_memory_diagnostics(self, stage, record):
        """Keep the evidence for a killed stage beside the other products."""

        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        path = self.outdir / "diagnostics" / f"{stage}_sigkill_{stamp}.json"
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            atomic_write_json(path, record)
        except OSError as exc:
            print(f"WARNING: could not write {path}: {exc}", flush=True)
            return None
        return path

    def stage_memory_requirement(self, stage, snapshot_value=None):
        """Memory range for a stage whose imaging needs are known in advance.

        Only self-calibration is covered: it always images serially with
        w-projection, from settings the orchestrator can read.
        """

        if stage != "selfcal" or CONFIG is None:
            return None

        band = self.resolved_receiver_band() or self.args.receiver_band
        settings = dict(
            (getattr(CONFIG, "SELFCAL_IMAGING_CONFIG", {}) or {}).get(band)
            or {}
        )
        for key, value in (self.profile_settings.get("selfcal") or {}).items():
            if key in {"imsize", "nterms", "wprojplanes"}:
                settings[key] = value
        if not settings.get("imsize"):
            return None

        budget = self.casa_memory_budget_mb(snapshot_value)
        try:
            return MEMORY_DIAGNOSTICS.wproject_imaging_memory(
                settings["imsize"],
                settings.get("nterms", 2),
                settings.get("wprojplanes", 128),
                self.stage_thread_limit(stage),
                budget_bytes=(budget * 1024 ** 2) if budget else None,
            )
        except (TypeError, ValueError):
            return None

    def check_stage_memory(self, stage, snapshot_value):
        """Stop before a stage that cannot fit, and explain a short allocation."""

        allocation = MEMORY_DIAGNOSTICS.allocation_warning(snapshot_value)
        if allocation:
            print(f"WARNING: {allocation}", flush=True)

        requirement = self.stage_memory_requirement(stage, snapshot_value)
        limit = MEMORY_DIAGNOSTICS.effective_limit_bytes(snapshot_value)
        if not requirement or not limit:
            return

        size = MEMORY_DIAGNOSTICS.format_bytes
        shape = (
            f"{requirement['imsize']} px, nterms={requirement['nterms']}, "
            f"{requirement['wprojplanes']} w-planes, "
            f"{requirement['threads']} threads"
        )
        print(
            f"Self-calibration imaging ({shape}) needs at least "
            f"{size(requirement['minimum_bytes'])}, and up to about "
            f"{size(requirement['largest_bytes'])} depending on how much "
            "memory CASA believes is free.",
            flush=True,
        )

        if limit < requirement["minimum_bytes"]:
            message = (
                f"Stage '{stage}' cannot run in the {size(limit)} this process "
                "may use. Its w-projection imaging needs at least "
                f"{size(requirement['minimum_bytes'])}: CASA builds the "
                "convolution functions on a "
                f"{requirement['convolution_size']}-pixel square grid, and "
                "the image products and gridding buffer come on top."
            )
            if allocation:
                message += " " + allocation
            message += " Nothing was run; request more memory and resubmit."
            if self.args.dry_run:
                print(f"WARNING: {message}", flush=True)
                return
            raise RuntimeError(message)

        planned = requirement.get("planned_bytes")
        if planned:
            print(
                "Told it may use "
                f"{size(requirement['planning_budget_bytes'])}, CASA should "
                f"compute {requirement['planned_planes_per_chunk']} w-planes "
                f"at a time and use about {size(planned)}.",
                flush=True,
            )

        if planned and limit < planned:
            # Both levers are named with their arithmetic, because which one
            # helps depends on the shape: the per-thread screens are charged
            # only while a chunk is that wide, so a smaller image that fits
            # every plane in one chunk can cost more, not less.
            print(
                f"WARNING: this process may use {size(limit)}, less than the "
                f"{size(planned)} self-calibration imaging is expected to "
                "take. Lower `selfcal_imsize` in config.yaml (the "
                f"{requirement['convolution_size']}-pixel convolution grid "
                "scales as its square), or run fewer threads (--threads, "
                f"currently {requirement['threads']}); each thread holds a "
                f"{size(16 * requirement['convolution_size'] ** 2)} screen "
                "while a chunk is being built.",
                flush=True,
            )
        elif not planned and limit < requirement["largest_bytes"]:
            print(
                f"WARNING: this process may use {size(limit)}, less than the "
                f"{size(requirement['largest_bytes'])} self-calibration "
                "imaging can take if CASA plans its w-projection work from a "
                "larger memory figure. Fewer threads (--threads) lower that "
                "figure without changing the image.",
                flush=True,
            )

    def run_per_target_stages(self, stages, target_stages, blocked_upstream):
        """Run the per-target stages for every target this run covers.

        Targets are independent: one that fails does not stop the others,
        because each keeps its own products and resume state. The run still
        fails afterwards, naming every target that did not complete, and
        resubmitting continues from what each of them has.
        """

        scopes = self.science_target_scopes(target_stages)
        if len(scopes) > 1:
            print(
                f"\n{len(scopes)} science targets to process: "
                + ", ".join(str(target) for target, _ in scopes)
            )

        failures = []
        for target, root in scopes:
            if len(scopes) > 1 or target:
                print("\n" + "=" * 70)
                print(f"Target: {target}")
                print(f"Products: {root}")
                print("=" * 70, flush=True)
            try:
                with self.target_scope(target, root):
                    unavailable_stages = set(blocked_upstream)
                    for stage in target_stages:
                        self.run_selected_stage(
                            stage,
                            stages,
                            unavailable_stages,
                        )
            except (RuntimeError, OSError, ValueError) as exc:
                if len(scopes) == 1:
                    raise
                failures.append((target, str(exc)))
                print(f"\nERROR: target {target} did not complete: {exc}")
                print(
                    "Its products and resume state are kept. Continuing "
                    "with the remaining targets.",
                    flush=True,
                )

        if failures:
            raise RuntimeError(
                f"{len(failures)} of {len(scopes)} targets did not complete: "
                + "; ".join(f"{target} ({reason})" for target, reason in failures)
            )

    def run_selected_stage(self, stage, stages, unavailable_stages):
        """Run one selected stage, honouring dependencies and resume state."""

        resolved_telescope = self.resolved_telescope()
        if resolved_telescope:
            TELESCOPE_PROFILE_MODULE.validate_stage_capabilities(
                resolved_telescope,
                self.args.data_state,
                stages,
                imager=self.args.imager,
                receiver_band=(
                    self.resolved_receiver_band() or
                    self.args.receiver_band
                ),
            )

        blocked_by = sorted(
            self.selected_stage_dependencies(stage, stages)
            .intersection(unavailable_stages)
        )
        if blocked_by:
            message = (
                f"Stage '{stage}' is blocked by failed/incomplete "
                f"upstream stage(s): {', '.join(blocked_by)}."
            )
            self.update_stage_state(
                stage,
                "blocked",
                error=message,
                blocked_by=blocked_by,
            )
            unavailable_stages.add(stage)
            print(f"WARNING: {message}")
            return

        if self.should_skip_stage(stage):
            if stage == "inspection":
                self.validate_downstream_target_selection(
                    stages,
                    defer_until_flagging=("flagging" in stages),
                )
            elif stage == "flagging":
                self.validate_downstream_target_selection(stages)
            return

        returncode = self.run_stage(stage)
        if returncode:
            unavailable_stages.add(stage)
        elif stage == "inspection":
            self.validate_downstream_target_selection(
                stages,
                defer_until_flagging=("flagging" in stages),
            )
        elif stage == "flagging":
            self.validate_downstream_target_selection(stages)

    def run_stage(self, stage):

        command = self.command_for_stage(stage)
        self.commands.append(command)

        print("\n")
        print("=" * 70)
        print(f"Stage: {stage}")
        print("=" * 70)
        print(command_to_string(command))

        limit = (
            1
            if stage == "ddcal"
            else self.stage_thread_limit(stage)
        )

        if self.stage_uses_mpi(stage):
            print(
                f"CASA MPI processes: {self.mpi_process_count()} "
                f"({self.mpi_threads_per_process()} numerical threads "
                "per process)."
            )
        else:
            print(
                f"Numerical threads: {limit} "
                "(OMP/BLAS/Accelerate)."
            )

        # Recorded before and during the stage, so that a SIGKILL can be
        # attributed to the job's own memory limit or to something outside it.
        memory_before = self.stage_memory_snapshot()
        if memory_before.get("node") or memory_before.get("job_limit"):
            print(MEMORY_DIAGNOSTICS.summary_line(memory_before), flush=True)
            warning = MEMORY_DIAGNOSTICS.low_node_memory_warning(memory_before)
            if warning:
                print(warning, flush=True)
        self.check_stage_memory(stage, memory_before)

        if self.args.dry_run:

            return 0

        self.update_stage_state(
            stage,
            "running",
            command=[str(part) for part in command],
            returncode=None,
        )

        memory_watch = MEMORY_DIAGNOSTICS.MemoryWatch(memory_before)

        started = time.monotonic()
        started_unix = time.time()

        try:

            process = subprocess.Popen(
                command,
                env=self.stage_environment(stage, memory_before),
            )

        except Exception as exc:

            self.update_stage_state(
                stage,
                "failed",
                returncode=None,
                error=str(exc),
            )
            raise

        interval = self.args.progress_interval

        try:

            if interval > 0:

                while True:

                    try:
                        returncode = process.wait(timeout=interval)
                        break

                    except subprocess.TimeoutExpired:

                        memory_watch.observe(
                            self.stage_memory_snapshot(),
                            time.monotonic() - started,
                        )

                        elapsed_minutes = (
                            time.monotonic() - started
                        ) / 60.0

                        print(
                            f"Stage '{stage}' is still running "
                            f"({elapsed_minutes:.1f} min elapsed).",
                            flush=True,
                        )

            else:

                returncode = process.wait()

        except KeyboardInterrupt:

            self.update_stage_state(
                stage,
                "interrupted",
                returncode=process.poll(),
                elapsed_seconds=time.monotonic() - started,
            )
            raise

        elapsed_seconds = time.monotonic() - started

        if returncode != 0:

            memory_record = None

            if returncode < 0:
                signal_number = -returncode
                try:
                    signal_name = signal.Signals(signal_number).name
                except ValueError:
                    signal_name = f"signal {signal_number}"
                message = (
                    f"Stage '{stage}' was terminated by {signal_name} "
                    f"(exit code {returncode})."
                )
                casa_log = self.casa_log_tail(started_unix)
                if signal_number != signal.SIGKILL:
                    # A crash can also be memory running out, reported by the
                    # allocator rather than by the kernel's killer, so the
                    # same evidence is worth having for any signal.
                    memory_after = self.stage_memory_snapshot()
                    context = MEMORY_DIAGNOSTICS.summary_at_exit(
                        memory_before,
                        memory_after,
                        memory_watch,
                    )
                    if context:
                        message += f" Memory at the time: {context}."
                    diagnostics_path = self.write_memory_diagnostics(
                        stage,
                        {
                            "stage": stage,
                            "returncode": returncode,
                            "signal": signal_name,
                            "elapsed_seconds": elapsed_seconds,
                            "memory_context": context,
                            "casa_log": casa_log,
                            "before": memory_before,
                            "after": memory_after,
                        },
                    )
                    if diagnostics_path:
                        message += f" Evidence: {diagnostics_path}."

                if signal_number == signal.SIGKILL:
                    memory_after = self.stage_memory_snapshot()
                    attribution = MEMORY_DIAGNOSTICS.attribute_kill(
                        memory_before,
                        memory_after,
                        memory_watch,
                    )
                    diagnostics_path = self.write_memory_diagnostics(
                        stage,
                        {
                            "stage": stage,
                            "returncode": returncode,
                            "elapsed_seconds": elapsed_seconds,
                            "attribution": attribution,
                            "casa_log": casa_log,
                            "before": memory_before,
                            "after": memory_after,
                        },
                    )
                    memory_record = {
                        "cause": attribution["cause"],
                        "diagnostics": (
                            str(diagnostics_path)
                            if diagnostics_path else None
                        ),
                    }
                    message += " " + attribution["explanation"]
                    allocation = MEMORY_DIAGNOSTICS.allocation_warning(
                        memory_before
                    )
                    if allocation:
                        message += " " + allocation
                    if diagnostics_path:
                        message += f" Memory evidence: {diagnostics_path}."
                    message += (
                        " SIGKILL is imposed externally, commonly by a "
                        "Slurm/cgroup out-of-memory limit, job cancellation, "
                        "or node failure. The pipeline state is resumable; "
                        "inspect `sacct -j $SLURM_JOB_ID --format="
                        "JobID,State,ExitCode,Elapsed,MaxRSS,ReqMem` and the "
                        "Slurm error log before resubmitting."
                    )

                if casa_log:
                    # Whatever CASA was doing when it died is the one account
                    # of the crash that the traceback cannot give, so it is
                    # printed where the failure is read rather than left for
                    # someone to find on the cluster.
                    print(
                        f"\nLast lines of {casa_log['path']} before the "
                        f"{signal_name}:",
                        flush=True,
                    )
                    for line in casa_log["lines"][-25:]:
                        print(f"  {line}", flush=True)
                    message += f" CASA log: {casa_log['path']}."
            elif returncode == STAGE_INTERRUPTED_EXIT_CODE:
                message = (
                    f"Stage '{stage}' was interrupted from outside the "
                    "pipeline: the job was cancelled, reached its time "
                    "limit, or its session ended. Completed checkpoints are "
                    "kept; resubmitting continues from them."
                )
                self.update_stage_state(
                    stage,
                    "interrupted",
                    returncode=returncode,
                    elapsed_seconds=elapsed_seconds,
                    error=message,
                )
                raise StageInterrupted(message)

            else:
                message = (
                    f"Stage '{stage}' failed with exit code "
                    f"{returncode}."
                )

            state_values = {}
            if memory_record:
                state_values["memory_kill"] = memory_record

            self.update_stage_state(
                stage,
                "failed",
                returncode=returncode,
                elapsed_seconds=elapsed_seconds,
                error=message,
                **state_values,
            )

            if self.args.continue_on_error:

                print(f"WARNING: {message}")
                return returncode

            raise RuntimeError(message)

        if stage in {"ingest", "inspection"}:
            # A new upstream report changes every dependent fingerprint.
            self._stage_fingerprints.clear()
            self._stage_recipe_fingerprints.clear()

        if stage == "inspection":
            # Auto detection becomes available only after the CASA process
            # writes its report. Persist it before recording completion.
            self.refresh_run_manifest_instrument_context()

        valid, reason = self.validate_stage_completion(stage)

        if not valid:

            message = (
                f"Stage '{stage}' exited successfully but its durable "
                f"products are incomplete: {reason}"
            )
            self.update_stage_state(
                stage,
                "failed",
                returncode=returncode,
                elapsed_seconds=elapsed_seconds,
                validation_error=reason,
                error=message,
            )

            if self.args.continue_on_error:

                print(f"WARNING: {message}")
                return 1

            raise RuntimeError(message)

        self.update_stage_state(
            stage,
            "completed",
            returncode=returncode,
            elapsed_seconds=elapsed_seconds,
            validation=reason,
            inferred=False,
        )

        elapsed_minutes = elapsed_seconds / 60.0
        print(
            f"Stage '{stage}' completed "
            f"({elapsed_minutes:.1f} min).",
            flush=True,
        )

        return returncode

    # --------------------------------------------------------

    def print_plan(self, stages):

        print("\n")
        print("=" * 70)
        print("Radio Continuum Pipeline Plan")
        print("=" * 70)
        print(f"Pipeline dir : {self.root}")
        print(f"Output dir   : {self.outdir}")
        if self.input_is_lta():
            print(f"Input LTA    : {self.input_data}")
            print(f"Converters   : {self.args.listscan} ; {self.args.gvfits}")
            print(f"Imported MS  : {self.imported_ms}")
        elif self.input_kind == "fits":
            print(f"Input FITS   : {self.input_data}")
            print(f"FITS format  : {self.args.fits_format}")
            print(f"Imported MS  : {self.imported_ms}")
        elif self.source_ms:
            print(f"Source MS    : {self.source_ms}")
            print(f"Working MS   : {self.working_ms}")
        else:
            print(f"Input MS     : {self.args.ms}")
        print(f"CASA         : {self.args.casa_python or self.args.casa}")
        print(f"Python       : {self.args.python}")
        print(f"Profile      : {self.args.profile}")
        print(
            f"Telescope    : {self.resolved_telescope() or self.args.telescope}"
        )
        print(
            f"Receiver band: {self.resolved_receiver_band() or 'auto'}"
        )
        print(f"Data state   : {self.args.data_state}")
        print(f"Imager       : {self.imaging_backend()}")
        print(f"Resume       : {self.args.resume}")
        print(f"Serial threads: {self.thread_limit()}")
        if self.imaging_backend() == "casa" and self.mpi_enabled():
            print(
                f"CASA MPI      : {self.mpi_process_count()} processes "
                "for final imaging"
            )
        else:
            suffix = ""

            if self.imaging_backend() == "wsclean":
                suffix = " (WSClean backend)"
            elif self.args.mpi == "auto":
                if platform.system() == "Darwin":
                    suffix = " (macOS serial/OpenMP fallback)"
                elif not os.environ.get("SLURM_JOB_ID"):
                    suffix = " (outside Slurm; serial/OpenMP fallback)"
                elif self.mpi_launcher() is None:
                    suffix = " (launcher unavailable; serial fallback)"
                elif self.mpi_rank_capacity() < 2:
                    suffix = " (Slurm exposes one MPI task; serial fallback)"

            print(f"CASA MPI      : off{suffix}")
        print(f"Stages       : {', '.join(stages)}")
        print("=" * 70)

    # --------------------------------------------------------

    def run_manifest_path(self):

        return self.outdir / "pipeline_run_manifest.json"

    # --------------------------------------------------------

    def refresh_run_manifest_instrument_context(self):
        """Persist telescope/band resolved by inspection into provenance."""

        if self.args.dry_run:
            return

        path = self.run_manifest_path()
        if not path.is_file():
            return

        telescope = self.resolved_telescope()
        receiver_band = self.resolved_receiver_band()
        if not telescope:
            return

        manifest = read_json(path, {})
        manifest["instrument_context"] = {
            "telescope": telescope,
            "receiver_band": receiver_band,
            "imaging_backend": self.imaging_backend(),
        }
        manifest["updated_utc"] = utc_now()
        atomic_write_json(path, manifest)

    # --------------------------------------------------------

    def stage_products_exist(self):

        return any(
            (self.outdir / stage).exists()
            for stage in PIPELINE_STAGES
        )

    # --------------------------------------------------------

    def configuration_digest(self):

        if self._configuration_digest is not None:
            return self._configuration_digest

        digest = hashlib.sha256()

        for name in [
            "config.py",
            "config.yaml",
            "telescope_profiles.py",
            "casa_compat.py",
            "fits_ingest.py",
            "inspect.py",
            "flagging.py",
            "flagging_selectors.py",
            "calibration.py",
            "selfcal.py",
            "imaging.py",
            "imaging_safety.py",
            "memory_diagnostics.py",
            "wsclean_imaging.py",
            "meerkat_beam.py",
            "ddcal.py",
            "qa_products.py",
            "pipeline.py",
        ]:

            path = self.root / name

            if path.exists():

                digest.update(path.read_bytes())

        self._configuration_digest = digest.hexdigest()
        return self._configuration_digest

    # --------------------------------------------------------

    def run_identity(self):

        source_ms = self.source_ms

        if source_ms is None and self.args.ms:

            source_ms = Path(self.args.ms).expanduser().resolve()

        source_ms_identity = self.source_ms_identity_for_provenance()

        return {
            "schema_version": 3,
            "input_data": (
                str(self.input_data)
                if self.input_data is not None
                else None
            ),
            "input_kind": self.input_kind,
            "fits_format": (
                self.args.fits_format
                if self.input_kind == "fits"
                else None
            ),
            "imported_ms": (
                str(self.imported_ms)
                if self.imported_ms is not None
                else None
            ),
            "gmrt_flagfiles": list(self.args.gmrt_flagfile),
            "source_ms": (
                str(source_ms)
                if source_ms is not None
                else None
            ),
            "source_ms_identity": source_ms_identity,
            "working_ms": (
                str(self.working_ms)
                if self.working_ms is not None
                else None
            ),
            "configuration_sha256": self.configuration_digest(),
            "execution_profile": self.args.profile,
            "instrument_request": {
                "telescope": self.args.telescope,
                "receiver_band": self.args.receiver_band,
                "data_state": self.args.data_state,
                "imager": self.args.imager,
                "lofar_beam_mode": self.args.lofar_beam_mode,
                "profile_schema_version": (
                    TELESCOPE_PROFILE_MODULE.PROFILE_SCHEMA_VERSION
                ),
            },
            "instrument_context": {
                "telescope": self.resolved_telescope(),
                "receiver_band": self.resolved_receiver_band(),
                "imaging_backend": self.imaging_backend(),
            },
            "thread_limit": self.thread_limit(),
            "mpi_mode": self.args.mpi,
            "mpi_processes": (
                self.mpi_process_count()
                if self.imaging_backend() == "casa" and self.mpi_enabled()
                else None
            ),
            "mpi_threads_per_process": (
                self.mpi_threads_per_process()
                if self.imaging_backend() == "casa" and self.mpi_enabled()
                else None
            ),
            "field_selection": {
                "flux_calibrator": self.args.flux_calibrator,
                "bandpass_calibrator": self.args.bandpass_calibrator,
                "secondary_flux_calibrator": (
                    self.args.secondary_flux_calibrator
                ),
                "phase_calibrator": self.args.phase_calibrator,
                "target_phase_map": mapping_value(
                    self.args.target_phase_map
                ),
                "target_fields": csv_value(self.args.target_fields),
                "reference_antenna": self.args.reference_antenna,
            },
            "downstream_selection": {
                "target": self.args.target,
                "mode": (
                    self._target_selection_origin or
                    (
                        "pending"
                        if self.args.target_selection == "auto"
                        else "explicit_required"
                    )
                ),
                "target_selection": self.args.target_selection,
                "auto_target_max_flagged_fraction": (
                    self.args.auto_target_max_flagged_fraction
                ),
            },
            "calibrator_model": {
                "flux_standard": self.args.flux_standard,
                "flux_model": self.args.flux_model,
                "flux_model_type": self.args.flux_model_type,
                "flux_model_sha256": (
                    sha256_file_or_directory(self.args.flux_model)
                    if self.args.flux_model else None
                ),
            },
        }

    # --------------------------------------------------------

    def validate_outdir_support(self, stages):

        if self.args.dry_run:

            return

        manifest_path = self.run_manifest_path()
        identity = self.run_identity()

        if manifest_path.exists():

            manifest = read_json(manifest_path, {})

            if manifest.get("schema_version", 1) == 1:
                telescope = self.resolved_telescope()
                if telescope not in {None, "ugmrt"}:
                    raise RuntimeError(
                        "Legacy pipeline products have no telescope "
                        "provenance and cannot be adopted for MeerKAT or "
                        "LOFAR. Use a new --outdir."
                    )

                legacy_fields = manifest.setdefault("field_selection", {})
                legacy_fields.setdefault("bandpass_calibrator", None)
                manifest["schema_version"] = 2
                manifest["instrument_request"] = identity[
                    "instrument_request"
                ]
                manifest["instrument_context"] = {
                    "telescope": telescope or "ugmrt",
                    "receiver_band": self.resolved_receiver_band(),
                    "imaging_backend": "casa",
                }
                manifest["calibrator_model"] = identity[
                    "calibrator_model"
                ]
                manifest.setdefault("migration_history", []).append(
                    {
                        "migrated_utc": utc_now(),
                        "from_schema": 1,
                        "assumption": (
                            "legacy uGMRT products; adoption forbidden for "
                            "MeerKAT/LOFAR"
                        ),
                    }
                )
                atomic_write_json(manifest_path, manifest)

            if manifest.get("schema_version") == 2:
                raise RuntimeError(
                    "This output directory predates immutable source-MS "
                    "identity tracking and cannot be adopted safely. Use a "
                    "new --outdir so products cannot be relabelled after an "
                    "input MS is replaced at the same pathname."
                )

            # Schema-2 manifests created before visibility-FITS support are
            # unambiguously Measurement Set inputs.
            manifest.setdefault("input_data", manifest.get("source_ms"))
            manifest.setdefault("input_kind", "ms")
            manifest.setdefault("fits_format", None)
            manifest.setdefault("imported_ms", None)
            manifest.setdefault("gmrt_flagfiles", [])

            # A historical development manifest briefly stored --target in
            # field_selection.  That option selects one already-calibrated
            # split target and must not invalidate shared inspection,
            # flagging, or calibration products.
            manifest.setdefault("field_selection", {}).pop("target", None)

            for key in [
                "schema_version",
                "field_selection",
            ]:

                if manifest.get(key) != identity.get(key):

                    raise RuntimeError(
                        "Output directory belongs to a different MS, "
                        "field selection, or pipeline configuration:\n"
                        f"{self.outdir}\n"
                        "Use a new --outdir for this run."
                    )

            # A stage-scoped downstream restart (for example
            # ``--from-stage imaging``) does not consume the raw source MS.
            # With no positional input, leave its immutable identity anchored
            # to the existing run manifest and validate the actual downstream
            # MS through the selfcal/imaging report chain.  When a source is
            # explicitly supplied it must still match exactly.
            if (
                identity.get("source_ms_identity") is not None and
                manifest.get("source_ms_identity") !=
                identity.get("source_ms_identity")
            ):
                raise RuntimeError(
                    "The explicitly supplied source Measurement Set does "
                    "not match the immutable source identity recorded for "
                    f"this output directory:\n{self.outdir}\n"
                    "Use the original source MS, or omit the positional "
                    "input for a report-driven downstream-only restart."
                )

            old_request = manifest.get("instrument_request", {})
            new_request = identity.get("instrument_request", {})

            if identity.get("input_data") is not None:
                for key in [
                    "input_data",
                    "input_kind",
                    "fits_format",
                    "imported_ms",
                    "gmrt_flagfiles",
                ]:
                    if manifest.get(key) != identity.get(key):
                        raise RuntimeError(
                            "The output directory belongs to a different "
                            "input data product or FITS import destination. "
                            "Use a new --outdir."
                        )
            if (
                old_request.get("data_state") and
                old_request.get("data_state") != new_request.get("data_state")
            ):
                raise RuntimeError(
                    "The output directory was created for a different input "
                    "calibration state. Use a new --outdir."
                )
            old_telescope = old_request.get("telescope")
            new_telescope = new_request.get("telescope")
            if (
                old_telescope not in (None, "", "auto") and
                new_telescope not in (None, "", "auto") and
                old_telescope != new_telescope
            ):
                raise RuntimeError(
                    "The output directory was created for a different "
                    "explicit telescope. Use a new --outdir."
                )
            old_band = old_request.get("receiver_band")
            new_band = new_request.get("receiver_band")
            if old_band and new_band and str(old_band).upper() != str(new_band).upper():
                raise RuntimeError(
                    "The output directory was created for a different "
                    "explicit receiver band. Use a new --outdir."
                )

            # Later-stage restarts may not receive the raw MS again.  When it
            # is supplied, however, it must agree with the original run.
            for key in ["source_ms", "working_ms"]:

                if (
                    identity.get(key) is not None and
                    manifest.get(key) != identity.get(key)
                ):

                    raise RuntimeError(
                        "Output directory belongs to a different MS, "
                        "field selection, or pipeline configuration:\n"
                        f"{self.outdir}\n"
                        "Use a new --outdir for this run."
                    )

            configuration_changed = (
                manifest.get("configuration_sha256") !=
                identity.get("configuration_sha256")
            )
            thread_changed = (
                manifest.get("thread_limit") !=
                identity.get("thread_limit")
            )
            mpi_changed = (
                manifest.get("mpi_mode") != identity.get("mpi_mode") or
                manifest.get("mpi_processes") !=
                identity.get("mpi_processes") or
                manifest.get("mpi_threads_per_process") !=
                identity.get("mpi_threads_per_process")
            )
            profile_changed = (
                manifest.get("execution_profile") !=
                identity.get("execution_profile")
            )

            if configuration_changed and not self.args.resume:

                raise RuntimeError(
                    "Pipeline code/configuration changed since this output "
                    "directory was created. Use --resume to retain validated "
                    "completed stages, or use a new --outdir for a clean run."
                )

            if configuration_changed:

                history = manifest.setdefault(
                    "configuration_history",
                    [],
                )
                previous_digest = manifest.get("configuration_sha256")

                if previous_digest and not any(
                    item.get("sha256") == previous_digest
                    for item in history
                    if isinstance(item, dict)
                ):

                    history.append({
                        "sha256": previous_digest,
                        "recorded_utc": manifest.get("created_utc"),
                    })

                history.append({
                    "sha256": identity["configuration_sha256"],
                    "recorded_utc": utc_now(),
                    "reason": "resume after pipeline code/configuration change",
                })

                print(
                    "WARNING: pipeline code/configuration changed. Resume "
                    "will validate every durable product. Completed stages "
                    "with an unchanged normalized science command are "
                    "adopted under the new implementation fingerprint; "
                    "genuine input or recipe changes still require a new "
                    "output directory."
                )

            if thread_changed:

                print(
                    "Numerical thread limit changed for this resume: "
                    f"{manifest.get('thread_limit')} -> "
                    f"{identity.get('thread_limit')}"
                )

            if mpi_changed:
                print(
                    "CASA MPI execution changed for this resume: "
                    f"{manifest.get('mpi_processes')} -> "
                    f"{identity.get('mpi_processes')} processes "
                    f"x {identity.get('mpi_threads_per_process')} threads "
                    f"(mode={identity.get('mpi_mode')})"
                )

            downstream_changed = (
                manifest.get("downstream_selection") !=
                identity.get("downstream_selection")
            )
            request_changed = old_request != new_request
            calibrator_model_changed = (
                manifest.get("calibrator_model") !=
                identity.get("calibrator_model")
            )

            if (
                configuration_changed or
                profile_changed or
                thread_changed or
                mpi_changed or
                downstream_changed or
                request_changed or
                calibrator_model_changed
            ):

                manifest["configuration_sha256"] = identity[
                    "configuration_sha256"
                ]
                manifest["execution_profile"] = identity[
                    "execution_profile"
                ]
                manifest["thread_limit"] = identity["thread_limit"]
                manifest["mpi_mode"] = identity["mpi_mode"]
                manifest["mpi_processes"] = identity["mpi_processes"]
                manifest["mpi_threads_per_process"] = identity[
                    "mpi_threads_per_process"
                ]
                manifest["downstream_selection"] = identity["downstream_selection"]
                manifest["instrument_request"] = new_request
                manifest["calibrator_model"] = identity["calibrator_model"]
                manifest["updated_utc"] = utc_now()
                atomic_write_json(manifest_path, manifest)

            return

        if self.stage_products_exist():

            inspection, reason = self.read_stage_report("inspection")
            reported_ms = (
                inspection.get("measurement_set")
                if inspection is not None
                else None
            )
            reported_path = self.resolve_product_path(reported_ms)
            allowed_paths = {
                Path(value).expanduser().resolve()
                for value in [
                    identity.get("source_ms"),
                    identity.get("working_ms"),
                ]
                if value
            }

            if (
                self.args.resume and
                reported_path is not None and
                reported_path.resolve() in allowed_paths
            ):

                identity["created_utc"] = utc_now()
                identity["adopted_existing_products"] = True
                atomic_write_json(manifest_path, identity)
                print(
                    "Adopted existing pipeline products after matching the "
                    "inspection-report MS to this run."
                )
                return

            raise RuntimeError(
                "Output directory contains untracked pipeline products:\n"
                f"{self.outdir}\n"
                f"Could not safely match them to this MS ({reason}). "
                "Use a new --outdir or restore its run manifest."
            )

        if identity["source_ms"] is None:

            raise RuntimeError(
                "Cannot create a new pipeline run without an input "
                "Measurement Set. Run the raw-data stages first, or use "
                "an existing output directory with its run manifest."
            )

        identity["created_utc"] = utc_now()
        atomic_write_json(manifest_path, identity)

    # --------------------------------------------------------

    def required_source_versions(self, stages):
        """The pipeline files this run depends on, with minimum versions."""

        required = dict(SOURCE_INTERFACE_VERSIONS)
        if "ingest" in stages and self.input_is_lta():
            required.update(LTA_SOURCE_INTERFACE_VERSIONS)
        telescope = self.resolved_telescope()
        # Before inspection the telescope may be unknown, in which case
        # either telescope that defaults to WSClean could turn up.
        wsclean_possible = "imaging" in stages and (
            self.args.imager == "wsclean" or (
                self.args.imager == "auto" and
                telescope in (None, "meerkat", "lofar")
            )
        )
        if wsclean_possible:
            required.update(WSCLEAN_SOURCE_INTERFACE_VERSIONS)
            if telescope in (None, "meerkat"):
                required.update(MEERKAT_WSCLEAN_SOURCE_INTERFACE_VERSIONS)
        return required

    def check_lta_tools(self):
        """Stop at the start if an LTA input cannot be converted.

        listscan and gvfits are NCRA programs, not part of CASA. The ingest
        stage runs them before anything else, so a missing one would stop
        the run at once anyway; checking here says what to do about it
        before a job is queued. A dry run only warns.
        """

        telescope = self.args.telescope
        if (
            telescope not in (None, "", "auto") and
            TELESCOPE_PROFILE_MODULE.normalize_telescope(telescope) != "ugmrt"
        ):
            raise RuntimeError(
                f"An LTA file is GMRT data; --telescope {telescope} cannot "
                "be used with it."
            )

        problems = []
        for name in ("listscan", "gvfits"):
            command = getattr(self.args, name, None) or name
            try:
                parts = shlex.split(str(command))
            except ValueError as exc:
                problems.append(f"{name} ({command}): {exc}")
                continue
            if not parts or shutil.which(parts[0]) is None:
                problems.append(
                    f"{name} ({command}) was not found or is not executable"
                )
        if not problems:
            return
        message = (
            "The LTA input cannot be converted: " + "; ".join(problems) +
            ". Linux builds of listscan and gvfits are in "
            "https://github.com/ruta-k/CAPTURE-CASA6. Download them, make "
            "them executable (chmod +x), and give their paths with "
            "--listscan and --gvfits (gmrt_listscan and gmrt_gvfits in "
            "config.yaml)."
        )
        if self.args.dry_run:
            print("WARNING: " + message)
            return
        raise FileNotFoundError(message)

    def check_wsclean_tools(self, stages):
        """Stop at the start if WSClean, or MeerKAT's beam, is unavailable.

        Imaging is the last stage. Discovering that WSClean is not installed,
        or that the Python running it has no katbeam, after calibration and
        self-calibration have used hours of a job is the worst time to learn
        it, so both are checked now. Nothing here changes the data.
        """

        if (
            "imaging" not in stages or
            self.imaging_backend() != "wsclean" or
            self.args.dry_run
        ):
            return

        wsclean = self.args.wsclean
        if shutil.which(wsclean) is None and not Path(wsclean).is_file():
            if self.resolved_telescope() == "lofar":
                # LOFAR has no CASA imaging path to fall back on.
                advice = (
                    "LOFAR imaging requires WSClean. Load or install it, or "
                    "pass --wsclean /path/to/wsclean."
                )
            else:
                advice = (
                    f"{self.resolved_telescope() or 'This'} imaging uses "
                    "WSClean by default. Load or install it (pass --wsclean "
                    "/path/to/wsclean), or image with --imager casa to use "
                    "CASA tclean instead."
                )
            raise FileNotFoundError(
                f"WSClean executable not found: {wsclean}. {advice}"
            )

        if self.resolved_telescope() != "meerkat":
            return
        try:
            completed = subprocess.run(
                [
                    self.args.python,
                    "-c",
                    "import katbeam, astropy.io.fits",
                ],
                capture_output=True,
                text=True,
                timeout=120,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise RuntimeError(
                f"Could not run {self.args.python} to check MeerKAT's "
                f"primary-beam dependencies: {exc}"
            ) from exc
        if completed.returncode != 0:
            detail = (completed.stderr or "").strip().splitlines()
            raise RuntimeError(
                "MeerKAT WSClean imaging applies the katbeam primary beam "
                f"and needs katbeam and astropy in {self.args.python} "
                "(--python): "
                + (detail[-1] if detail else "the import failed")
                + ". Install them there (pip install katbeam astropy), or "
                "image with --imager casa."
            )

    def preflight(self, stages):

        outdated = outdated_source_files(
            self.root, self.required_source_versions(stages)
        )
        if outdated:
            raise RuntimeError(
                "This pipeline.py needs newer copies of other pipeline files "
                f"in {self.root}: " + "; ".join(outdated) + ". Copy the "
                "matching versions of those files, then resubmit. Nothing was "
                "run."
            )

        explicit_telescope = (
            self.args.telescope
            if self.args.telescope not in (None, "", "auto")
            else None
        )
        if explicit_telescope:
            TELESCOPE_PROFILE_MODULE.validate_stage_capabilities(
                explicit_telescope,
                self.args.data_state,
                stages,
                imager=self.args.imager,
                receiver_band=self.args.receiver_band,
            )

        if self.args.gmrt_flagfile:
            if self.input_kind != "fits":
                raise RuntimeError(
                    "--gmrt-flagfile is valid only with visibility-FITS "
                    "ingestion, not an existing Measurement Set."
                )
            if explicit_telescope not in (None, "ugmrt"):
                raise RuntimeError(
                    "--gmrt-flagfile requires --telescope ugmrt (or "
                    "auto-detected GMRT/uGMRT FITS metadata)."
                )
            if self.args.fits_format == "fitsidi":
                raise RuntimeError(
                    "--gmrt-flagfile is supported only by uGMRT UVFITS "
                    "ingestion, not FITS-IDI."
                )
            if not self.args.dry_run:
                missing = [
                    path for path in self.args.gmrt_flagfile
                    if not Path(path).is_file()
                ]
                if missing:
                    raise FileNotFoundError(
                        "GMRT online flag file not found: " + ", ".join(missing)
                    )

        preflight_telescope = (
            TELESCOPE_PROFILE_MODULE.normalize_telescope(explicit_telescope)
            if explicit_telescope
            else self.resolved_telescope()
        )
        if preflight_telescope == "lofar" and self.input_kind == "fits":
            raise RuntimeError(
                "LOFAR visibility FITS is not a validated calibrated-input "
                "handoff. Run the raw observation through LINC/DP3/LoSoTo "
                "and provide its finalized, CASA-readable target "
                "Measurement Set; this pipeline does not infer LINC "
                "calibration or beam provenance from a FITS container."
            )
        if (
            preflight_telescope == "lofar" and
            "imaging" in stages and
            not self.args.lofar_beam_mode
        ):
            raise RuntimeError(
                "LOFAR imaging requires --lofar-beam-mode full, "
                "differential, or none, including when validated products "
                "might otherwise be reused."
            )
        if (
            preflight_telescope == "lofar" and
            self.args.lofar_beam_mode == "none" and
            "source_detection" in stages and
            str(
                SOURCE_DETECTION_TASK_CONFIG.get(
                    "image_preference",
                    "pbcor",
                )
            ).lower() == "pbcor"
        ):
            raise RuntimeError(
                "LOFAR --lofar-beam-mode none produces only an uncorrected "
                "quick-look image, but source detection is configured to "
                "require a PB-corrected image. Omit source_detection from "
                "--stages, disable it in config.yaml, or explicitly set "
                "source_detection_image_preference: restored for a "
                "non-publication quick-look catalogue."
            )

        raw_ms_stages = self.raw_ms_stages()

        if raw_ms_stages.intersection(stages) or "ingest" in stages:

            if not self.args.ms:

                raise RuntimeError(
                    "No input Measurement Set or visibility FITS is "
                    "configured. "
                    "Set PIPELINE_RUN_CONFIG['ms'] in config.py "
                    "or pass the input path on the command line."
                )

            if (
                self.input_kind == "fits" and
                self._prepared_working_from_manifest and
                "ingest" not in stages
            ):
                resumed_ms = Path(self.args.ms).expanduser().resolve()
                if not self.casa_table_exists(resumed_ms):
                    raise FileNotFoundError(
                        "The imported Measurement Set recorded by the run "
                        f"manifest is unavailable:\n{resumed_ms}"
                    )
            elif self.input_kind == "fits":
                if self.input_data is None:
                    raise RuntimeError("No visibility FITS input is available.")
                if not self.args.dry_run and not self.input_data.is_file():
                    raise FileNotFoundError(
                        (
                            "Input GMRT LTA file not found:\n"
                            if self.input_is_lta()
                            else "Input visibility FITS file not found:\n"
                        ) +
                        f"{self.input_data}"
                    )
                # A dry run still validates an existing file's headers.  This
                # prevents an ordinary restored-image FITS from ever being
                # planned as calibratable visibility data, while retaining
                # the useful ability to plan commands for a not-yet-mounted
                # input path. An LTA file has no FITS header: only listscan
                # reads it, so the converters are checked instead.
                if self.input_is_lta():
                    self.check_lta_tools()
                elif self.input_data.is_file():
                    try:
                        FITS_INGEST.detect_visibility_fits(
                            self.input_data,
                            self.args.fits_format,
                        )
                    except (ValueError, OSError) as exc:
                        raise ValueError(
                            "Input FITS is not supported interferometric "
                            f"visibility data: {exc}"
                        ) from exc
                if self.imported_ms == self.input_data:
                    raise RuntimeError(
                        "The imported Measurement Set path must differ from "
                        f"the FITS input: {self.input_data}"
                    )
                if self.imported_ms.is_symlink():
                    raise RuntimeError(
                        "The imported Measurement Set destination cannot be "
                        f"a symlink: {self.imported_ms}"
                    )
                ingest_products = (self.outdir / "ingest").resolve()
                if (
                    self.imported_ms == ingest_products or
                    self.is_within(ingest_products, self.imported_ms)
                ):
                    raise RuntimeError(
                        "The imported Measurement Set cannot contain the "
                        "ingest report/checkpoint directory. Choose a path "
                        "below that directory (the default) or elsewhere: "
                        f"{self.imported_ms}"
                    )
            else:
                input_ms = Path(self.args.ms)
                if not self.args.dry_run and not input_ms.exists():
                    raise FileNotFoundError(
                        "Input Measurement Set not found:\n"
                        f"{input_ms}"
                    )
                if not self.args.dry_run and not input_ms.is_dir():
                    raise RuntimeError(
                        "Input Measurement Set must be a directory:\n"
                        f"{input_ms}"
                    )
                if (
                    not self.args.dry_run and
                    not self.casa_table_exists(input_ms)
                ):
                    raise RuntimeError(
                        "Input Measurement Set is not a readable CASA table "
                        "(table.dat is missing):\n"
                        f"{input_ms}"
                    )

        if (
            self.args.working_ms and
            self.input_kind == "fits"
        ):
            raise RuntimeError(
                "--working-ms is not used for FITS input; choose the "
                "conversion destination with --imported-ms."
            )

        if (
            self.input_kind == "ms" and
            self.mutable_ms_stages().intersection(stages) and
            not self.args.copy_input_ms
        ):
            raise RuntimeError(
                "Raw-MS flagging/calibration cannot run with "
                "--no-copy-input-ms. Those stages mutate flags and "
                "CORRECTED_DATA, which would both alter the archival input "
                "and invalidate its resume identity. Use the default "
                "managed working copy (or --working-ms PATH)."
            )

        if (
            self.args.working_ms and
            not self.args.copy_input_ms
        ):
            raise RuntimeError(
                "--working-ms requires --copy-input-ms."
            )

        casa_stages = {
            "ingest",
            "inspection",
            "flagging",
            "calibration",
            "selfcal",
        }
        if self.imaging_backend() == "casa":
            casa_stages.add("imaging")

        if self.casa_python_command() and self.args.mpi == "on":
            raise RuntimeError(
                "--mpi on needs mpicasa, which a modular CASA (casa_python) "
                "does not have. Use --mpi off or auto."
            )

        if (
            casa_stages.intersection(stages) and
            not self.args.dry_run and
            self.casa_python_command()
        ):
            # Start CASA's Python once now: a wrong image or a Python without
            # casatasks stops the run before any stage, not hours in.
            version = casa_python_version(self.casa_python_command())
            print(f"CASA {version} via {self.args.casa_python}")

        elif (
            casa_stages.intersection(stages) and
            not self.args.dry_run
        ):

            casa_path = shutil.which(self.args.casa)

            if casa_path is None and not Path(self.args.casa).exists():

                module = getattr(self.args, "casa_module", None)
                raise FileNotFoundError(
                    "CASA executable not found:\n"
                    f"{self.args.casa}\n" + (
                        f"Loading the module {module} did not provide it; "
                        "check the name with 'module avail casa'."
                        if module else
                        "Set casa_module in config.yaml (for example "
                        "casa/6.7.0; 'module avail casa' lists them), or "
                        "casa_python for a modular CASA (for example "
                        "'singularity exec IMAGE.sif python3'), load CASA "
                        "into PATH, or pass --casa /path/to/casa."
                    )
                )

        python_stages = {
            "source_detection",
            "ddcal",
        }
        if self.imaging_backend() == "wsclean":
            python_stages.add("imaging")

        if (
            python_stages.intersection(stages) and
            not self.args.dry_run
        ):

            python_path = shutil.which(self.args.python)

            if python_path is None and not Path(self.args.python).exists():

                raise FileNotFoundError(
                    "Python executable for non-CASA stages not found:\n"
                    f"{self.args.python}\n"
                    "Set PIPELINE_RUN_CONFIG['python'] or pass --python."
                )

        self.check_wsclean_tools(stages)

        if (
            "ddcal" in stages and
            not self.args.dd_prepare_only and
            self.args.dd_engine in {"auto", "killms-ddfacet"} and
            not self.args.dry_run
        ):
            missing = [
                name
                for name in ["DDF.py", "MakeModel.py", "kMS.py"]
                if shutil.which(name) is None
            ]
            if missing:
                raise FileNotFoundError(
                    "DD execution was requested, but required killMS/"
                    "DDFacet commands are not in PATH: " +
                    ", ".join(missing) +
                    ". Load or activate the solver environment before "
                    "starting the pipeline (or in slurm_module_commands)."
                )

        if (
            "ddcal" in stages and
            not self.args.dd_prepare_only and
            not self.args.skip_dd_post_imaging and
            not self.args.dry_run
        ):
            casa_path = (
                "casa_python" if self.casa_python_command()
                else shutil.which(self.args.casa)
            )
            if casa_path is None and not Path(self.args.casa).exists():
                raise FileNotFoundError(
                    "CASA is required for final DD primary-beam correction. "
                    "Load CASA, pass --casa, or use "
                    "--skip-dd-post-imaging."
                )

        if self.args.mpi == "on" and platform.system() == "Darwin":
            raise RuntimeError(
                "CASA mpicasa is not supported on macOS. Use --mpi off "
                "for local Mac runs."
            )

        if self.args.mpi == "on" and self.mpi_launcher() is None:
            raise FileNotFoundError(
                "CASA MPI was forced with --mpi on, but mpicasa was not "
                "found. Load the matching CASA installation or pass "
                "--mpicasa /path/to/mpicasa."
            )

        if self.args.mpi == "on" or self.mpi_enabled():
            self.mpi_process_count()

        if not self.args.dry_run:

            self.outdir.mkdir(
                parents=True,
                exist_ok=True,
            )

    # --------------------------------------------------------

    def print_final_products(self, stages):

        if not {"imaging", "ddcal"}.intersection(stages):

            return

        scopes = self.science_target_scopes(stages) or [(None, self.outdir)]
        if len(scopes) > 1 or self.target_scope_name is None:
            for target, root in scopes:
                with self.target_scope(target, root):
                    self.print_target_final_products(stages)
            return

        self.print_target_final_products(stages)

    def print_target_final_products(self, stages):

        report = read_json(
            self.report_path("imaging"),
            {},
        )

        products = report.get(
            "products",
            {},
        )

        pbcor = products.get("pbcor_image")
        subbands = products.get("subbands", [])
        fits = products.get("fits", [])

        print("\n")
        print("=" * 70)
        title = "Final PB-Corrected Products"
        if self.target_scope_name:
            title += f" -- {self.target_scope_name}"
        print(title)
        print("=" * 70)

        if pbcor:

            print(f"Total intensity : {pbcor}")

        else:

            print(
                "Total intensity : not found in imaging report"
            )

        pbcor_subbands = [
            item.get("pbcor_image")
            for item in subbands
            if item.get("pbcor_image")
        ]

        if pbcor_subbands:

            print(
                f"In-band images  : {len(pbcor_subbands)} "
                "PB-corrected sub-band images"
            )

            for image in pbcor_subbands:

                print(f"                 {image}")

        if fits:

            print(
                f"FITS exports    : {len(fits)} files in "
                f"{self.outdir / 'imaging' / 'fits'}"
            )

        print("=" * 70)

        if "ddcal" in stages:
            dd_report = read_json(
                self.report_path("ddcal"),
                {},
            )
            dd_products = dd_report.get("products", {})

            print("\n")
            print("=" * 70)
            print("Direction-Dependent Calibration Products")
            print("=" * 70)
            print(
                f"Status / mode   : "
                f"{dd_report.get('status', 'unknown')} / "
                f"{dd_report.get('mode', 'unknown')}"
            )

            for label, key in [
                ("DD restored FITS", "dd_restored_fits"),
                ("DD PB-corrected", "dd_pbcor_fits"),
            ]:
                values = dd_products.get(key, [])
                if isinstance(values, str):
                    values = [values]
                if values:
                    for index, value in enumerate(values):
                        prefix = label if index == 0 else ""
                        print(f"{prefix:17s}: {value}")
                else:
                    print(f"{label:17s}: not produced")

            print(
                f"Gain QA         : "
                f"{dd_products.get('solution_qa', 'not produced')}"
            )
            print("=" * 70)

    # --------------------------------------------------------

    @staticmethod
    def selected_stage_dependencies(stage, selected_stages):
        """Return selected upstream stages whose failure blocks *stage*.

        ``--continue-on-error`` is diagnostic convenience, not permission to
        consume stale products after a current upstream failure.  Dependencies
        are conditional so a deliberately supplied downstream MS can still be
        processed when raw-data stages were not selected at all.
        """

        selected = set(selected_stages or ())

        if stage == "inspection":
            return {"ingest"}.intersection(selected)
        if stage == "flagging":
            return {"inspection"}.intersection(selected)
        if stage == "calibration":
            return {"inspection", "flagging"}.intersection(selected)
        if stage == "selfcal":
            preferred = "calibration" if "calibration" in selected else "inspection"
            return {preferred}.intersection(selected)
        if stage == "imaging":
            for preferred in ("selfcal", "calibration", "inspection"):
                if preferred in selected:
                    return {preferred}
            return set()
        if stage == "source_detection":
            return {"imaging"}.intersection(selected)
        if stage == "ddcal":
            for preferred in ("source_detection", "imaging"):
                if preferred in selected:
                    return {preferred}
            return set()
        return set()

    # --------------------------------------------------------

    def rehydrate_resume_status_context(self):
        """Restore the original read-only run context from its manifest.

        Status mode runs before the normal working-copy preparation on
        purpose.  Reconstructing both the immutable source and the effective
        mutable MS here makes stage selection and fingerprints identical to a
        normal resume without copying, locking, or writing any product.
        """

        manifest_path = self.run_manifest_path()
        if not manifest_path.is_file():
            if self.input_data is None:
                raise RuntimeError(
                    "--resume-status needs the original input path or an "
                    "existing pipeline_run_manifest.json."
                )
            return

        try:
            manifest = read_json(manifest_path, {})
        except (OSError, ValueError, TypeError) as exc:
            raise RuntimeError(
                "Cannot reconstruct --resume-status because the run "
                f"manifest is unreadable: {exc}"
            ) from exc
        if not isinstance(manifest, dict):
            raise RuntimeError(
                "Cannot reconstruct --resume-status from a non-object run "
                "manifest."
            )

        manifest_source_value = (
            manifest.get("input_data") or manifest.get("source_ms")
        )
        if not manifest_source_value:
            raise RuntimeError(
                "The run manifest does not record the immutable source "
                "input required by --resume-status."
            )
        manifest_source = Path(
            manifest_source_value
        ).expanduser().resolve()
        manifest_kind = str(
            manifest.get("input_kind") or
            self.classify_input_path(
                manifest_source,
                manifest.get("fits_format", "auto"),
            )
        ).strip().lower()
        if manifest_kind not in {"ms", "fits"}:
            raise RuntimeError(
                "The run manifest contains an unsupported input kind: "
                f"{manifest_kind!r}."
            )

        explicit_source = (
            None if self._input_from_config else self.input_data
        )
        no_explicit_input = explicit_source is None
        if (
            explicit_source is not None and
            explicit_source.resolve() != manifest_source
        ):
            raise RuntimeError(
                "The explicit --resume-status input does not match the "
                f"manifest source: {manifest_source}"
            )
        if (
            explicit_source is not None and
            self.input_kind and
            self.input_kind != manifest_kind
        ):
            raise RuntimeError(
                "The explicit --resume-status input kind does not match "
                "the run manifest."
            )

        self.input_data = manifest_source
        self.input_kind = manifest_kind

        if manifest_kind == "ms":
            self.source_ms = Path(
                manifest.get("source_ms") or manifest_source
            ).expanduser().resolve()
            recorded_source_identity = manifest.get("source_ms_identity")
            self._recorded_source_ms_identity = recorded_source_identity
            if explicit_source is not None:
                try:
                    live_source_identity = self.measurement_set_identity(
                        manifest_source
                    )
                except (OSError, ValueError) as exc:
                    raise RuntimeError(
                        "The manifest source Measurement Set cannot be "
                        f"validated for --resume-status: {exc}"
                    ) from exc
                if (
                    recorded_source_identity is not None and
                    recorded_source_identity != live_source_identity
                ):
                    raise RuntimeError(
                        "The source Measurement Set no longer matches the "
                        "immutable identity recorded in the run manifest."
                    )
            working_value = manifest.get("working_ms")
            if working_value:
                self.working_ms = Path(
                    working_value
                ).expanduser().resolve()
                self.args.ms = str(self.working_ms)
            else:
                self.args.ms = str(manifest_source)
        else:
            self.args.fits_format = self.normalize_fits_format(
                manifest.get("fits_format") or self.args.fits_format
            )
            imported_value = manifest.get("imported_ms")
            if not imported_value:
                raise RuntimeError(
                    "The FITS run manifest does not record its imported "
                    "Measurement Set."
                )
            self.imported_ms = Path(
                imported_value
            ).expanduser().resolve()
            self.args.imported_ms = str(self.imported_ms)
            self.args.ms = str(self.imported_ms)

        # With no positional input, restore every omitted manifest field that
        # affects stage selection or provenance. Explicit CLI options remain
        # authoritative and produce a mismatch if intentionally changed.
        if no_explicit_input:
            self.restore_manifest_science_context(manifest)

            if not self.args.stages:
                state = read_json(self.state_path(), {})
                recorded_stages = (
                    state.get("stages", {})
                    if isinstance(state, dict)
                    else {}
                )
                if isinstance(recorded_stages, dict) and recorded_stages:
                    ordered = [
                        stage for stage in PIPELINE_STAGES
                        if stage in recorded_stages
                    ]
                    if ordered:
                        self.args.stages = ",".join(ordered)

        self._stage_fingerprints.clear()
        self._stage_recipe_fingerprints.clear()

    # --------------------------------------------------------

    def run(self):

        if self.args.resume_status:
            self.rehydrate_resume_status_context()

        stages = self.selected_stages()

        if self.args.list_stages:

            for stage in PIPELINE_STAGES:

                suffix = (
                    " (visibility FITS only)"
                    if stage == "ingest"
                    else " (optional)"
                    if stage in {"source_detection", "ddcal"}
                    else ""
                )
                print(f"{stage}{suffix}")

            return

        self.discard_unused_configured_primary_input(stages)
        input_kind_before_rehydration = self.input_kind
        self.rehydrate_resume_run_context(stages)
        if self.input_kind != input_kind_before_rehydration:
            # FITS inputs may need the ingest stage reintroduced now that the
            # manifest has restored the input kind. MS stage windows are
            # unchanged, but recomputing keeps one deterministic code path.
            stages = self.selected_stages()
            self.discard_unused_configured_primary_input(stages)
            self.rehydrate_resume_run_context(stages)

        if self.args.resume_status:
            self.print_resume_status(stages)
            return

        self.load_casa_module()
        self.preflight(stages)
        self.acquire_run_lock()

        try:

            self.prepare_working_copy(stages)

            self.validate_outdir_support(stages)

            if "inspection" not in stages:
                self.validate_downstream_target_selection(
                    stages,
                    defer_until_flagging=("flagging" in stages),
                )

            self.print_plan(stages)

            unavailable_stages = set()

            shared_stages = [
                stage for stage in stages
                if stage not in PER_TARGET_STAGES
            ]
            target_stages = [
                stage for stage in stages
                if stage in PER_TARGET_STAGES
            ]

            for stage in shared_stages:
                self.run_selected_stage(stage, stages, unavailable_stages)

            if target_stages:
                self.run_per_target_stages(
                    stages,
                    target_stages,
                    set(unavailable_stages),
                )

            if self.args.dry_run:

                print("\nDry run complete. No stages were executed.")

            else:

                print("\nPipeline run complete.")
                self.print_final_products(stages)

        finally:

            self.release_run_lock()


# ============================================================
# CLI
# ============================================================

def build_parser():

    parser = TrackingArgumentParser(
        description=(
            "Run the uGMRT, MeerKAT, or LOFAR continuum pipeline stages."
        )
    )

    parser.add_argument(
        "ms",
        nargs="?",
        default=None,
        metavar="INPUT",
        help=(
            "Input Measurement Set directory or UVFITS/FITS-IDI visibility "
            "file. Visibility FITS is imported to a managed MS first."
        ),
    )

    parser.add_argument(
        "--outdir",
        default=RUN_CONFIG.get(
            "outdir",
            DEFAULT_OUTDIR,
        ),
        help="Pipeline output directory.",
    )

    parser.add_argument(
        "--fits-format",
        choices=["auto", "uvfits", "fitsidi", "fits-idi", "lta"],
        default=RUN_CONFIG.get("fits_format", "auto"),
        help=(
            "Visibility-FITS dialect. Auto detects random-groups UVFITS or "
            "FITS-IDI UV_DATA, and takes a name ending in .lta, .ltb or "
            ".lta_N as a GMRT LTA file; ordinary image FITS is rejected."
        ),
    )

    parser.add_argument(
        "--listscan",
        default=RUN_CONFIG.get("gmrt_listscan") or "listscan",
        metavar="PROGRAM",
        help=(
            "NCRA's listscan, for GMRT LTA input: a path, a name on PATH, "
            "or a whole command (e.g. one that runs it in a container)."
        ),
    )

    parser.add_argument(
        "--gvfits",
        default=RUN_CONFIG.get("gmrt_gvfits") or "gvfits",
        metavar="PROGRAM",
        help="NCRA's gvfits, for GMRT LTA input; given like --listscan.",
    )

    parser.add_argument(
        "--keep-lta-fits",
        action="store_true",
        default=bool(RUN_CONFIG.get("keep_lta_fits", False)),
        help=(
            "Keep the FITS converted from an LTA file in <outdir>/ingest/lta "
            "after the import (it is as large as the data)."
        ),
    )

    parser.add_argument(
        "--imported-ms",
        default=RUN_CONFIG.get("imported_ms"),
        metavar="PATH",
        help=(
            "Destination MS for visibility-FITS ingestion. Default: "
            "<outdir>/ingest/<input-stem>.ms."
        ),
    )

    parser.add_argument(
        "--gmrt-flagfile",
        action="append",
        default=list(RUN_CONFIG.get("gmrt_flagfiles") or []),
        metavar="PATH",
        help=(
            "uGMRT online flag file passed to CASA importgmrt. Repeat for "
            "multiple files. Valid only for uGMRT UVFITS ingestion."
        ),
    )

    parser.add_argument(
        "--telescope",
        choices=TELESCOPE_PROFILE_MODULE.TELESCOPE_CHOICES,
        default=RUN_CONFIG.get("telescope", "auto"),
        help=(
            "Telescope identity. Auto reads the MS metadata; an explicit "
            "value is verified against it."
        ),
    )

    parser.add_argument(
        "--receiver-band",
        default=RUN_CONFIG.get("receiver_band"),
        help=(
            "Receiver band override (for example B4, UHF, L, HBA). "
            "Required for ambiguous MeerKAT overlap sub-bands."
        ),
    )

    parser.add_argument(
        "--data-state",
        choices=TELESCOPE_PROFILE_MODULE.DATA_STATE_CHOICES,
        default=RUN_CONFIG.get("data_state", "raw"),
        help=(
            "Use 'calibrated' for a target MS produced by a suitable "
            "upstream workflow such as LOFAR LINC/DP3."
        ),
    )

    parser.add_argument(
        "--imager",
        choices=TELESCOPE_PROFILE_MODULE.IMAGER_CHOICES,
        default=RUN_CONFIG.get("imager", "auto"),
        help=(
            "Auto selects CASA tclean for uGMRT, and WSClean for MeerKAT "
            "and LOFAR. uGMRT is always imaged with CASA; MeerKAT may use "
            "--imager casa; LOFAR requires WSClean."
        ),
    )

    parser.add_argument(
        "--wsclean",
        default=RUN_CONFIG.get("wsclean", "wsclean"),
        help="WSClean executable used by the LOFAR imaging backend.",
    )

    parser.add_argument(
        "--wsclean-data-column",
        default=RUN_CONFIG.get("wsclean_data_column", "CORRECTED_DATA"),
        help="Data column read by WSClean.",
    )

    parser.add_argument(
        "--lofar-beam-mode",
        choices=["full", "differential", "none"],
        default=RUN_CONFIG.get("lofar_beam_mode"),
        help=(
            "Required for LOFAR imaging: full trusts MS beam metadata; "
            "differential is an expert forced override; none makes no "
            "PB-corrected product."
        ),
    )

    parser.add_argument(
        "--profile",
        choices=EXECUTION_PROFILE_NAMES,
        default=RUN_CONFIG.get("profile", "science"),
        help=(
            "Execution profile. 'laptop' is the fast, memory-bounded "
            "configuration for a 36-GB Mac; 'science' preserves full "
            "resolution."
        ),
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=RUN_CONFIG.get("threads"),
        metavar="N",
        help=(
            "CASA/OpenMP/BLAS threads for serial work. Default: all CPUs "
            "available to the local process or Slurm allocation."
        ),
    )

    parser.add_argument(
        "--flagging-threads",
        type=int,
        default=RUN_CONFIG.get("flagging_threads", 4),
        metavar="N",
        help=(
            "Memory-bounded numerical threads for CASA flagging. Default: 4; "
            "this changes execution parallelism, not science thresholds."
        ),
    )

    parser.add_argument(
        "--mpi",
        choices=["auto", "on", "off"],
        default=str(RUN_CONFIG.get("mpi", "auto")).lower(),
        help=(
            "CASA MPI mode for final-imaging tclean. 'auto' enables MPI "
            "only in Slurm when a matching mpicasa is found; macOS safely "
            "uses serial CASA/OpenMP."
        ),
    )

    parser.add_argument(
        "--mpi-processes",
        type=int,
        default=RUN_CONFIG.get("mpi_processes"),
        metavar="N",
        help=(
            "CASA MPI process count including rank 0 (default: up to 16 "
            "CPUs from the allocation)."
        ),
    )

    parser.add_argument(
        "--mpicasa",
        default=RUN_CONFIG.get("mpicasa", "mpicasa"),
        metavar="PATH",
        help=(
            "mpicasa launcher. The executable beside --casa is preferred "
            "when this is left as 'mpicasa'."
        ),
    )

    parser.add_argument(
        "--copy-input-ms",
        action="store_true",
        default=bool(
            RUN_CONFIG.get(
                "copy_input_ms",
                True,
            )
        ),
        help=(
            "Copy the input MS before mutable stages and run the "
            "pipeline on that working copy."
        ),
    )

    parser.add_argument(
        "--no-copy-input-ms",
        action="store_false",
        dest="copy_input_ms",
        help=(
            "Skip a copy only for non-mutating calibrated-input stages. "
            "Raw-MS flagging/calibration reject this option."
        ),
    )

    parser.add_argument(
        "--working-ms",
        default=RUN_CONFIG.get("working_ms"),
        metavar="PATH",
        help=(
            "Destination for --copy-input-ms. Default: "
            "<outdir>/working_ms/<input-name>.ms."
        ),
    )

    parser.add_argument(
        "--target",
        default=RUN_CONFIG.get("target"),
        help=(
            "Explicit target name used when resolving split targets and "
            "downstream products. This always overrides automatic selection."
        ),
    )

    parser.add_argument(
        "--targets",
        default=RUN_CONFIG.get("targets"),
        metavar="A,B",
        help=(
            "Comma-separated science targets to self-calibrate and image. "
            "Defaults to every calibrated target."
        ),
    )

    parser.add_argument(
        "--target-selection",
        choices=("all", "auto", "explicit"),
        default=RUN_CONFIG.get("target_selection", "all"),
        help=(
            "Policy when multiple calibrated targets are available and "
            "neither --target nor --targets is given. 'all' images every "
            "target (the default); 'auto' ranks validated post-flagging "
            "field statistics and images one; 'explicit' asks for a target."
        ),
    )

    parser.add_argument(
        "--auto-target-max-flagged-fraction",
        type=float,
        default=RUN_CONFIG.get(
            "auto_target_max_flagged_fraction",
            0.80,
        ),
        metavar="FRACTION",
        help=(
            "Reject an automatically considered target when its post-flag "
            "flagged fraction exceeds this value (default: 0.80)."
        ),
    )

    parser.add_argument(
        "--flux-calibrator",
        default=RUN_CONFIG.get("flux_calibrator"),
        help=(
            "Explicit primary flux-calibrator field or comma-separated "
            "list for calibration."
        ),
    )

    parser.add_argument(
        "--bandpass-calibrator",
        default=RUN_CONFIG.get("bandpass_calibrator"),
        help="Explicit delay/bandpass-calibrator field.",
    )

    parser.add_argument(
        "--flux-standard",
        default=RUN_CONFIG.get("flux_standard"),
        help="Override the telescope profile's CASA setjy flux standard.",
    )

    parser.add_argument(
        "--flux-model",
        default=RUN_CONFIG.get("flux_model"),
        help=(
            "Validated absolute-flux CASA calibrator model image or component "
            "list (for example a MeerKAT UHF multi-component model)."
        ),
    )

    parser.add_argument(
        "--flux-model-type",
        choices=("casa-image", "component-list"),
        default=RUN_CONFIG.get("flux_model_type", "casa-image"),
        help="Representation supplied by --flux-model.",
    )

    parser.add_argument(
        "--secondary-flux-calibrator",
        default=RUN_CONFIG.get("secondary_flux_calibrator"),
        help=(
            "Optional additional flux-calibrator field or comma-separated "
            "list appended after --flux-calibrator."
        ),
    )

    parser.add_argument(
        "--phase-calibrator",
        default=RUN_CONFIG.get("phase_calibrator"),
        help=(
            "Explicit phase-calibrator field or comma-separated list for "
            "calibration."
        ),
    )

    parser.add_argument(
        "--target-phase-map",
        default=RUN_CONFIG.get("target_phase_map"),
        metavar="TARGET:PHASE[,TARGET:PHASE...]",
        help=(
            "Explicit target-to-phase-calibrator mapping used when "
            "multiple phase calibrators are present."
        ),
    )

    parser.add_argument(
        "--target-fields",
        default=RUN_CONFIG.get("target_fields"),
        metavar="FIELD[,FIELD...]",
        help="Explicit comma-separated target field names for calibration.",
    )

    parser.add_argument(
        "--reference-antenna",
        default=RUN_CONFIG.get("reference_antenna"),
        help="Explicit reference antenna for calibration and self-calibration.",
    )

    parser.add_argument(
        "--stages",
        default=RUN_CONFIG.get("stages"),
        help=(
            "Comma-separated stage list. Example: "
            "inspection,flagging,calibration"
        ),
    )

    parser.add_argument(
        "--from-stage",
        "--resume-from",
        dest="from_stage",
        choices=PIPELINE_STAGES,
        default=None,
        help=(
            "Start/resume from this stage. With --resume, missing input "
            "context is reconstructed from the run manifest and validated "
            "stage reports."
        ),
    )

    parser.add_argument(
        "--to-stage",
        choices=PIPELINE_STAGES,
        default=None,
        help="Stop after this stage.",
    )

    parser.add_argument(
        "--include-ddcal",
        action="store_true",
        default=bool(
            RUN_CONFIG.get(
                "include_ddcal",
                False,
            )
        ),
        help="Include optional DD calibration stage.",
    )

    parser.add_argument(
        "--no-ddcal",
        action="store_false",
        dest="include_ddcal",
        help="Do not include the optional DD calibration stage.",
    )

    parser.add_argument(
        "--dd-mode",
        choices=[
            "auto",
            "single",
            "survey",
        ],
        default=RUN_CONFIG.get(
            "dd_mode",
            DD_TASK_CONFIG.get("mode", "auto"),
        ),
        help=(
            "DD strategy: compact single-object, multi-facet survey, "
            "or automatic selection."
        ),
    )

    parser.add_argument(
        "--dd-engine",
        choices=[
            "auto",
            "prepare",
            "cubical",
            "killms-ddfacet",
        ],
        default=RUN_CONFIG.get(
            "dd_engine",
            DD_TASK_CONFIG.get(
                "engine",
                "auto",
            ),
        ),
        help="DD-cal backend when ddcal is included.",
    )

    parser.add_argument(
        "--solver-command-template",
        default=RUN_CONFIG.get("solver_command_template"),
        help=(
            "Expert CubiCal DD solve command template forwarded to "
            "ddcal.py; use only with --dd-engine cubical after site "
            "validation."
        ),
    )

    parser.add_argument(
        "--apply-command-template",
        default=RUN_CONFIG.get("apply_command_template"),
        help=(
            "Expert CubiCal DD apply command template forwarded to ddcal.py."
        ),
    )

    parser.add_argument(
        "--execute-ddcal",
        action="store_false",
        dest="dd_prepare_only",
        help=(
            "Allow ddcal.py to execute external solver commands. "
            "Default is prepare-only."
        ),
    )

    parser.add_argument(
        "--prepare-ddcal",
        action="store_true",
        dest="dd_prepare_only",
        help="Prepare DD-cal products without executing external solvers.",
    )

    parser.set_defaults(
        dd_prepare_only=bool(
            RUN_CONFIG.get(
                "dd_prepare_only",
                not DD_TASK_CONFIG.get(
                    "execute",
                    False,
                ),
            )
        ),
    )

    parser.add_argument(
        "--skip-dd-post-imaging",
        action="store_true",
        default=bool(
            RUN_CONFIG.get(
                "skip_dd_post_imaging",
                not DD_TASK_CONFIG.get(
                    "post_dd_imaging",
                    True,
                ),
            )
        ),
        help="Pass --skip-imaging to ddcal.py.",
    )

    parser.add_argument(
        "--run-dd-post-imaging",
        action="store_false",
        dest="skip_dd_post_imaging",
        help="Run post-DD imaging when ddcal is included.",
    )

    parser.add_argument(
        "--subbands",
        type=int,
        default=None,
        help=(
            "Number of in-band images when sub-band imaging is enabled "
            "with --make-subbands. Defaults are band dependent."
        ),
    )

    parser.add_argument(
        "--no-subbands",
        action="store_const",
        const=True,
        default=None,
        dest="no_subbands",
        help="Do not make in-band sub-band images (the default).",
    )

    parser.add_argument(
        "--make-subbands",
        action="store_const",
        const=False,
        dest="no_subbands",
        help=(
            "Opt in to in-band sub-band imaging, which is off by default. "
            "Combine with --subbands N to choose the count."
        ),
    )

    parser.add_argument(
        "--selfcal-ms",
        default=RUN_CONFIG.get("selfcal_ms"),
        help="Explicit calibrated target MS for selfcal.py.",
    )

    parser.add_argument(
        "--imaging-ms",
        default=RUN_CONFIG.get("imaging_ms"),
        help="Explicit MS for imaging.py.",
    )

    parser.add_argument(
        "--dd-ms",
        default=RUN_CONFIG.get("dd_ms"),
        help="Explicit MS for ddcal.py.",
    )

    parser.add_argument(
        "--casa-module",
        default=RUN_CONFIG.get("casa_module") or None,
        metavar="MODULE",
        help=(
            "Environment module that provides CASA (for example casa/6.7.0 "
            "on ilifu); the run loads it itself, on this node or in a SLURM "
            "job."
        ),
    )

    parser.add_argument(
        "--casa-python",
        default=RUN_CONFIG.get("casa_python") or None,
        metavar="COMMAND",
        help=(
            "Command that starts a Python with casatasks and casatools (a "
            "modular CASA), for example 'singularity exec IMAGE.sif python3'. "
            "Stages then run as COMMAND script args instead of casa -c."
        ),
    )

    parser.add_argument(
        "--casa",
        default=RUN_CONFIG.get(
            "casa",
            "casa",
        ),
        help="CASA executable.",
    )

    parser.add_argument(
        "--python",
        default=(
            RUN_CONFIG.get("python") or
            sys.executable
        ),
        help="Python executable for non-CASA stages.",
    )

    parser.add_argument(
        "--ncpu",
        type=int,
        default=RUN_CONFIG.get(
            "ncpu",
            DEFAULT_NCPU,
        ),
        help="CPU count passed to DD-cal.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running them.",
    )

    parser.add_argument(
        "--resume",
        action="store_true",
        default=bool(RUN_CONFIG.get("resume", True)),
        help=(
            "Resume automatically (default): skip completed stages only "
            "after validating their durable products, and restart the "
            "failed stage from its internal checkpoints."
        ),
    )

    parser.add_argument(
        "--no-resume",
        action="store_false",
        dest="resume",
        help=(
            "Invoke every selected stage even when its final products "
            "already validate. Individual stage checkpoints remain active."
        ),
    )

    parser.add_argument(
        "--resume-status",
        action="store_true",
        help=(
            "Show recorded stage status and whether each stage's products "
            "are safe to reuse, then exit."
        ),
    )

    parser.add_argument(
        "--progress-interval",
        type=float,
        default=float(
            RUN_CONFIG.get(
                "progress_interval_s",
                60.0,
            )
        ),
        metavar="SECONDS",
        help=(
            "Print a heartbeat while a stage is running. "
            "Set to 0 to disable (default: 60)."
        ),
    )

    parser.add_argument(
        "--list-stages",
        action="store_true",
        help="Print known stages and exit.",
    )

    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        default=bool(
            RUN_CONFIG.get(
                "continue_on_error",
                False,
            )
        ),
        help="Continue after a failed stage.",
    )

    parser.add_argument(
        "--write-sbatch",
        nargs="?",
        const="__CONFIG_DEFAULT__",
        default=None,
        metavar="PATH",
        help=(
            "Write an ilifu/Slurm sbatch file from config.py "
            "and exit. If PATH is omitted, use "
            "ILIFU_SLURM_CONFIG['script_path']."
        ),
    )

    launch = parser.add_mutually_exclusive_group()
    launch.add_argument(
        "--background",
        action="store_true",
        help=(
            "Run on this machine (for example the Jupyter node), detached "
            "from the terminal and at low priority; output goes to "
            "<outdir>/pipeline_run.log."
        ),
    )
    launch.add_argument(
        "--submit",
        action="store_true",
        help=(
            "Write the SLURM job from config.yaml (or --write-sbatch PATH), "
            "check that a node of the partition can take it, and submit it "
            "with sbatch."
        ),
    )

    parser.add_argument(
        "--sbatch-dry-run",
        action="store_true",
        help=(
            "When used with --write-sbatch, make the batch job "
            "run pipeline.py --dry-run."
        ),
    )

    parser.add_argument(
        "--sbatch-extra-arg",
        action="append",
        default=[],
        help=(
            "Extra argument to append to the pipeline command "
            "inside the generated sbatch file. May be repeated. "
            "For values beginning with --, use "
            "--sbatch-extra-arg=--flag."
        ),
    )

    return parser


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":

    arguments = build_parser().parse_args()

    try:

        if arguments.submit:

            sys.exit(submit_to_slurm(arguments))

        if arguments.background:

            sys.exit(start_in_background(arguments, sys.argv[1:]))

        if arguments.write_sbatch is not None:

            write_sbatch_from_config(arguments)
            sys.exit(0)

        PipelineRunner(arguments).run()

    except Exception as exc:

        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )
        sys.exit(1)
