#!/usr/bin/env python3
#
# flagging.py
#
# Telescope-aware uGMRT/MeerKAT CASA flagging module
# (calibrated LOFAR inputs bypass this raw-data stage)
#
#--------------------------------------------
# - Restart checkpoints
# - Iterative tfcrop/rflag
# - SPW statistics
# - Antenna statistics
# - Baseline statistics
# - Channel statistics (framework)
# - Edge-channel flagging
# - Improved bad antenna detection
# - Automatic reference antenna recommendation
# - QA plot generation (PNG)
# - JSON/ASCII reports
#---------------------------------------------
# Author: Sushant Dutta
# Compatible with CASA 6+

import importlib.util
import json
import math
import os


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of failing inside a stage.
PIPELINE_INTERFACE_VERSION = 1
import sys
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

if (
    sys.path and
    os.path.abspath(sys.path[0]) == SCRIPT_DIR
):

    sys.path.pop(0)

import argparse


def load_casa_compat_module():

    compat_path = Path(SCRIPT_DIR) / "casa_compat.py"
    spec = importlib.util.spec_from_file_location(
        "ugmrt_casa_compat",
        str(compat_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


CASA_COMPAT = load_casa_compat_module()
flagdata = CASA_COMPAT.import_casa_task("flagdata")
flagmanager = CASA_COMPAT.import_casa_task("flagmanager")
from casatools import table

PLOTTING_IMPORT_ERROR = None

try:
    import matplotlib

    # CASA is normally non-interactive, but set the backend explicitly so a
    # batch run never depends on an X/Quartz display being available.
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
except Exception as exc:
    plt = None
    PLOTTING_IMPORT_ERROR = str(exc)


def load_config_module():

    config_path = Path(SCRIPT_DIR) / "config.py"

    if not config_path.exists():

        return None

    removed_script_dir = False

    if (
        sys.path and
        os.path.abspath(sys.path[0]) == SCRIPT_DIR
    ):

        sys.path.pop(0)
        removed_script_dir = True

    try:

        spec = importlib.util.spec_from_file_location(
            "ugmrt_pipeline_config",
            str(config_path),
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        return module

    finally:

        if removed_script_dir:

            sys.path.insert(0, SCRIPT_DIR)


def load_telescope_profiles_module():

    module_path = Path(SCRIPT_DIR) / "telescope_profiles.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_telescope_profiles",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_flagging_selectors_module():

    module_path = Path(SCRIPT_DIR) / "flagging_selectors.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_flagging_selectors",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def deep_merge(base, override):

    result = dict(base or {})
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


CONFIG = load_config_module()
RUN_CONFIG = (
    getattr(CONFIG, "PIPELINE_RUN_CONFIG", {})
    if CONFIG is not None
    else {}
)
FLAGGING_DEFAULTS = (
    getattr(CONFIG, "FLAGGING_CONFIG", {})
    if CONFIG is not None
    else {}
)
TELESCOPE_FLAGGING_DEFAULTS = (
    getattr(CONFIG, "TELESCOPE_FLAGGING_CONFIG", {})
    if CONFIG is not None
    else {}
)
TELESCOPE_PROFILE_MODULE = load_telescope_profiles_module()
FLAGGING_SELECTORS = load_flagging_selectors_module()

# The comments this pipeline gives the flag versions it saves. Finding one on
# an MS whose flagging has not started in this run means an earlier run
# already flagged it.
PIPELINE_FLAG_VERSION_COMMENTS = (
    "radio continuum pipeline flag state",
    "Before corrected-data tfcrop/rflag and round-2 solves",
    "Before self-calibration residual-RFI flagging",
)


def earlier_pipeline_flag_versions(ms):
    """Flag versions this pipeline saved on ``ms``, read from its list.

    CASA keeps saved versions in ``<ms>.flagversions/FLAG_VERSION_LIST``, one
    ``name : comment`` per line. Reading the list needs no CASA and changes
    nothing.
    """

    path = Path(str(ms) + ".flagversions") / "FLAG_VERSION_LIST"
    try:
        lines = path.read_text(errors="replace").splitlines()
    except OSError:
        return []
    found = []
    for line in lines:
        name, _, comment = line.partition(" : ")
        if any(marker in comment for marker in PIPELINE_FLAG_VERSION_COMMENTS):
            found.append(name.strip())
    return found


# Manual flags are observation-specific: data the log or an earlier run shows
# to be bad, which the automatic flaggers cannot recognise. On the Band 5 run
# of 2026-10-25 the first ten minutes of the 3C48 scan gave amplitude solves
# with SNR < 3 on 22 of 26 antennas (phases were fine), and the bandpass was
# solved over them.
MANUAL_FLAG_SELECTION_KEYS = frozenset({
    "field", "scan", "antenna", "timerange", "spw", "correlation",
    "uvrange", "observation", "intent", "array", "feed",
})
# A single command that selects more than this share of the observation is
# far more likely a typing error (a missing field or time range) than a
# considered flag.
MANUAL_FLAG_MAX_FRACTION = 0.25


def parse_manual_flag_command(command):
    """Read one flagdata command in CASA's list syntax into task arguments.

    ``mode='manual' field='3C48' timerange='...' reason='...'``. Only manual
    flags are accepted, and each must select something: a command naming no
    field, scan, antenna, time or channels would flag the whole observation.
    Returns the selection (flagdata arguments) and the stated reason.
    """

    import shlex

    try:
        tokens = shlex.split(str(command))
    except ValueError as exc:
        raise ValueError(f"Manual flag command {command!r}: {exc}") from None
    selection = {}
    reason = None
    for token in tokens:
        key, separator, value = token.partition("=")
        key = key.strip().lower()
        if not separator or not key:
            raise ValueError(
                f"Manual flag command {command!r}: {token!r} is not "
                "key='value'."
            )
        if key == "mode":
            if value.strip().lower() != "manual":
                raise ValueError(
                    f"Manual flag command {command!r}: only mode='manual' is "
                    "accepted here; the automatic flaggers have settings of "
                    "their own."
                )
            continue
        if key == "reason":
            reason = value
            continue
        if key not in MANUAL_FLAG_SELECTION_KEYS:
            raise ValueError(
                f"Manual flag command {command!r}: unknown selection "
                f"{key!r}; use {', '.join(sorted(MANUAL_FLAG_SELECTION_KEYS))}."
            )
        if key in selection:
            raise ValueError(
                f"Manual flag command {command!r}: {key!r} is given twice."
            )
        if not value.strip():
            raise ValueError(
                f"Manual flag command {command!r}: {key!r} is empty."
            )
        selection[key] = value.strip()
    if not selection:
        raise ValueError(
            f"Manual flag command {command!r} selects nothing, so it would "
            "flag the whole observation."
        )
    return selection, reason


def resolve_execution_profile(profile):
    """Resolve a profile without changing the science flagging defaults."""

    name = str(profile or RUN_CONFIG.get("profile", "science"))
    name = name.strip().lower()

    if CONFIG is not None and hasattr(CONFIG, "get_execution_profile"):
        return name, CONFIG.get_execution_profile(name)

    if name == "science":
        return name, {}

    raise ValueError(
        "Execution profiles require config.py; use the 'science' profile "
        "when running flagging.py as a standalone script without it."
    )


class CASAFlagger:

    def __init__(self,
                 ms,
                 outdir="Products",
                 max_iterations=None,
                 bad_ant_threshold=None,
                 edge_fraction=None,
                 profile=None,
                 provenance_fingerprint=None):

        self.ms = ms
        self.provenance_fingerprint = (
            str(provenance_fingerprint).strip()
            if provenance_fingerprint
            else None
        )
        self.profile, self.profile_settings = (
            resolve_execution_profile(profile)
        )
        profile_flagging = self.profile_settings.get("flagging", {})
        self.max_iterations = (
            max_iterations
            if max_iterations is not None
            else profile_flagging.get(
                "max_iterations",
                FLAGGING_DEFAULTS.get(
                    "max_iterations",
                    4,
                ),
            )
        )
        self.bad_ant_threshold = (
            bad_ant_threshold
            if bad_ant_threshold is not None
            else FLAGGING_DEFAULTS.get(
                "bad_ant_threshold",
                0.7,
            )
        )
        self.edge_fraction = edge_fraction

        # Optional, conservative unsupervised ML complement to CASA's
        # deterministic tfcrop/rflag passes. It operates on persistent
        # per-channel statistics and is deliberately bounded so it cannot
        # flag a large fraction of a spectral window automatically.
        self.ml_config = {
            "enabled": False,
            "algorithm": "isolation_forest",
            "contamination": 0.01,
            "max_flag_fraction_per_spw": 0.01,
            "min_channels": 64,
            "min_valid_samples_per_channel": 256,
            "sample_rows": 50000,
            "chunk_rows": 2000,
            "random_state": 42,
            "n_estimators": 200,
            "edge_guard_fraction": 0.03,
            "require_persistence_sigma": 4.0,
        }
        self.ml_config.update(FLAGGING_DEFAULTS.get("ml_rfi", {}))

        # load inspection

        inspection_file = (
            Path(outdir) /
            "inspection" /
            "inspection_report.json"
        )

        if inspection_file.exists():

            with open(inspection_file) as f:
                self.inspection = json.load(f)

        else:

            self.inspection = {}

        self.telescope = TELESCOPE_PROFILE_MODULE.telescope_name_from_report(
            self.inspection
        )
        if self.telescope is None:
            raise RuntimeError(
                "Inspection report has no supported telescope identity. "
                "Rerun inspection before flagging."
            )

        self.flagging_defaults = deep_merge(
            FLAGGING_DEFAULTS,
            TELESCOPE_FLAGGING_DEFAULTS.get(self.telescope, {}),
        )

        if max_iterations is None:
            self.max_iterations = profile_flagging.get(
                "max_iterations",
                self.flagging_defaults.get("max_iterations", 4),
            )
        if bad_ant_threshold is None:
            self.bad_ant_threshold = self.flagging_defaults.get(
                "bad_ant_threshold",
                0.7,
            )
        self.ml_config.update(
            self.flagging_defaults.get("ml_rfi", {})
        )

        #
        # Band-specific edge channel fractions
        #

        EDGE_FRACTION = {

            "b2": 0.07,
            "band2": 0.07,

            "b3": 0.05,
            "band3": 0.05,

            "b4": 0.05,
            "band4": 0.05,

            "b5": 0.03,
            "band5": 0.03,
        }

        band = self.inspection.get("band", {})

        if isinstance(band, dict):

            band = (
                band.get("name") or
                band.get("band") or
                ""
            )

        band_upper = str(band).upper()
        band_lower = band_upper.lower()
        self.band = band_upper
        self.known_frequency_masks = list(
            self.flagging_defaults.get(
                "known_frequency_masks",
                {},
            ).get(self.band, [])
        )
        self.precal_rflag_enabled = bool(
            self.flagging_defaults.get("precal_rflag", True)
        )

        self.edge_fraction = (
            edge_fraction
            if edge_fraction is not None
            else self.flagging_defaults.get(
                "edge_fraction",
                {},
            ).get(
                band_upper,
                EDGE_FRACTION.get(
                    band_lower,
                    0.05,
                )
            )
        )

        autoflag_config = self.flagging_defaults.get(
            "autoflag",
            {},
        )
        self.autoflag_config = dict(
            autoflag_config.get(
                "default",
                {},
            )
        )
        self.autoflag_config.update(
            autoflag_config.get(
                self.band,
                {},
            )
        )
        self.quack_interval_s = float(
            self.flagging_defaults.get(
                "quack_interval_s",
                0.0,
            )
        )
        # Parsed now, so a mistyped command stops the stage before any flag
        # is written.
        self.manual_flag_commands = [
            str(command)
            for command in (
                self.flagging_defaults.get("manual_commands") or []
            )
        ]
        for command in self.manual_flag_commands:
            parse_manual_flag_command(command)
        self.flag_autocorrelations_enabled = bool(
            self.flagging_defaults.get(
                "flag_autocorrelations",
                True,
            )
        )
        self.flag_shadowed_enabled = bool(
            self.flagging_defaults.get(
                "flag_shadowed",
                True,
            )
        )
        self.shadow_tolerance_m = float(
            self.flagging_defaults.get(
                "shadow_tolerance_m",
                0.0,
            )
        )
        self.extend_growtime = float(
            self.flagging_defaults.get(
                "extend_growtime",
                80.0,
            )
        )
        self.extend_growfreq = float(
            self.flagging_defaults.get(
                "extend_growfreq",
                80.0,
            )
        )
        self.extend_polarizations = bool(
            self.flagging_defaults.get(
                "extend_polarizations",
                True,
            )
        )
        self.convergence_new_fraction = float(
            self.flagging_defaults.get(
                "convergence_new_fraction",
                0.001,
            )
        )
        self.max_iteration_new_fraction = float(
            self.flagging_defaults.get(
                "max_iteration_new_fraction",
                0.20,
            )
        )
        self.max_total_flagged_fraction = float(
            self.flagging_defaults.get(
                "max_total_flagged_fraction",
                0.70,
            )
        )
        self.max_group_iteration_new_fraction = float(
            self.flagging_defaults.get(
                "max_group_iteration_new_fraction",
                0.35,
            )
        )
        self.max_group_total_flagged_fraction = float(
            self.flagging_defaults.get(
                "max_group_total_flagged_fraction",
                0.85,
            )
        )
        self.max_extension_new_fraction = float(
            self.flagging_defaults.get(
                "max_extension_new_fraction",
                0.10,
            )
        )
        self.max_group_extension_new_fraction = float(
            self.flagging_defaults.get(
                "max_group_extension_new_fraction",
                0.20,
            )
        )
        self._parallel_correlation = None

        for name, value in [
            (
                "max_iteration_new_fraction",
                self.max_iteration_new_fraction,
            ),
            (
                "max_total_flagged_fraction",
                self.max_total_flagged_fraction,
            ),
            (
                "max_extension_new_fraction",
                self.max_extension_new_fraction,
            ),
            (
                "max_group_iteration_new_fraction",
                self.max_group_iteration_new_fraction,
            ),
            (
                "max_group_total_flagged_fraction",
                self.max_group_total_flagged_fraction,
            ),
            (
                "max_group_extension_new_fraction",
                self.max_group_extension_new_fraction,
            ),
            (
                "convergence_new_fraction",
                self.convergence_new_fraction,
            ),
        ]:
            if not 0.0 <= value <= 1.0:
                raise ValueError(
                    f"{name} must be between 0 and 1; received {value}."
                )

        for name, value in [
            ("extend_growtime", self.extend_growtime),
            ("extend_growfreq", self.extend_growfreq),
        ]:
            if not 0.0 <= value <= 100.0:
                raise ValueError(
                    f"{name} must be between 0 and 100; received {value}."
                )

        self.max_iterations = max(1, int(self.max_iterations))

        self.root = Path(outdir) / "flagging"
        self.plotdir = self.root / "plots"
        self.checkdir = self.root / "checkpoints"

        self.root.mkdir(parents=True, exist_ok=True)
        self.plotdir.mkdir(parents=True, exist_ok=True)
        self.checkdir.mkdir(parents=True, exist_ok=True)

        self.donefile = self.root / "flagging.done"

    # --------------------------------------------------
    # Restart logic
    # --------------------------------------------------

    def done(self):
        return self.donefile.exists()

    def mark_done(self):
        self.donefile.touch()

    # --------------------------------------------------

    def step_done(self, name):

        return (self.checkdir / f"{name}.done").is_file()

    def mark_step_done(self, name):

        (self.checkdir / f"{name}.done").touch()

    def run_checkpointed_step(self,
                              name,
                              label,
                              operation,
                              *args,
                              **kwargs):
        """Run a mutating flagging step once, then persist its completion."""

        if self.step_done(name):

            print(
                f"[flagging] {label} already completed; reusing checkpoint.",
                flush=True,
            )
            return None

        result = self.run_step(
            label,
            operation,
            *args,
            **kwargs,
        )
        self.mark_step_done(name)

        return result

    # --------------------------------------------------

    def flagdata_chunk_plan(self):
        """Size flagdata time chunks so no single chunk exceeds memory.

        flagdata holds one time chunk per field and SPW in memory, and the
        pipeline used whole scans. For high-resolution data a single scan can
        exceed the node's allocation, and the kernel or Slurm then kills CASA
        with SIGKILL partway through a mutating pass. The plan keeps whole
        scans wherever they fit, bounds the chunk only where they do not,
        fails before any flag is written when even a short chunk cannot fit,
        and records the decision because ntime shapes TFCrop's time fits and
        the extension step's growtime.
        """

        cached = getattr(self, "_flagdata_chunk_plan", None)
        if cached is not None:
            return cached

        limit = FLAGGING_SELECTORS.detect_memory_limit_bytes(
            configured_gb=RUN_CONFIG.get("flagdata_memory_limit_gb"),
        )

        try:
            from casatools import msmetadata

            groups = FLAGGING_SELECTORS.summarize_scan_chunks(
                self.ms,
                msmetadata,
                table,
            )
        except Exception as exc:
            plan = {
                "ntime": "scan",
                "status": "estimation_failed",
                "reason": f"scan geometry could not be read: {exc}",
            }
        else:
            plan = FLAGGING_SELECTORS.plan_flagdata_time_chunk(
                groups,
                limit["limit_bytes"],
                memory_fraction=float(RUN_CONFIG.get(
                    "flagdata_memory_fraction",
                    FLAGGING_SELECTORS.FLAGDATA_MEMORY_FRACTION,
                )),
                bytes_per_element=float(RUN_CONFIG.get(
                    "flagdata_bytes_per_element",
                    FLAGGING_SELECTORS.FLAGDATA_BYTES_PER_ELEMENT,
                )),
            )

        plan["memory_limit"] = limit
        self._flagdata_chunk_plan = plan

        root = getattr(self, "root", None)
        if root is not None:
            self.save_json(plan, Path(root) / "flagdata_chunk_plan.json")

        gib = 1024 ** 3
        if plan["status"] == "insufficient_memory":
            raise RuntimeError(
                "flagdata cannot fit even "
                f"{plan['minimum_integrations']} integrations of the largest "
                "scan in memory: about "
                f"{plan['required_memory_bytes'] / gib:.0f} GiB is needed, "
                f"against a limit of {limit['limit_bytes'] / gib:.0f} GiB "
                f"({limit['limiting_source']}). Request more memory, or "
                "channel-average the Measurement Set before flagging. No "
                "flags were written."
            )

        if plan["status"] == "bounded":
            print(
                "[flagging] Whole-scan flagdata chunks would need about "
                f"{plan['largest_scan_chunk']['scan_chunk_bytes'] / gib:.0f} "
                f"GiB against a {plan['budget_bytes'] / gib:.0f} GiB budget; "
                f"using ntime={plan['ntime']} instead. Recorded in "
                "flagdata_chunk_plan.json.",
                flush=True,
            )
        elif plan["status"] == "estimation_failed":
            print(
                "WARNING: flagdata memory could not be estimated "
                f"({plan['reason']}); retaining whole-scan chunks.",
                flush=True,
            )

        return plan

    def flagdata_ntime(self):

        return self.flagdata_chunk_plan()["ntime"]

    def summary(self, **kwargs):
        return flagdata(
            vis=self.ms,
            mode="summary",
            action="calculate",
            ntime=self.flagdata_ntime(),
            **kwargs
        )

    def spw_summary(self, **kwargs):
        return flagdata(
            vis=self.ms,
            mode="summary",
            action="calculate",
            spwchan=True,
            ntime=self.flagdata_ntime(),
            **kwargs
        )

    def save_json(self, obj, outfile):

        outfile = Path(outfile)
        outfile.parent.mkdir(parents=True, exist_ok=True)
        temporary = outfile.with_name(
            f".{outfile.name}.tmp.{os.getpid()}"
        )

        with open(temporary, "w") as f:
            json.dump(
                obj,
                f,
                indent=4,
                default=str
            )
            f.flush()
            os.fsync(f.fileno())

        os.replace(temporary, outfile)

    @staticmethod
    def read_json(path, default):

        path = Path(path)

        if not path.exists():
            return default

        try:
            with open(path) as f:
                return json.load(f)
        except (OSError, ValueError, TypeError) as exc:
            print(
                f"WARNING: Could not read flagging QA input {path}: {exc}",
                flush=True,
            )
            return default

    @staticmethod
    def format_elapsed(seconds):

        if seconds < 60:
            return f"{seconds:.1f} s"

        return f"{seconds / 60.0:.1f} min"

    def run_step(self, label, operation, *args, **kwargs):
        """Run a CASA operation with visible start and completion messages."""

        print(f"[flagging] {label} started.", flush=True)
        started = time.monotonic()
        result = operation(*args, **kwargs)
        elapsed = time.monotonic() - started
        print(
            f"[flagging] {label} completed "
            f"({self.format_elapsed(elapsed)}).",
            flush=True,
        )

        return result

    # --------------------------------------------------
    # Flagging operations
    # --------------------------------------------------

    def refuse_previously_flagged_input(self):
        """Stop before stacking a second round of flagging on the first.

        A working copy reused by a run in a new --outdir (a fixed
        ``working_ms`` path, a FITS run's reused ``--imported-ms``, or
        ``copy_input_ms: false``) arrives carrying the earlier run's flags:
        its band cuts, its dead antennas and its interference flags. Flagging
        it again counts all of that as the input's own and flags more on
        top, and each further run erodes the data again.
        """

        versions = earlier_pipeline_flag_versions(self.ms)
        if not versions:
            return None
        listed = ", ".join(versions)
        if FLAGGING_DEFAULTS.get("allow_previously_flagged_input", False):
            print(
                "WARNING: flagging a Measurement Set that an earlier run of "
                f"this pipeline already flagged (flag versions: {listed}), "
                "because flagging_allow_previously_flagged_input is set.",
                flush=True,
            )
            return versions
        raise RuntimeError(
            f"{self.ms} already carries flags written by an earlier run of "
            f"this pipeline (flag versions: {listed}), and this run has not "
            "flagged it yet. Flagging it again would stack a second round of "
            "interference flagging on the first, and count the earlier "
            "run's band cuts and dead antennas as the input's own. This "
            "happens when a new --outdir reuses a working copy: a fixed "
            "`working_ms` path, a FITS run's `--imported-ms`, or "
            "`copy_input_ms: false`. If this MS is a working copy, start "
            "from a fresh one: rename it (with its .flagversions) out of the "
            "way, or give `working_ms` a new path, and the input is copied "
            "again. If it is your original data (copy_input_ms: false), do "
            "not delete it: re-import it from the raw GMRT data, or restore "
            "flags you saved before the first pipeline run. To flag it "
            "again on purpose, set flagging_allow_previously_flagged_input: "
            "true."
        )

    def checkpoint(self, name):

        flagmanager(
            vis=self.ms,
            mode="save",
            versionname=name,
            comment=f"radio continuum pipeline flag state: {name}",
            merge="replace",
        )

    def restore_checkpoint(self, name):

        flagmanager(
            vis=self.ms,
            mode="restore",
            versionname=name,
            merge="replace",
        )

    def parallel_hand_correlation(self):
        """Return the parallel-hand correlation selector present in the MS."""

        if self._parallel_correlation is not None:
            return self._parallel_correlation

        tool = table()
        opened = False

        try:
            tool.open(str(Path(self.ms) / "POLARIZATION"))
            opened = True
            values = np.asarray(
                tool.getcol("CORR_TYPE"),
                dtype=int,
            ).ravel()
            present = set(int(value) for value in values)

            circular = [
                name
                for code, name in [(5, "RR"), (8, "LL")]
                if code in present
            ]
            linear = [
                name
                for code, name in [(9, "XX"), (12, "YY")]
                if code in present
            ]

            if circular and linear:
                raise RuntimeError(
                    "The POLARIZATION table contains both circular- and "
                    "linear-feed parallel hands; split the polarization "
                    "setups before auto-flagging."
                )

            hands = circular or linear

            if len(hands) != 2:
                raise RuntimeError(
                    "Stokes-I auto-flagging requires both parallel hands "
                    "(RR/LL or XX/YY); the POLARIZATION table provides "
                    + (", ".join(hands) if hands else "neither")
                    + "."
                )

            self._parallel_correlation = ",".join(hands)

        except Exception as exc:
            # Falling back to every correlation is not a safe default. Cross
            # hands are noise-dominated for an unpolarized field, so their
            # tfcrop/rflag statistics differ from the parallel hands, and the
            # guarded extension propagates any flags found there back into
            # RR/LL with extendpols=True. calibration.py already fails closed
            # on exactly this condition; the pre-calibration stage must match.
            raise RuntimeError(
                "Cannot identify the parallel-hand correlations required for "
                f"Stokes-I auto-flagging: {exc}"
            ) from exc

        finally:
            if opened:
                tool.close()

        return self._parallel_correlation

    def absolute_parallel_hand_correlation(self):

        hands = self.parallel_hand_correlation()

        if not hands:
            raise RuntimeError(
                "No parallel-hand correlation selector is available for "
                "auto-flagging."
            )

        # CASA applies an expression operator such as ABS to every
        # correlation that follows it.  The valid two-hand form is therefore
        # ``ABS_RR,LL`` (or ``ABS_XX,YY``), not ``ABS_RR,ABS_LL``.
        return f"ABS_{hands}"

    def iteration_thresholds(self, iteration):
        """Return a reproducible, bounded threshold schedule."""

        tightening = max(
            0.0,
            float(
                self.autoflag_config.get(
                    "tighten_per_iteration",
                    0.50,
                )
            ),
        )
        offset = max(0, int(iteration) - 1) * tightening

        def tightened(key, default, floor_key, floor):
            initial = float(
                self.autoflag_config.get(
                    key,
                    default,
                )
            )
            minimum = float(
                self.autoflag_config.get(
                    floor_key,
                    floor,
                )
            )
            return max(minimum, initial - offset)

        return {
            "tfcrop_timecutoff": tightened(
                "tfcrop_timecutoff",
                4.5,
                "tfcrop_timecutoff_floor",
                4.0,
            ),
            "tfcrop_freqcutoff": tightened(
                "tfcrop_freqcutoff",
                4.0,
                "tfcrop_freqcutoff_floor",
                3.5,
            ),
            "tfcrop_maxnpieces": int(
                self.autoflag_config.get(
                    "tfcrop_maxnpieces",
                    7,
                )
            ),
            "rflag_timedevscale": tightened(
                "rflag_timedevscale",
                6.0,
                "rflag_timedevscale_floor",
                5.0,
            ),
            "rflag_freqdevscale": tightened(
                "rflag_freqdevscale",
                6.0,
                "rflag_freqdevscale_floor",
                5.0,
            ),
            "rflag_winsize": int(
                self.autoflag_config.get(
                    "rflag_winsize",
                    5,
                )
            ),
        }

    def flag_autocorrelations(self):

        if not self.flag_autocorrelations_enabled:
            return

        flagdata(
            vis=self.ms,
            mode="manual",
            ntime=self.flagdata_ntime(),
            autocorr=True,
            action="apply",
            flagbackup=False,
        )

    def flag_shadowed_data(self):

        if not self.flag_shadowed_enabled:
            return

        flagdata(
            vis=self.ms,
            mode="shadow",
            ntime=self.flagdata_ntime(),
            tolerance=self.shadow_tolerance_m,
            action="apply",
            flagbackup=False,
        )

    def quack(self):

        if self.quack_interval_s <= 0:
            print("Quack flagging disabled by telescope profile.")
            return

        flagdata(
            vis=self.ms,
            mode="quack",
            quackinterval=self.quack_interval_s,
            quackmode="beg",
            quackincrement=False,
            action="apply",
            flagbackup=False,
        )

    def clip(self):

        flagdata(
            vis=self.ms,
            mode="clip",
            ntime=self.flagdata_ntime(),
            datacolumn="DATA",
            clipzeros=True,
            correlation=self.absolute_parallel_hand_correlation(),
            action="apply",
            flagbackup=False,
        )

    def flag_edge_channels(self):
        """Remove receiver-invalid spectrum and per-SPW filter roll-off.

        These are two different defects and both must be excised.  The
        guaranteed-range mask removes digitized channels outside the range
        over which the receiver is specified (MeerKAT UHF/L, and uGMRT B4
        when a 400-MHz GWB setup reaches past 850 MHz).  It says nothing
        about the sub-band filter roll-off at the edges of each spectral
        window, which is present whether or not the observation lies inside
        the guaranteed range.  Running only the first check silently left
        every uGMRT B4 setup contained within 550-850 MHz with completely
        unflagged SPW edges.
        """

        band_profile = TELESCOPE_PROFILE_MODULE.get_band_profile(
            self.telescope,
            self.band,
        )

        science_range_report = None
        if (
            "guaranteed_freq_min" in band_profile and
            "guaranteed_freq_max" in band_profile
        ):
            science_range_report = (
                self.flag_outside_guaranteed_science_range(band_profile)
            )

        percentage_report = self.flag_percentage_edge_channels()

        # Start from the guaranteed-range report so its keys -- scheme,
        # science_range_MHz, per-SPW coverage -- stay at the top level where
        # existing readers and the validation contract expect them, then layer
        # the per-SPW roll-off mask on top rather than nesting both and
        # silently changing the published report's shape.
        combined = dict(science_range_report or {})
        combined.update({
            "telescope": self.telescope,
            "band": self.band,
            "guaranteed_science_range_mask": science_range_report,
            "percentage_edge_mask": percentage_report,
            "edge_flagging": (
                "performed"
                if (
                    (science_range_report or {}).get("spw_selection") or
                    (percentage_report or {}).get("spw_selection")
                )
                else "not_required"
            ),
            "edge_fraction": self.edge_fraction,
            "spw_string": ",".join(
                selection
                for selection in [
                    (science_range_report or {}).get("spw_selection", ""),
                    (percentage_report or {}).get("spw_selection", ""),
                ]
                if selection
            ),
        })
        self.save_json(
            combined,
            self.root / "edge_flagging_report.json",
        )
        return combined

    def flag_percentage_edge_channels(self):
        """Flag a band-dependent fraction of channels at both SPW edges."""

        if self.edge_fraction <= 0:
            print("Percentage edge-channel flagging disabled by profile.")
            return {
                "status": "disabled",
                "edge_fraction": self.edge_fraction,
                "spw_selection": "",
            }

        if "spws" not in self.inspection:

            print(
                "Inspection report missing. "
                "Skipping edge-channel flagging."
            )

            return {
                "status": "skipped_missing_inspection_spws",
                "edge_fraction": self.edge_fraction,
                "spw_selection": "",
            }

        spw_ranges = []

        for spw, info in self.inspection["spws"].items():

            nchan = info["nchan"]

            nedge = int(self.edge_fraction * nchan)

            if nedge < 1 or 2 * nedge >= nchan:
                continue

            spw_ranges.append(

                f"{spw}:0~{nedge-1};"
                f"{nchan-nedge}~{nchan-1}"
            )

        spw_string = ",".join(spw_ranges)

        if not spw_string:
            print(
                "No valid edge-channel ranges were generated; skipping.",
                flush=True,
            )
            return {
                "status": "no_valid_ranges",
                "edge_fraction": self.edge_fraction,
                "spw_selection": "",
            }

        flagdata(
            vis=self.ms,
            mode="manual",
            ntime=self.flagdata_ntime(),
            spw=spw_string,
            action="apply",
            flagbackup=False,
        )

        print(
            f"Flagged edge channels: {spw_string}"
        )

        return {
            "status": "applied",
            "edge_fraction": self.edge_fraction,
            "spw_selection": spw_string,
            # Legacy key kept so existing consumers keep working.
            "spw_string": spw_string,
        }

    def flag_outside_guaranteed_science_range(self, band_profile):

        allowed = {
            ("ugmrt", "B4"),
            ("meerkat", "MEERKAT_UHF"),
            ("meerkat", "MEERKAT_L"),
        }
        if (self.telescope, self.band) not in allowed:
            raise RuntimeError(
                "Guaranteed-range masking was requested for an unvalidated "
                f"telescope/band pair: {self.telescope}/{self.band}."
            )

        report = FLAGGING_SELECTORS.build_science_range_channel_masks(
            self.inspection.get("spws", {}),
            float(band_profile["guaranteed_freq_min"]) / 1e6,
            float(band_profile["guaranteed_freq_max"]) / 1e6,
        )
        report.update({
            "telescope": self.telescope,
            "band": self.band,
            "edge_flagging": (
                "performed" if report["spw_selection"] else "not_required"
            ),
            # Legacy report readers used spw_string for percentage edges.
            "spw_string": report["spw_selection"],
            "digitized_receiver_range_MHz": [
                float(band_profile["freq_min"]) / 1e6,
                float(band_profile["freq_max"]) / 1e6,
            ],
        })

        selection = report["spw_selection"]
        if selection:
            flagdata(
                vis=self.ms,
                mode="manual",
                ntime=self.flagdata_ntime(),
                spw=selection,
                action="apply",
                flagbackup=False,
            )
            print(
                f"Flagged channels outside the {self.telescope}/{self.band} "
                "guaranteed science "
                f"range: {selection}",
                flush=True,
            )
        else:
            print(
                f"All observed channels lie inside the {self.telescope}/"
                f"{self.band} guaranteed "
                "science range; no receiver-edge mask was applied.",
                flush=True,
            )

        self.save_json(
            report,
            self.root / "guaranteed_range_mask_report.json",
        )
        return report

    def flag_known_frequency_ranges(self):

        if not self.known_frequency_masks:
            print("No telescope-specific frequency masks are configured.")
            return

        # The bundled persistent masks are validated only for MeerKAT L.
        # Refuse to apply them if configuration leakage ever makes them
        # visible under another telescope or receiver band.
        if self.telescope != "meerkat" or self.band != "MEERKAT_L":
            report = {
                "telescope": self.telescope,
                "band": self.band,
                "configured_frequency_ranges": self.known_frequency_masks,
                "frequency_ranges": self.known_frequency_masks,
                "spw_selection": "",
                "status": "skipped_telescope_band_guard",
                "reason": (
                    "Static MeerKAT L masks never apply outside "
                    "meerkat/MEERKAT_L."
                ),
            }
            self.save_json(
                report,
                self.root / "known_frequency_mask_report.json",
            )
            print(report["reason"], flush=True)
            return report

        report = FLAGGING_SELECTORS.build_intersected_static_masks(
            self.inspection.get("spws", {}),
            self.known_frequency_masks,
        )
        report.update({
            "telescope": self.telescope,
            "band": self.band,
            # Compatibility alias for earlier report consumers.
            "frequency_ranges": self.known_frequency_masks,
        })
        selection = report["spw_selection"]
        if selection:
            flagdata(
                vis=self.ms,
                mode="manual",
                ntime=self.flagdata_ntime(),
                spw=selection,
                action="apply",
                flagbackup=False,
            )
        else:
            print(
                "Configured MeerKAT L masks do not overlap any observed "
                "channel; no static frequency flags were applied.",
                flush=True,
            )
        self.save_json(
            report,
            self.root / "known_frequency_mask_report.json",
        )
        return report

    @staticmethod
    def summary_counts(summary):
        """(flagged, total) visibilities of a CASA summary, or (0, 0)."""

        if not isinstance(summary, dict):
            return 0.0, 0.0
        values = []
        for key in ("flagged", "total"):
            value = summary.get(key, 0)
            if hasattr(value, "item"):
                value = value.item()
            try:
                values.append(float(value))
            except (TypeError, ValueError):
                values.append(0.0)
        return values[0], values[1]

    def flag_manual_commands(self, observation_total):
        """Apply the observation's manual flag commands, one by one.

        Each command's selection is counted before and after it is applied,
        so the report shows what it removed. A command that selects no data
        (a misspelt field, a time outside the observation) or more than a
        quarter of the observation stops the stage instead of being applied.
        """

        report_path = self.root / "manual_flag_report.json"
        if not self.manual_flag_commands:
            print("No manual flag commands are configured.", flush=True)
            report = {"status": "none", "commands": []}
            self.save_json(report, report_path)
            return report

        records = []
        for command in self.manual_flag_commands:
            selection, reason = parse_manual_flag_command(command)
            flagged, total = self.summary_counts(self.summary(**selection))
            if total <= 0:
                raise RuntimeError(
                    f"Manual flag command {command!r} selects no data in "
                    f"{self.ms}. Check its field, scan or time range against "
                    "the inspection report."
                )
            share = (
                total / float(observation_total)
                if observation_total else None
            )
            if share is not None and share > MANUAL_FLAG_MAX_FRACTION:
                raise RuntimeError(
                    f"Manual flag command {command!r} selects {share:.1%} of "
                    f"the observation (more than "
                    f"{MANUAL_FLAG_MAX_FRACTION:.0%}); it was not applied. "
                    "Narrow its selection."
                )
            flagdata(
                vis=self.ms,
                mode="manual",
                action="apply",
                flagbackup=False,
                **selection,
            )
            flagged_after, _ = self.summary_counts(self.summary(**selection))
            record = {
                "command": command,
                "selection": selection,
                "reason": reason,
                "selected_visibilities": total,
                "share_of_observation": share,
                "flagged_before": flagged / total,
                "flagged_after": flagged_after / total,
            }
            records.append(record)
            print(
                f"Manual flags: {command}: {total:.0f} visibilities "
                f"({(share or 0.0):.2%} of the observation), "
                f"{record['flagged_before']:.1%} flagged before, "
                f"{record['flagged_after']:.1%} after.",
                flush=True,
            )
        report = {"status": "applied", "commands": records}
        self.save_json(report, report_path)
        return report

    def tfcrop(self, iteration):

        thresholds = self.iteration_thresholds(iteration)

        flagdata(
            vis=self.ms,
            mode="tfcrop",
            datacolumn="DATA",
            correlation=self.absolute_parallel_hand_correlation(),
            ntime=self.flagdata_ntime(),
            combinescans=False,
            timecutoff=thresholds["tfcrop_timecutoff"],
            freqcutoff=thresholds["tfcrop_freqcutoff"],
            timefit="line",
            freqfit="poly",
            maxnpieces=thresholds["tfcrop_maxnpieces"],
            flagdimension="freqtime",
            usewindowstats="none",
            extendflags=False,
            action="apply",
            flagbackup=False,
        )

    def rflag(self, iteration):

        if not self.precal_rflag_enabled:
            print(
                "Pre-calibration RFlag is disabled by the telescope "
                "profile; it will run after a provisional bandpass."
            )
            return

        thresholds = self.iteration_thresholds(iteration)

        flagdata(
            vis=self.ms,
            mode="rflag",
            datacolumn="DATA",
            correlation=self.absolute_parallel_hand_correlation(),
            ntime=self.flagdata_ntime(),
            combinescans=False,
            winsize=thresholds["rflag_winsize"],
            timedevscale=thresholds["rflag_timedevscale"],
            freqdevscale=thresholds["rflag_freqdevscale"],
            extendflags=False,
            action="apply",
            flagbackup=False,
        )

    def extend(self):

        flagdata(
            vis=self.ms,
            mode="extend",
            ntime=self.flagdata_ntime(),
            combinescans=False,
            growtime=self.extend_growtime,
            growfreq=self.extend_growfreq,
            growaround=False,
            flagneartime=False,
            flagnearfreq=False,
            extendpols=self.extend_polarizations,
            action="apply",
            flagbackup=False,
        )

    @staticmethod
    def _robust_zscore(values):
        values = np.asarray(values, dtype=float)
        median = np.nanmedian(values)
        mad = np.nanmedian(np.abs(values - median))
        scale = 1.4826 * mad
        if not np.isfinite(scale) or scale <= 0:
            scale = np.nanstd(values)
        if not np.isfinite(scale) or scale <= 0:
            scale = 1.0
        return (values - median) / scale

    @staticmethod
    def _compress_channels(channels):
        channels = sorted(set(int(value) for value in channels))
        if not channels:
            return []
        ranges = []
        start = previous = channels[0]
        for channel in channels[1:]:
            if channel == previous + 1:
                previous = channel
                continue
            ranges.append((start, previous))
            start = previous = channel
        ranges.append((start, previous))
        return ranges

    def _ddid_to_spw(self):
        tool = table()
        try:
            tool.open(str(Path(self.ms) / "DATA_DESCRIPTION"))
            values = tool.getcol("SPECTRAL_WINDOW_ID")
            return {index: int(value) for index, value in enumerate(values)}
        finally:
            try:
                tool.close()
            except Exception:
                pass

    def _collect_ml_channel_features(self):
        """Collect bounded, streaming per-channel visibility features.

        The model sees persistent statistics only; it does not classify
        individual integrations. This makes it a conservative complement to
        tfcrop/rflag and avoids loading the complete MS into memory.
        """

        config = self.ml_config
        sample_rows = max(1, int(config.get("sample_rows", 50000)))
        chunk_rows = max(1, int(config.get("chunk_rows", 2000)))
        ddid_to_spw = self._ddid_to_spw()
        tool = table()
        aggregates = {}

        try:
            tool.open(str(self.ms))
            nrows = int(tool.nrows())
            if nrows <= 0:
                return {}
            stride = max(1, int(math.ceil(nrows / float(sample_rows))))

            for start in range(0, nrows, chunk_rows * stride):
                row_ids = np.arange(
                    start, min(nrows, start + chunk_rows * stride), stride,
                    dtype=int,
                )
                if row_ids.size == 0:
                    continue

                # CASA table getcol is contiguous, so read the enclosing block
                # and then sub-select the sampled rows in memory.
                block_start = int(row_ids[0])
                block_stop = int(row_ids[-1]) + 1
                block_nrow = block_stop - block_start
                local = row_ids - block_start
                data = tool.getcol("DATA", block_start, block_nrow)[:, :, local]
                flags = tool.getcol("FLAG", block_start, block_nrow)[:, :, local]
                ddids = tool.getcol(
                    "DATA_DESC_ID", block_start, block_nrow
                )[local]

                amplitudes = np.abs(data)
                finite = np.isfinite(amplitudes)
                usable = finite & (~flags) & (amplitudes > 0)

                for ddid in np.unique(ddids):
                    selector = ddids == ddid
                    if not np.any(selector):
                        continue
                    spw = ddid_to_spw.get(int(ddid), int(ddid))
                    nchan = amplitudes.shape[1]
                    record = aggregates.setdefault(
                        int(spw),
                        {
                            "sum_logamp": np.zeros(nchan, dtype=float),
                            "sum_sq_logamp": np.zeros(nchan, dtype=float),
                            "valid": np.zeros(nchan, dtype=np.int64),
                            "flagged": np.zeros(nchan, dtype=np.int64),
                            "total": np.zeros(nchan, dtype=np.int64),
                            "zeros": np.zeros(nchan, dtype=np.int64),
                        },
                    )
                    subset_amp = amplitudes[:, :, selector]
                    subset_usable = usable[:, :, selector]
                    subset_flags = flags[:, :, selector]
                    safe_log = np.zeros_like(subset_amp, dtype=float)
                    safe_log[subset_usable] = np.log10(subset_amp[subset_usable])
                    record["sum_logamp"] += np.sum(safe_log, axis=(0, 2))
                    record["sum_sq_logamp"] += np.sum(
                        safe_log * safe_log, axis=(0, 2)
                    )
                    record["valid"] += np.sum(subset_usable, axis=(0, 2))
                    record["flagged"] += np.sum(subset_flags, axis=(0, 2))
                    record["total"] += subset_flags.shape[0] * subset_flags.shape[2]
                    record["zeros"] += np.sum(
                        finite[:, :, selector] & (subset_amp == 0), axis=(0, 2)
                    )
        finally:
            try:
                tool.close()
            except Exception:
                pass

        return aggregates

    def ml_rfi_flagging(self):
        """Flag only persistent channel anomalies proposed by IsolationForest.

        This stage is intentionally conservative: it guards SPW edges,
        requires a robust high-amplitude persistence test, and caps the
        automatically flagged fraction in every SPW. If scikit-learn is not
        installed, the stage writes a report and leaves the data unchanged.
        """

        report_path = self.root / "ml_rfi_report.json"
        config = dict(self.ml_config)
        report = {
            "enabled": bool(config.get("enabled", False)),
            "algorithm": config.get("algorithm", "isolation_forest"),
            "config": config,
            "flagged_spw_channels": {},
            "status": "skipped",
        }

        if not report["enabled"]:
            report["reason"] = "Disabled in FLAGGING_CONFIG['ml_rfi']."
            self.save_json(report, report_path)
            print("[flagging] ML-assisted RFI flagging is disabled.", flush=True)
            return report

        try:
            from sklearn.ensemble import IsolationForest
        except Exception as exc:
            report["reason"] = f"scikit-learn unavailable: {exc}"
            self.save_json(report, report_path)
            print(
                "WARNING: ML-assisted RFI flagging skipped because "
                f"scikit-learn is unavailable: {exc}",
                flush=True,
            )
            return report

        aggregates = self._collect_ml_channel_features()
        selections = []
        diagnostics = {}

        for spw, record in sorted(aggregates.items()):
            valid = record["valid"].astype(float)
            total = np.maximum(record["total"].astype(float), 1.0)
            nchan = len(valid)
            minimum_channels = int(config.get("min_channels", 64))
            if nchan < minimum_channels:
                diagnostics[str(spw)] = {"status": "too_few_channels"}
                continue

            mean = np.divide(
                record["sum_logamp"], valid,
                out=np.full(nchan, np.nan), where=valid > 0,
            )
            second = np.divide(
                record["sum_sq_logamp"], valid,
                out=np.full(nchan, np.nan), where=valid > 0,
            )
            scatter = np.sqrt(np.maximum(second - mean * mean, 0.0))
            flag_fraction = record["flagged"] / total
            zero_fraction = record["zeros"] / total
            gradient = np.abs(np.gradient(mean))
            min_samples = int(
                config.get("min_valid_samples_per_channel", 256)
            )
            usable_channels = np.isfinite(mean) & (valid >= min_samples)

            edge_guard = max(1, int(
                float(config.get("edge_guard_fraction", 0.03)) * nchan
            ))
            usable_channels[:edge_guard] = False
            usable_channels[-edge_guard:] = False
            indices = np.flatnonzero(usable_channels)
            if len(indices) < minimum_channels:
                diagnostics[str(spw)] = {"status": "insufficient_samples"}
                continue

            matrix = np.column_stack([
                self._robust_zscore(mean[indices]),
                self._robust_zscore(scatter[indices]),
                self._robust_zscore(flag_fraction[indices]),
                self._robust_zscore(zero_fraction[indices]),
                self._robust_zscore(gradient[indices]),
            ])
            matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)
            contamination = min(0.1, max(0.001, float(
                config.get("contamination", 0.01)
            )))
            model = IsolationForest(
                n_estimators=int(config.get("n_estimators", 200)),
                contamination=contamination,
                random_state=int(config.get("random_state", 42)),
                n_jobs=1,
            )
            labels = model.fit_predict(matrix)
            scores = -model.score_samples(matrix)
            amplitude_z = self._robust_zscore(mean)
            persistence_sigma = float(
                config.get("require_persistence_sigma", 4.0)
            )
            proposed = indices[(labels < 0) & (amplitude_z[indices] >= persistence_sigma)]

            max_fraction = max(0.0, min(0.05, float(
                config.get("max_flag_fraction_per_spw", 0.01)
            )))
            max_channels = max(1, int(math.floor(max_fraction * nchan)))
            if proposed.size > max_channels:
                score_lookup = {int(index): float(score) for index, score in zip(indices, scores)}
                proposed = np.asarray(sorted(
                    proposed, key=lambda channel: score_lookup[int(channel)],
                    reverse=True,
                )[:max_channels], dtype=int)

            channel_ranges = self._compress_channels(proposed.tolist())
            for first, last in channel_ranges:
                channel_text = str(first) if first == last else f"{first}~{last}"
                selections.append(f"{spw}:{channel_text}")

            report["flagged_spw_channels"][str(spw)] = proposed.tolist()
            diagnostics[str(spw)] = {
                "status": "analysed",
                "nchan": nchan,
                "usable_channels": int(len(indices)),
                "proposed_channels": proposed.tolist(),
            }

        report["diagnostics"] = diagnostics
        report["spw_selection"] = ",".join(selections)

        if not selections:
            report["status"] = "completed_no_flags"
            self.save_json(report, report_path)
            print(
                "[flagging] ML-assisted RFI analysis found no persistent "
                "channel anomalies that passed the safeguards.",
                flush=True,
            )
            return report

        flagdata(
            vis=self.ms,
            mode="manual",
            ntime=self.flagdata_ntime(),
            spw=report["spw_selection"],
            flagbackup=False,
        )
        report["status"] = "flags_applied"
        self.save_json(report, report_path)
        print(
            "[flagging] ML-assisted RFI flags applied to: "
            + report["spw_selection"],
            flush=True,
        )
        return report

    # --------------------------------------------------
    # Statistics
    # --------------------------------------------------

    def inspected_antenna_names(self):

        antenna_file = (
            self.root.parent /
            "inspection" /
            "antennas.json"
        )

        if not antenna_file.exists():

            return []

        with open(antenna_file) as f:

            antennas = json.load(f)

        return [
            item.get("name")
            for item in antennas
            if item.get("name")
        ]

    def active_antenna_names(self):
        """Return antennas that occur in at least one MAIN-table row.

        The ANTENNA subtable can retain entries for antennas that were absent,
        fully removed, or excluded while producing a split target MS. Calling
        ``flagdata(antenna=...)`` for such entries makes CASA emit a SEVERE
        ``MSSelectionNullSelection`` before Python can catch the exception.
        """

        tool = table()
        main_open = False
        antenna_open = False

        try:
            tool.open(str(self.ms))
            main_open = True

            if tool.nrows() == 0:
                return set()

            # Read in bounded blocks and reduce each to its unique IDs. The
            # previous full-column read converted every row to a Python int
            # list, roughly 36 bytes per row per column, which is several GB
            # on a long, short-integration MeerKAT observation.
            nrows = int(tool.nrows())
            block = 5_000_000
            active_ids = set()
            for column in ("ANTENNA1", "ANTENNA2"):
                if column not in tool.colnames():
                    continue
                for start in range(0, nrows, block):
                    values = np.asarray(
                        tool.getcol(column, start, min(block, nrows - start))
                    )
                    active_ids.update(int(value) for value in np.unique(values))

            tool.close()
            main_open = False

            tool.open(str(Path(self.ms) / "ANTENNA"))
            antenna_open = True
            names = tool.getcol("NAME")
            if hasattr(names, "tolist"):
                names = names.tolist()

            active_names = {
                str(names[index])
                for index in active_ids
                if 0 <= index < len(names)
            }

            # Keep numeric selectors as well because some inspection reports
            # store antenna IDs rather than station names.
            active_names.update(str(index) for index in active_ids)
            return active_names

        except Exception as exc:
            print(
                "WARNING: Could not determine active antennas before QA: "
                f"{exc}. Falling back to the inspection list.",
                flush=True,
            )
            return None

        finally:
            if main_open or antenna_open:
                try:
                    tool.close()
                except Exception:
                    pass

    def antenna_stats(self):

        # ``antenna`` is a data-selection parameter in flagdata summary mode,
        # not a request for an antenna-axis breakdown. Pre-filter the
        # inspection list against antennas that actually occur in MAIN so CASA
        # is never asked to evaluate a zero-row selection.
        stats = {}
        inspected = self.inspected_antenna_names()
        active = self.active_antenna_names()

        if active is not None:
            skipped = [
                antenna for antenna in inspected
                if str(antenna) not in active
            ]
            inspected = [
                antenna for antenna in inspected
                if str(antenna) in active
            ]

            if skipped:
                print(
                    "Skipping antennas with no rows in the current MS: "
                    + ", ".join(map(str, skipped)),
                    flush=True,
                )

        for antenna in inspected:

            try:
                summary = self.summary(antenna=str(antenna))
            except Exception as exc:
                print(
                    f"WARNING: Could not summarise active antenna "
                    f"{antenna}: {exc}"
                )
                continue

            if not isinstance(summary, dict) or summary.get("total", 0) <= 0:
                print(
                    f"Skipping antenna {antenna}: summary contains no rows.",
                    flush=True,
                )
                continue

            stats[str(antenna)] = summary

        if not stats:

            print(
                "WARNING: No active per-antenna statistics are available; "
                "bad-antenna detection is skipped."
            )

        return stats

    def spw_stats(self):

        by_spw = {}

        for spw in sorted(
            (self.inspection.get("spws") or {}).keys(),
            key=lambda value: int(value),
        ):

            try:
                by_spw[str(spw)] = self.summary(
                    spw=str(spw),
                )
            except Exception as exc:
                print(
                    f"WARNING: Could not summarise spectral window "
                    f"{spw}: {exc}",
                    flush=True,
                )

        try:
            by_spw_channel = self.spw_summary()
        except Exception as exc:
            print(
                "WARNING: Could not calculate per-channel flagging "
                f"statistics: {exc}",
                flush=True,
            )
            by_spw_channel = {}

        return {
            "by_spw": by_spw,
            "by_spw_channel": by_spw_channel,
        }

    def field_stats(self):

        names = set()
        classified = self.inspection.get("fields", {})

        if isinstance(classified, dict):
            for value in classified.values():
                if isinstance(value, (list, tuple)):
                    names.update(
                        str(item)
                        for item in value
                        if str(item).strip()
                    )

        for key in (
            "flux_calibrators",
            "phase_calibrators",
            "target_fields",
        ):
            value = self.inspection.get(key, [])
            if isinstance(value, (list, tuple)):
                names.update(
                    str(item)
                    for item in value
                    if str(item).strip()
                )

        by_field = {}

        for name in sorted(names):
            try:
                by_field[name] = self.summary(
                    field=name,
                )
            except Exception as exc:
                print(
                    f"WARNING: Could not summarise field {name}: {exc}",
                    flush=True,
                )

        try:
            breakdown = self.summary(fieldcnt=True)
        except Exception as exc:
            print(
                "WARNING: Could not calculate per-field flagging "
                f"statistics: {exc}",
                flush=True,
            )
            breakdown = {}

        return {
            "by_field": by_field,
            "breakdown": breakdown,
        }

    def baseline_stats(self):

        try:
            return self.summary(basecnt=True)
        except Exception:
            return {}

    # --------------------------------------------------
    # QA logic
    # --------------------------------------------------

    def detect_bad_antennas(self, stats):

        bad = []

        #
        # Collect total visibility counts
        #

        totals = []

        for k, v in stats.items():

            if not isinstance(v, dict):
                continue

            total = v.get("total", 0)

            if total > 0:
                totals.append(total)

        #
        # No valid antennas
        #

        if len(totals) == 0:
            return [], 0

        #
        # Robust threshold
        #

        totals.sort()

        median_total = totals[len(totals)//2]

        for k, v in stats.items():

            if not isinstance(v, dict):
                continue

            total = v.get("total", 0)
            flagged = v.get("flagged", 0)

            if total <= 0:
                continue

            #
            # Ignore antennas with very few samples
            #

            if total < median_total / 2:
                continue

            frac = flagged / total

            if frac > self.bad_ant_threshold:
                bad.append(k)

        return bad, median_total

    def recommend_refant(self,
                      stats,
                      bad_ant=None):

        if bad_ant is None:
            bad_ant = []

        ranking = []

        for k, v in stats.items():

            if not isinstance(v, dict):
                continue

            total = v.get("total", 0)
            flagged = v.get("flagged", 0)

            if total <= 0:
                continue

            #
            # Exclude bad antennas
            #

            if k in bad_ant:
                continue

            frac = flagged / total

            ranking.append(
                (frac, -total, k)
            )

        ranking.sort()

        if not ranking:
            return "UNKNOWN"

        self.refant_ranking = [
            {"antenna": name, "flagged_fraction": frac, "total": -negative}
            for frac, negative, name in ranking
        ]

        # The least-flagged antenna is not automatically the best phase
        # reference. uGMRT arm antennas sit tens of km out, so they have lower
        # per-baseline S/N and stronger differential ionosphere than the
        # central-square antennas the telescope profile lists. Prefer a
        # profile antenna whenever its flag fraction is competitive; with
        # refantmode='strict' in calibration a poor reference now costs
        # flagged solutions rather than failing silently.
        best_fraction = ranking[0][0]
        tolerance = float(
            self.flagging_defaults.get(
                "refant_preference_tolerance",
                0.05,
            )
        )
        acceptable = {
            name
            for frac, _negative, name in ranking
            if frac <= best_fraction + tolerance
        }
        preferred = TELESCOPE_PROFILE_MODULE.get_telescope_profile(
            self.telescope
        )["preferred_refants"]

        fractions = {
            item["antenna"]: item["flagged_fraction"]
            for item in self.refant_ranking
        }

        for name in preferred:
            if name in acceptable:
                if name != ranking[0][2]:
                    print(
                        f"[flagging] Reference antenna {name} "
                        f"(flagged {fractions[name]:.4f}) preferred over "
                        f"least-flagged {ranking[0][2]} "
                        f"(flagged {best_fraction:.4f}); it is a "
                        f"{self.telescope} profile antenna and lies within "
                        f"the {tolerance:.4f} tolerance.",
                        flush=True,
                    )
                self.refant_selection = {
                    "recommended": name,
                    "method": "telescope_preferred_within_tolerance",
                    "least_flagged": ranking[0][2],
                    "tolerance": tolerance,
                }
                return name

        self.refant_selection = {
            "recommended": ranking[0][2],
            "method": "least_flagged",
            "least_flagged": ranking[0][2],
            "tolerance": tolerance,
        }
        return ranking[0][2]

    # --------------------------------------------------
    # Plots
    # --------------------------------------------------

    @staticmethod
    def summary_flagged_fraction(summary):
        """Return a finite flagged fraction from a CASA summary mapping."""

        if not isinstance(summary, dict):
            return None

        total = summary.get("total")
        flagged = summary.get("flagged")

        try:
            if hasattr(total, "item"):
                total = total.item()
            if hasattr(flagged, "item"):
                flagged = flagged.item()
            total = float(total)
            flagged = float(flagged)
        except (TypeError, ValueError):
            return None

        if total <= 0 or not math.isfinite(total):
            return None

        fraction = flagged / total

        if not math.isfinite(fraction):
            return None

        return fraction

    @classmethod
    def summary_group_fractions(cls, summary, category):
        """Return valid per-group fractions from a CASA summary category."""

        if not isinstance(summary, dict):
            return {}

        entries = summary.get(category, {})

        if not isinstance(entries, dict):
            return {}

        fractions = {}

        for key, value in entries.items():
            fraction = cls.summary_flagged_fraction(value)

            if fraction is not None:
                fractions[str(key)] = fraction

        return fractions

    def antenna_plot(self, stats):
        """Write per-antenna flag fractions, or return None when unavailable."""

        names = []
        fracs = []

        for k, v in (stats or {}).items():

            fraction = self.summary_flagged_fraction(v)

            if fraction is None:
                continue

            names.append(k)
            fracs.append(100.0 * fraction)

        if not names:
            return None

        figure, axis = plt.subplots(figsize=(max(10, len(names) * 0.35), 5))
        axis.bar(names, fracs, color="tab:blue")
        axis.set_ylabel("Flagged (%)")
        axis.set_xlabel("Antenna")
        axis.set_title("Flagged visibility fraction by antenna")
        for label in axis.get_xticklabels():
            label.set_rotation(90)
        axis.set_ylim(bottom=0)

        outfile = self.plotdir / "antenna_flag_fraction.png"

        figure.savefig(str(outfile), dpi=160, bbox_inches="tight")
        plt.close(figure)
        return outfile

    def category_plot(self,
                      stats,
                      category,
                      outfile_name,
                      title,
                      xlabel):
        """Write a flag-fraction bar plot for a CASA summary category."""

        if not isinstance(stats, dict):
            return None

        entries = stats.get(category, {})

        if not isinstance(entries, dict):
            return None

        labels = []
        fractions = []

        for label, summary in entries.items():

            fraction = self.summary_flagged_fraction(summary)

            if fraction is None:
                continue

            labels.append(str(label))
            fractions.append(100.0 * fraction)

        if not labels:
            return None

        figure, axis = plt.subplots(
            figsize=(max(8, len(labels) * 0.7), 5)
        )
        axis.bar(labels, fractions, color="tab:green")
        axis.set_ylabel("Flagged (%)")
        axis.set_xlabel(xlabel)
        axis.set_title(title)
        axis.set_ylim(bottom=0)

        for label in axis.get_xticklabels():
            label.set_rotation(45)
            label.set_horizontalalignment("right")

        outfile = self.plotdir / outfile_name
        figure.savefig(str(outfile), dpi=160, bbox_inches="tight")
        plt.close(figure)
        return outfile

    def spw_plot(self, stats):
        """Write the aggregate flag fraction of every spectral window."""

        return self.category_plot(
            stats,
            "by_spw",
            "spw_flag_fraction.png",
            "Flagged visibility fraction by spectral window",
            "Spectral window",
        )

    def field_plot(self, stats):
        """Write the flag fraction of every observing field."""

        return self.category_plot(
            stats,
            "by_field",
            "field_flag_fraction.png",
            "Flagged visibility fraction by field",
            "Field",
        )

    @classmethod
    def flagging_attribution(cls, pre, post):
        """Separate what the input arrived with from what this run flagged.

        A cumulative fraction mixes the two. The number that says whether
        this run over-flagged is the share of the data it received unflagged
        that it then flagged: (post - pre) / (1 - pre).
        """

        before = cls.summary_flagged_fraction(pre)
        after = cls.summary_flagged_fraction(post)
        if before is None or after is None or before >= 1.0:
            return None
        return {
            "input_flagged_fraction": before,
            "final_flagged_fraction": after,
            "this_run_fraction_of_received": max(0.0, after - before) / (
                1.0 - before
            ),
        }

    def flag_fraction_history_plot(self, pre, history, post):
        """Write the QA plot that is available even without antenna stats."""

        labels = []
        fractions = []

        pre_fraction = self.summary_flagged_fraction(pre)

        if pre_fraction is not None:
            labels.append("pre")
            fractions.append(100.0 * pre_fraction)

        for item in history or []:
            if not isinstance(item, dict):
                continue

            try:
                number = int(item.get("iteration"))
                fraction = float(item.get("cumulative_fraction"))
            except (TypeError, ValueError):
                continue

            if not math.isfinite(fraction):
                continue

            labels.append(f"iter {number}")
            fractions.append(100.0 * fraction)

        post_fraction = self.summary_flagged_fraction(post)

        if post_fraction is not None:
            labels.append("post")
            fractions.append(100.0 * post_fraction)

        if not labels:
            return None

        figure, axis = plt.subplots(figsize=(max(7, len(labels) * 1.25), 5))
        x_values = list(range(len(labels)))
        axis.plot(
            x_values,
            fractions,
            marker="o",
            linewidth=2,
            color="tab:orange",
        )
        axis.set_xticks(x_values)
        axis.set_xticklabels(labels)
        axis.set_ylabel("Cumulative flagged data (%)")
        axis.set_title("Flagging progression")
        axis.grid(axis="y", alpha=0.3)
        axis.set_ylim(bottom=0)
        attribution = self.flagging_attribution(pre, post)
        if attribution:
            axis.text(
                0.02, 0.04,
                f"{attribution['input_flagged_fraction'] * 100:.1f}% arrived "
                "already flagged; this run flagged "
                f"{attribution['this_run_fraction_of_received'] * 100:.1f}% "
                "of the data it received unflagged",
                transform=axis.transAxes, fontsize=8,
                bbox={"facecolor": "white", "alpha": 0.8, "linewidth": 0},
            )

        outfile = self.plotdir / "flagging_fraction_history.png"
        figure.savefig(str(outfile), dpi=160, bbox_inches="tight")
        plt.close(figure)
        return outfile

    def generate_qa_plots(self,
                          pre,
                          history,
                          post,
                          antenna_stats,
                          spw_stats=None,
                          field_stats=None):
        """Create all available flagging QA plots and record their status."""

        status = {
            "matplotlib_available": plt is not None,
            "plots": [],
            "warnings": [],
        }

        if plt is None:
            status["warnings"].append(
                "Matplotlib is unavailable: "
                f"{PLOTTING_IMPORT_ERROR or 'unknown import error'}"
            )
        else:
            for label, operation in [
                (
                    "flagging_fraction_history",
                    lambda: self.flag_fraction_history_plot(pre, history, post),
                ),
                (
                    "antenna_flag_fraction",
                    lambda: self.antenna_plot(antenna_stats),
                ),
                (
                    "spw_flag_fraction",
                    lambda: self.spw_plot(spw_stats),
                ),
                (
                    "field_flag_fraction",
                    lambda: self.field_plot(field_stats),
                ),
            ]:
                try:
                    outfile = operation()
                except Exception as exc:
                    status["warnings"].append(
                        f"Could not create {label}: {exc}"
                    )
                    continue

                if outfile is None:
                    status["warnings"].append(
                        f"{label} has no usable input statistics."
                    )
                else:
                    status["plots"].append(str(outfile))

        self.save_json(
            status,
            self.root / "qa_plot_status.json",
        )

        if status["plots"]:
            print(
                "[flagging] QA plots written: " +
                ", ".join(status["plots"]),
                flush=True,
            )
        else:
            print(
                "WARNING: No flagging QA plot was created. See "
                f"{self.root / 'qa_plot_status.json'} for the reason.",
                flush=True,
            )

        return status

    def regenerate_qa_plots(self):
        """Regenerate QA figures from saved reports without changing the MS."""

        history_record = self.read_json(
            self.root / "iteration_history.json",
            {},
        )
        history = (
            history_record.get("history", [])
            if isinstance(history_record, dict)
            else []
        )

        return self.generate_qa_plots(
            self.read_json(self.root / "preflag_summary.json", {}),
            history,
            self.read_json(self.root / "postflag_summary.json", {}),
            self.read_json(self.root / "antenna_statistics.json", {}),
            self.read_json(self.root / "spw_statistics.json", {}),
            self.read_json(self.root / "field_statistics.json", {}),
        )

    # --------------------------------------------------
    # Iterative flagging
    # --------------------------------------------------

    def iterative_cycles(self):

        baseline_file = (
            self.root /
            "iterative_baseline_summary.json"
        )
        baseline = self.read_json(
            baseline_file,
            {},
        )

        if not baseline:
            baseline = self.run_step(
                "Auto-flag baseline summary",
                self.summary,
            )
            self.save_json(
                baseline,
                baseline_file,
            )

        prev = self.summary_flagged_fraction(
            baseline
        )
        previous_summary = baseline

        if prev is None:
            raise RuntimeError(
                "Cannot run guarded auto-flagging because the baseline "
                "flag fraction is unavailable."
            )

        history = []

        for i in range(self.max_iterations):

            number = i + 1
            label = f"Iteration {number}/{self.max_iterations}"
            summary_file = self.root / f"iter_{number}_summary.json"
            summary = self.read_json(summary_file, {})
            thresholds = self.iteration_thresholds(number)

            # A saved summary is durable evidence that both algorithms and
            # the CASA flagmanager checkpoint completed in an older run.
            if summary:

                self.mark_step_done(f"iteration_{number}")
                print(
                    f"[flagging] {label} already completed; loading "
                    "its saved summary.",
                    flush=True,
                )

            else:
                before_step = f"iteration_{number}_before"
                before_version = f"before_iter_{number}"

                # If a previous process stopped between saving the transaction
                # checkpoint and writing the accepted/rejected summary, its
                # FLAG column may contain a partial tfcrop/rflag pass. Restore
                # the exact pre-pass flags and repeat the full pass.
                if self.step_done(before_step):
                    print(
                        f"[flagging] Recovering incomplete {label} from "
                        f"{before_version}.",
                        flush=True,
                    )
                    self.run_step(
                        f"{label}: restore pre-pass flags",
                        self.restore_checkpoint,
                        before_version,
                    )

                    for stale_step in [
                        f"iteration_{number}_tfcrop",
                        f"iteration_{number}_rflag",
                        f"iteration_{number}_flagmanager",
                    ]:
                        (
                            self.checkdir /
                            f"{stale_step}.done"
                        ).unlink(missing_ok=True)

                else:
                    self.run_checkpointed_step(
                        before_step,
                        f"{label}: save pre-pass flags",
                        self.checkpoint,
                        before_version,
                    )

                print(
                    "[flagging] "
                    f"{label} thresholds: "
                    f"tfcrop(time={thresholds['tfcrop_timecutoff']:.2f}, "
                    f"freq={thresholds['tfcrop_freqcutoff']:.2f}); "
                    f"rflag(time={thresholds['rflag_timedevscale']:.2f}, "
                    f"freq={thresholds['rflag_freqdevscale']:.2f}).",
                    flush=True,
                )

                self.run_checkpointed_step(
                    f"iteration_{number}_tfcrop",
                    f"{label}: tfcrop",
                    self.tfcrop,
                    number,
                )
                self.run_checkpointed_step(
                    f"iteration_{number}_rflag",
                    f"{label}: rflag",
                    self.rflag,
                    number,
                )

                attempted = self.run_step(
                    f"{label}: flag summary",
                    self.summary,
                )
                attempted_fraction = self.summary_flagged_fraction(
                    attempted
                )

                if attempted_fraction is None:
                    self.run_step(
                        f"{label}: restore after invalid summary",
                        self.restore_checkpoint,
                        before_version,
                    )
                    raise RuntimeError(
                        f"{label} produced no valid flag summary; "
                        "pre-pass flags were restored."
                    )

                attempted_delta = max(
                    0.0,
                    attempted_fraction - prev,
                )
                violations = []

                if (
                    attempted_delta >
                    self.max_iteration_new_fraction
                ):
                    violations.append(
                        "new flag fraction "
                        f"{attempted_delta:.4f} exceeds "
                        f"{self.max_iteration_new_fraction:.4f}"
                    )

                if (
                    attempted_fraction >
                    self.max_total_flagged_fraction
                ):
                    violations.append(
                        "cumulative flag fraction "
                        f"{attempted_fraction:.4f} exceeds "
                        f"{self.max_total_flagged_fraction:.4f}"
                    )

                for category in ("field", "spw"):
                    before_groups = self.summary_group_fractions(
                        previous_summary,
                        category,
                    )
                    after_groups = self.summary_group_fractions(
                        attempted,
                        category,
                    )

                    for group, after_fraction in after_groups.items():
                        before_fraction = before_groups.get(
                            group,
                            0.0,
                        )
                        group_delta = max(
                            0.0,
                            after_fraction - before_fraction,
                        )

                        if (
                            group_delta >
                            self.max_group_iteration_new_fraction
                        ):
                            violations.append(
                                f"{category} {group} added "
                                f"{group_delta:.4f}, exceeding "
                                f"{self.max_group_iteration_new_fraction:.4f}"
                            )

                        if (
                            after_fraction >
                            self.max_group_total_flagged_fraction
                        ):
                            violations.append(
                                f"{category} {group} total "
                                f"{after_fraction:.4f} exceeds "
                                f"{self.max_group_total_flagged_fraction:.4f}"
                            )

                if violations:
                    reason = "; ".join(violations)
                    print(
                        f"WARNING: Rejecting {label}: {reason}.",
                        flush=True,
                    )
                    self.run_step(
                        f"{label}: restore rejected flags",
                        self.restore_checkpoint,
                        before_version,
                    )
                    summary = self.run_step(
                        f"{label}: restored flag summary",
                        self.summary,
                    )
                    restored_fraction = (
                        self.summary_flagged_fraction(
                            summary
                        )
                    )
                    summary["_pipeline"] = {
                        "accepted": False,
                        "rejection_reason": reason,
                        "thresholds": thresholds,
                        "ntime": self.flagdata_ntime(),
                        "before_fraction": prev,
                        "attempted_fraction": attempted_fraction,
                        "attempted_new_fraction": attempted_delta,
                        "restored_fraction": restored_fraction,
                    }
                    self.save_json(
                        summary,
                        summary_file,
                    )
                    self.mark_step_done(
                        f"iteration_{number}"
                    )
                else:
                    self.run_checkpointed_step(
                        f"iteration_{number}_flagmanager",
                        f"{label}: accepted checkpoint",
                        self.checkpoint,
                        f"iter_{number}",
                    )
                    summary = attempted
                    summary["_pipeline"] = {
                        "accepted": True,
                        "thresholds": thresholds,
                        "ntime": self.flagdata_ntime(),
                        "before_fraction": prev,
                        "attempted_fraction": attempted_fraction,
                        "attempted_new_fraction": attempted_delta,
                    }
                    self.save_json(
                        summary,
                        summary_file,
                    )
                    self.mark_step_done(
                        f"iteration_{number}"
                    )

            metadata = summary.get(
                "_pipeline",
                {},
            )
            accepted = bool(
                metadata.get(
                    "accepted",
                    True,
                )
            )
            frac = self.summary_flagged_fraction(
                summary
            )

            if frac is None:
                raise RuntimeError(
                    f"{label} summary contains no valid flag fraction."
                )

            delta = (
                float(
                    metadata.get(
                        "attempted_new_fraction",
                        max(0.0, frac - prev),
                    )
                )
                if accepted
                else 0.0
            )
            record = {
                "iteration": number,
                "accepted": accepted,
                "cumulative_fraction": frac,
                "new_fraction": delta,
                "thresholds": metadata.get(
                    "thresholds",
                    thresholds,
                ),
                # Iterations completed before chunk planning existed ran on
                # whole scans.
                "ntime": metadata.get("ntime", "scan"),
            }

            for key in [
                "rejection_reason",
                "attempted_fraction",
                "attempted_new_fraction",
                "restored_fraction",
            ]:
                if key in metadata:
                    record[key] = metadata[key]

            history.append(record)

            print(
                f"[flagging] {label}: "
                f"cumulative={frac:.4f}, "
                f"new={delta:.4f}, "
                f"accepted={accepted}",
                flush=True,
            )

            if not accepted:
                print(
                    "[flagging] Auto-flagging stopped after rollback. "
                    "Review the saved attempted fractions before changing "
                    "the safety limits.",
                    flush=True,
                )
                break

            #
            # Stop when newly flagged data become negligible
            #

            if delta <= self.convergence_new_fraction:

                print(
                    "[flagging] Convergence reached.",
                    flush=True,
                )

                break

            prev = frac
            previous_summary = summary

        return history

    # --------------------------------------------------

    def guarded_extension(self):
        """Extend flags once, with automatic rollback on a large jump."""

        report_path = (
            self.root /
            "extension_report.json"
        )
        existing = self.read_json(
            report_path,
            {},
        )

        if existing and self.step_done("extend"):
            print(
                "[flagging] Guarded flag extension already completed.",
                flush=True,
            )
            return existing

        if self.step_done("extend") and not existing:
            migrated_summary = self.run_step(
                "Migrated extension flag summary",
                self.summary,
            )
            migrated_fraction = self.summary_flagged_fraction(
                migrated_summary
            )
            report = {
                "accepted": None,
                "rejection_reason": None,
                "migrated_from_legacy_checkpoint": True,
                "growtime": None,
                "growfreq": None,
                "extendpols": None,
                "before_fraction": None,
                "attempted_fraction": None,
                "attempted_new_fraction": None,
                "final_fraction": migrated_fraction,
            }
            self.save_json(
                report,
                report_path,
            )
            print(
                "WARNING: Reused a legacy extension checkpoint; its "
                "original growth parameters cannot be reconstructed.",
                flush=True,
            )
            return report

        before_version = "before_extend"
        before_step = "extend_before"

        if self.step_done(before_step) and not existing:
            print(
                "[flagging] Recovering an incomplete flag-extension "
                "transaction.",
                flush=True,
            )
            self.run_step(
                "Restore pre-extension flags",
                self.restore_checkpoint,
                before_version,
            )
            (
                self.checkdir /
                "extend_operation.done"
            ).unlink(missing_ok=True)
        else:
            self.run_checkpointed_step(
                before_step,
                "Pre-extension checkpoint",
                self.checkpoint,
                before_version,
            )

        before = self.run_step(
            "Pre-extension flag summary",
            self.summary,
        )
        before_fraction = self.summary_flagged_fraction(
            before
        )

        if before_fraction is None:
            raise RuntimeError(
                "Cannot guard flag extension because the pre-extension "
                "flag fraction is unavailable."
            )

        self.run_checkpointed_step(
            "extend_operation",
            "Flag extension",
            self.extend,
        )
        attempted = self.run_step(
            "Post-extension flag summary",
            self.summary,
        )
        attempted_fraction = self.summary_flagged_fraction(
            attempted
        )

        if attempted_fraction is None:
            self.restore_checkpoint(before_version)
            raise RuntimeError(
                "Flag extension produced no valid summary; pre-extension "
                "flags were restored."
            )

        delta = max(
            0.0,
            attempted_fraction - before_fraction,
        )
        violations = []

        if delta > self.max_extension_new_fraction:
            violations.append(
                f"extension added {delta:.4f}, exceeding "
                f"{self.max_extension_new_fraction:.4f}"
            )

        if attempted_fraction > self.max_total_flagged_fraction:
            violations.append(
                f"total fraction {attempted_fraction:.4f} exceeds "
                f"{self.max_total_flagged_fraction:.4f}"
            )

        for category in ("field", "spw"):
            before_groups = self.summary_group_fractions(
                before,
                category,
            )
            after_groups = self.summary_group_fractions(
                attempted,
                category,
            )

            for group, after_fraction in after_groups.items():
                before_group_fraction = before_groups.get(
                    group,
                    0.0,
                )
                group_delta = max(
                    0.0,
                    after_fraction - before_group_fraction,
                )

                if (
                    group_delta >
                    self.max_group_extension_new_fraction
                ):
                    violations.append(
                        f"{category} {group} extension added "
                        f"{group_delta:.4f}, exceeding "
                        f"{self.max_group_extension_new_fraction:.4f}"
                    )

                if (
                    after_fraction >
                    self.max_group_total_flagged_fraction
                ):
                    violations.append(
                        f"{category} {group} total "
                        f"{after_fraction:.4f} exceeds "
                        f"{self.max_group_total_flagged_fraction:.4f}"
                    )

        if violations:
            reason = "; ".join(violations)
            print(
                "WARNING: Rejecting flag extension and restoring its "
                f"checkpoint: {reason}.",
                flush=True,
            )
            self.run_step(
                "Restore rejected flag extension",
                self.restore_checkpoint,
                before_version,
            )
            final = self.run_step(
                "Restored extension summary",
                self.summary,
            )
            accepted = False
        else:
            final = attempted
            accepted = True
            reason = None

        report = {
            "accepted": accepted,
            "rejection_reason": reason,
            "growtime": self.extend_growtime,
            "growfreq": self.extend_growfreq,
            # growtime is a fraction of the time chunk, not of the scan, once
            # the chunk is bounded.
            "ntime": self.flagdata_ntime(),
            "extendpols": self.extend_polarizations,
            "before_fraction": before_fraction,
            "attempted_fraction": attempted_fraction,
            "attempted_new_fraction": delta,
            "final_fraction": self.summary_flagged_fraction(
                final
            ),
        }
        self.save_json(
            report,
            report_path,
        )
        self.mark_step_done("extend")
        return report

    # --------------------------------------------------

    def run(self):

        if self.done():

            post_saved = self.read_json(
                self.root / "postflag_summary.json",
                {},
            )

            if post_saved:

                expected_ms = str(
                    Path(self.ms).expanduser().resolve()
                )
                reported_ms = post_saved.get("measurement_set")
                if reported_ms and str(
                    Path(reported_ms).expanduser().resolve()
                ) != expected_ms:
                    raise RuntimeError(
                        "flagging.done belongs to a different Measurement "
                        f"Set: {reported_ms} != {expected_ms}"
                    )
                reported_telescope = str(
                    post_saved.get("telescope") or self.telescope
                ).strip().lower()
                if reported_telescope != str(self.telescope).strip().lower():
                    raise RuntimeError(
                        "flagging.done belongs to a different telescope: "
                        f"{reported_telescope} != {self.telescope}"
                    )
                reported_band = str(
                    post_saved.get("band") or self.band
                ).strip().upper()
                if reported_band != str(self.band).strip().upper():
                    raise RuntimeError(
                        "flagging.done belongs to a different receiver band: "
                        f"{reported_band} != {self.band}"
                    )

                previous_fingerprint = post_saved.get(
                    "provenance_fingerprint"
                )
                if previous_fingerprint != self.provenance_fingerprint:
                    history = post_saved.setdefault(
                        "provenance_migration_history",
                        [],
                    )
                    history.append({
                        "updated_unix_time": time.time(),
                        "previous_fingerprint": previous_fingerprint,
                        "current_fingerprint": self.provenance_fingerprint,
                        "reason": (
                            "completed flagging checkpoint reused after "
                            "the orchestrator validated an unchanged input "
                            "and science command"
                        ),
                    })

                # A completed flagging run is allowed to refresh only its
                # orchestration binding. No flagdata task is repeated here.
                post_saved["telescope"] = self.telescope
                post_saved["band"] = self.band
                post_saved["measurement_set"] = expected_ms
                post_saved["provenance_fingerprint"] = (
                    self.provenance_fingerprint
                )
                self.save_json(
                    post_saved,
                    self.root / "postflag_summary.json",
                )

                print(
                    "Flagging already completed; refreshed its validated "
                    "provenance and regenerating QA plots from saved "
                    "summaries."
                )
                self.regenerate_qa_plots()
                return

            print(
                "WARNING: flagging.done exists but the post-flag report is "
                "missing or invalid; rebuilding reports from checkpoints."
            )

        pre_file = self.root / "preflag_summary.json"
        pre = self.read_json(pre_file, {})

        if not pre:

            # Nothing has been saved yet in this run, so a version of ours
            # on the MS can only have come from an earlier one.
            self.refuse_previously_flagged_input()

            pre = self.run_step(
                "Pre-flag summary",
                self.summary,
            )
            self.save_json(pre, pre_file)

        # Migration support: an iteration summary from an older version is
        # proof that all initial deterministic flagging steps completed.
        if any(self.root.glob("iter_*_summary.json")):

            for name in [
                "before_flagging",
                "autocorrelations",
                "shadow",
                "quack",
                "clip",
                "known_frequency_masks",
                "edge_channels",
                "manual_commands",
            ]:

                self.mark_step_done(name)

        self.run_checkpointed_step(
            "before_flagging",
            "Pre-flag checkpoint",
            self.checkpoint,
            "before_flagging",
        )

        band_profile = TELESCOPE_PROFILE_MODULE.get_band_profile(
            self.telescope,
            self.band,
        )
        guaranteed_range = None
        if (
            self.telescope == "meerkat" and
            "guaranteed_freq_min" in band_profile and
            "guaranteed_freq_max" in band_profile
        ):
            guaranteed_range = [
                float(band_profile["guaranteed_freq_min"]) / 1e6,
                float(band_profile["guaranteed_freq_max"]) / 1e6,
            ]

        self.save_json(
            {
                "execution_profile": self.profile,
                "telescope": self.telescope,
                "band": self.band,
                "parallel_hand_correlation": (
                    self.parallel_hand_correlation()
                ),
                "max_iterations": self.max_iterations,
                "edge_fraction": self.edge_fraction,
                "edge_flagging_scheme": (
                    "outside_guaranteed_science_range"
                    if guaranteed_range is not None else
                    "percentage_per_spw"
                ),
                "guaranteed_science_range_MHz": guaranteed_range,
                "known_frequency_masks": self.known_frequency_masks,
                "manual_commands": self.manual_flag_commands,
                "precal_rflag": self.precal_rflag_enabled,
                "quack_interval_s": self.quack_interval_s,
                "flag_autocorrelations": (
                    self.flag_autocorrelations_enabled
                ),
                "flag_shadowed": self.flag_shadowed_enabled,
                "shadow_tolerance_m": self.shadow_tolerance_m,
                "autoflag": self.autoflag_config,
                "internal_autoflag_extension": False,
                "flagdata_chunk_plan": self.flagdata_chunk_plan(),
                "extend_growtime": self.extend_growtime,
                "extend_growfreq": self.extend_growfreq,
                "extend_polarizations": self.extend_polarizations,
                "convergence_new_fraction": (
                    self.convergence_new_fraction
                ),
                "max_iteration_new_fraction": (
                    self.max_iteration_new_fraction
                ),
                "max_total_flagged_fraction": (
                    self.max_total_flagged_fraction
                ),
                "max_group_iteration_new_fraction": (
                    self.max_group_iteration_new_fraction
                ),
                "max_group_total_flagged_fraction": (
                    self.max_group_total_flagged_fraction
                ),
                "max_extension_new_fraction": (
                    self.max_extension_new_fraction
                ),
                "max_group_extension_new_fraction": (
                    self.max_group_extension_new_fraction
                ),
                "ml_rfi": self.ml_config,
            },
            self.root / "flagging_parameters.json",
        )

        self.run_checkpointed_step(
            "autocorrelations",
            "Autocorrelation flagging",
            self.flag_autocorrelations,
        )
        self.run_checkpointed_step(
            "shadow",
            "Shadowed-data flagging",
            self.flag_shadowed_data,
        )
        self.run_checkpointed_step(
            "quack",
            "Quack flagging",
            self.quack,
        )
        self.run_checkpointed_step(
            "clip",
            "Zero-value clipping",
            self.clip,
        )
        self.run_checkpointed_step(
            "known_frequency_masks",
            "Telescope-specific frequency-mask flagging",
            self.flag_known_frequency_ranges,
        )
        self.run_checkpointed_step(
            "edge_channels",
            "Edge-channel flagging",
            self.flag_edge_channels,
        )
        self.run_checkpointed_step(
            "manual_commands",
            "Manual flag commands",
            self.flag_manual_commands,
            self.summary_counts(pre)[1],
        )

        history = self.iterative_cycles()

        self.run_checkpointed_step(
            "ml_rfi",
            "ML-assisted persistent-channel RFI flagging",
            self.ml_rfi_flagging,
        )

        extension_report = self.guarded_extension()

        antenna_file = self.root / "antenna_statistics.json"
        ant = self.read_json(antenna_file, {})

        if not ant:

            ant = self.run_step(
                "Antenna statistics",
                self.antenna_stats,
            )

        spw_file = self.root / "spw_statistics.json"
        spw = self.read_json(spw_file, None)

        if spw is None:

            spw = self.run_step(
                "Spectral-window statistics",
                self.spw_stats,
            )

        baseline_file = self.root / "baseline_statistics.json"
        base = self.read_json(baseline_file, None)

        if base is None:

            base = self.run_step(
                "Baseline statistics",
                self.baseline_stats,
            )

        field_file = self.root / "field_statistics.json"
        fields = self.read_json(field_file, None)

        if fields is None:

            fields = self.run_step(
                "Field statistics",
                self.field_stats,
            )

        bad_ant, median_total = (
            self.detect_bad_antennas(
                ant
            )
        )

        refant = self.recommend_refant(
            ant,
            bad_ant
        )

        self.save_json(
            ant,
            self.root /
            "antenna_statistics.json"
        )

        self.save_json(
            {
                "median_total": median_total,
                "bad_antennas": bad_ant,
                "recommended_refant": refant,
                "reference_antenna_selection": getattr(
                    self,
                    "refant_selection",
                    {},
                ),
                "reference_antenna_ranking": getattr(
                    self,
                    "refant_ranking",
                    [],
                )[:10],
            },
            self.root / "bad_antenna_report.json"
        )

        self.save_json(
            spw,
            self.root /
            "spw_statistics.json"
        )

        self.save_json(
            base,
            self.root /
            "baseline_statistics.json"
        )

        self.save_json(
            fields,
            self.root /
            "field_statistics.json"
        )

        with open(
            self.root /
            "bad_antennas.txt",
            "w"
        ) as f:

            f.write(
                "\n".join(bad_ant)
            )

        with open(
            self.root /
            "recommended_refant.txt",
            "w"
        ) as f:

            f.write(refant)

        post_file = self.root / "postflag_summary.json"
        post = self.run_step(
            "Post-flag summary",
            self.summary,
        )
        # CASA's flagdata summary is backend-neutral and contains no
        # instrument identity of its own.  Add immutable context so resume
        # validation cannot adopt a MeerKAT/LOFAR flag set under another
        # telescope or receiver profile.
        post["telescope"] = self.telescope
        post["band"] = self.band
        post["measurement_set"] = str(
            Path(self.ms).expanduser().resolve()
        )
        post["provenance_fingerprint"] = self.provenance_fingerprint
        self.save_json(post, post_file)

        self.save_json(
            {
                "execution_profile": self.profile,
                "max_iterations": self.max_iterations,
                "baseline_fraction": (
                    self.summary_flagged_fraction(
                        self.read_json(
                            self.root /
                            "iterative_baseline_summary.json",
                            {},
                        )
                    )
                ),
                "history": history,
                "extension": extension_report,
            },
            self.root /
            "iteration_history.json"
        )

        self.generate_qa_plots(
            pre,
            history,
            post,
            ant,
            spw,
            fields,
        )

        attribution = self.flagging_attribution(pre, post)
        if attribution:
            self.save_json(attribution, self.root / "flagging_attribution.json")
            print(
                f"Flagging: {attribution['input_flagged_fraction'] * 100:.1f}% "
                "of the visibilities arrived already flagged; this run "
                "flagged "
                f"{attribution['this_run_fraction_of_received'] * 100:.1f}% "
                "of the data it received unflagged, for "
                f"{attribution['final_flagged_fraction'] * 100:.1f}% in total.",
                flush=True,
            )

        self.run_checkpointed_step(
            "after_flagging",
            "Post-flag checkpoint",
            self.checkpoint,
            "after_flagging",
        )

        self.mark_done()

        print(
            "Flagging completed.",
            flush=True,
        )


def build_parser():

    parser = argparse.ArgumentParser(
        description="Run profile-aware CASA flagging."
    )

    parser.add_argument(
        "ms",
        help="Input Measurement Set.",
    )

    parser.add_argument(
        "--outdir",
        default=RUN_CONFIG.get(
            "outdir",
            "Products",
        ),
        help="Pipeline product directory.",
    )

    parser.add_argument(
        "--profile",
        choices=getattr(
            CONFIG,
            "EXECUTION_PROFILE_NAMES",
            ("science",),
        ),
        default=RUN_CONFIG.get("profile", "science"),
        help=(
            "Execution profile. 'laptop' caps iterative RFI flagging at "
            "three passes while retaining QA statistics."
        ),
    )

    parser.add_argument(
        "--max-iterations",
        type=int,
        default=None,
        help="Override iterative flagging passes.",
    )

    parser.add_argument(
        "--bad-ant-threshold",
        type=float,
        default=None,
        help="Override bad antenna flag-fraction threshold.",
    )

    parser.add_argument(
        "--edge-fraction",
        type=float,
        default=None,
        help="Override band edge fraction to flag.",
    )

    parser.add_argument(
        "--provenance-fingerprint",
        default=None,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--plots-only",
        action="store_true",
        help=(
            "Regenerate flagging QA plots from saved JSON reports without "
            "modifying the Measurement Set."
        ),
    )

    return parser



# An interrupt reaches a batch stage when the job is cancelled, when it runs
# into its time limit, or when an interactive session ends. CASA would print
# a bare KeyboardInterrupt traceback from wherever it happened to be, which
# reads like a crash. Exit code 130 is the conventional one for SIGINT, and
# the orchestrator reports it as an interruption rather than a failure.
def _interrupted(stage):
    print(
        f"\n{stage} was interrupted from outside the pipeline (SIGINT): the "
        "job was cancelled, reached its time limit, or its session ended. "
        "Completed checkpoints are kept; resubmitting continues from them.",
        flush=True,
    )
    sys.exit(130)

if __name__ == "__main__":

    args = build_parser().parse_args()

    kwargs = {}

    if args.max_iterations is not None:

        kwargs["max_iterations"] = args.max_iterations

    if args.bad_ant_threshold is not None:

        kwargs["bad_ant_threshold"] = args.bad_ant_threshold

    if args.edge_fraction is not None:

        kwargs["edge_fraction"] = args.edge_fraction

    flagger = CASAFlagger(
        args.ms,
        outdir=args.outdir,
        profile=args.profile,
        provenance_fingerprint=args.provenance_fingerprint,
        **kwargs,
    )

    try:
        if args.plots_only:
            flagger.regenerate_qa_plots()
        else:
            flagger.run()
    except KeyboardInterrupt:
        _interrupted("Flagging")
