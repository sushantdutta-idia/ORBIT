#!/usr/bin/env python3
#
# ============================================================
# selfcal.py
#
# CASA-native science-grade uGMRT/MeerKAT self-calibration pipeline
#
# LOFAR self-calibration is deliberately handled by external LOFAR workflows.
#
# Compatible with CASA 6+
#
# Author: Sushant Dutta
# ============================================================

import importlib.util
import os
import sys
import time
import json
import shutil
import math
import numpy as np
from pathlib import Path


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of failing partway through a
# stage. Increase it whenever pipeline.py starts to rely on something new here.
PIPELINE_INTERFACE_VERSION = 12


SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

if (
    sys.path and
    os.path.abspath(sys.path[0]) == SCRIPT_DIR
):

    sys.path.pop(0)

import argparse


def significance_floor(configured, diagnostics):
    """The larger of a configured sigma level and what noise alone reaches.

    Peaks are the true image extremes, so a fixed 5-sigma trigger would fire
    on noise in any large image. The noise-only level for the image's size
    (see casa_compat.noise_peak_snr) is the least a trigger may use.
    """

    floor = float(configured)
    noise = (diagnostics or {}).get("noise_only_peak_snr")
    return max(floor, float(noise)) if noise else floor


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
applycal = CASA_COMPAT.import_casa_task("applycal")
flagdata = CASA_COMPAT.import_casa_task("flagdata")
flagmanager = CASA_COMPAT.import_casa_task("flagmanager")
gaincal = CASA_COMPAT.import_casa_task("gaincal")
imhead = CASA_COMPAT.import_casa_task("imhead")
rmtables = CASA_COMPAT.import_casa_task("rmtables")
imstat = CASA_COMPAT.import_casa_task("imstat")
# plotms is deliberately disabled by default.  The ilifu CASA 6.7 build
# may expose a partial/broken plotting installation, so merely importing or
# invoking plotms can terminate an otherwise valid reduction.  QA products
# are generated directly from CASA tables by qa_products.py instead.
plotms = None
split = CASA_COMPAT.import_casa_task("split")
tclean = CASA_COMPAT.import_casa_task("tclean")
from casatools import table


tb = table()


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


def load_qa_products_module():

    module_path = Path(SCRIPT_DIR) / "qa_products.py"
    spec = importlib.util.spec_from_file_location(
        "ugmrt_qa_products",
        str(module_path),
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
        "radio_continuum_selfcal_identity",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# MODEL_DATA checks read at most this many values, in blocks of at most
# MODEL_CHECK_BLOCK_VALUES, spread across the whole Measurement Set.
MODEL_CHECK_MAX_VALUES = 1 << 26
MODEL_CHECK_BLOCK_VALUES = 1 << 23


CONFIG = load_config_module()
RUN_CONFIG = (
    getattr(CONFIG, "PIPELINE_RUN_CONFIG", {})
    if CONFIG is not None
    else {}
)
SELFCAL_TASK_DEFAULTS = (
    getattr(CONFIG, "SELFCAL_TASK_CONFIG", {})
    if CONFIG is not None
    else {}
)
QA_PRODUCTS = load_qa_products_module()
TELESCOPE_PROFILE_MODULE = load_telescope_profiles_module()
IMAGING_SAFETY = load_imaging_safety_module()


def resolve_execution_profile(profile):
    """Resolve a profile without mutating the shared configuration."""

    name = str(profile or RUN_CONFIG.get("profile", "science"))
    name = name.strip().lower()

    if CONFIG is not None and hasattr(CONFIG, "get_execution_profile"):
        return name, CONFIG.get_execution_profile(name)

    if name == "science":
        return name, {}

    raise ValueError(
        "Execution profiles require config.py; use the 'science' profile "
        "when running selfcal.py as a standalone script without it."
    )

# ============================================================
# BAND CONFIGURATION
# ============================================================

SELFCAL_CONFIG = {

    "B2": {
        "phase_solints": ["8min", "4min", "2min", "int"],
        "amp_solint": "inf",
        "min_snr_phase": 5.0,
        "min_snr_amp": 30.0,
    },

    "B3": {
        "phase_solints": ["8min", "4min", "2min", "int"],
        "amp_solint": "inf",
        "min_snr_phase": 5.0,
        "min_snr_amp": 30.0,
    },

    "B4": {
        "phase_solints": ["4min", "2min", "1min", "int"],
        "amp_solint": "inf",
        "min_snr_phase": 5.0,
        "min_snr_amp": 40.0,
    },

    "B5": {
        "phase_solints": ["4min", "2min", "1min"],
        "amp_solint": "scan",
        "min_snr_phase": 5.0,
        "min_snr_amp": 50.0,
    }
}

IMAGING_CONFIG = {

    "B2": {
        "cell": "1.0arcsec",
        "imsize": 8192,
        "robust": 0.5,
        "nterms": 2,
        "niter": 100000,
        "nsigma": 3.0,
        "wprojplanes": 128,
    },

    "B3": {
        "cell": "0.8arcsec",
        "imsize": 8192,
        "robust": 0.5,
        "nterms": 2,
        "niter": 100000,
        "nsigma": 3.0,
        "wprojplanes": 128,
    },

    "B4": {
        "cell": "0.5arcsec",
        "imsize": 8192,
        "robust": 0.5,
        "nterms": 2,
        "niter": 100000,
        "nsigma": 3.0,
        "wprojplanes": 128,
    },

    "B5": {
        "cell": "0.3arcsec",
        "imsize": 8192,
        "robust": 0.5,
        "nterms": 2,
        "niter": 100000,
        "nsigma": 3.0,
        "wprojplanes": 128,
    }
}


if CONFIG is not None:

    SELFCAL_CONFIG = getattr(
        CONFIG,
        "SELFCAL_CONFIG",
        SELFCAL_CONFIG,
    )
    IMAGING_CONFIG = getattr(
        CONFIG,
        "SELFCAL_IMAGING_CONFIG",
        IMAGING_CONFIG,
    )


# ============================================================
# SELF-CAL MANAGER
# ============================================================

# Most major cycles a masked (auto-multithresh) self-calibration pass may run;
# see SelfCalManager.cap_masked_major_cycles. Those of the uGMRT Band 4 run
# that set it converged in 13, 19 and 21 cycles, and its fourth was still
# oscillating after 22; 20 would have stopped only the 21-cycle pass, one
# cycle early. It lives here rather than in config.py so that a
# self-calibration already under way keeps its recipe, and with it its
# accepted loops, when this file is updated.
DEFAULT_CLEAN_MAX_MAJOR_CYCLES = 20

# The amp+phase self-calibration solve. It used the band's ``min_snr_amp``
# (40 for uGMRT Band 4) as gaincal's minsnr, one solution per polarization,
# and a single reference antenna. On the first field it was tried on, a
# 0.7 Jy field whose model visibility is about 8 mJy rms against 0.9 mJy of
# noise per baseline, scan and polarization, an antenna's amplitude solution
# has a signal-to-noise near 40: gaincal flagged all but eight of 81
# solutions, the reference antenna's among them in two of three scans, and
# with no alternate reference CASA then flagged every antenna in those
# scans. Solving the two polarizations together gains a factor of root two;
# a solution good to 10% is kept, and whether the set of them helps is
# decided on the image, not by this threshold. Like the major-cycle limit
# these default here so that a self-calibration under way keeps its recipe.
DEFAULT_AMP_SELFCAL_MINSNR = 10.0
DEFAULT_AMP_SELFCAL_GAINTYPE = "T"

# How applycal carries an amp+phase solution to the data. The solve holds each
# gain constant over its interval, which for the usual ``inf`` is a whole
# scan, so each visibility takes the solution of its own interval. Linear
# interpolation instead joins one scan's solution to the next scan's, hours
# away, and leaves uncorrected every visibility that lies between a good
# solution and a failed one: with a failed last scan that is the second half
# of the scan before it, whose gain would then step from the solved value to
# one at its midpoint.
AMP_SELFCAL_INTERP = "nearest"


class SelfCalManager:

    # --------------------------------------------------------
    # Constructor
    # --------------------------------------------------------

    def __init__(self,
                 ms,
                 outdir="Products",
                 upstream_outdir=None,
                 profile=None,
                 telescope=None,
                 receiver_band=None,
                 reference_antenna=None,
                 target=None,
                 provenance_fingerprint=None):

        self.ms = ms
        self.input_ms = str(Path(ms).expanduser().resolve())
        self.requested_telescope = telescope
        self.requested_band = receiver_band
        self.reference_antenna_override = reference_antenna
        self.target = target or Path(str(ms)).stem
        self.provenance_fingerprint = (
            str(provenance_fingerprint).strip()
            if provenance_fingerprint else None
        )
        self.profile, self.profile_settings = (
            resolve_execution_profile(profile)
        )

        self.task_config = {
            "minblperant": 4,
            "gain_minsnr": 3.0,
            "phase_improvement_stop": 0.02,
            "applymode": "calonly",
            "phase_gaintype": "G",
            "amp_gaintype": "G",
            "phase_refantmode": "flex",
            "phase_solmode": "",
            "phase_start_dynamic_range": 20.0,
            "phase_start_min_solution_snr": 10.0,
            "amp_start_dynamic_range": 100.0,
            "phase_retry_enabled": True,
            "phase_retry_solint": "inf",
            "phase_retry_combine": "scan",
            "phase_retry_minsnr": 3.0,
            "phase_retry_minblperant": 3,
            "phase_retry_gaintype": "T",
            "model_check_rows": 20000,
            "model_nonzero_epsilon": 1.0e-12,
            "max_caltable_unusable_fraction": 0.50,
            "max_caltable_group_unusable_fraction": 0.80,
            "calwt": False,
            "parang": False,
        }
        self.task_config.update(
            SELFCAL_TASK_DEFAULTS
        )

        # CASA stages are often launched from a Slurm submission directory,
        # rather than the pipeline directory.  Keep every self-calibration
        # product path absolute so gaincal() and applycal() always refer to
        # the same calibration table.
        self.outdir = Path(outdir).expanduser().resolve()
        # Inspection and calibration describe the observation, so with one
        # product directory per target they stay at the run's root.
        self.upstream_outdir = (
            Path(upstream_outdir).expanduser().resolve()
            if upstream_outdir else self.outdir
        )

        self.selfcaldir = (
            self.outdir /
            "selfcal"
        )

        self.plotdir = (
            self.selfcaldir /
            "plots"
        )

        self.imagedir = (
            self.selfcaldir /
            "images"
        )

        self.caltabledir = (
            self.selfcaldir /
            "caltables"
        )

        self.caltabledir.mkdir(
            parents=True,
            exist_ok=True
        )

        self.checkdir = (
            self.selfcaldir /
            "checkpoints"
        )

        self.reportdir = (
            self.selfcaldir /
            "reports"
        )

        self.msdir = (
            self.selfcaldir /
            "selfcal_ms"
        )

        self.msdir.mkdir(
            parents=True,
            exist_ok=True
        )

        self.current_ms = self.ms
        self.clean_diagnostics = {}
        self.caltable_qa = {}

        #
        # Create directories
        #

        for d in [

            self.selfcaldir,
            self.plotdir,
            self.imagedir,
            self.checkdir,
            self.reportdir

        ]:

            d.mkdir(
                parents=True,
                exist_ok=True
            )

        #
        # Load metadata
        #

        self.load_inspection_report()

        self.load_calibration_report()

        detected_telescope = (
            TELESCOPE_PROFILE_MODULE.telescope_name_from_report(
                self.inspection
            ) or
            TELESCOPE_PROFILE_MODULE.telescope_name_from_report(
                self.calibration
            )
        )
        if self.requested_telescope not in (None, "", "auto"):
            requested = TELESCOPE_PROFILE_MODULE.normalize_telescope(
                self.requested_telescope
            )
            if detected_telescope and requested != detected_telescope:
                raise RuntimeError(
                    f"Requested telescope {requested} conflicts with the "
                    f"inspection/calibration reports ({detected_telescope})."
                )
            self.telescope = requested
        else:
            self.telescope = detected_telescope
        if self.telescope is None:
            raise RuntimeError(
                "No supported telescope identity is present in the "
                "inspection/calibration reports."
            )
        if not TELESCOPE_PROFILE_MODULE.get_telescope_profile(
            self.telescope
        )["casa_selfcal_supported"]:
            raise RuntimeError(
                "CASA self-calibration is disabled for LOFAR because it "
                "does not solve LOFAR beam/clock/TEC terms. Use a DP3/LoSoTo "
                "or facet-calibrated target MS and start at imaging."
            )

        # Determine band

        self.choose_band_config()
        self.apply_execution_profile()

        self.reference_antenna = (
            self.reference_antenna_override or
            self.get_reference_antenna()
        )

        if (
            not self.reference_antenna or
            str(self.reference_antenna).upper() == "UNKNOWN"
        ):

            raise RuntimeError(
                "No usable reference antenna is recorded in the "
                "calibration products. Rerun calibration with --refant "
                "or fix the flagging recommendation first."
            )

    # --------------------------------------------------------
    # Inspection metadata
    # --------------------------------------------------------

    def load_inspection_report(self):

        inspection_file = (
            self.upstream_outdir /
            "inspection" /
            "inspection_report.json"
        )

        if inspection_file.exists():

            with open(inspection_file) as f:

                self.inspection = json.load(f)

        else:

            print(
                "WARNING: "
                "inspection_report.json not found."
            )

            self.inspection = {}

    # --------------------------------------------------------
    # Calibration metadata
    # --------------------------------------------------------

    def load_calibration_report(self):

        calibration_file = (
            self.upstream_outdir /
            "calibration" /
            "calibration_report.json"
        )

        if calibration_file.exists():

            with open(calibration_file) as f:

                self.calibration = json.load(f)

        else:

            print(
                "WARNING: "
                "calibration_report.json not found."
            )

            self.calibration = {} 

    # --------------------------------------------------------
    # Get reference antenna
    # --------------------------------------------------------

    def get_reference_antenna(self):

        # Calibration report preferred

        if "reference_antenna" in self.calibration:

                ref = self.calibration["reference_antenna"]

                if isinstance(ref, dict):

                    return ref.get("recommended", "")

                return ref

        # Inspection report fallback

        if "reference_antenna" in self.inspection:

            ref = self.inspection["reference_antenna"]

            if isinstance(ref, dict):

                return ref.get("recommended", "")

            return ref

        return ""

    # --------------------------------------------------------
    # Band selection
    # --------------------------------------------------------

    def choose_band_config(self):

        band = None

        #
        # Inspection report preferred
        #

        if "band" in self.inspection:

            band_entry = self.inspection["band"]

            if isinstance(band_entry, dict):

                band = (
                    band_entry.get("name") or
                    band_entry.get("band") or
                    ""
                )

            else:

                band = band_entry

            band = str(band).upper()

        #
        # Calibration report fallback
        #

        elif "band" in self.calibration:

            band_entry = self.calibration.get(
                "band",
                "",
            )

            if isinstance(band_entry, dict):

                band = (
                    band_entry.get("name") or
                    band_entry.get("band") or
                    ""
                )

            else:

                band = band_entry

            band = str(band).upper()

        #
        # Default
        #

        if self.requested_band:
            requested_band = TELESCOPE_PROFILE_MODULE.normalize_band(
                self.telescope,
                self.requested_band,
            )
            if band and band in SELFCAL_CONFIG and band != requested_band:
                raise RuntimeError(
                    f"Requested band {requested_band} conflicts with the "
                    f"stage reports ({band})."
                )
            band = requested_band

        if band not in SELFCAL_CONFIG:

            raise RuntimeError(
                f"No self-calibration profile exists for "
                f"{self.telescope}/{band or 'UNKNOWN'}. Rerun inspection "
                "with the correct --telescope and --receiver-band."
            )

        self.band = band

        self.band_config = dict(SELFCAL_CONFIG[band])
        self.band_config["phase_solints"] = list(
            self.band_config["phase_solints"]
        )

        self.image_config = dict(IMAGING_CONFIG[self.band])
        self.flat_imaging_overrides = [
            item
            for item in (
                getattr(CONFIG, "FLAT_IMAGING_OVERRIDES", [])
                if CONFIG is not None
                else []
            )
            if item.get("table") == "SELFCAL_IMAGING_CONFIG" and
            item.get("band") == self.band
        ]

        for override in self.flat_imaging_overrides:
            print(
                "WARNING: config.yaml key "
                f"{override['key']}={override['config_yaml_value']!r} "
                f"replaced the {self.band} self-calibration default "
                f"{override['band_default']!r}."
            )

        # Imaging field size (degrees)
        self.field_size = (
            self.inspection
            .get("field_size_deg", 1.0)
        )

    # --------------------------------------------------------
    # Execution profile
    # --------------------------------------------------------

    def apply_execution_profile(self):
        """Apply resource-aware overrides to this manager's local settings."""

        settings = self.profile_settings.get("selfcal", {})

        for key in [
            "cell",
            "imsize",
            "robust",
            "niter",
            "nterms",
            "nsigma",
            "wprojplanes",
        ]:
            if key in settings:
                self.image_config[key] = settings[key]

        phase_loops = settings.get("phase_loops")

        if phase_loops is not None:
            try:
                phase_loops = int(phase_loops)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "Execution-profile selfcal.phase_loops must be a "
                    "positive integer."
                ) from exc

            if phase_loops < 1:
                raise ValueError(
                    "Execution-profile selfcal.phase_loops must be at least 1."
                )

            self.band_config["phase_solints"] = (
                self.band_config["phase_solints"][:phase_loops]
            )

    # --------------------------------------------------------
    # Checkpoints
    # --------------------------------------------------------

    def checkpoint_exists(self,
                          name):

        return (
            self.checkdir /
            f"{name}.done"
        ).exists()

    # --------------------------------------------------------
    # Mark checkpoint
    # --------------------------------------------------------

    def mark_checkpoint(self,
                        name):

        checkpoint = (
            self.checkdir /
            f"{name}.done"
        )

        checkpoint.touch()

    # --------------------------------------------------------
    # Remove checkpoint
    # --------------------------------------------------------

    def clear_checkpoint(self,
                         name):

        checkpoint = (
            self.checkdir /
            f"{name}.done"
        )

        if checkpoint.exists():

            checkpoint.unlink()

    # --------------------------------------------------------

    def rejection_path(self, name):

        return self.checkdir / f"{name}.rejected.json"

    # --------------------------------------------------------

    def phase_rejected(self, loop):

        return self.rejection_path(
            f"apply_phase{loop}"
        ).exists()

    # --------------------------------------------------------

    def mark_phase_rejected(self,
                            loop,
                            record):

        self.save_json(
            record,
            self.rejection_path(f"apply_phase{loop}"),
        )

    # --------------------------------------------------------

    def integration_rows(self, measurement_set):
        """Integration times of a Measurement Set and the rows at each."""

        key = str(measurement_set)
        cache = getattr(self, "_integration_rows", None)
        if cache is None:
            cache = self._integration_rows = {}
        if key in cache:
            return cache[key]

        result = None
        opened = False
        try:
            tb.open(key)
            opened = True
            if "TIME" in tb.colnames() and int(tb.nrows()) > 0:
                times, rows = np.unique(
                    np.asarray(tb.getcol("TIME"), dtype=float).reshape(-1),
                    return_counts=True,
                )
                if times.size and np.all(np.isfinite(times)):
                    result = (times, rows.astype(float))
        except Exception:
            result = None
        finally:
            if opened:
                tb.close()

        if result is not None:
            cache[key] = result
        return result

    def solution_time_shares(self, solution_times, measurement_set=None):
        """Share of the observing time that each solution time stands for.

        Solution intervals are not equal. With one solution per scan, a scan
        that flagging has cut to a minute is one interval and so is a scan of
        twenty minutes. Each integration of the Measurement Set is given to
        the solution nearest to it in time, and a solution's share is the
        fraction of the rows it receives. Returns None when that cannot be
        measured; the caller then counts intervals.
        """

        measurement_set = measurement_set or getattr(self, "current_ms", None)
        if not measurement_set:
            return None

        try:
            solutions = np.asarray(solution_times, dtype=float).reshape(-1)
        except (TypeError, ValueError):
            return None
        if (
            solutions.size < 1 or
            not np.all(np.isfinite(solutions)) or
            np.any(np.diff(solutions) <= 0)
        ):
            return None

        integrations = self.integration_rows(measurement_set)
        if integrations is None:
            return None
        times, rows = integrations

        # A Measurement Set observed at another time is not the one these
        # solutions belong to.
        if solutions[0] > times[-1] or solutions[-1] < times[0]:
            return None

        edges = 0.5 * (solutions[1:] + solutions[:-1])
        shares = np.bincount(
            np.searchsorted(edges, times),
            weights=rows,
            minlength=solutions.size,
        )
        total = float(shares.sum())
        if total <= 0:
            return None
        return shares / total

    # --------------------------------------------------------
    def caltable_status(self, caltable):
        """Validate global and per-group self-calibration solution coverage."""

        table_path = Path(caltable).expanduser().resolve()
        qa = {
            "schema_version": 1,
            "table": str(table_path),
            "status": "failed",
            "reason": "",
        }

        def json_value(value):
            item = getattr(value, "item", None)
            if callable(item):
                try:
                    value = item()
                except (TypeError, ValueError):
                    pass
            if value is None or isinstance(value, (bool, int, float, str)):
                return value
            return str(value)

        def finish(status, reason=""):
            qa["status"] = status
            qa["reason"] = reason
            if not hasattr(self, "caltable_qa"):
                self.caltable_qa = {}
            self.caltable_qa[str(table_path)] = dict(qa)

            qa_dir = getattr(self, "caltabledir", table_path.parent)
            sidecar = Path(qa_dir) / f"{table_path.name}.qa.json"
            sidecar.parent.mkdir(parents=True, exist_ok=True)
            temporary = sidecar.with_name(
                f".{sidecar.name}.tmp.{os.getpid()}"
            )
            try:
                with open(temporary, "w") as handle:
                    json.dump(self.to_builtin(qa), handle, indent=2)
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(temporary, sidecar)
            finally:
                if temporary.exists():
                    temporary.unlink()

            return status == "passed", reason

        try:
            max_unusable = float(
                self.task_config.get(
                    "max_caltable_unusable_fraction",
                    0.50,
                )
            )
            max_group_unusable = float(
                self.task_config.get(
                    "max_caltable_group_unusable_fraction",
                    0.80,
                )
            )
            # A self-calibration solve is not a calibrator solve: an antenna
            # or an interval that cannot be solved keeps the calibration it
            # already has (applymode calonly), so a bounded number of them
            # is a property of the data, not a reason to discard every other
            # antenna's solution.
            tolerated_group_failures = {
                "ANTENNA1": float(
                    self.task_config.get(
                        "max_caltable_failed_antenna_fraction",
                        0.25,
                    )
                ),
                "TIME": float(
                    self.task_config.get(
                        "max_caltable_failed_time_fraction",
                        0.25,
                    )
                ),
            }
        except (TypeError, ValueError):
            return finish(
                "failed",
                "self-calibration solution-coverage limits must be numeric",
            )

        qa["max_unusable_solution_fraction"] = max_unusable
        qa["max_group_unusable_solution_fraction"] = max_group_unusable
        qa["max_failed_group_fraction"] = dict(tolerated_group_failures)
        if (
            not math.isfinite(max_unusable) or
            not 0.0 <= max_unusable <= 1.0 or
            not math.isfinite(max_group_unusable) or
            not 0.0 <= max_group_unusable <= 1.0 or
            any(
                not math.isfinite(value) or not 0.0 <= value <= 1.0
                for value in tolerated_group_failures.values()
            )
        ):
            return finish(
                "failed",
                "self-calibration solution-coverage limits must be between "
                "zero and one",
            )

        if not table_path.is_dir():
            return finish("failed", "the table directory was not created")

        if not (table_path / "table.dat").is_file():
            return finish(
                "failed",
                "the table directory has no table.dat file",
            )

        opened = False

        try:
            tb.open(str(table_path))
            opened = True
            rows = int(tb.nrows())
            columns = set(tb.colnames())
            qa["rows"] = rows

            if rows < 1:
                return finish("failed", "the table contains no solution rows")

            value_column = None
            for candidate in ("CPARAM", "FPARAM"):
                if candidate in columns:
                    value_column = candidate
                    break

            if value_column is None:
                return finish(
                    "failed",
                    "the table has no CPARAM or FPARAM column",
                )

            if "FLAG" not in columns:
                return finish("failed", "the table has no FLAG column")

            required_groups = ["ANTENNA1", "SPECTRAL_WINDOW_ID"]
            missing_groups = [
                name for name in required_groups if name not in columns
            ]
            if missing_groups:
                return finish(
                    "failed",
                    "the table lacks required solution-group column(s): "
                    + ", ".join(missing_groups),
                )

            values = np.asarray(tb.getcol(value_column))
            flags = np.asarray(tb.getcol("FLAG"), dtype=bool)
            qa["parameter_column"] = value_column
            qa["solution_shape"] = list(values.shape)

            if values.size < 1:
                return finish("failed", "the table contains no solutions")
            if values.shape != flags.shape:
                return finish(
                    "failed",
                    f"{value_column} and FLAG shapes differ: "
                    f"{values.shape} != {flags.shape}",
                )
            if values.ndim < 1 or values.shape[-1] != rows:
                return finish(
                    "failed",
                    "solution-array rows do not align with the calibration "
                    "table row count",
                )

            try:
                finite = np.isfinite(values)
                nonzero = np.abs(values) > 0.0
            except (TypeError, ValueError) as exc:
                return finish(
                    "failed",
                    f"the solution column is not numeric ({exc})",
                )

            # Antennas with no data to solve from (dead, quarantined, or
            # flagged away upstream) always return flagged solutions. They
            # are a fact about the Measurement Set, measured before the
            # solve, so their rows are left out of every coverage figure
            # instead of letting them condemn the rest of the array.
            excluded_antennas = sorted(
                int(value) for value in (
                    getattr(self, "unsolvable_antenna_ids", None) or ()
                )
            )
            keep_rows = np.ones(rows, dtype=bool)
            if excluded_antennas:
                antenna_rows = np.asarray(
                    tb.getcol("ANTENNA1")
                ).reshape(-1)
                if len(antenna_rows) == rows:
                    keep_rows = ~np.isin(antenna_rows, excluded_antennas)
            qa["excluded_unsolvable_antenna_ids"] = excluded_antennas
            qa["excluded_solution_rows"] = int(
                rows - np.count_nonzero(keep_rows)
            )
            if not np.any(keep_rows):
                return finish(
                    "failed",
                    "every solution row belongs to an antenna that holds "
                    "no usable data",
                )

            usable = ((~flags) & finite & nonzero)[..., keep_rows]
            unflagged_nonfinite = ((~flags) & (~finite))[..., keep_rows]
            unflagged_zero = (
                (~flags) & finite & (~nonzero)
            )[..., keep_rows]
            flags = flags[..., keep_rows]
            total = int(usable.size)
            usable_count = int(np.count_nonzero(usable))
            unusable_count = total - usable_count
            unusable_fraction = unusable_count / total
            qa.update({
                "total_solutions": total,
                "usable_solutions": usable_count,
                "unusable_solutions": unusable_count,
                "unusable_fraction": unusable_fraction,
                "flagged_solutions": int(np.count_nonzero(flags)),
                "unflagged_nonfinite_solutions": int(
                    np.count_nonzero(unflagged_nonfinite)
                ),
                "unflagged_zero_solutions": int(
                    np.count_nonzero(unflagged_zero)
                ),
            })

            group_metrics = {}
            for group_column in [
                "ANTENNA1",
                "SPECTRAL_WINDOW_ID",
                "FIELD_ID",
                "OBSERVATION_ID",
                "TIME",
            ]:
                if group_column not in columns:
                    continue
                group_values = np.asarray(
                    tb.getcol(group_column)
                ).reshape(-1)
                if len(group_values) != rows:
                    return finish(
                        "failed",
                        f"{group_column} has {len(group_values)} values for "
                        f"{rows} calibration-table rows",
                    )
                group_values = group_values[keep_rows]

                # Group with np.unique rather than a nested scan.  TIME on a
                # solint='int' table has as many distinct values as solution
                # intervals, and the previous pure-Python comparison was
                # quadratic in the row count: a routine self-cal table could
                # stall QA for hours.
                groups = []
                unique_values, inverse = np.unique(
                    group_values,
                    return_inverse=True,
                )
                inverse = np.asarray(inverse).reshape(-1)
                for index, raw_value in enumerate(unique_values):
                    row_mask = inverse == index
                    group_usable = usable[..., row_mask]
                    group_total = int(group_usable.size)
                    group_usable_count = int(
                        np.count_nonzero(group_usable)
                    )
                    group_unusable = group_total - group_usable_count
                    fraction = (
                        group_unusable / group_total
                        if group_total else 1.0
                    )
                    status = (
                        "passed"
                        if group_usable_count > 0 and
                        fraction <= max_group_unusable
                        else "failed"
                    )
                    groups.append({
                        "value": json_value(raw_value),
                        "rows": int(np.count_nonzero(row_mask)),
                        "total_solutions": group_total,
                        "usable_solutions": group_usable_count,
                        "unusable_solutions": group_unusable,
                        "unusable_fraction": fraction,
                        "status": status,
                    })
                group_metrics[group_column] = groups
            qa["group_metrics"] = group_metrics

        except Exception as exc:
            return finish(
                "failed",
                f"CASA could not inspect the table ({exc})",
            )

        finally:
            if opened:
                tb.close()

        if usable_count < 1:
            return finish(
                "failed",
                "all calibration solutions are flagged or invalid",
            )
        if qa["unflagged_nonfinite_solutions"]:
            return finish(
                "failed",
                "the table contains unflagged non-finite solutions",
            )
        if qa["unflagged_zero_solutions"]:
            return finish(
                "failed",
                "the table contains unflagged zero-valued solutions",
            )
        if unusable_fraction > max_unusable:
            return finish(
                "failed",
                "unusable calibration-solution fraction "
                f"{unusable_fraction:.6f} exceeds the configured maximum "
                f"{max_unusable:.6f}",
            )
        # What a failed interval costs is the data it covers, not its place
        # in a list: one solution per scan over scans of 9.5 min, 22 min and
        # 86 s (uGMRT Band 4, G0212-0335) made the last scan a third of the
        # intervals and 4% of the observation.
        time_groups = group_metrics.get("TIME") or []
        time_shares = self.solution_time_shares(
            [group["value"] for group in time_groups]
        )
        if time_shares is not None:
            for group, share in zip(time_groups, time_shares):
                group["observing_time_share"] = float(share)

        tolerated = {}
        for column_name, groups in group_metrics.items():
            failed = [
                group for group in groups if group["status"] == "failed"
            ]
            if not failed:
                continue
            allowed = tolerated_group_failures.get(column_name)
            failed_fraction = len(failed) / float(len(groups))
            extent = (
                f"{len(failed)} of {len(groups)} {column_name} groups failed"
            )
            if column_name == "TIME" and time_shares is not None:
                failed_fraction = float(sum(
                    group["observing_time_share"] for group in failed
                ))
                qa["failed_observing_time_share"] = failed_fraction
                extent += (
                    f", covering {failed_fraction:.1%} of the observing time"
                )
            if allowed is not None and failed_fraction <= allowed:
                tolerated[column_name] = [
                    group["value"] for group in failed
                ]
                continue
            group = failed[0]
            return finish(
                "failed",
                f"{column_name}={group['value']!r} has unusable "
                f"solution fraction {group['unusable_fraction']:.6f}, "
                "or no usable solution"
                + (
                    f" ({extent}; at most {allowed:.0%} are tolerated)"
                    if allowed is not None else ""
                ),
            )

        qa["tolerated_failed_groups"] = tolerated
        if tolerated:
            print(
                "Self-calibration table "
                f"{table_path.name}: solutions are missing for "
                + "; ".join(
                    f"{len(values)} {name} group(s) {values}"
                    for name, values in tolerated.items()
                )
                + ". Those data keep the calibration they already have."
            )

        return finish("passed")

    # --------------------------------------------------------

    def archive_incomplete_caltable(self, caltable):
        """Preserve an incomplete table so a fresh gaincal solve can run."""

        table_path = Path(caltable)

        if not table_path.exists():

            return None

        index = 1

        while True:

            archived = table_path.with_name(
                f"{table_path.name}.incomplete.{index}"
            )

            if not archived.exists():

                break

            index += 1

        table_path.rename(archived)

        print(
            "Preserved incomplete calibration table as "
            f"{archived}"
        )

        return archived

    # --------------------------------------------------------

    def prepare_caltable(self, caltable):
        """Make a calibration-table path safe for a fresh gaincal solve.

        A table without a matching solve checkpoint may belong to a previous
        run with different model data, flags, or solve settings. Preserve it
        for diagnosis, but never apply or silently reuse it.
        """

        valid, reason = self.caltable_status(caltable)

        if Path(caltable).exists():

            if valid:

                print(
                    "Existing calibration table has no matching solve "
                    "checkpoint and will not be reused."
                )

            else:

                print(
                    "Existing calibration table is unusable: "
                    f"{reason}."
                )

            self.archive_incomplete_caltable(caltable)

        return False

    # --------------------------------------------------------

    def archive_incomplete_product(self, path, label):
        """Move an uncheckpointed product aside, keeping it for inspection."""

        path = Path(path)
        if not (path.exists() or path.is_symlink()):
            return None

        index = 1
        while True:
            archived = path.with_name(f"{path.name}.incomplete.{index}")
            if not (archived.exists() or archived.is_symlink()):
                break
            index += 1

        path.rename(archived)
        if not hasattr(self, "archived_incomplete_products"):
            self.archived_incomplete_products = []
        self.archived_incomplete_products.append({
            "product": label,
            "original": str(path),
            "archived": str(archived),
        })
        print(
            f"Preserved incomplete {label} from an interrupted run as "
            f"{archived}"
        )
        return archived

    def archive_incomplete_images(self, image_base, label):
        """Move aside every product of one uncheckpointed tclean run.

        tclean reuses existing images of the same name by default, so a
        partial set from a killed run would otherwise seed the new one.
        """

        base = Path(image_base)
        if not base.parent.is_dir():
            return []

        archived = []
        for product in sorted(base.parent.glob(base.name + ".*")):
            if ".incomplete." in product.name:
                continue
            result = self.archive_incomplete_product(product, label)
            if result is not None:
                archived.append(result)
        return archived

    def recover_interrupted_phase_loop(self, loop):
        """Clear what a killed, uncheckpointed phase loop left behind.

        The loop works on a candidate copy of the last accepted MS and never
        modifies that MS, so redoing it from the accepted MS cannot apply a
        solution twice. Its scratch products are preserved, not deleted.
        """

        for product in (
            self.msdir / f"phase{loop}.candidate.ms",
            self.msdir / f"phase{loop}.ms",
        ):
            self.archive_incomplete_product(
                product,
                f"phase-selfcal loop {loop} Measurement Set",
            )

        self.archive_incomplete_images(
            self.imagedir / f"phase_{loop}",
            f"phase-selfcal loop {loop} image",
        )

    def discard_phase_candidate(self, loop):
        """Remove only the temporary MS used to evaluate one phase solve."""

        candidate_ms = self.msdir / f"phase{loop}.candidate.ms"

        if candidate_ms.exists():

            shutil.rmtree(candidate_ms)

    # --------------------------------------------------------

    def commit_phase_candidate(self, loop):
        """Commit an accepted temporary phase-selfcal MS.

        The temporary MS contains the newly corrected data. Splitting its
        CORRECTED_DATA into the committed MS prevents a rejected candidate
        from modifying the previously accepted MS.
        """

        candidate_ms = self.msdir / f"phase{loop}.candidate.ms"
        split_ms = self.msdir / f"phase{loop}.ms"

        if not candidate_ms.exists():

            raise RuntimeError(
                "The temporary phase-selfcal candidate is missing: "
                f"{candidate_ms}"
            )

        if split_ms.exists():

            raise RuntimeError(
                "An unaccepted phase-selfcal candidate already exists:\n"
                f"{split_ms}\n"
                "Review or remove the candidate before rerunning this "
                "self-calibration loop."
            )

        tb.open(str(candidate_ms))
        cols = tb.colnames()
        tb.close()

        if "CORRECTED_DATA" not in cols:

            raise RuntimeError(
                "Temporary phase-selfcal candidate has no "
                "CORRECTED_DATA column."
            )

        print(
            f"Committing {split_ms.name}"
        )

        split(
            vis=str(candidate_ms),
            outputvis=str(split_ms),
            datacolumn="corrected",
            keepflags=True,
        )

        self.discard_phase_candidate(loop)
        self.current_ms = str(split_ms)

        return self.current_ms

    # --------------------------------------------------------

    # Save JSON helper
    # --------------------------------------------------------

    def save_json(self,
                  data,
                  outfile):

        with open(outfile, "w") as f:

            json.dump(
                self.to_builtin(data),
                f,
                indent=4
            )

    # --------------------------------------------------------

    @classmethod
    def to_builtin(cls, value):
        """Convert CASA/NumPy containers and scalars to JSON-safe values."""

        if hasattr(value, "tolist"):

            value = value.tolist()

        if isinstance(value, dict):

            return {
                str(key): cls.to_builtin(item)
                for key, item in value.items()
            }

        if isinstance(value, (list, tuple)):

            return [
                cls.to_builtin(item)
                for item in value
            ]

        try:

            json.dumps(value)
            return value

        except TypeError:

            return str(value)

    # --------------------------------------------------------
    # Print summary
    # --------------------------------------------------------

    def print_summary(self):

        print("\n")
        print("=" * 70)
        print(
            f"{TELESCOPE_PROFILE_MODULE.get_telescope_profile(self.telescope)['label']} "
            "Self-Calibration Summary"
        )
        print("=" * 70)

        print(
            f"Measurement Set : {self.ms}"
        )

        print(
            f"Band            : {self.band}"
        )

        print(
            f"Profile         : {self.profile}"
        )

        print(
            f"Phase solints   : "
            f"{self.band_config['phase_solints']}"
        )

        print(
            f"Amp solint      : "
            f"{self.band_config['amp_solint']}"
        )

        print(
            f"Min SNR phase   : "
            f"{self.band_config['min_snr_phase']}"
        )

        print(
            f"Min SNR amp     : "
            f"{self.band_config['min_snr_amp']}"
        )

        print(
            f"Output          : "
            f"{self.selfcaldir}"
        )

        if hasattr(self, "final_metrics"):
            print("\n")

            print(
                f"Final RMS : "
                f"{self.final_metrics['rms_Jy_beam']:.3e}"
            )

            print(
                f"Peak Flux : "
                f"{self.final_metrics['peak_Jy_beam']:.3e}"
            )

            print(
                f"Dynamic Range : "
                f"{self.final_metrics['dynamic_range']:.0f}"
            )

        print("=" * 70)

    # --------------------------------------------------------
    # Adaptive automasking
    # --------------------------------------------------------

    def automask_kwargs(self, previous_image=None, relaxed=False):
        """Return conservative masking or one bounded empty-model recovery."""

        if relaxed:
            sidelobe = float(
                self.task_config.get(
                    "automask_retry_sidelobethreshold",
                    1.5,
                )
            )
            noise = float(
                self.task_config.get(
                    "automask_retry_noisethreshold",
                    3.5,
                )
            )
            minbeam = float(
                self.task_config.get(
                    "automask_retry_minbeamfrac",
                    0.05,
                )
            )
            growiterations = int(
                self.task_config.get(
                    "automask_retry_growiterations",
                    150,
                )
            )
        else:
            sidelobe = float(
                self.task_config.get(
                    "automask_sidelobethreshold",
                    2.5,
                )
            )
            noise = float(
                self.task_config.get(
                    "automask_noisethreshold",
                    4.25,
                )
            )
            minbeam = float(
                self.task_config.get(
                    "automask_minbeamfrac",
                    0.20,
                )
            )
            growiterations = int(
                self.task_config.get(
                    "automask_growiterations",
                    75,
                )
            )

        values = {
            "usemask": "auto-multithresh",
            "sidelobethreshold": max(1.0, sidelobe),
            "noisethreshold": max(3.0, noise),
            "lownoisethreshold": float(
                self.task_config.get(
                    "automask_lownoisethreshold",
                    1.5,
                )
            ),
            "minbeamfrac": max(0.01, minbeam),
            "growiterations": max(1, growiterations),
            "negativethreshold": float(
                self.task_config.get(
                    "automask_negativethreshold",
                    0.0,
                )
            ),
            "smoothfactor": float(
                self.task_config.get(
                    "automask_smoothfactor",
                    1.0,
                )
            ),
            "cutthreshold": float(
                self.task_config.get(
                    "automask_cutthreshold",
                    0.01,
                )
            ),
            "dogrowprune": bool(
                self.task_config.get(
                    "automask_dogrowprune",
                    True,
                )
            ),
            "minpercentchange": float(
                self.task_config.get(
                    "automask_minpercentchange",
                    -1.0,
                )
            ),
            "fastnoise": bool(
                self.task_config.get(
                    "automask_fastnoise",
                    False,
                )
            ),
        }
        print(
            "Auto-multithresh: "
            f"sidelobethreshold={values['sidelobethreshold']}, "
            f"noisethreshold={values['noisethreshold']}, "
            f"minbeamfrac={values['minbeamfrac']}, "
            f"fastnoise={values['fastnoise']}"
        )
        return values

    # --------------------------------------------------------

    def gridder_kwargs(self):
        """Return the wide-field gridder shared by every self-cal image.

        CASA's tclean default is ``wprojplanes=1``, which disables
        w-projection even when ``gridder='wproject'`` is selected.  Self-cal
        images must use the same w-plane depth as final imaging, otherwise
        the sky model that the gains are solved against, and the
        dynamic-range test that accepts or rejects each loop, are both
        measured on a non-coplanar-distorted image.
        """

        return {
            "gridder": "wproject",
            "wprojplanes": int(
                self.image_config.get("wprojplanes", 128)
            ),
            # Keep the outer field in the model rather than letting CASA's
            # internal beam model blank it; those sources still contribute
            # visibility power that the gain solve would otherwise absorb.
            "pblimit": float(self.task_config.get("pblimit", 0.05)),
        }

    # --------------------------------------------------------

    def clean_control_kwargs(self):
        """Stable minor-cycle controls shared by all self-cal images."""

        controls = self.minor_cycle_kwargs()
        # Multi-scale sizes, when the band configures them. Without, the
        # model is point components, which cannot hold extended emission:
        # on the uGMRT field that prompted the option it held 0.57 Jy of the
        # 0.78 Jy a multi-scale clean of the same data found.
        scales = self.image_config.get("scales")
        if scales:
            controls["scales"] = [int(item) for item in scales]
        return controls

    def minor_cycle_kwargs(self):

        return {
            "gain": float(
                self.task_config.get(
                    "clean_gain",
                    0.1,
                )
            ),
            "cyclefactor": float(
                self.task_config.get(
                    "clean_cyclefactor",
                    1.5,
                )
            ),
            "minpsffraction": float(
                self.task_config.get(
                    "clean_minpsffraction",
                    0.05,
                )
            ),
            "maxpsffraction": float(
                self.task_config.get(
                    "clean_maxpsffraction",
                    0.8,
                )
            ),
            "smallscalebias": float(
                self.task_config.get(
                    "clean_smallscalebias",
                    0.0,
                )
            ),
        }

    # --------------------------------------------------------

    def deconvolution_diagnostics(self, imagename):
        """Measure CLEAN model occupancy and robust residual S/N."""

        def first_image(suffixes):
            for suffix in suffixes:
                candidate = Path(str(imagename) + suffix)
                if candidate.is_dir() and (candidate / "table.dat").is_file():
                    return str(candidate)
            return None

        model_image = first_image([
            ".model.tt0",
            ".model",
        ])
        residual_image = first_image([
            ".residual.tt0",
            ".residual",
        ])
        diagnostics = {
            "model_image": model_image,
            "residual_image": residual_image,
            "model_nonzero": False,
            "model_peak_abs_Jy": None,
            "residual_rms_Jy_beam": None,
            "residual_peak_abs_Jy_beam": None,
            "residual_peak_snr": None,
        }

        if model_image:
            try:
                stats = imstat(imagename=model_image)
                model_peak = max(
                    abs(float(stats["max"][0])),
                    abs(float(stats["min"][0])),
                )
                diagnostics["model_peak_abs_Jy"] = model_peak
                diagnostics["model_nonzero"] = model_peak > 1.0e-12
                try:
                    diagnostics["model_flux_Jy"] = float(stats["sum"][0])
                except (KeyError, IndexError, TypeError, ValueError):
                    diagnostics["model_flux_Jy"] = None
                try:
                    # By Parseval, the mean square of the model's visibilities
                    # over the uv plane is the sum of its squared pixels, so
                    # this is the amplitude a gain solve has to work with.
                    # Components spread over several pixels make it an
                    # underestimate on short baselines, so it errs low.
                    diagnostics["model_rms_visibility_Jy"] = math.sqrt(
                        max(0.0, float(stats["sumsq"][0]))
                    )
                except (KeyError, IndexError, TypeError, ValueError):
                    diagnostics["model_rms_visibility_Jy"] = None
            except Exception as exc:
                diagnostics["model_error"] = str(exc)

        if residual_image:
            try:
                stats = CASA_COMPAT.peak_and_noise(imstat, residual_image)
                rms = abs(stats["rms"]) if stats["rms"] else None
                diagnostics["residual_rms_Jy_beam"] = rms
                diagnostics["residual_peak_abs_Jy_beam"] = stats["peak_abs"]
                diagnostics["residual_peak_position"] = (
                    stats["peak_abs_position"]
                )
                diagnostics["residual_peak_snr"] = stats["peak_abs_snr"]
                diagnostics["noise_only_peak_snr"] = (
                    stats.get("noise_only_peak_snr")
                )
            except Exception as exc:
                diagnostics["residual_error"] = str(exc)

        return diagnostics

    # --------------------------------------------------------


    def clear_stale_clean_mask(self, imagename):
        """Remove a mask left by an earlier pass before re-cleaning.

        A retry changes how the mask is made: from auto-multithresh to
        different thresholds, or to the primary beam. tclean restarted over an
        existing ``.mask`` from the previous pass mixes the two, which is a
        fragile state in CASA. Removing it makes each pass build its own.
        """

        mask = Path(str(imagename) + ".mask")
        if not (mask.is_dir() and (mask / "table.dat").is_file()):
            return False
        try:
            rmtables(tablenames=str(mask))
        except Exception as exc:
            print(f"WARNING: could not remove {mask.name}: {exc}")
            return False
        if mask.exists():
            try:
                shutil.rmtree(mask)
            except OSError:
                pass
        return True

    def primary_beam_image_ready(self, imagename):
        """Whether a primary-beam image exists to build a pb mask from.

        ``usemask='pb'`` reads it, and with ``calcpsf=False`` tclean does not
        make one. Asking for a mask that cannot be built is not worth the
        risk of finding out how CASA fails.
        """

        for suffix in (".pb.tt0", ".pb"):
            candidate = Path(str(imagename) + suffix)
            if candidate.is_dir() and (candidate / "table.dat").is_file():
                return True
        return False

    def saved_record(self, name):
        """A record an earlier run of this stage wrote to the report folder."""

        reportdir = getattr(self, "reportdir", None)
        if reportdir is None:
            return None
        return self.read_json(Path(reportdir) / name)

    def keep_record(self, record, name):
        """Write a record where a resumed stage will find it.

        A resumed stage restores its loops from checkpoints and repeats
        neither their images nor their flagging, so what those steps
        reported exists only in what the first run wrote down. The report
        must not be what stops a stage: a record that cannot be written is
        only missing from a later report.
        """

        reportdir = getattr(self, "reportdir", None)
        if reportdir is None:
            return
        try:
            self.save_json(record, Path(reportdir) / name)
        except (OSError, TypeError, ValueError):
            pass

    def recorded_clean_passes(self):
        """CASA's account of every clean pass, including earlier runs'."""

        passes = getattr(self, "clean_passes", None)
        if passes is None:
            saved = self.saved_record("clean_passes.json")
            passes = self.clean_passes = (
                saved if isinstance(saved, dict) else {}
            )
        return passes

    def tclean_pass(self, kwargs, label, name):
        """Run one tclean call and keep CASA's account of why it stopped."""

        summary = CASA_COMPAT.summarize_tclean(tclean(**kwargs))
        passes = self.recorded_clean_passes()
        passes.setdefault(label, []).append(dict(summary or {}, name=name))
        self.keep_record(passes, "clean_passes.json")
        if summary and summary.get("stop_reason"):
            print(
                f"CLEAN {name} of {label} {summary['stop_reason']}"
                + (
                    f" after {summary['iterdone']} of {summary['niter']} "
                    "iterations"
                    if summary.get("iterdone") is not None and
                    summary.get("niter") is not None else ""
                )
                + (
                    f" and {summary['nmajordone']} major cycles"
                    if summary.get("nmajordone") is not None else ""
                )
                + (
                    f"; PSF largest sidelobe {summary['max_psf_sidelobe']:.2f}"
                    if summary.get("max_psf_sidelobe") else ""
                )
                + "."
            )
        return summary

    def continue_to_convergence(self, kwargs, label, summary, diagnostics):
        """Carry on cleaning while CLEAN stopped only for want of iterations.

        An unfinished model is a worse model to self-calibrate against: the
        flux it lacks is absorbed into the gains. CLEAN that reaches its
        iteration limit with significant emission left continues from where
        it stopped, with its mask and model, for at most
        ``clean_max_continuations`` further runs of ``niter``.

        A pass stopped by ``clean_max_major_cycles`` is not continued: that
        limit exists for a masked pass that has stopped converging, and
        continuing it would undo the limit. The significance tested here is
        the whole image's, which in a field limited by calibration errors is
        an artefact beside the brightest source that no cleaning removes.
        """

        limit = int(self.task_config.get("clean_max_continuations", 3))
        configured_floor = self.task_config.get("clean_continue_min_peak_snr", 5.0)
        floor = significance_floor(configured_floor, diagnostics)
        rounds = 0
        carried = {
            key: diagnostics[key] for key in ("bootstrap",)
            if key in diagnostics
        }
        while (
            summary and summary.get("stopcode") == 1 and
            rounds < limit and
            (diagnostics.get("residual_peak_snr") or 0.0) >= floor
        ):
            rounds += 1
            print(
                f"CLEAN of {label} stopped at its iteration limit with a "
                f"{diagnostics['residual_peak_snr']:.1f}-sigma residual peak. "
                f"Continuing ({rounds} of at most {limit})."
            )
            more = dict(kwargs)
            more.update({"restart": True, "calcpsf": False, "calcres": False})
            summary = self.tclean_pass(more, label, f"continuation {rounds}")
            diagnostics = self.deconvolution_diagnostics(kwargs["imagename"])
            diagnostics.update(carried)

        diagnostics["continuations"] = rounds
        diagnostics["final_stop"] = summary
        return diagnostics

    def unmasked_clean_kwargs(self, kwargs, nsigma, niter):
        """tclean arguments for a clean bounded by nsigma instead of a mask.

        CASA blanks everything below ``pblimit`` of its own beam model, so
        with a positive pblimit "unmasked" still means "inside that beam",
        and the mask is the beam at that level. A negative pblimit switches
        CASA's blanking off; the clean then reaches the whole image, so a
        bright source in the outer primary beam or its first sidelobe is
        deconvolved with the rest instead of leaving its sidelobes across
        the field.
        """

        unmasked = dict(kwargs)
        for key in ("sidelobethreshold", "noisethreshold", "lownoisethreshold",
                    "negativethreshold", "minbeamfrac", "growiterations",
                    "dogrowprune", "minpercentchange", "fastnoise",
                    "smoothfactor", "cutthreshold", "verbose", "pbmask",
                    "mask"):
            unmasked.pop(key, None)

        configured = self.task_config.get("clean_bootstrap_pbmask")
        try:
            pblimit = float(kwargs.get("pblimit"))
        except (TypeError, ValueError):
            pblimit = None

        if configured is None and pblimit is not None and pblimit < 0.0:
            unmasked.update({"usemask": "user", "mask": ""})
        else:
            unmasked.update({
                "usemask": "pb",
                "pbmask": CASA_COMPAT.bootstrap_pbmask(
                    configured,
                    kwargs.get("pblimit"),
                ),
            })
        unmasked.update({"nsigma": float(nsigma), "niter": int(niter)})
        return unmasked

    def adapt_automask_to_psf(self, kwargs, label):
        """Lower the auto-multithresh sidelobe threshold to suit this PSF.

        The first pass on this image measured the PSF's largest sidelobe.
        Where it is high, the configured threshold asks for emission above
        the residual peak itself and no mask can open; it is capped so that
        the criterion sits at a fixed fraction of the peak instead.
        """

        if (
            kwargs.get("usemask") != "auto-multithresh" or
            "sidelobethreshold" not in kwargs
        ):
            return kwargs

        sidelobe = None
        for summary in reversed(
            (getattr(self, "clean_passes", None) or {}).get(label, [])
        ):
            if summary.get("max_psf_sidelobe"):
                sidelobe = summary["max_psf_sidelobe"]
                break

        fraction = self.task_config.get("automask_max_sidelobe_fraction", 0.6)
        value, lowered = CASA_COMPAT.adaptive_sidelobethreshold(
            kwargs["sidelobethreshold"],
            sidelobe,
            fraction,
        )
        if not lowered:
            return kwargs

        print(
            f"The PSF of {label} has a largest sidelobe of {sidelobe:.2f}. "
            f"sidelobethreshold={kwargs['sidelobethreshold']} would require "
            "emission above "
            f"{float(kwargs['sidelobethreshold']) * sidelobe:.2f} of the "
            f"residual peak before any mask opens; using {value:.2f} "
            f"({value * sidelobe:.2f} of the peak) for this image."
        )
        adapted = dict(kwargs)
        adapted["sidelobethreshold"] = value
        if not hasattr(self, "automask_adaptations"):
            self.automask_adaptations = {}
        self.automask_adaptations[label] = {
            "max_psf_sidelobe": sidelobe,
            "configured_sidelobethreshold": kwargs["sidelobethreshold"],
            "used_sidelobethreshold": value,
        }
        return adapted

    def cap_masked_major_cycles(self, kwargs):
        """Bound the masked pass by major cycles as well as by depth.

        On a short observation's point-spread function the masked pass does
        not settle onto its nsigma threshold. In the uGMRT Band 4 run that
        prompted this, the peak residual inside the mask at the start of
        each major cycle went 831, 592, 485, 367, 639, 383, 353, 597 ...
        uJy/beam: a deep minor cycle ends near the threshold, and the next
        major cycle finds a peak twice as high, which a few iterations
        remove. Those passes took 13 to 21 major cycles of 6 to 16 minutes,
        and CASA reported "Possible divergence" in the minor cycle of three
        of them. Nothing but the nsigma test and niter ends such a pass, and
        a pass that exhausts niter is then continued. ``nmajor`` makes the
        number of major cycles a limit of its own:
        ``clean_max_major_cycles``, ``DEFAULT_CLEAN_MAX_MAJOR_CYCLES`` unless
        configured, 0 for none.
        """

        value = self.task_config.get(
            "clean_max_major_cycles", DEFAULT_CLEAN_MAX_MAJOR_CYCLES
        )
        try:
            cap = int(value) if value is not None else 0
        except (TypeError, ValueError):
            cap = 0
        if cap <= 0 or kwargs.get("usemask") != "auto-multithresh":
            return kwargs
        capped = dict(kwargs)
        capped["nmajor"] = cap
        return capped

    def floor_clean(self, kwargs, label, diagnostics):
        """Finish with an nsigma-bounded clean when no mask could be opened.

        auto-multithresh does the deep clean inside the masks it draws. When
        it draws none, what is left is the residual of the shallow bootstrap,
        with every source still carrying the part of its point-spread
        function below that depth. Cleaning what remains above
        ``clean_floor_nsigma`` guarantees the image is never left shallower
        than a plain threshold clean would have made it.
        """

        floor = float(self.task_config.get("clean_floor_nsigma", 5.0))
        niter = int(self.task_config.get("clean_bootstrap_niter", 20000))
        if (
            floor <= 0 or niter <= 0 or
            not bool(self.task_config.get("clean_bootstrap_enabled", True)) or
            not diagnostics.get("model_nonzero") or
            not CASA_COMPAT.is_zero_mask_stop(diagnostics.get("final_stop"))
        ):
            return diagnostics

        peak_snr = diagnostics.get("residual_peak_snr")
        if peak_snr is None or peak_snr < floor:
            return diagnostics

        print(
            f"auto-multithresh opened no mask for {label}, and the residual "
            f"still holds a {peak_snr:.1f}-sigma peak. Cleaning what remains "
            f"above {floor:.1f} sigma, at most {niter} iterations."
        )
        self.clear_stale_clean_mask(kwargs["imagename"])
        floor_kwargs = self.unmasked_clean_kwargs(kwargs, floor, niter)
        floor_kwargs.update({
            "restart": True,
            "calcpsf": False,
            "calcres": False,
        })
        summary = self.tclean_pass(floor_kwargs, label, "unmasked floor")
        after = self.deconvolution_diagnostics(kwargs["imagename"])
        for key in ("bootstrap", "continuations"):
            if key in diagnostics:
                after[key] = diagnostics[key]
        after["final_stop"] = summary
        after["floor_clean"] = {
            "nsigma": floor,
            "residual_peak_snr_before": peak_snr,
            "residual_rms_before_Jy_beam": diagnostics.get(
                "residual_rms_Jy_beam"
            ),
        }
        return after

    def bootstrap_clean(self, kwargs, label):
        """Take the bright sources out before auto-multithresh measures noise.

        auto-multithresh opens a mask at ``max(sidelobethreshold *
        max_sidelobe * peak, noisethreshold * rms)``, both measured on the
        residual, which on the first pass is the dirty map. Where a field
        holds several comparable sources that map's rms is sidelobe
        confusion rather than noise, so the threshold sits at or above the
        peak and no mask opens. For self-calibration that is fatal twice
        over: there is no model to solve against, so the loop is skipped and
        the run continues on uncalibrated gains.

        Cleaning the plainly real emission first breaks it. ``nsigma`` is
        recomputed from the residual each major cycle, so the first cycle
        reaches only the brightest peaks and the threshold follows the noise
        down as their sidelobes come out. The pass is bounded by nsigma and a
        capped iteration count and writes no model column; auto-multithresh
        still does the deep, structure-following clean afterwards.
        """

        if not bool(self.task_config.get("clean_bootstrap_enabled", True)):
            return None

        nsigma = float(self.task_config.get("clean_bootstrap_nsigma", 8.0))
        niter = int(self.task_config.get("clean_bootstrap_niter", 20000))
        if niter <= 0 or nsigma <= 0:
            return None

        bootstrap_kwargs = self.unmasked_clean_kwargs(kwargs, nsigma, niter)
        bootstrap_kwargs["savemodel"] = "none"
        pbmask = bootstrap_kwargs.get("pbmask")

        print(
            f"Bootstrap clean of {label}: taking emission above "
            f"{nsigma:.1f} sigma "
            + (
                f"inside pb>{pbmask:.2f} " if pbmask is not None
                else "anywhere in the image "
            )
            + f"out of the dirty map first, at most {niter} iterations, so "
            "that auto-multithresh measures noise rather than sidelobes."
        )

        self.clear_stale_clean_mask(kwargs["imagename"])
        self.tclean_pass(bootstrap_kwargs, label, "bootstrap")
        diagnostics = self.deconvolution_diagnostics(kwargs["imagename"])
        rms = diagnostics.get("residual_rms_Jy_beam")
        print(
            f"Bootstrap clean of {label}: "
            + ("model present" if diagnostics.get("model_nonzero")
               else "nothing above the threshold")
            + (f", residual RMS {rms:.3e} Jy/beam" if rms else "")
        )
        return diagnostics

    def run_tclean_with_recovery(self, kwargs, label):
        """Run self-cal imaging and retry only an empty high-S/N model."""

        # An image made again replaces the account of the one it replaces.
        self.recorded_clean_passes().pop(label, None)

        bootstrap = self.bootstrap_clean(kwargs, label)
        if bootstrap is not None:
            # The bootstrap made the PSF and the first residual; the
            # auto-multithresh pass continues from them, adding its
            # components to the ones already found, and writes the model
            # column the gain solve reads.
            kwargs = dict(kwargs)
            kwargs.update({
                "restart": True,
                "calcpsf": False,
                "calcres": False,
            })
            self.clear_stale_clean_mask(kwargs["imagename"])

        kwargs = self.adapt_automask_to_psf(kwargs, label)
        kwargs = self.cap_masked_major_cycles(kwargs)
        summary = self.tclean_pass(kwargs, label, "auto-multithresh")
        diagnostics = self.deconvolution_diagnostics(
            kwargs["imagename"]
        )
        if bootstrap is not None:
            diagnostics["bootstrap"] = bootstrap
        diagnostics = self.continue_to_convergence(
            kwargs, label, summary, diagnostics
        )
        diagnostics = self.floor_clean(kwargs, label, diagnostics)
        diagnostics["recovery_attempted"] = False
        diagnostics["recovery_succeeded"] = False
        peak_snr = diagnostics.get("residual_peak_snr")
        retry_min_snr = significance_floor(
            self.task_config.get("automask_retry_min_peak_snr", 5.0),
            diagnostics,
        )
        rms = diagnostics.get("residual_rms_Jy_beam")
        peak = diagnostics.get("residual_peak_abs_Jy_beam")
        print(
            f"Deconvolution of {label}: model "
            + ("present" if diagnostics.get("model_nonzero") else "EMPTY")
            + (f", residual RMS {rms:.3e} Jy/beam" if rms else "")
            + (f", peak {peak:.3e} Jy/beam" if peak else "")
            + (f", peak/RMS {peak_snr:.1f}" if peak_snr else "")
        )
        should_retry = (
            bool(
                self.task_config.get(
                    "automask_retry_enabled",
                    True,
                )
            ) and
            not diagnostics.get("model_nonzero") and
            peak_snr is not None and
            peak_snr >= retry_min_snr
        )

        if should_retry:
            print(
                "WARNING: Self-cal CLEAN model is empty while the residual "
                f"contains a {peak_snr:.1f}-sigma feature. Retrying once "
                "with the bounded recovery mask."
            )
            self.clear_stale_clean_mask(kwargs["imagename"])
            retry_kwargs = dict(kwargs)
            retry_kwargs.update(
                self.automask_kwargs(relaxed=True)
            )
            retry_kwargs.update({
                "restart": True,
                "calcpsf": False,
                "calcres": True,
            })
            initial = dict(diagnostics)
            self.tclean_pass(retry_kwargs, label, "relaxed-mask retry")
            diagnostics = self.deconvolution_diagnostics(
                retry_kwargs["imagename"]
            )
            diagnostics["initial"] = initial
            diagnostics["recovery_attempted"] = True
            diagnostics["recovery_succeeded"] = bool(
                diagnostics.get("model_nonzero")
            )
            print(
                "The relaxed-mask retry "
                + (
                    "produced a model."
                    if diagnostics["recovery_succeeded"] else
                    "still produced no model."
                )
            )

        elif not diagnostics.get("model_nonzero"):
            print(
                f"WARNING: CLEAN produced no model for {label}. The residual "
                "peak is "
                + (
                    f"{peak_snr:.1f} sigma, below the {retry_min_snr:.1f} "
                    "sigma the bounded relaxed-mask retry needs"
                    if peak_snr is not None else
                    "not measurable"
                )
                + ", so there is nothing here to deconvolve. The image is "
                "the dirty map, and self-calibration has no model to solve "
                "against."
            )


        if (
            not diagnostics.get("model_nonzero") and
            bool(self.task_config.get("automask_fallback_enabled", True)) and
            peak_snr is not None and
            peak_snr >= significance_floor(
                self.task_config.get("automask_fallback_min_peak_snr", 5.0),
                diagnostics,
            ) and
            self.primary_beam_image_ready(kwargs["imagename"])
        ):
            # auto-multithresh never opened a mask. Clean what is plainly
            # there instead: compact emission inside the primary beam, to a
            # conservative depth with a capped iteration count, which cannot
            # follow sidelobes into the noise.
            fallback_nsigma = float(self.task_config.get("automask_fallback_nsigma", 6.0))
            fallback_niter = int(self.task_config.get("automask_fallback_niter", 5000))
            print(
                "No mask was opened for "
                f"{label}, and the residual holds a {peak_snr:.1f}-sigma "
                f"feature. Cleaning inside the primary beam to "
                f"{fallback_nsigma:.1f} sigma, at most {fallback_niter} "
                "iterations."
            )
            self.clear_stale_clean_mask(kwargs["imagename"])
            fallback_kwargs = dict(kwargs)
            fallback_kwargs.update({
                "usemask": "pb",
                "pbmask": float(self.task_config.get("automask_fallback_pbmask", 0.2)),
                "nsigma": fallback_nsigma,
                "niter": fallback_niter,
                "restart": True,
                "calcpsf": False,
                "calcres": True,
            })
            before_fallback = dict(diagnostics)
            self.tclean_pass(fallback_kwargs, label, "primary-beam fallback")
            diagnostics = self.deconvolution_diagnostics(
                fallback_kwargs["imagename"]
            )
            diagnostics["before_fallback"] = before_fallback
            diagnostics["fallback_attempted"] = True
            diagnostics["fallback_succeeded"] = bool(
                diagnostics.get("model_nonzero")
            )
            print(
                "The primary-beam clean "
                + (
                    "produced a model."
                    if diagnostics["fallback_succeeded"] else
                    "still produced no model; the data hold nothing this "
                    "stage can deconvolve."
                )
            )

        if not hasattr(self, "clean_diagnostics"):
            self.clean_diagnostics = {}
        self.clean_diagnostics[label] = diagnostics
        return diagnostics

    # --------------------------------------------------------
    # Self-calibration model and solve diagnostics
    # --------------------------------------------------------

    def model_column_status(self, ms):
        """Check that MODEL_DATA is present, finite and non-zero throughout.

        Reads are sized by value count and spread over the whole Measurement
        Set, first and last rows included. A kill while tclean writes
        MODEL_DATA can leave the early rows filled and later rows zero, and
        a check that stops at the first non-zero value near the start of the
        table would accept that model and solve the later times against
        zero. Every sampled block holding unflagged data must therefore
        contain a non-zero model value.
        """

        opened = False
        checked = 0
        usable = 0
        max_rows = max(1, int(self.task_config.get("model_check_rows", 20000)))
        epsilon = float(
            self.task_config.get("model_nonzero_epsilon", 1.0e-12)
        )

        try:
            tb.open(str(ms))
            opened = True
            columns = set(tb.colnames())

            if "MODEL_DATA" not in columns:
                return False, "MODEL_DATA is absent"

            nrows = int(tb.nrows())
            if nrows < 1:
                return False, "the Measurement Set contains zero rows"

            probe = np.asarray(
                tb.getcol("MODEL_DATA", startrow=0, nrow=1)
            )
            values_per_row = max(1, int(probe.size))
            row_budget = max(
                1,
                min(max_rows, MODEL_CHECK_MAX_VALUES // values_per_row),
            )
            rows_per_block = max(
                1,
                min(row_budget, MODEL_CHECK_BLOCK_VALUES // values_per_row),
            )

            if nrows <= row_budget:
                starts = list(range(0, nrows, rows_per_block))
            else:
                blocks = max(2, -(-row_budget // rows_per_block))
                last_start = nrows - rows_per_block
                starts = sorted({
                    round(index * last_start / (blocks - 1))
                    for index in range(blocks)
                })

            for startrow in starts:
                nrow = min(rows_per_block, nrows - startrow)
                rows = f"rows {startrow}-{startrow + nrow - 1}"
                model = np.asarray(
                    tb.getcol("MODEL_DATA", startrow=startrow, nrow=nrow)
                )
                finite = np.isfinite(model)
                unflagged = np.ones(model.shape, dtype=bool)

                if "FLAG" in columns:
                    flags = np.asarray(
                        tb.getcol("FLAG", startrow=startrow, nrow=nrow),
                        dtype=bool,
                    )
                    if flags.shape == model.shape:
                        unflagged = ~flags

                if np.any(unflagged & ~finite):
                    return False, (
                        f"MODEL_DATA has non-finite unflagged values in {rows}"
                    )

                valid = unflagged & finite
                block_checked = int(np.count_nonzero(valid))
                if block_checked == 0:
                    continue

                block_usable = int(
                    np.count_nonzero(valid & (np.abs(model) > epsilon))
                )
                if block_usable == 0:
                    return False, (
                        f"MODEL_DATA is zero over the unflagged data in "
                        f"{rows}; the model column appears incompletely "
                        "written"
                    )

                checked += block_checked
                usable += block_usable

        except Exception as exc:
            return False, f"MODEL_DATA could not be inspected ({exc})"

        finally:
            if opened:
                tb.close()

        if checked < 1:
            return False, "no unflagged model samples were found"

        if usable < 1:
            return False, "MODEL_DATA is zero over the inspected unflagged data"

        return True, ""

    def predict_model_column(self, ms, image_base):
        """Restore an existing MT-MFS image model into an MS MODEL_DATA column."""

        image_base = str(image_base)
        model_tt0 = Path(image_base + ".model.tt0")

        if not model_tt0.is_dir():
            raise RuntimeError(
                "Cannot rebuild MODEL_DATA because the model image is missing: "
                f"{model_tt0}"
            )

        print(
            "Rebuilding MODEL_DATA from "
            f"{model_tt0.name} before gaincal."
        )

        # CASA documents this restart/predict-only pattern for writing an
        # existing tclean model into MODEL_DATA without another minor cycle.
        tclean(
            vis=str(ms),
            imagename=image_base,
            specmode="mfs",
            deconvolver="mtmfs",
            nterms=self.image_config["nterms"],
            gridder="wproject",
            wprojplanes=int(self.image_config.get("wprojplanes", 128)),
            pblimit=float(self.task_config.get("pblimit", 0.05)),
            weighting="briggs",
            robust=self.image_config["robust"],
            imsize=self.image_config["imsize"],
            cell=self.image_config["cell"],
            niter=0,
            savemodel="modelcolumn",
            calcpsf=False,
            calcres=False,
            restart=True,
            interactive=False,
        )

        valid, reason = self.model_column_status(ms)
        if not valid:
            raise RuntimeError(
                "tclean did not produce a usable MODEL_DATA column in "
                f"{ms}: {reason}"
            )

    def ensure_model_column(self, ms, image_base):
        """Ensure gaincal sees a valid model, including after split/restart."""

        valid, reason = self.model_column_status(ms)
        if valid:
            return

        print(f"MODEL_DATA requires repair: {reason}.")
        self.predict_model_column(ms, image_base)

    def measure_antenna_data(self, ms):
        """Count unflagged data per antenna and find the unsolvable ones.

        A gain can only be solved for an antenna that has data. An antenna
        that is dead, was quarantined during calibration, or lost its
        calibration upstream arrives here fully or almost fully flagged, and
        gaincal can only return flagged solutions for it. Counting those as
        failures of the solve let one such antenna reject the solution of
        every other antenna, so no target with a flagged antenna was ever
        self-calibrated. They are identified from the data, before the
        solve, and left out of the table's coverage test.
        """

        record = {
            "measurement_set": str(ms),
            "status": "unavailable",
            "antennas": {},
            "unsolvable_ids": [],
            "unsolvable_names": [],
        }
        self.antenna_data = record
        self.unsolvable_antenna_ids = set()

        try:
            summary = flagdata(
                vis=str(ms),
                mode="summary",
                action="calculate",
            )
        except Exception as exc:
            record["reason"] = f"flag summary failed: {exc}"
            return record

        per_antenna = (
            summary.get("antenna") if isinstance(summary, dict) else None
        )
        if not isinstance(per_antenna, dict) or not per_antenna:
            record["reason"] = "the flag summary holds no antenna counts"
            return record

        opened = False
        try:
            tb.open(str(Path(str(ms)) / "ANTENNA"))
            opened = True
            names = [str(item) for item in tb.getcol("NAME")]
        except Exception as exc:
            record["reason"] = f"the ANTENNA table could not be read: {exc}"
            return record
        finally:
            if opened:
                tb.close()

        unflagged = {}
        for name, counts in per_antenna.items():
            try:
                total = float(counts.get("total", 0.0))
                flagged = float(counts.get("flagged", 0.0))
            except (AttributeError, TypeError, ValueError):
                continue
            if math.isfinite(total) and math.isfinite(flagged):
                unflagged[str(name)] = max(0.0, total - flagged)

        populated = [value for value in unflagged.values() if value > 0]
        if not populated:
            record["reason"] = "no antenna has unflagged data"
            return record

        median = float(np.median(populated))
        minimum_fraction = float(
            self.task_config.get("min_antenna_data_fraction", 0.20)
        )
        threshold = minimum_fraction * median

        unsolvable_ids = []
        for index, name in enumerate(names):
            # An antenna absent from the summary has no rows at all: CASA
            # still writes a flagged placeholder solution for it.
            count = unflagged.get(name, 0.0)
            record["antennas"][name] = {
                "id": index,
                "unflagged_samples": count,
                "solvable": bool(count > threshold),
            }
            if count <= threshold:
                unsolvable_ids.append(index)

        record.update({
            "status": "measured",
            "median_unflagged_samples": median,
            "min_antenna_data_fraction": minimum_fraction,
            "unsolvable_ids": unsolvable_ids,
            "unsolvable_names": [names[index] for index in unsolvable_ids],
        })
        self.unsolvable_antenna_ids = set(unsolvable_ids)

        solvable = len(names) - len(unsolvable_ids)
        print(
            f"Antennas with data to solve from: {solvable} of {len(names)}."
            + (
                " Left out of the solution-coverage test (less than "
                f"{minimum_fraction:.0%} of the median antenna's unflagged "
                "data): " + ", ".join(record["unsolvable_names"]) + "."
                if unsolvable_ids else ""
            )
        )
        return record

    def ensure_antenna_census(self, ms):
        """Measure which antennas hold data once, and keep the answer.

        Every validation of a self-calibration table, including the one a
        resumed run makes of tables it accepted earlier, must leave out the
        same antennas. The census is therefore written next to the reports
        and read back on resume instead of being measured again on a
        Measurement Set whose flags a later step may have changed.
        """

        if getattr(self, "antenna_data", None) is not None:
            return self.antenna_data

        reportdir = getattr(self, "reportdir", None)
        saved = Path(reportdir) / "antenna_data.json" if reportdir else None
        if saved is not None and saved.is_file():
            try:
                with open(saved) as handle:
                    record = json.load(handle)
            except (OSError, ValueError, TypeError):
                record = None
            if (
                isinstance(record, dict) and
                record.get("status") == "measured"
            ):
                self.antenna_data = record
                self.unsolvable_antenna_ids = {
                    int(value) for value in record.get("unsolvable_ids", [])
                }
                print(
                    "Antennas left out of the solution-coverage test (from "
                    f"{saved.name}): "
                    + (", ".join(record.get("unsolvable_names", [])) or "none")
                    + "."
                )
                return record

        record = self.measure_antenna_data(ms)
        if saved is not None and record.get("status") == "measured":
            try:
                self.save_json(record, saved)
            except Exception as exc:
                print(f"WARNING: could not write {saved.name}: {exc}")
        return record

    def flag_unsolvable_antenna_data(self, ms):
        """Flag what is left of the antennas self-calibration cannot solve.

        An antenna with a small remnant of data gets no self-calibration
        solutions, and with ``calonly`` its visibilities are passed through
        with the phases the calibrator transfer gave them. Where phase
        self-calibration is worth doing those phases are poor (on the uGMRT
        field that prompted this the first loop raised the peak by 70%), so
        the remnant, about 3% of the visibilities there, is imaged with a
        different and worse calibration than everything around it. It is
        flagged in the candidate Measurement Set of a loop, and so only
        takes effect if that loop is accepted.
        """

        if not bool(self.task_config.get("flag_unsolvable_antennas", True)):
            return None

        census = getattr(self, "antenna_data", None) or {}
        if census.get("status") != "measured":
            return None

        remnants = {
            name: float(entry.get("unflagged_samples") or 0.0)
            for name, entry in (census.get("antennas") or {}).items()
            if not entry.get("solvable") and
            float(entry.get("unflagged_samples") or 0.0) > 0
        }
        if not remnants:
            return None

        total = sum(
            float(entry.get("unflagged_samples") or 0.0)
            for entry in census["antennas"].values()
        )
        names = sorted(remnants)
        print(
            "Flagging the remaining data of antennas that cannot be "
            f"self-calibrated: {', '.join(names)}."
        )
        flagdata(
            vis=str(ms),
            mode="manual",
            antenna=",".join(names),
            flagbackup=False,
        )
        self.unsolvable_antenna_flagging = {
            "antennas": names,
            # Each visibility is counted at both of its antennas.
            "fraction_of_antenna_samples": (
                sum(remnants.values()) / total if total > 0 else None
            ),
        }
        self.keep_record(
            self.unsolvable_antenna_flagging,
            "unsolvable_antenna_flagging.json",
        )
        return self.unsolvable_antenna_flagging

    def phase_reference_antennas(self, ms):
        """Return a prioritized active-reference-antenna list for gaincal."""

        census = getattr(self, "antenna_data", None) or {}
        if census.get("status") == "measured":
            # Rank by unflagged data, not by row count: a flagged antenna
            # has as many rows as a healthy one and makes a poor reference.
            usable = {
                name: entry["unflagged_samples"]
                for name, entry in census["antennas"].items()
                if entry.get("solvable")
            }
            ranked = sorted(usable, key=usable.get, reverse=True)
            preferred = str(self.reference_antenna).strip()
            ordered = [preferred] if preferred in usable else []
            for name in ranked:
                if name not in ordered:
                    ordered.append(name)
            if ordered:
                return ",".join(ordered[:5])

        opened = False
        try:
            tb.open(str(ms))
            opened = True
            ant1 = np.asarray(tb.getcol("ANTENNA1"), dtype=int)
            ant2 = np.asarray(tb.getcol("ANTENNA2"), dtype=int)
        except Exception:
            return str(self.reference_antenna)
        finally:
            if opened:
                tb.close()

        counts = {}
        for antenna_id in np.concatenate((ant1, ant2)):
            antenna_id = int(antenna_id)
            counts[antenna_id] = counts.get(antenna_id, 0) + 1

        ant_table = Path(ms) / "ANTENNA"
        opened = False
        try:
            tb.open(str(ant_table))
            opened = True
            names = [str(item) for item in tb.getcol("NAME")]
        except Exception:
            return str(self.reference_antenna)
        finally:
            if opened:
                tb.close()

        ranked = [
            names[index]
            for index, _ in sorted(
                counts.items(), key=lambda item: item[1], reverse=True
            )
            if 0 <= index < len(names)
        ]

        preferred = str(self.reference_antenna).strip()
        ordered = []
        for name in [preferred] + ranked:
            if name and name not in ordered:
                ordered.append(name)

        return ",".join(ordered[:5])

    def solve_phase_caltable(self, loop, solint, caltable):
        """Run the requested phase solve and one conservative rescue solve."""

        primary_minsnr = float(
            self.band_config.get(
                "min_snr_phase", self.task_config["gain_minsnr"]
            )
        )
        attempts = [
            {
                "name": "configured",
                "solint": solint,
                "combine": "",
                "minsnr": primary_minsnr,
                "minblperant": int(self.task_config["minblperant"]),
                "gaintype": self.task_config["phase_gaintype"],
            }
        ]

        if self.task_config.get("phase_retry_enabled", True):
            rescue = {
                "name": "long-interval rescue",
                "solint": self.task_config.get("phase_retry_solint", "inf"),
                "combine": self.task_config.get(
                    "phase_retry_combine", "scan"
                ),
                "minsnr": float(
                    self.task_config.get("phase_retry_minsnr", 3.0)
                ),
                "minblperant": int(
                    self.task_config.get("phase_retry_minblperant", 3)
                ),
                "gaintype": self.task_config.get(
                    "phase_retry_gaintype", "T"
                ),
            }
            signature = tuple(rescue[key] for key in (
                "solint", "combine", "minsnr", "minblperant", "gaintype"
            ))
            primary_signature = tuple(attempts[0][key] for key in (
                "solint", "combine", "minsnr", "minblperant", "gaintype"
            ))
            if signature != primary_signature:
                attempts.append(rescue)

        # Flags do not change between phase loops (solutions are applied
        # with calonly), so one census of the input serves every loop.
        self.ensure_antenna_census(self.current_ms)

        refants = self.phase_reference_antennas(self.current_ms)
        records = []

        for attempt_number, attempt in enumerate(attempts, start=1):
            if Path(caltable).exists():
                self.archive_incomplete_caltable(caltable)

            print(
                "Running phase gain solve "
                f"{attempt_number}/{len(attempts)}: "
                f"solint={attempt['solint']}, combine={attempt['combine']!r}, "
                f"minsnr={attempt['minsnr']}, "
                f"minblperant={attempt['minblperant']}, "
                f"gaintype={attempt['gaintype']}."
            )

            try:
                gaincal(
                    vis=self.current_ms,
                    caltable=caltable,
                    calmode="p",
                    solint=attempt["solint"],
                    combine=attempt["combine"],
                    refant=refants,
                    refantmode=self.task_config.get(
                        "phase_refantmode", "flex"
                    ),
                    minsnr=attempt["minsnr"],
                    minblperant=attempt["minblperant"],
                    gaintype=attempt["gaintype"],
                    solmode=self.task_config.get("phase_solmode", ""),
                    solnorm=False,
                    parang=self.task_config["parang"],
                )
                valid, reason = self.caltable_status(caltable)
            except Exception as exc:
                valid, reason = False, f"gaincal raised {exc}"

            record = dict(attempt)
            record.update({"valid": bool(valid), "reason": reason})
            records.append(record)

            if valid:
                return True, record, records

            print(
                "Phase solve attempt did not yield usable solutions: "
                f"{reason}"
            )

        return False, None, records

    # --------------------------------------------------------
    # Initial Imaging
    # --------------------------------------------------------

    def initial_imaging(self):

        if self.checkpoint_exists(
            "initial_image"
        ):

            print(
                "Initial image exists."
            )

            self.initial_image = str(
                self.imagedir /
                "initial.image.tt0"
            )    

            self.clean_diagnostics[
                "initial"
            ] = self.deconvolution_diagnostics(
                self.imagedir / "initial"
            )

            return

        imagename = str(
            self.imagedir /
            "initial"
        )

        self.archive_incomplete_images(imagename, "initial selfcal image")

        print(
            "\nGenerating initial image..."
        )

        kwargs = {
            "vis": self.current_ms,
            "imagename": imagename,
            "specmode": "mfs",
            "deconvolver": "mtmfs",
            "nterms": self.image_config["nterms"],
            "weighting": "briggs",
            "robust": self.image_config["robust"],
            "imsize": self.image_config["imsize"],
            "cell": self.image_config["cell"],
            "niter": self.image_config["niter"],
            "nsigma": self.image_config["nsigma"],
            "savemodel": "modelcolumn",
            "calcpsf": True,
            "calcres": True,
            "interactive": False,
        }
        kwargs.update(self.gridder_kwargs())
        kwargs.update(self.automask_kwargs())
        kwargs.update(self.clean_control_kwargs())
        self.run_tclean_with_recovery(
            kwargs,
            "initial",
        )

        self.initial_image = (
            imagename + ".image.tt0"
        )

        self.mark_checkpoint(
            "initial_image"
        )

    # --------------------------------------------------------
    # Estimate Image Statistics
    # --------------------------------------------------------

    def estimate_snr(self):

        print(
            "\nEstimating image statistics..."
        )

        stats = CASA_COMPAT.peak_and_noise(imstat, self.initial_image)

        peak = stats["max"]

        rms = stats["rms"]

        dynamic_range = stats["dynamic_range"] or 0.0

        self.image_stats = {

            "peak_Jy_beam": peak,
            "rms_Jy_beam": rms,
            "dynamic_range": dynamic_range
        }

        print(
            f"Peak       : {peak:.4e}"
        )

        print(
            f"RMS         : {rms:.4e}"
        )

        print(
            f"Dynamic Range: {dynamic_range:.1f}"
        )

        return dynamic_range

    # --------------------------------------------------------
    # Phase-only self calibration
    # --------------------------------------------------------

    def restore_phase_selfcal_state(self):
        """Restore accepted phase-loop state after a later-stage failure."""

        self.phase_tables = []
        self.phase_history = []
        previous_dr = self.image_stats.get("dynamic_range", 0.0)

        for loop, solint in enumerate(
            self.selfcal_strategy["phase_solints"],
            start=1,
        ):

            checkpoint = f"apply_phase{loop}"
            caltable = str(
                (self.caltabledir / f"phase_{loop}.cal").resolve()
            )

            if self.checkpoint_exists(checkpoint):

                split_ms = self.msdir / f"phase{loop}.ms"
                valid, reason = self.caltable_status(caltable)

                if not split_ms.is_dir():

                    raise RuntimeError(
                        "Accepted phase-selfcal checkpoint is missing its "
                        f"Measurement Set: {split_ms}"
                    )

                if not valid:

                    raise RuntimeError(
                        "Accepted phase-selfcal checkpoint has no usable "
                        f"calibration table: {caltable} ({reason})"
                    )

                current_dr = previous_dr
                image = self.imagedir / f"phase_{loop}.image.tt0"

                if image.exists():

                    stats = CASA_COMPAT.peak_and_noise(imstat, image)

                    if stats["dynamic_range"]:

                        current_dr = stats["dynamic_range"]

                improvement = (
                    (current_dr - previous_dr) / previous_dr
                    if previous_dr > 0
                    else 0.0
                )
                self.current_ms = str(split_ms)
                self.phase_tables.append(caltable)
                self.phase_history.append({
                    "loop": loop,
                    "solint": solint,
                    "dynamic_range": current_dr,
                    "improvement": improvement,
                    "measurement_set": self.current_ms,
                    "caltable": caltable,
                    "accepted": True,
                    "restored_from_checkpoint": True,
                })
                previous_dr = current_dr
                continue

            rejection_files = [
                self.rejection_path(f"apply_phase{loop}"),
                self.rejection_path(f"solve_phase{loop}"),
            ]

            for rejection_file in rejection_files:

                if rejection_file.is_file():

                    try:

                        with open(rejection_file) as handle:

                            record = json.load(handle)

                    except (OSError, ValueError, TypeError):

                        record = {
                            "loop": loop,
                            "solint": solint,
                            "accepted": False,
                            "rejection_reason": (
                                "Saved rejection record could not be read."
                            ),
                        }

                    record["restored_from_checkpoint"] = True
                    self.phase_history.append(record)
                    break

            break

        print(
            "Phase selfcal already completed; restored "
            f"{len(self.phase_tables)} accepted phase loop(s)."
        )

    # --------------------------------------------------------

    def run_phase_selfcal(self):

        if not self.selfcal_strategy["phase_selfcal"]:

            stats = getattr(self, "image_stats", {}) or {}
            threshold = float(self.task_config.get(
                "phase_start_dynamic_range",
                max(20.0, 3.0 * self.band_config["min_snr_phase"]),
            ))
            evidence = getattr(self, "phase_start_evidence", {}) or {}
            snr = evidence.get("expected_solution_snr")
            print(
                "\nSkipping phase selfcal: the initial image reaches a "
                f"dynamic range of {stats.get('dynamic_range', 0.0):.1f} "
                f"(peak {stats.get('peak_Jy_beam', 0.0):.3e} Jy/beam, RMS "
                f"{stats.get('rms_Jy_beam', 0.0):.3e} Jy/beam), below "
                f"{threshold:.1f}, and "
                + (
                    "each antenna's phase solution at "
                    f"solint={evidence.get('first_solint')} would reach only "
                    f"S/N {snr:.1f}, below "
                    f"{evidence.get('min_solution_snr', 10.0):.0f}"
                    if snr is not None else
                    "the solution S/N could not be estimated (no CLEAN model "
                    "or no visibility-noise measurement)"
                )
                + ". Solving against a model this weak would fit noise, so "
                "the calibrated data are imaged as they are."
            )

            self.phase_tables = []

            return

        # Before any table is validated, on a fresh run or a resumed one.
        self.ensure_antenna_census(self.current_ms)

        if self.checkpoint_exists("phase_selfcal"):

            self.restore_phase_selfcal_state()
            return

        self.phase_tables = []

        previous_dr = self.image_stats["dynamic_range"]

        self.phase_history = []

        for loop, solint in enumerate(
            self.selfcal_strategy["phase_solints"],
            start=1
        ):

            print(
                f"\nPhase selfcal loop {loop}"
            )

            caltable = str(
                (self.caltabledir / f"phase_{loop}.cal").resolve()
            )

            phase_checkpoint = f"apply_phase{loop}"

            if self.phase_rejected(loop):

                print(
                    "A prior candidate solution was rejected; "
                    "keeping the last accepted self-calibration MS."
                )

                break

            if self.checkpoint_exists(phase_checkpoint):

                split_ms = self.msdir / f"phase{loop}.ms"

                if not split_ms.exists():

                    raise RuntimeError(
                        "Accepted phase-selfcal checkpoint is missing "
                        f"its Measurement Set: {split_ms}"
                    )

                valid, reason = self.caltable_status(caltable)

                if not valid:

                    raise RuntimeError(
                        "Accepted phase-selfcal checkpoint has no usable "
                        "calibration table:\n"
                        f"{caltable}\n"
                        f"Reason: {reason}\n"
                        "The accepted MS cannot be reproduced safely. "
                        "Remove this phase checkpoint and its phase MS, "
                        "then rerun selfcal."
                    )

                self.current_ms = str(split_ms)
                self.phase_tables.append(caltable)

                image = self.imagedir / f"phase_{loop}.image.tt0"

                if image.exists():

                    stats = CASA_COMPAT.peak_and_noise(imstat, image)

                    if stats["dynamic_range"]:

                        previous_dr = stats["dynamic_range"]

                continue

            # Nothing from an uncheckpointed attempt at this loop is trusted.
            self.recover_interrupted_phase_loop(loop)

            previous_ms = self.current_ms

            model_base = (
                self.imagedir / "initial"
                if loop == 1
                else self.imagedir / f"phase_{loop - 1}"
            )
            self.ensure_model_column(self.current_ms, model_base)

            valid, accepted_attempt, solve_attempts = (
                self.solve_phase_caltable(loop, solint, caltable)
            )

            if not valid:

                reasons = "; ".join(
                    f"{item['name']}: {item['reason']}"
                    for item in solve_attempts
                )
                rejected = {
                    "loop": loop,
                    "solint": solint,
                    "caltable": caltable,
                    "measurement_set": previous_ms,
                    "accepted": False,
                    "solve_attempts": solve_attempts,
                    "rejection_reason": (
                        "The target/model did not support a reliable phase "
                        f"solution after all configured attempts ({reasons})."
                    ),
                }
                self.phase_history.append(rejected)
                self.save_json(
                    rejected,
                    self.rejection_path(f"solve_phase{loop}"),
                )

                print(
                    "No reliable phase-selfcal solution is available after "
                    "the configured and rescue solves. The last accepted "
                    "Measurement Set will be retained; this is a normal "
                    "low-S/N stopping condition."
                )

                break

            self.apply_phase_selfcal(loop, caltable)

            current_dr = self.image_after_phase_selfcal(
                loop
            )

            if previous_dr > 0:
                improvement = (current_dr - previous_dr) / previous_dr
            else:
                improvement = 0.0

            print(
                f"Dynamic range improvement = "
                f"{100*improvement:.1f}%"
            )

            self.phase_history.append({

                "loop": loop,

                "solint": solint,

                "dynamic_range": current_dr,

                "improvement": improvement,

                "measurement_set": self.current_ms,

                "caltable": caltable,

                "solve_parameters": accepted_attempt,

                "accepted": False,

            })

            if improvement < self.task_config[
                "phase_improvement_stop"
            ]:

                rejected = self.phase_history[-1]
                rejected["rejection_reason"] = (
                    "Dynamic-range improvement did not meet the "
                    "configured acceptance threshold."
                )
                rejected["candidate_measurement_set"] = self.current_ms

                self.discard_phase_candidate(loop)
                self.current_ms = previous_ms
                rejected["measurement_set"] = self.current_ms
                self.mark_phase_rejected(loop, rejected)

                print(
                    "Rejecting phase-selfcal candidate and retaining "
                    "the last accepted Measurement Set."
                )

                break

            self.commit_phase_candidate(loop)
            self.phase_history[-1]["accepted"] = True
            self.phase_history[-1]["measurement_set"] = self.current_ms
            self.phase_tables.append(caltable)
            self.mark_checkpoint(phase_checkpoint)

            previous_dr = current_dr

        self.mark_checkpoint("phase_selfcal")

    # --------------------------------------------------------
    # Apply Phase Selfcal Solutions and Split
    # --------------------------------------------------------

    def apply_phase_selfcal(self, loop, caltable):

        phase_checkpoint = f"apply_phase{loop}"

        if self.checkpoint_exists(phase_checkpoint):

            split_ms = self.msdir / f"phase{loop}.ms"

            if split_ms.exists():
                self.current_ms = str(split_ms)

            else:

                raise RuntimeError(
                    "Phase-selfcal checkpoint is missing its "
                    f"Measurement Set: {split_ms}"
                )

            print(
                f"Apply phase selfcal already done for loop {loop}."
            )

            return str(split_ms)

        valid, reason = self.caltable_status(caltable)

        if not valid:

            raise RuntimeError(
                "Refusing to apply an unusable phase-selfcal table:\n"
                f"{caltable}\n"
                f"Reason: {reason}"
            )

        split_ms = self.msdir / f"phase{loop}.ms"
        candidate_ms = self.msdir / f"phase{loop}.candidate.ms"

        if split_ms.exists():

            raise RuntimeError(
                "An unaccepted phase-selfcal candidate already exists:\n"
                f"{split_ms}\n"
                "Review or remove the candidate before rerunning this "
                "self-calibration loop."
            )

        if candidate_ms.exists():

            raise RuntimeError(
                "A temporary phase-selfcal candidate already exists:\n"
                f"{candidate_ms}\n"
                "Review or remove it before rerunning this "
                "self-calibration loop."
            )

        # Apply candidate solutions to a temporary copy. This preserves the
        # previously accepted MS if the dynamic-range test rejects the loop.
        print(
            f"\nCreating temporary candidate {candidate_ms.name}"
        )

        split(
            vis=self.current_ms,
            outputvis=str(candidate_ms),
            datacolumn="all",
            keepflags=True,
        )

        print(
            f"Applying {Path(caltable).name} to {candidate_ms.name}"
        )

        try:
            # In the candidate only: if the loop is rejected these flags go
            # with it, and the data are imaged as calibration left them.
            self.flag_unsolvable_antenna_data(candidate_ms)
            applycal(
                vis=str(candidate_ms),
                gaintable=[caltable],
                interp=["linear"],
                applymode=self.task_config["applymode"],
                calwt=self.task_config["calwt"],
                parang=self.task_config["parang"],
            )
        except Exception:
            self.discard_phase_candidate(loop)
            raise

        opened = False
        try:
            tb.open(str(candidate_ms))
            opened = True
            cols = set(tb.colnames())
            rows = int(tb.nrows())
        finally:
            if opened:
                tb.close()

        if rows < 1 or "CORRECTED_DATA" not in cols:
            self.discard_phase_candidate(loop)
            raise RuntimeError(
                "applycal did not create a usable CORRECTED_DATA column "
                "in the temporary phase-selfcal candidate."
            )

        self.current_ms = str(candidate_ms)

        return self.current_ms

    # --------------------------------------------------------
    # Image after phase selfcal
    # --------------------------------------------------------

    def image_candidate(self, label):
        """Image the current (candidate) Measurement Set and measure it.

        Every self-calibration loop is judged on an image made the same way:
        same grid, weighting, masks and depth, with the model written to
        MODEL_DATA for the solve that follows if the loop is accepted.
        """

        imagename = str(self.imagedir / label)

        kwargs = {
            "vis": self.current_ms,
            "imagename": imagename,
            "specmode": "mfs",
            "deconvolver": "mtmfs",
            "nterms": self.image_config["nterms"],
            "weighting": "briggs",
            "robust": self.image_config["robust"],
            "imsize": self.image_config["imsize"],
            "cell": self.image_config["cell"],
            "niter": self.image_config["niter"],
            "nsigma": self.image_config["nsigma"],
            "savemodel": "modelcolumn",
            "calcpsf": True,
            "calcres": True,
            "interactive": False,
        }
        kwargs.update(self.gridder_kwargs())
        kwargs.update(
            self.automask_kwargs(
                self.initial_image
            )
        )
        kwargs.update(self.clean_control_kwargs())
        self.run_tclean_with_recovery(
            kwargs,
            label,
        )

        return CASA_COMPAT.peak_and_noise(
            imstat, imagename + ".image.tt0"
        )

    def image_after_phase_selfcal(
        self,
        loop
    ):

        stats = self.image_candidate(f"phase_{loop}")

        dynamic_range = stats["dynamic_range"] or 0.0

        print(
            f"Dynamic Range = "
            f"{dynamic_range:.0f}"
        )

        return dynamic_range

    # --------------------------------------------------------
    # Amplitude + Phase Selfcal
    # --------------------------------------------------------

    def run_amp_phase_selfcal(self):

        if not self.task_config.get("enable_amp_selfcal", False):

            print(
                "\nSkipping amp+phase selfcal (disabled by "
                "configuration; it can change the flux scale)."
            )

            return

        amp_solint = (
            self.selfcal_strategy.get("amp_solint") or
            self.band_config.get("amp_solint")
        )

        if not amp_solint:

            print(
                "\nSkipping amp+phase selfcal: no amplitude solution "
                "interval is configured for this band."
            )

            return

        if not getattr(self, "phase_tables", []):

            print(
                "\nSkipping amp+phase selfcal because no phase-only "
                "solution was accepted."
            )

            return

        checkpoint = "amp_selfcal"

        caltable = str(
            (self.caltabledir / "amp_phase.cal").resolve()
        )

        if self.checkpoint_exists(checkpoint):

            amp_ms = (
                self.msdir /
                "final_selfcal.ms"
            )

            if amp_ms.exists():

                self.current_ms = str(amp_ms)
                self.amp_table = caltable
                self.amp_history = self.read_json(
                    self.reportdir / "amp_selfcal.json"
                )

                return

            raise RuntimeError(
                "Amp selfcal checkpoint is missing its Measurement Set: "
                f"{amp_ms}"
            )

        rejection_file = self.rejection_path("amp_selfcal")

        if rejection_file.is_file():

            # Decided by an earlier attempt of this stage; the phase-only
            # Measurement Set it kept is the one in use.
            self.amp_history = self.read_json(rejection_file)
            print(
                "\nAmp+phase selfcal was rejected earlier; keeping the "
                "phase-only result."
            )

            return

        final_ms = self.msdir / "final_selfcal.ms"
        candidate_ms = self.msdir / "amp_phase.candidate.ms"

        # Left by a kill part-way through. The candidate is a copy and the
        # accepted phase Measurement Set was never written to, so redoing
        # the step from it is safe.
        self.archive_incomplete_product(
            final_ms,
            "amp+phase selfcal Measurement Set",
        )
        self.archive_incomplete_product(
            candidate_ms,
            "amp+phase selfcal candidate Measurement Set",
        )
        self.archive_incomplete_images(
            self.imagedir / "amp_phase",
            "amp+phase selfcal image",
        )

        print(
            "\nRunning amp+phase selfcal."
        )

        # commit_phase_candidate splits with datacolumn="corrected", which
        # does not carry MODEL_DATA.  Without this the amplitude solve would
        # run against CASA's implicit 1 Jy point source at the phase centre
        # and drive every antenna gain towards that fictitious model.
        last_loop = len(getattr(self, "phase_tables", []))
        reference_base = self.imagedir / (
            f"phase_{last_loop}" if last_loop else "initial"
        )
        self.ensure_model_column(
            self.current_ms,
            reference_base,
        )

        reference = CASA_COMPAT.peak_and_noise(
            imstat, str(reference_base) + ".image.tt0"
        )

        # Amplitude self-calibration absorbs sky-model error into the gains,
        # so it needs a well-modelled field. The test is made on the image
        # the solve starts from: a field can begin below the threshold and
        # pass it once its phases have been corrected.
        amp_start_dynamic_range = float(
            self.task_config.get("amp_start_dynamic_range", 100.0)
        )
        reference_dr = float(reference.get("dynamic_range") or 0.0)

        if reference_dr < amp_start_dynamic_range:

            self.amp_history = {
                "accepted": False,
                "attempted": False,
                "reference_dynamic_range": reference_dr,
                "rejection_reason": (
                    f"The phase-only image reaches a dynamic range of "
                    f"{reference_dr:.0f}, below the "
                    f"{amp_start_dynamic_range:.0f} a reliable amplitude "
                    "solution needs."
                ),
            }
            self.save_json(self.amp_history, rejection_file)

            print(
                "\nSkipping amp+phase selfcal: "
                f"{self.amp_history['rejection_reason']}"
            )

            return

        amp_minsnr = float(
            self.task_config.get("amp_minsnr", DEFAULT_AMP_SELFCAL_MINSNR)
        )
        amp_gaintype = str(
            self.task_config.get(
                "amp_solve_gaintype", DEFAULT_AMP_SELFCAL_GAINTYPE
            )
        )
        # The same priority list the phase loops use. With one reference
        # antenna and no alternate, CASA flags every antenna in any interval
        # where that antenna's own solution fails.
        refants = self.phase_reference_antennas(self.current_ms)
        amp_solve = {
            "solint": amp_solint,
            "minsnr": amp_minsnr,
            "gaintype": amp_gaintype,
            "refant": refants,
            "interp": AMP_SELFCAL_INTERP,
        }

        if not self.prepare_caltable(caltable):

            gaincal(

                vis=self.current_ms,

                caltable=caltable,

                calmode="ap",

                solint=amp_solint,

                refant=refants,

                refantmode=self.task_config.get(
                    "phase_refantmode", "flex"
                ),

                minsnr=amp_minsnr,

                minblperant=self.task_config[
                    "minblperant"
                ],

                gaintype=amp_gaintype,

                combine="",

                # Amplitudes are normalised to a median of one, so the solve
                # redistributes gain between antennas and times and leaves
                # the flux scale set by the calibrators where it was.
                solnorm=True,

                normtype="median",

                parang=self.task_config[
                    "parang"
                ]
            )

        valid, reason = self.caltable_status(caltable)

        if not valid:

            rejected = {
                "caltable": caltable,
                "measurement_set": self.current_ms,
                "accepted": False,
                "attempted": True,
                "solve_parameters": amp_solve,
                "reference_dynamic_range": reference_dr,
                "rejection_reason": (
                    "gaincal did not produce a usable CASA calibration "
                    f"table: {reason}"
                ),
            }
            self.save_json(
                rejected,
                self.rejection_path("solve_amp_phase"),
            )
            self.amp_history = rejected

            print(
                "WARNING: amp+phase selfcal solve produced no usable "
                f"calibration table ({reason}). Retaining the phase-only "
                "result."
            )

            return

        # The solutions go onto a copy, and the copy is imaged, before
        # anything is kept: an amplitude solve against an incomplete model
        # can lower the dynamic range or move the flux scale, and phase
        # loops are held to the same test.
        previous_ms = self.current_ms

        print(
            f"\nCreating temporary candidate {candidate_ms.name}"
        )

        split(
            vis=self.current_ms,
            outputvis=str(candidate_ms),
            datacolumn="all",
            keepflags=True,
        )

        try:
            applycal(
                vis=str(candidate_ms),
                gaintable=[caltable],
                interp=[AMP_SELFCAL_INTERP],
                applymode=self.task_config[
                    "applymode"
                ],
                calwt=self.task_config[
                    "calwt"
                ],
                parang=self.task_config[
                    "parang"
                ]
            )

            tb.open(str(candidate_ms))
            try:
                cols = tb.colnames()
            finally:
                tb.close()

            if "CORRECTED_DATA" not in cols:
                raise RuntimeError("CORRECTED_DATA column does not exist.")

            self.current_ms = str(candidate_ms)
            stats = self.image_candidate("amp_phase")

        except Exception:
            self.current_ms = previous_ms
            if candidate_ms.exists():
                shutil.rmtree(candidate_ms)
            raise

        record = self.judge_amp_selfcal(reference, stats)
        record.update({
            "solint": amp_solint,
            "solve_parameters": amp_solve,
            "caltable": caltable,
            "reference_image": str(reference_base) + ".image.tt0",
            "measurement_set": previous_ms,
        })
        self.amp_history = record

        print(
            "Amp+phase selfcal: dynamic range "
            f"{record['reference_dynamic_range']:.0f} -> "
            f"{record['dynamic_range']:.0f} "
            f"({100 * record['improvement']:+.1f}%), peak "
            f"{100 * record['peak_change']:+.1f}%."
        )

        if not record["accepted"]:

            shutil.rmtree(candidate_ms)
            self.current_ms = previous_ms
            self.save_json(record, rejection_file)

            print(
                "Rejecting the amp+phase candidate and retaining the "
                f"phase-only Measurement Set: {record['rejection_reason']}"
            )

            return

        split(

            vis=str(candidate_ms),

            outputvis=str(final_ms),

            datacolumn="corrected",

            keepflags=True
        )

        shutil.rmtree(candidate_ms)

        self.current_ms = str(final_ms)
        record["measurement_set"] = self.current_ms

        self.amp_table = caltable

        self.save_json(record, self.reportdir / "amp_selfcal.json")

        # A final image made before this solve, in an earlier attempt of the
        # stage, shows the phase-only data. final_imaging then moves it
        # aside and images the Measurement Set now in use.
        self.clear_checkpoint("final_image")

        self.mark_checkpoint(
            checkpoint
        )

    def judge_amp_selfcal(self, reference, stats):
        """Whether an amp+phase candidate image is better than what it replaces.

        Two conditions, both needed. The dynamic range has to rise by at
        least ``amp_improvement_stop`` (the phase loops' threshold unless set
        otherwise). And the peak may not move by more than
        ``amp_max_peak_change``: the solutions are normalised, so a peak that
        shifts by more than phase coherence explains means the solve has
        moved the flux scale to suit its model.
        """

        threshold = float(self.task_config.get(
            "amp_improvement_stop",
            self.task_config.get("phase_improvement_stop", 0.02),
        ))
        max_peak_change = float(
            self.task_config.get("amp_max_peak_change", 0.10)
        )

        reference_dr = float(reference.get("dynamic_range") or 0.0)
        reference_peak = float(reference.get("max") or 0.0)
        dynamic_range = float(stats.get("dynamic_range") or 0.0)
        peak = float(stats.get("max") or 0.0)

        improvement = (
            (dynamic_range - reference_dr) / reference_dr
            if reference_dr > 0 else 0.0
        )
        peak_change = (
            (peak - reference_peak) / reference_peak
            if reference_peak > 0 else 0.0
        )

        record = {
            "reference_dynamic_range": reference_dr,
            "reference_peak_Jy_beam": reference_peak,
            "dynamic_range": dynamic_range,
            "peak_Jy_beam": peak,
            "rms_Jy_beam": stats.get("rms"),
            "improvement": improvement,
            "peak_change": peak_change,
            "improvement_threshold": threshold,
            "max_peak_change": max_peak_change,
            "accepted": True,
        }

        if reference_dr <= 0 or reference_peak <= 0:
            record["accepted"] = False
            record["rejection_reason"] = (
                "The phase-only image could not be measured, so there is "
                "nothing to compare the candidate with."
            )
        elif improvement < threshold:
            record["accepted"] = False
            record["rejection_reason"] = (
                f"Dynamic range changed by {100 * improvement:+.1f}%, below "
                f"the {100 * threshold:.0f}% an amplitude solution has to "
                "earn."
            )
        elif abs(peak_change) > max_peak_change:
            record["accepted"] = False
            record["rejection_reason"] = (
                f"The peak moved by {100 * peak_change:+.1f}%, more than the "
                f"{100 * max_peak_change:.0f}% allowed: the flux scale would "
                "no longer be the calibrators'."
            )

        return record

    @staticmethod
    def read_json(path):
        try:
            with open(path) as handle:
                return json.load(handle)
        except (OSError, ValueError, TypeError):
            return None

    # --------------------------------------------------------
    # Final Imaging
    # --------------------------------------------------------

    def visibility_datacolumn(self, ms):
        """The column holding the best visibilities in this Measurement Set."""

        opened = False
        try:
            tb.open(str(ms))
            opened = True
            columns = set(tb.colnames())
        except Exception:
            return None
        finally:
            if opened:
                tb.close()

        if "CORRECTED_DATA" in columns:
            return "CORRECTED_DATA"
        return "DATA" if "DATA" in columns else None

    def measure_coherence(self, ms, rows_wanted=20000, block_rows=512):
        """Test whether the calibrated data add up coherently.

        Averaging the visibilities as complex numbers gives the flux at the
        phase centre; averaging their moduli gives, for a target field, the
        per-visibility noise, which is far larger. The ratio of the two is
        therefore small for any target field and says nothing on its own.
        What does say something is whether the vector average stands above
        its own uncertainty: with N samples of per-component noise sigma,
        that uncertainty is sigma/sqrt(N). A vector average many times larger
        can only come from visibilities whose phases agree, which is what
        calibrated data look like.

        The test is one-sided. A high significance shows the data are
        coherent; a low one can also mean a field with little flux at its
        phase centre, so it is reported as inconclusive rather than as a
        fault. Rows are sampled across the whole Measurement Set in bounded
        blocks.
        """

        column = self.visibility_datacolumn(ms)
        if column is None:
            return None

        opened = False
        try:
            tb.open(str(ms))
            opened = True
            nrows = int(tb.nrows())
            if nrows < 1:
                return None

            columns = set(tb.colnames())
            block_rows = max(1, min(int(block_rows), nrows))
            starts = list(range(0, nrows, block_rows))
            wanted_blocks = max(1, int(rows_wanted) // block_rows)
            if len(starts) > wanted_blocks:
                step = len(starts) / float(wanted_blocks)
                starts = sorted({
                    starts[min(len(starts) - 1, int(index * step))]
                    for index in range(wanted_blocks)
                })

            vector_sum = 0.0 + 0.0j
            scalar_sum = 0.0
            samples = 0
            autocorrelations = 0
            rows_read = 0
            correlations = None
            antennas = set()
            descriptions = set()
            scans = set()
            intervals = []
            weight_sum = 0.0
            weight_column = (
                "WEIGHT_SPECTRUM" if "WEIGHT_SPECTRUM" in columns
                else "WEIGHT" if "WEIGHT" in columns
                else None
            )

            for startrow in starts:
                nrow = min(block_rows, nrows - startrow)
                rows_read += nrow
                values = np.asarray(
                    tb.getcol(column, startrow=startrow, nrow=nrow)
                )
                if values.ndim >= 1:
                    correlations = int(values.shape[0])
                keep = np.isfinite(values)
                if "FLAG" in columns:
                    flags = np.asarray(
                        tb.getcol("FLAG", startrow=startrow, nrow=nrow),
                        dtype=bool,
                    )
                    if flags.shape == values.shape:
                        keep &= ~flags
                if "ANTENNA1" in columns and "ANTENNA2" in columns:
                    antenna1 = np.asarray(
                        tb.getcol("ANTENNA1", startrow=startrow, nrow=nrow),
                        dtype=int,
                    )
                    antenna2 = np.asarray(
                        tb.getcol("ANTENNA2", startrow=startrow, nrow=nrow),
                        dtype=int,
                    )
                    cross = antenna1 != antenna2
                    autocorrelations += int(np.count_nonzero(~cross))
                    keep &= cross.reshape((1,) * (values.ndim - 1) + (-1,))

                selected = values[keep]
                if selected.size == 0:
                    continue
                if weight_column:
                    try:
                        weights = np.asarray(
                            tb.getcol(weight_column, startrow=startrow,
                                      nrow=nrow),
                            dtype=float,
                        )
                        if weights.ndim == values.ndim - 1:
                            # WEIGHT is (correlation, row): every channel
                            # of a row shares it.
                            weights = np.broadcast_to(
                                weights[:, np.newaxis, :], values.shape
                            )
                        weight_sum += float(np.sum(weights[keep]))
                    except Exception:
                        weight_column = None
                vector_sum += complex(np.sum(selected))
                scalar_sum += float(np.sum(np.abs(selected)))
                samples += int(selected.size)

                # Which rows carried any unflagged data, for the shape of
                # the full dataset a gain solve will see.
                useful = keep.reshape(-1, keep.shape[-1]).any(axis=0)
                if "ANTENNA1" in columns and "ANTENNA2" in columns:
                    antennas.update(int(a) for a in antenna1[useful])
                    antennas.update(int(a) for a in antenna2[useful])
                for name, target in (
                    ("DATA_DESC_ID", descriptions),
                    ("SCAN_NUMBER", scans),
                ):
                    if name in columns:
                        target.update(
                            int(value) for value in np.asarray(
                                tb.getcol(name, startrow=startrow, nrow=nrow)
                            )[useful]
                        )
                if "INTERVAL" in columns:
                    intervals.extend(
                        float(value) for value in np.asarray(
                            tb.getcol("INTERVAL", startrow=startrow, nrow=nrow)
                        )[useful][:64]
                    )

            # Only the solvability estimate needs this; a failure here must
            # not cost the coherence measurement itself.
            try:
                integrations = self.count_integrations(nrows) \
                    if "TIME" in columns else None
            except Exception:
                integrations = None

        except Exception as exc:
            return {"status": "unavailable", "reason": str(exc)}
        finally:
            if opened:
                tb.close()

        if samples < 1 or scalar_sum <= 0.0:
            return {
                "status": "no_unflagged_data",
                "column": column,
                "samples": samples,
            }

        scalar_mean = scalar_sum / samples
        vector_mean = abs(vector_sum) / samples
        # For a complex Gaussian the mean modulus is sigma * sqrt(pi/2).
        component_sigma = scalar_mean / math.sqrt(math.pi / 2.0)
        vector_mean_error = component_sigma / math.sqrt(samples)

        # The rows read are a sample. The whole dataset holds proportionally
        # more unflagged samples, and its thermal noise is correspondingly
        # lower than the sample's uncertainty above.
        estimated_samples = samples * (nrows / float(max(1, rows_read)))
        interval = float(np.median(intervals)) if intervals else None
        on_source = (
            integrations * interval
            if integrations and interval else None
        )

        return {
            "status": "measured",
            "column": column,
            "samples": samples,
            "autocorrelations_skipped": autocorrelations,
            "scalar_mean_amplitude": scalar_mean,
            "vector_mean_amplitude": vector_mean,
            "component_sigma": component_sigma,
            "vector_mean_error": vector_mean_error,
            "significance": (
                vector_mean / vector_mean_error
                if vector_mean_error > 0 else None
            ),
            # Kept for reference. For a target field this is set by the
            # per-visibility noise, not by the calibration.
            "vector_to_scalar_ratio": abs(vector_sum) / scalar_sum,
            "rows": nrows,
            "rows_sampled": rows_read,
            "estimated_unflagged_samples": estimated_samples,
            # Naturally weighted, at the phase centre, from every unflagged
            # sample. Briggs weighting and w-projection raise it somewhat;
            # it is the floor the image is being compared against.
            "thermal_noise_Jy_beam": (
                component_sigma / math.sqrt(estimated_samples)
            ),
            "correlations": correlations,
            "antennas_with_data": len(antennas) or None,
            "spectral_windows": len(descriptions) or None,
            "scans": len(scans) or None,
            "integration_s": interval,
            "integrations": integrations,
            "on_source_s": on_source,
            # Only meaningful when the weights are inverse variances in
            # Jy^-2, as statwt writes them; the caller decides.
            "weight_column": weight_column,
            "estimated_weight_sum": (
                weight_sum * (nrows / float(max(1, rows_read)))
                if weight_column and weight_sum > 0 else None
            ),
        }

    @staticmethod
    def count_integrations(nrows, block_rows=1 << 20):
        """Distinct timestamps in the open MAIN table, read in blocks.

        Only the TIME column is read, a block at a time, so memory stays at
        one block however long the observation.
        """

        seen = np.empty(0)
        for startrow in range(0, nrows, block_rows):
            times = np.asarray(
                tb.getcol(
                    "TIME",
                    startrow=startrow,
                    nrow=min(block_rows, nrows - startrow),
                )
            )
            seen = np.union1d(seen, np.unique(times))
        return int(seen.size) or None

    def statistical_weights_applied(self):
        """Whether calibration wrote inverse-variance weights for this target."""

        record = (
            (getattr(self, "calibration", None) or {}).get(
                "statistical_weights"
            ) or {}
        ).get(getattr(self, "target", None)) or {}
        return record.get("status") == "applied"

    def adopt_statistical_thermal_noise(self, measurement):
        """Use the weights for the thermal floor once they mean Jy^-2.

        Without inverse-variance weights the floor is sigma/sqrt(N), which
        treats every sample as equally noisy. With statwt's weights the
        naturally weighted image noise is 1/sqrt(sum of weights) directly,
        and the channels the bandpass left noisy count for what they are
        worth. Both are kept; the second becomes the floor.
        """

        total = measurement.get("estimated_weight_sum")
        if not total or not self.statistical_weights_applied():
            return None
        measurement["thermal_noise_equal_weights_Jy_beam"] = measurement.get(
            "thermal_noise_Jy_beam"
        )
        measurement["thermal_noise_Jy_beam"] = 1.0 / math.sqrt(total)
        measurement["thermal_noise_source"] = "statistical_weights"
        return measurement["thermal_noise_Jy_beam"]

    def check_calibrated_coherence(self):
        """Report whether the calibrated target data are coherent.

        Diagnostic only: it never changes the data and never fails the stage.
        """

        try:
            measurement = self.measure_coherence(self.current_ms)
        except Exception as exc:
            measurement = {"status": "unavailable", "reason": str(exc)}

        self.coherence = measurement
        if not measurement:
            return None
        self.adopt_statistical_thermal_noise(measurement)

        print("\nCalibrated-data check")
        if measurement.get("status") != "measured":
            print(
                f"  not measured: {measurement.get('reason') or measurement.get('status')}"
            )
            return measurement

        floor = float(
            self.task_config.get("coherence_min_significance", 5.0)
        )
        significance = measurement.get("significance")
        print(
            f"  column {measurement['column']}, "
            f"{measurement['samples']} unflagged cross-hand samples"
        )
        print(
            "  per-visibility noise "
            f"{measurement['component_sigma']:.4e} Jy per component"
        )
        print(
            "  coherent flux at the phase centre "
            f"{measurement['vector_mean_amplitude']:.4e} Jy "
            f"+/- {measurement['vector_mean_error']:.2e} Jy"
            + (f" ({significance:.0f} sigma)" if significance else "")
        )
        thermal = measurement.get("thermal_noise_Jy_beam")
        if thermal:
            print(
                "  thermal noise of all unflagged data: about "
                f"{thermal * 1e6:.1f} uJy/beam (natural weighting; "
                "Briggs weighting raises it somewhat)"
                + (
                    ", from the statistical weights"
                    if measurement.get("thermal_noise_source") ==
                    "statistical_weights" else ""
                )
            )

        if significance is None or significance < floor:
            print(
                "  NOTE: the visibilities show no significant coherent flux "
                "at the phase centre. That happens both when calibration has "
                "not transferred and when the field simply has little flux "
                "at its centre, so this is inconclusive on its own. If the "
                "image also shows no sources, check the phase-calibrator "
                "solutions, the target-to-calibrator mapping and the "
                "bandpass.",
                flush=True,
            )
        else:
            print(
                "  The data are coherent: visibilities with random phases "
                "could not produce this.",
                flush=True,
            )
        return measurement

    # --------------------------------------------------------

    @staticmethod
    def summary_flagged_fraction(summary):
        """Flagged fraction from a flagdata summary, or None if unusable."""

        try:
            total = float(summary.get("total", 0.0))
            flagged = float(summary.get("flagged", 0.0))
        except (AttributeError, TypeError, ValueError):
            return None

        if (
            total <= 0 or
            flagged < 0 or
            flagged > total or
            not math.isfinite(total) or
            not math.isfinite(flagged)
        ):
            return None

        fraction = flagged / total
        return fraction if math.isfinite(fraction) else None

    # --------------------------------------------------------

    def residual_datacolumn(self, ms):
        """Which residual CASA can form for this Measurement Set.

        ``residual`` is CORRECTED minus MODEL. A Measurement Set split from an
        accepted self-calibration round carries its corrected visibilities in
        DATA and has no CORRECTED column; there the residual is DATA minus
        MODEL instead.
        """

        opened = False
        try:
            tb.open(str(ms))
            opened = True
            columns = set(tb.colnames())
        except Exception:
            return None
        finally:
            if opened:
                tb.close()

        if "MODEL_DATA" not in columns:
            return None
        return "residual" if "CORRECTED_DATA" in columns else "residual_data"

    # --------------------------------------------------------

    def sky_model_is_empty(self):
        """Whether the published image contains no clean components."""

        diagnostics = (getattr(self, "clean_diagnostics", {}) or {})
        for label in ("final", "initial"):
            entry = diagnostics.get(label)
            if isinstance(entry, dict) and "model_nonzero" in entry:
                return not entry["model_nonzero"]
        return False

    def warn_if_no_sky_model(self):
        """Say plainly when nothing was deconvolved, and what that means."""

        if not self.sky_model_is_empty():
            return False

        diagnostics = (
            (getattr(self, "clean_diagnostics", {}) or {}).get("final") or
            (getattr(self, "clean_diagnostics", {}) or {}).get("initial") or
            {}
        )
        rms = diagnostics.get("residual_rms_Jy_beam")
        peak = diagnostics.get("residual_peak_abs_Jy_beam")
        snr = diagnostics.get("residual_peak_snr")

        print("\n" + "=" * 70)
        print(f"NO SKY MODEL WAS DECONVOLVED FOR {self.target}")
        print("=" * 70)
        print(
            "CLEAN found nothing to clean, so the published image is the\n"
            "dirty map: sources are not deconvolved, sidelobes are not\n"
            "removed, and self-calibration had no model to solve against.\n"
            "Interference and sidelobe structure dominate such an image."
        )
        if rms:
            print(
                f"\nResidual RMS {rms:.3e} Jy/beam, peak "
                f"{peak:.3e} Jy/beam"
                + (f", peak/RMS {snr:.1f}" if snr else "")
            )
        print(
            "\nThis is a property of the data reaching imaging, not of the\n"
            "deconvolution settings. Check, in this order: the flagged\n"
            "fraction for this field, the calibration-table QA sidecars in\n"
            "calibration/, and whether the phase calibrator brackets this\n"
            "target in time. A field whose calibration did not transfer\n"
            "coherently produces exactly this result.",
            flush=True,
        )
        return True

    def record_skipped_residual_flagging(self, reason, diagnostics=None):
        """Record why the residual pass did not run."""

        record = {
            "status": "skipped",
            "reason": reason,
            "measurement_set": getattr(self, "current_ms", None),
        }
        if diagnostics is not None:
            record["deconvolution_diagnostics"] = diagnostics
        self.residual_flagging = record
        self.save_json(record, self.reportdir / "residual_flagging.json")

    # --------------------------------------------------------

    def latest_image_base(self):
        """The image whose model describes the current Measurement Set.

        That is the image of the last *accepted* loop. A rejected loop leaves
        its image behind too, made from a candidate whose solutions were
        thrown away; the data kept are those of the loop before it. Its
        model is not a model of them, and where the loop was rejected for
        making the image worse, subtracting it would leave real emission in
        the "residual" for the interference detector to flag.
        """

        # An accepted amp+phase solution is the last word on the data.
        if getattr(self, "amp_table", None):
            candidate = self.imagedir / "amp_phase"
            if Path(str(candidate) + ".model.tt0").is_dir():
                return candidate

        history = getattr(self, "phase_history", None) or []
        accepted = sorted(
            {
                int(item.get("loop", index))
                for index, item in enumerate(history, start=1)
                if isinstance(item, dict) and
                item.get("accepted", True) is not False
            },
            reverse=True,
        )
        for loop in accepted:
            candidate = self.imagedir / f"phase_{loop}"
            if Path(str(candidate) + ".model.tt0").is_dir():
                return candidate

        initial = self.imagedir / "initial"
        if Path(str(initial) + ".model.tt0").is_dir():
            return initial
        return None

    # --------------------------------------------------------

    def flag_residual_rfi(self):
        """Flag residual interference, without ever failing the stage.

        This is a quality step: the images are produced with or without it,
        so anything that goes wrong here is recorded and reported rather than
        raised. Its own guard against over-flagging stays inside
        ``attempt_residual_flagging``.
        """

        try:
            self.attempt_residual_flagging()
        except Exception as exc:
            record = {
                "status": "failed",
                "reason": f"{type(exc).__name__}: {exc}",
                "measurement_set": getattr(self, "current_ms", None),
            }
            self.residual_flagging = record
            try:
                self.save_json(
                    record,
                    self.reportdir / "residual_flagging.json",
                )
            except Exception:
                pass
            print(
                "WARNING: residual-RFI flagging could not run "
                f"({record['reason']}). The final image is made from the "
                "data as it stands."
            )

    def attempt_residual_flagging(self):
        """Flag interference left in the residuals before the final image.

        Interference that survives calibration is what leaves stripes and
        ripples in the image. Detection runs on the residual, so the field's
        own model is subtracted first and real emission is not what the
        detector sees. The pass is bounded: the flag state is saved first and
        restored in full when more than the configured fraction of the data is
        removed, which is what a wrong model or a mis-scaled solution causes.
        """

        checkpoint = "residual_flagging"

        if not bool(self.task_config.get("residual_flagging", True)):
            print("\nSkipping residual-RFI flagging (disabled).")
            return

        if self.checkpoint_exists(checkpoint):
            print("Residual-RFI flagging already done.")
            if getattr(self, "residual_flagging", None) is None:
                self.residual_flagging = self.saved_record(
                    "residual_flagging.json"
                )
            return

        image_base = self.latest_image_base()
        if image_base is None:
            print(
                "\nSkipping residual-RFI flagging: no self-calibration "
                "image is available to subtract."
            )
            self.record_skipped_residual_flagging(
                "no self-calibration image is available to subtract"
            )
            return

        # Detection on the residual only means anything once the field has
        # been deconvolved. With no clean components the "residual" is the
        # data itself, and the detector would see real emission.
        diagnostics = self.deconvolution_diagnostics(image_base)
        if not diagnostics.get("model_nonzero"):
            print(
                "\nWARNING: skipping residual-RFI flagging because "
                f"{Path(image_base).name} has no clean components. Nothing "
                "was deconvolved, so there is no sky model to subtract."
            )
            self.record_skipped_residual_flagging(
                f"{Path(image_base).name} has an empty clean model",
                diagnostics=diagnostics,
            )
            return

        print("\nFlagging residual RFI on the model-subtracted data.")
        self.ensure_model_column(self.current_ms, image_base)
        datacolumn = self.residual_datacolumn(self.current_ms)
        if datacolumn is None:
            print(
                "WARNING: no usable MODEL_DATA column is present, so "
                "residual-RFI flagging was skipped."
            )
            return

        before = flagdata(
            vis=self.current_ms,
            mode="summary",
            action="calculate",
        )
        before_fraction = self.summary_flagged_fraction(before)

        version = f"before_residual_flagging_{os.getpid()}_{time.time_ns()}"
        flagmanager(
            vis=self.current_ms,
            mode="save",
            versionname=version,
            comment="Before self-calibration residual-RFI flagging",
        )

        parameters = {
            "datacolumn": datacolumn,
            "timedevscale": float(
                self.task_config.get("residual_flag_timedevscale", 5.0)
            ),
            "freqdevscale": float(
                self.task_config.get("residual_flag_freqdevscale", 5.0)
            ),
            "winsize": int(self.task_config.get("residual_flag_winsize", 5)),
            "extend": bool(self.task_config.get("residual_flag_extend", True)),
            "growtime": float(
                self.task_config.get("residual_flag_growtime", 80.0)
            ),
            "growfreq": float(
                self.task_config.get("residual_flag_growfreq", 80.0)
            ),
            "image": str(image_base),
            "flag_version": version,
        }

        flagdata(
            vis=self.current_ms,
            mode="rflag",
            datacolumn=datacolumn,
            timedevscale=parameters["timedevscale"],
            freqdevscale=parameters["freqdevscale"],
            winsize=parameters["winsize"],
            # Each scan is detected on its own, and the detector never
            # extends flags itself; extension is a separate bounded pass.
            combinescans=False,
            extendflags=False,
            action="apply",
            flagbackup=False,
        )

        if parameters["extend"]:
            flagdata(
                vis=self.current_ms,
                mode="extend",
                growtime=parameters["growtime"],
                growfreq=parameters["growfreq"],
                growaround=False,
                flagneartime=False,
                flagnearfreq=False,
                combinescans=False,
                extendpols=True,
                action="apply",
                flagbackup=False,
            )

        after = flagdata(
            vis=self.current_ms,
            mode="summary",
            action="calculate",
        )
        after_fraction = self.summary_flagged_fraction(after)
        new_fraction = (
            None
            if before_fraction is None or after_fraction is None
            else max(0.0, after_fraction - before_fraction)
        )
        limit = float(
            self.task_config.get("residual_flag_max_new_fraction", 0.10)
        )

        record = {
            "status": "completed",
            "parameters": parameters,
            "before_fraction": before_fraction,
            "after_fraction": after_fraction,
            "new_fraction": new_fraction,
            "max_new_fraction": limit,
            "measurement_set": self.current_ms,
        }

        if new_fraction is None or new_fraction > limit:
            flagmanager(
                vis=self.current_ms,
                mode="restore",
                versionname=version,
            )
            record["status"] = "reverted"
            record["reason"] = (
                "the pass could not be measured"
                if new_fraction is None else
                f"it removed {new_fraction:.1%} of the data, more than the "
                f"{limit:.1%} limit"
            )
            print(
                "WARNING: residual-RFI flagging was undone because "
                f"{record['reason']}. The final image is made from the data "
                "as it was before this pass."
            )
        else:
            print(
                f"Residual-RFI flagging removed {new_fraction:.2%} of the "
                f"data ({after_fraction:.2%} flagged in total)."
            )

        self.residual_flagging = record
        self.save_json(record, self.reportdir / "residual_flagging.json")
        self.mark_checkpoint(checkpoint)

    # --------------------------------------------------------

    def final_imaging(self):

        checkpoint = "final_image"

        if self.checkpoint_exists(checkpoint):

            self.final_image = str(
                self.imagedir /
                "final.image.tt0"
            )

            self.clean_diagnostics[
                "final"
            ] = self.deconvolution_diagnostics(
                self.imagedir / "final"
            )

            return

        imagename = str(

            self.imagedir /
            "final"

        )

        self.archive_incomplete_images(imagename, "final selfcal image")

        print(
            "\nGenerating final image."
        )

        kwargs = {
            "vis": self.current_ms,
            "imagename": imagename,
            "specmode": "mfs",
            "deconvolver": "mtmfs",
            "nterms": self.image_config["nterms"],
            "weighting": "briggs",
            "robust": self.image_config["robust"],
            "imsize": self.image_config["imsize"],
            "cell": self.image_config["cell"],
            "niter": self.image_config["niter"],
            "nsigma": self.image_config["nsigma"],
            "savemodel": "modelcolumn",
            "calcpsf": True,
            "calcres": True,
            "interactive": False,
        }
        kwargs.update(self.gridder_kwargs())
        kwargs.update(
            self.automask_kwargs(
                getattr(
                    self,
                    "initial_image",
                    None,
                )
            )
        )
        kwargs.update(self.clean_control_kwargs())
        self.run_tclean_with_recovery(
            kwargs,
            "final",
        )

        self.final_image = (
            imagename +
            ".image.tt0"
        )

        self.mark_checkpoint(
            checkpoint
        )

    # --------------------------------------------------------
    # Generate QA plots
    # --------------------------------------------------------

    def generate_native_qa_plots(self):
        """Render self-cal QA figures without relying on CASA plotms."""

        status = {"plots": [], "warnings": []}
        plot_specs = []

        if hasattr(self, "initial_image"):
            plot_specs.append(
                (
                    "initial_image.png",
                    "Initial self-calibration image",
                    lambda outfile: QA_PRODUCTS.plot_casa_image(
                        self.initial_image, outfile, "Initial self-calibration image"
                    ),
                )
            )

        if hasattr(self, "final_image"):
            plot_specs.append(
                (
                    "final_image.png",
                    "Final self-calibration image",
                    lambda outfile: QA_PRODUCTS.plot_casa_image(
                        self.final_image, outfile, "Final self-calibration image"
                    ),
                )
            )

        for index, caltable in enumerate(getattr(self, "phase_tables", []), 1):
            plot_specs.append(
                (
                    f"phase{index}_phase_vs_time.png",
                    f"Self-cal phase solutions: iteration {index}",
                    lambda outfile, table=caltable, number=index: (
                        QA_PRODUCTS.plot_caltable_time_series(
                            table,
                            outfile,
                            "phase",
                            f"Self-cal phase solutions: iteration {number}",
                        )
                    ),
                )
            )

        amp_table = getattr(self, "amp_table", None)
        if amp_table:
            plot_specs.extend(
                [
                    (
                        "amp_gain_vs_time.png",
                        "Self-cal amplitude solutions",
                        lambda outfile: QA_PRODUCTS.plot_caltable_time_series(
                            amp_table,
                            outfile,
                            "amplitude",
                            "Self-cal amplitude solutions",
                        ),
                    ),
                    (
                        "amp_phase_vs_time.png",
                        "Self-cal amplitude-solve phases",
                        lambda outfile: QA_PRODUCTS.plot_caltable_time_series(
                            amp_table,
                            outfile,
                            "phase",
                            "Self-cal amplitude-solve phases",
                        ),
                    ),
                ]
            )

        for filename, label, plotter in plot_specs:
            try:
                outfile = plotter(self.plotdir / filename)
            except Exception as exc:
                status["warnings"].append(
                    f"Could not create {filename} ({label}): {exc}"
                )
                continue

            if outfile is None:
                status["warnings"].append(
                    f"No usable input is available for {filename} ({label})."
                )
            else:
                status["plots"].append(str(outfile))

        if not plot_specs:
            status["warnings"].append(
                "No self-calibration image or gain table is available yet."
            )

        status_path = self.reportdir / "qa_plot_status.json"
        QA_PRODUCTS.write_status(status_path, status)

        if status["plots"]:
            print(
                "Self-calibration QA plots written: " +
                ", ".join(status["plots"]),
                flush=True,
            )
        else:
            print(
                "WARNING: No self-calibration QA plot was created. See "
                f"{status_path}.",
                flush=True,
            )

        return status

    def measure_final_image(self):
        """Record the final self-calibration image's peak, noise and range.

        These used to be measured after the plotms figures. When the figures
        moved to the CASA-native path, ``generate_qa_plots`` began returning
        before it reached them, so every report said ``"final_metrics": {}``
        and ``selfcal_metrics.json`` was never written.
        """

        if not getattr(self, "final_image", None):
            return None

        try:
            stats = CASA_COMPAT.peak_and_noise(imstat, self.final_image)
        except Exception as exc:
            print(
                "WARNING: could not measure the final self-calibration "
                f"image {self.final_image}: {exc}"
            )
            return None

        metrics = {
            "peak_Jy_beam": stats["max"],
            "rms_Jy_beam": stats["rms"],
            "dynamic_range": stats["dynamic_range"],
        }
        self.save_json(metrics, self.reportdir / "selfcal_metrics.json")
        self.final_metrics = metrics
        return metrics

    def generate_qa_plots(self):

        print("\nGenerating QA products...")

        self.measure_final_image()

        # This CASA-native path works in the ilifu CASA 6.7 build, where
        # casatasks does not provide plotms.
        return self.generate_native_qa_plots()

    # --------------------------------------------------------
    # Determine Selfcal Strategy
    # --------------------------------------------------------

    @staticmethod
    def solint_seconds(solint, integration_s=None, scan_s=None):
        """A gaincal solint in seconds, or None if it cannot be known."""

        text = str(solint).strip().lower()
        if text == "int":
            return integration_s
        if text in {"inf", "scan"}:
            return scan_s
        for unit, factor in (("min", 60.0), ("ms", 1.0e-3), ("s", 1.0),
                             ("h", 3600.0)):
            if text.endswith(unit):
                try:
                    return float(text[:-len(unit)]) * factor
                except ValueError:
                    return None
        try:
            return float(text)
        except ValueError:
            return None

    def phase_solution_snr(self, solint):
        """Expected S/N of one antenna's phase solution at ``solint``.

        One antenna's gain is solved from every baseline to it, every channel
        of its spectral window and, for gaintype T, both parallel hands, over
        one solution interval. Each of those samples carries the model's
        visibility amplitude against the per-visibility noise, so

            S/N = A_model * sqrt(n) / sigma

        with A_model the model's rms visibility amplitude, sigma the
        per-component visibility noise, and n that sample count. This is the
        quantity gaincal's minsnr is applied to, measured before the solve.
        Returns None when any ingredient is missing.
        """

        coherence = getattr(self, "coherence", None) or {}
        initial = (getattr(self, "clean_diagnostics", None) or {}).get(
            "initial"
        ) or {}
        amplitude = initial.get("model_rms_visibility_Jy")
        sigma = coherence.get("component_sigma")
        total = coherence.get("estimated_unflagged_samples")
        antennas = coherence.get("antennas_with_data")
        on_source = coherence.get("on_source_s")
        if not all((amplitude, sigma, total, antennas, on_source)):
            return None
        if antennas < 2:
            return None

        scans = coherence.get("scans") or 1
        scan_s = on_source / float(scans)
        seconds = self.solint_seconds(
            solint, coherence.get("integration_s"), scan_s
        )
        if not seconds:
            return None
        # combine='' keeps each solution inside one scan.
        seconds = min(seconds, scan_s)

        correlations = coherence.get("correlations") or 2
        parallel = min(2, correlations)
        gaintype = str(self.task_config.get("phase_gaintype", "G")).upper()
        hands = parallel if gaintype == "T" else 1

        per_solution = (
            total
            * (2.0 / antennas)
            / float(coherence.get("spectral_windows") or 1)
            * (hands / float(correlations))
            * (seconds / on_source)
        )
        if per_solution <= 0:
            return None
        return amplitude * math.sqrt(per_solution) / sigma

    def phase_selfcal_is_warranted(self, dynamic_range, phase_solints):
        """Whether phase self-calibration should start.

        The image dynamic range is a sufficient reason, but it is the wrong
        statistic to require. A field whose dynamic range is held down by
        antenna gain errors, the spokes around every bright source, is
        exactly the field self-calibration exists to fix, and those same
        errors keep its dynamic range below any fixed threshold. What
        decides whether phase solutions can be trusted is their own S/N, so
        that is measured too, for the longest solution interval, and either
        reason is enough. Each loop is still accepted only if the image
        improves, and gaincal still rejects solutions below minsnr.
        """

        start_dr = float(
            self.task_config.get(
                "phase_start_dynamic_range",
                max(20.0, 3.0 * self.band_config["min_snr_phase"]),
            )
        )
        minimum_snr = float(
            self.task_config.get("phase_start_min_solution_snr", 10.0)
        )
        solint = phase_solints[0] if phase_solints else None
        snr = self.phase_solution_snr(solint) if solint else None
        self.phase_start_evidence = {
            "dynamic_range": dynamic_range,
            "start_dynamic_range": start_dr,
            "first_solint": solint,
            "expected_solution_snr": snr,
            "min_solution_snr": minimum_snr,
        }

        if dynamic_range >= start_dr:
            self.phase_start_evidence["reason"] = "dynamic_range"
            return True
        if snr is not None and minimum_snr > 0 and snr >= minimum_snr:
            self.phase_start_evidence["reason"] = "solution_snr"
            print(
                f"Image dynamic range {dynamic_range:.1f} is below "
                f"{start_dr:.1f}, but each antenna's phase solution at "
                f"solint={solint} is expected to reach S/N {snr:.0f} "
                f"(needs {minimum_snr:.0f}). A field held back by gain "
                "errors is what self-calibration is for, so it starts."
            )
            return True
        self.phase_start_evidence["reason"] = "neither"
        return False

    def check_weights_against_sky(self):
        """Whether the sky could have biased the statistical weights.

        statwt takes the scatter of the visibilities for noise. The sky adds
        to that scatter, most on long baselines where off-centre sources
        fringe fastest, so a field bright enough to matter would have its
        long baselines down-weighted. It matters only when the model's
        rms visibility amplitude is a sizeable fraction of the noise.
        """

        if not self.statistical_weights_applied():
            return None
        amplitude = (
            (getattr(self, "clean_diagnostics", None) or {}).get("initial")
            or {}
        ).get("model_rms_visibility_Jy")
        sigma = (getattr(self, "coherence", None) or {}).get(
            "component_sigma"
        )
        if not amplitude or not sigma:
            return None
        ratio = amplitude / sigma
        limit = float(self.task_config.get("statwt_max_sky_to_noise", 0.1))
        self.weights_sky_to_noise = ratio
        if ratio > limit:
            print(
                "WARNING: the sky model's rms visibility amplitude is "
                f"{ratio:.2f} of the per-visibility noise. statwt counted "
                "that signal as noise, so long baselines are under-weighted "
                "and the image resolution and noise suffer. Set "
                "ugmrt_statwt_enabled: false for this field, or reweight "
                "from residuals."
            )
        return ratio

    def determine_solints(self):

        self.check_weights_against_sky()

        dynamic_range = self.image_stats["dynamic_range"]

        phase_solints = []

        for solint in self.band_config[
            "phase_solints"
        ]:

            phase_solints.append(
                solint
            )

        amp_solint = None

        # An image dynamic range and a per-solution gaincal S/N are different
        # quantities.  min_snr_amp is the gaincal minsnr (symmetric with
        # min_snr_phase); the decision to attempt amplitude self-calibration
        # at all is governed by its own image-quality threshold.  Amplitude
        # self-calibration absorbs sky-model error into the gains, so it needs
        # a genuinely well-modelled field, not merely a detection.
        amp_start_dynamic_range = float(
            self.task_config.get("amp_start_dynamic_range", 100.0)
        )

        if dynamic_range >= amp_start_dynamic_range:

            amp_solint = (
                self.band_config[
                    "amp_solint"
                ]
            )

        self.selfcal_strategy = {

            "phase_solints":
                phase_solints,

            "amp_solint":
                amp_solint,

            "phase_selfcal":
                self.phase_selfcal_is_warranted(
                    dynamic_range, phase_solints
                ),

            # Whether it is attempted is decided after the phase loops, on
            # the image they produce; see run_amp_phase_selfcal.
            "amp_selfcal":
                bool(
                    self.task_config.get(
                        "enable_amp_selfcal",
                        False,
                    )
                )
        }

        print(
            "\nSelected strategy:"
        )

        print(
            f"Phase selfcal : "
            f"{self.selfcal_strategy['phase_selfcal']}"
        )

        print(
            "Phase start DR : "
            f"{float(self.task_config.get('phase_start_dynamic_range', 20.0)):.1f}"
        )

        print(
            f"Phase solints : "
            f"{phase_solints}"
        )

        print(
            f"Amp selfcal   : "
            f"{self.selfcal_strategy['amp_selfcal']}"
        )

        print(
            "Amp start DR  : "
            f"{amp_start_dynamic_range:.1f}"
        )

        print(
            f"Amp solint    : "
            f"{amp_solint}"
        )    
    
    # --------------------------------------------------------
    # Selfcal Report
    # --------------------------------------------------------

    def write_report(self):

        beam = None

        if hasattr(self, "final_image"):

            try:

                header = imhead(
                    imagename=self.final_image,
                    mode="summary",
                )

                if isinstance(header, dict):

                    beam = header.get("restoringbeam")

            except Exception as exc:

                print(
                    "WARNING: could not read the final-image restoring "
                    f"beam for the report: {exc}"
                )

        final_ms = getattr(
            self,
            "current_ms",
            self.ms,
        )
        input_ms_identity = IMAGING_SAFETY.path_stat_fingerprint(
            self.input_ms
        )
        final_ms_identity = IMAGING_SAFETY.path_stat_fingerprint(final_ms)

        report = {

            "provenance_fingerprint": getattr(
                self,
                "provenance_fingerprint",
                None,
            ),

            "input_measurement_set": self.input_ms,

            "input_measurement_set_stat_sha256": input_ms_identity,

            "measurement_set":
                self.ms,

            "band":
                self.band,

            "telescope":
                self.telescope,

            "target":
                self.target,

            "execution_profile": self.profile,

            "residual_flagging": getattr(self, "residual_flagging", None),

            # True when CLEAN found nothing: the image is the dirty map.
            "sky_model_empty": self.sky_model_is_empty(),

            # Vector- versus scalar-averaged amplitudes of the calibrated
            # data, which separate a faint field from incoherent data.
            "calibrated_coherence": getattr(self, "coherence", None),
            "antenna_data": getattr(self, "antenna_data", None),
            "phase_start_evidence": getattr(
                self, "phase_start_evidence", None
            ),
            "weights_sky_to_noise": getattr(
                self, "weights_sky_to_noise", None
            ),

            # Uncheckpointed products left by an interrupted run and moved
            # aside, not reused, before their step was redone.
            "archived_incomplete_products": list(
                getattr(self, "archived_incomplete_products", [])
            ),

            "qa_plot_status": str(
                self.reportdir / "qa_plot_status.json"
            ),

            "reference_antenna":
                self.reference_antenna,

            "phase_solints":
                self.band_config[
                    "phase_solints"
                ],

            "amp_solint":
                self.band_config[
                    "amp_solint"
                ],

            "min_snr_phase":
                self.band_config[
                    "min_snr_phase"
                ],

            "min_snr_amp":
                self.band_config[
                    "min_snr_amp"
                ],

            "image_statistics":
                getattr(
                    self,
                    "image_stats",
                    {}
                ),

            "selfcal_strategy":
                getattr(
                    self,
                    "selfcal_strategy",
                    {}
                ), 

            "final_ms":
                final_ms,

            "final_ms_stat_sha256": final_ms_identity,

            "amp_table":
                getattr(
                    self,
                    "amp_table",
                    None
                ),

            # Antennas whose remnant data were flagged because no
            # self-calibration solution could be made for them.
            "unsolvable_antenna_flagging":
                getattr(
                    self,
                    "unsolvable_antenna_flagging",
                    None
                ) or self.saved_record("unsolvable_antenna_flagging.json"),

            # The amp+phase candidate against the phase-only image it was
            # compared with, and why it was kept or discarded.
            "amp_history":
                getattr(
                    self,
                    "amp_history",
                    None
                ),

            "phase_tables":
                getattr(
                    self,
                    "phase_tables",
                    []
                ), 

            "phase_history":
                getattr(
                    self,
                    "phase_history",
                    []
                ),

            "calibration_table_qa":
                getattr(
                    self,
                    "caltable_qa",
                    {}
                ),

            "final_image":
                getattr(
                    self,
                    "final_image",
                    None
                ),

            "final_metrics":
                getattr(
                    self,
                    "final_metrics",
                    {}
                ), 

            "imaging_parameters": {
                "cell": self.image_config["cell"],
                "imsize": self.image_config["imsize"],
                "weighting": "briggs",
                "robust": self.image_config["robust"],
                "nterms": self.image_config["nterms"],
                "niter": self.image_config["niter"],
                "nsigma": self.image_config["nsigma"],
                "automask": self.automask_kwargs(),
                "clean_controls": self.clean_control_kwargs(),
                "flat_imaging_overrides": getattr(
                    self,
                    "flat_imaging_overrides",
                    [],
                ),
            },

            "clean_diagnostics": self.clean_diagnostics,
            "clean_passes": self.recorded_clean_passes(),
            "automask_adaptations": getattr(
                self, "automask_adaptations", {}
            ),

            "beam": beam,

            "checkpoints": {

                "initial_image":
                    self.checkpoint_exists(
                        "initial_image"
                    ),

                "phase_selfcal":
                    self.checkpoint_exists(
                        "phase_selfcal"
                    ),

                "amp_selfcal":
                    self.checkpoint_exists(
                        "amp_selfcal"
                    ),

                "final_image":
                    self.checkpoint_exists(
                        "final_image"
                    ),

                "selfcal":
                    self.checkpoint_exists(
                        "selfcal"
                    )
            }
        }

        outfile = (
            self.reportdir /
            "selfcal_report.json"
        )

        self.save_json(
            report,
            outfile
        )

        print(
            f"Selfcal report written:\n"
            f"{outfile}"
        )

    # --------------------------------------------------------
    # Main driver
    # --------------------------------------------------------

    def run(self):

        print("\n")
        print("=" * 70)
        print(
            f"Starting {TELESCOPE_PROFILE_MODULE.get_telescope_profile(self.telescope)['label']} "
            "Self-Calibration"
        )
        print("=" * 70)

        self.print_summary()

        self.check_calibrated_coherence()

        self.initial_imaging()

        self.estimate_snr()

        self.determine_solints()

        self.run_phase_selfcal()

        self.run_amp_phase_selfcal()

        self.flag_residual_rfi()

        self.final_imaging()

        self.warn_if_no_sky_model()

        self.generate_qa_plots()

        self.mark_checkpoint(
            "selfcal"
        )

        self.write_report()

        print("\n")
        print("=" * 70)
        print(
            f"{TELESCOPE_PROFILE_MODULE.get_telescope_profile(self.telescope)['label']} "
            "self-calibration completed successfully."
        )
        print("=" * 70)

# ============================================================
# CASA ENTRY POINT
# ============================================================


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

    parser = argparse.ArgumentParser(
        description="Run profile-aware CASA target self-calibration."
    )

    parser.add_argument(
        "ms",
        help="Calibrated target Measurement Set.",
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
        "--upstream-outdir",
        default=None,
        help=(
            "Directory holding the shared inspection and calibration "
            "reports. Defaults to --outdir; the orchestrator points it at "
            "the run's root when each target has its own product directory."
        ),
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
            "Execution profile. 'laptop' reduces image size and self-cal "
            "iterations for a 36-GB Mac."
        ),
    )

    parser.add_argument(
        "--telescope",
        choices=TELESCOPE_PROFILE_MODULE.TELESCOPE_CHOICES,
        default=RUN_CONFIG.get("telescope", "auto"),
        help="Telescope identity; verified against available stage reports.",
    )

    parser.add_argument(
        "--receiver-band",
        default=RUN_CONFIG.get("receiver_band"),
        help="Explicit telescope-scoped receiver band.",
    )

    parser.add_argument(
        "--reference-antenna",
        default=RUN_CONFIG.get("reference_antenna"),
        help="Reference antenna for a direct calibrated-MS selfcal run.",
    )

    parser.add_argument(
        "--target",
        default=RUN_CONFIG.get("target"),
        help="Target field label recorded in self-calibration provenance.",
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
            "Regenerate self-calibration QA plots from existing images and "
            "gain tables without modifying the Measurement Set."
        ),
    )

    args = parser.parse_args()

    manager = SelfCalManager(
        args.ms,
        outdir=args.outdir,
        upstream_outdir=args.upstream_outdir,
        profile=args.profile,
        telescope=args.telescope,
        receiver_band=args.receiver_band,
        reference_antenna=args.reference_antenna,
        target=args.target,
        provenance_fingerprint=args.provenance_fingerprint,
    )

    try:
        if args.plots_only:
            manager.generate_qa_plots()
        else:
            manager.run()
    except KeyboardInterrupt:
        _interrupted("Self-calibration")
