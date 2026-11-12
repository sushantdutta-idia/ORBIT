#!/usr/bin/env python3
# ============================================================
# config.py
#
# Central configuration for the multi-telescope radio continuum pipeline.
#
# This module is intentionally CASA-free so it can be imported by normal Python tooling, wrappers, tests, and pipeline launchers.
#
# Pipeline order:
#   fits_ingest.py (visibility FITS or GMRT LTA only) -> inspect.py -> flagging.py
#   -> calibration.py -> selfcal.py
#   -> imaging.py -> source_detection.py -> ddcal.py
#
# Author: Sushant Dutta
# ============================================================

import copy
import importlib.util
import math
import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of failing inside a stage.
PIPELINE_INTERFACE_VERSION = 5


def _load_telescope_profiles_module():

    module_path = Path(__file__).with_name("telescope_profiles.py")
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_telescope_profiles",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TELESCOPE_PROFILE_MODULE = _load_telescope_profiles_module()
TELESCOPE_PROFILES = TELESCOPE_PROFILE_MODULE.TELESCOPE_PROFILES
TELESCOPE_CHOICES = TELESCOPE_PROFILE_MODULE.TELESCOPE_CHOICES
DATA_STATE_CHOICES = TELESCOPE_PROFILE_MODULE.DATA_STATE_CHOICES
IMAGER_CHOICES = TELESCOPE_PROFILE_MODULE.IMAGER_CHOICES
PROFILE_SCHEMA_VERSION = TELESCOPE_PROFILE_MODULE.PROFILE_SCHEMA_VERSION


def _parse_config_yaml_scalar(raw_value):

    value = raw_value.strip()

    if value == "":
        return ""

    if (
        len(value) >= 2 and
        value[0] in {"'", '"'} and
        value[-1] == value[0]
    ):
        return value[1:-1]

    lower = value.lower()

    if lower in {"true", "yes", "on"}:
        return True

    if lower in {"false", "no", "off"}:
        return False

    if lower in {"none", "null", "~"}:
        return None

    if value.startswith("[") and value.endswith("]"):
        items = value[1:-1].strip()

        if not items:
            return []

        return [
            _parse_config_yaml_scalar(item)
            for item in items.split(",")
        ]

    try:
        return int(value)
    except ValueError:
        pass

    try:
        return float(value)
    except ValueError:
        return value


def parse_calibrator_flux_models(value):
    """Calibrators given their flux density outright, for CASA setjy.

    ``FIELD=FLUX_JY@REFFREQ_MHZ[,SPIX[,CURVATURE]]``, several separated by
    ``;``, or a mapping of field to ``{flux_jy, reffreq_mhz, spix}``. A
    flat string survives the fallback parser used when PyYAML is absent.
    Returns ``{field: {"flux_jy", "reffreq_mhz", "spix"}}``.
    """

    def number(text, what):
        try:
            result = float(text)
        except (TypeError, ValueError):
            raise ValueError(f"{what} {text!r} is not a number") from None
        if not math.isfinite(result):
            raise ValueError(f"{what} {text!r} is not finite")
        return result

    def entry(field, flux, reffreq, spix):
        field = str(field).strip()
        if not field:
            raise ValueError("a calibrator flux model names no field")
        flux = number(flux, f"the flux density of {field}")
        reffreq = number(reffreq, f"the reference frequency of {field}")
        if flux <= 0 or reffreq <= 0:
            raise ValueError(
                f"{field}: flux density and reference frequency must be "
                "positive"
            )
        if spix is None or spix == "":
            spix = []
        elif not isinstance(spix, (list, tuple)):
            spix = [spix]
        return field, {
            "flux_jy": flux,
            "reffreq_mhz": reffreq,
            "spix": [number(term, f"the spectral index of {field}")
                     for term in spix],
        }

    models = {}
    if not value:
        return models
    if isinstance(value, dict):
        for field, spec in value.items():
            spec = spec or {}
            name, model = entry(field, spec.get("flux_jy"),
                                spec.get("reffreq_mhz"), spec.get("spix"))
            models[name] = model
        return models
    for part in str(value).split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" not in part or "@" not in part:
            raise ValueError(
                f"calibrator flux model {part!r} is not "
                "FIELD=FLUX_JY@REFFREQ_MHZ[,SPIX...]"
            )
        field, rest = part.split("=", 1)
        flux, rest = rest.split("@", 1)
        terms = [term.strip() for term in rest.split(",")]
        name, model = entry(field, flux, terms[0], terms[1:])
        models[name] = model
    return models


def parse_manual_flag_commands(value):
    """config.yaml's ``flagging_manual_commands`` as a list of commands.

    Each is one flagdata command in CASA's list syntax, checked by the
    flagging stage. Write them as a one-line ``["...", "..."]`` list without
    commas inside a command: the fallback parser used when PyYAML is absent
    splits that list at commas.
    """

    if value is None or value == "":
        return []
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, (list, tuple)):
        raise ValueError(
            "flagging_manual_commands must be a list of flagdata commands, "
            f"not {type(value).__name__}"
        )
    return [str(command).strip() for command in value if str(command).strip()]


def _load_config_yaml(path=None):

    path = (
        Path(__file__).with_name("config.yaml")
        if path is None
        else Path(path)
    )

    if not path.exists():
        return {}

    try:
        import yaml
    except Exception:
        yaml = None

    if yaml is not None:
        with open(path) as f:
            loaded = yaml.safe_load(f) or {}

        if not isinstance(loaded, dict):
            raise ValueError(
                f"Expected a mapping in {path}"
            )

        return loaded

    loaded = {}

    for line in path.read_text().splitlines():
        stripped = line.strip()

        if not stripped or stripped.startswith("#"):
            continue

        if ":" not in stripped:
            continue

        key, value = stripped.split(":", 1)
        key = key.strip()

        if not key:
            continue

        loaded[key] = _parse_config_yaml_scalar(value)

    return loaded


CONFIG_YAML = _load_config_yaml()


# ============================================================
# GLOBAL PIPELINE SETTINGS
# ============================================================

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
DEFAULT_RAM_GB = 256
# PB corrections become increasingly noisy and beam-model dependent toward
# the field edge.  A 0.30 cutoff keeps the default science image inside the
# region where a factor >3.3 noise amplification is avoided.  Users needing
# a wider exploratory image can explicitly lower this value.
DEFAULT_PBMIN = 0.3
DEFAULT_REFANTS = [
    "C00",
    "C05",
    "C08",
    "C09",
    "C01",
    "C02",
]


# ============================================================
# EXECUTION PROFILES
# ============================================================

# The default ``science`` profile retains the full-resolution settings below.
# ``laptop`` is deliberately opt-in: it is a resource-aware quick-look and
# development profile for a 36-GB Mac, not a replacement for the full
# field-of-view science imaging configuration.
EXECUTION_PROFILES = {
    "science": {
        "description": (
            "Full-resolution science imaging; suited to a cluster or a "
            "large-memory workstation."
        ),
    },
    "laptop": {
        "description": (
            "Fast, memory-bounded local run for a 36-GB Mac; validate the "
            "result with the science profile before final analysis."
        ),
        # A 36-GB Apple-silicon Mac commonly has six performance cores. A
        # small cap avoids memory-pressure and BLAS/OpenMP oversubscription.
        "threads": 6,
        "flagging": {
            "max_iterations": 1,
        },
        "selfcal": {
            "imsize": 2048,
            "niter": 25000,
            "phase_loops": 2,
            # Match the reduced final-imaging w-plane depth below so a laptop
            # self-cal image keeps the same wide-field treatment it will get
            # at final imaging.
            "wprojplanes": 64,
        },
        "imaging": {
            "imsize": 2048,
            "niter": 50000,
            "scales": [0, 5, 15, 45],
            "wprojplanes": 64,
            "subbands": 0,
            "subband_niter": 10000,
            "make_subbands": False,
            "parallel": False,
        },
    },
}

EXECUTION_PROFILE_NAMES = tuple(EXECUTION_PROFILES)


def get_execution_profile(name="science"):
    """Return an isolated execution-profile mapping.

    Stage scripts may adjust their local copy, so returning a deep copy keeps
    repeated imports and different pipeline runs independent.
    """

    normalized = str(name or "science").strip().lower()

    if normalized not in EXECUTION_PROFILES:
        choices = ", ".join(EXECUTION_PROFILE_NAMES)
        raise ValueError(
            f"Unknown execution profile '{name}'. Choose one of: {choices}."
        )

    return copy.deepcopy(EXECUTION_PROFILES[normalized])


# ============================================================
# USER RUN CONFIGURATION
# ============================================================

# Edit this block before running pipeline.py. The command line can still override these values, but a normal science run should only require setting ms, target, and any imaging/DD preferences below. Then run:
#   python3 pipeline.py

PIPELINE_RUN_CONFIG = {
    # Raw input Measurement Set, visibility FITS file or GMRT LTA file.
    # FITS and LTA are imported into a managed MS before inspection.
    # "ms": "/data/observation.uvfits",
    # "ms": "/data/observation_gwb.lta",
    "ms": None,

    # Visibility-FITS dialect. ``auto`` distinguishes random-groups UVFITS
    # from FITS-IDI UV_DATA. Ordinary image FITS is rejected.
    "fits_format": "auto",
    # Optional destination for a FITS import. By default the pipeline uses
    # <outdir>/ingest/<input-stem>.ms.
    "imported_ms": None,
    # Optional importgmrt online flag files for uGMRT UVFITS.  Paths are
    # content-bound in ingest provenance and may also be repeated on the CLI.
    "gmrt_flagfiles": [],
    # A GMRT LTA input ("ms": "/data/observation.lta") is converted to UVFITS
    # with NCRA's listscan and gvfits before importgmrt. Each is a path, a
    # name on PATH, or a whole command. Linux builds of both are in
    # https://github.com/ruta-k/CAPTURE-CASA6.
    "gmrt_listscan": "listscan",
    "gmrt_gvfits": "gvfits",
    # The converted FITS is removed once the MS is validated unless kept.
    "keep_lta_fits": False,

    # Telescope identity is read from OBSERVATION/TELESCOPE_NAME.  Keep an
    # explicit override available for older or repaired Measurement Sets.
    "telescope": "auto",
    # Required when a MeerKAT sub-band lies wholly inside the UHF/L overlap.
    "receiver_band": None,
    # ``calibrated`` means a target MS already cross-calibrated by a suitable
    # upstream workflow (notably LINC/DP3 for LOFAR).
    "data_state": "raw",
    # Auto selects CASA for uGMRT/MeerKAT and WSClean+EveryBeam for LOFAR.
    "imager": "auto",
    "wsclean": "wsclean",
    "wsclean_data_column": "CORRECTED_DATA",
    # Required explicitly for LOFAR: ``full`` trusts the MS beam metadata,
    # ``differential`` forces phase-centre differential correction, and
    # ``none`` deliberately writes no PB-corrected product.
    "lofar_beam_mode": None,
    # Explicit override of the telescope profile's CASA setjy standard.
    "flux_standard": None,
    # Optional absolute-flux CASA image or component list predicted with ft.
    # MeerKAT UHF calibration requires a validated multi-component model.
    "flux_model": None,
    "flux_model_type": "casa-image",

    # Target field name. Leave as None to self-calibrate and image every
    # calibrated science target.
    "target": None,

    # Comma-separated subset of science targets, when only some of the
    # calibrated targets should be self-calibrated and imaged.
    "targets": None,

    # What to do when several calibrated targets are available and neither
    # target nor targets is set. "all" images each of them in its own product
    # directory; "auto" ranks post-flagging statistics and images one;
    # "explicit" stops and asks.
    "target_selection": "all",

    # Explicit field selections take precedence over inspection inference.
    # Set these when the MS lacks reliable observing intents.
    "flux_calibrator": None,
    "bandpass_calibrator": None,
    "secondary_flux_calibrator": None,
    "phase_calibrator": None,
    "target_phase_map": None,
    "target_fields": None,
    "reference_antenna": None,

    # All pipeline products, including PB-corrected images, are written below this directory.
    "outdir": DEFAULT_OUTDIR,

    # Copy the input MS before flagging or calibration. The
    # pipeline operates only on this working copy, preserving the original.
    "copy_input_ms": True,

    # Optional location for that working copy. By default it is created at
    # <outdir>/working_ms/<input-name>.ms.
    "working_ms": None,

    # None means the standard continuum sequence. Visibility FITS starts with
    # ingest; a Measurement Set starts with inspection. Both then continue
    # through flagging, calibration, selfcal, and imaging.
    # Use a list to customize, for example:
    # ["inspection", "flagging", "calibration", "selfcal", "imaging",
    #  "source_detection"]
    "stages": None,

    # Optional direction-dependent calibration stage.
    "include_ddcal": False,
    "dd_mode": "auto",
    "dd_engine": "auto",
    "solver_command_template": None,
    "apply_command_template": None,

    # Keep DD-cal prepare-only unless killMS/DDFacet is available (or a
    # site-validated CubiCal workflow is configured).
    "dd_prepare_only": True,
    # When executing killMS/DDFacet, also apply the CASA uGMRT primary-beam
    # model to the final DD restored image.
    "skip_dd_post_imaging": False,

    # Executables and resources.
    "casa": "casa",
    # Environment module that provides CASA (for example "casa/6.7.0" on
    # ilifu). The pipeline loads it itself, so a run started from a terminal
    # (python3 pipeline.py --background) or as a SLURM job (--submit) finds
    # CASA without a `module load` first. None: CASA must already be on PATH.
    "casa_module": None,
    # Command that starts a Python with casatasks and casatools (a modular
    # CASA), for example "singularity exec /path/casa-modular.sif python3".
    # Stages then run as `<command> script args` instead of `casa -c`. When
    # casa_module names a module that runs CASA from a container and has no
    # `casa` program (ilifu's casa/6.7.0), that container is used by itself.
    "casa_python": None,
    "python": None,
    "ncpu": DEFAULT_NCPU,

    # ``science`` preserves the full-resolution configuration. Use
    # ``laptop`` (or ``--profile laptop``) on a 36-GB Mac.
    "profile": "science",
    # Optional explicit cap for CASA/OpenMP/BLAS threads. ``None`` uses the
    # selected profile recommendation, then all CPUs available locally or in
    # the Slurm allocation. This changes execution only, not science settings.
    "threads": None,

    # Speed up final-imaging tclean on Linux clusters. Self-cal tclean keeps
    # serial-process/OpenMP execution because it writes MODEL_DATA. ``auto``
    # requires both Slurm and the matching mpicasa;
    # it falls back to serial/OpenMP CASA on a local machine and on macOS.
    "mpi": "auto",
    "mpi_processes": None,
    "mpicasa": "mpicasa",

    # Optional explicit stage inputs. Leave as None for automatic handoff through the pipeline reports.
    "selfcal_ms": None,
    "imaging_ms": None,
    "dd_ms": None,

    # Runtime behavior.
    # Validate and reuse completed stages after a job failure or timeout.
    "resume": True,
    "continue_on_error": False,

    # flagdata memory sizing. Whole-scan chunks are kept wherever they fit;
    # otherwise the time chunk is bounded and the choice is recorded in
    # flagging/flagdata_chunk_plan.json. None detects the limit from Slurm,
    # the cgroup and physical memory, taking the smallest. The fraction leaves
    # headroom for CASA itself; bytes per element is an estimate to calibrate
    # against `sacct --format=MaxRSS`. These live here rather than in
    # FLAGGING_CONFIG so a memory change never invalidates a flagging recipe.
    "flagdata_memory_limit_gb": None,
    "flagdata_memory_fraction": 0.5,

    # Fraction of the job's memory limit that CASA is told it may use.
    # casacore otherwise reads the node's free memory under Slurm and can
    # plan single allocations larger than the whole job allowance, which is
    # what w-projection does when it sizes its convolution-function chunks.
    # Set to 0 to leave CASA's own view alone.
    "casa_memory_fraction": 0.5,
    "flagdata_bytes_per_element": 12.0,

    # Print a pipeline heartbeat while a CASA stage is still running.
    # Set to 0 to disable it.
    "progress_interval_s": 60,
}

# ============================================================
# ILIFU / SLURM BATCH CONFIGURATION
# ============================================================

# Edit this block when running on ilifu through Slurm. Generate the sbatch file with:
#   python3 pipeline.py --write-sbatch run_radio_continuum_pipeline.sbatch
#
# Then submit with: sbatch run_radio_continuum_pipeline.sbatch
#
# Keep account, partition, qos, and module commands consistent with the current ilifu environment assigned to your project.

ILIFU_SLURM_CONFIG = {
    "job_name": "radio_continuum",

    # Leave blank if your ilifu allocation does not require it.
    "account": "",
    "partition": "",
    "qos": "",

    "nodes": 1,
    # Hybrid layout: MPI ranks accelerate final tclean, while the batch
    # process can use the full allocation for serial-process/OpenMP stages.
    "ntasks": min(32, DEFAULT_NCPU),
    "cpus_per_task": max(1, DEFAULT_NCPU // min(32, DEFAULT_NCPU)),
    "mem": f"{DEFAULT_RAM_GB}G",
    "time": "200:00:00",

    # Slurm logs are written relative to the submit directory.
    "log_dir": "slurm_logs",
    "output": "%x-%j.out",
    "error": "%x-%j.err",

    # Optional email notifications.
    "mail_user": "",
    "mail_type": "END,FAIL",

    # Run location. Empty means the directory containing pipeline.py.
    "workdir": "",

    # Commands emitted before the pipeline run. Keep this empty unless the
    # required CASA/Python environment needs an explicit module or activation
    # command. A bare "module purge" would discard a working CASA setup.
    "module_commands": [],
    "activate_commands": [],

    # Extra shell commands before/after the pipeline. Useful for
    # scratch setup, environment variables, or site-specific storage.
    "pre_commands": [],
    "post_commands": [],

    # Thread environment. pipeline.py switches MPI imaging stages to one
    # OpenMP/BLAS thread per rank to prevent oversubscription.
    "set_thread_env": True,

    # Additional command-line arguments passed to pipeline.py inside
    # the batch job. Normally leave empty and edit PIPELINE_RUN_CONFIG.
    "pipeline_args": [],

    # Default generated script name.
    "script_path": "run_radio_continuum_pipeline.sbatch",
}

# ============================================================
# BAND DEFINITIONS
# ============================================================

# Telescope-scoped frequency limits come from telescope_profiles.py.  The
# legacy UGMRT_BANDS name remains available for existing imports.
TELESCOPE_BANDS = {
    telescope: copy.deepcopy(profile["bands"])
    for telescope, profile in TELESCOPE_PROFILES.items()
}
UGMRT_BANDS = TELESCOPE_BANDS["ugmrt"]
MEERKAT_BANDS = TELESCOPE_BANDS["meerkat"]
LOFAR_BANDS = TELESCOPE_BANDS["lofar"]
ALL_BANDS = TELESCOPE_PROFILE_MODULE.all_band_profiles()
MEERKAT_S_BANDS = tuple(
    band for band in MEERKAT_BANDS
    if band.startswith("MEERKAT_S")
)


# ============================================================
# INSPECTION
# ============================================================

INSPECTION_CONFIG = {
    "flux_catalogue": "flux-cals.list",
    "phase_catalogue": "vla-cals.list",
    "preferred_refants": DEFAULT_REFANTS,
    "preferred_refants_by_telescope": {
        name: list(profile["preferred_refants"])
        for name, profile in TELESCOPE_PROFILES.items()
    },
}


# ============================================================
# FLAGGING
# ============================================================

FLAGGING_CONFIG = {
    # An MS that an earlier run of this pipeline already flagged is refused,
    # so that a reused working copy cannot have its flags stacked again.
    "allow_previously_flagged_input": False,
    # One conservative raw-DATA pass limits cumulative erosion. A second,
    # guarded TFCrop+RFlag pass occurs on CORRECTED_DATA after a provisional
    # solve, followed by a completely fresh calibration solve.
    "max_iterations": 1,
    "precal_rflag": False,
    "bad_ant_threshold": 0.7,
    # A telescope-profile reference antenna is preferred over the strictly
    # least-flagged one whenever its flagged fraction is within this margin.
    # Array geometry matters as much as cleanliness for a phase reference.
    "refant_preference_tolerance": 0.05,
    "edge_fraction": {
        "B2": 0.07,
        "B3": 0.05,
        "B4": 0.05,
        "B5": 0.03,
        "MEERKAT_UHF": 0.02,
        "MEERKAT_L": 0.01,
        **{band: 0.01 for band in MEERKAT_S_BANDS},
        # LOFAR sub-band edges and RFI are handled upstream by LINC/AOFlagger.
        # A calibrated LOFAR input is therefore not percentage-edge flagged.
        "LOFAR_LBA": 0.0,
        "LOFAR_HBA": 0.0,
    },
    # Persistent MeerKAT L-band ranges from the processMeerKAT defaults.
    # They are telescope/band scoped and are never applied to overlapping
    # uGMRT or LOFAR frequencies.
    "known_frequency_masks": {
        "MEERKAT_L": [
            "933~960MHz",
            "1163~1299MHz",
            "1524~1630MHz",
        ],
    },
    # A blind scan-start quack is not a telescope-independent science rule.
    # Enable it only after the observation log/waterfall demonstrates a
    # repeatable settling transient, then record the justified interval.
    "quack_interval_s": 0.0,
    "flag_autocorrelations": True,
    "flag_shadowed": True,
    "shadow_tolerance_m": 0.0,
    # Detection is done on parallel hands. Flags are propagated to all
    # correlations once, after the auto-flag passes have converged.
    "extend_growtime": 80.0,
    "extend_growfreq": 80.0,
    "extend_polarizations": True,
    "convergence_new_fraction": 0.001,
    "max_iteration_new_fraction": 0.20,
    "max_total_flagged_fraction": 0.70,
    # A single calibrator field or SPW can be damaged while the whole-MS
    # fraction still looks harmless. These group-level limits are checked
    # whenever CASA includes field/SPW breakdowns in its summary.
    "max_group_iteration_new_fraction": 0.35,
    "max_group_total_flagged_fraction": 0.85,
    "max_extension_new_fraction": 0.10,
    "max_group_extension_new_fraction": 0.20,
    "autoflag": {
        "default": {
            "tfcrop_timecutoff": 4.5,
            "tfcrop_freqcutoff": 4.0,
            "tfcrop_maxnpieces": 7,
            "rflag_timedevscale": 6.0,
            "rflag_freqdevscale": 6.0,
            "rflag_winsize": 5,
            # Start conservatively, then approach CASA's standard RFlag
            # 5-sigma scale only if the newly flagged fraction has not
            # converged after the preceding pass.
            "tighten_per_iteration": 0.50,
            "tfcrop_timecutoff_floor": 4.0,
            "tfcrop_freqcutoff_floor": 3.5,
            "rflag_timedevscale_floor": 5.0,
            "rflag_freqdevscale_floor": 5.0,
        },
        # Lower-frequency uGMRT bands usually contain more persistent
        # narrow-band RFI. The small change remains more conservative than
        # CASA's default RFlag scale of 5.
        "B2": {
            "rflag_freqdevscale": 5.5,
        },
        "B3": {
            "rflag_freqdevscale": 5.5,
        },
    },
    # The experimental channel-level IsolationForest complement is opt-in.
    # Deterministic CASA flags and their saved thresholds remain the default.
    "ml_rfi": {
        "enabled": False,
    },
}

# Overrides are merged after inspection has identified the instrument.  They
# intentionally affect flagging mechanics only; receiver-band thresholds stay
# in the band-indexed mappings above.
TELESCOPE_FLAGGING_CONFIG = {
    "ugmrt": {
        "precal_rflag": False,
    },
    "meerkat": {
        "quack_interval_s": 0.0,
        "flag_shadowed": True,
        # processMeerKAT defers median-based RFlag until the provisional
        # bandpass has flattened the data.
        "precal_rflag": False,
    },
    "lofar": {
        "max_iterations": 1,
        "quack_interval_s": 0.0,
        "flag_shadowed": False,
        "extend_polarizations": True,
        "autoflag": {
            "default": {
                "tfcrop_timecutoff": 6.0,
                "tfcrop_freqcutoff": 6.0,
                "rflag_timedevscale": 8.0,
                "rflag_freqdevscale": 8.0,
            },
        },
    },
}


# ============================================================
# CALIBRATION
# ============================================================

CAL_CONFIG = {
    "B2": {
        "delay_solint": "inf",
        "phase_solint": "60s",
        "amp_solint": "600s",
        "combine_delay": "scan",
        "hanningsmooth": True,
    },
    "B3": {
        "delay_solint": "inf",
        "phase_solint": "30s",
        "amp_solint": "300s",
        "combine_delay": "scan",
        "hanningsmooth": True,
    },
    "B4": {
        "delay_solint": "inf",
        "phase_solint": "20s",
        "amp_solint": "180s",
        "combine_delay": "scan",
        "hanningsmooth": True,
        # Widest channel a calibrated target is averaged to. The band's
        # 97.7 kHz channels become 390.6 kHz: for the 5120 x 1 arcsec field
        # and its 4.9 x 4.5 arcsec beam, bandwidth smearing lowers a source's
        # peak by 2% at the image edge and under 1% within 30 arcmin, and
        # each major cycle reads and grids a quarter of the data. Other
        # bands keep their native channels until their fields are reviewed.
        "target_max_channel_width_khz": 400.0,
    },
    "B5": {
        "delay_solint": "inf",
        "phase_solint": "10s",
        "amp_solint": "120s",
        "combine_delay": "",
        "hanningsmooth": False,
    },
}

CAL_CONFIG.update({
    "MEERKAT_UHF": {
        "delay_solint": "inf",
        "phase_solint": "16s",
        "amp_solint": "128s",
        "combine_delay": "",
        "hanningsmooth": False,
    },
    "MEERKAT_L": {
        "delay_solint": "inf",
        "phase_solint": "8s",
        "amp_solint": "64s",
        "combine_delay": "",
        "hanningsmooth": False,
    },
    **{
        band: {
            "delay_solint": "inf",
            "phase_solint": "8s",
            "amp_solint": "64s",
            "combine_delay": "",
            "hanningsmooth": False,
        }
        for band in MEERKAT_S_BANDS
    },
})

CALIBRATION_TASK_CONFIG = {
    "flux_standard": "Perley-Butler 2017",
    "delay_minsnr": 3.0,
    "gain_minsnr": 3.0,
    "minblperant": 4,
    "gain_solmode": "L1R",
    # The provisional pre-bandpass solve must be receptor-dependent. CASA's
    # dual-polarization T solve is only valid after a prior G or unnormalised
    # bandpass has made the receptor phases coherent.
    "prebp_gaintype": "G",
    "bandpass_fillgaps": 8,
    "bandpass_solint": "inf",
    "bandpass_solnorm": True,
    # Spectral order used by fluxscale to bootstrap a transfer calibrator.
    # It is clamped at run time to (number of spectral windows - 1), because
    # a wideband single-SPW setup cannot support a spectral fit.
    "fluxscale_fitorder": 1,
    # 'flex' lets CASA re-reference a solution interval to another antenna
    # when the reference antenna is flagged in it, holding the alternate's
    # phase continuous. 'strict' instead flags every antenna's solution in
    # that interval. For the parallel-hand, total-intensity calibration this
    # pipeline does, a change of reference adds the same phase to every
    # antenna at that time, which cancels on every baseline, so nothing is
    # lost by 'flex'. 'strict' is costly: in the uGMRT Band 4 run that
    # prompted this default, the reference antenna was flagged for the first
    # phase-calibrator scan, 'strict' flagged that scan's solutions for the
    # whole array, and linear interpolation then flagged 100% of the science
    # scan next to it, a third of that target's hour-angle coverage. Set
    # 'strict' only for work that needs one physical reference antenna
    # throughout, such as cross-hand polarization calibration.
    "refantmode": "flex",
    # Solution interval of the phase table transferred to science targets.
    # 'inf' means one solution per calibrator scan, which is averaged down
    # before being interpolated across the calibrator cycle.
    "target_phase_transfer_solint": "inf",
    # Fail closed when more than half of a solver table is flagged. This is
    # checked element-by-element across polarization, channel, and row axes.
    "max_caltable_flagged_fraction": 0.50,
    # A globally healthy table can still contain an unsolved antenna, SPW,
    # field, or observation. Guard each available identifier independently.
    "max_caltable_group_flagged_fraction": 0.80,
    # A single antenna can be completely absent from a delay solution even
    # when the rest of the array has a healthy K table. Permit one bounded,
    # recorded recovery: exclude that antenna and recompute delay.cal.
    "auto_quarantine_failed_delay_antennas": True,
    "max_auto_quarantine_delay_antennas": 1,
    "max_auto_quarantine_delay_fraction": 0.10,
    "min_antennas_after_delay_quarantine": 8,
    # Reject selected data when any applied table lacks the corresponding
    # solution instead of silently leaving those visibilities uncorrected.
    "applymode": "calflagstrict",
    # Propagate amplitude-calibration factors into the visibility weights.
    # Phase-only tables do not alter weights; amplitude/bandpass tables do.
    "calwt": True,
    "parang": False,
}

TELESCOPE_CALIBRATION_TASK_CONFIG = {
    "ugmrt": {
        "flux_standard": "Perley-Butler 2017",
        "phase_gaintype": "T",
        "amp_gaintype": "T",
        # uGMRT RFI is most safely found after a provisional bandpass has
        # flattened the spectrum.  Accepted corrected-data flags therefore
        # trigger a complete fresh solve, just as they do for MeerKAT.
        "postcal_flagging": True,
        "postcal_tfcrop_timecutoff": 5.0,
        "postcal_tfcrop_freqcutoff": 4.5,
        "postcal_rflag_timedevscale": 6.0,
        "postcal_rflag_freqdevscale": 6.0,
        "postcal_rflag_winsize": 5,
        "postcal_max_new_fraction": 0.10,
        "postcal_max_total_fraction": 0.70,
        "postcal_max_group_new_fraction": 0.15,
        "postcal_max_group_total_fraction": 0.85,
        # Interference leaves wings just under the detection threshold on
        # either side in time and frequency. Extending each accepted flag
        # over its neighbours removes them before the solutions are redone.
        # The pass never crosses a scan boundary.
        "postcal_extend": True,
        "postcal_extend_growtime": 80.0,
        "postcal_extend_growfreq": 80.0,
        # GMRT data carry one weight per row for every channel, while a
        # channel's calibrated noise follows its SEFD and residual
        # interference, which vary across Band 4. statwt measures that noise
        # and writes per-channel weights (Jy^-2) to each calibrated target.
        "statwt_enabled": True,
        "statwt_timebin": "120s",
        "statwt_chanbin": 16,
        "statwt_minsamp": 16,
        # applycal flags a target visibility when either of the two
        # calibrator solutions bracketing it failed, so one failed solution
        # removes an antenna from the science scans on both sides. For the
        # targets only, a failed solution with good ones of the same antenna
        # no further apart than this is bridged from them (or, at the start
        # or end of the run, taken from the nearest good one within half of
        # it). 120 minutes spans one missed visit of a calibrator observed
        # every hour. An antenna with fewer than this fraction of its
        # solutions good is treated as failing and left flagged.
        "target_transfer_bridge": True,
        "target_transfer_bridge_max_gap_minutes": 120.0,
        "target_transfer_bridge_min_valid_fraction": 0.5,
    },
    "meerkat": {
        "flux_standard": "Stevens-Reynolds 2016",
        # Stokes-I parallel-hand calibration in the linear-feed basis.  The
        # full-polarization XY/leakage workflow remains outside this pipeline.
        "phase_gaintype": "T",
        "amp_gaintype": "T",
        "postcal_flagging": True,
        "postcal_tfcrop_timecutoff": 4.0,
        "postcal_tfcrop_freqcutoff": 4.0,
        "postcal_rflag_timedevscale": 5.0,
        "postcal_rflag_freqdevscale": 5.0,
        "postcal_rflag_winsize": 5,
        "postcal_max_new_fraction": 0.20,
        "postcal_max_total_fraction": 0.75,
        "postcal_max_group_new_fraction": 0.30,
        "postcal_max_group_total_fraction": 0.90,
        # Interference leaves wings just under the detection threshold on
        # either side in time and frequency. Extending each accepted flag
        # over its neighbours removes them before the solutions are redone.
        # The pass never crosses a scan boundary.
        "postcal_extend": True,
        "postcal_extend_growtime": 80.0,
        "postcal_extend_growfreq": 80.0,
    },
    "lofar": {
        "flux_standard": "Scaife-Heald 2012",
    },
}


# ============================================================
# SELF-CALIBRATION
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
    },
}

SELFCAL_CONFIG.update({
    "MEERKAT_UHF": {
        "phase_solints": ["4min", "2min", "1min", "int"],
        "amp_solint": "inf",
        "min_snr_phase": 5.0,
        "min_snr_amp": 40.0,
    },
    "MEERKAT_L": {
        "phase_solints": ["4min", "2min", "1min", "int"],
        "amp_solint": "inf",
        "min_snr_phase": 5.0,
        "min_snr_amp": 40.0,
    },
    **{
        band: {
            "phase_solints": ["2min", "1min", "30s", "int"],
            "amp_solint": "inf",
            "min_snr_phase": 5.0,
            "min_snr_amp": 50.0,
        }
        for band in MEERKAT_S_BANDS
    },
    # These values describe downstream metadata only.  The capability gate
    # prevents the CASA self-calibration implementation from using them.
    "LOFAR_LBA": {
        "phase_solints": ["2min", "1min", "30s"],
        "amp_solint": "inf",
        "min_snr_phase": 8.0,
        "min_snr_amp": 50.0,
    },
    "LOFAR_HBA": {
        "phase_solints": ["4min", "2min", "1min"],
        "amp_solint": "inf",
        "min_snr_phase": 8.0,
        "min_snr_amp": 50.0,
    },
})

# Every band's field and pixel follow one rule (finding 167): four or five
# pixels across the beam of the longest baselines at the band's centre, and
# a field reaching out to where the primary beam at the band's lowest
# frequency has fallen to about a tenth, as Band 4's below does (0.12 at 550
# MHz). Band 2's 1.0" cells were 13 pixels to its 13" beam over a field
# whose edge still had 0.80 of the beam at 125 MHz, Band 3's 0.8" 8 pixels
# with 0.54, Band 5's 0.3" x 8192 0.26 at 1050 MHz. Bands 4 and 5 and
# MeerKAT L are as their real runs were made; the others are not yet run.
SELFCAL_IMAGING_CONFIG = {
    "B2": {
        "cell": "3.0arcsec",
        "imsize": 8192,
        "robust": 0.5,
        "nterms": 2,
        "niter": 100000,
        "nsigma": 3.0,
    },
    "B3": {
        "cell": "1.5arcsec",
        "imsize": 8192,
        "robust": 0.5,
        "nterms": 2,
        "niter": 100000,
        "nsigma": 3.0,
    },
    # 1.0" cells put 4 to 5 pixels across the 4" to 5" Band 4 beam, which is
    # what CLEAN needs; 0.5" spent four times the pixels on 9 per beam. The
    # field is what matters: 8192 x 0.5" reached only 34 arcmin, where the
    # measured beam is still 0.11 at 700 MHz and 0.28 at 550 MHz, so sources
    # beyond it were never in the model. 5120 x 1.0" reaches 43 arcmin
    # (beam 0.03 and 0.12) in 39% of the pixels.
    "B4": {
        "cell": "1.0arcsec",
        "imsize": 5120,
        "robust": 0.5,
        "nterms": 2,
        "niter": 100000,
        "nsigma": 3.0,
    },
    # As the Band 5 targets were imaged (2026-11, run_targets.py): 0.5" is
    # four pixels to the 2" beam, and the 51' field ends where the beam at
    # 1050 MHz is 0.12.
    "B5": {
        "cell": "0.5arcsec",
        "imsize": 6144,
        "robust": 0.5,
        "nterms": 2,
        "niter": 100000,
        "nsigma": 3.0,
    },
    # MeerKAT UHF and S by the same rule; L keeps its band profile's grid.
    "MEERKAT_UHF": {
        "cell": "2.0arcsec",
        "imsize": 7168,
        "robust": -0.5,
        "nterms": 2,
        "niter": 100000,
        "nsigma": 3.0,
    },
    **{
        band: {
            "cell": "0.6arcsec",
            "imsize": 8192,
            "robust": -0.5,
            "nterms": 2,
            "niter": 100000,
            "nsigma": 3.0,
        }
        for band in MEERKAT_S_BANDS
    },
}

SELFCAL_IMAGING_CONFIG.update({
    band: {
        "cell": values["cell"],
        "imsize": values["imsize"],
        "robust": -0.5 if band.startswith("MEERKAT_") else -0.25,
        "nterms": 2,
        "niter": 100000,
        "nsigma": 3.0,
    }
    for band, values in ALL_BANDS.items()
    if band not in SELFCAL_IMAGING_CONFIG
})

SELFCAL_TASK_CONFIG = {
    # Science-safe auto-multithresh defaults.  The initial mask is deliberately
    # resistant to sidelobes/noise; a bounded retry is allowed only when a
    # strong residual source is present but the model stayed empty.
    "automask_sidelobethreshold": 2.5,
    "automask_noisethreshold": 4.25,
    "automask_lownoisethreshold": 1.5,
    "automask_minbeamfrac": 0.20,
    "automask_growiterations": 75,
    "automask_negativethreshold": 0.0,
    "automask_smoothfactor": 1.0,
    "automask_cutthreshold": 0.01,
    "automask_dogrowprune": True,
    "automask_minpercentchange": -1.0,
    "automask_fastnoise": False,
    "automask_retry_enabled": True,
    # With poor uv coverage the sidelobe criterion of auto-multithresh can
    # sit above the image peak, so no mask opens and CLEAN produces no model
    # at all, however bright the field. The bounded relaxed-mask retry is
    # what recovers such a field, so it triggers on any feature CLEAN could
    # legitimately deconvolve rather than only on a strong one.
    "automask_retry_min_peak_snr": 5.0,
    "automask_retry_sidelobethreshold": 1.5,
    "automask_retry_noisethreshold": 3.5,
    "automask_retry_minbeamfrac": 0.05,
    "automask_retry_growiterations": 150,
    # A snapshot with sparse uv coverage has a point-spread function whose
    # sidelobes reach tens of percent of the peak. auto-multithresh then
    # opens no mask at all, because its sidelobe criterion sits at or above
    # the image peak, and CLEAN returns an empty model however bright the
    # sources are: the published image stays the dirty map, and
    # self-calibration has nothing to solve against. When that happens, the
    # compact emission inside the primary beam is cleaned directly, to a
    # conservative depth and a capped iteration count so the pass cannot
    # chase sidelobes down into the noise.
    "automask_fallback_enabled": True,
    "automask_fallback_min_peak_snr": 5.0,
    "automask_fallback_nsigma": 6.0,
    "automask_fallback_niter": 5000,
    "automask_fallback_pbmask": 0.2,
    # auto-multithresh measures its mask threshold on the residual, which on
    # the first pass is the dirty map. Where several comparable sources make
    # that map's rms sidelobe confusion rather than noise, the threshold sits
    # at or above the peak and no mask ever opens. A bounded first clean of
    # what is plainly there removes those sidelobes, so the rms means what it
    # says by the time auto-multithresh reads it. nsigma is re-measured each
    # major cycle, so this cannot run away: it reaches the brightest peaks
    # first and follows the noise down.
    "clean_bootstrap_enabled": True,
    # Shallow on purpose. The bootstrap only has to take out enough bright
    # emission for the residual rms to mean something; auto-multithresh then
    # masks every source above its noise threshold and cleans it to nsigma.
    # A 5-sigma bootstrap left nothing above that threshold, so the masked
    # pass found no mask and each source kept up to 5 sigma of its PSF.
    "clean_bootstrap_nsigma": 8.0,
    "clean_bootstrap_niter": 20000,
    # None: reach as far as the deconvolution does (pblimit), so a bright
    # source in the outer primary beam is taken out with the rest.
    "clean_bootstrap_pbmask": None,
    # auto-multithresh masks emission above sidelobethreshold x the PSF's
    # largest sidelobe x the residual peak. A short observation's PSF has
    # sidelobes of 0.4 or more, which puts that above the peak itself. The
    # threshold is lowered, per image, so the criterion is at most this
    # fraction of the peak.
    "automask_max_sidelobe_fraction": 0.6,
    # If auto-multithresh still opens no mask, what remains above this many
    # sigma is cleaned without one, so an image is never left shallower
    # than a plain threshold clean would make it.
    "clean_floor_nsigma": 5.0,
    # A CLEAN that stops at its iteration limit with emission above this
    # many sigma in the residual has not finished; it is continued, at most
    # this many further runs of niter.
    "clean_continue_min_peak_snr": 5.0,
    "clean_max_continuations": 3,
    # Fewer pixels than this across the restoring beam's minor axis leave
    # every point source between pixel centres, which CLEAN's components
    # cannot represent: errors around each bright source that no cleaning
    # removes.
    "min_pixels_per_beam": 3.0,
    # A residual peak at least this significant is located and explained.
    "residual_peak_report_min_snr": 10.0,
    # An image with no CLEAN model is the sky convolved with the array's
    # point-spread function. Publishing one as final is refused when the
    # residual holds a peak above this, because the emission is real and
    # something failed to reach it.
    "undeconvolved_max_peak_snr": 5.0,
    # A residual this many times the thermal floor measured by the
    # self-calibration stage is limited by errors in the data, not by noise.
    "noise_over_thermal_warn": 3.0,
    "clean_gain": 0.1,
    "clean_cyclefactor": 1.5,
    "clean_minpsffraction": 0.05,
    "clean_maxpsffraction": 0.8,
    "clean_smallscalebias": 0.0,
    # Self-calibration images define the sky model the gains are solved
    # against.  CASA blanks, and does not deconvolve, everything below
    # pblimit of its own beam model, which for uGMRT is an Airy disk for a
    # uniformly illuminated 45 m dish: it reaches 0.05 at 31 arcmin at
    # 700 MHz, where the measured beam is still 0.16, and 0.34 at 550 MHz.
    # Sources there are bright enough to matter and were left out of the
    # model, their power going into the gain solutions. A negative value switches the blanking off,
    # so the whole image is deconvolved.
    "pblimit": -0.01,
    "minblperant": 4,
    "gain_minsnr": 3.0,
    "phase_improvement_stop": 0.02,
    # Each phase loop solves an increment on data already corrected by the
    # loops before it. Where an increment cannot be solved (too few
    # baselines, S/N below minsnr), calonly leaves those data with the
    # calibration they already have; calflagstrict would flag them, and at
    # short solution intervals on a faint field that discards a large part
    # of the observation to gain nothing.
    "applymode": "calonly",
    # An antenna with less than this fraction of the median antenna's
    # unflagged data cannot be solved and is left out of the solution
    # coverage test, instead of failing the table for every other antenna.
    "min_antenna_data_fraction": 0.20,
    # Antennas (and solution intervals) with data whose solutions still
    # fail are tolerated up to these fractions; beyond them the table is
    # rejected. A spectral window, field or observation with no solutions is
    # never tolerated.
    "max_caltable_failed_antenna_fraction": 0.25,
    "max_caltable_failed_time_fraction": 0.25,
    # Amplitude self-calibration can change the flux scale. Enable it only
    # after phase-only self-calibration and the model have been reviewed.
    "enable_amp_selfcal": False,
    "phase_gaintype": "G",
    "amp_gaintype": "G",
    "phase_refantmode": "flex",
    "phase_solmode": "",
    # Do not confuse image dynamic range with gaincal per-solution S/N.
    # min_snr_phase / min_snr_amp in SELFCAL_CONFIG are gaincal minsnr values;
    # these two are the image-quality thresholds for attempting each mode.
    # Amplitude self-calibration absorbs sky-model error into the gains, so it
    # requires a genuinely well-modelled field rather than a mere detection.
    "phase_start_dynamic_range": 20.0,
    # Phase self-calibration also starts below that dynamic range when each
    # antenna's solution at the longest solint is expected to reach this S/N.
    # Dynamic range alone is the wrong gate: gain errors are what hold it
    # down. gaincal's own minsnr and the per-loop image test still apply.
    "phase_start_min_solution_snr": 10.0,
    "amp_start_dynamic_range": 100.0,
    # One conservative recovery attempt for weak/sparse target models.
    "phase_retry_enabled": True,
    "phase_retry_solint": "inf",
    # One solution per scan. Combining scans averages over the hours that
    # separate a target's visits, across which antenna phases are unrelated.
    "phase_retry_combine": "",
    "phase_retry_minsnr": 3.0,
    "phase_retry_minblperant": 3,
    "phase_retry_gaintype": "T",
    "model_check_rows": 20000,
    "model_nonzero_epsilon": 1.0e-12,

    # Residual RFI that survives calibration is what puts stripes and
    # ripples into the final image. It is flagged on CORRECTED minus MODEL,
    # so real source structure is in the model and is not flagged, using the
    # sky model from the last self-calibration image. The pass runs once,
    # after the last accepted solution and before the final image, and is
    # undone in full if it removes more than residual_flag_max_new_fraction
    # of the data.
    # The vector-averaged visibility is the flux at the phase centre, and
    # its uncertainty is the per-visibility noise divided by sqrt(N). Below
    # this many sigma there is no significant coherent flux, which is
    # inconclusive: calibration that did not transfer looks the same as a
    # field with little flux at its centre.
    "coherence_min_significance": 5.0,

    "residual_flagging": True,
    "residual_flag_timedevscale": 5.0,
    "residual_flag_freqdevscale": 5.0,
    "residual_flag_winsize": 5,
    "residual_flag_max_new_fraction": 0.10,
    # Narrow-band and short-lived interference leaves wings that sit just
    # below the detection threshold. Extending each flag over neighbouring
    # channels and times removes them; the pass never crosses a scan.
    "residual_flag_extend": True,
    "residual_flag_growtime": 80.0,
    "residual_flag_growfreq": 80.0,
    # A self-calibration solution table is accepted only when the global
    # solution coverage and every antenna/SPW group stay inside these
    # fail-closed limits.  Non-finite and zero-valued unflagged solutions
    # count as unusable, not merely as warnings.
    "max_caltable_unusable_fraction": 0.50,
    "max_caltable_group_unusable_fraction": 0.80,
    "calwt": False,
    "parang": False,
}


# ============================================================
# FINAL CONTINUUM IMAGING
# ============================================================

FINAL_IMAGING_CONFIG = {
    # Fields and pixels by the rule given at SELFCAL_IMAGING_CONFIG.
    "B2": {
        "cell": "3.0arcsec",
        "imsize": 8192,
        "robust": 0.5,
        "nterms": 2,
        "niter": 200000,
        "nsigma": 3.0,
        "scales": [0, 5, 15, 45, 135],
        "common_beam": "",
        "wprojplanes": 128,
        "subbands": 6,
        "subband_niter": 50000,
    },
    "B3": {
        "cell": "1.5arcsec",
        "imsize": 8192,
        "robust": 0.5,
        "nterms": 2,
        "niter": 200000,
        "nsigma": 3.0,
        "scales": [0, 5, 15, 45, 135],
        "common_beam": "",
        "wprojplanes": 128,
        "subbands": 6,
        "subband_niter": 50000,
    },
    # Same geometry as the self-calibration images (see there). Scales are
    # in pixels: the beam, three beams and nine beams. The 135-pixel scale
    # was dropped with the change to 1.0" cells, where it would be 2.25
    # arcmin; multi-scale cost grows as the square of the number of scales.
    "B4": {
        "cell": "1.0arcsec",
        "imsize": 5120,
        "robust": 0.5,
        "nterms": 2,
        "niter": 200000,
        "nsigma": 3.0,
        "scales": [0, 5, 15, 45],
        "common_beam": "",
        "wprojplanes": 128,
        "subbands": 4,
        "subband_niter": 50000,
    },
    "B5": {
        "cell": "0.5arcsec",
        "imsize": 6144,
        "robust": 0.5,
        "nterms": 2,
        "niter": 200000,
        "nsigma": 3.0,
        "scales": [0, 5, 15, 45, 135],
        "common_beam": "",
        "wprojplanes": 96,
        "subbands": 4,
        "subband_niter": 50000,
    },
}

FINAL_IMAGING_CONFIG.update({
    "MEERKAT_UHF": {
        "cell": "2.0arcsec",
        "imsize": 7168,
        "robust": -0.5,
        "nterms": 2,
        "niter": 150000,
        "nsigma": 3.0,
        "scales": [0, 5, 10, 20, 40],
        "common_beam": "",
        "wprojplanes": 256,
        "subbands": 6,
        "subband_niter": 40000,
        # WSClean joins this many output channels in one deconvolution and
        # fits a spectrum (nterms) across them. Without it WSClean makes one
        # flat-spectrum image and ignores nterms, which for a band spanning a
        # factor of two biases both the model and the restored fluxes. The
        # same count gives the sub-band images, each primary-beam corrected
        # at its own frequency.
        "wsclean_channels_out": 6,
    },
    "MEERKAT_L": {
        "cell": "1.5arcsec",
        "imsize": 6144,
        "robust": -0.5,
        "nterms": 2,
        "niter": 150000,
        "nsigma": 3.0,
        "scales": [0, 5, 10, 20, 40],
        "common_beam": "",
        "wprojplanes": 512,
        "subbands": 8,
        "subband_niter": 40000,
        "wsclean_channels_out": 8,
    },
    **{
        band: {
            "cell": "0.6arcsec",
            "imsize": 8192,
            "robust": -0.5,
            "nterms": 2,
            "niter": 150000,
            "nsigma": 3.0,
            "scales": [0, 5, 10, 20, 40],
            "common_beam": "",
            "wprojplanes": 512,
            "subbands": 8,
            "subband_niter": 40000,
        }
        for band in MEERKAT_S_BANDS
    },
    "LOFAR_LBA": {
        "cell": "4.5arcsec",
        "imsize": 8192,
        "robust": -0.25,
        "nterms": 2,
        "niter": 200000,
        "nsigma": 3.0,
        "scales": [0, 4, 12, 36, 108],
        "common_beam": "",
        "wprojplanes": 256,
        "subbands": 8,
        "subband_niter": 50000,
        # LOFAR continuum bands are wide in fractional terms (LBA and HBA
        # both reach tens of per cent), so the PSF, the station beam, and the
        # sky spectrum all vary appreciably across the band.  Deconvolving
        # the whole band as one channel biases both the model and the
        # restored flux scale.  Joined-channel imaging with a spectral fit is
        # the standard LOFAR recipe, and WSClean does apply the primary-beam
        # correction to the joined MFS product, publishing
        # <prefix>-MFS-image-pb.fits alongside the per-channel images.
        "wsclean_channels_out": 6,
    },
    "LOFAR_HBA": {
        "cell": "2.0arcsec",
        "imsize": 8192,
        "robust": -0.25,
        "nterms": 2,
        "niter": 200000,
        "nsigma": 3.0,
        "scales": [0, 4, 12, 36, 108],
        "common_beam": "",
        "wprojplanes": 256,
        "subbands": 8,
        "subband_niter": 50000,
        "wsclean_channels_out": 6,
    },
})

# Self-calibration images must use the same wide-field gridding as final
# imaging.  CASA's tclean default is wprojplanes=1, which silently disables
# w-projection even when gridder='wproject' is requested; the resulting
# non-coplanar distortion away from the phase centre would corrupt the sky
# model that the self-calibration gains are solved against.
for _band, _selfcal_imaging in SELFCAL_IMAGING_CONFIG.items():
    _selfcal_imaging.setdefault(
        "wprojplanes",
        FINAL_IMAGING_CONFIG.get(_band, {}).get("wprojplanes", 128),
    )
del _band, _selfcal_imaging


FINAL_IMAGING_TASK_CONFIG = {
    "pbmin": DEFAULT_PBMIN,
    # WSClean sizes itself from the machine's memory, not from a Slurm
    # cgroup, so inside a job with a lower limit it plans for memory it may
    # not use. When the job has a limit, WSClean is capped at this fraction
    # of it (-abs-mem). 0 disables the cap; with no job limit none is added.
    "wsclean_memory_fraction": 0.8,
    # How the MeerKAT MFS image is divided by the primary beam. A band that
    # spans a factor of two has a beam that changes by tens of per cent
    # across it. band_average divides by the beam averaged over the band,
    # which is exact for a flat-spectrum source; reference_frequency divides
    # by the beam at the image's reference frequency, as the CASA backend
    # does, and leaves a frequency-dependent error off-axis.
    "meerkat_mfs_beam_mode": "band_average",
    # tclean masks .image wherever CASA's OWN primary-beam model drops below
    # pblimit (CASA default 0.2).  Where this pipeline corrects with a
    # measured beam instead (the uGMRT polynomial, or katbeam), that CASA
    # model is not the correction model, so leaving pblimit at 0.2 lets an
    # unvalidated beam decide the published aperture.  The measured uGMRT
    # main lobe is 7-20% wider than a plain 45 m diffraction estimate, and
    # the margin between CASA's 0.2 contour and this pipeline's 0.3 contour
    # is only about 15% -- thin enough that a narrower CASA model would blank
    # part of the intended science region.  Pushing the internal mask well
    # outside the science cutoff makes impbcor's pbmin the single authority.
    # A negative value disables CASA's internal mask entirely, and with it
    # the limit on where CLEAN may place components. That is the default:
    # CASA's uGMRT beam model (a 45 m Airy disk) reaches 0.05 at 31 arcmin
    # at 700 MHz, where the measured beam is still 0.16 (0.34 at 550 MHz),
    # so a positive value left real sources in the outer beam
    # undeconvolved, with their sidelobes across the published region.
    "pblimit_measured_beam": -0.01,
    # Where CASA's own beam is also the correction model, keep CASA's default.
    "pblimit_casa_beam": 0.2,
    "weighting": "briggs",
    "specmode": "mfs",
    "deconvolver": "mtmfs",
    "default_gridder": "wproject",
    "widefield_gridder": "wproject",
    "dd_gridder": "awproject",
    "usemask": "auto-multithresh",
    # Begin with a conservative mask that rejects sidelobes and isolated
    # noise.  If the model remains empty despite a strong residual source,
    # imaging.py performs one bounded, explicitly reported recovery pass.
    "automask_sidelobethreshold": 2.5,
    "automask_noisethreshold": 4.25,
    "automask_lownoisethreshold": 1.5,
    "automask_minbeamfrac": 0.20,
    "automask_growiterations": 75,
    "automask_negativethreshold": 0.0,
    "automask_smoothfactor": 1.0,
    "automask_cutthreshold": 0.01,
    "automask_dogrowprune": True,
    "automask_minpercentchange": -1.0,
    # CASA's improved (rather than legacy fast) residual-noise estimate is
    # slower but safer for fields containing extended emission.
    "automask_fastnoise": False,
    "automask_retry_enabled": True,
    # With poor uv coverage the sidelobe criterion of auto-multithresh can
    # sit above the image peak, so no mask opens and CLEAN produces no model
    # at all, however bright the field. The bounded relaxed-mask retry is
    # what recovers such a field, so it triggers on any feature CLEAN could
    # legitimately deconvolve rather than only on a strong one.
    "automask_retry_min_peak_snr": 5.0,
    "automask_retry_sidelobethreshold": 1.5,
    "automask_retry_noisethreshold": 3.5,
    "automask_retry_minbeamfrac": 0.05,
    "automask_retry_growiterations": 150,
    # A snapshot with sparse uv coverage has a point-spread function whose
    # sidelobes reach tens of percent of the peak. auto-multithresh then
    # opens no mask at all, because its sidelobe criterion sits at or above
    # the image peak, and CLEAN returns an empty model however bright the
    # sources are: the published image stays the dirty map, and
    # self-calibration has nothing to solve against. When that happens, the
    # compact emission inside the primary beam is cleaned directly, to a
    # conservative depth and a capped iteration count so the pass cannot
    # chase sidelobes down into the noise.
    "automask_fallback_enabled": True,
    "automask_fallback_min_peak_snr": 5.0,
    "automask_fallback_nsigma": 6.0,
    "automask_fallback_niter": 5000,
    "automask_fallback_pbmask": 0.2,
    # Earlier major cycles reduce minor-cycle divergence around high uGMRT
    # PSF sidelobes without changing the requested clean depth.
    # auto-multithresh measures its mask threshold on the residual, which on
    # the first pass is the dirty map. Where several comparable sources make
    # that map's rms sidelobe confusion rather than noise, the threshold sits
    # at or above the peak and no mask ever opens. A bounded first clean of
    # what is plainly there removes those sidelobes, so the rms means what it
    # says by the time auto-multithresh reads it. nsigma is re-measured each
    # major cycle, so this cannot run away: it reaches the brightest peaks
    # first and follows the noise down.
    "clean_bootstrap_enabled": True,
    # Shallow on purpose. The bootstrap only has to take out enough bright
    # emission for the residual rms to mean something; auto-multithresh then
    # masks every source above its noise threshold and cleans it to nsigma.
    # A 5-sigma bootstrap left nothing above that threshold, so the masked
    # pass found no mask and each source kept up to 5 sigma of its PSF.
    "clean_bootstrap_nsigma": 8.0,
    "clean_bootstrap_niter": 20000,
    # None: reach as far as the deconvolution does (pblimit), so a bright
    # source in the outer primary beam is taken out with the rest.
    "clean_bootstrap_pbmask": None,
    # auto-multithresh masks emission above sidelobethreshold x the PSF's
    # largest sidelobe x the residual peak. A short observation's PSF has
    # sidelobes of 0.4 or more, which puts that above the peak itself. The
    # threshold is lowered, per image, so the criterion is at most this
    # fraction of the peak.
    "automask_max_sidelobe_fraction": 0.6,
    # If auto-multithresh still opens no mask, what remains above this many
    # sigma is cleaned without one, so an image is never left shallower
    # than a plain threshold clean would make it.
    "clean_floor_nsigma": 5.0,
    # A CLEAN that stops at its iteration limit with emission above this
    # many sigma in the residual has not finished; it is continued, at most
    # this many further runs of niter.
    "clean_continue_min_peak_snr": 5.0,
    "clean_max_continuations": 3,
    # Fewer pixels than this across the restoring beam's minor axis leave
    # every point source between pixel centres, which CLEAN's components
    # cannot represent: errors around each bright source that no cleaning
    # removes.
    "min_pixels_per_beam": 3.0,
    # A residual peak at least this significant is located and explained.
    "residual_peak_report_min_snr": 10.0,
    # An image with no CLEAN model is the sky convolved with the array's
    # point-spread function. Publishing one as final is refused when the
    # residual holds a peak above this, because the emission is real and
    # something failed to reach it.
    "undeconvolved_max_peak_snr": 5.0,
    # A residual this many times the thermal floor measured by the
    # self-calibration stage is limited by errors in the data, not by noise.
    "noise_over_thermal_warn": 3.0,
    "clean_gain": 0.1,
    "clean_cyclefactor": 1.5,
    "clean_minpsffraction": 0.05,
    "clean_maxpsffraction": 0.8,
    "clean_smallscalebias": 0.0,
    "savemodel": "none",
    # In-band sub-band imaging is opt-in. It multiplies imaging time and disk
    # use by roughly the number of sub-bands; enable it with `make_subbands:
    # true` in config.yaml or `--make-subbands` on the command line.
    "make_subbands": False,
    "restart": True,
    "parallel": False,
    "require_pbcor": True,
    # CASA lists standalone widebandpbcor as a known issue when the image and
    # MS live in different directories.  The stable default corrects Taylor-0
    # intensity with the tclean PB through impbcor.  Advanced users may opt
    # into the staged wideband path after validating it with their CASA build.
    "wideband_pbcor_enabled": False,
    # Never sample the fully flagged filter-rolloff channels when explicitly
    # enabling standalone widebandpbcor.
    "wideband_pb_sample_edge_guard": 0.10,
}

TELESCOPE_IMAGING_TASK_CONFIG = {
    "ugmrt": {
        "primary_beam_provider": "ugmrt_polynomial",
        "require_pbcor": True,
    },
    "meerkat": {
        "primary_beam_provider": "katbeam",
        "require_pbcor": True,
    },
    "lofar": {
        "primary_beam_provider": "wsclean_everybeam",
        "require_pbcor": True,
        "default_gridder": "wsclean",
    },
}


# ============================================================
# SOURCE DETECTION
# ============================================================

# PyBDSF runs under PIPELINE_RUN_CONFIG["python"], not inside CASA.  Science
# catalogues use the PB-corrected image so off-axis integrated/peak fluxes are
# on the corrected intensity scale.  Adaptive local-RMS estimation handles
# the position-dependent noise inside the PB cutoff.
SOURCE_DETECTION_CONFIG = {
    "enabled": True,
    "engine": "pybdsf",
    "strict": False,
    "image_preference": "pbcor",
    "thresh_isl": 4.0,
    "thresh_pix": 5.0,
    "adaptive_rms_box": True,
    "rms_box": None,
    "rms_map": None,
    "mean_map": "default",
    # Wavelet decomposition helps recover diffuse/extended emission but can
    # substantially increase runtime for very large images.  Set False for a
    # faster compact-source catalogue.
    "atrous_do": True,
    "flag_maxsize_bm": 50.0,
    "export_mean": False,
    "quiet": True,
}


# ============================================================
# DIRECTION-DEPENDENT CALIBRATION
# ============================================================

DD_CONFIG = {
    "B2": {
        "min_flux_jy": 0.30,
        "max_directions": 12,
        "merge_radius_arcmin": 6.0,
        "solint_time": "10min",
        "solint_freq": "8MHz",
        "facet_radius_arcmin": 25.0,
    },
    "B3": {
        "min_flux_jy": 0.20,
        "max_directions": 12,
        "merge_radius_arcmin": 5.0,
        "solint_time": "8min",
        "solint_freq": "8MHz",
        "facet_radius_arcmin": 20.0,
    },
    "B4": {
        "min_flux_jy": 0.10,
        "max_directions": 10,
        "merge_radius_arcmin": 4.0,
        "solint_time": "5min",
        "solint_freq": "16MHz",
        "facet_radius_arcmin": 15.0,
    },
    "B5": {
        "min_flux_jy": 0.08,
        "max_directions": 8,
        "merge_radius_arcmin": 3.0,
        "solint_time": "3min",
        "solint_freq": "16MHz",
        "facet_radius_arcmin": 10.0,
    },
}

DD_CONFIG.update({
    "MEERKAT_UHF": {
        "min_flux_jy": 0.08,
        "max_directions": 12,
        "merge_radius_arcmin": 4.0,
        "solint_time": "5min",
        "solint_freq": "16MHz",
        "facet_radius_arcmin": 20.0,
    },
    "MEERKAT_L": {
        "min_flux_jy": 0.05,
        "max_directions": 10,
        "merge_radius_arcmin": 3.0,
        "solint_time": "3min",
        "solint_freq": "32MHz",
        "facet_radius_arcmin": 12.0,
    },
    **{
        band: {
            "min_flux_jy": 0.03,
            "max_directions": 8,
            "merge_radius_arcmin": 2.0,
            "solint_time": "2min",
            "solint_freq": "32MHz",
            "facet_radius_arcmin": 8.0,
        }
        for band in MEERKAT_S_BANDS
    },
    "LOFAR_LBA": {
        "min_flux_jy": 0.50,
        "max_directions": 20,
        "merge_radius_arcmin": 8.0,
        "solint_time": "2min",
        "solint_freq": "2MHz",
        "facet_radius_arcmin": 35.0,
    },
    "LOFAR_HBA": {
        "min_flux_jy": 0.15,
        "max_directions": 20,
        "merge_radius_arcmin": 5.0,
        "solint_time": "2min",
        "solint_freq": "2MHz",
        "facet_radius_arcmin": 25.0,
    },
})

DD_TASK_CONFIG = {
    "mode": "auto",
    "engine": "auto",
    "data_column": "auto",
    "model_column": "MODEL_DATA",
    "single_max_directions": 3,
    "survey_min_directions": 4,
    "pybdsf_thresh_isl": 4.0,
    "pybdsf_thresh_pix": 5.0,
    "pybdsf_adaptive_rms_box": True,
    "pybdsf_atrous_do": True,
    "killms_pol_mode": "Scalar",
    "killms_weighting": "Natural",
    # None derives the number of solution bands from the inspected bandwidth
    # and the band-dependent solint_freq above.
    "killms_nchan_solutions": None,
    "killms_max_nchan_solutions": 16,
    "ddfacet_major_cycles": 5,
    "ddfacet_rms_factor": 3.0,
    "ddfacet_mask_sigma": 5.0,
    "ddfacet_min_facet_diameter_deg": 0.05,
    # Preserve the global sky-model flux scale while retaining relative
    # time/frequency/antenna DD amplitude corrections.
    "ddfacet_global_norm": "MeanAbs",
    "validate_solver_cli": True,
    "strict": True,
    "execute": False,
    "post_dd_imaging": True,
}


# Every per-band imaging value that a flat config.yaml key replaced.  A flat
# key applies to *all* bands, so it silently discards telescope-specific
# tuning (MeerKAT's negative Briggs robust, for example).  Recording the
# replacements lets the imaging stages report what was overridden instead of
# leaving the difference invisible.
FLAT_IMAGING_OVERRIDES = []


def _apply_config_yaml_overrides():

    if not CONFIG_YAML:
        return

    for key in [
        "ms",
        "fits_format",
        "imported_ms",
        "gmrt_flagfiles",
        "gmrt_listscan",
        "gmrt_gvfits",
        "keep_lta_fits",
        "telescope",
        "receiver_band",
        "data_state",
        "imager",
        "wsclean",
        "wsclean_data_column",
        "lofar_beam_mode",
        "flux_standard",
        "flux_model",
        "flux_model_type",
        "target",
        "targets",
        "target_selection",
        "flux_calibrator",
        "bandpass_calibrator",
        "secondary_flux_calibrator",
        "phase_calibrator",
        "target_phase_map",
        "target_fields",
        "reference_antenna",
        "outdir",
        "copy_input_ms",
        "working_ms",
        "stages",
        "include_ddcal",
        "dd_mode",
        "dd_engine",
        "solver_command_template",
        "apply_command_template",
        "dd_prepare_only",
        "skip_dd_post_imaging",
        "casa",
        "casa_module",
        "casa_python",
        "python",
        "ncpu",
        "profile",
        "threads",
        "mpi",
        "mpi_processes",
        "mpicasa",
        "selfcal_ms",
        "imaging_ms",
        "dd_ms",
        "resume",
        "continue_on_error",
        "progress_interval_s",
        "flagdata_memory_limit_gb",
        "flagdata_memory_fraction",
        "casa_memory_fraction",
        "flagdata_bytes_per_element",
    ]:
        if key in CONFIG_YAML:
            PIPELINE_RUN_CONFIG[key] = CONFIG_YAML[key]

    if "ram_gb" in CONFIG_YAML:
        PIPELINE_RUN_CONFIG["ram_gb"] = CONFIG_YAML["ram_gb"]
        ILIFU_SLURM_CONFIG["mem"] = (
            f"{CONFIG_YAML['ram_gb']}G"
        )

    if "ncpu" in CONFIG_YAML:
        ILIFU_SLURM_CONFIG["cpus_per_task"] = CONFIG_YAML["ncpu"]

    # Slurm settings can be configured without editing this Python file. Use
    # either ``slurm_partition: my-partition`` (recommended) or the matching
    # short key (``partition: my-partition``) in config.yaml. Nested
    # ``slurm:`` mappings are supported when PyYAML is available.
    slurm_values = CONFIG_YAML.get("slurm", {})

    if not isinstance(slurm_values, dict):

        slurm_values = {}

    slurm_keys = {
        "job_name",
        "account",
        "partition",
        "qos",
        "nodes",
        "ntasks",
        "cpus_per_task",
        "mem",
        "time",
        "log_dir",
        "output",
        "error",
        "mail_user",
        "mail_type",
        "workdir",
        "module_commands",
        "activate_commands",
        "pre_commands",
        "post_commands",
        "set_thread_env",
        "pipeline_args",
        "script_path",
    }

    for key in slurm_keys:

        prefixed = f"slurm_{key}"

        if prefixed in CONFIG_YAML:

            ILIFU_SLURM_CONFIG[key] = CONFIG_YAML[prefixed]

        elif key in CONFIG_YAML:

            ILIFU_SLURM_CONFIG[key] = CONFIG_YAML[key]

        elif key in slurm_values:

            ILIFU_SLURM_CONFIG[key] = slurm_values[key]

    flagging_yaml_keys = {
        "flagging_max_iterations": "max_iterations",
        "flagging_precal_rflag": "precal_rflag",
        "flagging_bad_ant_threshold": "bad_ant_threshold",
        "flagging_refant_preference_tolerance": (
            "refant_preference_tolerance"
        ),
        "flagging_quack_interval_s": "quack_interval_s",
        "flagging_flag_autocorrelations": "flag_autocorrelations",
        "flagging_flag_shadowed": "flag_shadowed",
        "flagging_shadow_tolerance_m": "shadow_tolerance_m",
        "flagging_extend_growtime": "extend_growtime",
        "flagging_extend_growfreq": "extend_growfreq",
        "flagging_extend_polarizations": "extend_polarizations",
        "flagging_convergence_new_fraction": "convergence_new_fraction",
        "flagging_max_iteration_new_fraction": (
            "max_iteration_new_fraction"
        ),
        "flagging_max_total_fraction": "max_total_flagged_fraction",
        "flagging_max_group_iteration_new_fraction": (
            "max_group_iteration_new_fraction"
        ),
        "flagging_max_group_total_fraction": (
            "max_group_total_flagged_fraction"
        ),
        "flagging_max_extension_new_fraction": (
            "max_extension_new_fraction"
        ),
        "flagging_max_group_extension_new_fraction": (
            "max_group_extension_new_fraction"
        ),
    }

    for yaml_key, config_key in flagging_yaml_keys.items():
        if yaml_key in CONFIG_YAML:
            FLAGGING_CONFIG[config_key] = CONFIG_YAML[yaml_key]

    # Observation-specific flagdata commands (mode='manual'), applied before
    # the automatic flaggers. Absent unless given, so the flagging recipe of
    # runs without any is unchanged.
    manual_flags = parse_manual_flag_commands(
        CONFIG_YAML.get("flagging_manual_commands")
    )
    if manual_flags:
        FLAGGING_CONFIG["manual_commands"] = manual_flags

    edge_fractions = FLAGGING_CONFIG.setdefault(
        "edge_fraction",
        {},
    )
    for band in ALL_BANDS:
        yaml_key = f"flagging_edge_fraction_{band.lower()}"
        if yaml_key in CONFIG_YAML:
            edge_fractions[band] = float(CONFIG_YAML[yaml_key])

    if "meerkat_l_known_frequency_masks" in CONFIG_YAML:
        masks = CONFIG_YAML["meerkat_l_known_frequency_masks"]
        if isinstance(masks, str):
            masks = [item.strip() for item in masks.split(",") if item.strip()]
        FLAGGING_CONFIG.setdefault("known_frequency_masks", {})[
            "MEERKAT_L"
        ] = list(masks)

    postcal_keys = {
        "postcal_tfcrop_timecutoff": "postcal_tfcrop_timecutoff",
        "postcal_tfcrop_freqcutoff": "postcal_tfcrop_freqcutoff",
        "postcal_rflag_timedevscale": "postcal_rflag_timedevscale",
        "postcal_rflag_freqdevscale": "postcal_rflag_freqdevscale",
        "postcal_rflag_winsize": "postcal_rflag_winsize",
        "postcal_max_new_fraction": "postcal_max_new_fraction",
        "postcal_max_total_fraction": "postcal_max_total_fraction",
        "postcal_max_group_new_fraction": (
            "postcal_max_group_new_fraction"
        ),
        "postcal_max_group_total_fraction": (
            "postcal_max_group_total_fraction"
        ),
    }
    for telescope in ("ugmrt", "meerkat", "lofar"):
        telescope_config = TELESCOPE_CALIBRATION_TASK_CONFIG[telescope]
        for suffix in ("statwt_enabled", "statwt_timebin", "statwt_chanbin",
                       "statwt_minsamp", "target_transfer_bridge",
                       "target_transfer_bridge_max_gap_minutes",
                       "target_transfer_bridge_min_valid_fraction"):
            yaml_key = f"{telescope}_{suffix}"
            if yaml_key in CONFIG_YAML:
                telescope_config[suffix] = CONFIG_YAML[yaml_key]
        if telescope == "lofar":
            continue
        for suffix, config_key in postcal_keys.items():
            yaml_key = f"{telescope}_{suffix}"
            if yaml_key in CONFIG_YAML:
                telescope_config[config_key] = CONFIG_YAML[yaml_key]

    autoflag_yaml_keys = {
        "flagging_tfcrop_timecutoff": "tfcrop_timecutoff",
        "flagging_tfcrop_freqcutoff": "tfcrop_freqcutoff",
        "flagging_tfcrop_maxnpieces": "tfcrop_maxnpieces",
        "flagging_rflag_timedevscale": "rflag_timedevscale",
        "flagging_rflag_freqdevscale": "rflag_freqdevscale",
        "flagging_rflag_winsize": "rflag_winsize",
        "flagging_tighten_per_iteration": "tighten_per_iteration",
    }
    default_autoflag = FLAGGING_CONFIG.setdefault(
        "autoflag",
        {},
    ).setdefault(
        "default",
        {},
    )

    for yaml_key, config_key in autoflag_yaml_keys.items():
        if yaml_key in CONFIG_YAML:
            default_autoflag[config_key] = CONFIG_YAML[yaml_key]

    if "flagging_allow_previously_flagged_input" in CONFIG_YAML:
        FLAGGING_CONFIG["allow_previously_flagged_input"] = bool(
            CONFIG_YAML["flagging_allow_previously_flagged_input"]
        )

    if "flagging_ml_enabled" in CONFIG_YAML:
        FLAGGING_CONFIG.setdefault(
            "ml_rfi",
            {},
        )["enabled"] = bool(CONFIG_YAML["flagging_ml_enabled"])

    imaging_keys = [
        "cell",
        "imsize",
        "robust",
        "niter",
        "nterms",
        "nsigma",
        "common_beam",
        "subbands",
        "subband_niter",
    ]

    for key in imaging_keys:
        if key not in CONFIG_YAML:
            continue

        value = CONFIG_YAML[key]

        for table_name, table in (
            ("FINAL_IMAGING_CONFIG", FINAL_IMAGING_CONFIG),
            ("SELFCAL_IMAGING_CONFIG", SELFCAL_IMAGING_CONFIG),
        ):
            if table is SELFCAL_IMAGING_CONFIG and key not in {
                "cell",
                "imsize",
                "robust",
                "niter",
                "nterms",
                "nsigma",
            }:
                continue

            for band, config in table.items():
                previous = config.get(key)
                if previous is not None and previous != value:
                    FLAT_IMAGING_OVERRIDES.append({
                        "table": table_name,
                        "band": band,
                        "key": key,
                        "band_default": previous,
                        "config_yaml_value": value,
                    })
                config[key] = value

    for key in [
        "pbmin",
        "pblimit_measured_beam",
        "pblimit_casa_beam",
        "restart",
        "parallel",
        "make_subbands",
        "wideband_pbcor_enabled",
        "wideband_pb_sample_edge_guard",
    ]:
        if key in CONFIG_YAML:
            FINAL_IMAGING_TASK_CONFIG[key] = CONFIG_YAML[key]

    if "calibration_calwt" in CONFIG_YAML:
        CALIBRATION_TASK_CONFIG["calwt"] = bool(
            CONFIG_YAML["calibration_calwt"]
        )

    if "calibration_prebp_gaintype" in CONFIG_YAML:
        prebp_gaintype = str(
            CONFIG_YAML["calibration_prebp_gaintype"]
        ).strip().upper()
        if prebp_gaintype != "G":
            raise ValueError(
                "calibration_prebp_gaintype must be G for the supported "
                "dual-polarization Stokes-I calibration chain"
            )
        CALIBRATION_TASK_CONFIG["prebp_gaintype"] = prebp_gaintype

    for yaml_key, config_key in {
        "calibration_fluxscale_fitorder": "fluxscale_fitorder",
        "calibration_refantmode": "refantmode",
        "calibration_target_phase_transfer_solint": (
            "target_phase_transfer_solint"
        ),
    }.items():
        if yaml_key not in CONFIG_YAML:
            continue
        value = CONFIG_YAML[yaml_key]
        if config_key == "fluxscale_fitorder":
            value = int(value)
        else:
            value = str(value).strip()
        if config_key == "refantmode" and value not in {"strict", "flex"}:
            raise ValueError(
                "calibration_refantmode must be 'strict' or 'flex'"
            )
        CALIBRATION_TASK_CONFIG[config_key] = value

    # One value for every band of the run's telescope; 0 keeps the native
    # channels. Whatever is set becomes part of the calibration recipe.
    if "calibration_target_max_channel_width_khz" in CONFIG_YAML:
        _max_width = float(
            CONFIG_YAML["calibration_target_max_channel_width_khz"]
        )
        if not math.isfinite(_max_width) or _max_width < 0:
            raise ValueError(
                "calibration_target_max_channel_width_khz must be a "
                "finite number of kHz, 0 or more"
            )
        for _band_config in CAL_CONFIG.values():
            _band_config["target_max_channel_width_khz"] = _max_width

    if "calibration_max_flagged_solution_fraction" in CONFIG_YAML:
        CALIBRATION_TASK_CONFIG["max_caltable_flagged_fraction"] = float(
            CONFIG_YAML["calibration_max_flagged_solution_fraction"]
        )

    if "calibration_max_group_flagged_solution_fraction" in CONFIG_YAML:
        CALIBRATION_TASK_CONFIG[
            "max_caltable_group_flagged_fraction"
        ] = float(
            CONFIG_YAML[
                "calibration_max_group_flagged_solution_fraction"
            ]
        )

    if "calibration_auto_quarantine_failed_delay_antennas" in CONFIG_YAML:
        CALIBRATION_TASK_CONFIG[
            "auto_quarantine_failed_delay_antennas"
        ] = bool(
            CONFIG_YAML[
                "calibration_auto_quarantine_failed_delay_antennas"
            ]
        )

    if "calibration_max_auto_quarantine_delay_antennas" in CONFIG_YAML:
        CALIBRATION_TASK_CONFIG[
            "max_auto_quarantine_delay_antennas"
        ] = int(
            CONFIG_YAML[
                "calibration_max_auto_quarantine_delay_antennas"
            ]
        )

    if "calibration_max_auto_quarantine_delay_fraction" in CONFIG_YAML:
        CALIBRATION_TASK_CONFIG[
            "max_auto_quarantine_delay_fraction"
        ] = float(
            CONFIG_YAML[
                "calibration_max_auto_quarantine_delay_fraction"
            ]
        )

    # A calibrator given its flux density outright, with CASA setjy's
    # manual model: in uGMRT Band 5 (2021-02-22) the transfer from 3C48
    # made 0323+055 2.30 Jy where RACS-mid, five months later, had 3.08 Jy
    # at 1250 MHz, and every target came out a quarter low with it. Named
    # as a secondary flux calibrator with a survey's flux, it becomes an
    # absolute reference of its own. Set only when given, so that runs
    # without it keep their recorded calibration recipe.
    if "calibration_calibrator_flux_models" in CONFIG_YAML:
        CALIBRATION_TASK_CONFIG["calibrator_flux_models"] = (
            parse_calibrator_flux_models(
                CONFIG_YAML["calibration_calibrator_flux_models"]
            )
        )

    if "calibration_min_antennas_after_delay_quarantine" in CONFIG_YAML:
        CALIBRATION_TASK_CONFIG[
            "min_antennas_after_delay_quarantine"
        ] = int(
            CONFIG_YAML[
                "calibration_min_antennas_after_delay_quarantine"
            ]
        )

    for yaml_key, config_key in {
        "selfcal_phase_start_dynamic_range": "phase_start_dynamic_range",
        "selfcal_phase_start_min_solution_snr": "phase_start_min_solution_snr",
        "selfcal_amp_start_dynamic_range": "amp_start_dynamic_range",
        "selfcal_pblimit": "pblimit",
        "selfcal_min_antenna_data_fraction": "min_antenna_data_fraction",
        "selfcal_max_failed_antenna_fraction": (
            "max_caltable_failed_antenna_fraction"
        ),
        "selfcal_max_failed_time_fraction": (
            "max_caltable_failed_time_fraction"
        ),
        "selfcal_amp_improvement_stop": "amp_improvement_stop",
        "selfcal_amp_max_peak_change": "amp_max_peak_change",
        # gaincal's minsnr for the amp+phase solve (default 10, in
        # selfcal.py). The bands' min_snr_amp is no longer used for it.
        "selfcal_amp_minsnr": "amp_minsnr",
    }.items():
        if yaml_key in CONFIG_YAML:
            SELFCAL_TASK_CONFIG[config_key] = float(CONFIG_YAML[yaml_key])

    # Amplitude self-calibration is off unless asked for. When on, its
    # solution is applied to a copy and kept only if the image improves and
    # the peak stays put (see SelfCalManager.judge_amp_selfcal).
    if "selfcal_enable_amp" in CONFIG_YAML:
        _enable_amp = CONFIG_YAML["selfcal_enable_amp"]
        if isinstance(_enable_amp, str):
            _enable_amp = _enable_amp.strip().lower() in {"true", "yes", "1"}
        SELFCAL_TASK_CONFIG["enable_amp_selfcal"] = bool(_enable_amp)

    if "selfcal_flag_unsolvable_antennas" in CONFIG_YAML:
        _flag_unsolvable = CONFIG_YAML["selfcal_flag_unsolvable_antennas"]
        if isinstance(_flag_unsolvable, str):
            _flag_unsolvable = (
                _flag_unsolvable.strip().lower() in {"true", "yes", "1"}
            )
        SELFCAL_TASK_CONFIG["flag_unsolvable_antennas"] = bool(
            _flag_unsolvable
        )

    # G solves each polarization's amplitude, T both together (default, in
    # selfcal.py), which a faint field needs for the signal-to-noise.
    if "selfcal_amp_gaintype" in CONFIG_YAML:
        _amp_gaintype = str(CONFIG_YAML["selfcal_amp_gaintype"]).strip().upper()
        if _amp_gaintype not in {"G", "T"}:
            raise ValueError("selfcal_amp_gaintype must be G or T")
        SELFCAL_TASK_CONFIG["amp_solve_gaintype"] = _amp_gaintype

    # Multi-scale sizes, in pixels, for the self-calibration images of every
    # band; without it they are cleaned with point components only.
    if "selfcal_scales" in CONFIG_YAML:
        _scales = CONFIG_YAML["selfcal_scales"]
        if isinstance(_scales, str):
            _scales = [
                item for item in _scales.strip("[] ").split(",") if item.strip()
            ]
        _scales = [int(float(item)) for item in (_scales or [])]
        if any(item < 0 for item in _scales):
            raise ValueError("selfcal_scales must be pixel sizes of 0 or more")
        for _band_config in SELFCAL_IMAGING_CONFIG.values():
            _band_config["scales"] = list(_scales)

    if "selfcal_applymode" in CONFIG_YAML:
        _applymode = str(CONFIG_YAML["selfcal_applymode"]).strip()
        if _applymode not in {"calonly", "calflag", "calflagstrict"}:
            raise ValueError(
                "selfcal_applymode must be calonly, calflag or "
                f"calflagstrict, not {CONFIG_YAML['selfcal_applymode']!r}"
            )
        SELFCAL_TASK_CONFIG["applymode"] = _applymode

    # The masked clean's major-cycle limit defaults in selfcal.py and
    # imaging.py (DEFAULT_CLEAN_MAX_MAJOR_CYCLES), not here, so that adding
    # it left the recipe of a stage already under way unchanged. Setting it
    # here does change the recipe, as any setting would.
    # clean_max_major_cycles sets both stages; the prefixed keys set one.
    for _yaml_key, _targets in (
        ("clean_max_major_cycles",
         (SELFCAL_TASK_CONFIG, FINAL_IMAGING_TASK_CONFIG)),
        ("selfcal_clean_max_major_cycles", (SELFCAL_TASK_CONFIG,)),
        ("imaging_clean_max_major_cycles", (FINAL_IMAGING_TASK_CONFIG,)),
    ):
        if _yaml_key not in CONFIG_YAML:
            continue
        _major_cycles = int(CONFIG_YAML[_yaml_key])
        if _major_cycles < 0:
            raise ValueError(
                f"{_yaml_key} must be 0 (no limit) or a positive number of "
                "major cycles"
            )
        for _target in _targets:
            _target["clean_max_major_cycles"] = _major_cycles

    # Self-calibration images only have to carry a sky model, so they can be
    # made smaller than the published image when memory is tight: the padded
    # convolution grid, and with it CASA's largest allocation, grows as the
    # square of imsize. Reducing these narrows the field the sky model
    # covers, so the published image is left at its own size.
    for yaml_key, config_key, cast in (
        ("selfcal_imsize", "imsize", int),
        ("selfcal_cell", "cell", str),
        ("selfcal_wprojplanes", "wprojplanes", int),
        ("selfcal_niter", "niter", int),
        ("selfcal_robust", "robust", float),
    ):
        if yaml_key not in CONFIG_YAML:
            continue
        try:
            value = cast(CONFIG_YAML[yaml_key])
        except (TypeError, ValueError):
            raise ValueError(
                f"{yaml_key} must be a {cast.__name__}, not "
                f"{CONFIG_YAML[yaml_key]!r}"
            )
        if cast is int and value <= 0:
            raise ValueError(f"{yaml_key} must be positive, not {value}")
        for _band, _band_settings in SELFCAL_IMAGING_CONFIG.items():
            previous = _band_settings.get(config_key)
            if previous is not None and previous != value:
                # Recorded like any other displaced band default, so the
                # self-calibration report says the model was built over a
                # smaller field than the band expects.
                FLAT_IMAGING_OVERRIDES.append({
                    "table": "SELFCAL_IMAGING_CONFIG",
                    "band": _band,
                    "key": config_key,
                    "band_default": previous,
                    "config_yaml_value": value,
                })
            _band_settings[config_key] = value

    if "selfcal_max_unusable_solution_fraction" in CONFIG_YAML:
        SELFCAL_TASK_CONFIG["max_caltable_unusable_fraction"] = float(
            CONFIG_YAML["selfcal_max_unusable_solution_fraction"]
        )

    if "selfcal_max_group_unusable_solution_fraction" in CONFIG_YAML:
        SELFCAL_TASK_CONFIG[
            "max_caltable_group_unusable_fraction"
        ] = float(
            CONFIG_YAML[
                "selfcal_max_group_unusable_solution_fraction"
            ]
        )

    shared_clean_keys = [
        "automask_sidelobethreshold",
        "automask_noisethreshold",
        "automask_lownoisethreshold",
        "automask_minbeamfrac",
        "automask_growiterations",
        "automask_negativethreshold",
        "automask_smoothfactor",
        "automask_cutthreshold",
        "automask_dogrowprune",
        "automask_minpercentchange",
        "automask_fastnoise",
        "automask_retry_enabled",
        "automask_retry_min_peak_snr",
        "automask_retry_sidelobethreshold",
        "automask_retry_noisethreshold",
        "automask_retry_minbeamfrac",
        "automask_retry_growiterations",
        "automask_fallback_enabled",
        "automask_fallback_min_peak_snr",
        "automask_fallback_nsigma",
        "automask_fallback_niter",
        "automask_fallback_pbmask",
        "clean_bootstrap_enabled",
        "clean_bootstrap_nsigma",
        "clean_bootstrap_niter",
        "clean_bootstrap_pbmask",
        "automask_max_sidelobe_fraction",
        "clean_floor_nsigma",
        "clean_continue_min_peak_snr",
        "clean_max_continuations",
        "min_pixels_per_beam",
        "residual_peak_report_min_snr",
        "undeconvolved_max_peak_snr",
        "noise_over_thermal_warn",
        "clean_gain",
        "clean_cyclefactor",
        "clean_minpsffraction",
        "clean_maxpsffraction",
        "clean_smallscalebias",
    ]

    for key in shared_clean_keys:
        if key in CONFIG_YAML:
            FINAL_IMAGING_TASK_CONFIG[key] = CONFIG_YAML[key]
            SELFCAL_TASK_CONFIG[key] = CONFIG_YAML[key]

    if "wsclean_memory_fraction" in CONFIG_YAML:
        try:
            fraction = float(CONFIG_YAML["wsclean_memory_fraction"])
        except (TypeError, ValueError):
            fraction = None
        if fraction is None or not 0.0 <= fraction <= 1.0:
            raise ValueError(
                "wsclean_memory_fraction must be between 0 and 1, not "
                f"{CONFIG_YAML['wsclean_memory_fraction']!r}"
            )
        FINAL_IMAGING_TASK_CONFIG["wsclean_memory_fraction"] = fraction

    if "meerkat_mfs_beam_mode" in CONFIG_YAML:
        mode = str(CONFIG_YAML["meerkat_mfs_beam_mode"]).strip().lower()
        if mode not in {"band_average", "reference_frequency"}:
            raise ValueError(
                "meerkat_mfs_beam_mode must be band_average or "
                f"reference_frequency, not {mode!r}"
            )
        FINAL_IMAGING_TASK_CONFIG["meerkat_mfs_beam_mode"] = mode

    if "meerkat_wsclean_channels_out" in CONFIG_YAML:
        try:
            channels = int(CONFIG_YAML["meerkat_wsclean_channels_out"])
        except (TypeError, ValueError):
            channels = -1
        if channels < 0:
            raise ValueError(
                "meerkat_wsclean_channels_out must be a non-negative "
                "integer (0 or 1 images the band as one channel)."
            )
        for band in ("MEERKAT_UHF", "MEERKAT_L"):
            FINAL_IMAGING_CONFIG[band]["wsclean_channels_out"] = channels

    source_detection_yaml_keys = {
        "source_detection_enabled": "enabled",
        "source_detection_strict": "strict",
        "source_detection_image_preference": "image_preference",
        "pybdsf_thresh_isl": "thresh_isl",
        "pybdsf_thresh_pix": "thresh_pix",
        "pybdsf_adaptive_rms_box": "adaptive_rms_box",
        "pybdsf_atrous_do": "atrous_do",
        "pybdsf_rms_box": "rms_box",
        "pybdsf_rms_map": "rms_map",
        "pybdsf_mean_map": "mean_map",
        "pybdsf_flag_maxsize_bm": "flag_maxsize_bm",
    }

    for yaml_key, config_key in source_detection_yaml_keys.items():
        if yaml_key in CONFIG_YAML:
            SOURCE_DETECTION_CONFIG[config_key] = CONFIG_YAML[yaml_key]

    dd_yaml_keys = {
        "dd_mode": "mode",
        "dd_engine": "engine",
        "dd_data_column": "data_column",
        "dd_single_max_directions": "single_max_directions",
        "dd_survey_min_directions": "survey_min_directions",
        "dd_killms_pol_mode": "killms_pol_mode",
        "dd_killms_weighting": "killms_weighting",
        "dd_killms_nchan_solutions": "killms_nchan_solutions",
        "dd_killms_max_nchan_solutions": "killms_max_nchan_solutions",
        "ddfacet_major_cycles": "ddfacet_major_cycles",
        "ddfacet_rms_factor": "ddfacet_rms_factor",
        "ddfacet_mask_sigma": "ddfacet_mask_sigma",
        "ddfacet_min_facet_diameter_deg": (
            "ddfacet_min_facet_diameter_deg"
        ),
        "ddfacet_global_norm": "ddfacet_global_norm",
        "dd_validate_solver_cli": "validate_solver_cli",
        "dd_strict": "strict",
    }

    for yaml_key, config_key in dd_yaml_keys.items():
        if yaml_key in CONFIG_YAML:
            DD_TASK_CONFIG[config_key] = CONFIG_YAML[yaml_key]

    dd_band_yaml_keys = {
        "dd_min_flux_jy": "min_flux_jy",
        "dd_max_directions": "max_directions",
        "dd_merge_radius_arcmin": "merge_radius_arcmin",
        "dd_facet_radius_arcmin": "facet_radius_arcmin",
        "dd_solint_time": "solint_time",
        "dd_solint_freq": "solint_freq",
    }

    for yaml_key, config_key in dd_band_yaml_keys.items():
        if yaml_key not in CONFIG_YAML:
            continue
        for band_config in DD_CONFIG.values():
            band_config[config_key] = CONFIG_YAML[yaml_key]


_apply_config_yaml_overrides()


# ============================================================
# PRODUCT TREE
# ============================================================

PRODUCT_DIRS = {
    "ingest": {
        "root": "ingest",
    },
    "inspection": {
        "root": "inspection",
        "plots": "inspection/plots",
    },
    "flagging": {
        "root": "flagging",
        "plots": "flagging/plots",
    },
    "calibration": {
        "root": "calibration",
        "plots": "calibration/plots",
        "calibrated": "calibration/calibrated",
    },
    "selfcal": {
        "root": "selfcal",
        "plots": "selfcal/plots",
        "images": "selfcal/images",
        "caltables": "selfcal/caltables",
        "checkpoints": "selfcal/checkpoints",
        "reports": "selfcal/reports",
        "selfcal_ms": "selfcal/selfcal_ms",
    },
    "imaging": {
        "root": "imaging",
        "images": "imaging/images",
        "subbands": "imaging/subbands",
        "fits": "imaging/fits",
        "plots": "imaging/plots",
        "reports": "imaging/reports",
        "checkpoints": "imaging/checkpoints",
    },
    "source_detection": {
        "root": "source_detection",
        "catalogs": "source_detection/catalogs",
        "images": "source_detection/images",
        "logs": "source_detection/logs",
        "reports": "source_detection/reports",
        "checkpoints": "source_detection/checkpoints",
    },
    "ddcal": {
        "root": "ddcal",
        "catalogs": "ddcal/catalogs",
        "directions": "ddcal/directions",
        "configs": "ddcal/configs",
        "logs": "ddcal/logs",
        "images": "ddcal/images",
        "plots": "ddcal/plots",
        "reports": "ddcal/reports",
        "checkpoints": "ddcal/checkpoints",
        "solutions": "ddcal/solutions",
    },
}

REPORT_FILES = {
    "ingest": "ingest/ingest_report.json",
    "inspection": "inspection/inspection_report.json",
    "flagging_pre": "flagging/preflag_summary.json",
    "flagging_post": "flagging/postflag_summary.json",
    "bad_antennas": "flagging/bad_antennas.txt",
    "recommended_refant": "flagging/recommended_refant.txt",
    "calibration": "calibration/calibration_report.json",
    "split_targets": "calibration/split_targets.json",
    "selfcal": "selfcal/reports/selfcal_report.json",
    "selfcal_metrics": "selfcal/reports/selfcal_metrics.json",
    "imaging": "imaging/reports/imaging_report.json",
    "imaging_metrics": "imaging/reports/imaging_metrics.json",
    "subband_imaging": "imaging/reports/subband_imaging_report.json",
    "source_detection": "source_detection/reports/source_detection_report.json",
    "ddcal": "ddcal/reports/ddcal_report.json",
}


# ============================================================
# COMPATIBILITY ALIASES
# ============================================================

# Older helpers used this generic name. Keep it as the final science-imaging defaults, while SELFCAL_IMAGING_CONFIG remains available for selfcal.py.

IMAGING_CONFIG = FINAL_IMAGING_CONFIG

# Convenience lookup for common beam values.
COMMON_BEAMS = {
    band: config["common_beam"]
    for band, config in FINAL_IMAGING_CONFIG.items()
}


# ============================================================
# SLURM / SBATCH HELPERS
# ============================================================

def _shell_join(parts):

    return " ".join(
        shlex.quote(str(part))
        for part in parts
        if part is not None and str(part) != ""
    )


def _sbatch_line(key, value):

    if value in (None, "", [], {}):

        return None

    return f"#SBATCH --{key}={value}"


def slurm_pipeline_command(pipeline_script=None,
                           dry_run=False,
                           extra_args=None):

    if pipeline_script is None:

        pipeline_script = (
            Path(__file__).resolve().with_name(
                "pipeline.py"
            )
        )

    python_exe = (
        PIPELINE_RUN_CONFIG.get("python") or
        "python3"
    )

    command = [
        python_exe,
        str(pipeline_script),
    ]

    command.append(
        "--resume"
        if PIPELINE_RUN_CONFIG.get("resume", True)
        else "--no-resume"
    )

    if dry_run:

        command.append("--dry-run")

    command.extend(
        ILIFU_SLURM_CONFIG.get(
            "pipeline_args",
            [],
        )
    )

    if extra_args:

        command.extend(extra_args)

    return _shell_join(command)


def render_sbatch(pipeline_script=None,
                  dry_run=False,
                  extra_args=None):

    slurm = dict(ILIFU_SLURM_CONFIG)

    if pipeline_script is None:

        pipeline_script = (
            Path(__file__).resolve().with_name(
                "pipeline.py"
            )
        )

    pipeline_script = Path(pipeline_script).resolve()

    workdir = slurm.get("workdir") or str(
        pipeline_script.parent
    )

    log_dir = slurm.get("log_dir", "slurm_logs")
    output = str(Path(log_dir) / slurm.get(
        "output",
        "%x-%j.out",
    ))
    error = str(Path(log_dir) / slurm.get(
        "error",
        "%x-%j.err",
    ))

    directive_values = {
        "job-name": slurm.get("job_name"),
        "account": slurm.get("account"),
        "partition": slurm.get("partition"),
        "qos": slurm.get("qos"),
        "nodes": slurm.get("nodes"),
        "ntasks": slurm.get("ntasks"),
        "cpus-per-task": slurm.get("cpus_per_task"),
        "mem": slurm.get("mem"),
        "time": slurm.get("time"),
        "output": output,
        "error": error,
        "mail-user": slurm.get("mail_user"),
        "mail-type": (
            slurm.get("mail_type")
            if slurm.get("mail_user")
            else ""
        ),
    }

    lines = [
        "#!/bin/bash",
        "#",
        "# Auto-generated radio continuum pipeline Slurm job.",
        "# Edit config.yaml, regenerate this file, then submit with sbatch.",
    ]

    for key, value in directive_values.items():

        line = _sbatch_line(key, value)

        if line:

            lines.append(line)

    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
            f"mkdir -p {_shell_join([log_dir])}",
            f"cd {_shell_join([workdir])}",
            "",
            "echo \"Starting radio continuum pipeline on $(date)\"",
            "echo \"Running on node(s): ${SLURM_NODELIST:-unknown}\"",
            "echo \"MPI tasks: ${SLURM_NTASKS:-1}\"",
            "echo \"CPUs per task: ${SLURM_CPUS_PER_TASK:-1}\"",
            "echo \"CPUs on node: ${SLURM_CPUS_ON_NODE:-unknown}\"",
            "echo \"Submit directory: ${SLURM_SUBMIT_DIR:-unknown}\"",
            "",
        ]
    )

    for command in slurm.get("module_commands", []):

        lines.append(command)

    if slurm.get("module_commands"):

        lines.append("")

    for command in slurm.get("activate_commands", []):

        lines.append(command)

    if slurm.get("activate_commands"):

        lines.append("")

    if slurm.get("set_thread_env", True):

        lines.extend(
            [
                "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}",
                "export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}",
                "export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}",
                "export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}",
                "export VECLIB_MAXIMUM_THREADS=${SLURM_CPUS_PER_TASK:-1}",
                "export OMP_DYNAMIC=FALSE",
                "export PYTHONUNBUFFERED=1",
                "",
            ]
        )

    for command in slurm.get("pre_commands", []):

        lines.append(command)

    if slurm.get("pre_commands"):

        lines.append("")

    lines.append(
        slurm_pipeline_command(
            pipeline_script=pipeline_script,
            dry_run=dry_run,
            extra_args=extra_args,
        )
    )
    lines.append("")

    for command in slurm.get("post_commands", []):

        lines.append(command)

    if slurm.get("post_commands"):

        lines.append("")

    lines.extend(
        [
            "echo \"Finished radio continuum pipeline on $(date)\"",
            "",
        ]
    )

    return "\n".join(lines)


def write_sbatch(path=None,
                 pipeline_script=None,
                 dry_run=False,
                 extra_args=None):

    if path is None:

        path = ILIFU_SLURM_CONFIG.get(
            "script_path",
            "run_radio_continuum_pipeline.sbatch",
        )

    path = Path(path)
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    content = render_sbatch(
        pipeline_script=pipeline_script,
        dry_run=dry_run,
        extra_args=extra_args,
    )

    path.write_text(content)
    path.chmod(0o755)

    return path


# ============================================================
# CONFIG DATACLASS
# ============================================================

@dataclass
class PipelineConfig:
    ms: str
    outdir: str = DEFAULT_OUTDIR
    ncpu: int = DEFAULT_NCPU
    ram_gb: int = DEFAULT_RAM_GB
    restart: bool = True
    pbmin: float = DEFAULT_PBMIN
    stages: list = field(
        default_factory=lambda: list(PIPELINE_STAGES)
    )
    preferred_refants: list = field(
        default_factory=lambda: list(DEFAULT_REFANTS)
    )

    def __post_init__(self):

        self.outdir = str(self.outdir)
        Path(self.outdir).mkdir(
            exist_ok=True,
            parents=True,
        )

    @property
    def outpath(self):

        return Path(self.outdir)

    def product_dir(self, stage, key="root"):

        return self.outpath / PRODUCT_DIRS[stage][key]

    def report_path(self, name):

        return self.outpath / REPORT_FILES[name]

    def ensure_product_tree(self):

        for stage_dirs in PRODUCT_DIRS.values():

            for relpath in stage_dirs.values():

                (
                    self.outpath /
                    relpath
                ).mkdir(
                    parents=True,
                    exist_ok=True,
                )

    def band_from_frequency(self, frequency_hz):

        return detect_band(frequency_hz)

    def band_defaults(self, band):

        return get_band_defaults(band)


# ============================================================
# LOOKUP HELPERS
# ============================================================

def normalize_band(band, telescope="ugmrt"):

    if not band:

        return None

    try:
        return TELESCOPE_PROFILE_MODULE.normalize_band(
            telescope,
            band,
        )
    except ValueError:
        compact = re.sub(r"[^A-Z0-9]+", "", str(band).upper())
        if str(telescope).lower() == "meerkat" and compact in {"S", "SBAND"}:
            # Preserve the actionable S0-S4 ambiguity error from the
            # telescope registry instead of degrading it to a generic None.
            raise
        return None


def detect_band(frequency_hz, telescope="ugmrt"):

    try:
        return TELESCOPE_PROFILE_MODULE.detect_band(
            telescope,
            frequency_hz,
        )["name"]
    except ValueError:
        return "UNKNOWN"


def get_band_defaults(band, telescope="ugmrt"):

    band = normalize_band(band, telescope=telescope)

    if band is None:

        raise ValueError(
            f"Unsupported {telescope} band: {band}"
        )

    defaults = dict(ALL_BANDS[band])
    defaults.update(
        {
            "flagging": FLAGGING_CONFIG,
            "calibration": CAL_CONFIG[band],
            "selfcal": SELFCAL_CONFIG[band],
            "selfcal_imaging": SELFCAL_IMAGING_CONFIG[band],
            "imaging": FINAL_IMAGING_CONFIG[band],
            "ddcal": DD_CONFIG[band],
            "common_beam": COMMON_BEAMS[band],
        }
    )

    return defaults


def stage_order(include_ddcal=True):

    stages = list(PIPELINE_STAGES)

    if not include_ddcal:

        stages.remove("ddcal")

    return stages
