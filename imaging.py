#!/usr/bin/env python3
# ============================================================
#
# imaging.py
#
# CASA-native final continuum imaging stage for the uGMRT and MeerKAT
# pipeline.
#
# Intended order:
#   inspect.py -> flagging.py -> calibration.py -> selfcal.py
#   -> imaging.py
#
# Compatible with CASA 6+
#
# Author: Sushant Dutta
# ============================================================

import os
import sys

SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

if (
    sys.path and
    os.path.abspath(sys.path[0]) == SCRIPT_DIR
):

    sys.path.pop(0)

import argparse
import hashlib
import importlib.util
import json
import re
import math
import shutil
import numpy as np
from datetime import datetime, timezone
from pathlib import Path


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of failing partway through a
# stage. Increase it whenever pipeline.py starts to rely on something new here.
PIPELINE_INTERFACE_VERSION = 7

try:
    from importlib import metadata as importlib_metadata
except ImportError:  # pragma: no cover - CASA Python older than 3.8
    importlib_metadata = None


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
spectral_window_ids = CASA_COMPAT.spectral_window_ids
exportfits = CASA_COMPAT.import_casa_task("exportfits")
imhead = CASA_COMPAT.import_casa_task("imhead")
impbcor = CASA_COMPAT.import_casa_task("impbcor")
imsmooth = CASA_COMPAT.import_casa_task("imsmooth")
imstat = CASA_COMPAT.import_casa_task("imstat")
rmtables = CASA_COMPAT.import_casa_task("rmtables")
tclean = CASA_COMPAT.import_casa_task("tclean")
widebandpbcor = CASA_COMPAT.import_casa_task(
    "widebandpbcor",
    required=False,
)
from casatools import image as image_tool, msmetadata, quanta as quanta_tool, table


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
        "radio_continuum_imaging_safety",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONFIG = load_config_module()
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
TELESCOPE_IMAGING_TASK_CONFIG = (
    getattr(CONFIG, "TELESCOPE_IMAGING_TASK_CONFIG", {})
    if CONFIG is not None
    else {}
)
QA_PRODUCTS = load_qa_products_module()
TELESCOPE_PROFILE_MODULE = load_telescope_profiles_module()
IMAGING_SAFETY = load_imaging_safety_module()
FLAT_IMAGING_OVERRIDES = (
    getattr(CONFIG, "FLAT_IMAGING_OVERRIDES", [])
    if CONFIG is not None
    else []
)


def flat_imaging_overrides_for(table_name, band):
    """Return flat config.yaml keys that displaced this band's own default."""

    return [
        item
        for item in FLAT_IMAGING_OVERRIDES
        if item.get("table") == table_name and item.get("band") == band
    ]


def resolve_execution_profile(profile):
    """Resolve a profile without mutating shared image defaults."""

    name = str(profile or RUN_CONFIG.get("profile", "science"))
    name = name.strip().lower()

    if CONFIG is not None and hasattr(CONFIG, "get_execution_profile"):
        return name, CONFIG.get_execution_profile(name)

    if name == "science":
        return name, {}

    raise ValueError(
        "Execution profiles require config.py; use the 'science' profile "
        "when running imaging.py as a standalone script without it."
    )


# ============================================================
# BAND CONFIGURATION
# ============================================================

IMAGING_CONFIG = {

    "B2": {
        "cell": "1.0arcsec",
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
        "cell": "0.8arcsec",
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

    "B4": {
        "cell": "0.5arcsec",
        "imsize": 8192,
        "robust": 0.5,
        "nterms": 2,
        "niter": 200000,
        "nsigma": 3.0,
        "scales": [0, 5, 15, 45, 135],
        "common_beam": "",
        "wprojplanes": 128,
        "subbands": 4,
        "subband_niter": 50000,
    },

    "B5": {
        "cell": "0.3arcsec",
        "imsize": 8192,
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


if CONFIG is not None:

    IMAGING_CONFIG = getattr(
        CONFIG,
        "FINAL_IMAGING_CONFIG",
        IMAGING_CONFIG,
    )


# ============================================================
# HELPERS
# ============================================================

def load_json(path, default=None):

    path = Path(path)

    if not path.exists():

        return default

    with open(path) as f:

        return json.load(f)


def to_builtin(value):

    if hasattr(value, "tolist"):

        value = value.tolist()

    if isinstance(value, dict):

        return {
            str(k): to_builtin(v)
            for k, v in value.items()
        }

    if isinstance(value, (list, tuple)):

        return [
            to_builtin(v)
            for v in value
        ]

    try:

        json.dumps(value)
        return value

    except TypeError:

        return str(value)


def casa_image_exists(path):

    if not path:

        return False

    path = Path(path)

    return bool(
        path.is_dir() and
        (path / "table.dat").is_file()
    )


def clean_name(name):

    if not name:

        return "target"

    name = str(name).strip()
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)

    return name.strip("_") or "target"


def now_utc():

    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def adaptive_automask_kwargs(config, previous_image=None, relaxed=False):
    """Return bounded uGMRT auto-multithresh parameters.

    The normal pass is deliberately resistant to sidelobes and isolated
    noise.  ``relaxed`` is reserved for one recovery pass after the pipeline
    has verified that a strong residual source exists but CLEAN wrote an
    empty model.
    """

    if relaxed:
        sidelobe = float(
            IMAGING_TASK_CONFIG.get(
                "automask_retry_sidelobethreshold",
                1.5,
            )
        )
        noise = float(
            IMAGING_TASK_CONFIG.get(
                "automask_retry_noisethreshold",
                3.5,
            )
        )
        minbeam = float(
            IMAGING_TASK_CONFIG.get(
                "automask_retry_minbeamfrac",
                0.05,
            )
        )
        growiterations = int(
            IMAGING_TASK_CONFIG.get(
                "automask_retry_growiterations",
                150,
            )
        )
    else:
        sidelobe = float(
            IMAGING_TASK_CONFIG.get(
                "automask_sidelobethreshold",
                2.5,
            )
        )
        noise = float(
            IMAGING_TASK_CONFIG.get(
                "automask_noisethreshold",
                4.25,
            )
        )
        minbeam = float(
            IMAGING_TASK_CONFIG.get(
                "automask_minbeamfrac",
                0.20,
            )
        )
        growiterations = int(
            IMAGING_TASK_CONFIG.get(
                "automask_growiterations",
                75,
            )
        )

    return {
        "usemask": IMAGING_TASK_CONFIG.get(
            "usemask",
            "auto-multithresh",
        ),
        "sidelobethreshold": max(1.0, sidelobe),
        "noisethreshold": max(3.0, noise),
        "lownoisethreshold": float(
            IMAGING_TASK_CONFIG.get(
                "automask_lownoisethreshold",
                1.5,
            )
        ),
        "minbeamfrac": max(0.01, minbeam),
        "growiterations": max(1, growiterations),
        "negativethreshold": float(
            IMAGING_TASK_CONFIG.get(
                "automask_negativethreshold",
                0.0,
            )
        ),
        "smoothfactor": float(
            IMAGING_TASK_CONFIG.get(
                "automask_smoothfactor",
                1.0,
            )
        ),
        "cutthreshold": float(
            IMAGING_TASK_CONFIG.get(
                "automask_cutthreshold",
                0.01,
            )
        ),
        "dogrowprune": bool(
            IMAGING_TASK_CONFIG.get(
                "automask_dogrowprune",
                True,
            )
        ),
        "minpercentchange": float(
            IMAGING_TASK_CONFIG.get(
                "automask_minpercentchange",
                -1.0,
            )
        ),
        "fastnoise": bool(
            IMAGING_TASK_CONFIG.get(
                "automask_fastnoise",
                False,
            )
        ),
    }


def clean_control_kwargs():
    """Numerical controls that keep minor cycles stable on uGMRT PSFs."""

    return {
        "gain": float(
            IMAGING_TASK_CONFIG.get(
                "clean_gain",
                0.1,
            )
        ),
        "cyclefactor": float(
            IMAGING_TASK_CONFIG.get(
                "clean_cyclefactor",
                1.5,
            )
        ),
        "minpsffraction": float(
            IMAGING_TASK_CONFIG.get(
                "clean_minpsffraction",
                0.05,
            )
        ),
        "maxpsffraction": float(
            IMAGING_TASK_CONFIG.get(
                "clean_maxpsffraction",
                0.8,
            )
        ),
        "smallscalebias": float(
            IMAGING_TASK_CONFIG.get(
                "clean_smallscalebias",
                0.0,
            )
        ),
    }


def wideband_pb_samples(vis, nterms, edge_guard=0.10):
    """Build valid interior SPW/channel samples for widebandpbcor.

    The flagging stage removes the filter roll-off at both edges of every
    spectral window.  Sampling channel zero or the final channel therefore
    violates widebandpbcor's requirement that data exist at every requested
    frequency.  Use interior, frequency-spanning samples instead.
    """

    edge_guard = float(edge_guard)

    if not 0.0 < edge_guard < 0.5:
        raise ValueError(
            "wideband_pb_sample_edge_guard must be greater than 0 and "
            f"less than 0.5; received {edge_guard}."
        )

    samples = []
    msmd = msmetadata()
    try:
        msmd.open(str(vis))
        for spw in spectral_window_ids(msmd):
            freqs = list(msmd.chanfreqs(int(spw)))
            if not freqs:
                continue
            last = len(freqs) - 1
            channel_samples = [
                int(round(last * edge_guard)),
                len(freqs) // 2,
                int(round(last * (1.0 - edge_guard))),
            ]
            for chan in sorted(set(channel_samples)):
                samples.append((float(freqs[chan]), int(spw), int(chan)))
    finally:
        try:
            msmd.close()
        except Exception:
            pass

    samples.sort(key=lambda row: row[0])
    minimum = max(int(nterms), 2)
    desired = max(minimum, int(nterms) + 1)
    if len(samples) < minimum:
        raise RuntimeError(
            f"Only {len(samples)} valid PB frequency samples were found; "
            f"widebandpbcor with nterms={nterms} requires at least {minimum}."
        )

    if len(samples) > desired:
        indices = [round(i * (len(samples) - 1) / (desired - 1)) for i in range(desired)]
        samples = [samples[i] for i in sorted(set(indices))]

    return (
        [row[1] for row in samples],
        [row[2] for row in samples],
        [1.0 for _ in samples],
    )


# ============================================================
# FINAL IMAGING MANAGER
# ============================================================

# Most major cycles the masked (auto-multithresh) pass of the final image may
# run; see FinalImagingManager.cap_masked_major_cycles. On the uGMRT Band 4
# target that set it, that pass took the peak residual inside its mask from
# 581 to 193 uJy/beam in eight cycles and then returned to 229-248 uJy/beam at
# the start of each of the next seven. A default here, not in config.py,
# leaves the imaging recipe of an existing output directory unchanged.
DEFAULT_CLEAN_MAX_MAJOR_CYCLES = 20


class FinalImagingManager:

    # --------------------------------------------------------
    # Constructor
    # --------------------------------------------------------

    def __init__(self,
                 ms=None,
                 outdir="Products",
                 target=None,
                 robust=None,
                 imsize=None,
                 cell=None,
                 niter=None,
                 nterms=None,
                 nsigma=None,
                 scales=None,
                 pbmin=None,
                 common_beam=None,
                 gridder=None,
                 datacolumn="auto",
                 parallel=None,
                 restart=None,
                 force=False,
                 make_subbands=None,
                 n_subbands=None,
                 subband_niter=None,
                 profile=None,
                 telescope=None,
                 receiver_band=None,
                 provenance_fingerprint=None,
                 upstream_outdir=None):

        self.explicit_ms = ms
        self.outdir = Path(outdir)
        # Inspection and calibration describe the observation, so with one
        # product directory per target they stay at the run's root.
        self.upstream_outdir = (
            Path(upstream_outdir) if upstream_outdir else self.outdir
        )
        self.target = target
        self.requested_telescope = telescope
        self.requested_band = receiver_band
        self.provenance_fingerprint = (
            str(provenance_fingerprint).strip()
            if provenance_fingerprint else None
        )
        self.profile, self.profile_settings = (
            resolve_execution_profile(profile)
        )
        profile_imaging = self.profile_settings.get("imaging", {})
        self.pbmin = (
            IMAGING_TASK_CONFIG.get(
                "pbmin",
                0.2,
            )
            if pbmin is None
            else pbmin
        )
        self.datacolumn = datacolumn
        self.parallel = (
            profile_imaging.get(
                "parallel",
                IMAGING_TASK_CONFIG.get(
                    "parallel",
                    False,
                ),
            )
            if parallel is None
            else parallel
        )
        self.restart = (
            IMAGING_TASK_CONFIG.get(
                "restart",
                True,
            )
            if restart is None
            else restart
        )
        self.force = force
        self.require_pbcor = bool(
            IMAGING_TASK_CONFIG.get(
                "require_pbcor",
                True,
            )
        )
        self.wideband_pbcor_enabled = bool(
            IMAGING_TASK_CONFIG.get(
                "wideband_pbcor_enabled",
                False,
            )
        )
        self.wideband_pb_sample_edge_guard = float(
            IMAGING_TASK_CONFIG.get(
                "wideband_pb_sample_edge_guard",
                0.10,
            )
        )

        if not 0.0 < self.wideband_pb_sample_edge_guard < 0.5:
            raise ValueError(
                "wideband_pb_sample_edge_guard must be greater than 0 "
                "and less than 0.5."
            )
        self.make_subbands = (
            profile_imaging.get(
                "make_subbands",
                IMAGING_TASK_CONFIG.get(
                    "make_subbands",
                    False,
                ),
            )
            if make_subbands is None
            else make_subbands
        )

        self.root = self.outdir / "imaging"
        self.imagedir = self.root / "images"
        self.subbanddir = self.root / "subbands"
        self.fitsdir = self.root / "fits"
        self.plotdir = self.root / "plots"
        self.reportdir = self.root / "reports"
        self.checkdir = self.root / "checkpoints"

        for d in [
            self.root,
            self.imagedir,
            self.subbanddir,
            self.fitsdir,
            self.plotdir,
            self.reportdir,
            self.checkdir,
        ]:

            d.mkdir(parents=True, exist_ok=True)

        self.inspection = load_json(
            self.upstream_outdir /
            "inspection" /
            "inspection_report.json",
            {},
        )

        self.calibration = load_json(
            self.upstream_outdir /
            "calibration" /
            "calibration_report.json",
            {},
        )

        self.selfcal = load_json(
            self.outdir /
            "selfcal" /
            "reports" /
            "selfcal_report.json",
            {},
        )

        detected_telescope = (
            TELESCOPE_PROFILE_MODULE.telescope_name_from_report(
                self.selfcal
            ) or
            TELESCOPE_PROFILE_MODULE.telescope_name_from_report(
                self.calibration
            ) or
            TELESCOPE_PROFILE_MODULE.telescope_name_from_report(
                self.inspection
            )
        )
        if self.requested_telescope not in (None, "", "auto"):
            requested = TELESCOPE_PROFILE_MODULE.normalize_telescope(
                self.requested_telescope
            )
            if detected_telescope and requested != detected_telescope:
                raise RuntimeError(
                    f"Requested telescope {requested} conflicts with the "
                    f"stage reports ({detected_telescope})."
                )
            self.telescope = requested
        else:
            self.telescope = detected_telescope
        if self.telescope is None:
            raise RuntimeError(
                "No supported telescope identity is present in the stage "
                "reports. Rerun inspection with --telescope."
            )

        if self.telescope == "lofar":
            raise RuntimeError(
                "LOFAR imaging must use wsclean_imaging.py so the station "
                "beam is handled by EveryBeam. Run pipeline.py with "
                "--imager wsclean (the LOFAR auto default)."
            )

        self.band = self.detect_band()
        self.config = dict(IMAGING_CONFIG[self.band])
        self.flat_imaging_overrides = flat_imaging_overrides_for(
            "FINAL_IMAGING_CONFIG",
            self.band,
        )

        for override in self.flat_imaging_overrides:
            print(
                "WARNING: config.yaml key "
                f"{override['key']}={override['config_yaml_value']!r} "
                f"replaced the {self.band} default "
                f"{override['band_default']!r}. A flat key overrides every "
                "band, including telescope-specific tuning."
            )

        band_profile = TELESCOPE_PROFILE_MODULE.get_band_profile(
            self.telescope,
            self.band,
        )
        telescope_task = TELESCOPE_IMAGING_TASK_CONFIG.get(
            self.telescope,
            {},
        )
        self.primary_beam_provider = band_profile.get(
            "primary_beam_provider",
            telescope_task.get("primary_beam_provider", "none"),
        )
        self.katbeam_model = band_profile.get("katbeam_model")
        self.ugmrt_pb_coefficients = band_profile.get(
            "primary_beam_coefficients"
        )
        self.ugmrt_pb_valid_frequency_mhz = band_profile.get(
            "primary_beam_valid_frequency_mhz"
        )
        self.require_pbcor = bool(
            telescope_task.get("require_pbcor", self.require_pbcor)
        )

        if self.primary_beam_provider != "casa":
            self.wideband_pbcor_enabled = False

        if self.primary_beam_provider == "katbeam":
            if make_subbands is True:
                raise ValueError(
                    "MeerKAT sub-band PB correction is not implemented in "
                    "the CASA backend. Omit --make-subbands or run "
                    "telescope-validated external sub-band imaging."
                )
            if make_subbands is None:
                self.make_subbands = False

        if self.primary_beam_provider == "none":
            if make_subbands is True:
                raise ValueError(
                    "Sub-band PB correction was requested, but this "
                    "telescope/band has no validated beam provider."
                )
            if make_subbands is None:
                self.make_subbands = False

        for key, value in profile_imaging.items():
            if key not in {"make_subbands", "parallel"}:
                self.config[key] = value

        if robust is not None:
            self.config["robust"] = robust

        if imsize is not None:
            self.config["imsize"] = imsize

        if cell is not None:
            self.config["cell"] = cell

        if niter is not None:
            self.config["niter"] = niter

        if nterms is not None:
            self.config["nterms"] = nterms

        if nsigma is not None:
            self.config["nsigma"] = nsigma

        if scales is not None:
            self.config["scales"] = scales

        if common_beam is not None:
            self.config["common_beam"] = common_beam

        if n_subbands is not None:
            self.config["subbands"] = n_subbands

        if subband_niter is not None:
            self.config["subband_niter"] = subband_niter

        self.gridder = (
            gridder
            if gridder is not None
            else self.default_gridder()
        )

        self.ms = self.select_measurement_set()
        self.target = self.select_target()
        self.datacolumn = self.select_datacolumn()

        image_root = clean_name(self.target)
        self.imagename = str(
            self.imagedir /
            f"{image_root}_final_continuum"
        )

        self.image_tt0 = (
            self.imagename +
            ".image.tt0"
        )

        self.pbcor_image = None
        self.pbcor_method = None
        self.pbcor_status = {
            "standalone_wideband_enabled": self.wideband_pbcor_enabled,
            "standalone_wideband_attempted": False,
            "standalone_wideband_preflight": None,
            "standalone_wideband_error": None,
            "reused_existing_product": False,
        }
        self.science_reffreq = None
        self.common_beam_image = None
        self.common_beam_pbcor_image = None
        self.subband_products = []
        self.clean_diagnostics = {}
        self.stale_pbcor_workdirectories = []
        self.stale_pbcor_products = []

        if self.force:

            self.clear_all_checkpoints()
            self.remove_existing_products()

    # --------------------------------------------------------
    # Reports and metadata
    # --------------------------------------------------------

    def detect_band(self):

        candidates = [
            self.selfcal.get("band"),
            self.calibration.get("band"),
            self.inspection.get("band"),
        ]

        reported_band = None

        for candidate in candidates:

            if isinstance(candidate, dict):

                candidate = (
                    candidate.get("name") or
                    candidate.get("band")
                )

            if not candidate:

                continue

            band = str(candidate).upper()
            band = band.replace("BAND-", "B")
            band = band.replace("BAND", "B")

            if band in IMAGING_CONFIG:
                reported_band = band
                break

        if self.requested_band:
            requested_band = TELESCOPE_PROFILE_MODULE.normalize_band(
                self.telescope,
                self.requested_band,
            )
            if reported_band and requested_band != reported_band:
                raise RuntimeError(
                    f"Requested band {requested_band} conflicts with the "
                    f"stage reports ({reported_band})."
                )
            return requested_band

        if reported_band:
            return reported_band

        raise RuntimeError(
            f"Could not determine a configured imaging band for "
            f"{self.telescope}. Rerun inspection with the correct "
            "--receiver-band."
        )

    # --------------------------------------------------------

    def select_target(self):

        if self.target:

            return self.target

        targets = (
            self.calibration.get("targets") or
            self.inspection.get("target_fields") or
            self.inspection.get("fields", {}).get("targets") or
            []
        )

        if isinstance(targets, str):

            return targets

        if len(targets) == 1:

            return targets[0]

        if len(targets) > 1:

            raise RuntimeError(
                "Multiple target fields are available. Pass --target "
                "explicitly before imaging."
            )

        return Path(self.ms).stem

    # --------------------------------------------------------

    def select_measurement_set(self):

        if self.explicit_ms:

            ms = Path(self.explicit_ms)

            if not ms.exists():

                raise FileNotFoundError(
                    f"Measurement Set not found:\n{ms}"
                )

            return str(ms)

        final_ms = self.selfcal.get("final_ms")

        if (
            final_ms and
            Path(final_ms).exists() and
            (
                not self.target or
                self.selfcal.get("target") == self.target
            )
        ):

            print(
                "Using final self-calibrated MS from "
                "selfcal_report.json."
            )

            return final_ms

        if final_ms and Path(final_ms).exists() and self.target:
            print(
                "Ignoring a self-calibrated MS for a different target "
                f"({self.selfcal.get('target')!r}); resolving {self.target!r} "
                "from the calibrated split products instead."
            )

        split_targets = load_json(
            self.outdir /
            "calibration" /
            "split_targets.json",
            {},
        )

        if split_targets:

            if self.target and self.target in split_targets:

                return split_targets[self.target]

            if len(split_targets) != 1:

                raise RuntimeError(
                    "Multiple calibrated target MSs are available. "
                    "Pass the desired MS or --target explicitly."
                )

            only_target = next(iter(split_targets))

            return split_targets[only_target]

        raise FileNotFoundError(
            "Could not choose an input Measurement Set. "
            "Pass an MS explicitly or run selfcal.py first."
        )

    # --------------------------------------------------------

    def select_datacolumn(self):

        if self.datacolumn != "auto":

            return self.datacolumn

        tb = table()
        opened = False

        try:

            tb.open(self.ms)
            opened = True
            columns = set(tb.colnames())

        finally:

            if opened:

                tb.close()

        ms_path = Path(self.ms)

        if (
            ms_path.name == "final_selfcal.ms" or
            ms_path.parent.name == "selfcal_ms"
        ):

            return "data"

        if "CORRECTED_DATA" in columns:

            return "corrected"

        return "data"

    # --------------------------------------------------------

    def default_gridder(self):

        return IMAGING_TASK_CONFIG.get(
            "default_gridder",
            "wproject",
        )

    # --------------------------------------------------------
    # Checkpoints
    # --------------------------------------------------------

    def checkpoint_exists(self, stage):

        return (
            self.checkdir /
            f"{stage}.done"
        ).exists()

    # --------------------------------------------------------

    def mark_checkpoint(self, stage):

        (
            self.checkdir /
            f"{stage}.done"
        ).touch()

    # --------------------------------------------------------

    def clear_all_checkpoints(self):

        for checkpoint in self.checkdir.glob("*.done"):

            checkpoint.unlink()

    # --------------------------------------------------------

    def clear_pbcor_dependent_checkpoints(self):
        """Invalidate every product derived from the corrected image."""

        for stage in IMAGING_SAFETY.PBCOR_DEPENDENT_CHECKPOINTS:
            checkpoint = self.checkdir / f"{stage}.done"
            if checkpoint.exists():
                checkpoint.unlink()

    # --------------------------------------------------------
    # Product helpers
    # --------------------------------------------------------

    def save_json(self, data, outfile):
        """Atomically replace JSON so a killed Slurm job can resume."""

        outfile = Path(outfile)
        outfile.parent.mkdir(parents=True, exist_ok=True)
        temporary = outfile.with_name(
            f".{outfile.name}.tmp.{os.getpid()}"
        )

        with open(temporary, "w") as f:

            json.dump(
                to_builtin(data),
                f,
                indent=4,
            )
            f.flush()
            os.fsync(f.fileno())

        os.replace(temporary, outfile)

    # --------------------------------------------------------

    def product_candidates(self):

        suffixes = [
            "image.tt0",
            "image.tt1",
            "image.tt0.pbcor",
            "image.tt0.pbcor.provenance.json",
            "pbcor.image.tt0",
            "pbcor.image.tt1",
            "residual.tt0",
            "model.tt0",
            "mask",
            "pb",
            "pb.tt0",
            "psf.tt0",
            "sumwt.tt0",
            "alpha",
            "alpha.error",
            "pbcor.alpha",
            "pbcor.alpha.error",
            "commonbeam.image",
            "commonbeam.image.pbcor",
            "katbeam.pb.tt0",
            "katbeam.pb.tt0.provenance.json",
            "ugmrt.pb.tt0",
            "ugmrt.pb.tt0.provenance.json",
        ]

        return [
            Path(f"{self.imagename}.{suffix}")
            for suffix in suffixes
        ]

    # --------------------------------------------------------

    def existing_products(self):

        products = [
            str(path)
            for path in self.product_candidates()
            if path.exists()
        ]

        for product in self.subband_product_candidates():

            if product.exists():

                products.append(str(product))

        seen = set()
        unique = []

        for product in products:

            if product in seen:

                continue

            seen.add(product)
            unique.append(product)

        return unique

    # --------------------------------------------------------

    def subband_product_candidates(self):

        candidates = []

        for suffix in [
            "image",
            "image.tt0",
            "image.tt1",
            "image.pbcor",
            "image.tt0.pbcor",
            "pbcor.image",
            "pbcor.image.tt0",
            "pb",
            "pb.tt0",
            "residual",
            "residual.tt0",
            "model",
            "model.tt0",
            "mask",
            "psf",
            "psf.tt0",
            "sumwt",
            "sumwt.tt0",
            "ugmrt.pb.tt0",
            "ugmrt.pb.tt0.provenance.json",
            "image.pbcor.provenance.json",
            "image.tt0.pbcor.provenance.json",
        ]:

            candidates.extend(
                self.subbanddir.glob(
                    f"*_subband_*.{suffix}"
                )
            )

        return candidates

    # --------------------------------------------------------

    def remove_existing_products(self):

        products = (
            self.product_candidates() +
            self.subband_product_candidates()
        )

        image_root = Path(self.imagename).expanduser().resolve()
        for pattern in IMAGING_SAFETY.image_root_partial_patterns(image_root):
            products.extend(image_root.parent.glob(pattern))

        for product in dict.fromkeys(Path(item) for item in products):

            if product.exists():

                print(
                    f"Removing existing imaging product: {product}"
                )

                if product.is_dir():
                    if casa_image_exists(product):
                        rmtables(
                            tablenames=str(product)
                        )
                    else:
                        shutil.rmtree(product)
                else:
                    product.unlink()

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    def print_summary(self):

        print("\n")
        print("=" * 70)
        print(
            f"{TELESCOPE_PROFILE_MODULE.get_telescope_profile(self.telescope)['label']} "
            "Final Continuum Imaging"
        )
        print("=" * 70)

        print(
            f"Measurement Set : {self.ms}"
        )

        print(
            f"Data column     : {self.datacolumn}"
        )

        print(
            f"Target          : {self.target}"
        )

        print(
            f"Band            : {self.band}"
        )

        print(
            f"PB provider     : {self.primary_beam_provider}"
        )

        print(
            f"Profile         : {self.profile}"
        )

        print(
            f"Image root      : {self.imagename}"
        )

        print(
            f"Gridder         : {self.gridder}"
        )

        print(
            f"Cell / imsize   : "
            f"{self.config['cell']} / "
            f"{self.config['imsize']}"
        )

        print(
            f"Robust          : {self.config['robust']}"
        )

        print(
            f"nterms          : {self.config['nterms']}"
        )

        print(
            f"niter / nsigma  : "
            f"{self.config['niter']} / "
            f"{self.config['nsigma']}"
        )

        print(
            f"PB cutoff       : {self.pbmin}"
        )

        print(
            f"Sub-band images : "
            f"{self.config.get('subbands', 0) if self.make_subbands else 0}"
        )

        print("=" * 70)

    # --------------------------------------------------------
    # Imaging
    # --------------------------------------------------------

    def deconvolution_diagnostics(self, imagename):
        """Measure whether CLEAN made a model and what remains in residuals."""

        def first_image(suffixes):
            for suffix in suffixes:
                candidate = Path(str(imagename) + suffix)
                if casa_image_exists(candidate):
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

    def tclean_pass(self, kwargs, label, name):
        """Run one tclean call and keep CASA's account of why it stopped.

        tclean returns its stop code, the iterations and major cycles it
        ran and the PSF's largest sidelobe. Without them a residual full of
        flux cannot be told apart from one that CLEAN reached and could not
        model, which need opposite fixes.
        """

        summary = CASA_COMPAT.summarize_tclean(tclean(**kwargs))
        passes = getattr(self, "clean_passes", None)
        if passes is None:
            passes = self.clean_passes = {}
        passes.setdefault(label, []).append(dict(summary or {}, name=name))
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

        A CLEAN that reaches its iteration limit with significant emission
        still in the residual has not finished, and what it leaves behind
        keeps its full point-spread function in the image. It continues from
        where it stopped, with its mask and model, for at most
        ``clean_max_continuations`` further runs of ``niter``.

        A pass stopped by ``clean_max_major_cycles`` is not continued: that
        limit exists for a masked pass that has stopped converging, and the
        significance tested here is the whole image's, which in a field
        limited by calibration errors is an artefact no cleaning removes.
        """

        limit = int(IMAGING_TASK_CONFIG.get("clean_max_continuations", 3))
        configured_floor = IMAGING_TASK_CONFIG.get("clean_continue_min_peak_snr", 5.0)
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
        if (
            summary and summary.get("stopcode") == 1 and
            (diagnostics.get("residual_peak_snr") or 0.0) >= floor
        ):
            print(
                f"WARNING: CLEAN of {label} is still stopping at its "
                "iteration limit with a "
                f"{diagnostics['residual_peak_snr']:.1f}-sigma residual peak "
                f"after {rounds} continuations. Raise `niter` or "
                "`clean_max_continuations`; the image is not fully "
                "deconvolved."
            )
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

        configured = IMAGING_TASK_CONFIG.get("clean_bootstrap_pbmask")
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

        fraction = IMAGING_TASK_CONFIG.get("automask_max_sidelobe_fraction", 0.6)
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

        On a short observation's point-spread function the masked pass
        oscillates about its nsigma threshold: a minor cycle cleans to the
        cycle threshold and the next major cycle finds the peak back where
        it started. Nothing but the nsigma test and niter ends such a pass.
        ``nmajor`` makes the number of major cycles a limit of its own:
        ``clean_max_major_cycles``, ``DEFAULT_CLEAN_MAX_MAJOR_CYCLES`` unless
        configured, 0 for none. Self-calibration has the same limit.
        """

        value = IMAGING_TASK_CONFIG.get(
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

        floor = float(IMAGING_TASK_CONFIG.get("clean_floor_nsigma", 5.0))
        niter = int(IMAGING_TASK_CONFIG.get("clean_bootstrap_niter", 20000))
        if (
            floor <= 0 or niter <= 0 or
            not bool(IMAGING_TASK_CONFIG.get("clean_bootstrap_enabled", True)) or
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
        max_sidelobe * peak, noisethreshold * rms)``. Both terms are measured
        on the residual, which on the first pass is the dirty map. In a field
        with several comparable sources that map's rms is sidelobe confusion,
        not noise, so the threshold sits at or above the peak, no mask opens,
        nothing is deconvolved, and what gets published is the dirty map with
        every source wearing its full point-spread function.

        Cleaning the plainly real emission first breaks the circle. ``nsigma``
        is recomputed from the residual at each major cycle, so the first
        cycle reaches only the brightest peaks, and the threshold follows the
        noise down as their sidelobes come out. The pass is bounded by nsigma
        and a capped iteration count and writes no model column;
        auto-multithresh still does the deep, structure-following clean
        afterwards, on a residual whose rms now means what it says.
        """

        if not bool(IMAGING_TASK_CONFIG.get("clean_bootstrap_enabled", True)):
            return None

        nsigma = float(IMAGING_TASK_CONFIG.get("clean_bootstrap_nsigma", 8.0))
        niter = int(IMAGING_TASK_CONFIG.get("clean_bootstrap_niter", 20000))
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
        """Run tclean and recover once from an empty high-S/N automask."""

        bootstrap = self.bootstrap_clean(kwargs, label)
        if bootstrap is not None:
            # The bootstrap made the PSF and the first residual; the
            # auto-multithresh pass continues from them and its components
            # add to the ones already found.
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
            IMAGING_TASK_CONFIG.get("automask_retry_min_peak_snr", 8.0),
            diagnostics,
        )

        should_retry = (
            bool(
                IMAGING_TASK_CONFIG.get(
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
                "WARNING: CLEAN model is empty while the residual contains "
                f"a {peak_snr:.1f}-sigma feature. Retrying once with the "
                "bounded recovery mask."
            )
            self.clear_stale_clean_mask(kwargs["imagename"])
            retry_kwargs = dict(kwargs)
            retry_kwargs.update(
                adaptive_automask_kwargs(
                    self.config,
                    relaxed=True,
                )
            )
            retry_kwargs.update({
                "restart": True,
                "calcpsf": False,
                "calcres": True,
            })
            diagnostics["recovery_attempted"] = True
            diagnostics["initial"] = dict(diagnostics)

            self.tclean_pass(retry_kwargs, label, "relaxed-mask retry")
            recovered = self.deconvolution_diagnostics(
                retry_kwargs["imagename"]
            )
            recovered["recovery_attempted"] = True
            recovered["recovery_succeeded"] = bool(
                recovered.get("model_nonzero")
            )
            recovered["initial"] = diagnostics["initial"]
            diagnostics = recovered

            if not diagnostics["recovery_succeeded"]:
                print(
                    "WARNING: The bounded automask retry still produced no "
                    "CLEAN model. The residual is retained rather than "
                    "allowing an unbounded clean into noise."
                )

        elif not diagnostics.get("model_nonzero"):
            if peak_snr is None:
                reason = "residual S/N could not be measured"
            else:
                reason = (
                    f"residual peak S/N={peak_snr:.1f} is below the "
                    f"{retry_min_snr:.1f} recovery threshold"
                )
            print(
                "WARNING: CLEAN produced an empty model and no relaxed "
                f"retry was allowed because {reason}."
            )

        if (
            not diagnostics.get("model_nonzero") and
            bool(IMAGING_TASK_CONFIG.get("automask_fallback_enabled", True)) and
            peak_snr is not None and
            peak_snr >= significance_floor(
                IMAGING_TASK_CONFIG.get("automask_fallback_min_peak_snr", 5.0),
                diagnostics,
            ) and
            self.primary_beam_image_ready(kwargs["imagename"])
        ):
            # auto-multithresh never opened a mask. Clean what is plainly
            # there instead: compact emission inside the primary beam, to a
            # conservative depth with a capped iteration count, which cannot
            # follow sidelobes into the noise.
            fallback_nsigma = float(IMAGING_TASK_CONFIG.get("automask_fallback_nsigma", 6.0))
            fallback_niter = int(IMAGING_TASK_CONFIG.get("automask_fallback_niter", 5000))
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
                "pbmask": float(IMAGING_TASK_CONFIG.get("automask_fallback_pbmask", 0.2)),
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

        self.clean_diagnostics[label] = diagnostics
        return diagnostics

    # --------------------------------------------------------

    def run_tclean(self):

        checkpoint = "tclean"

        if (
            self.checkpoint_exists(checkpoint) and
            casa_image_exists(self.image_tt0)
        ):

            print(
                "Final continuum image already exists."
            )

            self.clean_diagnostics[
                "final_continuum"
            ] = self.deconvolution_diagnostics(
                self.imagename
            )

            # A kept image is published unchanged, so it is held to the same
            # standard as one made now: an earlier run's dirty map must not
            # pass simply because it already exists.
            self.verify_deconvolution("final_continuum")

            return

        print("\nRunning final continuum tclean.")

        kwargs = {
            "vis": self.ms,
            "imagename": self.imagename,
            "datacolumn": self.datacolumn,
            "specmode": IMAGING_TASK_CONFIG.get(
                "specmode",
                "mfs",
            ),
            "deconvolver": IMAGING_TASK_CONFIG.get(
                "deconvolver",
                "mtmfs",
            ),
            "nterms": self.config["nterms"],
            "scales": self.config["scales"],
            "gridder": self.gridder,
            "weighting": IMAGING_TASK_CONFIG.get(
                "weighting",
                "briggs",
            ),
            "robust": self.config["robust"],
            "imsize": self.config["imsize"],
            "cell": self.config["cell"],
            "niter": self.config["niter"],
            "nsigma": self.config["nsigma"],
            "interactive": False,
            "restart": self.restart,
            "calcpsf": True,
            "calcres": True,
            "savemodel": "none",
            "parallel": self.parallel,
            "pblimit": self.pblimit_value(),
        }

        kwargs.update(adaptive_automask_kwargs(self.config))
        kwargs.update(clean_control_kwargs())

        self.science_reffreq = self.science_reference_frequency()

        if self.science_reffreq:
            reffreq_mhz = self.science_reffreq["reffreq_hz"] / 1e6
            kwargs["reffreq"] = f"{reffreq_mhz:.6f}MHz"
            print(
                "Setting the image reference frequency to "
                f"{reffreq_mhz:.3f} MHz, the centre of the science-valid "
                "band, rather than "
                f"{self.science_reffreq['tclean_default_reffreq_hz'] / 1e6:.3f}"
                " MHz from the full digitized span. The primary beam is "
                "evaluated at the image reference frequency."
            )

        print(
            "Auto-multithresh: "
            f"sidelobethreshold={kwargs['sidelobethreshold']}, "
            f"noisethreshold={kwargs['noisethreshold']}, "
            f"minbeamfrac={kwargs['minbeamfrac']}, "
            f"fastnoise={kwargs['fastnoise']}"
        )

        if self.gridder == "wproject":

            kwargs["wprojplanes"] = (
                self.config["wprojplanes"]
            )

        if self.gridder == "awproject":

            kwargs["wbawp"] = True
            kwargs["conjbeams"] = True
            kwargs["computepastep"] = 360.0
            kwargs["rotatepastep"] = 360.0

        self.run_tclean_with_recovery(
            kwargs,
            "final_continuum",
        )

        if not casa_image_exists(self.image_tt0):

            raise RuntimeError(
                "tclean completed without producing a usable restored "
                f"image: {self.image_tt0}"
            )

        self.verify_deconvolution("final_continuum")

        self.mark_checkpoint(checkpoint)

    @staticmethod
    def angle_in_arcsec(value, unit):
        factors = {"rad": 180.0 / math.pi * 3600.0, "deg": 3600.0,
                   "arcmin": 60.0, "arcsec": 1.0, "mas": 1.0e-3}
        factor = factors.get(str(unit).strip())
        return None if factor is None else float(value) * factor

    def direction_geometry(self, image):
        """Reference pixel and cell size (arcsec) of an image's sky axes."""

        try:
            header = imhead(imagename=str(image), mode="summary")
            refpix = [float(v) for v in list(header["refpix"])[:2]]
            units = list(header["axisunits"])[:2]
            cells = [
                abs(self.angle_in_arcsec(value, unit))
                for value, unit in zip(list(header["incr"])[:2], units)
            ]
        except Exception:
            return None
        if any(cell is None or cell <= 0 for cell in cells):
            return None
        return {"refpix": refpix, "cell_arcsec": cells, "header": header}

    def pixel_value(self, image, x, y):
        if not casa_image_exists(image):
            return None
        try:
            return float(
                imstat(imagename=str(image), box=f"{x},{y},{x},{y}")["max"][0]
            )
        except Exception:
            return None

    def describe_residual_peak(self, label):
        """Say where the largest residual is and what that implies.

        Three causes leave a strong residual, and they need different fixes:
        CLEAN was never allowed there (outside its mask), CLEAN stopped
        before reaching it, or CLEAN reached it and could not model it
        (calibration or direction-dependent errors, or a spectrum the
        Taylor terms cannot follow). Position, primary-beam level and mask
        membership tell them apart.
        """

        diagnostics = self.clean_diagnostics.get(label) or {}
        position = diagnostics.get("residual_peak_position")
        peak = diagnostics.get("residual_peak_abs_Jy_beam")
        snr = diagnostics.get("residual_peak_snr")
        floor = float(
            IMAGING_TASK_CONFIG.get("residual_peak_report_min_snr", 10.0)
        )
        if not position or not peak or snr is None or snr < floor:
            return None

        x, y = position
        geometry = self.direction_geometry(self.imagename + ".residual.tt0")
        offset = None
        if geometry:
            dx = (x - geometry["refpix"][0]) * geometry["cell_arcsec"][0]
            dy = (y - geometry["refpix"][1]) * geometry["cell_arcsec"][1]
            offset = math.hypot(dx, dy) / 60.0
        pb = (
            self.pixel_value(self.imagename + ".pb.tt0", x, y)
            if casa_image_exists(self.imagename + ".pb.tt0")
            else self.pixel_value(self.imagename + ".pb", x, y)
        )
        if pb is not None and pb >= 0.999:
            pb = None  # a flat CASA beam says nothing about position
        mask = self.pixel_value(self.imagename + ".mask", x, y)
        in_mask = None if mask is None else mask > 0.5
        stop = diagnostics.get("final_stop") or {}

        record = {
            "pixel": [int(x), int(y)],
            "peak_Jy_beam": peak,
            "snr": snr,
            "offset_arcmin": offset,
            "casa_primary_beam": pb,
            "inside_clean_mask": in_mask,
        }

        where = (
            f"{peak * 1e3:.2f} mJy/beam ({snr:.0f} sigma) at pixel "
            f"({x}, {y})"
            + (f", {offset:.1f}' from the pointing centre" if offset else "")
            + (f", primary beam {pb:.2f}" if pb is not None else "")
        )
        if in_mask is False:
            record["cause"] = "outside_mask"
            explanation = (
                "outside the clean mask, so it was never deconvolved and its "
                "sidelobes remain across the whole image. auto-multithresh "
                "did not mask it; the bootstrap reaches it only if it lies "
                "inside `clean_bootstrap_pbmask`."
            )
        elif stop.get("converged") is False:
            record["cause"] = "not_converged"
            explanation = (
                f"inside the clean mask, but CLEAN {stop.get('stop_reason')} "
                "before reaching it."
            )
        elif stop.get("converged") is None:
            record["cause"] = "unknown"
            explanation = (
                "inside the clean mask; CASA did not report why CLEAN "
                "stopped, so whether it finished is not known."
            )
        else:
            record["cause"] = "not_modelled"
            explanation = (
                "inside the clean mask and CLEAN finished, so this is "
                "emission CLEAN could not model: residual calibration or "
                "direction-dependent errors, or a spectrum the Taylor terms "
                "cannot follow."
            )
            band = (getattr(self, "science_reffreq", None) or {}).get(
                "valid_hz"
            )
            nterms = int(self.config.get("nterms", 2) or 2)
            if pb is not None and pb < 0.9 and band and nterms < 3:
                low, high = float(band[0]), float(band[1])
                if 0 < low < high:
                    # A beam whose width scales as 1/frequency multiplies
                    # the spectrum by pb**((nu/nu0)**2): curvature
                    # beta = 2 ln(pb) in log-log, which nterms=2 cannot hold.
                    # Its worst misfit over the band is (2/3)|beta|h^2.
                    half = 0.5 * math.log(high / low)
                    misfit = (2.0 / 3.0) * abs(2.0 * math.log(pb)) * half ** 2
                    record["beam_spectral_misfit"] = misfit
                    explanation += (
                        f" At primary beam {pb:.2f} across "
                        f"{low / 1e6:.0f}-{high / 1e6:.0f} MHz, the beam alone "
                        "bends the apparent spectrum enough that nterms="
                        f"{nterms} misfits it by up to {misfit * 100:.0f}% of "
                        "the source's apparent flux; nterms=3 models that "
                        "curvature."
                    )

        diagnostics["residual_peak"] = record
        print(f"Largest residual in {label}: {where}, {explanation}")
        return record

    def check_pixel_sampling(self, label):
        """Warn when the pixels are too coarse for the restoring beam."""

        geometry = self.direction_geometry(self.image_tt0)
        if not geometry:
            return None
        beam = geometry["header"].get("restoringbeam")
        if not isinstance(beam, dict) or "minor" not in beam:
            return None
        try:
            minor = self.angle_in_arcsec(
                beam["minor"]["value"], beam["minor"]["unit"]
            )
        except (KeyError, TypeError, ValueError):
            return None
        cell = max(geometry["cell_arcsec"])
        if not minor or not cell:
            return None

        pixels = minor / cell
        minimum = float(IMAGING_TASK_CONFIG.get("min_pixels_per_beam", 3.0))
        diagnostics = self.clean_diagnostics.setdefault(label, {})
        diagnostics["pixels_per_beam_minor"] = pixels
        print(
            f"Pixel sampling of {label}: {pixels:.1f} pixels across the "
            f"restoring beam's {minor:.2f}\" minor axis ({cell:.2f}\" cells)."
        )
        if pixels < minimum:
            print(
                f"WARNING: {pixels:.1f} pixels per beam is too coarse. A point "
                "source between pixel centres cannot be represented by CLEAN "
                "components, which leaves errors around every bright source "
                "that no amount of cleaning removes. Set `cell` to at most "
                f"{minor / 4.0:.2f}\" (a quarter of the beam) in config.yaml."
            )
        return pixels

    def report_noise_against_thermal(self, label, rms):
        """Say how far the residual noise sits above what the data allow.

        The self-calibration stage measures the per-visibility noise and the
        number of unflagged samples, which together give the thermal noise
        of a naturally weighted image. A residual within a factor of about
        two of it is noise-limited; Briggs weighting and w-projection
        account for that much. Well above it, the image is limited by
        errors in the data (calibration, interference, or sources the model
        does not hold), which no amount of cleaning removes.
        """

        coherence = (
            (getattr(self, "selfcal", None) or {}).get("calibrated_coherence")
            or {}
        )
        thermal = coherence.get("thermal_noise_Jy_beam")
        if not rms or not thermal:
            return None

        ratio = rms / thermal
        diagnostics = self.clean_diagnostics.setdefault(label, {})
        diagnostics["thermal_noise_Jy_beam"] = thermal
        diagnostics["noise_over_thermal"] = ratio
        print(
            f"Noise of {label}: {rms * 1e6:.1f} uJy/beam against a thermal "
            f"floor of {thermal * 1e6:.1f} uJy/beam, {ratio:.1f} times it."
        )
        if ratio > float(
            IMAGING_TASK_CONFIG.get("noise_over_thermal_warn", 3.0)
        ):
            evidence = (
                (getattr(self, "selfcal", None) or {}).get(
                    "phase_start_evidence"
                ) or {}
            )
            ran = bool(
                ((getattr(self, "selfcal", None) or {}).get(
                    "selfcal_strategy"
                ) or {}).get("phase_selfcal")
            )
            print(
                f"WARNING: {label} is limited by errors in the data, not by "
                "noise: artefacts around bright sources are what that looks "
                "like, and cleaning cannot remove them. "
                + (
                    "Phase self-calibration ran; compare the loops in the "
                    "self-calibration report."
                    if ran else
                    "Phase self-calibration did not run"
                    + (
                        f" (expected solution S/N "
                        f"{evidence['expected_solution_snr']:.1f})"
                        if evidence.get("expected_solution_snr") is not None
                        else ""
                    )
                    + ", and it is the correction for exactly this."
                )
            )
        return ratio

    def verify_deconvolution(self, label):
        """Refuse to publish a dirty map as if it were a final image.

        An undeconvolved image is not a faint version of the right answer: it
        is the sky convolved with the array's point-spread function, so every
        source carries sidelobes at the level of the other sources and no
        flux measured from it means anything. If the field really holds
        nothing to deconvolve that is a result and the image stands; if it
        holds a significant peak and the model is empty, every clean path
        this stage has failed and the product must not be presented as final.
        """

        for check in (self.check_pixel_sampling, self.describe_residual_peak):
            try:
                check(label)
            except Exception as exc:
                print(f"WARNING: {check.__name__} could not run: {exc}")

        diagnostics = self.clean_diagnostics.get(label) or {}
        peak_snr = diagnostics.get("residual_peak_snr")
        flux = diagnostics.get("model_flux_Jy")
        rms = diagnostics.get("residual_rms_Jy_beam")

        if diagnostics.get("model_nonzero"):
            print(
                f"Deconvolution of {label}: "
                + (f"{flux:.4f} Jy in the CLEAN model, " if flux else "")
                + (f"residual RMS {rms:.3e} Jy/beam, " if rms else "")
                + (f"residual peak/RMS {peak_snr:.1f}." if peak_snr else "")
            )
            self.report_noise_against_thermal(label, rms)
            return

        minimum = significance_floor(
            IMAGING_TASK_CONFIG.get("undeconvolved_max_peak_snr", 5.0),
            diagnostics,
        )
        if peak_snr is None or peak_snr < minimum:
            print(
                f"WARNING: nothing was deconvolved for {label}. The residual "
                "holds no feature above "
                f"{minimum:.1f} sigma"
                + (f" (peak/RMS {peak_snr:.1f})" if peak_snr else "")
                + ", so the field is consistent with empty and the image is "
                "published as the dirty map it is. Treat any apparent "
                "structure in it as the point-spread function."
            )
            return

        raise RuntimeError(
            f"Refusing to publish {label} as a final image: CLEAN produced "
            "no model, so the image is the dirty map. Its residual holds a "
            f"{peak_snr:.1f}-sigma peak"
            + (f" of {diagnostics['residual_peak_abs_Jy_beam']:.3e} Jy/beam"
               if diagnostics.get("residual_peak_abs_Jy_beam") else "")
            + ", so there is real emission that every clean pass failed to "
            "reach; every source in that image carries the full "
            "point-spread function and no flux measured from it is "
            "meaningful. Check the CASA log for what auto-multithresh chose "
            "as its mask threshold, and lower `automask_noisethreshold` or "
            "raise `clean_bootstrap_niter` in config.yaml. Set "
            "`undeconvolved_max_peak_snr` above "
            f"{peak_snr:.1f} to publish it anyway."
        )

    # --------------------------------------------------------
    # Primary-beam correction
    # --------------------------------------------------------

    def archive_incomplete_casa_product(self, product, label):
        """Move an uncommitted pipeline CASA table aside before rebuilding."""

        product = Path(product).expanduser().resolve()
        if not product.exists():
            return None

        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        archived = product.with_name(
            f"{product.name}.incomplete.{stamp}"
        )
        counter = 1
        while archived.exists():
            archived = product.with_name(
                f"{product.name}.incomplete.{stamp}.{counter}"
            )
            counter += 1
        product.rename(archived)
        if hasattr(self, "stale_pbcor_products"):
            self.stale_pbcor_products.append(str(archived))
        print(
            f"Preserved incomplete {label} before a safe rebuild:\n"
            f"{archived}"
        )
        return str(archived)

    # --------------------------------------------------------

    def archive_pbcor_workdirectory(self, imagename):
        """Move an incomplete widebandpbcor workspace out of retry's way."""

        image_root = Path(imagename).expanduser().resolve()
        workdirectory = image_root.with_name(
            image_root.name + ".pbcor.workdirectory"
        )

        if not workdirectory.exists():
            return None

        stamp = datetime.now(timezone.utc).strftime(
            "%Y%m%dT%H%M%SZ"
        )
        archived = workdirectory.with_name(
            f"{workdirectory.name}.stale.{stamp}"
        )
        counter = 1

        while archived.exists():
            archived = workdirectory.with_name(
                f"{workdirectory.name}.stale.{stamp}.{counter}"
            )
            counter += 1

        workdirectory.rename(archived)
        self.stale_pbcor_workdirectories.append(str(archived))
        print(
            "Preserved incomplete widebandpbcor workspace as:\n"
            f"{archived}"
        )
        return str(archived)

    # --------------------------------------------------------

    def archive_widebandpbcor_outputs(self, imagename):
        """Preserve partial widebandpbcor outputs before retry/fallback."""

        image_root = Path(imagename).expanduser().resolve()
        archived_outputs = []

        for suffix in [
            ".pbcor.image.tt0",
            ".pbcor.image.tt1",
            ".pbcor.alpha",
            ".pbcor.alpha.error",
        ]:
            product = Path(str(image_root) + suffix)

            if not product.exists():
                continue

            stamp = datetime.now(timezone.utc).strftime(
                "%Y%m%dT%H%M%SZ"
            )
            archived = product.with_name(
                f"{product.name}.stale.{stamp}"
            )
            counter = 1

            while archived.exists():
                archived = product.with_name(
                    f"{product.name}.stale.{stamp}.{counter}"
                )
                counter += 1

            product.rename(archived)
            self.stale_pbcor_products.append(str(archived))
            archived_outputs.append(str(archived))
            print(
                "Preserved incomplete widebandpbcor product as:\n"
                f"{archived}"
            )

        return archived_outputs

    # --------------------------------------------------------

    def wideband_pbcor_preflight(self, imagename, nterms):
        """Validate every known widebandpbcor prerequisite without calling it."""

        image_root = Path(imagename).expanduser().resolve()
        vis = Path(self.ms).expanduser().resolve()
        required = []

        for term in range(int(nterms)):
            required.extend([
                Path(f"{image_root}.image.tt{term}"),
                Path(f"{image_root}.residual.tt{term}"),
            ])

        missing = [
            str(path)
            for path in required
            if not casa_image_exists(path)
        ]

        if missing:
            return {
                "ready": False,
                "reason": (
                    "required Taylor image products are missing or invalid"
                ),
                "missing": missing,
            }

        if not vis.is_dir():
            return {
                "ready": False,
                "reason": f"Measurement Set is unavailable: {vis}",
            }

        if not os.access(str(image_root.parent), os.W_OK):
            return {
                "ready": False,
                "reason": (
                    "image directory is not writable: "
                    f"{image_root.parent}"
                ),
            }

        try:
            samples = wideband_pb_samples(
                vis,
                nterms,
                edge_guard=self.wideband_pb_sample_edge_guard,
            )
        except Exception as exc:
            return {
                "ready": False,
                "reason": f"invalid PB frequency samples: {exc}",
            }

        spwlist, chanlist, weightlist = samples

        if not (
            len(spwlist) ==
            len(chanlist) ==
            len(weightlist) and
            len(spwlist) >= int(nterms)
        ):
            return {
                "ready": False,
                "reason": (
                    "SPW/channel/weight sample lists are inconsistent or "
                    "shorter than nterms"
                ),
            }

        return {
            "ready": True,
            "reason": None,
            "spwlist": spwlist,
            "chanlist": chanlist,
            "weightlist": weightlist,
            "edge_guard": self.wideband_pb_sample_edge_guard,
        }

    # --------------------------------------------------------

    def call_widebandpbcor(self, imagename, nterms, preflight=None):
        """Call widebandpbcor with relative image and MS names in one cwd."""

        image_root = Path(imagename).expanduser().resolve()
        vis = Path(self.ms).expanduser().resolve()
        preflight = (
            self.wideband_pbcor_preflight(image_root, nterms)
            if preflight is None
            else dict(preflight)
        )

        if not preflight.get("ready"):
            raise RuntimeError(
                "widebandpbcor preflight rejected the call: "
                f"{preflight.get('reason', 'unknown reason')}"
            )

        spwlist = list(preflight["spwlist"])
        chanlist = list(preflight["chanlist"])
        weightlist = list(preflight["weightlist"])
        original_directory = Path.cwd()
        staged_vis = image_root.parent / (
            f".{image_root.name}.widebandpbcor-vis.{os.getpid()}.ms"
        )
        staged_vis_created = False

        self.archive_pbcor_workdirectory(image_root)
        self.archive_widebandpbcor_outputs(image_root)

        print(
            "widebandpbcor samples: "
            f"spw={spwlist}, channel={chanlist}"
        )

        try:
            # widebandpbcor is documented to fail when the MS and image are in
            # different directories. A process-unique read-only symlink puts
            # both names in the image directory. Passing only basenames also
            # prevents CASA from concatenating an absolute path underneath
            # its .pbcor.workdirectory.
            if staged_vis.exists() or staged_vis.is_symlink():
                raise RuntimeError(
                    "Refusing to replace an existing widebandpbcor staging "
                    f"path: {staged_vis}"
                )

            os.symlink(
                str(vis),
                str(staged_vis),
                target_is_directory=True,
            )
            staged_vis_created = True
            os.chdir(str(image_root.parent))
            widebandpbcor(
                vis=staged_vis.name,
                imagename=image_root.name,
                nterms=int(nterms),
                action="pbcor",
                pbmin=self.pbmin,
                spwlist=spwlist,
                chanlist=chanlist,
                weightlist=weightlist,
            )

        finally:
            os.chdir(str(original_directory))

            if staged_vis_created and staged_vis.is_symlink():
                staged_vis.unlink()

    # --------------------------------------------------------

    def find_pb_image(self, imagename):

        for candidate in [
            Path(imagename + ".pb.tt0"),
            Path(imagename + ".pb"),
        ]:

            if casa_image_exists(candidate):
                return str(candidate)

        return None

    # --------------------------------------------------------

    def call_impbcor_atomic(self, imagename, pbimage, outfile):
        """Publish impbcor output without discarding the prior valid image."""

        outfile = Path(outfile).expanduser().resolve()
        temporary = outfile.with_name(
            f".{outfile.name}.partial.{os.getpid()}"
        )

        if temporary.exists():
            raise RuntimeError(
                f"Temporary impbcor output already exists: {temporary}"
            )

        impbcor(
            imagename=str(Path(imagename).expanduser().resolve()),
            pbimage=str(Path(pbimage).expanduser().resolve()),
            outfile=str(temporary),
            cutoff=self.pbmin,
            overwrite=True,
        )

        if not casa_image_exists(temporary):
            raise RuntimeError(
                "impbcor completed without a usable CASA image: "
                f"{temporary}"
            )

        archived = None
        if outfile.exists():
            stamp = datetime.now(timezone.utc).strftime(
                "%Y%m%dT%H%M%SZ"
            )
            archived = outfile.with_name(
                f"{outfile.name}.stale.{stamp}"
            )
            counter = 1

            while archived.exists():
                archived = outfile.with_name(
                    f"{outfile.name}.stale.{stamp}.{counter}"
                )
                counter += 1

            # The previous corrected image remains in its canonical location
            # until the replacement has completed and passed validation.
            outfile.rename(archived)

        try:
            os.replace(temporary, outfile)
        except Exception:
            if (
                archived is not None and
                archived.exists() and
                not outfile.exists()
            ):
                archived.rename(outfile)
            raise

        if archived is not None:
            self.stale_pbcor_products.append(str(archived))
            print(
                "Preserved previous PB-corrected image as:\n"
                f"{archived}"
            )

        return str(outfile)

    # --------------------------------------------------------

    def find_pbcor_image(self):

        candidates = [
            Path(
                self.imagename +
                ".image.tt0.pbcor"
            ),
            Path(
                self.imagename +
                ".pbcor.image.tt0"
            ),
        ]

        for candidate in candidates:

            if casa_image_exists(candidate):

                return str(candidate)

        return None

    # --------------------------------------------------------

    def image_geometry_descriptor(self, image_path, include_record=False):
        """Read live CASA geometry used to validate a cached beam image."""

        image_path = Path(image_path).expanduser().resolve()
        if not casa_image_exists(image_path):
            raise RuntimeError(f"CASA image is missing or invalid: {image_path}")

        tool = image_tool()
        coords = None
        quanta = None
        try:
            tool.open(str(image_path))
            shape = [int(value) for value in tool.shape()]
            coords = tool.coordsys()
            direction = coords.findcoordinate("direction")
            pixel_axes = [int(value) for value in direction.get("pixel", [])]
            if len(pixel_axes) != 2:
                raise RuntimeError(
                    "The image does not expose exactly two direction axes."
                )
            increments = np.asarray(
                coords.increment(type="direction", format="n")["numeric"],
                dtype=float,
            )
            reference_pixels = np.asarray(
                coords.referencepixel(type="direction")["numeric"],
                dtype=float,
            )
            direction_values = np.asarray(
                coords.referencevalue(type="direction", format="n")["numeric"],
                dtype=float,
            )
            frequency_values = np.asarray(
                coords.referencevalue(type="spectral", format="n")["numeric"],
                dtype=float,
            )
            if frequency_values.size == 0:
                raise RuntimeError("The image has no spectral reference value.")
            spectral_units = coords.units(type="spectral")
            if isinstance(spectral_units, str):
                spectral_unit = spectral_units
            else:
                spectral_units = list(spectral_units)
                spectral_unit = spectral_units[0] if spectral_units else ""

            quanta = quanta_tool()
            frequency_hz = IMAGING_SAFETY.normalize_frequency_hz(
                frequency_values[0],
                spectral_unit,
                quanta.convert,
            )

            # Keep the original CASA record for image.fromarray(). A separate
            # JSON-safe copy is used only for deterministic provenance hashing.
            coordinate_record = coords.torecord()
            canonical_coordinate_record = to_builtin(coordinate_record)
        finally:
            if quanta is not None:
                try:
                    done = getattr(quanta, "done", None)
                    if callable(done):
                        done()
                except Exception:
                    pass
            if coords is not None:
                try:
                    coords.done()
                except Exception:
                    pass
            try:
                tool.close()
            except Exception:
                pass

        canonical_coordinates = json.dumps(
            canonical_coordinate_record,
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        ).encode("utf-8")
        descriptor = {
            "shape": shape,
            "pixel_axes": pixel_axes,
            "direction_increment_rad": increments.tolist(),
            "direction_reference_pixel": reference_pixels.tolist(),
            "direction_reference_value_rad": direction_values.tolist(),
            "frequency_hz": frequency_hz,
            "coordsys_sha256": hashlib.sha256(
                canonical_coordinates
            ).hexdigest(),
        }
        result = {
            "geometry": descriptor,
            "stat_sha256": IMAGING_SAFETY.path_stat_fingerprint(image_path),
        }
        if include_record:
            result["coordinate_record"] = coordinate_record
        return result

    def pointing_centre_offset(self, source_geometry):
        """Check that the image centre is where the antennas were pointing.

        The primary-beam models are evaluated about the image reference
        pixel.  That is the pointing direction for an ordinary single-field
        image, but a re-phased Measurement Set moves the image centre away
        from the antenna boresight, and the beam would then be applied off
        centre -- biasing every off-axis flux density with nothing in the
        products to show it.
        """

        record = {"status": "not_verifiable", "offset_arcmin": None}

        try:
            tool = table()
            opened = False
            try:
                tool.open(str(Path(self.ms) / "FIELD"))
                opened = True
                columns = set(tool.colnames())
                column = (
                    "REFERENCE_DIR"
                    if "REFERENCE_DIR" in columns
                    else "PHASE_DIR"
                )
                values = np.asarray(tool.getcol(column), dtype=float)
            finally:
                if opened:
                    tool.close()

            # CASA returns (2, npoly, nfields) for a direction column.
            flattened = values.reshape(2, -1)
            pointings = [
                (float(flattened[0, index]), float(flattened[1, index]))
                for index in range(flattened.shape[1])
            ]
            offset = IMAGING_SAFETY.nearest_pointing_offset_arcmin(
                source_geometry["direction_reference_value_rad"],
                pointings,
            )
        except Exception as exc:
            record["reason"] = str(exc)
            return record

        record.update({
            "status": "verified",
            "column": column,
            "pointing_count": len(pointings),
            "offset_arcmin": offset,
        })

        # A fractional primary-beam error of order offset/FWHM appears at the
        # half-power radius. Ten arcseconds is well below that for every
        # supported band, so anything larger deserves a look.
        if offset * 60.0 > 10.0:
            record["status"] = "offset"
            print(
                "WARNING: the image reference direction is "
                f"{offset:.3f} arcmin from the nearest {column} in the "
                "Measurement Set. The primary-beam model is centred on the "
                "image, so an offset of this size biases off-axis flux "
                "densities. Confirm the Measurement Set was not re-phased "
                "away from the antenna pointing."
            )

        return record

    # --------------------------------------------------------

    def science_frequency_range_hz(self):
        """Return the observed range clipped to the science-valid band.

        Some receivers digitize wider than they are specified (uGMRT B4 to
        950 MHz against an 850 MHz limit; MeerKAT UHF and L likewise). The
        flagging stage masks the out-of-range channels, but that does not
        move anything CASA computes from the *selected* frequency span.
        """

        profile = TELESCOPE_PROFILE_MODULE.get_band_profile(
            self.telescope,
            self.band,
        )
        low = profile.get("science_freq_min")
        high = profile.get("science_freq_max")
        if low is None or high is None:
            return None

        observed = []
        for info in (self.inspection.get("spws", {}) or {}).values():
            if not isinstance(info, dict):
                continue
            try:
                observed.append(float(info["fmin_MHz"]) * 1e6)
                observed.append(float(info["fmax_MHz"]) * 1e6)
            except (KeyError, TypeError, ValueError):
                continue

        if not observed:
            return None

        observed_low = min(observed)
        observed_high = max(observed)
        valid_low = max(observed_low, float(low))
        valid_high = min(observed_high, float(high))

        if valid_high <= valid_low:
            return None

        return {
            "observed_hz": [observed_low, observed_high],
            "science_hz": [float(low), float(high)],
            "valid_hz": [valid_low, valid_high],
            "extends_outside_science_range": bool(
                observed_low < float(low) or observed_high > float(high)
            ),
        }

    # --------------------------------------------------------

    def science_reference_frequency(self):
        """Reference frequency at the centre of the usable band.

        With reffreq unset, tclean puts the image reference frequency at the
        middle of the selected span, which includes channels the flagging
        stage removed as outside the receiver's specified range. The measured
        primary beam is then evaluated at a frequency the surviving data do
        not represent. For a uGMRT B4 observation covering the full 550-950
        MHz digitizer output the header frequency is 750 MHz while the usable
        data are centred near 700 MHz, and the resulting primary-beam error
        reaches roughly 6% at 15 arcmin and 18% at the PB=0.3 contour --
        past the 5% radial-trend limit the validation contract sets.
        """

        coverage = self.science_frequency_range_hz()
        if not coverage or not coverage["extends_outside_science_range"]:
            return None

        valid_low, valid_high = coverage["valid_hz"]
        observed_low, observed_high = coverage["observed_hz"]

        return {
            **coverage,
            "reffreq_hz": 0.5 * (valid_low + valid_high),
            "tclean_default_reffreq_hz": 0.5 * (observed_low + observed_high),
        }

    # --------------------------------------------------------

    def pblimit_value(self):
        """Return the tclean pblimit appropriate to the beam actually used.

        pblimit masks the restored image using CASA's internal beam model.
        When the flux correction is done with a measured beam supplied by
        this pipeline, CASA's model must not be what decides the published
        aperture; the impbcor pbmin cutoff is.  When CASA's own beam is the
        correction model, CASA's default remains correct.
        """

        if self.primary_beam_provider in {"ugmrt_polynomial", "katbeam"}:
            return float(
                IMAGING_TASK_CONFIG.get("pblimit_measured_beam", 0.05)
            )

        return float(IMAGING_TASK_CONFIG.get("pblimit_casa_beam", 0.2))

    # --------------------------------------------------------

    def effective_aperture(self, beam_2d, increments, reference_pixels):
        """Report the primary-beam level the image actually reaches.

        ``pbmin`` only limits the published product when the image is large
        enough to contain the ``pbmin`` contour.  When it is not, the field
        is truncated by ``imsize * cell`` at a *higher* primary-beam gain, and
        describing the result as "PB >= pbmin" would overstate the field of
        view.  Sources between the image edge and the main lobe also remain
        undeconvolved, so their sidelobes stay in the map.
        """

        beam = np.asarray(beam_2d, dtype=float)
        nx, ny = beam.shape
        x_arcmin = (
            np.arange(nx, dtype=float) - reference_pixels[0]
        ) * increments[0] * 180.0 / math.pi * 60.0
        y_arcmin = (
            np.arange(ny, dtype=float) - reference_pixels[1]
        ) * increments[1] * 180.0 / math.pi * 60.0
        radius = np.hypot(
            *np.meshgrid(x_arcmin, y_arcmin, indexing="ij")
        )

        finite = np.isfinite(beam)
        inside = finite & (beam >= float(self.pbmin))
        edge_levels = [
            float(beam[0, ny // 2]),
            float(beam[nx - 1, ny // 2]),
            float(beam[nx // 2, 0]),
            float(beam[nx // 2, ny - 1]),
        ]
        corner_levels = [
            float(beam[0, 0]),
            float(beam[0, ny - 1]),
            float(beam[nx - 1, 0]),
            float(beam[nx - 1, ny - 1]),
        ]
        # The pbmin contour is enclosed only when no pixel on the image
        # boundary still sits above the cutoff.
        boundary = np.zeros(beam.shape, dtype=bool)
        boundary[0, :] = boundary[-1, :] = True
        boundary[:, 0] = boundary[:, -1] = True
        encloses = not bool(np.any(inside & boundary))

        report = {
            "pbmin": float(self.pbmin),
            "image_half_width_arcmin": float(
                max(abs(x_arcmin[0]), abs(x_arcmin[-1]))
            ),
            "image_half_height_arcmin": float(
                max(abs(y_arcmin[0]), abs(y_arcmin[-1]))
            ),
            "pb_at_edge_midpoints": edge_levels,
            "pb_at_corners": corner_levels,
            "max_pb_on_image_boundary": (
                float(np.nanmax(beam[boundary])) if np.any(boundary) else None
            ),
            "max_radius_above_pbmin_arcmin": (
                float(np.nanmax(radius[inside])) if np.any(inside) else None
            ),
            "image_encloses_pbmin_contour": encloses,
            "effective_pb_cutoff": (
                float(self.pbmin)
                if encloses
                else float(np.nanmax(beam[boundary]))
            ),
        }

        if not encloses:
            print(
                "WARNING: the configured image does not reach the PB="
                f"{self.pbmin:g} cutoff. The field is truncated by "
                "imsize*cell at a primary-beam gain of about "
                f"{report['effective_pb_cutoff']:.3f}, not by pbmin. "
                "Do not describe this product as PB>=%.2f-limited; sources "
                "between the image edge and the main lobe are also left "
                "undeconvolved." % float(self.pbmin)
            )

        return report

    # --------------------------------------------------------

    def build_ugmrt_primary_beam_image(
        self,
        source_image=None,
        image_root=None,
    ):
        """Create or validate the measured 2023 uGMRT polynomial PB.

        The NCRA prescription is an eighth-order power-beam polynomial in
        ``x = radius_arcmin * frequency_GHz``.  CASA's default PB repository
        documents VLA/ALMA models, so an unknown-telescope Airy disk is not an
        acceptable substitute for flux correction of uGMRT images.
        """

        coefficients = self.ugmrt_pb_coefficients
        if not (
            isinstance(coefficients, (list, tuple)) and
            len(coefficients) == 4
        ):
            raise RuntimeError(
                f"No validated uGMRT polynomial coefficients are configured "
                f"for {self.band}."
            )
        coefficients = [float(value) for value in coefficients]

        source = Path(
            source_image or self.image_tt0
        ).expanduser().resolve()
        source_descriptor = self.image_geometry_descriptor(
            source,
            include_record=True,
        )
        source_geometry = source_descriptor["geometry"]
        frequency_mhz = source_geometry["frequency_hz"] / 1e6
        valid_range = self.ugmrt_pb_valid_frequency_mhz
        if not (
            isinstance(valid_range, (list, tuple)) and
            len(valid_range) == 2 and
            float(valid_range[0]) <= frequency_mhz <= float(valid_range[1])
        ):
            raise RuntimeError(
                "The restored-image reference frequency "
                f"({frequency_mhz:.3f} MHz) lies outside the measured "
                f"uGMRT {self.band} primary-beam model range "
                f"{valid_range}."
            )

        beam_root = str(image_root or self.imagename)
        pb_path = Path(beam_root + ".ugmrt.pb.tt0").resolve()
        provenance_path = Path(str(pb_path) + ".provenance.json")
        expected = {
            "schema_version": 1,
            "provider": "ugmrt_polynomial",
            "model": "NCRA-uGMRT-PB-2023-11-29",
            "band": self.band,
            "coefficients": coefficients,
            "valid_frequency_mhz": [
                float(valid_range[0]),
                float(valid_range[1]),
            ],
            "frequency_mhz": frequency_mhz,
            "source_image": str(source),
            "source_geometry": source_geometry,
            "imsize": int(self.config["imsize"]),
            "cell": str(self.config["cell"]),
        }
        radial_cutoff = IMAGING_SAFETY.first_outward_beam_cutoff(
            coefficients,
            float(self.pbmin),
        )
        expected.update({
            "pbmin": float(self.pbmin),
            "main_lobe_radial_argument_cutoff": radial_cutoff,
            "main_lobe_radius_arcmin": (
                radial_cutoff / (frequency_mhz / 1000.0)
            ),
        })

        # The sidecar is the commit marker.  Older code could publish the
        # CASA beam directory and then fail an over-strict coordinate-record
        # hash comparison before writing provenance.  Preserve that
        # uncommitted artifact and rebuild it from the current restored image.
        if pb_path.exists() and not provenance_path.is_file():
            self.archive_incomplete_casa_product(
                pb_path,
                "uGMRT primary-beam image",
            )

        if casa_image_exists(pb_path):
            stale = self.stale_beam_reason(
                load_json(provenance_path, {}),
                expected,
                source,
                source_geometry,
                pb_path,
            )
            if stale is None:
                return str(pb_path)
            print(f"Rebuilding the uGMRT primary beam: {stale}.")
            self.archive_incomplete_casa_product(
                pb_path,
                "stale uGMRT primary-beam image",
            )
        if pb_path.exists():
            raise RuntimeError(f"Incomplete uGMRT PB image exists: {pb_path}")

        shape = source_geometry["shape"]
        pixel_axes = source_geometry["pixel_axes"]
        increments = np.asarray(
            source_geometry["direction_increment_rad"],
            dtype=float,
        )
        reference_pixels = np.asarray(
            source_geometry["direction_reference_pixel"],
            dtype=float,
        )
        nx = shape[pixel_axes[0]]
        ny = shape[pixel_axes[1]]
        x_arcmin = (
            np.arange(nx, dtype=float) - reference_pixels[0]
        ) * increments[0] * 180.0 / math.pi * 60.0
        y_arcmin = (
            np.arange(ny, dtype=float) - reference_pixels[1]
        ) * increments[1] * 180.0 / math.pi * 60.0
        x_grid, y_grid = np.meshgrid(x_arcmin, y_arcmin, indexing="ij")
        radial_argument = np.hypot(x_grid, y_grid) * frequency_mhz / 1000.0
        beam_2d = IMAGING_SAFETY.ugmrt_primary_beam_response(
            radial_argument,
            coefficients,
        )
        beam_2d = np.where(
            np.isfinite(beam_2d) &
            (beam_2d > 0.0) &
            (radial_argument <= radial_cutoff),
            beam_2d,
            0.0,
        )
        if (
            beam_2d.shape != (nx, ny) or
            not np.isfinite(beam_2d).all() or
            float(np.nanmax(beam_2d)) > 1.05
        ):
            raise RuntimeError("The uGMRT PB polynomial produced invalid pixels.")

        aperture = self.effective_aperture(
            beam_2d,
            increments,
            reference_pixels,
        )
        pointing = self.pointing_centre_offset(source_geometry)

        reshape = [1] * len(shape)
        reshape[pixel_axes[0]] = nx
        reshape[pixel_axes[1]] = ny
        beam_pixels = np.broadcast_to(
            beam_2d.reshape(reshape),
            shape,
        ).copy()

        temporary = pb_path.with_name(
            f".{pb_path.name}.partial.{os.getpid()}"
        )
        if temporary.exists():
            raise RuntimeError(f"Temporary uGMRT PB image exists: {temporary}")
        output = image_tool()
        try:
            output.fromarray(
                outfile=str(temporary),
                pixels=beam_pixels,
                csys=source_descriptor["coordinate_record"],
                overwrite=False,
            )
        finally:
            try:
                output.done()
            except Exception:
                pass

        if not casa_image_exists(temporary):
            raise RuntimeError(
                "uGMRT PB generation completed without a usable CASA image."
            )
        temporary_descriptor = self.image_geometry_descriptor(temporary)
        geometry_reason = IMAGING_SAFETY.image_geometry_mismatch(
            source_geometry,
            temporary_descriptor["geometry"],
        )
        if geometry_reason:
            self.archive_incomplete_casa_product(
                temporary,
                "uGMRT primary-beam temporary image",
            )
            raise RuntimeError(
                "New uGMRT PB geometry does not match the restored image: "
                f"{geometry_reason}."
            )
        os.replace(temporary, pb_path)

        pb_descriptor = self.image_geometry_descriptor(pb_path)
        provenance = dict(expected)
        provenance.update({
            "source_stat_sha256": source_descriptor["stat_sha256"],
            "source_content_identity": (
                IMAGING_SAFETY.casa_table_content_identity(source)
            ),
            "pb_image": str(pb_path),
            "pb_stat_sha256": pb_descriptor["stat_sha256"],
            "pb_content_identity": (
                IMAGING_SAFETY.casa_table_content_identity(pb_path)
            ),
            "pb_geometry": pb_descriptor["geometry"],
            "effective_aperture": aperture,
            "pointing_centre_offset": pointing,
        })
        self.save_json(provenance, provenance_path)
        return str(pb_path)

    def stale_beam_reason(
        self,
        provenance,
        expected,
        source,
        source_geometry,
        pb_path,
    ):
        """Why a cached primary-beam image cannot be reused, or None.

        A beam is computed from the restored image and a fixed model, so a
        stale one is rebuilt rather than treated as an error.
        """

        mismatched = sorted(
            key for key, value in expected.items()
            if provenance.get(key) != value
        )
        if mismatched:
            return "its provenance differs in " + ", ".join(mismatched)
        try:
            matches, detail = IMAGING_SAFETY.artifact_matches_provenance(
                provenance,
                "source",
                source,
            )
            if not matches:
                return f"the restored image no longer matches it ({detail})"
            pb_descriptor = self.image_geometry_descriptor(pb_path)
            geometry_reason = IMAGING_SAFETY.image_geometry_mismatch(
                source_geometry,
                pb_descriptor["geometry"],
            )
            if geometry_reason:
                return (
                    "its geometry does not match the restored image "
                    f"({geometry_reason})"
                )
            matches, detail = IMAGING_SAFETY.artifact_matches_provenance(
                provenance,
                "pb",
                pb_path,
            )
        except (OSError, RuntimeError, TypeError, ValueError) as exc:
            return f"it could not be verified ({exc})"
        if not matches:
            return f"the beam image no longer matches its record ({detail})"
        return None

    def corrected_image_provenance(self, source, beam, corrected, **values):
        """Provenance binding a PB-corrected image to its inputs by content."""

        record = dict(values)
        for prefix, path in (
            ("source", source),
            ("pb", beam),
            ("pbcor", corrected),
        ):
            record[f"{prefix}_image"] = str(path)
            record[f"{prefix}_stat_sha256"] = (
                IMAGING_SAFETY.path_stat_fingerprint(path)
            )
            record[f"{prefix}_content_identity"] = (
                IMAGING_SAFETY.casa_table_content_identity(path)
            )
        return record

    def run_ugmrt_primary_beam_correction(self):
        pb_image = self.build_ugmrt_primary_beam_image()
        self.clear_pbcor_dependent_checkpoints()
        corrected_image = self.call_impbcor_atomic(
            self.image_tt0,
            pb_image,
            self.image_tt0 + ".pbcor",
        )
        source = Path(self.image_tt0).expanduser().resolve()
        beam_path = Path(pb_image).expanduser().resolve()
        corrected_path = Path(corrected_image).expanduser().resolve()
        provenance_path = Path(str(corrected_path) + ".provenance.json")
        provenance = self.corrected_image_provenance(
            source,
            beam_path,
            corrected_path,
            schema_version=1,
            provider="ugmrt_polynomial",
            method="ugmrt_polynomial_impbcor_tt0",
            model="NCRA-uGMRT-PB-2023-11-29",
            band=self.band,
            pbmin=float(self.pbmin),
        )
        self.save_json(provenance, provenance_path)
        self.pbcor_image = self.find_pbcor_image()
        if not self.pbcor_image:
            raise RuntimeError("uGMRT PB correction produced no usable image.")
        self.pbcor_method = "ugmrt_polynomial_impbcor_tt0"
        self.primary_beam_image = pb_image
        self.pbcor_status.update({
            "ugmrt_pb_model": "NCRA-uGMRT-PB-2023-11-29",
            "ugmrt_pb_coefficients": [
                float(value) for value in self.ugmrt_pb_coefficients
            ],
            "ugmrt_pb_valid_frequency_mhz": [
                float(value) for value in self.ugmrt_pb_valid_frequency_mhz
            ],
            "ugmrt_pb_image": pb_image,
            "ugmrt_pb_main_lobe_radial_argument_cutoff": (
                load_json(
                    Path(str(pb_image) + ".provenance.json"),
                    {},
                ).get("main_lobe_radial_argument_cutoff")
            ),
            "pbcor_provenance": str(provenance_path),
            "effective_aperture": load_json(
                Path(str(pb_image) + ".provenance.json"),
                {},
            ).get("effective_aperture"),
            "pointing_centre_offset": load_json(
                Path(str(pb_image) + ".provenance.json"),
                {},
            ).get("pointing_centre_offset"),
            "intensity_taylor_term": 0,
            "spectral_terms_corrected": False,
        })

    def build_katbeam_image(self):
        """Create or validate a measured MeerKAT PB CASA image."""

        if not self.katbeam_model:
            raise RuntimeError(
                f"No katbeam model is configured for {self.band}. "
                "Produce and validate an external MeerKAT beam before PB "
                "correction."
            )

        try:
            import katbeam
            JimBeam = katbeam.JimBeam
        except Exception as exc:
            raise RuntimeError(
                "MeerKAT PB correction requires the katbeam package in the "
                "CASA Python environment. Install katbeam or run inside the "
                "processMeerKAT environment."
            ) from exc

        katbeam_version = getattr(katbeam, "__version__", None)
        if not katbeam_version and importlib_metadata is not None:
            try:
                katbeam_version = importlib_metadata.version("katbeam")
            except Exception:
                katbeam_version = None
        if not katbeam_version:
            raise RuntimeError(
                "The installed katbeam package has no discoverable version. "
                "Publication PB provenance requires a versioned package; "
                "install a packaged katbeam release before imaging."
            )
        katbeam_version = str(katbeam_version)

        source = Path(self.image_tt0).expanduser().resolve()
        source_descriptor = self.image_geometry_descriptor(
            source,
            include_record=True,
        )
        source_geometry = source_descriptor["geometry"]
        frequency_mhz = source_geometry["frequency_hz"] / 1e6
        self.katbeam_frequency_mhz = frequency_mhz

        pb_path = Path(self.imagename + ".katbeam.pb.tt0").resolve()
        provenance_path = Path(str(pb_path) + ".provenance.json")
        expected = {
            "schema_version": 2,
            "provider": "katbeam",
            "katbeam_version": katbeam_version,
            "model": self.katbeam_model,
            "frequency_mhz": frequency_mhz,
            "source_image": str(source),
            "source_geometry": source_geometry,
            "imsize": int(self.config["imsize"]),
            "cell": str(self.config["cell"]),
        }

        if pb_path.exists() and not provenance_path.is_file():
            self.archive_incomplete_casa_product(
                pb_path,
                "katbeam primary-beam image",
            )

        if casa_image_exists(pb_path):
            stale = self.stale_beam_reason(
                load_json(provenance_path, {}),
                expected,
                source,
                source_geometry,
                pb_path,
            )
            if stale is None:
                return str(pb_path)
            print(f"Rebuilding the katbeam primary beam: {stale}.")
            self.archive_incomplete_casa_product(
                pb_path,
                "stale katbeam primary-beam image",
            )
        if pb_path.exists():
            raise RuntimeError(f"Incomplete katbeam image already exists: {pb_path}")

        shape = source_geometry["shape"]
        pixel_axes = source_geometry["pixel_axes"]
        increments = np.asarray(
            source_geometry["direction_increment_rad"],
            dtype=float,
        )
        reference_pixels = np.asarray(
            source_geometry["direction_reference_pixel"],
            dtype=float,
        )
        nx = shape[pixel_axes[0]]
        ny = shape[pixel_axes[1]]
        x_deg = (
            np.arange(nx, dtype=float) - reference_pixels[0]
        ) * increments[0] * 180.0 / math.pi
        y_deg = (
            np.arange(ny, dtype=float) - reference_pixels[1]
        ) * increments[1] * 180.0 / math.pi
        x_grid, y_grid = np.meshgrid(x_deg, y_deg, indexing="ij")
        beam_2d = np.asarray(
            JimBeam(self.katbeam_model).I(x_grid, y_grid, frequency_mhz),
            dtype=float,
        )
        # The uGMRT path requires every pixel to be finite and the peak to be
        # normalized.  katbeam must meet the same bar: `.any()` would accept a
        # plane that is almost entirely NaN, and those NaNs would propagate
        # straight through impbcor into the published intensity image.
        if beam_2d.shape != (nx, ny):
            raise RuntimeError("katbeam returned an invalid primary-beam plane.")
        if not np.isfinite(beam_2d).all():
            raise RuntimeError(
                f"katbeam returned {int(np.count_nonzero(~np.isfinite(beam_2d)))} "
                "non-finite primary-beam pixels; refusing to divide the "
                "restored image by a beam containing NaN or Inf."
            )
        beam_peak = float(np.nanmax(beam_2d))
        if not 0.95 <= beam_peak <= 1.05:
            raise RuntimeError(
                f"katbeam returned a beam whose peak is {beam_peak:.4f} "
                "rather than unity; a non-normalized beam would rescale the "
                "published flux-density scale."
            )

        aperture = self.effective_aperture(
            beam_2d,
            increments,
            reference_pixels,
        )
        pointing = self.pointing_centre_offset(source_geometry)

        reshape = [1] * len(shape)
        reshape[pixel_axes[0]] = nx
        reshape[pixel_axes[1]] = ny
        beam_pixels = np.broadcast_to(
            beam_2d.reshape(reshape),
            shape,
        ).copy()

        temporary = pb_path.with_name(
            f".{pb_path.name}.partial.{os.getpid()}"
        )
        if temporary.exists():
            raise RuntimeError(f"Temporary katbeam image exists: {temporary}")
        output = image_tool()
        try:
            output.fromarray(
                outfile=str(temporary),
                pixels=beam_pixels,
                csys=source_descriptor["coordinate_record"],
                overwrite=False,
            )
        finally:
            try:
                output.done()
            except Exception:
                pass

        if not casa_image_exists(temporary):
            raise RuntimeError(
                "katbeam calculation completed without a usable CASA image."
            )
        temporary_descriptor = self.image_geometry_descriptor(temporary)
        geometry_reason = IMAGING_SAFETY.image_geometry_mismatch(
            source_geometry,
            temporary_descriptor["geometry"],
        )
        if geometry_reason:
            self.archive_incomplete_casa_product(
                temporary,
                "katbeam primary-beam temporary image",
            )
            raise RuntimeError(
                "New katbeam PB geometry does not match the restored image: "
                f"{geometry_reason}."
            )
        os.replace(temporary, pb_path)

        pb_descriptor = self.image_geometry_descriptor(pb_path)
        provenance = dict(expected)
        provenance.update({
            "source_stat_sha256": source_descriptor["stat_sha256"],
            "source_content_identity": (
                IMAGING_SAFETY.casa_table_content_identity(source)
            ),
            "pb_image": str(pb_path),
            "pb_stat_sha256": pb_descriptor["stat_sha256"],
            "pb_content_identity": (
                IMAGING_SAFETY.casa_table_content_identity(pb_path)
            ),
            "pb_geometry": pb_descriptor["geometry"],
            "effective_aperture": aperture,
            "pointing_centre_offset": pointing,
        })
        # The sidecar is the cache commit marker and is written only after
        # the CASA table has been atomically published and re-opened.
        self.save_json(provenance, provenance_path)
        return str(pb_path)

    # --------------------------------------------------------

    def run_katbeam_correction(self):

        pb_image = self.build_katbeam_image()
        # A corrected image depends on the restored pixels and cutoff, not
        # only on PB geometry. Recompute it from the validated reusable PB and
        # invalidate every downstream checkpoint before publication.
        self.clear_pbcor_dependent_checkpoints()
        corrected_image = self.call_impbcor_atomic(
            self.image_tt0,
            pb_image,
            self.image_tt0 + ".pbcor",
        )
        source_path = Path(self.image_tt0).expanduser().resolve()
        beam_path = Path(pb_image).expanduser().resolve()
        corrected_path = Path(corrected_image).expanduser().resolve()
        corrected_provenance_path = Path(
            str(corrected_path) + ".provenance.json"
        )
        corrected_provenance = self.corrected_image_provenance(
            source_path,
            beam_path,
            corrected_path,
            schema_version=1,
            provider="katbeam",
            method="katbeam_impbcor_tt0",
            model=self.katbeam_model,
            pbmin=float(self.pbmin),
        )
        self.save_json(corrected_provenance, corrected_provenance_path)
        self.pbcor_image = self.find_pbcor_image()
        if not self.pbcor_image:
            raise RuntimeError(
                "katbeam PB correction did not create the expected image."
            )
        self.pbcor_method = "katbeam_impbcor_tt0"
        self.primary_beam_image = pb_image
        self.pbcor_status["katbeam_model"] = self.katbeam_model
        self.pbcor_status["katbeam_pb_image"] = pb_image
        self.pbcor_status["pbcor_provenance"] = str(
            corrected_provenance_path
        )
        self.pbcor_status["katbeam_frequency_mhz"] = getattr(
            self,
            "katbeam_frequency_mhz",
            None,
        )
        self.pbcor_status["effective_aperture"] = load_json(
            Path(str(pb_image) + ".provenance.json"),
            {},
        ).get("effective_aperture")
        self.pbcor_status["pointing_centre_offset"] = load_json(
            Path(str(pb_image) + ".provenance.json"),
            {},
        ).get("pointing_centre_offset")
        self.pbcor_status["intensity_taylor_term"] = 0
        self.pbcor_status["spectral_terms_corrected"] = False

    # --------------------------------------------------------

    def run_primary_beam_correction(self):

        checkpoint = "pbcor"

        if self.primary_beam_provider == "ugmrt_polynomial":
            print(
                "\nRunning uGMRT measured-polynomial primary-beam "
                "correction from the current restored image."
            )
            self.run_ugmrt_primary_beam_correction()
            self.mark_checkpoint(checkpoint)
            return

        if self.primary_beam_provider == "katbeam":
            print(
                "\nRunning MeerKAT katbeam primary-beam correction from "
                "the current restored image."
            )
            self.run_katbeam_correction()
            self.mark_checkpoint(checkpoint)
            return

        if self.primary_beam_provider == "none":
            reason = (
                f"No validated primary-beam provider is configured for "
                f"{self.telescope}/{self.band}; retaining the flat-noise "
                "restored image without PB correction."
            )
            self.pbcor_status["unavailable_reason"] = reason
            self.pbcor_method = "not_available"
            if self.require_pbcor:
                raise RuntimeError(reason)
            print(f"WARNING: {reason}")
            self.mark_checkpoint(checkpoint)
            return

        if self.primary_beam_provider != "casa":
            raise RuntimeError(
                f"Primary-beam provider {self.primary_beam_provider!r} is "
                "not implemented by the CASA imaging backend."
            )

        existing = self.find_pbcor_image()

        if existing:

            print(
                "Primary-beam corrected image already exists; reusing it."
            )

            self.pbcor_image = existing
            self.pbcor_method = (
                "widebandpbcor"
                if ".pbcor.image." in existing
                else "impbcor"
            )
            self.pbcor_status["reused_existing_product"] = True

            if not self.checkpoint_exists(checkpoint):
                self.mark_checkpoint(checkpoint)

            return

        print("\nRunning primary-beam correction.")

        wideband_error = None

        if (
            self.wideband_pbcor_enabled and
            self.config["nterms"] > 1 and
            CASA_COMPAT.task_available(widebandpbcor)
        ):

            preflight = self.wideband_pbcor_preflight(
                self.imagename,
                self.config["nterms"],
            )
            self.pbcor_status[
                "standalone_wideband_preflight"
            ] = preflight

            if not preflight.get("ready"):
                print(
                    "[imaging] Standalone wideband PB correction was "
                    "skipped before calling CASA: "
                    f"{preflight.get('reason', 'preflight failed')}."
                )
            else:
                self.pbcor_status[
                    "standalone_wideband_attempted"
                ] = True

            try:
                if preflight.get("ready"):
                    self.call_widebandpbcor(
                        self.imagename,
                        self.config["nterms"],
                        preflight=preflight,
                    )
                    self.pbcor_method = "widebandpbcor"

            except Exception as exc:
                wideband_error = exc
                self.pbcor_status[
                    "standalone_wideband_error"
                ] = str(exc)
                self.archive_pbcor_workdirectory(self.imagename)
                self.archive_widebandpbcor_outputs(self.imagename)
                print(
                    "WARNING: widebandpbcor failed; falling back to tt0 "
                    "impbcor so a valid intensity image is still produced: "
                    f"{exc}"
                )

        elif (
            self.wideband_pbcor_enabled and
            self.config["nterms"] > 1
        ):
            self.pbcor_status[
                "standalone_wideband_preflight"
            ] = {
                "ready": False,
                "reason": (
                    "widebandpbcor is not available in this CASA version"
                ),
            }
            print(
                "WARNING: widebandpbcor is not available in this "
                "CASA version. Falling back to impbcor on tt0."
            )

        elif self.config["nterms"] > 1:
            self.pbcor_status[
                "standalone_wideband_preflight"
            ] = {
                "ready": False,
                "reason": (
                    "disabled by wideband_pbcor_enabled=false because CASA "
                    "documents a cross-directory widebandpbcor failure"
                ),
            }
            print(
                "[imaging] Using atomic TT0 primary-beam correction. "
                "Standalone widebandpbcor is disabled by configuration."
            )

        self.pbcor_image = self.find_pbcor_image()

        if (
            self.pbcor_method == "widebandpbcor" and
            not self.pbcor_image
        ):
            wideband_error = RuntimeError(
                "widebandpbcor returned without a usable PB-corrected "
                "Taylor-0 image"
            )
            self.archive_pbcor_workdirectory(self.imagename)
            self.archive_widebandpbcor_outputs(self.imagename)
            print(
                "WARNING: widebandpbcor returned without its expected "
                "output; falling back to tt0 impbcor."
            )

        if not self.pbcor_image:
            pb_image = self.find_pb_image(self.imagename)

            if pb_image is None:
                message = (
                    "No primary-beam image was found. Cannot produce "
                    "a flux-corrected continuum image."
                )

                if self.require_pbcor:
                    raise RuntimeError(message)

                print(f"WARNING: {message}")
                return

            self.call_impbcor_atomic(
                self.image_tt0,
                pb_image,
                self.image_tt0 + ".pbcor",
            )
            self.pbcor_method = (
                "impbcor_tt0_fallback"
                if wideband_error is not None
                else "impbcor_tt0"
            )

        self.pbcor_image = self.find_pbcor_image()

        if not self.pbcor_image:

            message = (
                "PB correction completed, but the expected "
                "PB-corrected image was not found."
            )

            if self.require_pbcor:

                raise RuntimeError(message)

            print(f"WARNING: {message}")
            return

        if wideband_error is not None:
            print(
                "WARNING: Continuum-intensity PB correction succeeded, but "
                "wideband PB-corrected Taylor/spectral-index products were "
                "not created."
            )

        self.mark_checkpoint(checkpoint)

    # --------------------------------------------------------
    # In-band / sub-band imaging
    # --------------------------------------------------------

    def channel_ranges(self, channels):

        channels = sorted(set(channels))

        if not channels:

            return []

        ranges = []
        start = channels[0]
        previous = channels[0]

        for channel in channels[1:]:

            if channel == previous + 1:

                previous = channel
                continue

            ranges.append(
                (start, previous)
            )
            start = channel
            previous = channel

        ranges.append(
            (start, previous)
        )

        return [
            str(start)
            if start == end
            else f"{start}~{end}"
            for start, end in ranges
        ]

    # --------------------------------------------------------

    def make_subband_plan(self):

        n_subbands = int(
            self.config.get(
                "subbands",
                0
            )
        )

        if n_subbands <= 0:

            return []

        msmd = msmetadata()
        opened = False
        channels = []

        try:

            msmd.open(self.ms)
            opened = True

            for spw in spectral_window_ids(msmd):

                freqs = [
                    float(freq)
                    for freq in msmd.chanfreqs(spw)
                ]

                for chan, freq in enumerate(freqs):

                    channels.append(
                        {
                            "spw": int(spw),
                            "chan": int(chan),
                            "freq_Hz": freq,
                        }
                    )

        finally:

            if opened:

                msmd.close()

        if not channels:

            return []

        # A sub-band chunk built from the full digitized span can land wholly
        # in the range the flagging stage removed. That produces an image with
        # no unflagged data and, on the uGMRT path, a hard failure when its
        # reference frequency falls outside the measured beam model. Restrict
        # the plan to the science-valid band before dividing it up.
        coverage = self.science_frequency_range_hz()

        if coverage:
            valid_low, valid_high = coverage["valid_hz"]
            retained = [
                item
                for item in channels
                if valid_low <= item["freq_Hz"] <= valid_high
            ]
            if len(retained) != len(channels):
                print(
                    f"Sub-band planning uses {len(retained)} of "
                    f"{len(channels)} channels; the remainder lie outside the "
                    f"{valid_low / 1e6:.1f}-{valid_high / 1e6:.1f} MHz "
                    "science-valid range."
                )
            channels = retained

        if not channels:

            print(
                "WARNING: no channel lies inside the science-valid range; "
                "skipping in-band imaging."
            )
            return []

        channels.sort(
            key=lambda item: item["freq_Hz"]
        )

        n_subbands = min(
            n_subbands,
            len(channels)
        )

        plan = []

        for index in range(n_subbands):

            start = int(
                math.floor(
                    index *
                    len(channels) /
                    n_subbands
                )
            )
            end = int(
                math.floor(
                    (index + 1) *
                    len(channels) /
                    n_subbands
                )
            )

            chunk = channels[start:end]

            if not chunk:

                continue

            by_spw = {}

            for item in chunk:

                by_spw.setdefault(
                    item["spw"],
                    []
                ).append(
                    item["chan"]
                )

            spw_parts = []

            for spw in sorted(by_spw):

                ranges = ";".join(
                    self.channel_ranges(
                        by_spw[spw]
                    )
                )

                spw_parts.append(
                    f"{spw}:{ranges}"
                )

            freqs = [
                item["freq_Hz"]
                for item in chunk
            ]

            plan.append(
                {
                    "index": index + 1,
                    "name": f"subband_{index + 1:02d}",
                    "spw": ",".join(spw_parts),
                    "nchan": len(chunk),
                    "freq_min_MHz": min(freqs) / 1e6,
                    "freq_max_MHz": max(freqs) / 1e6,
                    "centre_MHz": (
                        sum(freqs) /
                        len(freqs) /
                        1e6
                    ),
                }
            )

        return plan

    # --------------------------------------------------------

    def find_subband_pb(self, imagename):

        candidates = [
            Path(imagename + ".pb"),
            Path(imagename + ".pb.tt0"),
        ]

        for candidate in candidates:

            if casa_image_exists(candidate):

                return str(candidate)

        return None

    # --------------------------------------------------------

    def find_subband_image(self, imagename):

        candidates = [
            Path(imagename + ".image.tt0"),
            Path(imagename + ".image"),
        ]

        for candidate in candidates:

            if casa_image_exists(candidate):

                return str(candidate)

        return None

    # --------------------------------------------------------

    def find_subband_pbcor(self, imagename):

        candidates = [
            Path(imagename + ".image.tt0.pbcor"),
            Path(imagename + ".pbcor.image.tt0"),
            Path(imagename + ".image.pbcor"),
            Path(imagename + ".pbcor.image"),
        ]

        for candidate in candidates:

            if casa_image_exists(candidate):

                return str(candidate)

        return None

    # --------------------------------------------------------

    def ugmrt_subband_product_valid(self, record):
        """Bind a cached sub-band correction to its image and measured PB."""

        if not isinstance(record, dict):
            return False
        image = record.get("image")
        imagename = record.get("imagename")
        corrected = record.get("pbcor_image")
        if not (
            image and imagename and corrected and
            casa_image_exists(image) and
            casa_image_exists(corrected) and
            record.get("pbcor_method") == "ugmrt_polynomial_impbcor"
        ):
            return False

        beam = self.build_ugmrt_primary_beam_image(
            source_image=image,
            image_root=imagename,
        )
        if Path(str(record.get("pb_image"))).resolve() != Path(beam).resolve():
            return False

        provenance_path = Path(str(corrected) + ".provenance.json")
        provenance = load_json(provenance_path, {})
        expected = {
            "schema_version": 1,
            "provider": "ugmrt_polynomial",
            "method": "ugmrt_polynomial_impbcor",
            "model": "NCRA-uGMRT-PB-2023-11-29",
            "band": self.band,
            "source_image": str(Path(image).expanduser().resolve()),
            "pb_image": str(Path(beam).expanduser().resolve()),
            "pbcor_image": str(Path(corrected).expanduser().resolve()),
            "pbmin": float(self.pbmin),
        }
        if any(provenance.get(key) != value for key, value in expected.items()):
            return False
        try:
            return all(
                IMAGING_SAFETY.artifact_matches_provenance(
                    provenance,
                    prefix,
                    path,
                )[0]
                for prefix, path in (
                    ("source", image),
                    ("pb", beam),
                    ("pbcor", corrected),
                )
            )
        except (OSError, TypeError, ValueError):
            return False

    # --------------------------------------------------------

    def run_subband_imaging(self):

        checkpoint = "subband_imaging"

        if not self.make_subbands:

            self.subband_products = []
            return

        report_file = (
            self.reportdir /
            "subband_imaging_report.json"
        )

        if self.checkpoint_exists(checkpoint) and report_file.exists():

            try:

                report = load_json(
                    report_file,
                    {}
                )

            except (OSError, ValueError, TypeError):

                report = {}

            products = report.get("subbands", [])

            if (
                products and
                all(
                    isinstance(item, dict) and
                    casa_image_exists(item.get("image")) and
                    (
                        not self.require_pbcor or
                        casa_image_exists(item.get("pbcor_image"))
                    )
                    for item in products
                )
            ):

                if self.primary_beam_provider == "ugmrt_polynomial":
                    cached_valid = all(
                        self.ugmrt_subband_product_valid(item)
                        for item in products
                    )
                    if not cached_valid:
                        print(
                            "WARNING: cached uGMRT sub-band PB correction "
                            "is not bound to the measured polynomial and "
                            "live image; rebuilding PB corrections."
                        )
                    else:
                        self.subband_products = products
                        print("Sub-band images already exist.")
                        return
                else:
                    self.subband_products = products
                    print("Sub-band images already exist.")
                    return

            print(
                "WARNING: sub-band checkpoint is stale; rebuilding its "
                "report and missing images."
            )
            (self.checkdir / f"{checkpoint}.done").unlink()

        plan = self.make_subband_plan()

        if not plan:

            print(
                "WARNING: No sub-band plan could be created. "
                "Skipping in-band imaging."
            )
            self.subband_products = []
            return

        print(
            f"\nGenerating {len(plan)} PB-corrected "
            "in-band images."
        )

        products = []

        for item in plan:

            imagename = str(
                self.subbanddir /
                (
                    f"{clean_name(self.target)}_"
                    f"{item['name']}"
                )
            )

            image = self.find_subband_image(
                imagename
            )

            if not image:

                kwargs = {
                    "vis": self.ms,
                    "imagename": imagename,
                    "datacolumn": self.datacolumn,
                    "spw": item["spw"],
                    "specmode": IMAGING_TASK_CONFIG.get(
                        "specmode",
                        "mfs",
                    ),
                    "deconvolver": IMAGING_TASK_CONFIG.get(
                        "deconvolver",
                        "mtmfs",
                    ),
                    "nterms": 1,
                    "scales": self.config["scales"],
                    "gridder": self.gridder,
                    "weighting": IMAGING_TASK_CONFIG.get(
                        "weighting",
                        "briggs",
                    ),
                    "robust": self.config["robust"],
                    "imsize": self.config["imsize"],
                    "cell": self.config["cell"],
                    "niter": self.config.get(
                        "subband_niter",
                        self.config["niter"],
                    ),
                    "nsigma": self.config["nsigma"],
                    "interactive": False,
                    "restart": self.restart,
                    "calcpsf": True,
                    "calcres": True,
                    "savemodel": "none",
                    "parallel": self.parallel,
                    "pblimit": self.pblimit_value(),
                }

                kwargs.update(adaptive_automask_kwargs(self.config))
                kwargs.update(clean_control_kwargs())

                if self.gridder == "wproject":

                    kwargs["wprojplanes"] = (
                        self.config["wprojplanes"]
                    )

                if self.gridder == "awproject":

                    kwargs["wbawp"] = True
                    kwargs["conjbeams"] = True
                    kwargs["computepastep"] = 360.0
                    kwargs["rotatepastep"] = 360.0

                print(
                    f"  tclean {item['name']} "
                    f"({item['centre_MHz']:.1f} MHz)"
                )

                self.run_tclean_with_recovery(
                    kwargs,
                    item["name"],
                )

                image = self.find_subband_image(
                    imagename
                )

            pbcor_image = self.find_subband_pbcor(
                imagename
            )

            if (
                self.primary_beam_provider == "ugmrt_polynomial" and
                image
            ):
                # CASA's generic tclean PB repository is not a validated
                # uGMRT beam model.  Every sub-band therefore uses the same
                # measured NCRA polynomial as the Taylor-0 continuum image,
                # evaluated at that sub-band image's own reference frequency.
                pb_image = self.build_ugmrt_primary_beam_image(
                    source_image=image,
                    image_root=imagename,
                )
                pbcor_image = None
            else:
                pb_image = self.find_subband_pb(
                    imagename
                )
            pbcor_method = (
                "existing"
                if pbcor_image
                else None
            )

            if (
                not pbcor_image and
                pb_image and
                image
            ):

                fallback = image + ".pbcor"

                # A one-term sub-band needs only an ordinary impbcor call.
                # The beam passed above is telescope-scoped: the measured
                # polynomial for uGMRT, or tclean's PB only for a profile that
                # explicitly declares CASA as its validated provider.
                self.archive_pbcor_workdirectory(imagename)

                print(
                    f"  impbcor {Path(fallback).name}"
                )

                fallback = self.call_impbcor_atomic(
                    image,
                    pb_image,
                    fallback,
                )

                if casa_image_exists(fallback):

                    pbcor_image = fallback
                    if self.primary_beam_provider == "ugmrt_polynomial":
                        pbcor_method = "ugmrt_polynomial_impbcor"
                        provenance_path = Path(
                            str(pbcor_image) + ".provenance.json"
                        )
                        self.save_json(
                            self.corrected_image_provenance(
                                Path(image).expanduser().resolve(),
                                Path(pb_image).expanduser().resolve(),
                                Path(pbcor_image).expanduser().resolve(),
                                schema_version=1,
                                provider="ugmrt_polynomial",
                                method="ugmrt_polynomial_impbcor",
                                model="NCRA-uGMRT-PB-2023-11-29",
                                band=self.band,
                                pbmin=float(self.pbmin),
                            ),
                            provenance_path,
                        )
                    else:
                        pbcor_method = "impbcor_tclean_pb"

            if not pbcor_image:
                message = (
                    "No PB-corrected product could be created for "
                    f"{item['name']}."
                )

                if self.require_pbcor:
                    print(f"ERROR: {message}")
                else:
                    print(f"WARNING: {message}")

            clean_diagnostics = self.clean_diagnostics.get(
                item["name"]
            )

            if clean_diagnostics is None:
                clean_diagnostics = self.deconvolution_diagnostics(
                    imagename
                )

            record = dict(item)
            record.update(
                {
                    "imagename": imagename,
                    "image": image,
                    "pb_image": pb_image,
                    "pbcor_image": pbcor_image,
                    "pbcor_method": pbcor_method,
                    "clean_diagnostics": clean_diagnostics,
                }
            )

            products.append(record)

            # Persist after each sub-band. A timeout or node failure can then
            # reuse every completed tclean/PB-correction product.
            self.save_json(
                {
                    "complete": False,
                    "completed_count": len(products),
                    "planned_count": len(plan),
                    "subbands": products,
                },
                report_file,
            )

            if not pbcor_image and self.require_pbcor:
                raise RuntimeError(message)

        self.subband_products = products

        self.save_json(
            {
                "complete": True,
                "completed_count": len(products),
                "planned_count": len(plan),
                "subbands": products
            },
            report_file,
        )

        self.mark_checkpoint(checkpoint)

    # --------------------------------------------------------
    # Common beam products
    # --------------------------------------------------------

    def make_common_beam_image(self):

        checkpoint = "common_beam"
        common_beam = self.config.get("common_beam")

        if not common_beam:

            return

        # Convolve the flat-noise restored image and correct afterwards.
        # Smoothing an already primary-beam-corrected map mixes pixels whose
        # noise has been amplified by different beam gains, so the convolved
        # product no longer has a well-defined noise or weighting.
        input_image = self.image_tt0

        outfile = (
            self.imagename +
            ".commonbeam.image"
        )
        corrected_outfile = outfile + ".pbcor"

        if (
            self.checkpoint_exists(checkpoint) and
            casa_image_exists(outfile)
        ):

            self.common_beam_image = outfile
            if casa_image_exists(corrected_outfile):
                self.common_beam_pbcor_image = corrected_outfile
            return

        stale_output = Path(outfile)
        if stale_output.exists():
            print(
                "Removing stale common-beam image after PB-corrected "
                f"input replacement: {stale_output}"
            )
            if casa_image_exists(stale_output):
                rmtables(tablenames=str(stale_output))
            elif stale_output.is_dir():
                shutil.rmtree(stale_output)
            else:
                stale_output.unlink()

        print(
            f"\nSmoothing final image to {common_beam}."
        )

        imsmooth(
            imagename=input_image,
            outfile=outfile,
            targetres=True,
            major=common_beam,
            minor=common_beam,
            pa="0deg",
        )

        if not casa_image_exists(outfile):

            raise RuntimeError(
                "imsmooth completed without producing a usable common-beam "
                f"image: {outfile}"
            )

        self.common_beam_image = outfile

        pb_image = getattr(self, "primary_beam_image", None)

        if pb_image and casa_image_exists(pb_image):
            print(
                "Primary-beam correcting the common-beam image with the "
                "same measured beam used for the continuum product."
            )
            self.common_beam_pbcor_image = self.call_impbcor_atomic(
                outfile,
                pb_image,
                corrected_outfile,
            )
        elif self.require_pbcor:
            print(
                "WARNING: no measured primary beam is available, so the "
                "common-beam product is published flat-noise (uncorrected)."
            )

        self.mark_checkpoint(checkpoint)

    # --------------------------------------------------------
    # Statistics
    # --------------------------------------------------------

    def image_metrics(self, imagename):

        if not imagename or not Path(imagename).exists():

            return {}

        stats = CASA_COMPAT.peak_and_noise(imstat, imagename)

        rms = float(stats["rms"])
        peak = float(stats["max"])
        minimum = float(stats["min"])

        if rms > 0:

            dynamic_range = peak / rms

        else:

            dynamic_range = None

        try:

            header = imhead(
                imagename=imagename,
                mode="summary",
            )

            beam = (
                header.get("restoringbeam")
                if isinstance(header, dict)
                else None
            )

        except Exception:

            beam = None

        return {
            "image": imagename,
            "peak_Jy_beam": peak,
            "min_Jy_beam": minimum,
            "rms_Jy_beam": rms,
            "dynamic_range": dynamic_range,
            "beam": beam,
        }

    # --------------------------------------------------------

    def measure_images(self):

        checkpoint = "metrics"

        if self.checkpoint_exists(checkpoint):

            try:

                report = load_json(
                    self.reportdir /
                    "imaging_metrics.json",
                    {},
                )

            except (OSError, ValueError, TypeError):

                report = {}

            if (
                isinstance(report.get("restored"), dict) and
                isinstance(report.get("residual"), dict)
            ):

                self.metrics = report

                return

            print(
                "WARNING: metrics checkpoint is stale; measuring images "
                "again."
            )
            (self.checkdir / f"{checkpoint}.done").unlink()

        print("\nMeasuring final image statistics.")

        self.metrics = {
            "restored": self.image_metrics(
                self.image_tt0
            ),
            "residual": self.image_metrics(
                self.imagename + ".residual.tt0"
            ),
            "pbcor": self.image_metrics(
                self.pbcor_image
            ),
            "common_beam": self.image_metrics(
                self.common_beam_image
            ),
            "common_beam_pbcor": self.image_metrics(
                self.common_beam_pbcor_image
            ),
            "rms_reference": {
                "map_sensitivity": "residual",
                "note": (
                    "Quote the residual (or restored) RMS as the map "
                    "sensitivity. RMS and dynamic range measured on a "
                    "primary-beam-corrected image are inflated by beam "
                    "amplification toward the field edge and are not a "
                    "sensitivity."
                ),
            },
            "subbands": [
                self.image_metrics(
                    item.get("pbcor_image") or
                    item.get("image")
                )
                for item in self.subband_products
            ],
        }

        self.save_json(
            self.metrics,
            self.reportdir /
            "imaging_metrics.json",
        )

        self.mark_checkpoint(checkpoint)

    # --------------------------------------------------------
    # QA plots
    # --------------------------------------------------------

    def generate_qa_plots(self, force=False):
        """Render final-image QA PNGs using a memory-bounded CASA read."""

        checkpoint = "qa_plots"
        status_path = self.reportdir / "qa_plot_status.json"

        if (
            self.checkpoint_exists(checkpoint) and
            status_path.exists() and
            not force
        ):

            try:

                self.qa_plot_status = load_json(status_path, {})

            except (OSError, ValueError, TypeError):

                self.qa_plot_status = {}

            if isinstance(self.qa_plot_status.get("plots"), list):

                print("Imaging QA plots already generated.")
                return self.qa_plot_status

            print(
                "WARNING: QA checkpoint is stale; regenerating plots."
            )
            (self.checkdir / f"{checkpoint}.done").unlink()

        status = {"plots": [], "warnings": []}
        image_specs = [
            ("restored_image.png", self.image_tt0, "Restored continuum image"),
            (
                "residual_image.png",
                self.imagename + ".residual.tt0",
                "Final continuum residual image",
            ),
            (
                "pbcor_image.png",
                self.pbcor_image,
                "Primary-beam-corrected continuum image",
            ),
            (
                "common_beam_image.png",
                self.common_beam_image,
                "Common-beam continuum image",
            ),
            (
                "common_beam_pbcor_image.png",
                self.common_beam_pbcor_image,
                "Common-beam PB-corrected continuum image",
            ),
        ]

        for filename, imagename, title in image_specs:
            if not imagename or not Path(imagename).exists():
                continue

            try:
                outfile = QA_PRODUCTS.plot_casa_image(
                    imagename,
                    self.plotdir / filename,
                    title,
                )
            except Exception as exc:
                status["warnings"].append(
                    f"Could not create {filename}: {exc}"
                )
                continue

            if outfile is None:
                status["warnings"].append(
                    f"No usable image plane is available for {filename}."
                )
            else:
                status["plots"].append(str(outfile))

        if not status["plots"] and not status["warnings"]:
            status["warnings"].append(
                "No final imaging product is available for QA rendering."
            )

        QA_PRODUCTS.write_status(status_path, status)
        self.qa_plot_status = status
        self.mark_checkpoint(checkpoint)

        if status["plots"]:
            print(
                "Imaging QA plots written: " + ", ".join(status["plots"]),
                flush=True,
            )
        else:
            print(
                "WARNING: No imaging QA plot was created. See "
                f"{status_path}.",
                flush=True,
            )

        return status

    # --------------------------------------------------------
    # FITS export
    # --------------------------------------------------------

    def fits_name(self, casa_image):

        name = Path(casa_image).name

        return str(
            self.fitsdir /
            f"{name}.fits"
        )

    # --------------------------------------------------------

    def export_image(self, casa_image):

        # Auxiliary files such as katbeam provenance JSON live beside CASA
        # image tables and are part of force/reuse accounting, but must never
        # be passed to exportfits.
        if not casa_image or not casa_image_exists(casa_image):

            return None

        fitsimage = self.fits_name(casa_image)

        print(
            f"  exportfits -> {Path(fitsimage).name}"
        )

        exportfits(
            imagename=casa_image,
            fitsimage=fitsimage,
            overwrite=True,
            dropdeg=True,
        )

        return fitsimage

    # --------------------------------------------------------

    def export_fits_products(self):

        checkpoint = "fits"

        if self.checkpoint_exists(checkpoint):

            self.fits_products = [
                str(path)
                for path in self.fitsdir.glob("*.fits")
            ]

            if self.fits_products:

                print(
                    "FITS products already exported."
                )

                return

            print(
                "WARNING: FITS checkpoint is stale; exporting images again."
            )
            (self.checkdir / f"{checkpoint}.done").unlink()

        print("\nExporting CASA images to FITS.")

        products = []
        seen = set()

        for image in [
            str(path)
            for path in self.product_candidates()
            if path.exists()
        ]:

            if (
                image.endswith(".mask") or
                image.endswith(".sumwt.tt0")
            ):

                continue

            if image in seen:

                continue

            seen.add(image)

            fitsimage = self.export_image(image)

            if fitsimage:

                products.append(fitsimage)

        for item in self.subband_products:

            for key in [
                "image",
                "pbcor_image",
            ]:

                image = item.get(key)

                if (
                    not image or
                    image in seen or
                    not Path(image).exists()
                ):

                    continue

                seen.add(image)
                fitsimage = self.export_image(image)

                if fitsimage:

                    products.append(fitsimage)
                    item[f"{key}_fits"] = fitsimage

        self.fits_products = products

        self.mark_checkpoint(checkpoint)

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    def write_report(self):

        fits_products = getattr(
            self,
            "fits_products",
            [],
        )
        fits_content_identities = {
            str(Path(path).expanduser().resolve()): (
                IMAGING_SAFETY.file_content_identity(path)
            )
            for path in fits_products
        }
        measurement_set_stat_sha256 = (
            IMAGING_SAFETY.path_stat_fingerprint(self.ms)
        )

        report = {
            "created_utc": now_utc(),
            "schema_version": 2,
            "provenance_fingerprint": getattr(
                self,
                "provenance_fingerprint",
                None,
            ),
            "measurement_set": self.ms,
            "measurement_set_stat_sha256": measurement_set_stat_sha256,
            "datacolumn": self.datacolumn,
            "target": self.target,
            "telescope": self.telescope,
            "band": self.band,
            "execution_profile": self.profile,
            "imaging_backend": "casa",
            "qa_plot_status": str(
                self.reportdir / "qa_plot_status.json"
            ),
            "imaging_parameters": {
                "imagename": self.imagename,
                "gridder": self.gridder,
                "cell": self.config["cell"],
                "imsize": self.config["imsize"],
                "weighting": IMAGING_TASK_CONFIG.get(
                    "weighting",
                    "briggs",
                ),
                "robust": self.config["robust"],
                "nterms": self.config["nterms"],
                "niter": self.config["niter"],
                "nsigma": self.config["nsigma"],
                "scales": self.config["scales"],
                "pbmin": self.pbmin,
                "pblimit": self.pblimit_value(),
                "science_reference_frequency": self.science_reffreq,
                "primary_beam_provider": self.primary_beam_provider,
                "katbeam_model": self.katbeam_model,
                "katbeam_frequency_mhz": getattr(
                    self,
                    "katbeam_frequency_mhz",
                    None,
                ),
                "wideband_pbcor_enabled": (
                    self.wideband_pbcor_enabled
                ),
                "wideband_pb_sample_edge_guard": (
                    self.wideband_pb_sample_edge_guard
                ),
                "common_beam": self.config.get(
                    "common_beam"
                ),
                "make_subbands": self.make_subbands,
                "subbands": self.config.get(
                    "subbands",
                    0
                ),
                "subband_niter": self.config.get(
                    "subband_niter",
                    self.config["niter"],
                ),
                "parallel": self.parallel,
                "automask": adaptive_automask_kwargs(
                    self.config
                ),
                "clean_controls": clean_control_kwargs(),
                "flat_imaging_overrides": self.flat_imaging_overrides,
            },
            "products": {
                "restored_image": self.image_tt0,
                "pbcor_image": self.pbcor_image,
                "pbcor_method": self.pbcor_method,
                "pbcor_status": self.pbcor_status,
                "common_beam_image": self.common_beam_image,
                "common_beam_pbcor_image": self.common_beam_pbcor_image,
                "subbands": self.subband_products,
                "casa_images": self.existing_products(),
                "fits": fits_products,
                "fits_content_identities": fits_content_identities,
            },
            "preserved_stale_pbcor_workdirectories": (
                self.stale_pbcor_workdirectories
            ),
            "preserved_stale_pbcor_products": (
                self.stale_pbcor_products
            ),
            "clean_diagnostics": self.clean_diagnostics,
            "clean_passes": getattr(self, "clean_passes", {}),
            "automask_adaptations": getattr(
                self, "automask_adaptations", {}
            ),
            "metrics": getattr(
                self,
                "metrics",
                {},
            ),
            "input_reports": {
                "inspection": str(
                    self.outdir /
                    "inspection" /
                    "inspection_report.json"
                ),
                "calibration": str(
                    self.outdir /
                    "calibration" /
                    "calibration_report.json"
                ),
                "selfcal": str(
                    self.outdir /
                    "selfcal" /
                    "reports" /
                    "selfcal_report.json"
                ),
            },
            "checkpoints": {
                "tclean": self.checkpoint_exists(
                    "tclean"
                ),
                "pbcor": self.checkpoint_exists(
                    "pbcor"
                ),
                "common_beam": self.checkpoint_exists(
                    "common_beam"
                ),
                "subband_imaging": self.checkpoint_exists(
                    "subband_imaging"
                ),
                "metrics": self.checkpoint_exists(
                    "metrics"
                ),
                "qa_plots": self.checkpoint_exists(
                    "qa_plots"
                ),
                "fits": self.checkpoint_exists(
                    "fits"
                ),
                "imaging": self.checkpoint_exists(
                    "imaging"
                ),
            },
        }

        outfile = (
            self.reportdir /
            "imaging_report.json"
        )

        self.save_json(report, outfile)

        print(
            f"\nImaging report written:\n{outfile}"
        )

    # --------------------------------------------------------
    # Driver
    # --------------------------------------------------------

    def run(self):

        print("\n")
        print("=" * 70)
        instrument_label = {
            "ugmrt": "uGMRT",
            "meerkat": "MeerKAT",
        }.get(self.telescope, str(self.telescope).upper())
        print(f"Starting {instrument_label} Final Imaging")
        print("=" * 70)

        self.print_summary()

        self.run_tclean()

        self.run_primary_beam_correction()

        self.run_subband_imaging()

        self.make_common_beam_image()

        self.measure_images()

        self.generate_qa_plots()

        self.export_fits_products()

        self.mark_checkpoint("imaging")

        self.write_report()

        print("\n")
        print("=" * 70)
        print(f"{instrument_label} final imaging completed successfully.")
        print("=" * 70)


# ============================================================
# CLI
# ============================================================

def parse_scales(value):

    if value is None:

        return None

    return [
        int(v.strip())
        for v in value.split(",")
        if v.strip()
    ]


def build_parser():

    parser = argparse.ArgumentParser(
        description=(
            "Final continuum imaging stage for the "
            "uGMRT/MeerKAT CASA pipeline."
        )
    )

    parser.add_argument(
        "ms",
        nargs="?",
        default=None,
        help=(
            "Input Measurement Set. If omitted, imaging.py "
            "uses Products/selfcal/reports/selfcal_report.json."
        ),
    )

    parser.add_argument(
        "--outdir",
        default=RUN_CONFIG.get(
            "outdir",
            "Products",
        ),
        help="Pipeline output directory.",
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
        "--target",
        default=RUN_CONFIG.get("target"),
        help="Target name for product naming.",
    )

    parser.add_argument(
        "--provenance-fingerprint",
        default=None,
        help=argparse.SUPPRESS,
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
            "Execution profile. 'laptop' uses a 2048-pixel image and "
            "disables sub-band imaging for a 36-GB Mac."
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
        "--robust",
        type=float,
        default=None,
        help="Briggs robust value.",
    )

    parser.add_argument(
        "--imsize",
        type=int,
        default=None,
        help="Image size in pixels.",
    )

    parser.add_argument(
        "--cell",
        default=None,
        help="Cell size, for example 0.5arcsec.",
    )

    parser.add_argument(
        "--niter",
        type=int,
        default=None,
        help="Maximum clean iterations.",
    )

    parser.add_argument(
        "--nterms",
        type=int,
        default=None,
        help="MT-MFS Taylor terms.",
    )

    parser.add_argument(
        "--nsigma",
        type=float,
        default=None,
        help="Auto-multithresh stopping level.",
    )

    parser.add_argument(
        "--scales",
        type=parse_scales,
        default=None,
        help="Comma-separated multiscale pixels.",
    )

    parser.add_argument(
        "--pbmin",
        type=float,
        default=IMAGING_TASK_CONFIG.get(
            "pbmin",
            0.2,
        ),
        help="Primary-beam correction cutoff.",
    )

    parser.add_argument(
        "--subbands",
        type=int,
        default=None,
        help=(
            "Number of in-band frequency chunks to image. "
            "Defaults are band dependent."
        ),
    )

    parser.add_argument(
        "--subband-niter",
        type=int,
        default=None,
        help="Maximum clean iterations for each sub-band image.",
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
            "Opt in to PB-corrected in-band sub-band images. Off by default; "
            "uGMRT only."
        ),
    )

    parser.add_argument(
        "--common-beam",
        default=None,
        help=(
            "Common circular beam for smoothed product. "
            "Use 'none' to disable."
        ),
    )

    parser.add_argument(
        "--gridder",
        choices=[
            "standard",
            "wproject",
            "awproject",
        ],
        default=None,
        help="CASA tclean gridder.",
    )

    parser.add_argument(
        "--datacolumn",
        choices=[
            "auto",
            "data",
            "corrected",
        ],
        default="auto",
        help="CASA tclean data column.",
    )

    parser.add_argument(
        "--parallel",
        action="store_const",
        const=True,
        default=None,
        dest="parallel",
        help="Run tclean with parallel=True.",
    )

    parser.add_argument(
        "--no-parallel",
        action="store_const",
        const=False,
        dest="parallel",
        help="Run tclean with parallel=False.",
    )

    parser.add_argument(
        "--no-restart",
        action="store_false",
        dest="restart",
        default=bool(
            IMAGING_TASK_CONFIG.get(
                "restart",
                True,
            )
        ),
        help="Run tclean with restart=False.",
    )

    parser.add_argument(
        "--restart",
        action="store_true",
        dest="restart",
        help="Run tclean with restart=True.",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Clear imaging checkpoints and remove existing "
            "CASA image products before running."
        ),
    )

    parser.add_argument(
        "--plots-only",
        action="store_true",
        help=(
            "Regenerate final-image QA plots from existing CASA images "
            "without running tclean."
        ),
    )

    return parser


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

    args = build_parser().parse_args()

    common_beam = args.common_beam

    if (
        common_beam is not None and
        common_beam.lower() == "none"
    ):

        common_beam = ""

    manager = FinalImagingManager(
        ms=args.ms,
        outdir=args.outdir,
        upstream_outdir=args.upstream_outdir,
        target=args.target,
        robust=args.robust,
        imsize=args.imsize,
        cell=args.cell,
        niter=args.niter,
        nterms=args.nterms,
        nsigma=args.nsigma,
        scales=args.scales,
        pbmin=args.pbmin,
        common_beam=common_beam,
        gridder=args.gridder,
        datacolumn=args.datacolumn,
        parallel=args.parallel,
        restart=args.restart,
        force=args.force,
        make_subbands=(
            None
            if args.no_subbands is None
            else not args.no_subbands
        ),
        n_subbands=args.subbands,
        subband_niter=args.subband_niter,
        telescope=args.telescope,
        receiver_band=args.receiver_band,
        profile=args.profile,
        provenance_fingerprint=args.provenance_fingerprint,
    )

    try:
        if args.plots_only:
            manager.generate_qa_plots(force=True)
        else:
            manager.run()
    except KeyboardInterrupt:
        _interrupted("Imaging")
