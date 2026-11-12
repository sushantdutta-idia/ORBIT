#!/usr/bin/env python3
"""Backend-neutral WSClean continuum imaging with LOFAR/EveryBeam support."""

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Avoid importing the local inspect.py when argparse/dataclasses load the
# standard-library module of the same name.
if sys.path and os.path.abspath(sys.path[0]) == SCRIPT_DIR:
    sys.path.pop(0)

import argparse
import importlib.util
import json
import math
import re
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of failing inside a stage.
PIPELINE_INTERFACE_VERSION = 1

# Beam handling by telescope. LOFAR's beam is applied by WSClean itself
# (EveryBeam). MeerKAT's measured katbeam model is applied to WSClean's output
# afterwards (meerkat_beam.py). uGMRT is never imaged by this module.
LOFAR_BEAM_MODES = ("full", "differential", "none")
MEERKAT_BEAM_MODES = ("katbeam", "none")


def load_local_module(filename, module_name):
    path = Path(SCRIPT_DIR) / filename
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONFIG = load_local_module("config.py", "radio_continuum_pipeline_config")
TELESCOPES = load_local_module(
    "telescope_profiles.py",
    "radio_continuum_telescope_profiles",
)
RUN_CONFIG = getattr(CONFIG, "PIPELINE_RUN_CONFIG", {})
IMAGING_CONFIG = getattr(CONFIG, "FINAL_IMAGING_CONFIG", {})
IMAGING_TASK_CONFIG = getattr(CONFIG, "FINAL_IMAGING_TASK_CONFIG", {})
IMAGING_SAFETY = load_local_module(
    "imaging_safety.py",
    "radio_continuum_imaging_safety",
)

_MEERKAT_BEAM = None


def meerkat_beam_module():
    """Load the MeerKAT beam module only when MeerKAT is being imaged.

    It needs numpy and astropy; a LOFAR run should not.
    """

    global _MEERKAT_BEAM
    if _MEERKAT_BEAM is None:
        _MEERKAT_BEAM = load_local_module(
            "meerkat_beam.py", "radio_continuum_meerkat_beam"
        )
    return _MEERKAT_BEAM


def load_json(path, default=None):
    path = Path(path)
    if not path.exists():
        return default
    with open(path) as handle:
        return json.load(handle)


def atomic_write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with open(temporary, "w") as handle:
        json.dump(value, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def command_string(command):
    return " ".join(shlex.quote(str(part)) for part in command)


def wsclean_scale(value):
    text = str(value).strip().lower()
    replacements = {
        "arcseconds": "asec",
        "arcsecond": "asec",
        "arcsec": "asec",
        "arcminutes": "amin",
        "arcminute": "amin",
        "arcmin": "amin",
    }
    for old, new in replacements.items():
        if text.endswith(old):
            return text[: -len(old)] + new
    return text


def same_path(left, right):
    try:
        return Path(left).expanduser().resolve() == (
            Path(right).expanduser().resolve()
        )
    except (OSError, TypeError, ValueError):
        return False


def ms_declares_column(ms, name):
    """Whether a Measurement Set's MAIN table declares ``name``.

    A CASA table stores its column names as plain strings in ``table.dat``, so
    this needs no CASA and changes nothing. None if the table cannot be read.
    """

    try:
        return name.encode("ascii") in (Path(ms) / "table.dat").read_bytes()
    except OSError:
        return None


def comparable_command(command):
    """A WSClean command without the memory cap, which depends on the job.

    The cap follows the memory a job happened to be given; a finished image
    must not be discarded because the next job had a different allocation.
    """

    cleaned = []
    skip = False
    for part in command:
        if skip:
            skip = False
            continue
        if part in {"-abs-mem", "-mem"}:
            skip = True
            continue
        cleaned.append(part)
    return cleaned


def report_telescope(*reports):
    for report in reports:
        value = TELESCOPES.telescope_name_from_report(report or {})
        if value:
            return value
    return None


def report_band(*reports):
    for report in reports:
        value = (report or {}).get("band")
        if isinstance(value, dict):
            value = value.get("name") or value.get("band")
        if value and str(value).upper() in IMAGING_CONFIG:
            return str(value).upper()
    return None


class WSCleanImagingManager:

    def __init__(
        self,
        ms,
        outdir="Products",
        target=None,
        telescope=None,
        receiver_band=None,
        executable="wsclean",
        threads=None,
        profile="science",
        beam_mode=None,
        datacolumn="CORRECTED_DATA",
        extra_args=None,
        dry_run=False,
        provenance_fingerprint=None,
        upstream_outdir=None,
    ):
        self.ms = str(Path(ms).expanduser().resolve())
        self.outdir = Path(outdir).expanduser().resolve()
        # Reports describing the observation as a whole (inspection,
        # calibration, the split targets) live at the run's root; this
        # target's own products and self-calibration report live in outdir.
        self.upstream_outdir = (
            Path(upstream_outdir).expanduser().resolve()
            if upstream_outdir else self.outdir
        )
        self.target_override = target
        self.requested_telescope = telescope
        self.requested_band = receiver_band
        self.executable = executable
        self.threads = int(threads or max(1, os.cpu_count() or 1))
        self.profile = str(profile or "science").lower()
        self.requested_beam_mode = beam_mode
        self.requested_datacolumn = str(datacolumn or "CORRECTED_DATA")
        self.extra_args = list(extra_args or [])
        self.dry_run = bool(dry_run)
        self.provenance_fingerprint = (
            str(provenance_fingerprint).strip()
            if provenance_fingerprint else None
        )

        self.root = self.outdir / "imaging"
        self.imagedir = self.root / "wsclean"
        self.fitsdir = self.root / "fits"
        self.reportdir = self.root / "reports"
        self.checkdir = self.root / "checkpoints"
        for directory in (
            self.root,
            self.imagedir,
            self.fitsdir,
            self.reportdir,
            self.checkdir,
        ):
            directory.mkdir(parents=True, exist_ok=True)

        self.inspection_path = (
            self.upstream_outdir / "inspection" / "inspection_report.json"
        )
        self.inspection = load_json(self.inspection_path, {})
        self.calibration = load_json(
            self.upstream_outdir / "calibration" / "calibration_report.json",
            {},
        )
        self.split_targets = load_json(
            self.upstream_outdir / "calibration" / "split_targets.json",
            {},
        ) or {}
        self.selfcal = load_json(
            self.outdir / "selfcal" / "reports" / "selfcal_report.json",
            {},
        )

        detected_telescope = report_telescope(
            self.selfcal,
            self.calibration,
            self.inspection,
        )
        if self.requested_telescope not in (None, "", "auto"):
            requested = TELESCOPES.normalize_telescope(self.requested_telescope)
            if detected_telescope and requested != detected_telescope:
                raise RuntimeError(
                    f"Requested telescope {requested} conflicts with the "
                    f"inspection report ({detected_telescope})."
                )
            self.telescope = requested
        else:
            self.telescope = detected_telescope

        if not self.telescope:
            raise RuntimeError(
                "No telescope identity is available. Run inspection first or "
                "pass --telescope explicitly."
            )
        if self.telescope == "ugmrt":
            raise RuntimeError(
                "uGMRT data are imaged with CASA tclean, not WSClean: its "
                "measured primary beam, self-calibration and "
                "direction-dependent handoffs are built on CASA images. Use "
                "--imager casa (or leave --imager auto)."
            )
        if self.telescope not in {"lofar", "meerkat"}:
            raise RuntimeError(
                f"WSClean imaging is not configured for {self.telescope}."
            )

        self.beam_mode = self.resolve_beam_mode(self.requested_beam_mode)

        inspected_ms = self.inspection.get("measurement_set")
        self.inspected_input = bool(inspected_ms) and same_path(
            inspected_ms, self.ms
        )
        self.derived_target = None
        if not self.inspected_input:
            if self.telescope == "lofar":
                if not inspected_ms:
                    raise RuntimeError(
                        "LOFAR WSClean imaging requires an inspection report "
                        "bound to the input MS. Run the inspection stage on "
                        "this exact calibrated target MS first."
                    )
                inspected_path = Path(inspected_ms).expanduser().resolve()
                raise RuntimeError(
                    "The inspection report belongs to a different Measurement "
                    f"Set ({inspected_path}) than the requested WSClean input "
                    f"({self.ms}). Use a matching --outdir or rerun "
                    "inspection."
                )
            # MeerKAT is imaged from a product of this pipeline: the
            # self-calibrated MS, or the split target. Inspection saw the raw
            # observation, so it cannot be bound to this MS, but the MS can
            # be bound to the pipeline's own reports.
            if not self.inspection:
                raise RuntimeError(
                    "No inspection report is available. Run the inspection "
                    "stage, or pass --upstream-outdir for a run that has one."
                )
            self.derived_target = self.pipeline_product_target()
            if self.derived_target is None:
                raise RuntimeError(
                    f"{self.ms} is neither the Measurement Set the "
                    "inspection stage saw nor a product recorded by this "
                    "pipeline (the self-calibrated MS in selfcal_report.json "
                    "or an entry of split_targets.json), so its field and "
                    "data column cannot be trusted. Image a calibrated "
                    "target MS from this run, or inspect this MS first."
                )

        self.datacolumn = self.resolve_datacolumn()

        self.band = report_band(
            self.selfcal,
            self.calibration,
            self.inspection,
        )
        if self.requested_band:
            requested_band = TELESCOPES.normalize_band(
                self.telescope,
                self.requested_band,
            )
            if self.band and requested_band != self.band:
                raise RuntimeError(
                    f"Requested band {requested_band} conflicts with the "
                    f"inspection report ({self.band})."
                )
            self.band = requested_band

        if not self.band or self.band not in IMAGING_CONFIG:
            raise RuntimeError(
                "No configured WSClean imaging band is available. Rerun "
                "inspection with --receiver-band."
            )

        self.band_profile = TELESCOPES.get_band_profile(
            self.telescope, self.band
        )
        if self.telescope == "meerkat" and self.beam_mode == "katbeam":
            if (
                self.band_profile.get("primary_beam_provider") != "katbeam" or
                not self.band_profile.get("katbeam_model")
            ):
                raise RuntimeError(
                    f"No validated MeerKAT primary-beam model exists for "
                    f"{self.band}, so WSClean cannot produce a "
                    "primary-beam-corrected image for it."
                )

        self.config = dict(IMAGING_CONFIG[self.band])
        execution_profile = CONFIG.get_execution_profile(self.profile)
        for key, value in execution_profile.get("imaging", {}).items():
            if key not in {"parallel", "make_subbands"}:
                self.config[key] = value

        self.target = self.select_target()
        self.field_id = self.resolve_field_id()
        clean_target = "".join(
            character if character.isalnum() or character in "-_" else "_"
            for character in self.target
        ).strip("_") or "target"
        self.prefix = self.imagedir / f"{clean_target}_final_continuum"
        self.report_path = self.reportdir / "imaging_report.json"

    # --------------------------------------------------------

    def resolve_beam_mode(self, requested):
        if self.telescope == "lofar":
            if requested in (None, ""):
                raise ValueError(
                    "LOFAR beam state must be explicit: use --beam-mode full "
                    "to trust the MS metadata, --beam-mode differential only "
                    "when a phase-centre beam correction is known to need "
                    "forcing, or --beam-mode none for an uncorrected "
                    "quick-look image."
                )
            mode = str(requested).lower()
            if mode not in LOFAR_BEAM_MODES:
                raise ValueError(
                    "--beam-mode must be full, differential, or none."
                )
            return mode

        mode = "katbeam" if requested in (None, "") else str(requested).lower()
        if mode not in MEERKAT_BEAM_MODES:
            raise ValueError(
                "MeerKAT images are corrected with the measured katbeam "
                "primary beam. --beam-mode full and differential are LOFAR "
                "EveryBeam modes; use katbeam (the default) or none for an "
                "uncorrected quick-look image."
            )
        return mode

    # --------------------------------------------------------

    def pipeline_product_target(self):
        """The target a pipeline-made MS belongs to, or None if not one."""

        final_ms = self.selfcal.get("final_ms")
        if final_ms and same_path(final_ms, self.ms):
            target = self.selfcal.get("target")
            if target:
                return str(target)
        for target, product in (self.split_targets or {}).items():
            if same_path(product, self.ms):
                return str(target)
        return None

    # --------------------------------------------------------

    def resolve_datacolumn(self):
        """The visibility column to image.

        ``auto`` follows the rule the CASA backend uses, so both backends
        read the same data: the self-calibrated MS is a ``split`` product
        whose calibrated visibilities are in DATA; any other MS uses
        CORRECTED_DATA when it has one.
        """

        requested = self.requested_datacolumn
        recorded = {
            str(column).upper()
            for column in (self.inspection.get("data_columns") or [])
        } if self.inspected_input else set()

        if requested.lower() == "auto":
            ms = Path(self.ms)
            if ms.name == "final_selfcal.ms" or ms.parent.name == "selfcal_ms":
                return "DATA"
            if recorded:
                return "CORRECTED_DATA" if "CORRECTED_DATA" in recorded else "DATA"
            present = ms_declares_column(self.ms, "CORRECTED_DATA")
            return "CORRECTED_DATA" if present else "DATA"

        if self.inspected_input:
            if not recorded:
                raise RuntimeError(
                    "The inspection report has no Measurement Set data-column "
                    "inventory. Rerun inspection before WSClean imaging."
                )
            if requested.upper() not in recorded:
                raise RuntimeError(
                    f"Requested WSClean data column {requested!r} is not "
                    "present in the inspected MS. Available columns: "
                    + ", ".join(sorted(recorded))
                )
        elif ms_declares_column(self.ms, requested) is False:
            raise RuntimeError(
                f"Requested WSClean data column {requested!r} is not "
                f"present in {self.ms}."
            )
        return requested

    def select_target(self):
        if self.derived_target is not None:
            if (
                self.target_override and
                str(self.target_override) != self.derived_target
            ):
                raise RuntimeError(
                    f"--target {self.target_override!r} does not match the "
                    "target this Measurement Set belongs to "
                    f"({self.derived_target!r})."
                )
            return self.derived_target
        if self.target_override:
            return str(self.target_override)
        candidates = (
            self.calibration.get("targets") or
            self.inspection.get("target_fields") or
            self.inspection.get("fields", {}).get("targets") or
            []
        )
        candidates = list(dict.fromkeys(str(value) for value in candidates))
        if len(candidates) == 1:
            return candidates[0]
        if len(candidates) > 1:
            raise RuntimeError(
                "Multiple targets are present. Pass --target explicitly: "
                + ", ".join(candidates)
            )
        fields = load_json(
            self.outdir / "inspection" / "fields.json",
            [],
        ) or []
        usable = [
            item for item in fields
            if isinstance(item, dict) and "id" in item and "name" in item
        ]
        if len(usable) == 1:
            return str(usable[0]["name"])
        return Path(self.ms).stem

    def resolve_field_id(self):
        """Resolve the selected target to one MS field, never just a label."""

        if not self.inspected_input:
            # A pipeline-made target MS holds that one target and nothing
            # else. Field ids in the inspection report describe the raw
            # observation and would select the wrong field, or none.
            return None

        fields = load_json(
            self.outdir / "inspection" / "fields.json",
            [],
        ) or []
        usable = [
            item
            for item in fields
            if isinstance(item, dict) and "id" in item and "name" in item
        ]
        matches = [
            item
            for item in usable
            if str(item["name"]) == str(self.target)
        ]
        if len(matches) == 1:
            return int(matches[0]["id"])
        if len(matches) > 1:
            raise RuntimeError(
                f"Target name {self.target!r} maps to multiple MS fields. "
                "Split the target to a single-field MS before WSClean."
            )
        if self.target_override:
            scope = (
                "the multiple MS fields"
                if len(usable) > 1
                else "the inspected MS field"
            )
            raise RuntimeError(
                f"Explicit target {self.target!r} could not be mapped to "
                f"{scope}. Pass an exact field name."
            )
        if len(usable) == 1:
            return int(usable[0]["id"])
        if len(usable) > 1:
            raise RuntimeError(
                f"Target {self.target!r} could not be mapped to one of the "
                "multiple MS fields. Run inspection and pass an exact field "
                "name, or provide a target-only calibrated MS."
            )
        raise RuntimeError(
            "The inspection fields inventory is missing or empty. Rerun "
            "inspection on the WSClean input MS before imaging."
        )

    def memory_cap_gb(self):
        """A cap on WSClean's memory, when the job has a limit below the node.

        WSClean sizes itself from the machine's memory, not from a Slurm
        cgroup, so inside a job with a lower limit it plans for memory it may
        not use and is killed. Told the job's limit it stays within it. With
        no job limit the node is the limit and nothing is added.
        """

        if hasattr(self, "_memory_cap"):
            return self._memory_cap
        self._memory_cap = None
        try:
            fraction = float(
                IMAGING_TASK_CONFIG.get("wsclean_memory_fraction", 0.8)
            )
        except (TypeError, ValueError):
            fraction = 0.0
        if fraction <= 0 or any(
            argument in {"-mem", "-abs-mem"} for argument in self.extra_args
        ):
            return None
        try:
            memory = load_local_module(
                "memory_diagnostics.py", "radio_continuum_memory_diagnostics"
            )
            job = (memory.snapshot() or {}).get("job_limit")
        except Exception:
            return None
        if job and job.get("limit_bytes"):
            self._memory_cap = max(
                1, int(job["limit_bytes"] * min(fraction, 1.0) / 1024 ** 3)
            )
        return self._memory_cap

    def build_command(self):
        config = self.config
        command = [
            self.executable,
            "-name", str(self.prefix),
            "-size", str(config["imsize"]), str(config["imsize"]),
            "-scale", wsclean_scale(config["cell"]),
            "-weight", "briggs", str(config.get("robust", -0.25)),
            "-pol", "I",
            "-data-column", self.datacolumn,
            "-niter", str(config.get("niter", 200000)),
            "-auto-mask", str(
                IMAGING_TASK_CONFIG.get("automask_noisethreshold", 4.25)
            ),
            "-auto-threshold", str(config.get("nsigma", 3.0)),
            "-gain", str(IMAGING_TASK_CONFIG.get("clean_gain", 0.1)),
            "-mgain", "0.8",
            "-j", str(self.threads),
        ]

        cap = self.memory_cap_gb()
        if cap:
            command.extend(["-abs-mem", str(cap)])

        if self.field_id is not None:
            command.extend(["-field", str(self.field_id)])

        scales = config.get("scales") or []
        if scales:
            command.extend(
                ["-multiscale", "-multiscale-scales", ",".join(map(str, scales))]
            )

        channels_out = int(config.get("wsclean_channels_out", 0) or 0)
        nterms = int(config.get("nterms", 1))
        if channels_out > 1:
            command.extend(
                [
                    "-join-channels",
                    "-channels-out", str(channels_out),
                    "-fit-spectral-pol", str(max(1, nterms)),
                ]
            )

        if self.telescope == "lofar" and self.beam_mode != "none":
            command.extend(
                [
                    "-apply-primary-beam",
                    "-primary-beam-limit",
                    str(IMAGING_TASK_CONFIG.get("pbmin", 0.3)),
                ]
            )
            if self.beam_mode == "differential":
                command.append("-use-differential-lofar-beam")

        command.extend(self.extra_args)
        command.append(self.ms)
        return command

    def output_products(self):
        paths = sorted(self.imagedir.glob(f"{self.prefix.name}*.fits"))

        def suffix(path):
            # Only what follows WSClean's prefix says what a file is; a
            # target called "model_galaxy" must not hide its own image.
            return path.name[len(self.prefix.name):].lower()

        def is_pbcor(path):
            name = suffix(path)
            return "image-pb" in name or "image_pb" in name

        def is_restored(path):
            name = suffix(path)
            return (
                "image" in name and
                not is_pbcor(path) and
                all(token not in name for token in (
                    "residual", "model", "dirty", "psf", "beam",
                ))
            )

        restored = [path for path in paths if is_restored(path)]
        pbcor = [path for path in paths if is_pbcor(path)]

        def prefer_mfs(values):
            return next(
                (path for path in values if "mfs" in path.name.lower()),
                values[0] if values else None,
            )

        return {
            "all": [str(path) for path in paths],
            "restored": str(prefer_mfs(restored)) if restored else None,
            "pbcor": str(prefer_mfs(pbcor)) if pbcor else None,
        }

    def meerkat_katbeam_applies(self):
        return self.telescope == "meerkat" and self.beam_mode == "katbeam"

    def beam_correction_settings(self):
        """The MeerKAT beam settings that shape the corrected products."""

        science = (
            self.band_profile.get("science_freq_min"),
            self.band_profile.get("science_freq_max"),
        )
        return {
            "model": self.band_profile.get("katbeam_model"),
            "pbmin": float(IMAGING_TASK_CONFIG.get("pbmin", 0.3)),
            "mfs_beam_mode": str(
                IMAGING_TASK_CONFIG.get("meerkat_mfs_beam_mode",
                                        "band_average")
            ),
            "science_limits_hz": (
                [float(science[0]), float(science[1])]
                if all(value is not None for value in science) else None
            ),
            "masked_ranges_hz": self.known_frequency_masks_hz(),
        }

    def known_frequency_masks_hz(self):
        """The frequency ranges the flagging stage masked, in Hz.

        Each sub-band is beam-corrected at the centre of what survived
        flagging, so this is what the flagging stage recorded applying, not
        the configured defaults. Without the report the masks are unknown
        and none is assumed.
        """

        report = load_json(
            self.upstream_outdir / "flagging" /
            "known_frequency_mask_report.json",
            {},
        )
        ranges = []
        items = report.get("applied_intersections") if isinstance(
            report, dict) else None
        for item in items or []:
            try:
                low, high = sorted(
                    float(value) * 1e6 for value in item["intersection_MHz"]
                )
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(low) and math.isfinite(high) and high > low:
                ranges.append([low, high])
        return sorted(ranges)

    def reusable_report(self):
        report = load_json(self.report_path, {})
        if report.get("imaging_backend") != "wsclean":
            return None
        if report.get("provenance_fingerprint") != self.provenance_fingerprint:
            return None
        if report.get("measurement_set") != self.ms:
            return None
        try:
            current_ms_identity = IMAGING_SAFETY.path_stat_fingerprint(self.ms)
        except (OSError, ValueError):
            return None
        if report.get("measurement_set_stat_sha256") != current_ms_identity:
            return None
        if report.get("telescope") != self.telescope or report.get("band") != self.band:
            return None
        if report.get("datacolumn") != self.datacolumn:
            return None
        if report.get("target") != self.target:
            return None
        if report.get("execution_profile") != self.profile:
            return None
        parameters = report.get("imaging_parameters", {})
        if parameters.get("beam_mode") != self.beam_mode:
            return None
        # The memory cap follows whatever the job was given; it does not
        # change the image, so it must not decide whether one is reused.
        if comparable_command(parameters.get("command") or []) != (
            comparable_command(self.build_command())
        ):
            return None
        products = report.get("products", {})
        restored = products.get("restored_image")
        pbcor = products.get("pbcor_image")
        if not restored or not Path(restored).is_file():
            return None
        if self.beam_mode != "none" and (not pbcor or not Path(pbcor).is_file()):
            return None
        if self.meerkat_katbeam_applies():
            settings = self.beam_correction_settings()
            for key in ("pbmin", "mfs_beam_mode"):
                if parameters.get(key) != settings[key]:
                    return None
            status = products.get("pbcor_status") or {}
            provenance = status.get("pbcor_provenance")
            if (
                status.get("katbeam_model") != settings["model"] or
                not provenance or
                not Path(provenance).is_file()
            ):
                return None
            # A correction made before finding 165 has no recorded weighting
            # and is off-axis biased; one made with other masks describes
            # other data. Either is redone, without imaging again.
            if (
                not status.get("mfs_weighting") or
                status.get("masked_ranges_hz") != settings["masked_ranges_hz"]
            ):
                return None
        valid_fits, _ = IMAGING_SAFETY.validate_file_content_identities(
            products.get("fits", []),
            products.get("fits_content_identities"),
        )
        if not valid_fits:
            return None
        return report

    def preflight(self):
        if not Path(self.ms).is_dir():
            raise FileNotFoundError(f"Measurement Set not found: {self.ms}")
        executable = shutil.which(self.executable)
        if executable is None and not Path(self.executable).is_file():
            raise FileNotFoundError(
                f"WSClean executable not found: {self.executable}. Install "
                "WSClean (with EveryBeam for LOFAR primary-beam correction) "
                "or pass --wsclean /path/to/wsclean."
            )
        if self.meerkat_katbeam_applies():
            # Fail now, in seconds, rather than after the image is made.
            module = meerkat_beam_module()
            module.load_katbeam()
            module._astropy_fits()

    def finished_marker_path(self):
        return self.checkdir / "wsclean_finished.json"

    def record_finished_wsclean(self, products):
        """Remember that WSClean itself has finished, and what it made.

        The MeerKAT beam correction runs after WSClean. If it fails, a retry
        should redo the correction, not hours of imaging. The marker binds the
        raw WSClean outputs to the command, the data and the Measurement Set
        that made them.
        """

        raw = [
            path for path in products["all"]
            if not any(token in Path(path).name for token in (
                "-pb.fits", "-katbeam-pb.fits",
            ))
        ]
        atomic_write_json(self.finished_marker_path(), {
            "schema_version": 1,
            "command": comparable_command(self.build_command()),
            "measurement_set": self.ms,
            "measurement_set_stat_sha256": (
                IMAGING_SAFETY.path_stat_fingerprint(self.ms)
            ),
            "datacolumn": self.datacolumn,
            "identities": {
                str(Path(path).expanduser().resolve()): (
                    IMAGING_SAFETY.file_content_identity(path)
                )
                for path in raw
            },
        })

    def finished_wsclean_products(self):
        """The products of an earlier WSClean run, if provably the same run."""

        marker = load_json(self.finished_marker_path(), {})
        if not isinstance(marker, dict) or marker.get("schema_version") != 1:
            return None
        try:
            current = IMAGING_SAFETY.path_stat_fingerprint(self.ms)
        except (OSError, ValueError):
            return None
        if (
            marker.get("measurement_set") != self.ms or
            marker.get("measurement_set_stat_sha256") != current or
            marker.get("datacolumn") != self.datacolumn or
            marker.get("command") != comparable_command(self.build_command())
        ):
            return None
        identities = marker.get("identities")
        valid, _ = IMAGING_SAFETY.validate_file_content_identities(
            list(identities or {}), identities
        )
        if not valid:
            return None
        products = self.output_products()
        return products if products["restored"] else None

    def apply_meerkat_beam(self, products):
        """Divide the restored images by the measured katbeam primary beam."""

        module = meerkat_beam_module()
        settings = self.beam_correction_settings()
        mfs, subbands = module.classify_outputs(
            [Path(path) for path in products["all"]], self.prefix.name
        )
        if mfs is None:
            raise RuntimeError(
                "WSClean did not produce a restored MFS image to correct."
            )
        print(
            f"Applying the MeerKAT katbeam primary beam ({settings['model']}, "
            f"{settings['mfs_beam_mode']} for the MFS image) to {mfs.name} "
            f"and {len(subbands)} sub-band image(s).",
            flush=True,
        )
        record = module.correct_products(
            mfs,
            subbands,
            model=settings["model"],
            pbmin=settings["pbmin"],
            science_limits_hz=settings["science_limits_hz"],
            identity=IMAGING_SAFETY.file_content_identity,
            mode=settings["mfs_beam_mode"],
            masked_ranges_hz=settings["masked_ranges_hz"],
        )
        provenance_path = module.write_provenance(
            Path(record["pbcor_image"] + ".provenance.json"), record
        )
        record["provenance_path"] = str(provenance_path)
        skipped = [
            item for item in record["subbands"] if item.get("skipped")
        ]
        for item in skipped:
            print(
                f"Sub-band {Path(item['source_image']).name} "
                f"({item['nominal_frequency_mhz']:.1f} MHz) was not "
                f"corrected: {item['skipped']}.",
                flush=True,
            )
        print(
            "MFS beam weighted by "
            + (
                "WSClean's imaging weights"
                if record["mfs_weighting"] == "wsclean_image_weight" else
                record["mfs_weighting"].replace("_", " ")
            )
            + f"; the MFS image describes {record['mfs_effective_frequency_mhz']:.1f} "
            f"MHz (its header says {record['frequency_mhz']:.1f} MHz).",
            flush=True,
        )

        bias = record["single_frequency_flux_bias"]
        if "error" in bias:
            print(
                "Could not estimate the error a single-frequency beam would "
                f"carry ({bias['error']}); the correction itself is "
                "unaffected.",
                flush=True,
            )
            bias = {}
        for level, item in sorted(bias.items(), reverse=True):
            print(
                f"A single-frequency beam at {record['frequency_mhz']:.0f} MHz "
                f"would leave a flat-spectrum source at beam level {level} "
                f"({item['radius_deg']:.2f} deg out) "
                f"{item['flux_bias_of_single_frequency_correction'] * 100:+.1f}% "
                "off in the MFS image; the "
                f"{record['mfs_beam_mode'].replace('_', ' ')} correction "
                "used here "
                + (
                    "removes that."
                    if record["mfs_beam_mode"] == "band_average" else
                    "carries it."
                ),
                flush=True,
            )
        print(
            "Spectral index is not corrected: a source whose spectrum is not "
            "flat is attenuated by a beam weighted differently across the "
            "band. The sub-band PB-corrected images are corrected at their "
            "own frequencies.",
            flush=True,
        )
        return record

    def run(self):
        existing = self.reusable_report()
        if existing and not self.dry_run:
            print(f"Reusing validated WSClean image: {existing['products']['restored_image']}")
            return existing

        command = self.build_command()
        print(command_string(command), flush=True)

        if self.dry_run:
            return {"command": command, "dry_run": True}

        self.preflight()

        products = (
            self.finished_wsclean_products()
            if self.meerkat_katbeam_applies() else None
        )
        if products is not None:
            print(
                "WSClean already finished for this data and command; "
                "redoing only the primary-beam correction.",
                flush=True,
            )
        else:
            completed = subprocess.run(command, check=False)
            if completed.returncode != 0:
                raise RuntimeError(
                    f"WSClean exited with status {completed.returncode}."
                    + (
                        " For LOFAR PB imaging, confirm this build includes "
                        "EveryBeam."
                        if self.telescope == "lofar" else ""
                    )
                )

            products = self.output_products()
            if not products["restored"]:
                raise RuntimeError("WSClean did not produce a restored continuum FITS image.")
            if self.meerkat_katbeam_applies():
                self.record_finished_wsclean(products)

        correction = None
        if self.meerkat_katbeam_applies():
            correction = self.apply_meerkat_beam(products)
            products = self.output_products()
        if self.beam_mode != "none" and not products["pbcor"]:
            raise RuntimeError(
                "WSClean did not produce a PB-corrected image. Confirm "
                "EveryBeam support and the LOFAR beam metadata in the MS."
                if self.telescope == "lofar" else
                "The MeerKAT primary-beam correction produced no "
                "PB-corrected image."
            )

        fits_content_identities = {
            str(Path(path).expanduser().resolve()): (
                IMAGING_SAFETY.file_content_identity(path)
            )
            for path in products["all"]
        }

        if self.telescope == "lofar":
            provider = "wsclean_everybeam" if self.beam_mode != "none" else "none"
            pbcor_method = (
                f"wsclean_everybeam_{self.beam_mode}"
                if self.beam_mode != "none" else "not_requested"
            )
        else:
            provider = "katbeam" if self.beam_mode != "none" else "none"
            pbcor_method = (
                correction["method"] if correction else "not_requested"
            )

        parameters = {
            "cell": self.config["cell"],
            "imsize": self.config["imsize"],
            "robust": self.config.get("robust"),
            "niter": self.config.get("niter"),
            "nterms": self.config.get("nterms"),
            "beam_mode": self.beam_mode,
            "primary_beam_provider": provider,
            "command": command,
        }
        subband_products = []
        pbcor_status = None
        if correction:
            settings = self.beam_correction_settings()
            parameters.update({
                "pbmin": settings["pbmin"],
                "mfs_beam_mode": correction["mfs_beam_mode"],
                "channels_out": int(
                    self.config.get("wsclean_channels_out", 0) or 0
                ),
            })
            pbcor_status = {
                "katbeam_model": correction["model"],
                "katbeam_version": correction["katbeam_version"],
                "katbeam_frequency_mhz": correction["frequency_mhz"],
                "effective_frequencies_mhz": (
                    correction["effective_frequencies_mhz"]
                ),
                "mfs_beam_mode": correction["mfs_beam_mode"],
                "mfs_weighting": correction["mfs_weighting"],
                "mfs_effective_frequency_mhz": (
                    correction["mfs_effective_frequency_mhz"]
                ),
                "masked_ranges_hz": correction["masked_ranges_hz"],
                "spectral_terms_corrected": False,
                "pbmin": correction["pbmin"],
                "katbeam_pb_image": correction["pb_image"],
                "pbcor_provenance": correction["provenance_path"],
                "single_frequency_flux_bias": (
                    correction["single_frequency_flux_bias"]
                ),
            }
            subband_products = [
                {
                    "restored_image": item["source_image"],
                    "pbcor_image": item.get("pbcor_image"),
                    "pb_image": item.get("pb_image"),
                    "frequency_mhz": item["frequency_mhz"],
                    "covered_fraction": item.get("covered_fraction"),
                    "skipped": item.get("skipped"),
                }
                for item in correction["subbands"]
            ]

        report = {
            "schema_version": 2,
            "provenance_fingerprint": self.provenance_fingerprint,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "measurement_set": self.ms,
            "measurement_set_stat_sha256": (
                IMAGING_SAFETY.path_stat_fingerprint(self.ms)
            ),
            "datacolumn": self.datacolumn,
            "target": self.target,
            "field_id": self.field_id,
            "telescope": self.telescope,
            "band": self.band,
            "execution_profile": self.profile,
            "imaging_backend": "wsclean",
            "imaging_parameters": parameters,
            "products": {
                "restored_image": products["restored"],
                "pbcor_image": products["pbcor"],
                "pbcor_method": pbcor_method,
                "pbcor_status": pbcor_status,
                "casa_images": [],
                "fits": products["all"],
                "fits_content_identities": fits_content_identities,
                "subbands": subband_products,
            },
            "quality_scope": (
                "direction-independent LOFAR image; deep wide-field science "
                "may require telescope-specific facet calibration"
                if self.telescope == "lofar"
                else "direction-independent continuum image with the "
                "measured katbeam primary beam applied; spectral index is "
                "not corrected"
            ),
        }
        atomic_write_json(self.report_path, report)
        (self.checkdir / "imaging.done").touch()
        return report


def build_parser():
    parser = argparse.ArgumentParser(
        description="Image a calibrated continuum MS with WSClean/EveryBeam."
    )
    parser.add_argument("ms", help="Calibrated target Measurement Set.")
    parser.add_argument("--outdir", default=RUN_CONFIG.get("outdir", "Products"))
    parser.add_argument(
        "--upstream-outdir",
        default=None,
        help=(
            "Run root holding the inspection and calibration reports, when "
            "--outdir is one target's own directory."
        ),
    )
    parser.add_argument("--target", default=RUN_CONFIG.get("target"))
    parser.add_argument("--provenance-fingerprint", default=None,
                        help=argparse.SUPPRESS)
    parser.add_argument(
        "--telescope",
        choices=TELESCOPES.TELESCOPE_CHOICES,
        default=RUN_CONFIG.get("telescope", "auto"),
    )
    parser.add_argument("--receiver-band", default=RUN_CONFIG.get("receiver_band"))
    parser.add_argument("--wsclean", default=RUN_CONFIG.get("wsclean", "wsclean"))
    parser.add_argument("--threads", type=int, default=RUN_CONFIG.get("threads"))
    parser.add_argument(
        "--profile",
        choices=getattr(CONFIG, "EXECUTION_PROFILE_NAMES", ("science",)),
        default=RUN_CONFIG.get("profile", "science"),
    )
    parser.add_argument(
        "--beam-mode",
        choices=["full", "differential", "none", "katbeam"],
        default=RUN_CONFIG.get("lofar_beam_mode"),
        help=(
            "LOFAR: required beam handling; full trusts the MS metadata, "
            "differential forcibly enables the expert differential-beam "
            "switch, none makes no PB-corrected product. MeerKAT: katbeam "
            "(the default) applies the measured katbeam primary beam to the "
            "result; none makes no PB-corrected product."
        ),
    )
    parser.add_argument(
        "--datacolumn",
        default="CORRECTED_DATA",
        help=(
            "Measurement Set data column read by WSClean, or auto to choose "
            "it the way the CASA backend does."
        ),
    )
    parser.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="Additional WSClean argument; may be repeated.",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser


if __name__ == "__main__":
    args = build_parser().parse_args()
    manager = WSCleanImagingManager(
        args.ms,
        outdir=args.outdir,
        target=args.target,
        telescope=args.telescope,
        receiver_band=args.receiver_band,
        executable=args.wsclean,
        threads=args.threads,
        profile=args.profile,
        beam_mode=args.beam_mode,
        datacolumn=args.datacolumn,
        extra_args=args.extra_arg,
        dry_run=args.dry_run,
        provenance_fingerprint=args.provenance_fingerprint,
        upstream_outdir=args.upstream_outdir,
    )
    manager.run()
