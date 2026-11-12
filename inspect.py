#!/usr/bin/env python3
#============================================
# inspect.py
#
# CASA-native uGMRT, MeerKAT, and LOFAR MS inspection module.
# Compatible with CASA 6+
# Author: Sushant Dutta
# Example Usage :
#   python3 inspect.py <MS_Data>
#============================================

import json
import os
import re
import sys
import importlib.util
import sysconfig
from pathlib import Path

SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

_IMPORTED_AS_STDLIB_INSPECT = (__name__ == "inspect")

if _IMPORTED_AS_STDLIB_INSPECT:

    _stdlib_inspect_path = (
        Path(sysconfig.get_path("stdlib")) /
        "inspect.py"
    )
    _stdlib_spec = importlib.util.spec_from_file_location(
        "_stdlib_inspect",
        str(_stdlib_inspect_path),
    )
    _stdlib_inspect = importlib.util.module_from_spec(
        _stdlib_spec
    )
    _stdlib_spec.loader.exec_module(_stdlib_inspect)

    for _name, _value in _stdlib_inspect.__dict__.items():
        if _name not in {
            "__builtins__",
            "__cached__",
            "__doc__",
            "__file__",
            "__loader__",
            "__name__",
            "__package__",
            "__spec__",
        }:
            globals()[_name] = _value

if (
    not _IMPORTED_AS_STDLIB_INSPECT and
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


if not _IMPORTED_AS_STDLIB_INSPECT:
    from casatools import msmetadata, table

    spectral_window_ids = (
        load_casa_compat_module().spectral_window_ids
    )


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


CONFIG = (
    None
    if _IMPORTED_AS_STDLIB_INSPECT
    else load_config_module()
)
RUN_CONFIG = (
    getattr(CONFIG, "PIPELINE_RUN_CONFIG", {})
    if CONFIG is not None
    else {}
)
INSPECTION_DEFAULTS = (
    getattr(CONFIG, "INSPECTION_CONFIG", {})
    if CONFIG is not None
    else {}
)
QA_PRODUCTS = (
    None
    if _IMPORTED_AS_STDLIB_INSPECT
    else load_qa_products_module()
)
TELESCOPE_PROFILE_MODULE = (
    None
    if _IMPORTED_AS_STDLIB_INSPECT
    else load_telescope_profiles_module()
)


CATALOGUE_STATUS = {}


def load_catalogue(filename):

    paths = [
        Path(__file__).parent / filename,
        Path(filename)
    ]

    for p in paths:

        if p.exists():

            with open(p) as f:
                entries = set(
                    line.strip().upper()
                    for line in f
                    if line.strip()
                )

            CATALOGUE_STATUS[str(filename)] = {
                "status": "loaded",
                "path": str(p.resolve()),
                "entries": len(entries),
            }
            return entries

    # A missing catalogue silently reclassifies every calibrator as a science
    # target, so record it rather than leaving only a console warning.
    CATALOGUE_STATUS[str(filename)] = {
        "status": "missing",
        "path": None,
        "entries": 0,
        "searched": [str(item) for item in paths],
    }
    print(f"WARNING: Cannot find {filename}")
    return set()

class MSInspector:

    def __init__(
        self,
        ms,
        outdir="Products",
        telescope="auto",
        receiver_band=None,
        data_state="raw",
    ):

        self.ms = ms
        self.outdir = Path(outdir)
        self.requested_telescope = telescope or "auto"
        self.requested_band = receiver_band
        self.data_state = str(data_state or "raw").strip().lower()

        if self.data_state not in TELESCOPE_PROFILE_MODULE.DATA_STATE_CHOICES:
            raise ValueError(
                "data_state must be 'raw' or 'calibrated'."
            )

        self.inspect_dir = self.outdir / "inspection"
        self.plotdir = self.inspect_dir / "plots"
        self.inspect_dir.mkdir(parents=True, exist_ok=True)
        self.plotdir.mkdir(parents=True, exist_ok=True)

        self.msmd = msmetadata()
        self.tb = table()

        self.flux_catalogue = load_catalogue(
            INSPECTION_DEFAULTS.get(
                "flux_catalogue",
                "flux-cals.list",
            )
        )
        self.phase_catalogue = load_catalogue(
            INSPECTION_DEFAULTS.get(
                "phase_catalogue",
                "vla-cals.list",
            )
        )

    def normalize_name(self, name):

        name = name.upper().strip()
        name = re.sub(r"^J(?=\d{4}[+-]\d+)", "", name)
        name = re.sub(r"^PHASE[_\-]*", "", name)
        name = re.sub(r"^FLUX[_\-]*", "", name)
        name = re.sub(r"^BP[ _\-]*", "", name)

        return name
    
    def is_flux_calibrator(self, name):

        n = self.normalize_name(name)

        return (
            n in self.flux_catalogue or
            str(name).upper().strip() in self.flux_catalogue
        )
    
    def is_phase_calibrator(self, name):

        n = self.normalize_name(name)

        return (
            n in self.phase_catalogue or
            str(name).upper().strip() in self.phase_catalogue
        )
        

    def open(self):
        self.msmd.open(self.ms)

    def close(self):
        self.msmd.close()

    def observatory_names(self):

        names = []

        try:
            names.extend(self.msmd.observatorynames())
        except Exception:
            pass

        if names:
            return sorted({str(name) for name in names if str(name).strip()})

        opened = False

        try:
            self.tb.open(str(Path(self.ms) / "OBSERVATION"))
            opened = True
            if "TELESCOPE_NAME" in self.tb.colnames():
                values = self.tb.getcol("TELESCOPE_NAME")
                names.extend(str(value) for value in values)
        finally:
            if opened:
                self.tb.close()

        return sorted({str(name) for name in names if str(name).strip()})

    def detect_telescope(self):

        return TELESCOPE_PROFILE_MODULE.resolve_telescope(
            self.requested_telescope,
            observatory_names=self.observatory_names(),
            antenna_names=self.msmd.antennanames(),
        )

    def detect_band(self, telescope=None):

        spws = spectral_window_ids(self.msmd)
        freqs = []

        for spw in spws:
            freqs.extend(self.msmd.chanfreqs(spw))

        if not freqs:
            raise RuntimeError(
                "No channel frequencies were found in the Measurement Set."
            )

        fmin = min(freqs)
        fmax = max(freqs)
        centre = 0.5 * (fmin + fmax)

        telescope = telescope or self.detect_telescope()

        return TELESCOPE_PROFILE_MODULE.detect_band(
            telescope["name"],
            centre,
            fmin_hz=fmin,
            fmax_hz=fmax,
            requested_band=self.requested_band,
        )

    def correlation_info(self):

        code_names = {
            1: "I", 2: "Q", 3: "U", 4: "V",
            5: "RR", 6: "RL", 7: "LR", 8: "LL",
            9: "XX", 10: "XY", 11: "YX", 12: "YY",
        }
        setups = {}
        opened = False

        try:
            self.tb.open(str(Path(self.ms) / "POLARIZATION"))
            opened = True
            values = self.tb.getcol("CORR_TYPE")

            # CASA commonly returns shape (ncorr, nsetup).
            if getattr(values, "ndim", 1) == 1:
                values = [values]
            else:
                values = values.T

            for index, setup in enumerate(values):
                setups[str(index)] = [
                    code_names.get(int(code), f"UNKNOWN_{int(code)}")
                    for code in setup
                ]
        except Exception as exc:
            return {
                "setups": {},
                "basis": "unknown",
                "warning": str(exc),
            }
        finally:
            if opened:
                self.tb.close()

        names = {name for setup in setups.values() for name in setup}
        basis = (
            "linear" if names.intersection({"XX", "XY", "YX", "YY"})
            else "circular" if names.intersection({"RR", "RL", "LR", "LL"})
            else "stokes" if names.intersection({"I", "Q", "U", "V"})
            else "unknown"
        )
        return {"setups": setups, "basis": basis}

    def data_columns(self):

        opened = False
        try:
            self.tb.open(str(self.ms))
            opened = True
            return list(self.tb.colnames())
        finally:
            if opened:
                self.tb.close()

    def fields(self):
        names = self.msmd.fieldnames()
        return [{"id": i, "name": n} for i, n in enumerate(names)]

    def classify_fields(self):

        flux = []
        bandpass = []
        phase = []
        targets = []

        # Preferred method: CASA intents

        try:

            names = self.msmd.fieldnames()

            flux_ids = set()

            try:
                flux_ids.update(
                    self.msmd.fieldsforintent("*FLUX*")
                )
            except:
                pass

            bandpass_ids = set()

            try:
                bandpass_ids.update(
                    self.msmd.fieldsforintent("*BANDPASS*")
                )
            except Exception:
                pass

            phase_ids = set()

            try:
                phase_ids.update(
                    self.msmd.fieldsforintent("*PHASE*")
                )
            except:
                pass

            try:
                phase_ids.update(
                    self.msmd.fieldsforintent("*GAIN*")
                )
            except:
                pass

            # CALIBRATE_AMPLI marks the amplitude calibrator. On its own,
            # or when no field carries CALIBRATE_FLUX, that is the
            # flux-scale calibrator (the ALMA and VLA convention), and
            # treating it as a phase intent once put a flux calibrator into
            # the phase list. But MeerKAT's katdal export tags every gain
            # calibrator CALIBRATE_PHASE,CALIBRATE_AMPLI beside a separate
            # CALIBRATE_FLUX field, and on 1708779077 (2026-11-11) that made
            # J2341-3506 a second flux standard, which the MeerKAT flux
            # policy refuses at setjy. A field with both intents next to a
            # distinct flux field is therefore the gain calibrator.
            ampli_ids = set()
            try:
                ampli_ids.update(
                    self.msmd.fieldsforintent("*AMPLI*")
                )
            except Exception:
                pass
            for field_id in ampli_ids:
                if flux_ids and field_id in phase_ids:
                    continue
                flux_ids.add(field_id)

            # A field carrying a FLUX or BANDPASS intent is not the phase
            # calibrator, provided a distinct phase calibrator remains. Keep
            # the overlap when nothing else is available, because a snapshot
            # observation may legitimately use one field for every role.
            distinct_phase_ids = set(phase_ids) - (
                set(flux_ids) | set(bandpass_ids)
            )
            if distinct_phase_ids:
                phase_ids = distinct_phase_ids

            used_ids = set(flux_ids) | set(bandpass_ids) | set(phase_ids)

            targets = [
                names[i]
                for i in range(len(names))
                if i not in used_ids
            ]

            flux = [names[i] for i in flux_ids]
            bandpass = [names[i] for i in bandpass_ids]
            phase = [names[i] for i in phase_ids]

            if len(flux) or len(bandpass) or len(phase):

                # Intent metadata is often incomplete in uGMRT MSs. Retain
                # the trustworthy intent assignments, then classify only
                # otherwise-unassigned fields with the local catalogues.
                for name in names:

                    if name in flux or name in bandpass or name in phase:

                        continue

                    if self.is_flux_calibrator(name):

                        flux.append(name)

                    elif self.is_phase_calibrator(name):

                        phase.append(name)

                used = set(flux) | set(bandpass) | set(phase)
                targets = [name for name in names if name not in used]

                return {
                    "flux_calibrators": sorted(set(flux)),
                    "bandpass_calibrators": sorted(set(bandpass or flux)),
                    "phase_calibrators": sorted(set(phase)),
                    "targets": sorted(set(targets))
                }

        except Exception:
            pass

        # Fallback: catalogue matching

        # No intent match was usable. Discard the provisional "everything is
        # a target" list created above and rebuild every class from the
        # catalogue results.
        flux = []
        bandpass = []
        phase = []
        targets = []

        for f in self.fields():

            name = f["name"]

            if self.is_flux_calibrator(name):

                flux.append(name)

            elif self.is_phase_calibrator(name):

                phase.append(name)

            else:

                targets.append(name)

        return {
            "flux_calibrators": sorted(set(flux)),
            "bandpass_calibrators": sorted(set(flux)),
            "phase_calibrators": sorted(set(phase)),
            "targets": sorted(set(targets))
        }

    def field_intents(self):
        """Record the raw observing intent of every field.

        Classification is inferred from these plus the local catalogues, so
        the inputs to that inference belong in the report. Without them a
        surprising calibrator assignment cannot be diagnosed after the fact.
        """

        intents = {}

        for field_id, name in enumerate(self.msmd.fieldnames()):
            try:
                values = sorted(
                    str(item)
                    for item in self.msmd.intentsforfield(field_id)
                )
            except Exception:
                values = []
            intents[str(name)] = values

        return intents

    # --------------------------------------------------------

    def scan_time_coverage(self):
        """Summarise when each field was observed, and how it interleaves.

        The calibration stage interpolates the phase calibrator's
        scan-averaged solutions onto the target. That is only interpolation
        while calibrator scans bracket the target in time, and its accuracy
        degrades with the gap between visits. Neither fact was measurable
        from any product, although the validation contract asks for both.
        """

        coverage = {
            "status": "unavailable",
            "fields": {},
            "target_phase_coverage": {},
        }

        try:
            names = list(self.msmd.fieldnames())
            scan_spans = {}

            for scan in self.msmd.scannumbers():
                try:
                    times = [
                        float(value)
                        for value in self.msmd.timesforscan(int(scan))
                    ]
                    fields = [
                        int(value)
                        for value in self.msmd.fieldsforscan(int(scan))
                    ]
                except Exception:
                    continue
                if not times or not fields:
                    continue
                scan_spans[int(scan)] = {
                    "fields": fields,
                    "start": min(times),
                    "end": max(times),
                }

            if not scan_spans:
                coverage["reason"] = "no usable scan metadata"
                return coverage

            by_field = {}
            for scan, span in scan_spans.items():
                for field_id in span["fields"]:
                    entry = by_field.setdefault(
                        field_id,
                        {"scans": [], "start": None, "end": None},
                    )
                    entry["scans"].append(scan)
                    entry["start"] = (
                        span["start"]
                        if entry["start"] is None
                        else min(entry["start"], span["start"])
                    )
                    entry["end"] = (
                        span["end"]
                        if entry["end"] is None
                        else max(entry["end"], span["end"])
                    )

            for field_id, entry in sorted(by_field.items()):
                if field_id >= len(names):
                    continue
                coverage["fields"][str(names[field_id])] = {
                    "field_id": field_id,
                    "scan_count": len(entry["scans"]),
                    "scans": sorted(entry["scans"]),
                    "time_start_s": entry["start"],
                    "time_end_s": entry["end"],
                    "elapsed_minutes": (
                        (entry["end"] - entry["start"]) / 60.0
                    ),
                }

            coverage["scan_spans"] = {
                str(scan): span for scan, span in sorted(scan_spans.items())
            }
            coverage["status"] = "measured"

        except Exception as exc:
            coverage["reason"] = str(exc)

        return coverage

    # --------------------------------------------------------

    def target_phase_coverage(self, coverage, field_info):
        """Check each target is bracketed in time by its phase calibrator."""

        result = {}
        scan_spans = coverage.get("scan_spans", {})
        fields = coverage.get("fields", {})

        if coverage.get("status") != "measured":
            return result

        names = list(self.msmd.fieldnames())
        name_by_id = {index: str(name) for index, name in enumerate(names)}

        def centres_for(field_name):
            centres = []
            for span in scan_spans.values():
                for field_id in span.get("fields", []):
                    if name_by_id.get(int(field_id)) == field_name:
                        centres.append(
                            0.5 * (span["start"] + span["end"])
                        )
                        break
            return sorted(centres)

        for target in field_info.get("targets", []):
            target_span = fields.get(str(target))
            if not target_span:
                continue

            start = target_span["time_start_s"]
            end = target_span["time_end_s"]
            entry = {}

            for calibrator in field_info.get("phase_calibrators", []):
                centres = centres_for(str(calibrator))
                if not centres:
                    continue

                before = [value for value in centres if value <= start]
                after = [value for value in centres if value >= end]

                # Largest interval between consecutive calibrator visits that
                # overlaps the time the target was being observed.
                gaps = [
                    (second - first) / 60.0
                    for first, second in zip(centres, centres[1:])
                    if second >= start and first <= end
                ]

                entry[str(calibrator)] = {
                    "calibrator_scans": len(centres),
                    "brackets_before": bool(before),
                    "brackets_after": bool(after),
                    "bracketed": bool(before and after),
                    "max_visit_gap_minutes": max(gaps) if gaps else None,
                    "leading_extrapolation_minutes": (
                        None if before else (min(centres) - start) / 60.0
                    ),
                    "trailing_extrapolation_minutes": (
                        None if after else (end - max(centres)) / 60.0
                    ),
                }

            if entry:
                result[str(target)] = entry

        return result

    # --------------------------------------------------------

    def antenna_info(self):
        return [
            {"id": antenna_id, "name": name}
            for antenna_id, name in enumerate(self.msmd.antennanames())
        ]

    def spw_info(self):

        info = {}

        for spw in spectral_window_ids(self.msmd):

            freqs = list(self.msmd.chanfreqs(spw))

            if not freqs:

                continue

            info[str(spw)] = {
                "id": int(spw),
                "nchan": len(freqs),
                # Preserve channel order so later hard masks can be emitted
                # as exact CASA channel-index selectors even for descending
                # or non-uniform spectral windows.
                "channel_frequencies_MHz": [
                    float(frequency) / 1e6
                    for frequency in freqs
                ],
                "fmin_MHz": min(freqs) / 1e6,
                "fmax_MHz": max(freqs) / 1e6,
                "centre_MHz": (
                    0.5 *
                    (min(freqs) + max(freqs)) /
                    1e6
                ),
                "chan_width_MHz": (
                    abs(freqs[1] - freqs[0]) / 1e6
                    if len(freqs) > 1
                    else None
                )
            }

        return info

    def recommend_refant(self):

        telescope = self.detect_telescope()["name"]
        preferred = INSPECTION_DEFAULTS.get(
            "preferred_refants_by_telescope",
            {},
        ).get(
            telescope,
            TELESCOPE_PROFILE_MODULE.get_telescope_profile(
                telescope
            )["preferred_refants"],
        )

        ants = self.msmd.antennanames()

        for ant in preferred:
            if ant in ants:
                return ant

        if len(ants) == 0:
            return None

        return ants[0]

    @staticmethod
    def coverage_warnings(coverage):
        """Warn when a target is not bracketed by its phase calibrator."""

        warnings = []
        per_target = (coverage or {}).get("target_phase_coverage", {})

        for target, calibrators in sorted(per_target.items()):
            if not calibrators:
                continue

            # Only complain when no configured calibrator brackets the target.
            if not any(item["bracketed"] for item in calibrators.values()):
                details = ", ".join(
                    f"{name}(before={item['brackets_before']}, "
                    f"after={item['brackets_after']})"
                    for name, item in sorted(calibrators.items())
                )
                warnings.append(
                    f"Target {target} is not bracketed in time by any phase "
                    f"calibrator ({details}); the transferred phase solution "
                    "is an extrapolation at one or both ends"
                )

            for name, item in sorted(calibrators.items()):
                gap = item.get("max_visit_gap_minutes")
                if gap is not None and gap > 30.0:
                    warnings.append(
                        f"Phase calibrator {name} revisits target {target} "
                        f"only every {gap:.1f} min at worst; check that the "
                        "interpolated phase remains valid over that interval"
                    )

        return warnings

    def calibration_warnings(self, field_info, coverage=None):

        warnings = []

        if self.data_state == "calibrated":
            if len(field_info["targets"]) == 0:
                warnings.append("No target field found")
            return warnings

        for filename, status in sorted(CATALOGUE_STATUS.items()):
            if status.get("status") != "loaded":
                warnings.append(
                    f"Calibrator catalogue {filename} was not found; every "
                    "unmatched field is classified as a science target"
                )

        if len(field_info["flux_calibrators"]) == 0:
            warnings.append("No flux calibrator found")

        if len(field_info["flux_calibrators"]) > 1:
            warnings.append(
                f"Multiple flux calibrators found "
                f"({len(field_info['flux_calibrators'])}): "
                + ", ".join(field_info["flux_calibrators"])
            )

        if len(field_info["phase_calibrators"]) == 0:
            warnings.append("No phase calibrator found")

        if len(field_info["phase_calibrators"]) > 1:
            warnings.append(
                f"Multiple phase calibrators found "
                f"({len(field_info['phase_calibrators'])}): "
                + ", ".join(field_info["phase_calibrators"])
            )

        if len(field_info["targets"]) == 0:
            warnings.append("No target field found")

        if len(field_info["targets"]) > 1:
            warnings.append(
                f"Multiple target fields found "
                f"({len(field_info['targets'])}): "
                + ", ".join(field_info["targets"])
            )

        warnings.extend(self.coverage_warnings(coverage))

        return warnings

    def calibration_strategy(self, field_info):

        flux = field_info["flux_calibrators"]
        bandpass = field_info.get("bandpass_calibrators", [])
        phase = field_info["phase_calibrators"]

        primary_flux = (
            flux[0]
            if len(flux)
            else None
        )

        secondary_flux = (
            flux[1]
            if len(flux) > 1
            else None
        )

        phase_cal = (
            phase[0]
            if len(phase)
            else None
        )

        return {
            "primary_flux_calibrator": primary_flux,
            "bandpass_calibrator": (
                bandpass[0] if len(bandpass) == 1 else primary_flux
            ),
            "secondary_flux_calibrator": secondary_flux,
            "phase_calibrator": phase_cal,
            "target_fields": field_info["targets"],
        }

    def summary(self):

        field_info = self.classify_fields()
        telescope = self.detect_telescope()
        band = self.detect_band(telescope)
        profile = TELESCOPE_PROFILE_MODULE.get_telescope_profile(
            telescope["name"]
        )
        refant = self.recommend_refant()
        coverage = self.scan_time_coverage()
        coverage["target_phase_coverage"] = self.target_phase_coverage(
            coverage,
            field_info,
        )

        return {

            "schema_version": 3,

            "measurement_set": self.ms,

            "data_state": self.data_state,

            "telescope": telescope,

            "band": band,

            "instrument": {
                "profile_schema_version": (
                    TELESCOPE_PROFILE_MODULE.PROFILE_SCHEMA_VERSION
                ),
                "telescope": telescope["name"],
                "band": band["name"],
                "calibration_backend": profile[
                    "raw_calibration_backend"
                ],
                "imaging_backend": profile["default_imager"],
                "primary_beam_provider": band.get(
                    "primary_beam_provider"
                ),
            },

            "correlations": self.correlation_info(),

            "data_columns": self.data_columns(),

            # Direct access for later pipeline stages

            "flux_calibrators":
                field_info["flux_calibrators"],

            "bandpass_calibrators":
                field_info.get("bandpass_calibrators", []),

            "phase_calibrators":
                field_info["phase_calibrators"],

            "target_fields":
                field_info["targets"],

            "fields": field_info,

            "field_intents": self.field_intents(),

            "calibrator_catalogues": dict(CATALOGUE_STATUS),

            "time_coverage": coverage,

            "spws": self.spw_info(),

            "calibration_strategy":
                self.calibration_strategy(field_info),

            # QA summary

            "calibration_summary": {

                "nscans":
                    len(self.msmd.scannumbers()),

                "n_flux_calibrators":
                    len(field_info["flux_calibrators"]),

                "n_phase_calibrators":
                    len(field_info["phase_calibrators"]),

                "n_bandpass_calibrators":
                    len(field_info.get("bandpass_calibrators", [])),

                "n_targets":
                    len(field_info["targets"]),
            },

            "warnings":
                self.calibration_warnings(field_info, coverage),

            "reference_antenna":
                {
                    "recommended": refant,
                    "method": "preferred_central_antenna"
                },

            "nantennas":
                len(self.msmd.antennanames()),

            "nfields":
                len(self.msmd.fieldnames()),

            "nspws":
                len(spectral_window_ids(self.msmd))
        }

    def write_products(self):

        summary = self.summary()

        with open(self.inspect_dir / "inspection_report.json", "w") as f:
            json.dump(summary, f, indent=4)

        with open(self.inspect_dir / "observation_summary.txt", "w") as f:
            f.write(json.dumps(summary, indent=4))

        with open(self.inspect_dir / "fields.json", "w") as f:
            json.dump(self.fields(), f, indent=4)

        with open(self.inspect_dir / "antennas.json", "w") as f:
            json.dump(self.antenna_info(), f, indent=4)

        self.generate_qa_plots(summary)

    def generate_qa_plots(self, summary):
        """Write a non-interactive spectral-window coverage diagnostic."""

        status = {"plots": [], "warnings": []}

        try:
            outfile = QA_PRODUCTS.plot_spw_coverage(
                summary.get("spws", {}),
                self.plotdir / "spw_frequency_coverage.png",
            )
            if outfile is None:
                status["warnings"].append(
                    "No usable spectral-window frequency information."
                )
            else:
                status["plots"].append(str(outfile))
        except Exception as exc:
            status["warnings"].append(
                f"Could not create spectral-window coverage plot: {exc}"
            )

        QA_PRODUCTS.write_status(
            self.inspect_dir / "qa_plot_status.json",
            status,
        )

        if status["plots"]:
            print(
                "Inspection QA plot written: " + ", ".join(status["plots"]),
                flush=True,
            )
        else:
            print(
                "WARNING: No inspection QA plot was created. See "
                f"{self.inspect_dir / 'qa_plot_status.json'}.",
                flush=True,
            )

    def run(self):

        self.open()

        try:
            self.write_products()

        finally:
            self.close()

def build_parser():

    parser = argparse.ArgumentParser(
        description="Inspect a uGMRT, MeerKAT, or LOFAR Measurement Set."
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
        "--telescope",
        choices=TELESCOPE_PROFILE_MODULE.TELESCOPE_CHOICES,
        default=RUN_CONFIG.get("telescope", "auto"),
        help=(
            "Telescope override. The default reads "
            "OBSERVATION/TELESCOPE_NAME and verifies antenna names."
        ),
    )

    parser.add_argument(
        "--receiver-band",
        default=RUN_CONFIG.get("receiver_band"),
        help=(
            "Explicit receiver band, required for a MeerKAT sub-band "
            "wholly inside the UHF/L frequency overlap."
        ),
    )

    parser.add_argument(
        "--data-state",
        choices=TELESCOPE_PROFILE_MODULE.DATA_STATE_CHOICES,
        default=RUN_CONFIG.get("data_state", "raw"),
        help="Whether the input is raw or already cross-calibrated.",
    )

    return parser


if __name__ == "__main__":

    args = build_parser().parse_args()

    inspector = MSInspector(
        args.ms,
        outdir=args.outdir,
        telescope=args.telescope,
        receiver_band=args.receiver_band,
        data_state=args.data_state,
    )
    inspector.run()

    print("Inspection complete.")
