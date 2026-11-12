#!/usr/bin/env python3
"""CASA-free telescope identity, receiver-band, and capability registry.

The original pipeline inferred a uGMRT band from frequency alone.  That is
unsafe once MeerKAT and LOFAR are accepted because their receiver ranges
overlap.  This module therefore resolves the telescope first, then resolves a
band inside that telescope's namespace.  It has no CASA imports so launchers,
tests, and batch-file generators can use the same rules as the CASA stages.
"""

from __future__ import annotations

# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of silently keeping an old
# default (this release moved MeerKAT's default imager to WSClean).
PIPELINE_INTERFACE_VERSION = 1

import copy
import re
from typing import Iterable, Mapping, Optional


PROFILE_SCHEMA_VERSION = 4

TELESCOPE_NAMES = ("ugmrt", "meerkat", "lofar")
TELESCOPE_CHOICES = ("auto",) + TELESCOPE_NAMES
DATA_STATE_CHOICES = ("raw", "calibrated")
IMAGER_CHOICES = ("auto", "casa", "wsclean")


TELESCOPE_PROFILES = {
    "ugmrt": {
        "label": "uGMRT",
        "metadata_patterns": (
            r"\bugmrt\b",
            r"\bgmrt\b",
            r"giant\s+metrewave\s+radio\s+telescope",
        ),
        # The array has a central square (C00-C14) plus three arms: east
        # (E02-E06), west (W01-W06) and south (S01-S06). Omitting S dropped
        # six of the thirty antennas from telescope inference, which matters
        # for a sub-array selection dominated by the south arm.
        "antenna_pattern": r"^(?:C\d{2}|[EWS]\d{2})$",
        "preferred_refants": ["C00", "C05", "C08", "C09", "C01", "C02"],
        "raw_calibration_backend": "casa",
        "default_imager": "casa",
        "casa_selfcal_supported": True,
        "raw_pipeline_supported": True,
        "bands": {
            "B2": {
                "aliases": ["2", "BAND2", "BAND-2"],
                "label": "uGMRT Band 2",
                "freq_min": 125e6,
                "freq_max": 250e6,
                "cell": "1.0arcsec",
                "imsize": 8192,
                "primary_beam_provider": "ugmrt_polynomial",
                "primary_beam_coefficients": [
                    -3.089, 39.314, -23.011, 5.037,
                ],
                "primary_beam_valid_frequency_mhz": [125.0, 250.0],
            },
            "B3": {
                "aliases": ["3", "BAND3", "BAND-3"],
                "label": "uGMRT Band 3",
                "freq_min": 250e6,
                "freq_max": 500e6,
                "cell": "0.8arcsec",
                "imsize": 8192,
                "primary_beam_provider": "ugmrt_polynomial",
                "primary_beam_coefficients": [
                    -3.129, 38.816, -21.608, 4.483,
                ],
                "primary_beam_valid_frequency_mhz": [250.0, 500.0],
            },
            "B4": {
                "aliases": ["4", "BAND4", "BAND-4"],
                "label": "uGMRT Band 4",
                "freq_min": 550e6,
                # A 400-MHz GWB setup can materialize channels through
                # 950 MHz even though the Band-4 front end is specified only
                # through 850 MHz.  Keep the digitized range distinct from
                # the unattended science-valid range: inspection may accept
                # the table, but flagging must mask every channel above the
                # receiver limit before calibration.
                "freq_max": 950e6,
                "science_freq_min": 550e6,
                "science_freq_max": 850e6,
                "guaranteed_freq_min": 550e6,
                "guaranteed_freq_max": 850e6,
                "cell": "0.5arcsec",
                "imsize": 8192,
                "primary_beam_provider": "ugmrt_polynomial",
                "primary_beam_coefficients": [
                    -3.263, 42.618, -25.580, 5.823,
                ],
                "primary_beam_valid_frequency_mhz": [550.0, 850.0],
            },
            "B5": {
                "aliases": ["5", "BAND5", "BAND-5"],
                "label": "uGMRT Band 5",
                "freq_min": 1000e6,
                "freq_max": 1460e6,
                "cell": "0.3arcsec",
                "imsize": 8192,
                "primary_beam_provider": "ugmrt_polynomial",
                "primary_beam_coefficients": [
                    -2.614, 27.594, -13.268, 2.395,
                ],
                "primary_beam_valid_frequency_mhz": [1050.0, 1450.0],
            },
        },
    },
    "meerkat": {
        "label": "MeerKAT",
        "metadata_patterns": (r"\bmeerkat\b", r"\bmkat\b"),
        "antenna_pattern": r"^m\d{3}$",
        # The flagging stage replaces this initial choice with the least
        # flagged suitable antenna.  These are merely deterministic fallbacks.
        "preferred_refants": ["m059", "m037", "m031", "m027", "m020", "m001"],
        "raw_calibration_backend": "casa",
        # WSClean images MeerKAT's wide band with joined-channel deconvolution
        # and its measured katbeam primary beam is applied to the result
        # (meerkat_beam.py). CASA tclean stays available with --imager casa.
        "default_imager": "wsclean",
        "casa_selfcal_supported": True,
        "raw_pipeline_supported": True,
        "bands": {
            "MEERKAT_UHF": {
                "aliases": ["UHF", "UHF-BAND", "UHF_BAND"],
                "label": "MeerKAT UHF band",
                "freq_min": 544e6,
                "freq_max": 1088e6,
                # Guaranteed/science-processing range.  The wider limits
                # above describe the digitized receiver output, not the band
                # over which unattended continuum calibration is guaranteed.
                "science_freq_min": 580e6,
                "science_freq_max": 1015e6,
                "guaranteed_freq_min": 580e6,
                "guaranteed_freq_max": 1015e6,
                "cell": "2.0arcsec",
                "imsize": 6144,
                "primary_beam_provider": "katbeam",
                "katbeam_model": "MKAT-AA-UHF-JIM-2020",
            },
            "MEERKAT_L": {
                "aliases": ["L", "LBAND", "L-BAND", "L_BAND"],
                "label": "MeerKAT L band",
                "freq_min": 856e6,
                "freq_max": 1712e6,
                "science_freq_min": 900e6,
                "science_freq_max": 1670e6,
                "guaranteed_freq_min": 900e6,
                "guaranteed_freq_max": 1670e6,
                "cell": "1.5arcsec",
                "imsize": 6144,
                "primary_beam_provider": "katbeam",
                "katbeam_model": "MKAT-AA-L-JIM-2020",
            },
            # The S receiver exposes five overlapping 875-MHz tunings. They
            # are not interchangeable: narrow observations in an overlap
            # require the tuning identity from metadata or an explicit CLI
            # selection.
            "MEERKAT_S0": {
                "aliases": ["S0", "S0-BAND", "S0_BAND", "S-BAND-0"],
                "label": "MeerKAT S0 band",
                "freq_min": 1750e6,
                "freq_max": 2625e6,
                "cell": "0.6arcsec",
                "imsize": 6144,
                # katbeam 0.1 publishes measured models only for UHF and L.
                "primary_beam_provider": "none",
            },
            "MEERKAT_S1": {
                "aliases": ["S1", "S1-BAND", "S1_BAND", "S-BAND-1"],
                "label": "MeerKAT S1 band",
                "freq_min": 1968.75e6,
                "freq_max": 2843.75e6,
                "cell": "0.6arcsec",
                "imsize": 6144,
                "primary_beam_provider": "none",
            },
            "MEERKAT_S2": {
                "aliases": ["S2", "S2-BAND", "S2_BAND", "S-BAND-2"],
                "label": "MeerKAT S2 band",
                "freq_min": 2187.5e6,
                "freq_max": 3062.5e6,
                "cell": "0.6arcsec",
                "imsize": 6144,
                "primary_beam_provider": "none",
            },
            "MEERKAT_S3": {
                "aliases": ["S3", "S3-BAND", "S3_BAND", "S-BAND-3"],
                "label": "MeerKAT S3 band",
                "freq_min": 2406.25e6,
                "freq_max": 3281.25e6,
                "cell": "0.6arcsec",
                "imsize": 6144,
                "primary_beam_provider": "none",
            },
            "MEERKAT_S4": {
                "aliases": ["S4", "S4-BAND", "S4_BAND", "S-BAND-4"],
                "label": "MeerKAT S4 band",
                "freq_min": 2625e6,
                "freq_max": 3500e6,
                "cell": "0.6arcsec",
                "imsize": 6144,
                "primary_beam_provider": "none",
            },
        },
    },
    "lofar": {
        "label": "LOFAR/ILT",
        "metadata_patterns": (
            r"\blofar\b",
            r"international\s+lofar\s+telescope",
            r"\bilt\b",
        ),
        "antenna_pattern": r"^(?:CS|RS|DE|FR|IE|LV|PL|SE|UK)\d{3}",
        "preferred_refants": [
            "CS001HBA0",
            "CS001HBA1",
            "CS002HBA0",
            "CS002HBA1",
            "CS001LBA",
            "ST001",
        ],
        "raw_calibration_backend": "linc",
        "default_imager": "wsclean",
        "casa_selfcal_supported": False,
        "raw_pipeline_supported": False,
        "bands": {
            "LOFAR_LBA": {
                "aliases": ["LBA", "LOW-BAND", "LOW_BAND"],
                "label": "LOFAR low-band antennas (LBA)",
                "freq_min": 10e6,
                "freq_max": 90e6,
                "cell": "4.5arcsec",
                "imsize": 8192,
                "primary_beam_provider": "wsclean_everybeam",
            },
            "LOFAR_HBA": {
                "aliases": ["HBA", "HIGH-BAND", "HIGH_BAND"],
                "label": "LOFAR high-band antennas (HBA)",
                "freq_min": 110e6,
                "freq_max": 240e6,
                "cell": "2.0arcsec",
                "imsize": 8192,
                "primary_beam_provider": "wsclean_everybeam",
            },
        },
    },
}


def _clean(value: object) -> str:
    return str(value or "").strip()


def normalize_telescope(value: object, *, allow_auto: bool = False) -> str:
    """Return a canonical telescope name or raise an actionable error."""

    raw = _clean(value).lower().replace("_", " ")
    compact = re.sub(r"[^a-z0-9]+", "", raw)

    aliases = {
        "auto": "auto",
        "ugmrt": "ugmrt",
        "gmrt": "ugmrt",
        "giantmetrewaveradiotelescope": "ugmrt",
        "meerkat": "meerkat",
        "mkat": "meerkat",
        "lofar": "lofar",
        "ilt": "lofar",
        "internationallofartelescope": "lofar",
    }
    normalized = aliases.get(compact)

    if normalized == "auto" and not allow_auto:
        normalized = None

    if normalized is None:
        choices = ", ".join(TELESCOPE_CHOICES if allow_auto else TELESCOPE_NAMES)
        raise ValueError(
            f"Unknown telescope '{value}'. Choose one of: {choices}."
        )

    return normalized


def _matching_telescopes(values: Iterable[object]) -> set:
    matches = set()

    for value in values or ():
        text = _clean(value).lower()
        if not text:
            continue
        for telescope, profile in TELESCOPE_PROFILES.items():
            if any(re.search(pattern, text, re.IGNORECASE) for pattern in profile["metadata_patterns"]):
                matches.add(telescope)

    return matches


def infer_telescope_from_antennas(antenna_names: Iterable[object]) -> Optional[str]:
    names = [_clean(name) for name in antenna_names or () if _clean(name)]
    if not names:
        return None

    scores = {}
    for telescope, profile in TELESCOPE_PROFILES.items():
        pattern = re.compile(profile["antenna_pattern"], re.IGNORECASE)
        scores[telescope] = sum(bool(pattern.search(name)) for name in names)

    best = max(scores.values(), default=0)
    winners = [name for name, score in scores.items() if score == best and score > 0]

    # Require a meaningful fraction rather than classifying from one odd
    # antenna name in a heterogeneous or synthetic MS.
    if len(winners) != 1 or best < max(1, (len(names) + 1) // 2):
        return None

    return winners[0]


def resolve_telescope(
    requested: object = "auto",
    *,
    observatory_names: Iterable[object] = (),
    antenna_names: Iterable[object] = (),
) -> dict:
    """Resolve telescope identity, rejecting metadata/override conflicts."""

    requested_name = normalize_telescope(requested or "auto", allow_auto=True)
    metadata_values = [_clean(value) for value in observatory_names or () if _clean(value)]
    metadata_matches = _matching_telescopes(metadata_values)

    if len(metadata_matches) > 1:
        raise ValueError(
            "Measurement Set metadata names multiple supported telescopes: "
            + ", ".join(sorted(metadata_matches))
            + ". Split the data or repair OBSERVATION/TELESCOPE_NAME."
        )

    metadata_name = next(iter(metadata_matches), None)
    antenna_name = infer_telescope_from_antennas(antenna_names)
    if metadata_name and antenna_name and metadata_name != antenna_name:
        raise ValueError(
            "Measurement Set OBSERVATION metadata identifies "
            f"{metadata_name}, but its antenna names identify {antenna_name}. "
            "Repair or split the MS before reduction."
        )
    detected = metadata_name or antenna_name

    if requested_name != "auto":
        if detected and detected != requested_name:
            raise ValueError(
                f"--telescope {requested_name} conflicts with Measurement Set "
                f"metadata identifying {detected}."
            )
        resolved = requested_name
        method = "explicit" if not detected else "explicit_verified"
    else:
        if not detected:
            raise ValueError(
                "Could not identify the telescope from "
                "OBSERVATION/TELESCOPE_NAME or antenna names. Pass "
                "--telescope ugmrt, --telescope meerkat, or --telescope lofar."
            )
        resolved = detected
        method = "observation_metadata" if metadata_name else "antenna_names"

    return {
        "name": resolved,
        "label": TELESCOPE_PROFILES[resolved]["label"],
        "requested": requested_name,
        "method": method,
        "observatory_names": metadata_values,
    }


def get_telescope_profile(telescope: object) -> dict:
    name = normalize_telescope(telescope)
    return copy.deepcopy(TELESCOPE_PROFILES[name])


def telescope_name_from_report(report: Mapping, default: Optional[str] = None) -> Optional[str]:
    value = (report or {}).get("telescope")
    if isinstance(value, Mapping):
        value = value.get("name") or value.get("telescope")
    if not value:
        return default
    try:
        return normalize_telescope(value)
    except ValueError:
        return default


def _band_alias_map(telescope: str) -> dict:
    aliases = {}
    for name, profile in TELESCOPE_PROFILES[telescope]["bands"].items():
        aliases[re.sub(r"[^A-Z0-9]+", "", name.upper())] = name
        for alias in profile.get("aliases", []):
            aliases[re.sub(r"[^A-Z0-9]+", "", str(alias).upper())] = name
    return aliases


def normalize_band(telescope: object, band: object) -> str:
    telescope_name = normalize_telescope(telescope)
    compact = re.sub(r"[^A-Z0-9]+", "", _clean(band).upper())
    if telescope_name == "meerkat" and compact in {"S", "SBAND"}:
        raise ValueError(
            "MeerKAT 'S' is ambiguous because the receiver has five "
            "overlapping tunings. Pass --receiver-band S0, S1, S2, S3, "
            "or S4 explicitly."
        )
    normalized = _band_alias_map(telescope_name).get(compact)
    if normalized is None:
        choices = ", ".join(TELESCOPE_PROFILES[telescope_name]["bands"])
        raise ValueError(
            f"Unknown {telescope_name} receiver band '{band}'. Choose one of: {choices}."
        )
    return normalized


def get_band_profile(telescope: object, band: object) -> dict:
    telescope_name = normalize_telescope(telescope)
    band_name = normalize_band(telescope_name, band)
    profile = copy.deepcopy(TELESCOPE_PROFILES[telescope_name]["bands"][band_name])
    profile.update({"name": band_name, "band": band_name, "telescope": telescope_name})
    return profile


def detect_band(
    telescope: object,
    centre_hz: float,
    *,
    fmin_hz: Optional[float] = None,
    fmax_hz: Optional[float] = None,
    requested_band: object = None,
) -> dict:
    """Resolve a receiver band within one telescope's namespace.

    MeerKAT UHF and L overlap.  A narrow sub-band wholly inside that overlap
    cannot be identified from frequency and deliberately requires
    ``--receiver-band``.
    """

    telescope_name = normalize_telescope(telescope)
    centre = float(centre_hz)
    low = centre if fmin_hz is None else float(fmin_hz)
    high = centre if fmax_hz is None else float(fmax_hz)
    if low > high:
        low, high = high, low

    bands = TELESCOPE_PROFILES[telescope_name]["bands"]

    if requested_band:
        band_name = normalize_band(telescope_name, requested_band)
        profile = bands[band_name]
        if not (
            profile["freq_min"] <= centre <= profile["freq_max"] and
            profile["freq_min"] <= low and
            high <= profile["freq_max"]
        ):
            raise ValueError(
                f"Requested band {band_name} does not contain the observed "
                f"range {low / 1e6:.3f}-{high / 1e6:.3f} MHz."
            )
        candidates = [band_name]
        method = "explicit"
    else:
        centre_candidates = [
            name
            for name, profile in bands.items()
            if profile["freq_min"] <= centre <= profile["freq_max"]
        ]
        method = "frequency"

        if not centre_candidates:
            raise ValueError(
                f"{centre / 1e6:.3f} MHz is outside the configured "
                f"{TELESCOPE_PROFILES[telescope_name]['label']} receiver bands."
            )

        # A centre frequency alone is not enough: concatenated, mislabeled,
        # or otherwise malformed data can extend beyond the selected
        # receiver. Require one profile to contain the complete observed
        # range before assigning it.
        candidates = [
            name
            for name in centre_candidates
            if bands[name]["freq_min"] <= low and high <= bands[name]["freq_max"]
        ]
        if not candidates:
            raise ValueError(
                f"No configured {TELESCOPE_PROFILES[telescope_name]['label']} "
                f"receiver band contains the complete observed range "
                f"{low / 1e6:.3f}-{high / 1e6:.3f} MHz. Split the MS or "
                "repair its spectral metadata before reduction."
            )
        if len(candidates) > 1:
            raise ValueError(
                f"The observed range {low / 1e6:.3f}-{high / 1e6:.3f} MHz "
                f"is ambiguous between {', '.join(candidates)}. Pass "
                "--receiver-band explicitly."
            )

    band_name = candidates[0]
    result = get_band_profile(telescope_name, band_name)
    result.update(
        {
            "method": method,
            "fmin_MHz": low / 1e6,
            "fmax_MHz": high / 1e6,
            "centre_MHz": centre / 1e6,
        }
    )
    return result


def default_imager(telescope: object) -> str:
    return get_telescope_profile(telescope)["default_imager"]


def validate_stage_capabilities(
    telescope: object,
    data_state: object,
    stages: Iterable[str],
    *,
    imager: object = "auto",
    receiver_band: object = None,
) -> None:
    """Fail closed for combinations the implementation cannot calibrate."""

    telescope_name = normalize_telescope(telescope)
    state = _clean(data_state).lower()
    if state not in DATA_STATE_CHOICES:
        raise ValueError(
            f"Unknown data state '{data_state}'. Choose raw or calibrated."
        )

    selected = set(stages or ())
    profile = TELESCOPE_PROFILES[telescope_name]

    band_value = receiver_band
    if isinstance(receiver_band, Mapping):
        band_value = receiver_band.get("name") or receiver_band.get("band")
    band_text = _clean(band_value)
    resolved_band = (
        normalize_band(telescope_name, band_text)
        if band_text and band_text.lower() != "auto"
        else None
    )

    if state == "calibrated" and "calibration" in selected:
        raise ValueError(
            "--data-state calibrated cannot be combined with the raw "
            "cross-calibration stage. Use the calibrated target MS directly "
            "for imaging, or start a new raw-data run."
        )

    if telescope_name == "lofar" and "flagging" in selected:
        raise ValueError(
            "The CASA flagging stage is disabled for LOFAR data. Apply the "
            "validated LINC/DP3/AOFlagger strategy upstream, then provide its "
            "output MS without selecting the pipeline flagging stage."
        )

    if (
        telescope_name == "meerkat" and
        resolved_band is not None and
        resolved_band.startswith("MEERKAT_S") and
        selected.intersection({"imaging", "source_detection"})
    ):
        blocked = ", ".join(
            sorted(selected.intersection({"imaging", "source_detection"}))
        )
        raise ValueError(
            f"MeerKAT {resolved_band[len('MEERKAT_'):]}-band has no "
            "validated primary-beam model in this "
            f"pipeline, so the {blocked} stage(s) cannot produce reliable "
            "publication products. Stop after calibration/self-calibration "
            "or supply a separately validated S-band imaging workflow."
        )

    lofar_raw_science_stages = selected.difference({"ingest", "inspection"})
    if (
        telescope_name == "lofar" and
        state == "raw" and
        lofar_raw_science_stages
    ):
        raise ValueError(
            "Raw LOFAR calibration/science reduction is not supported by "
            "these stages "
            f"({', '.join(sorted(lofar_raw_science_stages))}). Run the "
            "official LINC/DP3 calibration workflow first, then pass its "
            "target MS with --telescope lofar --data-state calibrated. An "
            "inspection-only raw-data run remains available."
        )

    if telescope_name == "lofar" and "selfcal" in selected:
        raise ValueError(
            "The CASA self-calibration stage is disabled for LOFAR because it "
            "does not solve the required LOFAR beam/clock/TEC terms. Use a "
            "DP3/LoSoTo or facet-calibrated target MS and start at imaging."
        )

    if telescope_name != "ugmrt" and "ddcal" in selected:
        raise ValueError(
            "The bundled killMS/DDFacet DD recipe and its CASA primary-beam "
            "handoff are currently validated only for uGMRT. Use a "
            "telescope-specific MeerKAT 3GC or LOFAR facet-calibration "
            "workflow instead of --include-ddcal."
        )

    resolved_imager = default_imager(telescope_name) if imager in (None, "", "auto") else str(imager)
    if telescope_name == "lofar" and "imaging" in selected and resolved_imager != "wsclean":
        raise ValueError(
            "LOFAR imaging requires WSClean/EveryBeam in this pipeline. Use "
            "--imager wsclean (or leave --imager auto)."
        )
    if (
        telescope_name == "ugmrt" and
        "imaging" in selected and
        resolved_imager != "casa"
    ):
        # uGMRT is imaged only with CASA tclean. Its measured NCRA primary
        # beam is applied by the CASA path, its self-calibration and
        # direction-dependent handoffs are built on CASA images, and nothing
        # in the WSClean path reproduces them.
        raise ValueError(
            f"{profile['label']} publication imaging is validated only with "
            "the CASA backend (tclean) in this pipeline. Use --imager casa "
            "(or leave --imager auto); the WSClean path has no uGMRT "
            "primary-beam correction."
        )
    if (
        telescope_name == "meerkat" and
        "imaging" in selected and
        resolved_imager not in {"casa", "wsclean"}
    ):
        raise ValueError(
            f"Unknown imager {resolved_imager!r} for {profile['label']}. "
            "Use --imager wsclean (the default) or --imager casa."
        )


def all_band_profiles() -> dict:
    result = {}
    for telescope, profile in TELESCOPE_PROFILES.items():
        for band, values in profile["bands"].items():
            item = copy.deepcopy(values)
            item["telescope"] = telescope
            result[band] = item
    return result
