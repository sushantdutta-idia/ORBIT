#!/usr/bin/env python3
# ============================================================
#
# ddcal.py
#
# Direction-dependent calibration workflow for the uGMRT continuum pipeline.
#
# Intended order:
#   inspect.py -> flagging.py -> calibration.py -> selfcal.py
#   -> imaging.py -> source_detection.py -> ddcal.py
#
# This module prepares DD directions from the final continuum
# image, writes source/direction products, writes solver command
# files, optionally executes the selected DD solver, optionally
# re-images, and writes a restartable JSON report.
#
# External tools are intentionally optional:
#   - PyBDSF for source extraction
#   - CubiCal/gocubical for DD gain solving
#   - KillMS + DDFacet for facet/DD solving and imaging
#   - CASA for post-DD imaging with imaging.py
#
# Author: Sushant Dutta
# ============================================================

import os
import sys


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of failing inside a stage.
PIPELINE_INTERFACE_VERSION = 1

SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

if (
    sys.path and
    os.path.abspath(sys.path[0]) == SCRIPT_DIR
):

    sys.path.pop(0)

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


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


def load_imaging_safety_module():

    module_path = Path(SCRIPT_DIR) / "imaging_safety.py"
    spec = importlib.util.spec_from_file_location(
        "ugmrt_ddcal_imaging_safety",
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
DD_TASK_CONFIG = (
    getattr(CONFIG, "DD_TASK_CONFIG", {})
    if CONFIG is not None
    else {}
)
QA_PRODUCTS = load_qa_products_module()
IMAGING_SAFETY = load_imaging_safety_module()


# ============================================================
# DEFAULTS
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


if CONFIG is not None:

    DD_CONFIG = getattr(
        CONFIG,
        "DD_CONFIG",
        DD_CONFIG,
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


def save_json(data, path):

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")

    with open(temporary, "w") as f:

        json.dump(
            to_builtin(data),
            f,
            indent=4,
        )
        f.write("\n")

    os.replace(temporary, path)


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


def clean_name(name):

    name = str(name or "target").strip()
    cleaned = []

    for char in name:

        if char.isalnum() or char in "._-":

            cleaned.append(char)

        else:

            cleaned.append("_")

    return "".join(cleaned).strip("_") or "target"


def normalize_key(key):

    return "".join(
        char.lower()
        for char in str(key)
        if char.isalnum()
    )


def get_first(row, candidates):

    normalized = {
        normalize_key(k): v
        for k, v in row.items()
    }

    for candidate in candidates:

        key = normalize_key(candidate)

        if key in normalized:

            value = normalized[key]

            if value not in ("", None):

                return value

    return None


def parse_float(value):

    if value is None:

        return None

    try:

        return float(str(value).strip())

    except Exception:

        return None


def parse_angle(value, is_ra=False):

    if value is None:

        return None

    text = str(value).strip()

    if not text:

        return None

    try:

        return float(text)

    except ValueError:

        pass

    cleaned = (
        text.replace("h", ":")
        .replace("m", ":")
        .replace("s", "")
        .replace("d", ":")
        .replace("'", ":")
        .replace('"', "")
    )

    parts = [
        p
        for p in cleaned.replace(" ", ":").split(":")
        if p
    ]

    if len(parts) == 0:

        return None

    sign = 1.0

    if parts[0].startswith("-"):

        sign = -1.0
        parts[0] = parts[0][1:]

    elif parts[0].startswith("+"):

        parts[0] = parts[0][1:]

    try:

        numbers = [float(p) for p in parts]

    except ValueError:

        return None

    while len(numbers) < 3:

        numbers.append(0.0)

    deg = (
        abs(numbers[0]) +
        numbers[1] / 60.0 +
        numbers[2] / 3600.0
    )

    if is_ra:

        return deg * 15.0

    return sign * deg


def angular_separation_deg(ra1, dec1, ra2, dec2):

    r1 = math.radians(ra1)
    d1 = math.radians(dec1)
    r2 = math.radians(ra2)
    d2 = math.radians(dec2)

    cos_sep = (
        math.sin(d1) * math.sin(d2) +
        math.cos(d1) * math.cos(d2) * math.cos(r1 - r2)
    )

    cos_sep = max(-1.0, min(1.0, cos_sep))

    return math.degrees(math.acos(cos_sep))


def which_any(names):

    for name in names:

        path = shutil.which(name)

        if path:

            return path

    return None


def command_to_string(command):

    return " ".join(
        shlex.quote(str(part))
        for part in command
    )


def split_command(command):

    if isinstance(command, list):

        return [
            str(part)
            for part in command
        ]

    return shlex.split(str(command))


def now_utc():

    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def parse_quantity(value, unit):
    """Return a positive numeric value expressed in the requested unit."""

    if isinstance(value, (int, float)):
        return float(value)

    text = str(value or "").strip().lower().replace(" ", "")
    match = re.fullmatch(
        r"([0-9]+(?:\.[0-9]+)?)([a-z]+)?",
        text,
    )

    if not match:
        raise ValueError(f"Cannot parse quantity: {value!r}")

    number = float(match.group(1))
    suffix = match.group(2) or unit

    conversions = {
        "min": {
            "s": 1.0 / 60.0,
            "sec": 1.0 / 60.0,
            "secs": 1.0 / 60.0,
            "second": 1.0 / 60.0,
            "seconds": 1.0 / 60.0,
            "m": 1.0,
            "min": 1.0,
            "mins": 1.0,
            "minute": 1.0,
            "minutes": 1.0,
            "h": 60.0,
            "hr": 60.0,
            "hour": 60.0,
            "hours": 60.0,
        },
        "mhz": {
            "hz": 1.0e-6,
            "khz": 1.0e-3,
            "mhz": 1.0,
            "ghz": 1.0e3,
        },
        "arcsec": {
            "arcsec": 1.0,
            "asec": 1.0,
            "arcmin": 60.0,
            "amin": 60.0,
            "deg": 3600.0,
        },
    }

    try:
        result = number * conversions[unit][suffix]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported unit in quantity {value!r}; expected {unit}."
        ) from exc

    if result <= 0:
        raise ValueError(f"Quantity must be positive: {value!r}")

    return result


# ============================================================
# SOURCE CATALOGUE
# ============================================================

class Source:

    def __init__(self,
                 name,
                 ra_deg,
                 dec_deg,
                 flux_jy,
                 peak_jy=None):

        self.name = name
        self.ra_deg = float(ra_deg)
        self.dec_deg = float(dec_deg)
        self.flux_jy = float(flux_jy)
        self.peak_jy = (
            None
            if peak_jy is None
            else float(peak_jy)
        )

    def to_dict(self):

        return {
            "name": self.name,
            "ra_deg": self.ra_deg,
            "dec_deg": self.dec_deg,
            "flux_jy": self.flux_jy,
            "peak_jy": self.peak_jy,
        }


class Direction:

    def __init__(self,
                 name,
                 ra_deg,
                 dec_deg,
                 flux_jy,
                 nsources=1,
                 sources=None):

        self.name = name
        self.ra_deg = float(ra_deg)
        self.dec_deg = float(dec_deg)
        self.flux_jy = float(flux_jy)
        self.nsources = int(nsources)
        self.sources = sources or []

    def to_dict(self):

        return {
            "name": self.name,
            "ra_deg": self.ra_deg,
            "dec_deg": self.dec_deg,
            "flux_jy": self.flux_jy,
            "nsources": self.nsources,
            "sources": self.sources,
        }


def read_source_catalogue(path):

    path = Path(path)

    if not path.exists():

        raise FileNotFoundError(
            f"Source catalogue not found:\n{path}"
        )

    if path.suffix.lower() not in [".csv", ".txt"]:

        raise RuntimeError(
            "DDCal catalogue parsing currently expects CSV/TXT."
            "Use PyBDSF CSV output or pass --directions for a pre-made direction list."
        )

    sources = []

    with open(path, newline="") as f:

        sample = f.read(4096)
        f.seek(0)

        dialect = csv.Sniffer().sniff(
            sample,
            delimiters=", \t",
        )

        reader = csv.DictReader(
            (
                line
                for line in f
                if line.strip() and not line.lstrip().startswith("#")
            ),
            dialect=dialect,
        )

        for index, row in enumerate(reader, start=1):

            ra = parse_angle(
                get_first(
                    row,
                    [
                        "RA",
                        "RA_deg",
                        "RAJ2000",
                        "Source_RA",
                        "ra",
                    ],
                ),
                is_ra=True,
            )

            dec = parse_angle(
                get_first(
                    row,
                    [
                        "DEC",
                        "Dec",
                        "DEC_deg",
                        "DEJ2000",
                        "Source_DEC",
                        "dec",
                    ],
                ),
                is_ra=False,
            )

            flux_value = get_first(
                row,
                [
                    "Total_flux",
                    "Total_flux_Jy",
                    "Flux",
                    "Flux_Jy",
                    "I",
                    "S",
                    "Peak_flux",
                ],
            )

            peak_value = get_first(
                row,
                [
                    "Peak_flux",
                    "Peak_flux_Jy",
                    "Peak",
                ],
            )

            flux = parse_float(flux_value)
            peak = parse_float(peak_value)

            if flux is None and peak is not None:

                flux = peak

            if ra is None or dec is None or flux is None:

                continue

            flux_key = normalize_key(
                next(
                    (
                        k
                        for k, v in row.items()
                        if str(v) == str(flux_value)
                    ),
                    "",
                )
            )

            if "mj" in flux_key:

                flux = flux / 1000.0

            if peak is not None:

                peak_key = normalize_key(
                    next(
                        (
                            k
                            for k, v in row.items()
                            if str(v) == str(peak_value)
                        ),
                        "",
                    )
                )

                if "mj" in peak_key:

                    peak = peak / 1000.0

            name = get_first(
                row,
                [
                    "Source_id",
                    "Source_name",
                    "Name",
                    "ID",
                ],
            )

            if not name:

                name = f"src{index:04d}"

            sources.append(
                Source(
                    clean_name(name),
                    ra,
                    dec,
                    flux,
                    peak_jy=peak,
                )
            )

    return sources


def write_source_catalogue(sources, path):

    with open(path, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "name",
                "ra_deg",
                "dec_deg",
                "flux_jy",
                "peak_jy",
            ],
        )

        writer.writeheader()

        for source in sources:

            writer.writerow(source.to_dict())


def write_directions_csv(directions, path):

    with open(path, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "name",
                "ra_deg",
                "dec_deg",
                "flux_jy",
                "nsources",
                "sources",
            ],
        )

        writer.writeheader()

        for direction in directions:

            row = direction.to_dict()
            row["sources"] = ";".join(direction.sources)
            writer.writerow(row)


def read_directions_csv(path):

    directions = []

    with open(path, newline="") as f:

        reader = csv.DictReader(
            line
            for line in f
            if line.strip() and not line.lstrip().startswith("#")
        )

        for index, row in enumerate(reader, start=1):

            ra = parse_angle(
                get_first(row, ["ra_deg", "RA", "RAJ2000"]),
                is_ra=True,
            )
            dec = parse_angle(
                get_first(row, ["dec_deg", "DEC", "DEJ2000"]),
                is_ra=False,
            )
            flux = parse_float(
                get_first(row, ["flux_jy", "Flux", "Total_flux"])
            )

            if ra is None or dec is None:

                continue

            if flux is None:

                flux = 0.0

            name = get_first(row, ["name", "Name", "id", "ID"])

            if not name:

                name = f"dir{index:03d}"

            sources = get_first(row, ["sources", "Sources"]) or ""

            directions.append(
                Direction(
                    clean_name(name),
                    ra,
                    dec,
                    flux,
                    nsources=int(
                        parse_float(
                            get_first(row, ["nsources", "n"])
                        ) or 1
                    ),
                    sources=[
                        s
                        for s in sources.split(";")
                        if s
                    ],
                )
            )

    return directions


def cluster_sources_into_directions(sources,
                                    min_flux_jy,
                                    max_directions,
                                    merge_radius_arcmin,
                                    uncovered=None):

    bright = [
        source
        for source in sources
        if source.flux_jy >= min_flux_jy
    ]

    bright.sort(
        key=lambda item: item.flux_jy,
        reverse=True,
    )

    directions = []
    merge_radius_deg = merge_radius_arcmin / 60.0

    for source in bright:

        matched = None

        for direction in directions:

            sep = angular_separation_deg(
                source.ra_deg,
                source.dec_deg,
                direction.ra_deg,
                direction.dec_deg,
            )

            if sep <= merge_radius_deg:

                matched = direction
                break

        if matched:

            total_flux = matched.flux_jy + source.flux_jy

            # Average right ascension on the unit circle.  A plain weighted
            # mean of degrees places the centroid of a pair straddling RA=0
            # near 180 degrees, i.e. on the opposite side of the sky.
            east = (
                matched.flux_jy * math.cos(math.radians(matched.ra_deg)) +
                source.flux_jy * math.cos(math.radians(source.ra_deg))
            )
            north = (
                matched.flux_jy * math.sin(math.radians(matched.ra_deg)) +
                source.flux_jy * math.sin(math.radians(source.ra_deg))
            )
            matched.ra_deg = math.degrees(math.atan2(north, east)) % 360.0

            matched.dec_deg = (
                matched.dec_deg * matched.flux_jy +
                source.dec_deg * source.flux_jy
            ) / total_flux

            matched.flux_jy = total_flux
            matched.nsources += 1
            matched.sources.append(source.name)

        else:

            if len(directions) >= max_directions:

                # This source is bright enough to need a direction but the
                # direction budget is exhausted.  It will not be DD-corrected
                # and its artefacts stay in the image, so record it instead of
                # dropping it silently.
                if uncovered is not None:
                    uncovered.append(source)

                continue

            directions.append(
                Direction(
                    f"dir{len(directions) + 1:03d}",
                    source.ra_deg,
                    source.dec_deg,
                    source.flux_jy,
                    nsources=1,
                    sources=[source.name],
                )
            )

    directions.sort(
        key=lambda item: item.flux_jy,
        reverse=True,
    )

    for index, direction in enumerate(directions, start=1):

        direction.name = f"dir{index:03d}"

    return directions


def write_ds9_regions(directions,
                      path,
                      radius_arcmin=5.0,
                      include_text=True):

    with open(path, "w") as f:

        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color=green dashlist=8 3 width=2 ")
        f.write("font=\"helvetica 10 normal roman\" ")
        f.write("select=1 highlite=1 dash=0 fixed=0 edit=1 ")
        f.write("move=1 delete=1 include=1 source=1\n")
        f.write("fk5\n")

        for direction in directions:

            label = (
                f" # text={{{direction.name} "
                f"{direction.flux_jy:.3f} Jy}}"
                if include_text
                else ""
            )

            f.write(
                "circle("
                f"{direction.ra_deg:.8f},"
                f"{direction.dec_deg:.8f},"
                f"{radius_arcmin:.3f}'"
                f"){label}\n"
            )


def write_tigger_lsm(directions, path):

    with open(path, "w") as f:

        f.write("#format:name ra_d dec_d i q u v spi freq0\n")

        for direction in directions:

            f.write(
                f"{direction.name} "
                f"{direction.ra_deg:.10f} "
                f"{direction.dec_deg:.10f} "
                f"{direction.flux_jy:.10e} "
                "0 0 0 0 1.0e9\n"
            )


# ============================================================
# DD CALIBRATION MANAGER
# ============================================================

class DDCalManager:

    def __init__(self,
                 ms=None,
                 image=None,
                 catalog=None,
                 directions=None,
                 outdir="Products",
                 target=None,
                 mode=None,
                 engine=None,
                 execute=None,
                 strict=None,
                 ncpu=None,
                 data_column=None,
                 model_column=None,
                 min_flux_jy=None,
                 max_directions=None,
                 merge_radius_arcmin=None,
                 facet_radius_arcmin=None,
                 solint_time=None,
                 solint_freq=None,
                 solver_command_template=None,
                 apply_command_template=None,
                 killms_command_template=None,
                 ddfacet_command_template=None,
                 skip_source_extraction=False,
                 skip_solver=False,
                 skip_imaging=None,
                 imaging_script=None,
                 casa=None,
                 provenance_fingerprint=None):

        self.explicit_ms = ms
        self.explicit_image = image
        self.explicit_catalog = catalog
        self.explicit_directions = directions
        self.provenance_fingerprint = (
            str(provenance_fingerprint).strip()
            if provenance_fingerprint
            else None
        )
        self.outdir = Path(outdir)
        self.target = target
        self.requested_mode = (
            mode or
            DD_TASK_CONFIG.get("mode", "auto")
        ).lower()
        self.mode = None
        self.requested_engine = (
            engine or
            DD_TASK_CONFIG.get(
                "engine",
                "auto",
            )
        ).lower()
        self.engine = self.requested_engine
        self.execute = (
            DD_TASK_CONFIG.get(
                "execute",
                True,
            )
            if execute is None
            else execute
        )
        self.strict = (
            DD_TASK_CONFIG.get(
                "strict",
                False,
            )
            if strict is None
            else strict
        )
        self.ncpu = (
            RUN_CONFIG.get(
                "ncpu",
                32,
            )
            if ncpu is None
            else ncpu
        )
        self.data_column = (
            data_column or
            DD_TASK_CONFIG.get(
                "data_column",
                "auto",
            )
        )
        self.model_column = (
            model_column or
            DD_TASK_CONFIG.get(
                "model_column",
                "MODEL_DATA",
            )
        )
        self.solver_command_template = solver_command_template
        self.apply_command_template = apply_command_template
        self.killms_command_template = killms_command_template
        self.ddfacet_command_template = ddfacet_command_template
        self.skip_source_extraction = skip_source_extraction
        self.skip_solver = skip_solver
        self.skip_imaging = (
            not DD_TASK_CONFIG.get(
                "post_dd_imaging",
                True,
            )
            if skip_imaging is None
            else skip_imaging
        )
        self.imaging_script = (
            Path(imaging_script)
            if imaging_script
            else Path(__file__).parent / "imaging.py"
        )
        self.casa = (
            casa or
            RUN_CONFIG.get(
                "casa",
                "casa",
            )
        )

        self.root = self.outdir / "ddcal"
        self.catalogdir = self.root / "catalogs"
        self.directiondir = self.root / "directions"
        self.configdir = self.root / "configs"
        self.logdir = self.root / "logs"
        self.imagedir = self.root / "images"
        self.plotdir = self.root / "plots"
        self.reportdir = self.root / "reports"
        self.checkdir = self.root / "checkpoints"
        self.soldir = self.root / "solutions"

        for directory in [
            self.root,
            self.catalogdir,
            self.directiondir,
            self.configdir,
            self.logdir,
            self.imagedir,
            self.plotdir,
            self.reportdir,
            self.checkdir,
            self.soldir,
        ]:

            directory.mkdir(parents=True, exist_ok=True)

        self.inspection = load_json(
            self.outdir /
            "inspection" /
            "inspection_report.json",
            {},
        )

        self.calibration = load_json(
            self.outdir /
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

        self.imaging = load_json(
            self.outdir /
            "imaging" /
            "reports" /
            "imaging_report.json",
            {},
        )

        self.source_detection = load_json(
            self.outdir /
            "source_detection" /
            "reports" /
            "source_detection_report.json",
            {},
        )

        if str(self.data_column).lower() == "auto":
            self.data_column = (
                self.imaging.get("datacolumn") or
                "CORRECTED_DATA"
            )

        column_aliases = {
            "data": "DATA",
            "corrected": "CORRECTED_DATA",
            "corrected_data": "CORRECTED_DATA",
            "model": "MODEL_DATA",
            "model_data": "MODEL_DATA",
        }
        self.data_column = column_aliases.get(
            str(self.data_column).strip().lower(),
            str(self.data_column).strip(),
        )

        self.band = self.detect_band()
        self.config = dict(
            DD_CONFIG.get(
                self.band,
                DD_CONFIG["B4"],
            )
        )

        if min_flux_jy is not None:
            self.config["min_flux_jy"] = min_flux_jy

        if max_directions is not None:
            self.config["max_directions"] = max_directions

        if merge_radius_arcmin is not None:
            self.config["merge_radius_arcmin"] = merge_radius_arcmin

        if facet_radius_arcmin is not None:
            self.config["facet_radius_arcmin"] = facet_radius_arcmin

        if solint_time is not None:
            self.config["solint_time"] = solint_time

        if solint_freq is not None:
            self.config["solint_freq"] = solint_freq

        self.ms = self.select_measurement_set()
        self.target = self.select_target()
        self.image = self.select_image()

        self.catalog_csv = (
            self.catalogdir /
            f"{clean_name(self.target)}_pybdsf_sources.csv"
        )

        self.catalog_fits = (
            self.catalogdir /
            f"{clean_name(self.target)}_pybdsf_sources.fits"
        )

        self.catalog_region = (
            self.catalogdir /
            f"{clean_name(self.target)}_pybdsf_sources.reg"
        )

        self.directions_csv = (
            self.directiondir /
            "dd_directions.csv"
        )

        self.directions_region = (
            self.directiondir /
            "dd_directions.reg"
        )

        self.solver_directions_region = (
            self.directiondir /
            "dd_solver_directions.reg"
        )

        self.facet_region = (
            self.directiondir /
            "dd_facets.reg"
        )

        self.skymodel = (
            self.directiondir /
            "dd_directions.lsm"
        )

        self.command_script = (
            self.configdir /
            "ddcal_commands.sh"
        )

        self.command_records = []
        self.status = {}
        self.workflow_status = "initializing"
        self.failure_reason = None
        self.dd_image_products = []
        self.pbcor_products = []
        self.dd_metrics = {}

    # --------------------------------------------------------
    # Metadata
    # --------------------------------------------------------

    def detect_band(self):

        candidates = [
            self.imaging.get("band"),
            self.selfcal.get("band"),
            self.calibration.get("band"),
            self.inspection.get("band"),
        ]

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

            if band in DD_CONFIG:

                return band

        print(
            "WARNING: Could not determine band. "
            "Using B4 DD defaults."
        )

        return "B4"

    # --------------------------------------------------------

    def select_measurement_set(self):

        candidates = [
            self.explicit_ms,
            self.imaging.get("measurement_set"),
            self.selfcal.get("final_ms"),
        ]

        for candidate in candidates:

            if candidate and Path(candidate).exists():

                return str(candidate)

        split_targets = load_json(
            self.outdir /
            "calibration" /
            "split_targets.json",
            {},
        )

        if split_targets:

            if self.target and self.target in split_targets:

                candidate = split_targets[self.target]

                if Path(candidate).exists():

                    return str(candidate)

            for candidate in split_targets.values():

                if Path(candidate).exists():

                    return str(candidate)

        raise FileNotFoundError(
            "Could not choose an input Measurement Set. "
            "Pass --ms or run selfcal.py/imaging.py first."
        )

    # --------------------------------------------------------

    def select_target(self):

        if self.target:

            return self.target

        candidates = [
            self.imaging.get("target"),
            self.selfcal.get("target"),
        ]

        targets = (
            self.calibration.get("targets") or
            self.inspection.get("target_fields") or
            self.inspection.get("fields", {}).get("targets") or
            []
        )

        if isinstance(targets, str):

            candidates.append(targets)

        elif targets:

            candidates.append(targets[0])

        for candidate in candidates:

            if candidate:

                return str(candidate)

        return Path(self.ms).stem

    # --------------------------------------------------------

    def select_image(self):

        if self.explicit_image:

            image = Path(self.explicit_image)

            if not image.exists():

                raise FileNotFoundError(
                    f"Image not found:\n{image}"
                )

            return str(image)

        products = self.imaging.get("products", {})
        fits_products = products.get("fits", [])

        preferred = [
            "pbcor.image.tt0",
            "image.tt0.pbcor",
            "commonbeam.image",
            "image.tt0",
        ]

        for token in preferred:

            for candidate in fits_products:

                if token in Path(candidate).name and Path(candidate).exists():

                    return str(candidate)

        for key in [
            "pbcor_image",
            "common_beam_image",
            "restored_image",
        ]:

            candidate = products.get(key)

            if candidate and Path(candidate).exists():

                return str(candidate)

        return None

    # --------------------------------------------------------
    # Checkpoints
    # --------------------------------------------------------

    def checkpoint_exists(self, name):

        return (
            self.checkdir /
            f"{name}.done"
        ).exists()

    # --------------------------------------------------------

    def checkpoint_provenance_path(self, name):

        return self.checkdir / f"{name}.provenance.json"

    # --------------------------------------------------------

    def mark_checkpoint(self, name):

        (
            self.checkdir /
            f"{name}.done"
        ).touch()

    # --------------------------------------------------------

    def file_fingerprint(self, path):

        if not path or not Path(path).is_file():
            return None

        return IMAGING_SAFETY.file_content_identity(path)

    # --------------------------------------------------------

    def measurement_set_stat_identity(self):

        return IMAGING_SAFETY.path_stat_fingerprint(self.ms)

    # --------------------------------------------------------

    def checkpoint_context(self):

        return {
            "schema_version": 1,
            "solver_signature": self.solver_signature(),
            "measurement_set": str(Path(self.ms).expanduser().resolve()),
            "measurement_set_stat_sha256": (
                self.measurement_set_stat_identity()
            ),
            "source_catalog_content_identity": self.file_fingerprint(
                self.catalog_csv
            ),
        }

    # --------------------------------------------------------

    @staticmethod
    def _identity_path(value):

        return str(Path(value).expanduser().resolve())

    # --------------------------------------------------------

    def product_identity(self, path):

        product = Path(path).expanduser().resolve()
        if product.is_file():
            return {
                "kind": "file",
                **IMAGING_SAFETY.file_content_identity(product),
            }
        if product.is_dir():
            return {
                "kind": "directory",
                "path": str(product),
                "stat_sha256": IMAGING_SAFETY.path_stat_fingerprint(
                    product
                ),
            }
        raise ValueError(f"DD checkpoint product is missing: {product}")

    # --------------------------------------------------------

    def product_identities(self, paths):

        unique = []
        for value in paths or ():
            if not value:
                continue
            path = Path(value).expanduser().resolve()
            if path not in unique:
                unique.append(path)

        return {
            str(path): self.product_identity(path)
            for path in unique
        }

    # --------------------------------------------------------

    def mark_bound_checkpoint(self, name, products=(), inputs=()):
        """Publish a checkpoint only after its inputs and outputs are bound."""

        record = {
            **self.checkpoint_context(),
            "checkpoint": name,
            "created_utc": now_utc(),
            "input_identities": self.product_identities(inputs),
            "product_identities": self.product_identities(products),
        }
        save_json(record, self.checkpoint_provenance_path(name))
        self.mark_checkpoint(name)
        return record

    # --------------------------------------------------------

    def bound_checkpoint_status(self, name, products=(), inputs=()):

        if not self.checkpoint_exists(name):
            return False, "completion marker is missing"

        provenance_path = self.checkpoint_provenance_path(name)
        record = load_json(provenance_path, {})
        if not isinstance(record, dict):
            return False, "checkpoint provenance is unreadable"

        try:
            current_context = self.checkpoint_context()
        except (OSError, TypeError, ValueError) as exc:
            return False, f"checkpoint input identity is unavailable: {exc}"

        for key, value in current_context.items():
            if record.get(key) != value:
                return False, f"checkpoint input provenance changed for {key}"

        try:
            current_inputs = self.product_identities(inputs)
            current_products = self.product_identities(products)
        except (OSError, TypeError, ValueError) as exc:
            return False, f"checkpoint product identity is invalid: {exc}"

        if record.get("input_identities") != current_inputs:
            return False, "checkpoint input-product identities changed"
        if record.get("product_identities") != current_products:
            return False, "checkpoint product identities changed"

        return True, "checkpoint inputs and products match"

    # --------------------------------------------------------

    def _archive_destination(self, path):

        token = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")
        candidate = path.with_name(f"{path.name}.stale.{token}.{os.getpid()}")
        index = 1
        while candidate.exists():
            candidate = path.with_name(
                f"{path.name}.stale.{token}.{os.getpid()}.{index}"
            )
            index += 1
        return candidate

    # --------------------------------------------------------

    def archive_owned_solver_product(self, path):
        """Preserve, rather than delete, one stale pipeline-owned output."""

        product = Path(path).expanduser().resolve()
        if not product.exists():
            return None

        owned_root = self.root.expanduser().resolve()
        inside_dd_root = False
        try:
            product.relative_to(owned_root)
            inside_dd_root = True
        except ValueError:
            pass

        solution_file = Path(
            self.solver_paths()["solution_file"]
        ).expanduser().resolve()
        if not inside_dd_root and product != solution_file:
            raise RuntimeError(
                "Refusing to archive a stale DD product outside the "
                f"pipeline-owned DD tree: {product}"
            )

        archived = self._archive_destination(product)
        product.rename(archived)
        return archived

    # --------------------------------------------------------

    def archive_checkpoint(self, name):

        archived = []
        for path in [
            self.checkdir / f"{name}.done",
            self.checkpoint_provenance_path(name),
        ]:
            if path.exists():
                destination = self._archive_destination(path)
                path.rename(destination)
                archived.append(destination)
        return archived

    # --------------------------------------------------------

    def invalidate_solver_state(self, reason, products=()):
        """Invalidate DD solver descendants without deleting evidence."""

        print(
            "WARNING: DD solver provenance is stale; preserving current "
            f"products and recomputing them ({reason})."
        )

        names = {
            "solver",
            "solution_qa",
            "post_dd_imaging",
            "ddcal",
        }
        names.update(
            path.name[:-len(".done")]
            for path in self.checkdir.glob("solver_*.done")
        )
        for name in sorted(names):
            self.archive_checkpoint(name)

        for product in products or ():
            self.archive_owned_solver_product(product)

        return False

    # --------------------------------------------------------

    def invalidate_post_dd_state(self, reason, products=()):

        print(
            "WARNING: DD primary-beam provenance is stale; preserving "
            f"current products and recomputing them ({reason})."
        )
        for name in ["post_dd_imaging", "ddcal"]:
            self.archive_checkpoint(name)
        for product in products or ():
            self.archive_owned_solver_product(product)

        return False

    # --------------------------------------------------------

    def direction_input_signature(self):

        payload = {
            "catalog": self.file_fingerprint(self.catalog_csv),
            "explicit_directions": (
                self.file_fingerprint(self.explicit_directions)
                if self.explicit_directions
                else None
            ),
            "min_flux_jy": self.config["min_flux_jy"],
            "max_directions": self.config["max_directions"],
            "merge_radius_arcmin": self.config[
                "merge_radius_arcmin"
            ],
        }
        encoded = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")

        return hashlib.sha256(encoded).hexdigest(), payload

    # --------------------------------------------------------

    def mark_direction_checkpoint(self):

        signature, inputs = self.direction_input_signature()
        save_json(
            {
                "created_utc": now_utc(),
                "signature": signature,
                "inputs": inputs,
            },
            self.checkdir / "directions_inputs.json",
        )
        self.mark_checkpoint("directions")

    # --------------------------------------------------------

    def direction_checkpoint_is_valid(self):

        if (
            not self.checkpoint_exists("directions") or
            not self.directions_csv.exists()
        ):
            return False

        record = load_json(
            self.checkdir / "directions_inputs.json",
            {},
        )
        signature, _inputs = self.direction_input_signature()

        return record.get("signature") == signature

    # --------------------------------------------------------
    # Command execution
    # --------------------------------------------------------

    def run_command(self,
                    name,
                    command,
                    required=False):

        command = split_command(command)
        command_text = command_to_string(command)
        log_file = self.logdir / f"{name}.log"

        record = {
            "name": name,
            "command": command_text,
            "log": str(log_file),
            "executed": False,
            "returncode": None,
            "status": "prepared",
        }

        self.command_records.append(record)

        if not self.execute:

            print(
                f"Prepared command [{name}]: {command_text}"
            )

            return record

        exe = shutil.which(command[0]) or command[0]

        if not Path(exe).exists() and shutil.which(command[0]) is None:

            message = (
                f"Command not found for stage {name}: "
                f"{command[0]}"
            )

            record["status"] = "missing"

            print(f"WARNING: {message}")

            if required or self.strict:

                raise RuntimeError(message)

            return record

        print(
            f"Running [{name}]: {command_text}"
        )

        with open(log_file, "w") as log:

            proc = subprocess.run(
                command,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
            )

        record["executed"] = True
        record["returncode"] = proc.returncode
        record["status"] = (
            "ok"
            if proc.returncode == 0
            else "failed"
        )

        if proc.returncode != 0 and (required or self.strict):

            raise RuntimeError(
                f"Stage {name} failed. See log:\n{log_file}"
            )

        return record

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    def print_summary(self):

        print("\n")
        print("=" * 70)
        print("uGMRT Direction-Dependent Calibration")
        print("=" * 70)
        print(f"MS              : {self.ms}")
        print(f"Target          : {self.target}")
        print(f"Band            : {self.band}")
        print(f"Reference image : {self.image}")
        print(f"DD mode         : {self.requested_mode}")
        print(f"Engine          : {self.requested_engine}")
        print(f"Execute         : {self.execute}")
        print(f"Output          : {self.root}")
        print(f"Min flux        : {self.config['min_flux_jy']} Jy")
        print(f"Max directions  : {self.config['max_directions']}")
        print(
            f"Merge radius    : "
            f"{self.config['merge_radius_arcmin']} arcmin"
        )
        print(
            f"Solint          : "
            f"{self.config['solint_time']} / "
            f"{self.config['solint_freq']}"
        )
        print("=" * 70)

    # --------------------------------------------------------
    # PyBDSF source extraction
    # --------------------------------------------------------

    def write_pybdsf_script(self):

        script = (
            self.configdir /
            "run_pybdsf.py"
        )

        image = self.image
        csv_path = str(self.catalog_csv)
        fits_path = str(self.catalog_fits)
        reg_path = str(self.catalog_region)

        lines = [
            "import bdsf",
            f"image = {image!r}",
            f"csv_path = {csv_path!r}",
            f"fits_path = {fits_path!r}",
            f"reg_path = {reg_path!r}",
            "img = bdsf.process_image(",
            "    image,",
            "    adaptive_rms_box=True,",
            "    atrous_do=True,",
            "    thresh_isl=4.0,",
            "    thresh_pix=5.0,",
            ")",
            "img.write_catalog(",
            "    outfile=csv_path,",
            "    catalog_type='srl',",
            "    format='csv',",
            "    clobber=True,",
            ")",
            "img.write_catalog(",
            "    outfile=fits_path,",
            "    catalog_type='srl',",
            "    format='fits',",
            "    clobber=True,",
            ")",
            "img.write_catalog(",
            "    outfile=reg_path,",
            "    catalog_type='srl',",
            "    format='ds9',",
            "    clobber=True,",
            ")",
        ]

        script.write_text("\n".join(lines) + "\n")

        return script

    # --------------------------------------------------------

    def extract_sources(self):

        if self.explicit_catalog:

            print(
                f"Using supplied source catalogue: "
                f"{self.explicit_catalog}"
            )

            self.catalog_csv = Path(self.explicit_catalog)
            self.mark_checkpoint("source_catalog")
            return

        detected_catalog = (
            self.source_detection.get("products", {})
            .get("source_catalog_csv")
        )

        if (
            self.source_detection.get("status") == "success" and
            detected_catalog and
            Path(detected_catalog).is_file()
        ):

            self.catalog_csv = Path(detected_catalog)
            print(
                "Reusing source-detection catalogue: "
                f"{self.catalog_csv}"
            )
            self.mark_checkpoint("source_catalog")
            return

        if self.checkpoint_exists("source_catalog") and self.catalog_csv.exists():

            print("Source catalogue already exists.")
            return

        if self.skip_source_extraction:

            print("Skipping source extraction by request.")
            return

        if not self.image:

            message = (
                "No final image was found for PyBDSF source "
                "extraction. Pass --image or --catalog."
            )

            if self.strict:

                raise RuntimeError(message)

            print(f"WARNING: {message}")
            return

        pybdsf_script = self.write_pybdsf_script()

        command = [
            sys.executable,
            str(pybdsf_script),
        ]

        record = self.run_command(
            "pybdsf",
            command,
            required=False,
        )

        if self.catalog_csv.exists():

            self.mark_checkpoint("source_catalog")

        elif record["status"] in ["missing", "failed"] and self.strict:

            raise RuntimeError(
                "PyBDSF source extraction did not create a CSV "
                f"catalogue: {self.catalog_csv}"
            )

    # --------------------------------------------------------
    # Direction selection
    # --------------------------------------------------------

    def make_directions(self):

        if self.explicit_directions:

            directions = read_directions_csv(
                self.explicit_directions
            )

            write_directions_csv(
                directions,
                self.directions_csv,
            )

            self.directions = directions
            self.select_operating_mode()
            self.write_direction_products()
            self.mark_direction_checkpoint()
            return

        if self.direction_checkpoint_is_valid():

            print("DD directions already exist.")
            self.directions = read_directions_csv(
                self.directions_csv
            )
            self.select_operating_mode()
            self.write_direction_products()
            return

        if self.checkpoint_exists("directions"):
            print(
                "WARNING: Direction inputs or thresholds changed; "
                "rebuilding DD directions."
            )
            (self.checkdir / "directions.done").unlink()

        if not self.catalog_csv.exists():

            message = (
                "No source catalogue is available for direction "
                "selection. Pass --catalog/--directions or install "
                "PyBDSF and run source extraction."
            )

            if self.strict:

                raise RuntimeError(message)

            print(f"WARNING: {message}")
            self.directions = []
            self.select_operating_mode()
            return

        sources = read_source_catalogue(
            self.catalog_csv
        )

        sources.sort(
            key=lambda item: item.flux_jy,
            reverse=True,
        )

        normalized_catalog = (
            self.catalogdir /
            "normalized_sources.csv"
        )

        write_source_catalogue(
            sources,
            normalized_catalog,
        )

        uncovered = []
        directions = cluster_sources_into_directions(
            sources,
            min_flux_jy=self.config["min_flux_jy"],
            max_directions=self.config["max_directions"],
            merge_radius_arcmin=self.config["merge_radius_arcmin"],
            uncovered=uncovered,
        )
        self.record_direction_selection(sources, directions, uncovered)

        if len(directions) == 0:

            message = (
                "No DD directions passed the flux threshold. "
                "Lower --min-flux-jy or supply --directions."
            )

            if self.strict:

                raise RuntimeError(message)

            print(f"WARNING: {message}")

        self.directions = directions
        self.select_operating_mode()
        self.write_direction_products()
        self.mark_direction_checkpoint()

    # --------------------------------------------------------

    def catalogue_flux_scale(self):
        """Report whether the input catalogue fluxes are apparent or intrinsic.

        ``min_flux_jy`` is a proxy for the signal available to the DD solver,
        which sees *apparent* flux.  A catalogue built from the PB-corrected
        image carries intrinsic flux, so an off-axis source at PB=0.35 enters
        direction selection looking about three times brighter than the solver
        will actually find it.
        """

        image = str(self.source_detection.get("input_image") or "")
        preference = (
            self.source_detection.get("parameters", {}) or {}
        ).get("image_preference")
        selection = self.source_detection.get("image_selection") or {}

        # Prefer the resolved role recorded by source detection over the
        # requested preference: the two can differ, and the filename is only
        # a guess.
        recorded = selection.get("is_primary_beam_corrected")
        if isinstance(recorded, bool):
            intrinsic = recorded
        else:
            intrinsic = (
                preference == "pbcor" or
                "pbcor" in Path(image).name.lower()
            )

        return {
            "catalogue_image": image or None,
            "image_preference": preference,
            "resolved_role": selection.get("resolved_role"),
            "flux_scale": (
                "primary-beam-corrected (intrinsic)"
                if intrinsic
                else "apparent"
            ),
            "min_flux_jy_is_apparent": not intrinsic,
            "note": (
                "min_flux_jy is compared against intrinsic flux; the solver "
                "sees the apparent flux, which is lower off axis by the "
                "primary-beam gain at each source position."
                if intrinsic
                else "min_flux_jy is compared against apparent flux."
            ),
        }

    # --------------------------------------------------------

    def record_direction_selection(self, sources, directions, uncovered):
        """Persist which bright sources did and did not get a direction."""

        threshold = float(self.config["min_flux_jy"])
        bright = [item for item in sources if item.flux_jy >= threshold]
        uncovered_csv = self.catalogdir / "uncovered_bright_sources.csv"

        if uncovered:
            write_source_catalogue(uncovered, uncovered_csv)
            print(
                f"WARNING: {len(uncovered)} source(s) above the "
                f"{threshold} Jy direction threshold received no DD "
                f"direction because max_directions="
                f"{self.config['max_directions']} was reached. They will not "
                "be direction-corrected and their artefacts remain in the "
                f"image. See {uncovered_csv}."
            )
        elif uncovered_csv.exists():
            uncovered_csv.unlink()

        self.direction_selection = {
            "min_flux_jy": threshold,
            "max_directions": int(self.config["max_directions"]),
            "merge_radius_arcmin": float(
                self.config["merge_radius_arcmin"]
            ),
            "catalogue_sources": len(sources),
            "sources_above_threshold": len(bright),
            "directions_selected": len(directions),
            "uncovered_bright_sources": len(uncovered),
            "uncovered_bright_sources_csv": (
                str(uncovered_csv) if uncovered else None
            ),
            "uncovered_bright_source_names": [
                item.name for item in uncovered
            ],
            "flux_scale": self.catalogue_flux_scale(),
        }
        save_json(
            self.direction_selection,
            self.directiondir / "direction_selection.json",
        )
        return self.direction_selection

    # --------------------------------------------------------

    def target_field_count(self):

        targets = (
            self.inspection.get("target_fields") or
            self.inspection.get("fields", {}).get("targets") or
            []
        )

        if isinstance(targets, str):
            targets = [targets]

        return len(set(str(item) for item in targets if item))

    # --------------------------------------------------------

    def select_operating_mode(self):
        """Choose and apply a conservative single-field or survey strategy."""

        valid = {"auto", "single", "survey"}

        if self.requested_mode not in valid:
            raise RuntimeError(
                f"Unsupported DD mode {self.requested_mode!r}; "
                f"choose from {sorted(valid)}."
            )

        if self.requested_mode == "auto":
            minimum = int(
                DD_TASK_CONFIG.get("survey_min_directions", 4)
            )
            multi_field = self.target_field_count() > 1
            self.mode = (
                "survey"
                if multi_field or len(self.directions) >= minimum
                else "single"
            )
        else:
            self.mode = self.requested_mode

        available = len(self.directions)

        if self.mode == "single":
            self.solver_direction_count = min(
                available,
                max(
                    1,
                    int(
                        DD_TASK_CONFIG.get(
                            "single_max_directions",
                            3,
                        )
                    ),
                ),
            )
        else:
            self.solver_direction_count = available

        if self.requested_engine == "auto":
            # DDFacet applies killMS DDE solutions during deconvolution.  This
            # is valid for both one/few-direction target fields and multi-facet
            # surveys, whereas CubiCal does not write generally DD-corrected
            # visibilities suitable for a later generic CASA imaging pass.
            self.engine = "killms-ddfacet"
        else:
            self.engine = self.requested_engine

        print(
            f"DD strategy      : {self.mode} "
            f"({self.solver_direction_count}/{available} directions, "
            f"{self.engine})"
        )

    # --------------------------------------------------------

    def write_direction_products(self):

        write_directions_csv(
            self.directions,
            self.directions_csv,
        )

        write_ds9_regions(
            self.directions,
            self.directions_region,
            radius_arcmin=2.0,
        )

        solver_count = getattr(
            self,
            "solver_direction_count",
            len(self.directions),
        )
        write_ds9_regions(
            self.directions[:solver_count],
            self.solver_directions_region,
            radius_arcmin=2.0,
        )

        write_ds9_regions(
            self.directions,
            self.facet_region,
            radius_arcmin=self.config["facet_radius_arcmin"],
        )

        write_tigger_lsm(
            self.directions,
            self.skymodel,
        )

        print(
            f"Selected {len(self.directions)} DD directions."
        )

    # --------------------------------------------------------
    # QA plots
    # --------------------------------------------------------

    def generate_qa_plots(self):
        """Render selected DD directions without requiring a display server."""

        directions = getattr(self, "directions", [])
        status = {"plots": [], "warnings": []}

        try:
            outfile = QA_PRODUCTS.plot_dd_directions(
                directions,
                self.plotdir / "dd_directions.png",
            )
        except Exception as exc:
            status["warnings"].append(
                f"Could not create DD direction plot: {exc}"
            )
        else:
            if outfile is None:
                status["warnings"].append(
                    "No usable DD directions are available for plotting."
                )
            else:
                status["plots"].append(str(outfile))

        status_path = self.root / "qa_plot_status.json"
        QA_PRODUCTS.write_status(status_path, status)

        if status["plots"]:
            print(
                "DD-calibration QA plots written: " + ", ".join(status["plots"]),
                flush=True,
            )
        else:
            print(
                "WARNING: No DD-calibration QA plot was created. See "
                f"{status_path}.",
                flush=True,
            )

        return status

    # --------------------------------------------------------
    # Solver command generation
    # --------------------------------------------------------

    def imaging_parameters(self):

        parameters = self.imaging.get("imaging_parameters", {})

        return {
            "imsize": int(
                parameters.get(
                    "imsize",
                    DD_TASK_CONFIG.get("ddfacet_imsize", 8192),
                )
            ),
            "cell_arcsec": parse_quantity(
                parameters.get(
                    "cell",
                    DD_TASK_CONFIG.get("ddfacet_cell", "1arcsec"),
                ),
                "arcsec",
            ),
            "robust": float(
                parameters.get(
                    "robust",
                    DD_TASK_CONFIG.get("ddfacet_robust", 0.5),
                )
            ),
        }

    # --------------------------------------------------------

    def frequency_span_mhz(self):

        spws = self.inspection.get("spws", {})
        lows = []
        highs = []

        for item in spws.values():
            try:
                lows.append(float(item["fmin_MHz"]))
                highs.append(float(item["fmax_MHz"]))
            except (KeyError, TypeError, ValueError):
                continue

        if not lows or not highs:
            return None

        return max(highs) - min(lows)

    # --------------------------------------------------------

    def killms_nchan_solutions(self):

        configured = DD_TASK_CONFIG.get("killms_nchan_solutions")

        if configured is not None:
            return max(1, int(configured))

        bandwidth = self.frequency_span_mhz()

        if not bandwidth:
            return 1

        interval = parse_quantity(
            self.config["solint_freq"],
            "mhz",
        )
        maximum = max(
            1,
            int(DD_TASK_CONFIG.get("killms_max_nchan_solutions", 16)),
        )

        return min(
            maximum,
            max(1, int(math.ceil(bandwidth / interval))),
        )

    # --------------------------------------------------------

    def solver_paths(self):

        target = clean_name(self.target)
        full_signature = self.solver_signature()
        signature = full_signature[:10]
        mode = self.mode or self.requested_mode
        stem = f"{target}_{mode}_{signature}"
        di_base = self.imagedir / f"{stem}_DI"
        dd_base = self.imagedir / f"{stem}_DD"
        # killMS materializes ``killMS.<OutSolsName>.sols.npz`` inside the
        # input MS.  Use a reserved full-signature namespace so the only
        # non-table file ignored by the generic MS stat fingerprint is exact,
        # pipeline-owned, and separately content-hashed by DD provenance.
        solution_name = f"radio_continuum_dd_{full_signature}"

        return {
            "mslist": self.configdir / "mslist.txt",
            "di_base": di_base,
            "dd_base": dd_base,
            "dico_model": Path(f"{di_base}.DicoModel"),
            "cluster_nodes": Path(
                f"{self.solver_directions_region}.ClusterCat.npy"
            ),
            "solution_name": solution_name,
            "solution_file": (
                Path(self.ms) /
                f"killMS.{solution_name}.sols.npz"
            ),
            "cache_dir": self.root / "cache",
        }

    # --------------------------------------------------------

    def solver_signature(self):

        payload = {
            "schema": 3,
            "measurement_set": str(Path(self.ms).resolve()),
            "measurement_set_stat_sha256": (
                self.measurement_set_stat_identity()
            ),
            "source_catalog_content_identity": self.file_fingerprint(
                self.catalog_csv
            ),
            "target": self.target,
            "band": self.band,
            "mode": self.mode or self.requested_mode,
            "engine": self.engine,
            "data_column": self.data_column,
            "configuration": self.config,
            "solver_configuration": {
                key: DD_TASK_CONFIG.get(key)
                for key in [
                    "killms_pol_mode",
                    "killms_weighting",
                    "killms_nchan_solutions",
                    "killms_max_nchan_solutions",
                    "ddfacet_major_cycles",
                    "ddfacet_rms_factor",
                    "ddfacet_mask_sigma",
                    "ddfacet_min_facet_diameter_deg",
                    "ddfacet_global_norm",
                ]
            },
            "templates": {
                "solver": self.solver_command_template,
                "apply": self.apply_command_template,
                "killms": self.killms_command_template,
                "ddfacet_initial": self.ddfacet_command_template,
                "make_model": DD_TASK_CONFIG.get(
                    "make_model_command_template"
                ),
                "ddfacet_final": DD_TASK_CONFIG.get(
                    "ddfacet_final_command_template"
                ),
            },
            "derived_nchan_solutions": (
                self.killms_nchan_solutions()
            ),
            "imaging": self.imaging_parameters(),
            "directions": [
                item.to_dict()
                for item in getattr(self, "directions", [])
            ],
            "solver_direction_count": getattr(
                self,
                "solver_direction_count",
                0,
            ),
        }
        encoded = json.dumps(
            to_builtin(payload),
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")

        return hashlib.sha256(encoded).hexdigest()

    # --------------------------------------------------------

    def template_context(self):

        paths = self.solver_paths()
        imaging = self.imaging_parameters()
        facet_diameter_deg = max(
            0.1,
            2.0 * float(self.config["facet_radius_arcmin"]) / 60.0,
        )

        return {
            "ms": self.ms,
            "target": clean_name(self.target),
            "band": self.band,
            "engine": self.engine,
            "mode": self.mode,
            "data_column": self.data_column,
            "model_column": self.model_column,
            "directions": str(self.directions_csv),
            "directions_reg": str(self.directions_region),
            "solver_directions_reg": str(
                self.solver_directions_region
            ),
            "facets": str(self.facet_region),
            "skymodel": str(self.skymodel),
            "solutions_dir": str(self.soldir),
            "output_dir": str(self.root),
            "config_dir": str(self.configdir),
            "log_dir": str(self.logdir),
            "image_dir": str(self.imagedir),
            "ncpu": self.ncpu,
            "solint_time": self.config["solint_time"],
            "solint_freq": self.config["solint_freq"],
            "solint_time_min": parse_quantity(
                self.config["solint_time"],
                "min",
            ),
            "nchan_solutions": self.killms_nchan_solutions(),
            "n_directions": max(
                1,
                getattr(
                    self,
                    "solver_direction_count",
                    len(self.directions),
                ),
            ),
            "n_facets": max(
                3,
                getattr(
                    self,
                    "solver_direction_count",
                    len(self.directions),
                ) + 1,
            ),
            "imsize": imaging["imsize"],
            "cell_arcsec": imaging["cell_arcsec"],
            "robust": imaging["robust"],
            "facet_diameter_deg": facet_diameter_deg,
            "mslist": str(paths["mslist"]),
            "di_base": str(paths["di_base"]),
            "dd_base": str(paths["dd_base"]),
            "dico_model": str(paths["dico_model"]),
            "cluster_nodes": str(paths["cluster_nodes"]),
            "solution_name": paths["solution_name"],
            "solution_file": str(paths["solution_file"]),
            "cache_dir": str(paths["cache_dir"]),
            "solver_signature": self.solver_signature(),
            "out_name": str(
                self.soldir /
                f"{clean_name(self.target)}_dd"
            ),
        }

    # --------------------------------------------------------

    def format_template(self, template):

        return template.format(
            **self.template_context()
        )

    # --------------------------------------------------------

    def write_ms_list(self):

        paths = self.solver_paths()
        paths["mslist"].write_text(
            str(Path(self.ms).resolve()) + "\n"
        )
        paths["cache_dir"].mkdir(parents=True, exist_ok=True)

        return paths["mslist"]

    # --------------------------------------------------------

    def ddfacet_common_command(self, output_name):

        exe = shutil.which("DDF.py") or "DDF.py"
        context = self.template_context()

        return [
            exe,
            f"--Output-Name={output_name}",
            f"--Data-MS={context['mslist']}",
            f"--Data-ColName={self.data_column}",
            f"--Parallel-NCPU={self.ncpu}",
            "--Image-Mode=Clean",
            "--Deconv-Mode=SSD",
            "--Deconv-PeakFactor=0.001",
            f"--Deconv-MaxMajorIter={int(DD_TASK_CONFIG.get('ddfacet_major_cycles', 5))}",
            f"--Deconv-RMSFactor={float(DD_TASK_CONFIG.get('ddfacet_rms_factor', 3.0))}",
            f"--Weight-Robust={context['robust']}",
            f"--Image-NPix={context['imsize']}",
            f"--Image-Cell={context['cell_arcsec']}",
            f"--Facets-NFacets={context['n_facets']}",
            "--Mask-Auto=1",
            f"--Mask-SigTh={float(DD_TASK_CONFIG.get('ddfacet_mask_sigma', 5.0))}",
            f"--Cache-Dir={context['cache_dir']}",
            "--Cache-Reset=0",
        ]

    # --------------------------------------------------------

    def default_ddfacet_initial_command(self):

        return self.ddfacet_common_command(
            self.solver_paths()["di_base"]
        )

    # --------------------------------------------------------

    def default_make_model_command(self):

        exe = shutil.which("MakeModel.py") or "MakeModel.py"

        return [
            exe,
            "--BaseImageName",
            str(self.solver_paths()["di_base"]),
            "--NCluster",
            "0",
            "--ds9PreClusterFile",
            str(self.solver_directions_region),
        ]

    # --------------------------------------------------------

    def default_killms_command(self):

        exe = shutil.which("kMS.py") or "kMS.py"
        context = self.template_context()

        return [
            exe,
            "--MSName",
            context["mslist"],
            "--SolverType",
            "KAFCA",
            "--PolMode",
            str(DD_TASK_CONFIG.get("killms_pol_mode", "Scalar")),
            "--BaseImageName",
            context["di_base"],
            "--dt",
            str(context["solint_time_min"]),
            "--NCPU",
            str(self.ncpu),
            "--OutSolsName",
            context["solution_name"],
            "--NChanSols",
            str(context["nchan_solutions"]),
            "--InCol",
            self.data_column,
            "--OutCol",
            self.data_column,
            "--Weighting",
            str(DD_TASK_CONFIG.get("killms_weighting", "Natural")),
            "--NodesFile",
            context["cluster_nodes"],
            "--MaxFacetSize",
            str(context["facet_diameter_deg"]),
        ]

    # --------------------------------------------------------

    def default_ddfacet_final_command(self):

        paths = self.solver_paths()
        command = self.ddfacet_common_command(paths["dd_base"])

        command.extend(
            [
                f"--DDESolutions-DDSols={paths['solution_name']}",
                (
                    "--DDESolutions-GlobalNorm=" +
                    str(
                        DD_TASK_CONFIG.get(
                            "ddfacet_global_norm",
                            "MeanAbs",
                        )
                    )
                ),
                (
                    "--DDESolutions-JonesMode=" +
                    str(
                        DD_TASK_CONFIG.get(
                            "killms_pol_mode",
                            "Scalar",
                        )
                    )
                ),
                f"--Predict-InitDicoModel={paths['dico_model']}",
                (
                    "--Facets-DiamMax=" +
                    str(self.template_context()["facet_diameter_deg"])
                ),
                (
                    "--Facets-DiamMin=" +
                    str(
                        float(
                            DD_TASK_CONFIG.get(
                                "ddfacet_min_facet_diameter_deg",
                                0.05,
                            )
                        )
                    )
                ),
            ]
        )

        return command

    # --------------------------------------------------------

    def write_solver_config(self):

        context = self.template_context()

        config = {
            "schema_version": 2,
            "created_utc": now_utc(),
            "status": "prepared" if not self.execute else "ready",
            "solver_signature": self.solver_signature(),
            "measurement_set_stat_sha256": (
                self.measurement_set_stat_identity()
            ),
            "source_catalog_content_identity": self.file_fingerprint(
                self.catalog_csv
            ),
            "mode": self.mode,
            "requested_mode": self.requested_mode,
            "engine": self.engine,
            "requested_engine": self.requested_engine,
            "measurement_set": self.ms,
            "target": self.target,
            "band": self.band,
            "data_column": self.data_column,
            "model_column": self.model_column,
            "directions": str(self.directions_csv),
            "skymodel": str(self.skymodel),
            "solutions_dir": str(self.soldir),
            "ncpu": self.ncpu,
            "solint_time": self.config["solint_time"],
            "solint_freq": self.config["solint_freq"],
            "template_context": context,
            "notes": [
                "The built-in killMS/DDFacet sequence follows the upstream "
                "DI image -> MakeModel clustering -> KAFCA solve -> "
                "DDESolutions DDFacet image workflow.",
                "CubiCal execution requires site-supplied templates and a "
                "direction-tagged sky model. CubiCal cannot generally write "
                "a visibility column corrected by direction-dependent gains.",
                "The direction products in Products/ddcal/directions "
                "are the stable interface between source extraction "
                "and the DD solver.",
            ],
        }

        outfile = self.configdir / "ddcal_config.json"

        save_json(config, outfile)

        return outfile

    # --------------------------------------------------------

    def solver_commands(self):

        if self.engine == "prepare":

            return []

        if self.engine == "cubical":

            if not self.solver_command_template:
                return []

            commands = [
                {
                    "name": "cubical_dd_workflow",
                    "command": self.format_template(
                        self.solver_command_template
                    ),
                    "required": True,
                    "expected": [],
                },
            ]

            if self.apply_command_template:
                commands.append(
                    {
                        "name": "cubical_dd_image_or_apply",
                        "command": self.format_template(
                            self.apply_command_template
                        ),
                        "required": True,
                        "expected": [],
                    }
                )

            return commands

        if self.engine == "killms-ddfacet":

            paths = self.solver_paths()
            initial = (
                self.format_template(
                    self.ddfacet_command_template
                )
                if self.ddfacet_command_template
                else self.default_ddfacet_initial_command()
            )
            make_model_template = DD_TASK_CONFIG.get(
                "make_model_command_template"
            )
            make_model = (
                self.format_template(make_model_template)
                if make_model_template
                else self.default_make_model_command()
            )
            kms = (
                self.format_template(
                    self.killms_command_template
                )
                if self.killms_command_template
                else self.default_killms_command()
            )
            final_template = DD_TASK_CONFIG.get(
                "ddfacet_final_command_template"
            )
            final = (
                self.format_template(final_template)
                if final_template
                else self.default_ddfacet_final_command()
            )

            return [
                {
                    "name": "ddfacet_initial",
                    "command": initial,
                    "required": True,
                    "expected": [str(paths["dico_model"])],
                },
                {
                    "name": "make_model_clusters",
                    "command": make_model,
                    "required": True,
                    "expected": [str(paths["cluster_nodes"])],
                },
                {
                    "name": "killms_solve",
                    "command": kms,
                    "required": True,
                    "expected": [str(paths["solution_file"])],
                },
                {
                    "name": "ddfacet_dd_image",
                    "command": final,
                    "required": True,
                    "expected": [],
                },
            ]

        raise RuntimeError(
            f"Unsupported DD engine: {self.engine}"
        )

    # --------------------------------------------------------

    def write_command_script(self):

        commands = self.solver_commands()

        with open(self.command_script, "w") as f:

            f.write("#!/bin/sh\n")
            f.write("set -e\n\n")
            f.write("# Auto-generated by ddcal.py\n")
            f.write(f"# Created: {now_utc()}\n\n")

            for spec in commands:

                f.write(f"echo 'Running {spec['name']}'\n")
                f.write(
                    command_to_string(
                        split_command(spec["command"])
                    ) + "\n\n"
                )

        os.chmod(self.command_script, 0o755)

        return self.command_script

    # --------------------------------------------------------

    def solver_step_products(self, spec):

        products = [
            str(Path(item).expanduser().resolve())
            for item in spec.get("expected", [])
            if Path(item).exists()
        ]
        if spec.get("name") == "ddfacet_dd_image":
            products.extend(self.find_dd_image_products())
        return list(dict.fromkeys(products))

    # --------------------------------------------------------

    def current_solver_products(self, commands=None):

        commands = self.solver_commands() if commands is None else commands
        products = []
        for spec in commands:
            products.extend(self.solver_step_products(spec))
        return list(dict.fromkeys(products))

    # --------------------------------------------------------

    def current_post_dd_products(self):

        products = list(self.current_solver_products())
        expected_pbcor = (
            self.imagedir /
            (
                f"{clean_name(self.target)}_{self.mode}_"
                f"{self.solver_signature()[:10]}_DD.pbcor.fits"
            )
        )
        for path in [
            expected_pbcor,
            self.reportdir / "dd_image_metrics.json",
        ]:
            if path.exists():
                products.append(str(path.resolve()))
        return list(dict.fromkeys(products))

    # --------------------------------------------------------

    def run_solver(self):

        if self.skip_solver:

            print("Skipping DD solver by request.")
            return

        if len(getattr(self, "directions", [])) == 0:

            message = (
                "No DD directions are available. Solver stage "
                "will not run."
            )

            if self.strict:

                raise RuntimeError(message)

            print(f"WARNING: {message}")
            return

        if self.execute:

            self.validate_execution_configuration()

        if self.engine == "killms-ddfacet":
            self.write_ms_list()

        self.write_solver_config()
        self.write_command_script()

        commands = self.solver_commands()
        solver_products = self.current_solver_products(commands)

        # A marker without matching input/output provenance is never trusted,
        # including markers written by an older pipeline schema.  Preserve
        # the signature-scoped outputs and restart the complete dependent
        # sequence so no later solve can consume a stale DI model or gain.
        for spec in commands:
            checkpoint = f"solver_{spec['name']}"
            if not self.checkpoint_exists(checkpoint):
                continue
            valid, reason = self.bound_checkpoint_status(
                checkpoint,
                self.solver_step_products(spec),
            )
            if not valid:
                self.invalidate_solver_state(
                    f"{checkpoint}: {reason}",
                    self.current_post_dd_products(),
                )
                solver_products = []
                break

        if self.checkpoint_exists("solver"):
            solver_products = self.current_solver_products(commands)
            valid, reason = self.bound_checkpoint_status(
                "solver",
                solver_products,
            )
            if valid:
                if self.engine == "killms-ddfacet":
                    self.validate_killms_solutions(
                        self.solver_paths()["solution_file"]
                    )
                    self.dd_image_products = self.find_dd_image_products()
                print("DD solver checkpoint and products are valid.")
                return
            self.invalidate_solver_state(
                f"solver: {reason}",
                self.current_post_dd_products(),
            )

        if not commands:
            message = (
                f"No executable {self.engine} DD workflow is configured."
            )
            if self.execute:
                raise RuntimeError(message)
            print(f"WARNING: {message}")
            return

        for spec in commands:

            checkpoint = f"solver_{spec['name']}"
            reusable, _reuse_reason = self.bound_checkpoint_status(
                checkpoint,
                self.solver_step_products(spec),
            )

            if reusable:
                print(
                    f"Reusing completed DD step: {spec['name']}"
                )
                self.command_records.append(
                    {
                        "name": spec["name"],
                        "command": command_to_string(
                            split_command(spec["command"])
                        ),
                        "executed": False,
                        "returncode": 0,
                        "status": "reused",
                    }
                )
                if (
                    spec["name"] == "killms_solve"
                ):
                    self.validate_killms_solutions(
                        self.solver_paths()["solution_file"]
                    )
                continue

            record = self.run_command(
                spec["name"],
                spec["command"],
                required=spec["required"],
            )

            if not self.execute:
                continue

            if record["status"] != "ok":
                raise RuntimeError(
                    f"DD solver step did not complete: {spec['name']}"
                )

            if not self.expected_products_exist(spec["expected"]):
                raise RuntimeError(
                    f"DD solver step {spec['name']} returned success but "
                    "did not create its required output products."
                )

            if (
                spec["name"] == "killms_solve" and
                self.engine == "killms-ddfacet"
            ):
                self.validate_killms_solutions(
                    self.solver_paths()["solution_file"]
                )

            if spec["name"] == "ddfacet_dd_image":
                self.dd_image_products = self.find_dd_image_products()
                if not self.dd_image_products:
                    raise RuntimeError(
                        "DDFacet returned success but no final DD restored "
                        "FITS image was found."
                    )

            self.mark_bound_checkpoint(
                checkpoint,
                self.solver_step_products(spec),
            )

        if self.execute:
            self.mark_bound_checkpoint(
                "solver",
                self.current_solver_products(commands),
            )

    # --------------------------------------------------------

    def expected_products_exist(self, products):

        return all(Path(item).exists() for item in products)

    # --------------------------------------------------------

    def numeric_solution_arrays(self, value, name="root"):

        try:
            import numpy as np
        except ImportError as exc:
            raise RuntimeError(
                "NumPy is required to validate killMS solutions."
            ) from exc

        if isinstance(value, dict):
            for key, item in value.items():
                yield from self.numeric_solution_arrays(
                    item,
                    f"{name}.{key}",
                )
            return

        if not isinstance(value, np.ndarray):
            return

        if value.dtype.names:
            for field in value.dtype.names:
                yield from self.numeric_solution_arrays(
                    value[field],
                    f"{name}.{field}",
                )
            return

        if value.dtype == object:
            for index, item in enumerate(value.flat):
                if index >= 10000:
                    break
                if isinstance(item, np.ndarray):
                    yield from self.numeric_solution_arrays(
                        item,
                        f"{name}[{index}]",
                    )
                elif isinstance(item, dict):
                    yield from self.numeric_solution_arrays(
                        item,
                        f"{name}[{index}]",
                    )
            return

        if np.issubdtype(value.dtype, np.number) and value.size:
            yield name, value

    # --------------------------------------------------------

    def validate_killms_solutions(self, solution_file):

        import numpy as np

        solution_file = Path(solution_file)
        report_path = self.soldir / "killms_solution_qa.json"
        plot_path = self.plotdir / "killms_gain_amplitudes.png"

        with np.load(solution_file, allow_pickle=True) as archive:
            arrays = []
            for key in archive.files:
                arrays.extend(
                    self.numeric_solution_arrays(
                        archive[key],
                        key,
                    )
                )

        gain_arrays = [
            (name, array)
            for name, array in arrays
            if np.iscomplexobj(array) or
            "gain" in name.lower() or
            name.lower().endswith(".g")
        ]

        if not gain_arrays:
            raise RuntimeError(
                "The killMS solution archive contains no recognizable gain "
                "arrays; refusing to accept an unvalidated solution."
            )

        samples = []
        array_stats = []
        total = 0
        usable = 0

        for name, array in gain_arrays:
            flat = np.asarray(array).reshape(-1)
            if flat.size > 200000:
                step = max(1, flat.size // 200000)
                flat = flat[::step]

            amplitude = np.abs(flat.astype(np.complex128, copy=False))
            finite = np.isfinite(amplitude)
            nonzero = amplitude > 0
            valid = finite & nonzero
            total += int(amplitude.size)
            usable += int(np.count_nonzero(valid))
            good = amplitude[valid]

            if good.size:
                samples.append(good)
                array_stats.append(
                    {
                        "name": name,
                        "samples": int(amplitude.size),
                        "usable": int(good.size),
                        "median_amplitude": float(np.median(good)),
                        "p05_amplitude": float(np.percentile(good, 5)),
                        "p95_amplitude": float(np.percentile(good, 95)),
                    }
                )

        usable_fraction = float(usable / total) if total else 0.0
        minimum = float(
            DD_TASK_CONFIG.get(
                "killms_min_usable_solution_fraction",
                0.5,
            )
        )

        if not samples or usable_fraction < minimum:
            raise RuntimeError(
                "killMS gain validation failed: usable finite/non-zero "
                f"solution fraction {usable_fraction:.3f} is below "
                f"{minimum:.3f}."
            )

        combined = np.concatenate(samples)
        median = float(np.median(combined))
        warnings = []

        if median < 0.05 or median > 20.0:
            warnings.append(
                "Median gain amplitude is far from unity; inspect flux-scale "
                "degeneracy and the sky model before using the DD image."
            )

        qa = {
            "created_utc": now_utc(),
            "solution_file": str(solution_file),
            "usable_fraction": usable_fraction,
            "minimum_usable_fraction": minimum,
            "median_amplitude": median,
            "p05_amplitude": float(np.percentile(combined, 5)),
            "p95_amplitude": float(np.percentile(combined, 95)),
            "arrays": array_stats,
            "warnings": warnings,
        }
        save_json(qa, report_path)

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            positive = combined[combined > 0]
            figure, axis = plt.subplots(figsize=(7.0, 4.5))
            axis.hist(
                np.log10(positive),
                bins=80,
                color="#315b8a",
                alpha=0.85,
            )
            axis.axvline(
                math.log10(median),
                color="#b22222",
                linestyle="--",
                label=f"median={median:.3g}",
            )
            axis.set_xlabel("log10 |DD gain|")
            axis.set_ylabel("Solution samples")
            axis.set_title("killMS direction-dependent gain amplitudes")
            axis.legend()
            figure.tight_layout()
            figure.savefig(plot_path, dpi=160)
            plt.close(figure)
        except Exception as exc:
            qa["warnings"].append(
                f"Could not render killMS gain QA plot: {exc}"
            )
            save_json(qa, report_path)

        qa_products = [solution_file, report_path]
        if plot_path.is_file():
            qa_products.append(plot_path)
        self.mark_bound_checkpoint("solution_qa", qa_products)

        return qa

    # --------------------------------------------------------

    def find_dd_image_products(self):

        base = self.solver_paths()["dd_base"].name
        candidates = sorted(
            self.imagedir.glob(f"{base}*.fits")
        )
        candidates = [
            item
            for item in candidates
            if ".pbcor." not in item.name.lower()
        ]
        restored = [
            item
            for item in candidates
            if "restor" in item.name.lower()
        ]

        restored.sort(
            key=lambda item: (
                "app.restored" not in item.name.lower(),
                item.name,
            )
        )

        return [
            str(item)
            for item in (restored or candidates)
            if item.is_file()
        ]

    # --------------------------------------------------------

    def validate_execution_configuration(self):

        """Validate solver dependencies before starting a costly DD run."""

        self.validate_ms_data_column()

        if self.engine == "prepare":

            raise RuntimeError(
                "DD engine 'prepare' cannot execute a solver. Select a "
                "solver engine and provide validated command templates."
            )

        required_templates = []

        if self.engine == "cubical":

            if not self.solver_command_template:
                raise RuntimeError(
                    "CubiCal DD execution needs a site-validated workflow "
                    "template with a direction-tagged Tigger model. CubiCal "
                    "does not generally produce a visibility column corrected "
                    "for DD gains. Use --engine killms-ddfacet for the built-in "
                    "end-to-end workflow, or provide "
                    "--solver-command-template."
                )

        elif self.engine == "killms-ddfacet":

            tools = {
                "DDF.py": shutil.which("DDF.py"),
                "MakeModel.py": shutil.which("MakeModel.py"),
                "kMS.py": shutil.which("kMS.py"),
            }
            missing = [
                name
                for name, path in tools.items()
                if path is None
            ]

            if missing:
                raise RuntimeError(
                    "killMS/DDFacet execution requires these commands in "
                    "PATH: " + ", ".join(missing) + ". On ilifu, load the "
                    "site module or activate the solver environment from "
                    "slurm_module_commands before starting the pipeline."
                )

            if bool(
                DD_TASK_CONFIG.get("validate_solver_cli", True)
            ):
                self.validate_solver_help(tools)

    # --------------------------------------------------------

    def validate_ms_data_column(self):

        try:
            from casacore.tables import table as casacore_table
        except ImportError:
            print(
                "WARNING: python-casacore is unavailable to ddcal.py; "
                f"the {self.data_column} column will be validated by the "
                "external solver."
            )
            return

        handle = None
        try:
            handle = casacore_table(
                self.ms,
                readonly=True,
                ack=False,
            )
            columns = set(handle.colnames())
        finally:
            if handle is not None:
                handle.close()

        if self.data_column not in columns:
            raise RuntimeError(
                f"DD input column {self.data_column!r} is absent from "
                f"{self.ms}. Available columns include: "
                + ", ".join(sorted(columns))
            )

    # --------------------------------------------------------

    def validate_solver_help(self, tools):

        checks = {
            "DDF.py": [
                "Output-Name",
                "Data-MS",
                "Data-ColName",
                "Parallel-NCPU",
                "Image-Mode",
                "Deconv-Mode",
                "Deconv-PeakFactor",
                "Deconv-MaxMajorIter",
                "Deconv-RMSFactor",
                "Weight-Robust",
                "Image-NPix",
                "Image-Cell",
                "Facets-NFacets",
                "Mask-Auto",
                "Mask-SigTh",
                "Cache-Dir",
                "Cache-Reset",
                "DDESolutions-DDSols",
                "DDESolutions-GlobalNorm",
                "DDESolutions-JonesMode",
                "Predict-InitDicoModel",
                "Facets-DiamMax",
                "Facets-DiamMin",
            ],
            "MakeModel.py": [
                "BaseImageName",
                "NCluster",
                "ds9PreClusterFile",
            ],
            "kMS.py": [
                "MSName",
                "SolverType",
                "PolMode",
                "BaseImageName",
                "dt",
                "NCPU",
                "OutSolsName",
                "NChanSols",
                "InCol",
                "OutCol",
                "Weighting",
                "NodesFile",
                "MaxFacetSize",
            ],
        }

        for name, tokens in checks.items():
            try:
                result = subprocess.run(
                    [tools[name], "-h"],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=30,
                )
            except Exception as exc:
                raise RuntimeError(
                    f"Could not inspect {name} command-line options: {exc}"
                ) from exc

            output = result.stdout or ""
            absent = [
                token
                for token in tokens
                if token.lower() not in output.lower()
            ]

            if absent:
                raise RuntimeError(
                    f"The installed {name} does not advertise required "
                    f"options: {', '.join(absent)}. Use a compatible "
                    "killMS/DDFacet environment or site command templates."
                )

    # --------------------------------------------------------
    # Post-DD imaging
    # --------------------------------------------------------

    MEASURED_PB_SUFFIXES = (".ugmrt.pb.tt0", ".katbeam.pb.tt0")

    def select_primary_beam_image(self):
        """Return the beam final imaging actually divided by.

        The DD product and the main continuum image must sit on one flux
        scale.  CASA's generic tclean ``.pb`` is not a validated uGMRT beam
        and is not what imaging.py corrected with, so selecting it here would
        publish two images with different off-axis scales.  The previous
        implementation matched on a ``.pb.tt0`` suffix, which the measured
        beam also ends with, and returned whichever came first in the
        recorded product inventory -- always CASA's.
        """

        products = self.imaging.get("products", {})
        parameters = self.imaging.get("imaging_parameters", {}) or {}
        status = products.get("pbcor_status", {}) or {}

        # The imaging stage records the measured beam it built, by provider.
        for key in ("ugmrt_pb_image", "katbeam_pb_image"):
            recorded = status.get(key)
            if recorded and Path(str(recorded)).exists():
                return Path(str(recorded))

        for item in products.get("casa_images", []):
            path = Path(str(item))
            if path.name.endswith(self.MEASURED_PB_SUFFIXES) and path.exists():
                return path

        imagename = parameters.get("imagename")
        if imagename:
            for suffix in self.MEASURED_PB_SUFFIXES:
                candidate = Path(str(imagename) + suffix)
                if candidate.exists():
                    return candidate

        # CASA's own beam is acceptable only where it is also the model that
        # corrected the continuum image; otherwise fail closed.
        provider = parameters.get("primary_beam_provider")
        method = str(products.get("pbcor_method") or "")
        if provider == "casa" or method.startswith(
            ("impbcor_tt0", "widebandpbcor")
        ):
            candidates = [
                Path(str(item))
                for item in products.get("casa_images", [])
                if Path(str(item)).name.endswith((".pb.tt0", ".pb"))
            ]
            if imagename:
                candidates.extend([
                    Path(str(imagename) + ".pb.tt0"),
                    Path(str(imagename) + ".pb"),
                ])
            for candidate in candidates:
                if candidate.exists():
                    return candidate

        return None

    # --------------------------------------------------------

    def primary_beam_provider(self, pb_image):
        """Name the provider of the beam selected for DD correction."""

        name = Path(str(pb_image)).name
        if name.endswith(".ugmrt.pb.tt0"):
            return "ugmrt_polynomial"
        if name.endswith(".katbeam.pb.tt0"):
            return "katbeam"
        return "casa_tclean_pb"

    # --------------------------------------------------------

    def write_dd_pbcor_script(self, dd_fits, pb_image):

        token = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")
        target = clean_name(self.target)
        signature = self.solver_signature()[:10]
        stem = f"{target}_{self.mode}_{signature}_DD"
        imported = self.imagedir / f".{target}_dd_import_{token}.image"
        regridded_pb = self.imagedir / f".{target}_pb_regrid_{token}.image"
        corrected = self.imagedir / f".{target}_dd_pbcor_{token}.image"
        output_fits = self.imagedir / f"{stem}.pbcor.fits"
        metrics_path = self.reportdir / "dd_image_metrics.json"
        script = self.configdir / "run_dd_primary_beam_correction.py"
        cutoff = float(
            self.imaging.get(
                "imaging_parameters",
                {},
            ).get(
                "pbmin",
                DD_TASK_CONFIG.get("pbmin", 0.3),
            )
        )

        provider = self.primary_beam_provider(pb_image)

        lines = [
            "import json",
            "from casatasks import exportfits, importfits, impbcor, imregrid",
            "from casatasks import imstat",
            f"dd_fits = {str(dd_fits)!r}",
            f"pb_image = {str(pb_image)!r}",
            f"pb_provider = {provider!r}",
            f"imported = {str(imported)!r}",
            f"regridded_pb = {str(regridded_pb)!r}",
            f"corrected = {str(corrected)!r}",
            f"output_fits = {str(output_fits)!r}",
            f"metrics_path = {str(metrics_path)!r}",
            "importfits(",
            "    fitsimage=dd_fits,",
            "    imagename=imported,",
            "    overwrite=True,",
            ")",
            "imregrid(",
            "    imagename=pb_image,",
            "    template=imported,",
            "    output=regridded_pb,",
            "    overwrite=True,",
            ")",
            "impbcor(",
            "    imagename=imported,",
            "    pbimage=regridded_pb,",
            "    outfile=corrected,",
            "    mode='divide',",
            f"    cutoff={cutoff!r},",
            "    overwrite=True,",
            ")",
            "exportfits(",
            "    imagename=corrected,",
            "    fitsimage=output_fits,",
            "    overwrite=True,",
            "    dropdeg=True,",
            ")",
            "stats = imstat(",
            "    imagename=corrected,",
            "    algorithm='chauvenet',",
            "    maxiter=5,",
            ")",
            # Clipping removes the sources before it reports max and min,
            # so the extremes come from the classic algorithm.
            "extremes = imstat(imagename=corrected)",
            "rms = float(stats['rms'][0])",
            "peak = float(extremes['max'][0])",
            "minimum = float(extremes['min'][0])",
            # Noise measured on a PB-corrected image is dominated by the
            # amplified field edge and is not the map sensitivity. Measure the
            # apparent (flat-noise) image too; that is the number to quote.
            "apparent = imstat(",
            "    imagename=imported,",
            "    algorithm='chauvenet',",
            "    maxiter=5,",
            ")",
            "apparent_rms = float(apparent['rms'][0])",
            "apparent_peak = float(imstat(imagename=imported)['max'][0])",
            "metrics = {",
            "    'image': output_fits,",
            "    'peak_Jy_beam': peak,",
            "    'min_Jy_beam': minimum,",
            "    'rms_Jy_beam': rms,",
            "    'dynamic_range': peak / rms if rms > 0 else None,",
            "    'apparent_image': dd_fits,",
            "    'apparent_rms_Jy_beam': apparent_rms,",
            "    'apparent_peak_Jy_beam': apparent_peak,",
            "    'apparent_dynamic_range': (",
            "        apparent_peak / apparent_rms if apparent_rms > 0 else None",
            "    ),",
            "    'rms_note': (",
            "        'Quote apparent_rms_Jy_beam as the map sensitivity; "
            "rms_Jy_beam is inflated by primary-beam amplification at the "
            "field edge.'",
            "    ),",
            "    'pb_image': pb_image,",
            "    'pb_provider': pb_provider,",
            f"    'pb_cutoff': {cutoff!r},",
            "}",
            "with open(metrics_path, 'w') as handle:",
            "    json.dump(metrics, handle, indent=4)",
            "    handle.write('\\n')",
        ]

        script.write_text("\n".join(lines) + "\n")

        return script, output_fits, metrics_path

    # --------------------------------------------------------

    def run_post_dd_imaging(self):

        if self.engine != "killms-ddfacet":
            if not self.skip_imaging:
                print(
                    "WARNING: Generic CASA post-imaging is disabled for "
                    f"{self.engine}: no generally valid DD-corrected "
                    "visibility column is available."
                )
            return

        if self.skip_imaging:

            print(
                "Skipping DD primary-beam correction by request. "
                "The DDFacet restored image remains available."
            )
            return

        if not self.execute:

            print(
                "DD primary-beam correction will run after solver execution."
            )
            return

        if not self.checkpoint_exists("solver"):
            raise RuntimeError(
                "Cannot primary-beam-correct the DD image before the solver "
                "and DDFacet imaging sequence has completed."
            )

        self.dd_image_products = self.find_dd_image_products()
        if not self.dd_image_products:
            raise RuntimeError(
                "No final DDFacet restored FITS image is available for "
                "primary-beam correction."
            )

        pb_image = self.select_primary_beam_image()
        if pb_image is None:
            raise RuntimeError(
                "No validated primary-beam image is available for DD "
                "correction. The DD product must use the same measured beam "
                "that final imaging corrected with (recorded as "
                "products.pbcor_status.ugmrt_pb_image), so that both "
                "published images share one off-axis flux scale. Run the "
                "imaging stage successfully before DD calibration, or use "
                "--skip-imaging and correct the DDFacet result separately "
                "with a validated beam model."
            )

        expected = (
            self.imagedir /
            (
                f"{clean_name(self.target)}_{self.mode}_"
                f"{self.solver_signature()[:10]}_DD.pbcor.fits"
            )
        )
        metrics_file = self.reportdir / "dd_image_metrics.json"
        post_inputs = [self.dd_image_products[0], pb_image]
        post_products = [expected, metrics_file]

        if self.checkpoint_exists("post_dd_imaging"):
            valid, reason = self.bound_checkpoint_status(
                "post_dd_imaging",
                post_products,
                inputs=post_inputs,
            )
            if valid:
                print("DD primary-beam correction checkpoint is valid.")
                self.pbcor_products = [str(expected)]
                self.dd_metrics = load_json(metrics_file, {})
                return

            self.invalidate_post_dd_state(reason, post_products)

        script, output_fits, metrics_path = self.write_dd_pbcor_script(
            self.dd_image_products[0],
            pb_image,
        )
        command = [
            self.casa,
            "--nologger",
            "--nogui",
            "-c",
            str(script),
        ]

        record = self.run_command(
            "dd_primary_beam_correction",
            command,
            required=True,
        )

        if (
            record.get("status") == "ok" and
            output_fits.is_file() and
            metrics_path.is_file()
        ):

            self.mark_bound_checkpoint(
                "post_dd_imaging",
                [output_fits, metrics_path],
                inputs=post_inputs,
            )
            self.pbcor_products = [str(output_fits)]
            self.dd_metrics = load_json(metrics_path, {})

        else:
            raise RuntimeError(
                "DD primary-beam correction did not create "
                f"{output_fits} and its metrics report."
            )

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    @staticmethod
    def materialized_product_files(value):

        if isinstance(value, dict):
            for item in value.values():
                yield from DDCalManager.materialized_product_files(item)
            return
        if isinstance(value, (list, tuple)):
            for item in value:
                yield from DDCalManager.materialized_product_files(item)
            return
        if value:
            path = Path(str(value)).expanduser()
            if path.is_file():
                yield path.resolve()

    # --------------------------------------------------------

    def write_report(self):

        try:
            solver_signature = self.solver_signature()
        except Exception:
            solver_signature = None

        try:
            discovered_dd_images = self.find_dd_image_products()
        except Exception:
            discovered_dd_images = []

        try:
            measurement_set_stat_sha256 = (
                self.measurement_set_stat_identity()
            )
        except (OSError, TypeError, ValueError):
            measurement_set_stat_sha256 = None

        try:
            solution_file = self.solver_paths()["solution_file"]
        except Exception:
            solution_file = None

        products = {
            "source_catalog_csv": str(self.catalog_csv),
            "source_catalog_fits": str(self.catalog_fits),
            "source_catalog_region": str(self.catalog_region),
            "directions_csv": str(self.directions_csv),
            "directions_region": str(self.directions_region),
            "solver_directions_region": str(
                self.solver_directions_region
            ),
            "facet_region": str(self.facet_region),
            "qa_plot_status": str(
                self.root / "qa_plot_status.json"
            ),
            "skymodel": str(self.skymodel),
            "solver_config": str(
                self.configdir / "ddcal_config.json"
            ),
            "command_script": str(self.command_script),
            "directions_input_provenance": str(
                self.checkdir / "directions_inputs.json"
            ),
            "direction_selection": str(
                self.directiondir / "direction_selection.json"
            ),
            "uncovered_bright_sources_csv": (
                getattr(self, "direction_selection", {}) or {}
            ).get("uncovered_bright_sources_csv"),
            "solver_checkpoint_provenance": str(
                self.checkpoint_provenance_path("solver")
            ),
            "post_dd_checkpoint_provenance": str(
                self.checkpoint_provenance_path("post_dd_imaging")
            ),
            "ddcal_checkpoint_provenance": str(
                self.checkpoint_provenance_path("ddcal")
            ),
            "solutions_dir": str(self.soldir),
            "solver_solution_file": (
                str(solution_file) if solution_file is not None else None
            ),
            "solution_qa": str(
                self.soldir / "killms_solution_qa.json"
            ),
            "solution_qa_plot": str(
                self.plotdir / "killms_gain_amplitudes.png"
            ),
            "dd_restored_fits": (
                self.dd_image_products or discovered_dd_images
            ),
            "dd_pbcor_fits": self.pbcor_products,
            "dd_image_metrics": str(
                self.reportdir / "dd_image_metrics.json"
            ),
        }
        published_files = sorted({
            str(path)
            for path in self.materialized_product_files(products)
        })
        product_content_identities = {
            path: IMAGING_SAFETY.file_content_identity(path)
            for path in published_files
        }

        report = {
            "schema_version": 2,
            "created_utc": now_utc(),
            "status": self.workflow_status,
            "reason": self.failure_reason,
            "provenance_fingerprint": self.provenance_fingerprint,
            "measurement_set": self.ms,
            "measurement_set_stat_sha256": measurement_set_stat_sha256,
            "source_catalog_content_identity": self.file_fingerprint(
                self.catalog_csv
            ),
            "target": self.target,
            "band": self.band,
            "reference_image": self.image,
            "engine": self.engine,
            "requested_engine": self.requested_engine,
            "mode": self.mode,
            "requested_mode": self.requested_mode,
            "solver_signature": solver_signature,
            "execute": self.execute,
            "strict": self.strict,
            "configuration": self.config,
            "direction_selection": getattr(
                self,
                "direction_selection",
                {},
            ),
            "solver_direction_count": getattr(
                self,
                "solver_direction_count",
                0,
            ),
            "products": products,
            "published_files": published_files,
            "product_content_identities": product_content_identities,
            "metrics": self.dd_metrics,
            "directions": [
                direction.to_dict()
                for direction in getattr(
                    self,
                    "directions",
                    [],
                )
            ],
            "commands": self.command_records,
            "checkpoints": {
                "source_catalog": self.checkpoint_exists(
                    "source_catalog"
                ),
                "directions": self.checkpoint_exists(
                    "directions"
                ),
                "solver": self.checkpoint_exists(
                    "solver"
                ),
                "post_dd_imaging": self.checkpoint_exists(
                    "post_dd_imaging"
                ),
                "ddcal": self.checkpoint_exists(
                    "ddcal"
                ),
            },
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
                "imaging": str(
                    self.outdir /
                    "imaging" /
                    "reports" /
                    "imaging_report.json"
                ),
                "source_detection": str(
                    self.outdir /
                    "source_detection" /
                    "reports" /
                    "source_detection_report.json"
                ),
            },
        }

        outfile = (
            self.reportdir /
            "ddcal_report.json"
        )

        save_json(report, outfile)

        print(
            f"\nDD calibration report written:\n{outfile}"
        )

    # --------------------------------------------------------
    # Driver
    # --------------------------------------------------------

    def run(self):

        print("\n")
        print("=" * 70)
        print("Starting uGMRT DD Calibration Workflow")
        print("=" * 70)

        try:
            self.print_summary()

            self.extract_sources()

            self.make_directions()

            self.generate_qa_plots()

            self.run_solver()

            self.run_post_dd_imaging()

            if self.execute and not self.skip_solver:
                if not self.checkpoint_exists("solver"):
                    raise RuntimeError(
                        "DD execution ended without a valid solver checkpoint."
                    )
                self.dd_image_products = self.find_dd_image_products()
                if (
                    self.engine == "killms-ddfacet" and
                    not self.dd_image_products
                ):
                    raise RuntimeError(
                        "DD execution ended without a DDFacet restored image."
                    )
                if (
                    self.engine == "killms-ddfacet" and
                    not self.skip_imaging and
                    not self.checkpoint_exists("post_dd_imaging")
                ):
                    raise RuntimeError(
                        "DD execution ended without the requested "
                        "primary-beam-corrected image."
                    )
                self.workflow_status = "success"
            elif self.execute and self.skip_solver:
                self.workflow_status = "skipped"
            else:
                self.workflow_status = "prepared"

            final_products = self.current_solver_products()
            final_products.extend(self.pbcor_products)
            metrics_path = self.reportdir / "dd_image_metrics.json"
            if metrics_path.is_file():
                final_products.append(str(metrics_path.resolve()))
            self.mark_bound_checkpoint("ddcal", final_products)
            self.write_report()

        except Exception as exc:
            self.workflow_status = "failed"
            self.failure_reason = str(exc)
            self.archive_checkpoint("ddcal")
            self.write_report()
            raise

        print("\n")
        print("=" * 70)
        print(
            "uGMRT DD calibration workflow completed "
            f"with status: {self.workflow_status}."
        )
        print("=" * 70)


# ============================================================
# BACKWARD-COMPATIBLE SMALL HELPERS
# ============================================================

def run(cmd):

    command = split_command(cmd)
    subprocess.run(command, check=True)


def run_pybdsf(image, outfile):

    outfile = Path(outfile)
    outfile.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    script = outfile.parent / "run_pybdsf.py"
    fits_out = outfile.with_suffix(".fits")
    reg_out = outfile.with_suffix(".reg")

    lines = [
        "import bdsf",
        f"image = {str(image)!r}",
        f"csv_path = {str(outfile)!r}",
        f"fits_path = {str(fits_out)!r}",
        f"reg_path = {str(reg_out)!r}",
        "img = bdsf.process_image(",
        "    image,",
        "    adaptive_rms_box=True,",
        "    atrous_do=True,",
        "    thresh_isl=4.0,",
        "    thresh_pix=5.0,",
        ")",
        "img.write_catalog(",
        "    outfile=csv_path,",
        "    catalog_type='srl',",
        "    format='csv',",
        "    clobber=True,",
        ")",
        "img.write_catalog(",
        "    outfile=fits_path,",
        "    catalog_type='srl',",
        "    format='fits',",
        "    clobber=True,",
        ")",
        "img.write_catalog(",
        "    outfile=reg_path,",
        "    catalog_type='srl',",
        "    format='ds9',",
        "    clobber=True,",
        ")",
    ]

    script.write_text("\n".join(lines) + "\n")

    run([
        sys.executable,
        str(script),
    ])


def run_cubical(ms, ncpu=64):

    exe = which_any(["gocubical", "cubical"])

    if exe is None:

        print("CubiCal not found.")
        return

    command = [
        exe,
        str(ms),
        "--sol-jones",
        "G",
        "--g-type",
        "complex-2x2",
        "--dist-ncpu",
        str(ncpu),
        "--out-casa-gaintables",
    ]

    run(command)


# ============================================================
# CLI
# ============================================================

def build_parser():

    parser = argparse.ArgumentParser(
        description=(
            "Direction-dependent calibration workflow for the uGMRT continuum pipeline."
        )
    )

    parser.add_argument(
        "ms",
        nargs="?",
        default=None,
        help=(
            "Input Measurement Set. If omitted, ddcal.py uses the imaging/selfcal reports."
        ),
    )

    parser.add_argument(
        "--image",
        default=None,
        help=(
            "Final continuum image/FITS used for PyBDSF source extraction."
        ),
    )

    parser.add_argument(
        "--catalog",
        default=None,
        help=(
            "Existing source catalogue CSV/TXT. Must contain RA, DEC, and flux columns."
        ),
    )

    parser.add_argument(
        "--directions",
        default=None,
        help=(
            "Existing DD direction CSV. Bypasses source clustering."
        ),
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
        "--target",
        default=RUN_CONFIG.get("target"),
        help="Target/product name.",
    )

    parser.add_argument(
        "--mode",
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
            "DD strategy: a compact direction set for one object, a "
            "multi-facet survey strategy, or automatic selection."
        ),
    )

    parser.add_argument(
        "--engine",
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
        help="DD solver backend.",
    )

    parser.add_argument(
        "--prepare-only",
        action="store_true",
        default=bool(
            RUN_CONFIG.get(
                "dd_prepare_only",
                not DD_TASK_CONFIG.get(
                    "execute",
                    True,
                ),
            )
        ),
        help=(
            "Prepare catalogues/directions/configs only; do not execute external commands."
        ),
    )

    parser.add_argument(
        "--execute",
        action="store_false",
        dest="prepare_only",
        help="Execute external DD commands after preparing products.",
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        default=bool(
            DD_TASK_CONFIG.get(
                "strict",
                False,
            )
        ),
        help="Raise errors if optional external tools fail.",
    )

    parser.add_argument(
        "--non-strict",
        action="store_false",
        dest="strict",
        help="Continue when optional external DD tools fail.",
    )

    parser.add_argument(
        "--ncpu",
        type=int,
        default=RUN_CONFIG.get(
            "ncpu",
            32,
        ),
        help="Number of CPUs for external DD tools.",
    )

    parser.add_argument(
        "--data-column",
        default=DD_TASK_CONFIG.get(
            "data_column",
            "auto",
        ),
        help=(
            "Input data column for DD solving. 'auto' reuses the column "
            "recorded by the final imaging stage."
        ),
    )

    parser.add_argument(
        "--model-column",
        default=DD_TASK_CONFIG.get(
            "model_column",
            "MODEL_DATA",
        ),
        help="Model column/list used by DD solver.",
    )

    parser.add_argument(
        "--min-flux-jy",
        type=float,
        default=None,
        help="Minimum source flux density for DD directions.",
    )

    parser.add_argument(
        "--max-directions",
        type=int,
        default=None,
        help="Maximum number of DD directions.",
    )

    parser.add_argument(
        "--merge-radius-arcmin",
        type=float,
        default=None,
        help="Merge bright sources within this radius.",
    )

    parser.add_argument(
        "--facet-radius-arcmin",
        type=float,
        default=None,
        help="DS9 facet radius for QA/inspection.",
    )

    parser.add_argument(
        "--solint-time",
        default=None,
        help="Time solution interval for DD gains.",
    )

    parser.add_argument(
        "--solint-freq",
        default=None,
        help="Frequency solution interval for DD gains.",
    )

    parser.add_argument(
        "--solver-command-template",
        default=DD_TASK_CONFIG.get(
            "solver_command_template"
        ),
        help=(
            "Custom CubiCal solve command template. Available "
            "fields include {ms}, {directions}, {skymodel}, "
            "{solutions_dir}, {ncpu}, {data_column}, "
            "{model_column}, {solint_time}, and {solint_freq}."
        ),
    )

    parser.add_argument(
        "--apply-command-template",
        default=DD_TASK_CONFIG.get(
            "apply_command_template"
        ),
        help="Custom CubiCal apply command template.",
    )

    parser.add_argument(
        "--killms-command-template",
        default=DD_TASK_CONFIG.get(
            "killms_command_template"
        ),
        help="Custom KillMS command template.",
    )

    parser.add_argument(
        "--ddfacet-command-template",
        default=DD_TASK_CONFIG.get(
            "ddfacet_command_template"
        ),
        help="Custom DDFacet command template.",
    )

    parser.add_argument(
        "--skip-source-extraction",
        action="store_true",
        help="Do not run PyBDSF.",
    )

    parser.add_argument(
        "--skip-solver",
        action="store_true",
        help="Do not run DD solver commands.",
    )

    parser.add_argument(
        "--skip-imaging",
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
        help="Do not run post-DD imaging.",
    )

    parser.add_argument(
        "--post-imaging",
        action="store_false",
        dest="skip_imaging",
        help=(
            "For killMS/DDFacet, regrid the CASA uGMRT primary beam and "
            "write a PB-corrected DD FITS image."
        ),
    )

    parser.add_argument(
        "--imaging-script",
        default=None,
        help="Path to imaging.py for post-DD imaging.",
    )

    parser.add_argument(
        "--casa",
        default=RUN_CONFIG.get(
            "casa",
            "casa",
        ),
        help="CASA executable for post-DD imaging.",
    )

    parser.add_argument(
        "--provenance-fingerprint",
        default=None,
        help=(
            "Orchestrator-provided DD input/recipe fingerprint recorded in "
            "the final report."
        ),
    )

    return parser


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":

    args = build_parser().parse_args()

    DDCalManager(
        ms=args.ms,
        image=args.image,
        catalog=args.catalog,
        directions=args.directions,
        outdir=args.outdir,
        target=args.target,
        mode=args.mode,
        engine=args.engine,
        execute=not args.prepare_only,
        strict=args.strict,
        ncpu=args.ncpu,
        data_column=args.data_column,
        model_column=args.model_column,
        min_flux_jy=args.min_flux_jy,
        max_directions=args.max_directions,
        merge_radius_arcmin=args.merge_radius_arcmin,
        facet_radius_arcmin=args.facet_radius_arcmin,
        solint_time=args.solint_time,
        solint_freq=args.solint_freq,
        solver_command_template=args.solver_command_template,
        apply_command_template=args.apply_command_template,
        killms_command_template=args.killms_command_template,
        ddfacet_command_template=args.ddfacet_command_template,
        skip_source_extraction=args.skip_source_extraction,
        skip_solver=args.skip_solver,
        skip_imaging=args.skip_imaging,
        imaging_script=args.imaging_script,
        casa=args.casa,
        provenance_fingerprint=args.provenance_fingerprint,
    ).run()
