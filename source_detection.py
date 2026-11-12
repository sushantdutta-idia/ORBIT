#!/usr/bin/env python3
# ============================================================
# source_detection.py
#
# Optional PyBDSF source-finding stage for the radio continuum
# pipeline. It consumes FITS products from either imaging backend,
# creates source/Gaussian catalogues and exports diagnostic images.
#
# The PB-corrected Stokes-I image is used for science catalogues so
# off-axis flux densities are corrected. Adaptive local-RMS estimation
# handles its position-dependent noise inside the configured PB cutoff.
#
# Author: Sushant Dutta
# ============================================================

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def load_config_module():
    path = SCRIPT_DIR / "config.py"
    if not path.exists():
        return None

    spec = importlib.util.spec_from_file_location(
        "ugmrt_pipeline_source_detection_config",
        str(path),
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
TASK_CONFIG = (
    getattr(CONFIG, "SOURCE_DETECTION_CONFIG", {})
    if CONFIG is not None
    else {}
)
DEFAULT_OUTDIR = (
    getattr(CONFIG, "DEFAULT_OUTDIR", "Products")
    if CONFIG is not None
    else "Products"
)


def now_utc():
    return datetime.now(timezone.utc).isoformat()


def atomic_write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")

    with open(temporary, "w") as handle:
        json.dump(data, handle, indent=2, default=str)
        handle.flush()
        os.fsync(handle.fileno())

    os.replace(temporary, path)


def file_content_identity(path, chunk_size=8 * 1024 * 1024):
    product = Path(path).expanduser().resolve()
    if not product.is_file():
        raise ValueError(f"Published source-detection file is missing: {product}")
    digest = hashlib.sha256()
    with product.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    stat = product.stat()
    return {
        "path": str(product),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
        "sha256": digest.hexdigest(),
    }


def validate_product_identities(products, identities):
    materialized = {
        str(Path(value).expanduser().resolve())
        for value in (products or {}).values()
        if value and Path(value).is_file()
    }
    if not isinstance(identities, dict) or set(identities) != materialized:
        return False
    try:
        return all(
            identities[path] == file_content_identity(path)
            for path in materialized
        )
    except (OSError, ValueError):
        return False


def read_json(path, default=None):
    path = Path(path)
    if not path.is_file():
        return default

    with open(path) as handle:
        return json.load(handle)


def resolve_path(value, base):
    if not value:
        return None

    path = Path(value).expanduser()
    if path.is_absolute():
        return path

    candidates = [
        Path.cwd() / path,
        SCRIPT_DIR / path,
        Path(base) / path,
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    return (Path.cwd() / path).resolve()


def clean_name(value):
    value = str(value or "source").strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_.-") or "source"


def normalized_fits_products(report, outdir):
    products = report.get("products", {}) if isinstance(report, dict) else {}
    fits_products = products.get("fits", []) or []
    resolved = []

    for item in fits_products:
        path = resolve_path(item, outdir)
        if path is not None and path.is_file():
            resolved.append(path)

    return resolved


def expected_fits_for_casa_image(casa_image, outdir):
    if not casa_image:
        return None

    return Path(outdir) / "imaging" / "fits" / (
        f"{Path(str(casa_image)).name}.fits"
    )


def score_fits_candidate(path, preference):
    name = path.name.lower()
    score = 0

    if preference == "pbcor":
        score += 100 if "pbcor" in name else -20
    else:
        score += -60 if "pbcor" in name else 20

    if ".image.tt0.fits" in name:
        score += 60
    elif ".image.fits" in name:
        score += 45
    elif "image" in name:
        score += 20

    for unwanted in [
        "residual",
        "model",
        "mask",
        "sumwt",
        "alpha",
        "beta",
        "psf",
        "pb.tt0",
    ]:
        if unwanted in name:
            score -= 100

    if "subband" in name or re.search(r"spw\d+", name):
        score -= 25

    return score


def select_input_image(explicit_image, report, outdir, preference):
    if explicit_image:
        image = resolve_path(explicit_image, outdir)
        if image is None or not image.is_file():
            raise FileNotFoundError(
                f"Source-detection image not found: {explicit_image}"
            )
        return image

    products = report.get("products", {}) if isinstance(report, dict) else {}
    key = "pbcor_image" if preference == "pbcor" else "restored_image"
    role_product = products.get(key)

    if preference == "pbcor" and not role_product:
        raise RuntimeError(
            "PB-corrected source detection was requested, but the imaging "
            "report records no PB-corrected product. Select --image-"
            "preference restored only if an attenuated, flat-noise image is "
            "intentional."
        )

    # Backend-neutral reports may point directly at a FITS image (WSClean),
    # while CASA reports point at an image table that is exported below the
    # imaging/fits directory.  Prefer the explicit product role before using
    # filename heuristics so a LOFAR channel or uncorrected MFS image cannot
    # accidentally outrank the recorded PB-corrected continuum product.
    selected_product = resolve_path(role_product, outdir)

    if selected_product is not None and selected_product.is_file():
        return selected_product.resolve()

    expected = expected_fits_for_casa_image(role_product, outdir)

    if expected is not None and expected.is_file():
        return expected.resolve()

    if role_product:
        # The imaging report names a product for this role, so the role is
        # decided; only its materialized FITS is missing.  Falling through to
        # filename scoring here would silently return a *different* image --
        # for a pbcor request, the uncorrected restored image scores well
        # enough to win, and the catalogue would carry primary-beam
        # attenuated off-axis flux densities while the report still claimed
        # image_preference=pbcor.  Fail closed instead.
        expected_name = Path(str(role_product)).name + ".fits"
        for candidate in normalized_fits_products(report, outdir):
            if candidate.name == expected_name:
                return candidate.resolve()

        raise FileNotFoundError(
            f"The imaging report records {key}={role_product!r}, but no "
            "materialized FITS export of it was found (looked for "
            f"{expected_name} under {Path(outdir) / 'imaging' / 'fits'}). "
            "Re-run the imaging stage's FITS export, or pass --image "
            "explicitly. Refusing to substitute a different image."
        )

    fits_products = normalized_fits_products(report, outdir)
    if not fits_products:
        fits_products = sorted(
            (Path(outdir) / "imaging" / "fits").glob("*.fits")
        )

    if not fits_products:
        raise FileNotFoundError(
            "No imaging FITS products were found. Run imaging.py before "
            "source_detection.py or pass --image explicitly."
        )

    # The heuristic path is now reachable only when the report names no
    # product for this role at all.  Keep it from ever contradicting the
    # requested preference.
    wants_pbcor = preference == "pbcor"
    fits_products = [
        path
        for path in fits_products
        if ("pbcor" in path.name.lower()) == wants_pbcor
    ]

    if not fits_products:
        raise RuntimeError(
            "No imaging FITS product matches the requested "
            f"--image-preference {preference}. Pass --image explicitly."
        )

    ranked = sorted(
        fits_products,
        key=lambda path: score_fits_candidate(path, preference),
        reverse=True,
    )
    selected = ranked[0]

    if score_fits_candidate(selected, preference) < 0:
        raise RuntimeError(
            "Could not identify a suitable science image among "
            "the imaging FITS products. Pass --image explicitly."
        )

    return selected.resolve()


def image_selection_provenance(image_path, report, outdir, preference):
    """State which imaging role the selected image actually fulfils.

    The requested preference is an intent; this records the outcome, so a
    catalogue can never claim a PB-corrected flux scale without the image to
    support it.
    """

    products = report.get("products", {}) if isinstance(report, dict) else {}
    resolved = Path(image_path).expanduser().resolve()

    roles = {}
    for key in ("pbcor_image", "restored_image"):
        product = products.get(key)
        if not product:
            continue
        candidates = set()
        direct = resolve_path(product, outdir)
        if direct is not None:
            candidates.add(direct.resolve())
        expected = expected_fits_for_casa_image(product, outdir)
        if expected is not None:
            candidates.add(expected.resolve())
        roles[key] = candidates

    matched = [key for key, paths in roles.items() if resolved in paths]
    role = matched[0] if matched else None

    return {
        "requested_preference": preference,
        "resolved_role": role,
        "is_primary_beam_corrected": (
            True if role == "pbcor_image"
            else False if role == "restored_image"
            else None
        ),
        "matches_requested_preference": (
            None
            if role is None
            else role == ("pbcor_image" if preference == "pbcor"
                          else "restored_image")
        ),
    }


def count_catalog_rows(path, id_column):
    path = Path(path)
    if not path.is_file():
        return None

    with open(path, newline="") as handle:
        lines = [
            line for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]

    if not lines:
        return 0

    reader = csv.DictReader(lines)
    if reader.fieldnames and id_column in reader.fieldnames:
        return sum(1 for row in reader if row.get(id_column, "").strip())

    # PyBDSF variants may emit a non-comment metadata line before the
    # actual CSV header. Locate the header explicitly before giving up.
    for index, line in enumerate(lines):
        if id_column in [item.strip() for item in line.split(",")]:
            reader = csv.DictReader(lines[index:])
            return sum(
                1 for row in reader
                if row.get(id_column, "").strip()
            )

    return None


class SourceDetectionManager:

    def __init__(self, args):
        self.args = args
        self.outdir = Path(args.outdir).expanduser().resolve()
        self.root = self.outdir / "source_detection"
        self.catalogdir = self.root / "catalogs"
        self.imagedir = self.root / "images"
        self.logdir = self.root / "logs"
        self.reportdir = self.root / "reports"
        self.checkdir = self.root / "checkpoints"
        self.report_path = self.reportdir / "source_detection_report.json"
        self.checkpoint = self.checkdir / "source_detection.done"
        self.imaging_report_path = (
            self.outdir / "imaging" / "reports" / "imaging_report.json"
        )
        self.imaging_report = read_json(self.imaging_report_path, {})

        for directory in [
            self.catalogdir,
            self.imagedir,
            self.logdir,
            self.reportdir,
            self.checkdir,
        ]:
            directory.mkdir(parents=True, exist_ok=True)

        self.config = dict(TASK_CONFIG)
        self.config.update({
            "enabled": args.enabled,
            "strict": args.strict,
            "image_preference": args.image_preference,
            "thresh_isl": args.thresh_isl,
            "thresh_pix": args.thresh_pix,
            "adaptive_rms_box": args.adaptive_rms_box,
            "atrous_do": args.atrous_do,
            "quiet": args.quiet,
        })

        if args.rms_box is not None:
            self.config["rms_box"] = tuple(args.rms_box)

        self.warnings = []
        self.products = {}

    def previous_success_is_valid(self):
        if self.args.overwrite or not self.checkpoint.is_file():
            return False

        report = read_json(self.report_path, {})
        if report.get("status") != "success":
            return False

        products = report.get("products", {})
        required = [
            products.get("source_catalog_csv"),
            products.get("source_catalog_fits"),
            products.get("source_regions"),
        ]

        return all(
            value and Path(value).is_file()
            for value in required
        ) and validate_product_identities(
            products,
            report.get("product_content_identities"),
        )

    def write_report(self, status, reason=None, image=None, engine_version=None):
        source_count = count_catalog_rows(
            self.products.get("source_catalog_csv", ""),
            "Source_id",
        )
        gaussian_count = count_catalog_rows(
            self.products.get("gaussian_catalog_csv", ""),
            "Gaus_id",
        )

        product_content_identities = {
            str(Path(value).expanduser().resolve()): file_content_identity(value)
            for value in self.products.values()
            if value and Path(value).is_file()
        }
        report = {
            "created_utc": now_utc(),
            "status": status,
            "reason": reason,
            "engine": "pybdsf",
            "engine_version": engine_version,
            "input_image": str(image) if image is not None else None,
            "image_selection": (
                image_selection_provenance(
                    image,
                    self.imaging_report,
                    self.outdir,
                    self.config.get("image_preference", "pbcor"),
                )
                if image is not None
                else None
            ),
            "target": (
                self.args.target or self.imaging_report.get("target")
            ),
            "telescope": self.imaging_report.get("telescope"),
            "band": self.imaging_report.get("band"),
            "imaging_backend": self.imaging_report.get(
                "imaging_backend",
                "casa",
            ),
            "primary_beam": {
                "method": self.imaging_report.get("products", {}).get(
                    "pbcor_method"
                ),
                "status": self.imaging_report.get("products", {}).get(
                    "pbcor_status"
                ),
            },
            "parameters": {
                key: self.config.get(key)
                for key in [
                    "image_preference",
                    "thresh_isl",
                    "thresh_pix",
                    "adaptive_rms_box",
                    "atrous_do",
                    "rms_box",
                    "rms_map",
                    "mean_map",
                    "flag_maxsize_bm",
                    "quiet",
                ]
            },
            "counts": {
                "sources": source_count,
                "gaussians": gaussian_count,
            },
            "products": self.products,
            "product_content_identities": product_content_identities,
            "warnings": self.warnings,
            "input_reports": {
                "imaging": str(self.imaging_report_path),
            },
        }

        atomic_write_json(self.report_path, report)
        return report

    def mark_done(self):
        self.checkpoint.write_text(now_utc() + "\n")

    def optional_catalog(self, image, outfile, catalog_type, fmt, **kwargs):
        try:
            image.write_catalog(
                outfile=str(outfile),
                catalog_type=catalog_type,
                format=fmt,
                clobber=True,
                **kwargs,
            )
        except Exception as exc:
            self.warnings.append(
                f"Could not write {catalog_type}/{fmt} catalogue "
                f"{outfile.name}: {exc}"
            )
            return None

        return str(outfile.resolve()) if outfile.is_file() else None

    def optional_image(self, image, outfile, img_type):
        try:
            image.export_image(
                outfile=str(outfile),
                img_type=img_type,
                img_format="fits",
                clobber=True,
            )
        except Exception as exc:
            self.warnings.append(
                f"Could not export PyBDSF {img_type} image "
                f"{outfile.name}: {exc}"
            )
            return None

        return str(outfile.resolve()) if outfile.is_file() else None

    def run(self):
        if not self.config.get("enabled", True):
            reason = "Source detection is disabled in configuration."
            print(reason)
            self.write_report("disabled", reason=reason)
            self.mark_done()
            return

        if self.previous_success_is_valid():
            print("PyBDSF source-detection products already exist.")
            return

        imaging_report = self.imaging_report
        if (
            self.args.target and
            imaging_report and
            imaging_report.get("target") != self.args.target
        ):
            raise RuntimeError(
                f"The imaging report belongs to target "
                f"{imaging_report.get('target')!r}, not requested source "
                f"detection target {self.args.target!r}."
            )
        image_path = select_input_image(
            self.args.image,
            imaging_report,
            self.outdir,
            self.config.get("image_preference", "pbcor"),
        )

        try:
            import bdsf
        except Exception as exc:
            reason = (
                "PyBDSF is unavailable in the selected pipeline Python "
                f"environment: {exc}"
            )
            self.write_report("skipped", reason=reason, image=image_path)
            self.mark_done()

            if self.config.get("strict", False):
                raise RuntimeError(reason) from exc

            print(f"WARNING: {reason}")
            print("Source detection was skipped; imaging products are unchanged.")
            return

        version = getattr(bdsf, "__version__", None)
        stem = clean_name(
            self.args.target
            or imaging_report.get("target")
            or image_path.stem
        )

        process_kwargs = {
            "interactive": False,
            "quiet": bool(self.config.get("quiet", True)),
            "adaptive_rms_box": bool(
                self.config.get("adaptive_rms_box", True)
            ),
            "atrous_do": bool(self.config.get("atrous_do", True)),
            "thresh_isl": float(self.config.get("thresh_isl", 4.0)),
            "thresh_pix": float(self.config.get("thresh_pix", 5.0)),
            "outdir": str(self.logdir),
        }

        for key in [
            "rms_box",
            "rms_map",
            "mean_map",
            "flag_maxsize_bm",
            "frequency",
            "beam",
        ]:
            value = self.config.get(key)
            if value is not None:
                process_kwargs[key] = value

        print("=" * 70)
        print("PyBDSF source detection")
        print(f"Input image     : {image_path}")
        print(f"Image selection : {self.config.get('image_preference')}")
        print(f"Thresholds      : {process_kwargs['thresh_isl']} / "
              f"{process_kwargs['thresh_pix']} sigma")
        print(f"Adaptive RMS    : {process_kwargs['adaptive_rms_box']}")
        print(f"Wavelets        : {process_kwargs['atrous_do']}")
        print(f"Output          : {self.root}")
        print("=" * 70)

        selection = image_selection_provenance(
            image_path,
            imaging_report,
            self.outdir,
            self.config.get("image_preference", "pbcor"),
        )
        if selection["matches_requested_preference"] is False:
            self.warnings.append(
                f"The selected image fulfils the {selection['resolved_role']} "
                f"role, not the requested {self.config.get('image_preference')} "
                "preference; catalogue flux densities are on that image's "
                "scale."
            )
            print(f"WARNING: {self.warnings[-1]}")
        elif selection["resolved_role"] is None:
            self.warnings.append(
                "The selected image does not correspond to any product role "
                "recorded by the imaging report; its flux scale is not "
                "established by pipeline provenance."
            )
            print(f"WARNING: {self.warnings[-1]}")

        try:
            result = bdsf.process_image(str(image_path), **process_kwargs)

            source_csv = self.catalogdir / f"{stem}_sources.csv"
            source_fits = self.catalogdir / f"{stem}_sources.fits"
            source_reg = self.catalogdir / f"{stem}_sources.reg"
            gaul_csv = self.catalogdir / f"{stem}_gaussians.csv"
            gaul_fits = self.catalogdir / f"{stem}_gaussians.fits"
            gaul_bbs = self.catalogdir / f"{stem}_gaussians.bbs"

            self.products["source_catalog_csv"] = self.optional_catalog(
                result, source_csv, "srl", "csv"
            )
            self.products["source_catalog_fits"] = self.optional_catalog(
                result, source_fits, "srl", "fits"
            )
            self.products["source_regions"] = self.optional_catalog(
                result, source_reg, "srl", "ds9"
            )
            self.products["gaussian_catalog_csv"] = self.optional_catalog(
                result, gaul_csv, "gaul", "csv"
            )
            self.products["gaussian_catalog_fits"] = self.optional_catalog(
                result, gaul_fits, "gaul", "fits"
            )
            self.products["bbs_sky_model"] = self.optional_catalog(
                result,
                gaul_bbs,
                "gaul",
                "bbs",
                bbs_patches="source",
            )

            image_exports = {
                "model_image": ("gaus_model", f"{stem}_pybdsf_model.fits"),
                "residual_image": (
                    "gaus_resid",
                    f"{stem}_pybdsf_residual.fits",
                ),
                "rms_image": ("rms", f"{stem}_pybdsf_rms.fits"),
                "island_mask": (
                    "island_mask",
                    f"{stem}_pybdsf_islands.fits",
                ),
            }

            if self.config.get("export_mean", False):
                image_exports["mean_image"] = (
                    "mean",
                    f"{stem}_pybdsf_mean.fits",
                )

            for key, (img_type, filename) in image_exports.items():
                self.products[key] = self.optional_image(
                    result,
                    self.imagedir / filename,
                    img_type,
                )

            required = [
                self.products.get("source_catalog_csv"),
                self.products.get("source_catalog_fits"),
                self.products.get("source_regions"),
            ]

            if not all(required):
                raise RuntimeError(
                    "PyBDSF completed, but one or more required source "
                    "catalogue products were not written."
                )

            report = self.write_report(
                "success",
                image=image_path,
                engine_version=version,
            )
            self.mark_done()

        except Exception as exc:
            self.write_report(
                "failed",
                reason=str(exc),
                image=image_path,
                engine_version=version,
            )
            if self.checkpoint.exists():
                self.checkpoint.unlink()
            raise

        counts = report.get("counts", {})
        print("\nSource detection complete.")
        print(f"Sources          : {counts.get('sources')}")
        print(f"Gaussian comps.  : {counts.get('gaussians')}")
        print(f"Catalogue        : {self.products['source_catalog_csv']}")
        print(f"DS9 regions      : {self.products['source_regions']}")
        print(f"Report           : {self.report_path}")


def bool_default(key, fallback):
    return bool(TASK_CONFIG.get(key, fallback))


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Run PyBDSF source detection on the final radio-continuum "
            "FITS image."
        )
    )

    parser.add_argument(
        "image",
        nargs="?",
        default=None,
        help=(
            "Optional explicit FITS image. By default the image is resolved "
            "from imaging/reports/imaging_report.json."
        ),
    )
    parser.add_argument(
        "--outdir",
        default=RUN_CONFIG.get("outdir", DEFAULT_OUTDIR),
        help="Pipeline output directory.",
    )
    parser.add_argument(
        "--target",
        default=RUN_CONFIG.get("target"),
        help="Target label used in output filenames.",
    )
    parser.add_argument(
        "--image-preference",
        choices=["restored", "pbcor"],
        default=TASK_CONFIG.get("image_preference", "pbcor"),
        help=(
            "Prefer the restored or PB-corrected FITS image when --image "
            "is omitted. PB-corrected is recommended for science fluxes; "
            "adaptive local-RMS estimation handles its varying noise."
        ),
    )
    parser.add_argument(
        "--thresh-isl",
        type=float,
        default=float(TASK_CONFIG.get("thresh_isl", 4.0)),
        help="PyBDSF island threshold in local-RMS units.",
    )
    parser.add_argument(
        "--thresh-pix",
        type=float,
        default=float(TASK_CONFIG.get("thresh_pix", 5.0)),
        help="PyBDSF pixel threshold in local-RMS units.",
    )
    parser.add_argument(
        "--rms-box",
        type=int,
        nargs=2,
        metavar=("BOX", "STEP"),
        default=TASK_CONFIG.get("rms_box"),
        help="Optional PyBDSF RMS box and step sizes in pixels.",
    )
    parser.add_argument(
        "--adaptive-rms-box",
        action="store_true",
        dest="adaptive_rms_box",
        help="Use PyBDSF adaptive RMS boxes.",
    )
    parser.add_argument(
        "--no-adaptive-rms-box",
        action="store_false",
        dest="adaptive_rms_box",
        help="Disable PyBDSF adaptive RMS boxes.",
    )
    parser.set_defaults(
        adaptive_rms_box=bool_default("adaptive_rms_box", True)
    )
    parser.add_argument(
        "--atrous",
        action="store_true",
        dest="atrous_do",
        help="Enable wavelet fitting for extended emission.",
    )
    parser.add_argument(
        "--no-atrous",
        action="store_false",
        dest="atrous_do",
        help="Disable wavelet fitting for a faster compact-source catalogue.",
    )
    parser.set_defaults(atrous_do=bool_default("atrous_do", True))
    parser.add_argument(
        "--strict",
        action="store_true",
        dest="strict",
        help="Fail the pipeline when PyBDSF is unavailable.",
    )
    parser.add_argument(
        "--no-strict",
        action="store_false",
        dest="strict",
        help="Skip safely when PyBDSF is unavailable.",
    )
    parser.set_defaults(strict=bool_default("strict", False))
    parser.add_argument(
        "--enable",
        action="store_true",
        dest="enabled",
        help="Enable source detection.",
    )
    parser.add_argument(
        "--disable",
        action="store_false",
        dest="enabled",
        help="Write a disabled report without running PyBDSF.",
    )
    parser.set_defaults(enabled=bool_default("enabled", True))
    parser.add_argument(
        "--quiet",
        action="store_true",
        dest="quiet",
        help="Reduce PyBDSF terminal output.",
    )
    parser.add_argument(
        "--verbose",
        action="store_false",
        dest="quiet",
        help="Show PyBDSF terminal output.",
    )
    parser.set_defaults(quiet=bool_default("quiet", True))
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-run and overwrite existing source-detection products.",
    )

    return parser


if __name__ == "__main__":
    arguments = build_parser().parse_args()
    try:
        SourceDetectionManager(arguments).run()
    except Exception as exc:
        print(f"ERROR: source detection failed: {exc}", file=sys.stderr)
        raise
