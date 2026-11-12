#!/usr/bin/env python3
"""Correct a WSClean image of a uGMRT or MeerKAT field for the primary beam.

For MeerKAT (``--telescope meerkat``) the measured katbeam model is applied
by ``meerkat_beam.py``, as the pipeline's own imaging stage applies it: each
sub-band at the centre of its unflagged frequencies, and the MFS image by
the mean of the sub-bands' beams weighted as WSClean weighted them (finding
165). The frequency ranges the flagging stage masked are read from the
products (``--products``). It writes ``<name>-MFS-image-pb.fits``, the
beam as ``<name>-MFS-katbeam-pb.fits``, the same for each sub-band, and a
provenance record, and needs numpy, astropy and katbeam. What follows is
the uGMRT correction.

WSClean's images are as the telescope saw the sky: a source away from the
pointing centre is fainter than it is. The correction divides by the
primary beam. For a band as wide as Band 4 there is not one beam: the
image is the weighted mean of its sub-band images, each seen through the
beam of its own frequency, so the beam of the image is the same weighted
mean of those beams. That is what is divided by here, with the weights
WSClean itself recorded in the sub-band images, and the NCRA polynomial
for the band that the pipeline already holds
(``telescope_profiles.py``, ``imaging_safety.ugmrt_primary_beam_response``).

The polynomial describes the main lobe only. Outside the radius at which
the narrowest beam, that of the highest sub-band, has fallen to the
cut-off, and wherever the mean beam is below the cut-off, the corrected
image is blank.

The pointing centre is taken to be the image centre, which is where
WSClean puts the phase centre and where these targets were observed.

It needs numpy and nothing else. It writes ``<name>-MFS-image.pbcor.fits``
and ``<name>-MFS-pb.fits`` beside the image and changes nothing else.

Example::

    python3 correct_primary_beam.py \\
        --name run/target/G0212-0335/images/final/G0212-0335/wsclean_robust_0_I
"""

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# This folder holds the pipeline's own inspect.py; see diagnose_noise_floor.
if sys.path and os.path.abspath(sys.path[0]) == SCRIPT_DIR:
    sys.path.pop(0)

import argparse
import importlib.util
import json
from pathlib import Path


BLOCK = 2880


def load_sibling(filename, name):
    spec = importlib.util.spec_from_file_location(
        name, Path(SCRIPT_DIR) / filename
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SAFETY = load_sibling("imaging_safety.py", "correct_primary_beam_safety")
PROFILES = load_sibling("telescope_profiles.py", "correct_primary_beam_profiles")


def band_coefficients(band):
    """The four coefficients of a uGMRT band's beam, and where they hold."""

    try:
        profile = PROFILES.TELESCOPE_PROFILES["ugmrt"]["bands"][band]
    except KeyError:
        raise SystemExit(f"No uGMRT band '{band}' in telescope_profiles.py.")
    return (
        [float(value) for value in profile["primary_beam_coefficients"]],
        [float(value) for value in profile["primary_beam_valid_frequency_mhz"]],
    )


def read_fits(path):
    """A WSClean image as (header bytes, header cards, plane [y, x], tail)."""

    import numpy as np

    raw = Path(path).read_bytes()
    cards, offset, done = {}, 0, False
    while not done:
        block = raw[offset:offset + BLOCK]
        if len(block) < BLOCK:
            raise ValueError(f"{path} is not a complete FITS file")
        offset += BLOCK
        for index in range(0, BLOCK, 80):
            card = block[index:index + 80].decode("ascii", "replace")
            key = card[:8].strip()
            if key == "END":
                done = True
                break
            if card[8:10] == "= ":
                cards[key] = card[10:].split("/")[0].strip().strip("'").strip()
    if int(cards.get("BITPIX", 0)) != -32:
        raise ValueError(f"{path} is not a 32-bit floating-point image")
    nx, ny = int(cards["NAXIS1"]), int(cards["NAXIS2"])
    plane = np.frombuffer(
        raw, dtype=">f4", count=nx * ny, offset=offset
    ).reshape(ny, nx).astype(float)
    return raw[:offset], cards, plane, raw[offset + 4 * nx * ny:]


def write_fits(path, header, plane, tail):
    with open(path, "wb") as handle:
        handle.write(header)
        handle.write(plane.astype(">f4").tobytes())
        handle.write(tail)


def sub_bands(name):
    """Frequency (Hz) and weight of each sub-band image of a WSClean run.

    The weight is WSClean's ``WSCIMGWG``; an image without it counts as
    one, and a sub-band that holds no data has weight zero and is left out.
    """

    name = Path(name)
    found = []
    for path in sorted(
        name.parent.glob(name.name + "-[0-9][0-9][0-9][0-9]-image.fits")
    ):
        raw = path.read_bytes()[:10 * BLOCK]
        cards = {}
        for index in range(0, len(raw), 80):
            card = raw[index:index + 80].decode("ascii", "replace")
            key = card[:8].strip()
            if key == "END":
                break
            if card[8:10] == "= ":
                cards[key] = card[10:].split("/")[0].strip().strip("'").strip()
        try:
            weight = float(cards.get("WSCIMGWG", 1.0))
        except ValueError:
            weight = 1.0
        if weight > 0:
            found.append({
                "image": str(path),
                "frequency_hz": float(cards["CRVAL3"]),
                "weight": weight,
            })
    return found


def mean_beam(shape, cell_arcmin, bands, coefficients, pbmin):
    """The weighted mean of the sub-bands' beams, and where it holds.

    Returns the beam (zero where it does not hold), the radius in arcmin
    inside which the polynomial is used, and the weighted mean frequency.
    """

    import numpy as np

    ny, nx = shape
    yy, xx = np.ogrid[:ny, :nx]
    # WSClean's phase centre is pixel (nx/2, ny/2), counting from zero.
    radius = np.hypot(xx - nx // 2, yy - ny // 2) * cell_arcmin
    total = sum(band["weight"] for band in bands)
    beam = np.zeros(shape)
    for band in bands:
        beam += band["weight"] / total * SAFETY.ugmrt_primary_beam_response(
            radius * band["frequency_hz"] / 1e9, coefficients
        )
    highest = max(band["frequency_hz"] for band in bands) / 1e9
    limit = SAFETY.first_outward_beam_cutoff(coefficients, pbmin) / highest
    holds = (radius <= limit) & (beam >= pbmin)
    frequency = sum(
        band["weight"] * band["frequency_hz"] for band in bands
    ) / total
    return np.where(holds, beam, 0.0), float(limit), float(frequency)


def flagging_masks_hz(products):
    """The frequency ranges the flagging stage masked, as (low, high) in Hz.

    What the stage recorded applying, as the pipeline's imaging stage reads
    it; without the record no mask is assumed.
    """

    if not products:
        return []
    path = Path(products) / "flagging" / "known_frequency_mask_report.json"
    try:
        with open(path) as handle:
            report = json.load(handle)
    except (OSError, ValueError):
        return []
    ranges = []
    for item in (report.get("applied_intersections") or []
                 if isinstance(report, dict) else []):
        try:
            low, high = sorted(
                float(value) * 1e6 for value in item["intersection_MHz"]
            )
        except (KeyError, TypeError, ValueError):
            continue
        if high > low:
            ranges.append([low, high])
    return sorted(ranges)


def content_identity(path):
    """Path, size and SHA-256 of a file, for the provenance record."""

    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    path = Path(path)
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size,
            "sha256": digest.hexdigest()}


def meerkat_band(band):
    """The katbeam model and science range of a MeerKAT band."""

    bands = PROFILES.TELESCOPE_PROFILES["meerkat"]["bands"]
    profile = bands.get(str(band or "").upper())
    if profile is None:
        raise SystemExit(
            f"No MeerKAT band '{band}' in telescope_profiles.py; give one of "
            f"{', '.join(sorted(bands))}."
        )
    model = profile.get("katbeam_model")
    if not model:
        raise SystemExit(
            f"MeerKAT band {band} has no measured katbeam model, so its "
            "images cannot be corrected for the primary beam here."
        )
    science = (profile.get("science_freq_min"), profile.get("science_freq_max"))
    limits = (
        (float(science[0]), float(science[1]))
        if None not in science else None
    )
    return model, limits


_KATBEAM_MODULE = None


def katbeam_module():
    """``meerkat_beam.py``, loaded once, and only for a MeerKAT image."""

    global _KATBEAM_MODULE
    if _KATBEAM_MODULE is None:
        _KATBEAM_MODULE = load_sibling(
            "meerkat_beam.py", "correct_primary_beam_katbeam")
    return _KATBEAM_MODULE


def run_meerkat(args):
    """The katbeam correction of a WSClean image of a MeerKAT field."""

    beam = katbeam_module()
    name = Path(args.name)
    image = name.with_name(name.name + "-MFS-image.fits")
    if not image.is_file():
        raise SystemExit(f"No image at {image}.")
    corrected = beam.pbcor_name(image)
    report_path = name.with_name(name.name + "-MFS-pbcor.json")
    if corrected.exists() and not args.overwrite:
        raise SystemExit(
            f"{corrected} already exists. Pass --overwrite to make it again."
        )
    model, limits = meerkat_band(args.band)
    mfs, subbands = beam.classify_outputs(
        sorted(name.parent.glob(name.name + "-*.fits")), name.name
    )
    if not subbands:
        raise SystemExit(
            f"No sub-band images named {name}-NNNN-image.fits: the beam of "
            "the image is the mean of its sub-bands' beams and cannot be "
            "had without them."
        )
    masks = flagging_masks_hz(args.products)
    record = beam.correct_products(
        mfs, subbands, model=model, pbmin=float(args.pbmin),
        science_limits_hz=limits, identity=content_identity,
        mode="band_average", masked_ranges_hz=masks,
    )
    provenance = beam.write_provenance(
        Path(record["pbcor_image"] + ".provenance.json"), record
    )
    skipped = [
        {"image": item["source_image"], "reason": item["skipped"]}
        for item in record["subbands"] if item.get("skipped")
    ]
    report = {
        "image": str(image),
        "corrected": record["pbcor_image"],
        "beam": record["pb_image"],
        "telescope": "meerkat",
        "band": str(args.band).upper(),
        "katbeam_model": model,
        "katbeam_version": record["katbeam_version"],
        "pbmin": float(args.pbmin),
        "sub_bands": len(subbands) - len(skipped),
        "skipped_sub_bands": skipped,
        "mfs_weighting": record["mfs_weighting"],
        "mean_frequency_mhz": record["mfs_effective_frequency_mhz"],
        "masked_ranges_hz": record["masked_ranges_hz"],
        "pixels_inside_pbmin": record["statistics"]["pixels_inside_pbmin"],
        "provenance": str(provenance),
    }
    with open(report_path, "w") as handle:
        json.dump(report, handle, indent=2)
    print(
        f"katbeam {model} primary beam applied: {report['sub_bands']} "
        f"sub-bands weighted as WSClean weighted them"
        + (f", {len(skipped)} without data left out" if skipped else "")
        + f", the image describing {report['mean_frequency_mhz']:.0f} MHz; "
        f"blank where the beam is below {float(args.pbmin):g}: {corrected}"
    )
    return report


def run(args):
    import numpy as np

    if getattr(args, "telescope", "ugmrt") == "meerkat":
        return run_meerkat(args)

    name = Path(args.name)
    image = name.with_name(name.name + "-MFS-image.fits")
    if not image.is_file():
        raise SystemExit(f"No image at {image}.")
    corrected = name.with_name(name.name + "-MFS-image.pbcor.fits")
    beam_path = name.with_name(name.name + "-MFS-pb.fits")
    report_path = name.with_name(name.name + "-MFS-pbcor.json")
    if corrected.exists() and not args.overwrite:
        raise SystemExit(
            f"{corrected} already exists. Pass --overwrite to make it again."
        )

    coefficients, valid = band_coefficients(args.band)
    bands = sub_bands(name)
    if not bands:
        raise SystemExit(
            f"No sub-band images named {name}-NNNN-image.fits: the beam of "
            "the image is the mean of its sub-bands' beams and cannot be "
            "had without them."
        )
    outside = [
        band for band in bands
        if not valid[0] <= band["frequency_hz"] / 1e6 <= valid[1]
    ]
    if outside:
        raise SystemExit(
            f"{len(outside)} sub-bands lie outside {valid[0]:g} to "
            f"{valid[1]:g} MHz, where the beam of band {args.band} was "
            "measured."
        )

    header, cards, plane, tail = read_fits(image)
    cell_arcmin = abs(float(cards["CDELT2"])) * 60.0
    beam, limit, frequency = mean_beam(
        plane.shape, cell_arcmin, bands, coefficients, float(args.pbmin)
    )
    holds = beam > 0
    result = np.full(plane.shape, np.nan)
    result[holds] = plane[holds] / beam[holds]
    write_fits(corrected, header, result, tail)
    write_fits(beam_path, header, np.where(holds, beam, np.nan), tail)

    report = {
        "image": str(image),
        "corrected": str(corrected),
        "beam": str(beam_path),
        "band": args.band,
        "coefficients": coefficients,
        "pbmin": float(args.pbmin),
        "sub_bands": len(bands),
        "mean_frequency_mhz": frequency / 1e6,
        "radius_arcmin": limit,
        "corrected_fraction_of_image": float(holds.mean()),
        "beam_at_edge": float(beam[holds].min()),
    }
    with open(report_path, "w") as handle:
        json.dump(report, handle, indent=2)
    print(
        f"Primary-beam corrected out to {limit:.1f} arcmin (beam down to "
        f"{report['beam_at_edge']:.2f}), the mean of {len(bands)} sub-bands' "
        f"beams about {frequency / 1e6:.0f} MHz: {corrected}"
    )
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--name", required=True,
                        help="the name WSClean was given, as in "
                        "<name>-MFS-image.fits")
    parser.add_argument("--telescope", default="ugmrt",
                        choices=("ugmrt", "meerkat"))
    parser.add_argument("--band",
                        help="the band, as in telescope_profiles.py "
                        "(default B4 for uGMRT; MEERKAT_L or MEERKAT_UHF "
                        "for MeerKAT)")
    parser.add_argument("--products",
                        help="the pipeline's output folder, read for the "
                        "frequency ranges its flagging masked (MeerKAT)")
    parser.add_argument("--pbmin", type=float, default=0.3,
                        help="the corrected image is blank where the beam "
                        "is below this (default 0.3, the pipeline's)")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    if args.band is None:
        if args.telescope == "meerkat":
            parser.error("--band is needed for MeerKAT: MEERKAT_L or "
                         "MEERKAT_UHF")
        args.band = "B4"
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
