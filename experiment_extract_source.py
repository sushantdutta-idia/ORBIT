#!/usr/bin/env python3
"""Leave one source in the data: subtract the rest of the field around it.

The test of a direction-dependent error. Every solve so far found one gain
per antenna for the whole field. If the antennas do not all see an
off-axis source alike (they point a little differently, and on the flank
of the primary beam a little is a lot), that one gain is right for none of
the sources out there, and each keeps holes beside it that grow with
frequency and with distance from the pointing centre. On G0212-0335 that
is what is left (finding 123).

Gains towards one source can only be solved where that source is alone.
This script makes such a Measurement Set:

1. from the model images of a WSClean run, one per sub-band, it makes a
   copy with a disc around the source set to zero: the rest of the field;
2. WSClean predicts that into the model column of a copy of the data;
3. CASA subtracts it, keeps the result, and moves the phase centre to the
   source.

What comes out holds the source, the noise, and whatever of the other
sources their model did not describe. ``selfcal_recipe.py`` on it, with a
small image, is then self-calibration in that one direction; if its gains
differ from antenna to antenna and its image loses the holes, the error is
direction-dependent.

``--peel`` is the way back, once ``selfcal_recipe.py`` has calibrated the
source alone. The source's own model, made from its calibrated data, is
taken out of the data the recipe started from, as those data see it: each
visibility of the model times the gain the recipe found towards the
source, which is the recipe's input over its output. What is left goes
back to the field's phase centre, and the model of the rest of the field
is added to it again. The result is the field without that source, taken
out with the gains the antennas had towards it, which is what a
subtraction with the field's one gain per antenna cannot do; everything
else keeps the field's own calibration and weights. The next source is
then taken from that, with ``--also-without`` naming the ones already
gone.

WSClean and CASA are run with the commands given, as in
``selfcal_recipe.py``; this script needs numpy for the model images. The
Measurement Set it is given is copied, never written to.

Example::

    python3 experiment_extract_source.py \\
        --ms wsclean_test/G0212-0335_scan15_ph_r2_r3_r4_r5_bpa.ms \\
        --models wsclean_test_scan15_bpa_ch12/G0212-0335/wsclean_robust_0_I \\
        --pixel 3162 4227 --name source1 --channel-range 51 768
    python3 experiment_extract_source.py --peel recipe_source1/G0212-0335 \
        --models wsclean_test_scan15_bpa_ch12/G0212-0335/wsclean_robust_0_I \
        --pixel 3162 4227 --name source1 --channel-range 51 768
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
import math
import shlex
import shutil
from pathlib import Path


CONTAINERS = "/software/astro/containers"
DEFAULT_WSCLEAN = f"singularity exec {CONTAINERS}/oxkat-latest.sif wsclean"
DEFAULT_CASA = (
    f"singularity exec {CONTAINERS}/casa-modular-py3.12-latest.sif python3"
)
BLOCK = 2880


def load_imager():
    spec = importlib.util.spec_from_file_location(
        "extract_source_imager",
        Path(SCRIPT_DIR) / "experiment_wsclean_image.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Model images
# ---------------------------------------------------------------------------

def read_fits(path):
    """A WSClean image as (header bytes, header cards, plane [y, x], tail).

    The plane is a copy that may be changed; ``tail`` is the padding after
    the data, so that the three parts written back give a file of the same
    size.
    """

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
    count = nx * ny
    plane = np.frombuffer(
        raw, dtype=">f4", count=count, offset=offset
    ).reshape(ny, nx).copy()
    return raw[:offset], cards, plane, raw[offset + 4 * count:]


def write_fits(path, header, plane, tail):
    with open(path, "wb") as handle:
        handle.write(header)
        handle.write(plane.astype(">f4").tobytes())
        handle.write(tail)


def model_images(prefix):
    """The per-sub-band model images of a WSClean run, in order."""

    prefix = Path(prefix)
    found = sorted(prefix.parent.glob(prefix.name + "-[0-9][0-9][0-9][0-9]-model.fits"))
    if not found:
        raise SystemExit(
            f"No model images named {prefix}-NNNN-model.fits. --models is "
            "the name WSClean was given, without the sub-band number."
        )
    return found


def pixel_direction(cards, x, y):
    """Right ascension and declination, in degrees, of a zero-based pixel.

    WSClean's images are in the sine projection about the phase centre.
    """

    ra0 = math.radians(float(cards["CRVAL1"]))
    dec0 = math.radians(float(cards["CRVAL2"]))
    l = math.radians((x + 1 - float(cards["CRPIX1"])) * float(cards["CDELT1"]))
    m = math.radians((y + 1 - float(cards["CRPIX2"])) * float(cards["CDELT2"]))
    n = math.sqrt(max(0.0, 1.0 - l * l - m * m))
    ra = ra0 + math.atan2(l, n * math.cos(dec0) - m * math.sin(dec0))
    dec = math.asin(m * math.cos(dec0) + n * math.sin(dec0))
    return math.degrees(ra) % 360.0, math.degrees(dec)


def without_disc(plane, x, y, radius):
    """A copy of the plane with a disc set to zero, and the flux taken out."""

    import numpy as np

    ny, nx = plane.shape
    y0, y1 = max(0, y - radius), min(ny, y + radius + 1)
    x0, x1 = max(0, x - radius), min(nx, x + radius + 1)
    yy, xx = np.ogrid[y0:y1, x0:x1]
    disc = (xx - x) ** 2 + (yy - y) ** 2 <= radius * radius
    rest = plane.copy()
    block = rest[y0:y1, x0:x1]
    taken = float(block[disc].sum())
    block[disc] = 0.0
    return rest, taken


def split_models(images, folder, x, y, radius, others=()):
    """Write the rest-of-the-field models; return what each sub-band lost.

    ``others`` are sources already taken out of the data: their discs are
    left out of the rest as well, so that nothing is subtracted, or added
    back, where there is no longer anything.
    """

    cards, table = None, []
    for index, image in enumerate(images):
        header, these, plane, tail = read_fits(image)
        if cards is None:
            cards = these
        ny, nx = plane.shape
        # Each source already out has its own disc; one given as a pixel
        # alone has this source's.
        discs = [
            (int(other[0]), int(other[1]),
             int(other[2]) if len(other) > 2 else int(radius))
            for other in others
        ]
        for px, py, _ in [(x, y, radius), *discs]:
            if not (0 <= px < nx and 0 <= py < ny):
                raise SystemExit(
                    f"Pixel {px}, {py} is outside the {nx} x {ny} model "
                    "images."
                )
        rest, taken = without_disc(plane, x, y, radius)
        for px, py, pr in discs:
            rest, _ = without_disc(rest, px, py, pr)
        write_fits(folder / f"rest-{index:04d}-model.fits", header, rest, tail)
        table.append({
            "sub_band": index,
            "source_model_Jy": taken,
            "rest_model_Jy": float(rest.sum()),
        })
    return cards, table


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def predict_command(args, ms, prefix, count, capabilities, imager):
    """WSClean writing the rest-of-the-field model into MODEL_DATA."""

    command = shlex.split(args.wsclean) + [
        "-predict",
        "-name", str(prefix),
        "-channels-out", str(int(count)),
        "-j", str(int(args.threads)),
        "-abs-mem", f"{float(args.memory_gb):g}",
        *imager.gridder_arguments(capabilities),
    ]
    if args.channel_range:
        command += [
            "-channel-range",
            str(int(args.channel_range[0])), str(int(args.channel_range[1])),
        ]
    return command + [str(ms)]


def casa_command(args, field_ms, out_ms, direction=None,
                 subtracted_ms=None, original_ms=None):
    """This script again, under CASA.

    With a direction: subtract the model column, keep the result, and move
    the phase centre there. Without one: add the model column and keep the
    result. With ``original_ms``, the data ``field_ms`` was calibrated
    from: the model is subtracted from those, as they see it.
    """

    command = shlex.split(args.casa) + [
        str(Path(SCRIPT_DIR) / "experiment_extract_source.py"),
        "--casa-part", "subtract-shift" if direction else "add",
        "--field-ms", str(field_ms),
        "--out-ms", str(out_ms),
    ]
    if direction:
        command += [
            "--subtracted-ms", str(subtracted_ms),
            "--phasecenter",
            f"J2000 {direction[0]:.8f}deg {direction[1]:.8f}deg",
        ]
        if original_ms is not None:
            command += ["--original-ms", str(original_ms)]
    return command


def subtract_as_calibrated(calibrated, original, out, chunk=20000):
    """The original data less the model as those data see it.

    ``calibrated`` holds a calibrated copy of ``original``, row for row,
    with the model in its model column. Calibration divides each
    visibility by a gain, so the original over the calibrated is that
    gain, exactly: what the model looks like through it is model times
    gain. Peeling with the calibrated data instead handed the whole field
    back with the gains solved towards one source: on G0216-0604 in Band 5
    a 10 mJy/beam source at 10' had amplitudes of 0.44 to 1.41, moving by
    56% within a scan, and the field's Stokes V noise fell from 16.5 to
    13.4 uJy/beam with them (2026-10-23). Where the gain cannot be had
    (the calibrated visibility is zero) the model is taken as it is.
    Returns the part of the visibilities without a gain and the median
    gain amplitude.
    """

    import numpy as np
    from casatools import table

    calibrated, original, out = Path(calibrated), Path(original), Path(out)
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(original, out)
    source, target = table(), table()
    source.open(str(calibrated))
    target.open(str(out), nomodify=False)
    try:
        rows = source.nrows()
        if target.nrows() != rows:
            raise SystemExit(
                f"{calibrated} has {rows} rows and {original} "
                f"{target.nrows()}: they are not the same data."
            )
        for column in ("ANTENNA1", "ANTENNA2"):
            if not np.array_equal(source.getcol(column),
                                  target.getcol(column)):
                raise SystemExit(
                    f"{calibrated} and {original} differ in {column}: they "
                    "are not the same data, row for row."
                )
        if not np.allclose(source.getcol("TIME"), target.getcol("TIME"),
                           rtol=0, atol=1e-3):
            raise SystemExit(
                f"{calibrated} and {original} differ in TIME: they are not "
                "the same data, row for row."
            )
        undefined, total, sizes = 0, 0, []
        for start in range(0, rows, chunk):
            count = min(chunk, rows - start)
            seen = source.getcol("DATA", start, count)
            model = source.getcol("MODEL_DATA", start, count)
            data = target.getcol("DATA", start, count)
            if seen.shape != data.shape:
                raise SystemExit(
                    f"{calibrated} has data of shape {seen.shape[:2]} and "
                    f"{original} {data.shape[:2]}: they are not the same "
                    "data."
                )
            with np.errstate(divide="ignore", invalid="ignore"):
                gain = data / seen
            missing = (seen == 0) | ~np.isfinite(gain)
            gain[missing] = 1.0
            target.putcol("DATA", data - model * gain, start, count)
            undefined += int(missing.sum())
            total += missing.size
            present = np.abs(gain[~missing & (data != 0)])
            if present.size:
                sizes.append(present[::max(1, present.size // 2000)])
    finally:
        source.close()
        target.close()
    median = float(np.median(np.concatenate(sizes))) if sizes else None
    return (undefined / total if total else 0.0), median


def casa_part(args):
    """The steps that need CASA, on a Measurement Set this run owns."""

    from casatasks import split, uvsub

    field, out = Path(args.field_ms), Path(args.out_ms)
    if out.exists():
        shutil.rmtree(out)
    if args.casa_part == "add":
        print(f"Adding the model column of {field} to its data ...",
              flush=True)
        uvsub(vis=str(field), reverse=True)
        split(vis=str(field), outputvis=str(out), datacolumn="corrected")
    else:
        from casatasks import phaseshift

        subtracted = Path(args.subtracted_ms)
        if subtracted.exists():
            shutil.rmtree(subtracted)
        if args.original_ms:
            print(f"Subtracting the model column of {field} from "
                  f"{args.original_ms}, through the gains between them ...",
                  flush=True)
            undefined, median = subtract_as_calibrated(
                field, args.original_ms, subtracted)
            print(f"Gains towards the source: median amplitude "
                  f"{median if median is None else round(median, 3)}; "
                  f"{100 * undefined:.2f}% of the visibilities had none and "
                  "took the model as it is.", flush=True)
        else:
            print(f"Subtracting the model column of {field} from its data "
                  "...", flush=True)
            uvsub(vis=str(field))
            split(vis=str(field), outputvis=str(subtracted),
                  datacolumn="corrected")
        print(f"Moving the phase centre to {args.phasecenter} ...",
              flush=True)
        phaseshift(
            vis=str(subtracted), outputvis=str(out),
            phasecenter=args.phasecenter,
        )
    if not out.is_dir():
        raise SystemExit("No Measurement Set was written.")
    print(f"Written: {out}")
    return 0


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

def sources_out(args):
    """``--also-without``: each source already out, with its disc."""

    out = []
    for values in args.also_without or []:
        if len(values) not in (2, 3):
            raise SystemExit(
                "--also-without takes X Y and, if its disc was not this "
                f"one's, RADIUS; got {' '.join(str(v) for v in values)}."
            )
        x, y = int(values[0]), int(values[1])
        out.append((x, y, int(values[2]) if len(values) == 3
                    else int(args.radius)))
    return out


def tools(args, imager):
    """Stop now if either tool does not answer; return what WSClean has."""

    for name, command in (("WSClean", args.wsclean), ("CASA", args.casa)):
        words = shlex.split(command)
        if not words or shutil.which(words[0]) is None:
            raise SystemExit(
                f"{name}: '{words[0] if words else ''}' is not on PATH. This "
                "script is started from the terminal that has singularity."
            )
    capabilities = imager.wsclean_capabilities(shlex.split(args.wsclean))
    if capabilities["version"] is None and not capabilities["options"]:
        raise SystemExit(
            f"'{args.wsclean} --version' did not answer as WSClean."
        )
    return capabilities


def recipe_products(recipe):
    """What ``selfcal_recipe.py`` left for a source: data and final model."""

    recipe = Path(recipe)
    try:
        with open(recipe / "recipe_report.json") as handle:
            report = json.load(handle)
    except (OSError, ValueError) as exc:
        raise SystemExit(
            f"No recipe report in {recipe} ({exc}). --peel is the folder "
            "selfcal_recipe.py worked in, with the target's name."
        ) from exc
    if "final:image" not in report.get("done", []):
        raise SystemExit(
            f"The recipe in {recipe} has not finished: its final image is "
            "the model that is subtracted."
        )
    calibrated = Path(report["calibrated_measurement_set"])
    if not calibrated.is_dir():
        raise SystemExit(f"The calibrated data {calibrated} are not there.")
    first = sorted((recipe / "images" / "final").glob("*/*_I-0000-model.fits"))
    if not first:
        raise SystemExit(
            f"No final model images under {recipe / 'images' / 'final'}."
        )
    prefix = first[0].with_name(first[0].name[:-len("-0000-model.fits")])
    # The data the recipe started from, if they are still there.
    original = Path(report["input"]) if report.get("input") else None
    if original is not None and not original.is_dir():
        original = None
    return calibrated, prefix, original


def run_peel(args):
    """Take a calibrated source out, and give the field back without it."""

    imager = load_imager()
    calibrated, source_prefix, original = recipe_products(args.peel)
    if original is None:
        print(
            "WARNING: the data the recipe started from are gone; the field "
            "is given back with the gains solved towards the source.",
            flush=True,
        )
    folder = Path(args.workdir) / args.name
    if imager.inside(folder, args.outdir):
        raise SystemExit(
            f"{folder} is inside the pipeline's products ({args.outdir}). "
            "Give --workdir a folder of its own."
        )
    images = model_images(args.models)
    source_images = model_images(source_prefix)
    if len(source_images) != len(images):
        raise SystemExit(
            f"The source's model has {len(source_images)} sub-bands and the "
            f"field's {len(images)}: they have to be the same."
        )
    x, y = int(args.pixel[0]), int(args.pixel[1])
    others = sources_out(args)
    work_ms = folder / "calibrated.ms"
    at_source_ms = folder / "residual_at_source.ms"
    residual_ms = folder / "residual.ms"
    out_ms = folder / f"{args.name}_peeled.ms"
    report_path = folder / "peel.json"
    if out_ms.exists() and not args.overwrite:
        raise SystemExit(
            f"{out_ms} already exists. Pass --overwrite to make it again."
        )
    capabilities = tools(args, imager)

    folder.mkdir(parents=True, exist_ok=True)
    with imager.claimed(folder):
        for path in (work_ms, at_source_ms, residual_ms, out_ms):
            if path.exists():
                shutil.rmtree(path)
        for path in folder.glob("rest-*-model.fits"):
            path.unlink()

        print(f"Copying {calibrated} to {work_ms} ...", flush=True)
        shutil.copytree(calibrated, work_ms)

        print("\nWSClean: the source's own model into the model column ...",
              flush=True)
        imager.run_command(
            predict_command(
                args, work_ms, source_prefix, len(images), capabilities,
                imager,
            ),
            folder / "predict_source.log",
        )

        cards, table = split_models(
            images, folder, x, y, int(args.radius), others
        )
        centre = (float(cards["CRVAL1"]) % 360.0, float(cards["CRVAL2"]))
        print("\nCASA: subtracting the source and going back to the field's "
              "centre ...", flush=True)
        imager.run_command(
            casa_command(args, work_ms, residual_ms, centre, at_source_ms,
                         original),
            folder / "casa_subtract.log",
        )
        if not residual_ms.is_dir():
            raise SystemExit(
                "No Measurement Set was written; see "
                f"{folder / 'casa_subtract.log'}."
            )

        print("\nWSClean: the rest of the field into the model column ...",
              flush=True)
        imager.run_command(
            predict_command(
                args, residual_ms, folder / "rest", len(images), capabilities,
                imager,
            ),
            folder / "predict_rest.log",
        )
        print("\nCASA: adding the rest of the field back ...", flush=True)
        imager.run_command(
            casa_command(args, residual_ms, out_ms),
            folder / "casa_add.log",
        )
        if not out_ms.is_dir():
            raise SystemExit(
                f"No Measurement Set was written; see {folder / 'casa_add.log'}."
            )
        if not args.keep_intermediate:
            for path in (work_ms, at_source_ms, residual_ms):
                shutil.rmtree(path, ignore_errors=True)
            for path in folder.glob("rest-*-model.fits"):
                path.unlink()

        report = {
            "recipe": str(args.peel),
            "calibrated": str(calibrated),
            "source_model": str(source_prefix),
            "field_models": [str(image) for image in images],
            "pixel": [x, y],
            "already_out": [list(other) for other in others],
            "radius_pixels": int(args.radius),
            "field_centre_deg": list(centre),
            "sub_bands": table,
            "measurement_set": str(out_ms),
            # Whose calibration the rest of the field is given back with.
            "calibration_kept": (
                "the field's" if original is not None
                else "the source's (its recipe's input was gone)"
            ),
            "original": None if original is None else str(original),
        }
        with open(report_path, "w") as handle:
            json.dump(report, handle, indent=2)
    print(f"\nThe field without the source: {out_ms}\n"
          f"Report written: {report_path}")
    return report


def run(args):
    imager = load_imager()
    source = Path(str(args.ms).rstrip("/"))
    if not source.is_dir():
        raise SystemExit(f"No Measurement Set at {source}.")
    folder = Path(args.workdir) / args.name
    if imager.inside(folder, args.outdir):
        raise SystemExit(
            f"{folder} is inside the pipeline's products ({args.outdir}). "
            "Give --workdir a folder of its own."
        )
    images = model_images(args.models)
    x, y = int(args.pixel[0]), int(args.pixel[1])
    field_ms = folder / "field.ms"
    subtracted_ms = folder / "subtracted.ms"
    out_ms = folder / f"{args.name}.ms"
    report_path = folder / "extract.json"

    if out_ms.exists() and not args.overwrite:
        raise SystemExit(
            f"{out_ms} already exists. Pass --overwrite to make it again, "
            "or --name for another source."
        )
    capabilities = tools(args, imager)
    others = sources_out(args)

    folder.mkdir(parents=True, exist_ok=True)
    with imager.claimed(folder):
        for path in (field_ms, subtracted_ms, out_ms):
            if path.exists():
                shutil.rmtree(path)
        for path in folder.glob("rest-*-model.fits"):
            path.unlink()

        print(
            f"Taking a disc of {args.radius} pixels about {x}, {y} out of "
            f"{len(images)} model images ...",
            flush=True,
        )
        cards, table = split_models(
            images, folder, x, y, int(args.radius), others
        )
        direction = pixel_direction(cards, x, y)
        kept = sum(row["source_model_Jy"] for row in table) / len(table)
        rest = sum(row["rest_model_Jy"] for row in table) / len(table)
        print(
            f"  the source's model holds {kept * 1e3:.1f} mJy on average, "
            f"the rest of the field {rest * 1e3:.0f} mJy; the source is at "
            f"RA {direction[0]:.5f} Dec {direction[1]:+.5f} deg.",
            flush=True,
        )

        print(f"\nCopying {source} to {field_ms} ...", flush=True)
        shutil.copytree(source, field_ms)

        print("\nWSClean: the rest of the field into the model column ...",
              flush=True)
        imager.run_command(
            predict_command(
                args, field_ms, folder / "rest", len(images), capabilities,
                imager,
            ),
            folder / "predict.log",
        )

        print("\nCASA: subtracting it and centring on the source ...",
              flush=True)
        imager.run_command(
            casa_command(args, field_ms, out_ms, direction, subtracted_ms),
            folder / "casa.log",
        )
        if not out_ms.is_dir():
            raise SystemExit(
                f"No Measurement Set was written; see {folder / 'casa.log'}."
            )
        if not args.keep_intermediate:
            # The copy with its model and corrected columns, and the
            # subtracted data before the shift, have served.
            for path in (field_ms, subtracted_ms):
                shutil.rmtree(path, ignore_errors=True)
            for path in folder.glob("rest-*-model.fits"):
                path.unlink()

        report = {
            "input": str(source),
            "models": [str(image) for image in images],
            "pixel": [x, y],
            "already_out": [list(other) for other in others],
            "radius_pixels": int(args.radius),
            "ra_deg": direction[0],
            "dec_deg": direction[1],
            "sub_bands": table,
            "measurement_set": str(out_ms),
        }
        with open(report_path, "w") as handle:
            json.dump(report, handle, indent=2)
    print(f"\nThe source alone: {out_ms}\nReport written: {report_path}")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--casa-part", choices=("subtract-shift", "add"),
                        help=argparse.SUPPRESS)
    parser.add_argument("--field-ms", help=argparse.SUPPRESS)
    parser.add_argument("--subtracted-ms", help=argparse.SUPPRESS)
    parser.add_argument("--out-ms", help=argparse.SUPPRESS)
    parser.add_argument("--phasecenter", help=argparse.SUPPRESS)
    parser.add_argument("--original-ms", help=argparse.SUPPRESS)
    parser.add_argument("--ms", help="the calibrated data; copied, never "
                        "written to")
    parser.add_argument("--models", help="the name WSClean was given for "
                        "the run whose model images to use, as in "
                        "<name>-0000-model.fits")
    parser.add_argument("--pixel", type=int, nargs=2, metavar=("X", "Y"),
                        help="the source, in zero-based pixels of those "
                        "images")
    parser.add_argument("--also-without", type=int, nargs="+",
                        action="append", metavar="X Y [RADIUS]",
                        help="a source already taken out of these data, "
                        "and the disc it was taken out with (this one's "
                        "--radius if not given); may be given several "
                        "times")
    parser.add_argument("--peel", metavar="RECIPE_FOLDER",
                        help="instead of leaving the source alone, take it "
                        "out: the folder (with the target's name) in which "
                        "selfcal_recipe.py calibrated it, and write the "
                        "field without it")
    parser.add_argument("--radius", type=int, default=90,
                        help="pixels about the source that are left in the "
                        "data (default 90)")
    parser.add_argument("--name", default="source",
                        help="names the folder and the Measurement Set")
    parser.add_argument("--workdir", default="extract_test")
    parser.add_argument("--outdir", default="uGMRT_Products",
                        help="the pipeline's products, which are not touched")
    parser.add_argument("--channel-range", type=int, nargs=2,
                        metavar=("FIRST", "END"),
                        help="the channels the models were made from, END "
                        "exclusive")
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--memory-gb", type=float, default=20.0)
    parser.add_argument("--wsclean", default=DEFAULT_WSCLEAN)
    parser.add_argument("--casa", default=DEFAULT_CASA)
    parser.add_argument("--keep-intermediate", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    if args.casa_part:
        return casa_part(args)
    for needed in ("models", "pixel") + (() if args.peel else ("ms",)):
        if getattr(args, needed) is None:
            parser.error(f"--{needed} is required")
    if args.peel:
        run_peel(args)
    else:
        run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
