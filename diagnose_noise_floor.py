#!/usr/bin/env python3
"""Measure what limits a target's image: the data, or what was made of them.

A self-calibrated image that sits well above its thermal noise can be limited
by two different things, and they call for different work:

* the visibilities themselves are noisier than their weights say (interference
  left in the data, a bad scan or antenna), in which case no calibration or
  deconvolution reaches the thermal figure and the flagging has to change;
* the visibilities are sound, and the excess is sky that was not removed
  properly: sidelobes of sources calibrated or deconvolved imperfectly.

Three measurements separate them, and they take minutes, not hours.

1. Stokes V images. The radio sky is almost unpolarized in circular
   polarization, so RR - LL holds the noise of the data and next to nothing
   of the sky. What it does hold is leakage: the R and L beams of an antenna
   do not coincide, so a bright source far from the pointing centre appears
   in V at several per cent of its Stokes I, and in a dirty image its
   sidelobes raise the rms. The dirty image is therefore made first, then
   cleaned shallowly to take those out; the residual's rms is the floor a
   perfect Stokes I image of these data could reach with that weighting.

2. The point-spread function for each weighting asked for: its width, its
   largest sidelobe and the rms of its far sidelobes. Together with the
   Stokes V floor this is what a change of ``robust`` buys and costs, since
   errors on a source spread through the image at the level of the
   sidelobes while the noise per beam rises as the weighting moves towards
   uniform.

3. The noise of the existing Stokes I image, measured in a grid of boxes.
   Where the image is quiet far from bright sources and noisy around them,
   the limit is those sources (direction-dependent errors: peel them, or
   model them better). Where it is uniformly high, the limit is everywhere.

4. With ``--model-residual``, the Stokes I residual of the model the
   pipeline already made: one major cycle that subtracts the model from the
   data and images what is left, for each weighting and for each scan. No
   CLEAN, so it takes minutes. Set against the Stokes V of the same data it
   says how much of each scan is not noise, which is where calibration
   errors show: a scan that carries them has a residual far above its own
   floor while the others sit near theirs.

Nothing here modifies the Measurement Set or any pipeline product. The images
and the report go to ``--workdir``.

Examples::

    python3 diagnose_noise_floor.py --outdir uGMRT_Products --target G0212-0335
    python3 diagnose_noise_floor.py --outdir uGMRT_Products \\
        --target G0212-0335 --robust 0.5 0 -0.3 --skip-scans
    python3 diagnose_noise_floor.py --outdir uGMRT_Products \\
        --target G0212-0335 --robust 0 0.5 --model-residual
"""

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# This folder holds the pipeline's own inspect.py. Run from here, anything
# that imports the standard library's inspect (numpy and dataclasses do)
# would get that file instead.
if sys.path and os.path.abspath(sys.path[0]) == SCRIPT_DIR:
    sys.path.pop(0)

import argparse
import json
import math
import shutil
from pathlib import Path


ARCSEC_PER_RADIAN = 206264.80624709636


# ---------------------------------------------------------------------------
# Statistics (numpy only)
# ---------------------------------------------------------------------------

def robust_rms(values):
    """Noise of an image region from its median absolute deviation.

    Sources and artefacts are outliers of a noise distribution, and a plain
    rms counts them as noise. Blanked pixels (not finite, or exactly zero as
    CASA writes them outside a mask) are left out. Returns None when fewer
    than ten pixels remain.
    """

    import numpy as np

    data = np.asarray(values, dtype=float).reshape(-1)
    data = data[np.isfinite(data) & (data != 0.0)]
    if data.size < 10:
        return None
    median = np.median(data)
    return float(1.4826 * np.median(np.abs(data - median)))


def plain_rms(values):
    import numpy as np

    data = np.asarray(values, dtype=float).reshape(-1)
    data = data[np.isfinite(data) & (data != 0.0)]
    if data.size < 1:
        return None
    return float(np.sqrt(np.mean(data * data)))


def rms_grid(plane, cells=10):
    """Robust rms of a 2-D image in ``cells`` x ``cells`` boxes.

    Returned as rows of the first image axis (right ascension pixel), each a
    list along the second (declination pixel).
    """

    import numpy as np

    plane = np.asarray(plane, dtype=float)
    if plane.ndim != 2:
        raise ValueError("rms_grid needs a two-dimensional image plane")
    nx, ny = plane.shape
    x_edges = np.linspace(0, nx, cells + 1).astype(int)
    y_edges = np.linspace(0, ny, cells + 1).astype(int)
    grid = []
    for i in range(cells):
        row = []
        for j in range(cells):
            box = plane[x_edges[i]:x_edges[i + 1], y_edges[j]:y_edges[j + 1]]
            row.append(robust_rms(box))
        grid.append(row)
    return grid


def summarize_grid(grid, shape):
    """Where an image is quiet, where it is not, and by how much."""

    import numpy as np

    cells = len(grid)
    entries = [
        (value, i, j)
        for i, row in enumerate(grid)
        for j, value in enumerate(row)
        if value is not None and value > 0
    ]
    if not entries:
        return None
    values = np.asarray([entry[0] for entry in entries])
    quiet = float(np.percentile(values, 10))
    median = float(np.median(values))

    def centre(i, j):
        return [
            int((i + 0.5) * shape[0] / cells),
            int((j + 0.5) * shape[1] / cells),
        ]

    noisiest = sorted(entries, reverse=True)[:5]
    return {
        "boxes": len(entries),
        "quietest_tenth_rms": quiet,
        "median_rms": median,
        "noisiest_rms": float(values.max()),
        "noisiest_over_quietest": float(values.max() / quiet) if quiet else None,
        "median_over_quietest": median / quiet if quiet else None,
        "noisiest_boxes": [
            {"rms": float(value), "centre_pixel": centre(i, j)}
            for value, i, j in noisiest
        ],
    }


def brightest_pixel(plane):
    import numpy as np

    plane = np.asarray(plane, dtype=float)
    finite = np.where(np.isfinite(plane), plane, -np.inf)
    index = np.unravel_index(int(np.argmax(finite)), plane.shape)
    return {
        "pixel": [int(index[0]), int(index[1])],
        "value": float(plane[index]),
    }


def psf_statistics(plane, cell_arcsec=1.0):
    """Width and sidelobes of a point-spread function.

    The width is that of a circle with the area the beam has above half its
    peak. Sidelobes are everything beyond one and a half widths from the
    peak, where a Gaussian of that width has fallen to 0.2%: the largest of
    them is what spreads a source's errors over its neighbourhood. The far
    sidelobes, beyond twenty widths, spread them over the whole image, and
    their rms is the factor by which an undeconvolved source of flux S
    raises the noise everywhere (by S times that rms).
    """

    import numpy as np

    plane = np.asarray(plane, dtype=float)
    peak_index = np.unravel_index(int(np.argmax(plane)), plane.shape)
    peak = float(plane[peak_index])
    if not peak > 0:
        return None
    normalised = plane / peak

    x = np.arange(plane.shape[0])[:, None] - peak_index[0]
    y = np.arange(plane.shape[1])[None, :] - peak_index[1]
    radius = np.hypot(x, y)

    # Only the lobe the peak belongs to: its extent is searched outwards.
    window = radius <= 100
    half = int(np.count_nonzero(window & (normalised >= 0.5)))
    width = 2.0 * math.sqrt(half / math.pi)
    outside = radius > 1.5 * width
    far = radius > 20.0 * width
    return {
        "width_pixels": width,
        "width_arcsec": width * float(cell_arcsec),
        "largest_sidelobe": (
            float(np.abs(normalised[outside]).max()) if outside.any() else None
        ),
        "far_sidelobe_rms": (
            float(np.sqrt(np.mean(normalised[far] ** 2))) if far.any() else None
        ),
    }


def weighting_floor(entry):
    """The data's noise for a weighting: cleaned Stokes V, else dirty."""

    cleaned = (entry.get("stokes_v_cleaned") or {}).get("robust_rms_Jy_beam")
    dirty = (entry.get("stokes_v_dirty") or {}).get("robust_rms_Jy_beam")
    return cleaned if cleaned else dirty


def scan_table(v_scans, i_scans):
    """Per scan, the noise and the Stokes I residual side by side."""

    table = {}
    for scan in sorted(
        set(v_scans or {}) | set(i_scans or {}),
        key=lambda value: (len(value), value),
    ):
        if scan == "robust":
            continue
        table[scan] = {
            "stokes_v": ((v_scans or {}).get(scan) or {}).get(
                "robust_rms_Jy_beam"
            ),
            "stokes_i": ((i_scans or {}).get(scan) or {}).get(
                "robust_rms_Jy_beam"
            ),
        }
    return table


def reading(weightings, image_summary, goal, scans=None):
    """Say in words what the measurements imply for the goal."""

    lines = []
    for entry in weightings:
        cleaned = (entry.get("stokes_v_cleaned") or {}).get(
            "robust_rms_Jy_beam"
        )
        floor = weighting_floor(entry)
        if floor is None:
            continue
        kind = "cleaned" if cleaned else "dirty, so an upper limit"
        psf = entry.get("psf") or {}
        detail = ""
        if psf.get("width_arcsec"):
            detail = (
                f" Beam {psf['width_arcsec']:.1f} arcsec, largest sidelobe "
                f"{psf['largest_sidelobe']:.2f}, far-sidelobe rms "
                f"{psf['far_sidelobe_rms']:.4f}."
            )
        lines.append(
            f"robust {entry['robust']:g}: the data reach "
            f"{floor * 1e6:.1f} uJy/beam in Stokes V ({kind}), "
            + (
                "at or below" if floor <= goal else "above"
            )
            + f" the {goal * 1e6:.0f} uJy/beam goal."
            + detail
        )
    for entry in weightings:
        residual = (entry.get("stokes_i_model_residual") or {}).get(
            "robust_rms_Jy_beam"
        )
        floor = weighting_floor(entry)
        if residual and floor:
            lines.append(
                f"robust {entry['robust']:g}: the Stokes I residual of the "
                f"existing model is {residual * 1e6:.1f} uJy/beam, "
                f"{residual / floor:.1f} times the data's noise."
            )
    for scan, values in (scans or {}).items():
        noise, residual = values.get("stokes_v"), values.get("stokes_i")
        if noise and residual:
            lines.append(
                f"scan {scan}: Stokes I residual {residual * 1e6:.1f} "
                f"against {noise * 1e6:.1f} uJy/beam of noise "
                f"({residual / noise:.1f} times)."
            )
    if image_summary:
        ratio = image_summary.get("noisiest_over_quietest") or 0.0
        quiet = image_summary["quietest_tenth_rms"]
        if ratio >= 2.0:
            lines.append(
                f"The Stokes I image is {ratio:.1f} times noisier in its "
                f"worst box than in its quietest tenth ({quiet * 1e6:.1f} "
                "uJy/beam): the excess is concentrated, which is what "
                "errors towards particular sources look like."
            )
        else:
            lines.append(
                "The Stokes I noise is nearly uniform across the image "
                f"(quietest tenth {quiet * 1e6:.1f} uJy/beam, worst box "
                f"{ratio:.1f} times that): the excess is not tied to one "
                "region."
            )
    return lines


# ---------------------------------------------------------------------------
# Locating the products
# ---------------------------------------------------------------------------

def resolve_products(outdir, target):
    """The self-calibrated Measurement Set and image of a target."""

    selfcal = Path(outdir) / "targets" / target / "selfcal"
    ms_dir = selfcal / "selfcal_ms"
    candidates = [ms_dir / "final_selfcal.ms"] + sorted(
        ms_dir.glob("phase[0-9]*.ms"),
        key=lambda path: int("".join(ch for ch in path.stem if ch.isdigit())),
        reverse=True,
    )
    ms = next((path for path in candidates if path.is_dir()), None)

    image = None
    for name in ("final", "amp_phase"):
        candidate = selfcal / "images" / f"{name}.image.tt0"
        if candidate.is_dir():
            image = candidate
            break
    return ms, image


# ---------------------------------------------------------------------------
# CASA
# ---------------------------------------------------------------------------

def read_plane(path):
    """First plane of a CASA image as a 2-D array, with its cell in arcsec."""

    import numpy as np
    from casatools import image

    tool = image()
    tool.open(str(path))
    try:
        data = np.asarray(tool.getchunk(), dtype=float)
        increment = tool.summary()["incr"]
    finally:
        tool.close()
    while data.ndim > 2:
        data = data[..., 0]
    cell = abs(float(increment[0])) * ARCSEC_PER_RADIAN
    return data, cell


def scan_numbers(ms):
    from casatools import msmetadata

    tool = msmetadata()
    tool.open(str(ms))
    try:
        return [int(scan) for scan in tool.scannumbers()]
    finally:
        tool.close()


def image_noise(path):
    plane, _ = read_plane(path)
    return {
        "image": str(path),
        "robust_rms_Jy_beam": robust_rms(plane),
        "rms_Jy_beam": plain_rms(plane),
        "max_abs_Jy_beam": float(abs(plane).max()),
    }


def imaging_arguments(ms, base, imsize, cell, robust, selection):
    return dict(
        vis=str(ms),
        imagename=str(base),
        datacolumn="data",
        specmode="mfs",
        gridder="standard",
        deconvolver="hogbom",
        imsize=int(imsize),
        cell=f"{cell:.4f}arcsec",
        weighting="briggs",
        robust=float(robust),
        stokes="V",
        pblimit=-0.01,
        interactive=False,
        **selection,
    )


def dirty_stokes_v(ms, base, imsize, cell, robust, **selection):
    """Make a dirty Stokes V image and return its noise."""

    from casatasks import tclean

    base = Path(base)
    for old in base.parent.glob(base.name + ".*"):
        shutil.rmtree(old, ignore_errors=True)

    tclean(
        niter=0,
        **imaging_arguments(ms, base, imsize, cell, robust, selection),
    )
    return image_noise(str(base) + ".residual")


def clean_stokes_v(ms, base, imsize, cell, robust, niter, nsigma, nmajor):
    """Clean the leaked sources out of a dirty Stokes V image.

    Continues from the dirty image's own products. It is shallow by design:
    it has to remove the few peaks whose sidelobes fill the image, not
    reach the noise.
    """

    from casatasks import tclean

    tclean(
        niter=int(niter),
        nsigma=float(nsigma),
        nmajor=int(nmajor),
        restart=True,
        calcpsf=False,
        calcres=False,
        **imaging_arguments(ms, base, imsize, cell, robust, {}),
    )
    result = image_noise(str(base) + ".residual")
    result.update({"niter": int(niter), "nsigma": float(nsigma)})
    return result


def model_images(image):
    """The Taylor-term model images that belong to a restored image."""

    image = str(image)
    if ".image." not in image:
        return []
    terms = []
    for term in range(8):
        candidate = Path(
            image.rsplit(".image.", 1)[0] + f".model.tt{term}"
        )
        if not candidate.is_dir():
            break
        terms.append(candidate)
    return terms


def stokes_i_model_residual(ms, base, imsize, cell, robust, model,
                            wprojplanes, grid_cells, **selection):
    """Image what a model leaves in the data, without cleaning.

    tclean with ``niter=0`` and a starting model makes one major cycle: it
    predicts the model, subtracts it from the visibilities and images the
    remainder. With ``savemodel='none'`` nothing is written to the
    Measurement Set. The gridder is the one the model was made with, since
    predicting a wide field without the w-term would leave every off-axis
    source in the residual.
    """

    from casatasks import tclean

    base = Path(base)
    for old in base.parent.glob(base.name + ".*"):
        shutil.rmtree(old, ignore_errors=True)

    tclean(
        vis=str(ms),
        imagename=str(base),
        datacolumn="data",
        specmode="mfs",
        gridder="wproject",
        wprojplanes=int(wprojplanes),
        deconvolver="mtmfs",
        nterms=len(model),
        imsize=int(imsize),
        cell=f"{cell:.4f}arcsec",
        weighting="briggs",
        robust=float(robust),
        stokes="I",
        niter=0,
        startmodel=[str(term) for term in model],
        pblimit=-0.01,
        savemodel="none",
        interactive=False,
        **selection,
    )
    path = str(base) + ".residual.tt0"
    plane, _ = read_plane(path)
    return {
        "image": path,
        "robust_rms_Jy_beam": robust_rms(plane),
        "rms_Jy_beam": plain_rms(plane),
        "max_abs_Jy_beam": float(abs(plane).max()),
        "grid": summarize_grid(rms_grid(plane, grid_cells), plane.shape),
    }


def measure_weighting(ms, workdir, imsize, cell, robust, args):
    """Stokes V floor and point-spread function for one ``robust``."""

    label = f"robust_{robust:g}".replace("-", "m").replace(".", "p")
    base = Path(workdir) / f"stokes_v_{label}"
    entry = {"robust": float(robust)}

    print(f"robust {robust:g}: dirty Stokes V image ...", flush=True)
    entry["stokes_v_dirty"] = dirty_stokes_v(ms, base, imsize, cell, robust)
    print(
        "  dirty: robust rms "
        f"{entry['stokes_v_dirty']['robust_rms_Jy_beam'] * 1e6:.1f} uJy/beam, "
        "largest |V| "
        f"{entry['stokes_v_dirty']['max_abs_Jy_beam'] * 1e3:.2f} mJy/beam",
        flush=True,
    )

    try:
        psf, _ = read_plane(str(base) + ".psf")
        entry["psf"] = psf_statistics(psf, cell)
        del psf
    except Exception as exc:
        entry["psf"] = {"error": f"{type(exc).__name__}: {exc}"}

    if args.v_niter > 0:
        print(f"robust {robust:g}: shallow clean of Stokes V ...", flush=True)
        try:
            entry["stokes_v_cleaned"] = clean_stokes_v(
                ms, base, imsize, cell, robust,
                args.v_niter, args.v_nsigma, args.v_nmajor,
            )
            print(
                "  cleaned: robust rms "
                f"{entry['stokes_v_cleaned']['robust_rms_Jy_beam'] * 1e6:.1f}"
                " uJy/beam",
                flush=True,
            )
        except Exception as exc:
            entry["stokes_v_cleaned"] = {
                "error": f"{type(exc).__name__}: {exc}"
            }
            print(f"  clean failed: {entry['stokes_v_cleaned']['error']}",
                  flush=True)

    if args.model:
        print(
            f"robust {robust:g}: Stokes I residual of the existing model ...",
            flush=True,
        )
        try:
            entry["stokes_i_model_residual"] = stokes_i_model_residual(
                ms, Path(workdir) / f"stokes_i_residual_{label}", imsize,
                cell, robust, args.model, args.wprojplanes, args.grid,
            )
            print(
                "  residual: robust rms "
                f"{entry['stokes_i_model_residual']['robust_rms_Jy_beam'] * 1e6:.1f}"
                " uJy/beam",
                flush=True,
            )
        except Exception as exc:
            entry["stokes_i_model_residual"] = {
                "error": f"{type(exc).__name__}: {exc}"
            }
            print(
                "  failed: "
                f"{entry['stokes_i_model_residual']['error']}",
                flush=True,
            )
    return entry


# ---------------------------------------------------------------------------

def run(args):
    ms, image = resolve_products(args.outdir, args.target)
    if args.ms:
        ms = Path(args.ms)
    if args.image:
        image = Path(args.image)
    if ms is None or not Path(ms).is_dir():
        raise SystemExit(
            f"No self-calibrated Measurement Set found for {args.target} "
            f"under {args.outdir}. Pass --ms."
        )

    workdir = Path(args.workdir) / args.target
    workdir.mkdir(parents=True, exist_ok=True)
    outfile = workdir / "noise_floor.json"
    report = {
        "target": args.target,
        "measurement_set": str(ms),
        "goal_Jy_beam": args.goal * 1e-6,
    }

    def save():
        with open(outfile, "w") as handle:
            json.dump(_plain(report), handle, indent=2, allow_nan=False)

    imsize, cell = args.imsize, args.cell
    if image is not None and Path(image).is_dir():
        plane, image_cell = read_plane(image)
        imsize = imsize or plane.shape[0]
        cell = cell or image_cell
        grid = rms_grid(plane, args.grid)
        report["stokes_i"] = {
            "image": str(image),
            "robust_rms_Jy_beam": robust_rms(plane),
            "brightest": brightest_pixel(plane),
            "grid": summarize_grid(grid, plane.shape),
            "grid_rms_Jy_beam": grid,
        }
        residual = Path(str(image).replace(".image.", ".residual."))
        if residual.is_dir() and residual != Path(image):
            res_plane, _ = read_plane(residual)
            report["stokes_i_residual"] = {
                "image": str(residual),
                "robust_rms_Jy_beam": robust_rms(res_plane),
                "grid": summarize_grid(
                    rms_grid(res_plane, args.grid), res_plane.shape
                ),
            }
        print(
            f"Stokes I image {Path(image).name}: robust rms "
            f"{report['stokes_i']['robust_rms_Jy_beam'] * 1e6:.1f} uJy/beam",
            flush=True,
        )
    if not imsize or not cell:
        raise SystemExit(
            "No Stokes I image to take the grid from. Pass --imsize and "
            "--cell (arcsec)."
        )
    report["imsize"], report["cell_arcsec"] = int(imsize), float(cell)

    args.model = []
    if args.model_residual:
        args.model = model_images(image) if image is not None else []
        if not args.model:
            raise SystemExit(
                "--model-residual needs the model images beside the Stokes "
                f"I image ({image}), and none were found."
            )
        report["model"] = [str(term) for term in args.model]

    # Written after each step: a run that is stopped keeps what it measured.
    report["weightings"] = []
    for robust in args.robust:
        report["weightings"].append(
            measure_weighting(ms, workdir, imsize, cell, robust, args)
        )
        save()

    if not args.skip_scans:
        report["stokes_v_scans"] = {"robust": float(args.robust[0])}
        for scan in scan_numbers(ms):
            print(f"Dirty Stokes V image of scan {scan} ...", flush=True)
            try:
                result = dirty_stokes_v(
                    ms, workdir / f"stokes_v_scan{scan}", imsize, cell,
                    args.robust[0], scan=str(scan),
                )
            except Exception as exc:
                result = {"error": f"{type(exc).__name__}: {exc}"}
            report["stokes_v_scans"][str(scan)] = result
            save()
            if result.get("robust_rms_Jy_beam"):
                print(
                    f"  robust rms "
                    f"{result['robust_rms_Jy_beam'] * 1e6:.1f} uJy/beam",
                    flush=True,
                )
            elif result.get("error"):
                print(f"  failed: {result['error']}", flush=True)

            if not args.model:
                continue
            residuals = report.setdefault(
                "stokes_i_model_residual_scans",
                {"robust": float(args.robust[0])},
            )
            print(f"Stokes I model residual of scan {scan} ...", flush=True)
            try:
                result = stokes_i_model_residual(
                    ms, workdir / f"stokes_i_residual_scan{scan}", imsize,
                    cell, args.robust[0], args.model, args.wprojplanes,
                    args.grid, scan=str(scan),
                )
            except Exception as exc:
                result = {"error": f"{type(exc).__name__}: {exc}"}
            residuals[str(scan)] = result
            save()
            if result.get("robust_rms_Jy_beam"):
                print(
                    f"  robust rms "
                    f"{result['robust_rms_Jy_beam'] * 1e6:.1f} uJy/beam",
                    flush=True,
                )
            elif result.get("error"):
                print(f"  failed: {result['error']}", flush=True)

    report["scans"] = scan_table(
        report.get("stokes_v_scans"),
        report.get("stokes_i_model_residual_scans"),
    )
    report["reading"] = reading(
        report["weightings"],
        (report.get("stokes_i") or {}).get("grid"),
        args.goal * 1e-6,
        report["scans"],
    )
    save()
    print("\n" + "\n".join(report["reading"]))
    print(f"\nReport written: {outfile}")
    return report


def _plain(value):
    """A report as JSON can hold it: no numpy types, no NaN."""

    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    item = getattr(value, "item", None)
    if callable(item):
        value = item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if value is None or isinstance(value, (int, float, str, bool)):
        return value
    return str(value)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--outdir", default="uGMRT_Products")
    parser.add_argument("--target", required=True)
    parser.add_argument("--ms", help="Measurement Set, if not the target's "
                        "final self-calibrated one")
    parser.add_argument("--image", help="Stokes I image to grid, if not the "
                        "target's final self-calibration image")
    parser.add_argument("--workdir", default="noise_floor_test")
    parser.add_argument("--robust", type=float, nargs="+", default=[0.5],
                        help="one or more Briggs robust values to measure "
                        "(default 0.5)")
    parser.add_argument("--imsize", type=int)
    parser.add_argument("--cell", type=float, help="arcsec")
    parser.add_argument("--grid", type=int, default=10,
                        help="boxes per image side (default 10)")
    parser.add_argument("--goal", type=float, default=20.0,
                        help="target noise in uJy/beam (default 20)")
    parser.add_argument("--skip-scans", action="store_true",
                        help="no per-scan images")
    parser.add_argument("--v-niter", type=int, default=5000,
                        help="iterations of the shallow Stokes V clean; 0 "
                        "leaves the image dirty (default 5000)")
    parser.add_argument("--model-residual", action="store_true",
                        help="also image what the existing Stokes I model "
                        "leaves in the data, for each weighting and scan")
    parser.add_argument("--wprojplanes", type=int, default=128,
                        help="w-projection planes for the model residual "
                        "(default 128, the pipeline's)")
    parser.add_argument("--v-nsigma", type=float, default=5.0)
    parser.add_argument("--v-nmajor", type=int, default=4)
    run(parser.parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
