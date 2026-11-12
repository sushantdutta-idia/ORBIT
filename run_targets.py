#!/usr/bin/env python3
"""Take every target from its calibrated data to a primary-beam-corrected image.

One command for all the targets the pipeline calibrated. First, for every
target, **the scans**: ``select_scans.py`` measures each scan from the data
alone and leaves out those that would cost the image more than they add
(finding 130); all the targets' scans are printed before any image is
started. Then for each target, in turn:

1. **the field**: ``selfcal_recipe.py`` on the kept scans, with images of
   the whole field: phases per scan and then within the scans, delays,
   amplitudes, the gains across the band, and an image;
2. **peeling**: ``peel_sources.py`` takes out the bright sources away from
   the pointing centre, found in that image, each with the gains the
   antennas had towards it;
3. **the target**: ``selfcal_recipe.py`` again on what is left, which the
   target now dominates, and the image of the whole field;
4. **the primary beam**: ``correct_primary_beam.py`` on that image;
5. with ``--extra-robust``, **images of other weightings**: the same
   calibrated data imaged again, cleaned as deeply and corrected for the
   primary beam, beside the robust-0 image, never instead of it.

This is the sequence worked out on G0212-0335 (findings 118 to 128). What
was done there by hand and is decided here by rule: the cores-only models
go down to a fraction of each field's own brightest pixel, not to a number
of mJy/beam; the sources to peel are those of at least ``--peel-above``
mJy/beam more than ``--keep-radius`` arcmin from the centre; and how much
is solved towards each follows from how bright it is.

WSClean deconvolves in twelve frequency sub-bands, without which an image
of a band this wide cannot be made (finding 123). Those are kept out of
sight: the product of a target is one image of the whole band.

The targets and their data come from the pipeline's own record,
``<products>/calibration/split_targets.json``; nothing in the products is
written to. Once one target has gone through, a target that fails does
not stop the others; before that a failure stops the run, since what
stopped the first would as likely stop the next. Given again, the same
command goes on where it stopped: each step keeps its own record.

The telescope and band are those the pipeline recorded in the products.
For a MeerKAT field the image is the band's own grid (6144 pixels of 1.5"
at L), the channels those inside the band's guaranteed range, the
reference antennas the one calibration used and then MeerKAT's preferred
ones, and the primary beam the measured katbeam model, applied as the
pipeline's imaging stage applies it (``correct_primary_beam.py
--telescope meerkat``). MeerKAT's data hold XX and YY only, with no
Stokes V, so its images are made without a noise image of the data. Run
it with a Python that has numpy, astropy and katbeam.

Example::

    python3 run_targets.py --channel-range 51 768
    python3 run_targets.py --channel-range 51 768 --scans G0212-0335=15
    ~/venvs/orbit-imaging/bin/python3 run_targets.py \\
        --products MeerKAT_Products_1708779077 --workdir MeerKAT_targets
"""

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# This folder holds the pipeline's own inspect.py; see diagnose_noise_floor.
if sys.path and os.path.abspath(sys.path[0]) == SCRIPT_DIR:
    sys.path.pop(0)

import argparse
import contextlib
import importlib.util
import json
import shlex
import socket
import time
from pathlib import Path


REPORT_NAME = "targets_report.json"
# Each target's own record, written only by the run that holds the target;
# the summary above is made from these.
TARGET_REPORT = "target_report.json"
STAGES = ("scans", "field", "peel", "target", "primary_beam")
# Where the recipe on what is left after peeling works. Named for how it
# is solved, against the whole field only (finding 132): a folder "target"
# left by an earlier run, whose first rounds were against cores-only
# models, is neither used nor removed, and the stage is done again here.
TARGET_FOLDER = "target_whole_field"


def load_imager():
    spec = importlib.util.spec_from_file_location(
        "run_targets_imager", Path(SCRIPT_DIR) / "experiment_wsclean_image.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


IMAGER = load_imager()


def load_recipe():
    spec = importlib.util.spec_from_file_location(
        "run_targets_recipe", Path(SCRIPT_DIR) / "selfcal_recipe.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_profiles():
    spec = importlib.util.spec_from_file_location(
        "run_targets_profiles", Path(SCRIPT_DIR) / "telescope_profiles.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PROFILES = load_profiles()
TELESCOPES = ("ugmrt", "meerkat")


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())


def script(name):
    return [sys.executable, "-u", str(Path(SCRIPT_DIR) / name)]


def read_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def write_json(path, data):
    """Write whole or not at all, even with other runs writing the same."""

    path = Path(path)
    scratch = path.with_name(
        f".{path.name}.{socket.gethostname()}.{os.getpid()}.tmp")
    with open(scratch, "w") as handle:
        json.dump(data, handle, indent=2)
    scratch.replace(path)


def products_context(products):
    """The telescope, band and reference antenna the pipeline recorded."""

    products = Path(products)
    inspection = read_json(products / "inspection" / "inspection_report.json")
    calibration = read_json(
        products / "calibration" / "calibration_report.json")
    inspection = inspection if isinstance(inspection, dict) else {}
    calibration = calibration if isinstance(calibration, dict) else {}

    def name_of(value):
        if isinstance(value, dict):
            value = value.get("name")
        return str(value).strip() if value else None

    refant = calibration.get("reference_antenna")
    if not refant:
        recorded = inspection.get("reference_antenna")
        refant = (recorded.get("recommended")
                  if isinstance(recorded, dict) else recorded)
    telescope = (name_of(inspection.get("telescope"))
                 or name_of(calibration.get("telescope")))
    return {
        "telescope": telescope.lower() if telescope else None,
        "band": (name_of(inspection.get("band"))
                 or name_of(calibration.get("band"))),
        "refant": str(refant).strip() if refant else None,
    }


# The field of each band's images and its pixel. Four or five pixels go
# across the beam the longest baselines give at the band's centre (25 km for
# the GMRT, 7.7 km for MeerKAT), and the field reaches out to where the
# primary beam, at the band's lowest frequency, has fallen to about a tenth:
# the rule that set Band 4's 5120 x 1" (its beam is 0.12 at the edge at 550
# MHz). uGMRT Bands 4 and 5 and MeerKAT L are as their runs were made
# (2026-10 and 11); the others follow the rule and have not been run.
TARGET_GRIDS = {
    ("ugmrt", "B2"): (8192, 3.0),
    ("ugmrt", "B3"): (8192, 1.5),
    ("ugmrt", "B4"): (5120, 1.0),
    ("ugmrt", "B5"): (6144, 0.5),
    ("meerkat", "MEERKAT_UHF"): (7168, 2.0),
    ("meerkat", "MEERKAT_L"): (6144, 1.5),
    **{("meerkat", f"MEERKAT_S{index}"): (8192, 0.6) for index in range(5)},
}


def selection_ranges(selection):
    """``"0:0~52;974~1023,0:0~9"`` as ``{"0": [(0, 52), (974, 1023), (0, 9)]}``.

    None if any part cannot be read.
    """

    found = {}
    for part in str(selection or "").split(","):
        part = part.strip()
        if not part:
            continue
        spw, _, channels = part.partition(":")
        if not channels:
            return None
        for piece in channels.split(";"):
            low, _, high = piece.partition("~")
            try:
                first, last = int(low), int(high or low)
            except ValueError:
                return None
            found.setdefault(spw.strip(), []).append((first, last))
    return found


def recorded_channel_range(products):
    """The channels of a calibrated target that hold data, END exclusive.

    What the flagging stage masked at the band's edges, in the target's own
    channels. Those masks are a share of each end of the band, and for
    MeerKAT and uGMRT Band 4 also what lies outside the guaranteed range.
    The calibration stage may have averaged several native channels into
    each of the target's (``calibration/target_channel_average.json``), and
    only channels made wholly of data are kept. Band 5 gives 122-3974, the
    range its run was given by hand. None for data in more than one
    spectral window, or without the flagging stage's record.
    """

    products = Path(products)
    inspection = read_json(products / "inspection" / "inspection_report.json")
    spws = inspection.get("spws") if isinstance(inspection, dict) else None
    guaranteed = read_json(
        products / "flagging" / "guaranteed_range_mask_report.json")
    guaranteed = guaranteed if isinstance(guaranteed, dict) else None
    total = spw = None
    if isinstance(spws, dict) and len(spws) == 1:
        (spw, record), = spws.items()
        try:
            total = int(record["nchan"])
        except (KeyError, TypeError, ValueError):
            total = None
    if total is None and guaranteed:
        items = guaranteed.get("spws") or []
        if len(items) == 1:
            spw = items[0].get("spw")
            try:
                total = int(items[0]["total_channels"])
            except (KeyError, TypeError, ValueError):
                total = None
    if total is None or spw is None:
        return None

    masked, recorded = set(), False
    edges = read_json(products / "flagging" / "edge_flagging_report.json")
    if isinstance(edges, dict):
        ranges = selection_ranges(
            edges.get("spw_string") or edges.get("spw_selection") or "")
        if ranges is None or set(ranges) - {str(spw)}:
            return None
        for first, last in ranges.get(str(spw), []):
            masked.update(range(first, last + 1))
        recorded = True
    if guaranteed:
        for item in guaranteed.get("spws") or []:
            if str(item.get("spw")) != str(spw):
                return None
            try:
                for first, last in item.get("channel_ranges") or []:
                    masked.update(range(int(first), int(last) + 1))
            except (TypeError, ValueError):
                return None
        recorded = True
    if not recorded:
        return None
    kept = [channel for channel in range(total) if channel not in masked]
    if not kept:
        return None
    first, end = kept[0], kept[-1] + 1
    width = 1
    average = read_json(
        products / "calibration" / "target_channel_average.json")
    if isinstance(average, dict):
        try:
            width = max(1, int(average.get("width", 1)))
        except (TypeError, ValueError):
            return None
    first, end = -(-first // width), end // width
    return (first, end) if end > first else None


def has_primary_beam(telescope, profile):
    """Whether a band has a measured primary-beam model in the pipeline."""

    if telescope == "ugmrt":
        return bool(profile.get("primary_beam_coefficients"))
    return bool(profile.get("katbeam_model"))


def resolve_telescope(args, error):
    """Fill in what depends on the telescope and band, from the products.

    ``error`` is called with a message when the products and the options
    disagree, or a setting cannot be found. Products that record nothing
    are taken as uGMRT Band 4, as this chain always did.
    """

    context = products_context(args.products)
    recorded = context["telescope"]
    if args.telescope and recorded and args.telescope != recorded:
        error(f"--telescope {args.telescope}, but the products "
              f"{args.products} are of {recorded}")
    telescope = args.telescope or recorded or "ugmrt"
    if telescope not in TELESCOPES:
        error(f"the products {args.products} are of {telescope}; this "
              "chain runs uGMRT and MeerKAT fields")
    args.telescope = telescope

    band = str(
        args.band or context["band"] or ("B4" if telescope == "ugmrt" else "")
    ).upper()
    if args.band and context["band"] and band != context["band"].upper():
        error(f"--band {args.band}, but the products {args.products} are of "
              f"{context['band']}")
    profile = PROFILES.TELESCOPE_PROFILES[telescope]["bands"].get(band)
    if profile is None:
        error(f"no {telescope} band '{band or 'unknown'}' in "
              "telescope_profiles.py; give --band if the products do not say")
    args.band = band
    grid = TARGET_GRIDS.get((telescope, band))
    if grid is None and not (args.imsize and args.cell):
        error(f"no field geometry for {telescope} {band}: give --imsize and "
              "--cell")
    args.imsize = args.imsize or grid[0]
    args.cell = args.cell or grid[1]
    if args.channel_range is None:
        found = recorded_channel_range(args.products)
        if found is None:
            error("--channel-range: the products hold no record of the "
                  "channels their flagging left at the band's edges")
        args.channel_range = list(found)
    if not has_primary_beam(telescope, profile) and not args.no_primary_beam:
        error(f"{telescope} {band} has no measured primary-beam model here "
              "(katbeam publishes MeerKAT UHF and L only): pass "
              "--no-primary-beam to image it without the correction")
    if telescope == "ugmrt":
        # The recipe's own reference antennas are uGMRT's, as every uGMRT
        # run of this chain has used them.
        return args
    if not args.refant:
        preferred = (
            PROFILES.TELESCOPE_PROFILES["meerkat"].get("preferred_refants")
            or [])
        order = [context["refant"]] + list(preferred)
        args.refant = ",".join(dict.fromkeys(name for name in order if name))
    return args


def targets_of(products):
    """The targets the pipeline calibrated, and where their data are."""

    record = Path(products) / "calibration" / "split_targets.json"
    found = read_json(record)
    if not isinstance(found, dict) or not found:
        raise SystemExit(
            f"No calibrated targets recorded in {record}. --products is the "
            "pipeline's output folder."
        )
    return {str(name): Path(path) for name, path in found.items()}


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def averaged_channel_range(args):
    """The channel range in the data the images are made from.

    ``--channel-range`` counts the calibrated data's channels. With
    ``--average-channels N`` the images are made from a copy with every N
    channels averaged; its range keeps only channels made wholly of data.
    """

    first, end = (int(value) for value in args.channel_range)
    average = int(getattr(args, "average_channels", 1) or 1)
    if average == 1:
        return first, end
    first, end = -(-first // average), end // average
    if end <= first:
        raise SystemExit(
            f"--channel-range {args.channel_range[0]} {args.channel_range[1]}"
            f" leaves no whole channel after averaging every {average}."
        )
    return first, end


def shared(args):
    first, end = averaged_channel_range(args)
    options = [
        "--outdir", str(args.products),
        "--channel-range", str(first), str(end),
        "--wsclean", args.wsclean,
    ]
    if args.casa:
        options += ["--casa", args.casa]
    return options


def forced_scans(values, available):
    """``--scans TARGET=3,15 ...`` as the scans to use for each target."""

    forced = {}
    for value in values or []:
        target, _, scans = value.partition("=")
        try:
            numbers = [int(scan) for scan in scans.split(",")]
        except ValueError:
            numbers = []
        if not target or not numbers:
            raise SystemExit(
                f"--scans {value}: give a target and its scans, as in "
                "G0212-0335=15 or G0216-0604=5,18."
            )
        if target not in available:
            raise SystemExit(
                f"--scans {value}: {target} is not among the calibrated "
                f"targets: {', '.join(available)}."
            )
        forced[target] = numbers
    return forced


def scans_rule(args, target):
    """What decides a target's scans; a change of it measures them again."""

    rule = {
        "least_share": float(args.least_share),
        "most_excess": float(args.most_excess),
        "scans": args.forced.get(target),
    }
    average = int(getattr(args, "average_channels", 1) or 1)
    if average > 1:
        rule["average_channels"] = average
    return rule


def scans_command(args, target, ms, folder):
    """Which scans of a target are worth imaging, from the data alone."""

    command = script("select_scans.py") + [
        "--target", target, "--ms", str(ms),
        "--workdir", str(folder / "scans"),
        "--least-share", f"{float(args.least_share):g}",
        "--most-excess", f"{float(args.most_excess):g}",
        "--outdir", str(args.products),
        "--channel-range",
        str(int(args.channel_range[0])), str(int(args.channel_range[1])),
    ]
    if args.forced.get(target):
        command += ["--scans", *(str(scan) for scan in args.forced[target])]
    average = int(getattr(args, "average_channels", 1) or 1)
    if average > 1:
        command += ["--average-channels", str(average)]
    if args.casa:
        command += ["--casa", args.casa]
    return command


def resources(args):
    """Threads, memory and parallel cleaning, as each script is told them."""

    options = [
        "--threads", str(int(args.threads)),
        "--memory-gb", f"{float(args.memory_gb):g}",
    ]
    if args.parallel_deconvolution is not None:
        options += [
            "--parallel-deconvolution", str(int(args.parallel_deconvolution)),
        ]
    return options


def refant_options(args):
    """The reference antennas, where the recipe's own (uGMRT's) will not do."""

    refant = getattr(args, "refant", None)
    return ["--refant", refant] if refant else []


def recipe_options(args):
    return [
        "--cell", f"{float(args.cell):g}",
        "--channels-out", str(int(args.channels_out)),
        "--core-fraction", f"{float(args.core_fraction):g}",
        *refant_options(args),
        *resources(args),
        "--tidy",
    ]


def time_phase_rounds(args):
    """The interval of the phases within a scan, or None for no such rounds.

    A faint field cannot carry them: G0231-0533 in Band 5 (one 12-minute
    scan, brightest source 3-4 mJy/beam) lost a third of its peak in its
    1-minute rounds, which every brighter field gained from.
    """

    value = str(args.time_phase_solint or "").strip()
    return None if value.lower() in ("", "none", "off", "0") else value


def field_command(args, target, ms, folder):
    """The whole field, from the kept scans, not yet self-calibrated."""

    interval = time_phase_rounds(args)
    return script("selfcal_recipe.py") + [
        "--target", target, "--ms", str(ms),
        "--workdir", str(folder / "field"),
        "--imsize", str(int(args.imsize)),
        *(["--time-phase-solint", interval] if interval else []),
        # The pipeline's calibrated data have had no amplitude
        # self-calibration; what is left of the field after peeling has.
        *(["--phase-only"] if args.phase_only else ["--amplitude-cores"]),
        *recipe_options(args), *shared(args),
    ]


def peel_command(args, target, ms, models, image, folder):
    return script("peel_sources.py") + [
        "--target", target, "--ms", str(ms), "--models", str(models),
        "--auto-sources", str(image),
        "--peel-above", f"{float(args.peel_above):g}",
        "--keep-radius", f"{float(args.keep_radius):g}",
        "--most", str(int(args.most)),
        "--workdir", str(folder / "peel"),
        "--field-imsize", str(int(args.imsize)),
        "--cell", f"{float(args.cell):g}",
        "--channels-out", str(int(args.channels_out)),
        *refant_options(args),
        *resources(args),
        "--no-image", "--tidy",
        *shared(args),
    ]


def target_command(args, target, ms, folder):
    """What is left once the bright off-axis sources are out.

    The models need only the inner field, where the sources that now
    matter are; the image is of the whole field. Nothing a solve cannot
    solve is flagged: the target alone is fainter than the field was.

    No round is against a cores-only model. These data are calibrated,
    and with the bright sources out such a model holds a tenth of the
    field's flux: on G0212-0335 the phases solved against it moved the
    central square by 8 to 14 degrees and doubled the holes beside the
    sources near the axis (finding 132).
    """

    return script("selfcal_recipe.py") + [
        "--target", target, "--ms", str(ms),
        "--workdir", str(folder / TARGET_FOLDER),
        "--imsize", str(int(args.inner_imsize)),
        "--final-imsize", str(int(args.imsize)),
        "--whole-field-only",
        "--bandpass-channels", "60", "--keep-unsolved",
        *(["--phase-only"]
          if args.phase_only or getattr(args, "phase_only_after_peel", False)
          else []),
        *recipe_options(args), *shared(args),
    ]


def robust_label(robust, taper=0):
    """As the images are named: robust_0, robust_0p5, robust_0_taper15."""

    return IMAGER.image_label(robust, taper)


# Tapered images are made at robust 0, on a grid of about five pixels to
# the beam over the same field, so that the cleaning scales (in pixels)
# grow with the beam.
TAPER_ROBUST = 0.0
TAPER_PIXELS_PER_BEAM = 5.0
# Untapered, robust 2 put nearly all the weight on the central square: on
# all nine Band 5 targets (2026-10-26) the clean diverged, to 0.7 mJy/beam
# and beyond, with fitted beams such as 14.4" x 1.9" (finding 157).
MOST_UNTAPERED_ROBUST = 1.0


def taper_grid(args, taper):
    """Cell (arcsec) and image size for a tapered image.

    The field is the main image's, times ``--taper-field-scale``: a taper
    leans on the short baselines, which see the sources beyond the main
    field, in the primary beam's first sidelobe, at their brightest. Left
    outside the image they are not cleaned, and at 15" their sidelobes
    left Stokes I 1.2-2 times as noisy as Stokes V in the quietest parts
    of all nine Band 5 fields (2026-10-28).
    """

    cell = max(float(args.cell), round(float(taper) / TAPER_PIXELS_PER_BEAM, 3))
    scale = float(getattr(args, "taper_field_scale", 1.0) or 1.0)
    imsize = int(round(
        float(args.imsize) * float(args.cell) * scale / cell / 2.0)) * 2
    return cell, max(imsize, 64)


def extra_folder_label(args, robust, taper):
    """The folder of an extra image: its weighting, and a wider field."""

    label = robust_label(robust, taper)
    scale = float(getattr(args, "taper_field_scale", 1.0) or 1.0)
    if taper and scale != 1.0:
        label += f"_field{scale:g}".replace(".", "p")
    return label


# The deconvolution of the recipe's final image (selfcal_recipe.py's own
# defaults), so that an image of another weighting differs in that alone.
FINAL_SCALES = (0, 5, 15, 30)
FINAL_AUTO_THRESHOLD = 1.0
FINAL_NMITER = 20


def extra_image_command(args, target, ms, workdir, robust, taper=0):
    """The calibrated data imaged again with another weighting.

    Robust 0 gives the finest beam the arms allow and pays for it in
    noise: on the GMRT, whose central square holds most of the
    visibilities, it weights them down. In Band 5 the robust-0 images of
    two 12-minute scans are at 1.5 times the noise natural weights allow
    (G0222-0321: 21.2 against 13.7 uJy/beam). A weighting nearer natural
    takes part of that back for a broader beam.

    With ``taper`` (arcsec) the long baselines are weighted down to that
    beam, on a coarser grid of the same field: the targets' own extended
    emission, 4-10 mJy over 65-136" in NVSS (G0216-0604, G0220-0318,
    G0222-0321, G0231-0533), is a few µJy in each 2" beam and shows only
    so.
    """

    first, end = averaged_channel_range(args)
    cell, imsize = (
        taper_grid(args, taper) if taper else (float(args.cell), args.imsize)
    )
    return script("experiment_wsclean_image.py") + [
        "--ms", str(ms), "--workdir", str(workdir), "--target", target,
        "--outdir", str(args.products),
        "--imsize", str(int(imsize)),
        "--cell", f"{float(cell):g}",
        "--robust", f"{float(robust):g}",
        *(["--taper", f"{float(taper):g}"] if taper else []),
        "--channels-out", str(int(args.channels_out)),
        "--channel-range", str(first), str(end),
        "--wsclean", args.wsclean,
        *resources(args),
        "--scales", *[str(scale) for scale in FINAL_SCALES],
        "--no-local-rms",
        "--auto-threshold", f"{FINAL_AUTO_THRESHOLD:g}",
        "--nmiter", str(FINAL_NMITER),
        "--overwrite",
    ]


def is_meerkat(args):
    return getattr(args, "telescope", "ugmrt") == "meerkat"


def beam_command(args, name):
    if is_meerkat(args):
        # The katbeam correction, with the masks the flagging stage applied.
        return script("correct_primary_beam.py") + [
            "--name", str(name), "--telescope", "meerkat",
            "--band", args.band, "--products", str(args.products),
            "--pbmin", f"{float(args.pbmin):g}", "--overwrite",
        ]
    return script("correct_primary_beam.py") + [
        "--name", str(name), "--band", args.band,
        "--pbmin", f"{float(args.pbmin):g}", "--overwrite",
    ]


def corrected_image(args, name):
    """Where the beam step leaves the corrected image of ``name``.

    The uGMRT step names it ``-MFS-image.pbcor.fits``; the katbeam step
    follows WSClean's ``-MFS-image-pb.fits``, as the pipeline's does.
    """

    suffix = "-MFS-image-pb.fits" if is_meerkat(args) else "-MFS-image.pbcor.fits"
    return Path(str(name) + suffix)


# ---------------------------------------------------------------------------
# What a step left
# ---------------------------------------------------------------------------

def scans_result(folder, args, target, source):
    """The scans kept for a target, if decided for these data by this rule."""

    report = read_json(folder / "scans" / "scans_report.json") or {}
    if (
        report.get("input") != str(source)
        or report.get("rule") != scans_rule(args, target)
        or not report.get("measurement_set")
        or not Path(report["measurement_set"]).is_dir()
    ):
        return None
    return {
        "measurement_set": report["measurement_set"],
        "kept": report.get("kept"),
        "left_out": report.get("left_out"),
        "reading": report.get("reading") or [],
    }


def recipe_input(folder, target):
    """The data a recipe's folder was started on, if it was started."""

    return (read_json(folder / target / "recipe_report.json") or {}).get(
        "input")


def recipe_result(folder, target):
    """Calibrated data, the name of the final image, and its numbers."""

    report = read_json(folder / target / "recipe_report.json") or {}
    if "final:image" not in report.get("done", []):
        return None
    first = sorted(
        (folder / target / "images" / "final").glob("*/*_I-MFS-image.fits")
    )
    if not first:
        return None
    name = first[0].with_name(first[0].name[:-len("-MFS-image.fits")])
    final = report.get("final") or {}
    return {
        "measurement_set": report.get("calibrated_measurement_set"),
        "name": str(name),
        "image": str(first[0]),
        "image_rms_Jy_beam": final.get("image_rms_Jy_beam"),
        "quietest_tenth_Jy_beam": final.get("quietest_tenth_Jy_beam"),
        "data_noise_Jy_beam": final.get("data_noise_Jy_beam"),
        "peak_Jy_beam": final.get("peak_Jy_beam"),
        "most_negative_Jy_beam": final.get("most_negative_Jy_beam"),
        "rounds": report.get("reading") or [],
    }


def restoring_beam(path):
    """BMAJ and BMIN of a FITS image, in arcsec, or None."""

    found = {}
    try:
        with open(path, "rb") as handle:
            while True:
                block = handle.read(2880)
                if len(block) < 2880:
                    return None
                for start in range(0, 2880, 80):
                    card = block[start:start + 80].decode("ascii", "replace")
                    key = card[:8].strip()
                    if key == "END":
                        if "BMAJ" in found and "BMIN" in found:
                            return [round(found["BMAJ"] * 3600, 2),
                                    round(found["BMIN"] * 3600, 2)]
                        return None
                    if key in ("BMAJ", "BMIN") and card[8:10] == "= ":
                        found[key] = float(card[10:].split("/")[0])
    except (OSError, ValueError):
        return None


def extra_result(name, robust, corrected, taper=0):
    """An image of another weighting: where it is, its noise, its beam."""

    name = Path(name)
    label = robust_label(robust, taper)
    report = read_json(name.parent / f"wsclean_{label}.json") or {}
    image = report.get("stokes_i") or {}
    return {
        "robust": float(robust),
        "taper_arcsec": float(taper) if taper else None,
        "image": str(name) + "-MFS-image.fits",
        "corrected": str(corrected) if corrected else None,
        "image_rms_Jy_beam": image.get("robust_rms_Jy_beam"),
        "quietest_tenth_Jy_beam": (image.get("grid") or {}).get(
            "quietest_tenth_rms"),
        "data_noise_Jy_beam": IMAGER.data_noise_of(report).get(
            "robust_rms_Jy_beam"),
        "beam_arcsec": restoring_beam(str(name) + "-MFS-image.fits"),
    }


def extra_images(args):
    """(robust, taper) of every extra image asked for, in order."""

    return (
        [(float(robust), 0.0) for robust in args.extra_robust or []] +
        [(TAPER_ROBUST, float(taper)) for taper in args.extra_taper or []]
    )


def micro(value):
    return "?" if value is None else f"{value * 1e6:.1f}"


def summary_line(target, entry):
    """A target in a line: what it came to, or where it stopped."""

    if entry.get("failed"):
        return f"{target}: stopped at '{entry['failed']['stage']}'."
    result = entry.get("result")
    if not result:
        return f"{target}: not finished."
    left_out = (entry.get("scans") or {}).get("left_out") or []
    noise = result.get("data_noise_Jy_beam")
    return (
        f"{target}: {micro(result.get('image_rms_Jy_beam'))} uJy/beam "
        f"({micro(result.get('quietest_tenth_Jy_beam'))} in the quietest "
        "tenth)"
        + (f" against {micro(noise)} for the data" if noise is not None
           else "")
        + f"; {entry.get('sources_out', 0)} sources taken out; "
        + (
            "without scan" + ("s " if len(left_out) > 1 else " ")
            + ", ".join(str(scan) for scan in left_out) + "; "
            if left_out else ""
        )
        + (
            "not corrected for the primary beam"
            if "corrected" in entry and entry["corrected"] is None
            else f"{entry.get('corrected', '?')}"
        )
        + "".join(
            (f"; taper {extra['taper_arcsec']:g}\": "
             if extra.get("taper_arcsec") else
             f"; robust {extra['robust']:g}: ")
            + f"{micro(extra.get('image_rms_Jy_beam'))} uJy/beam "
            f"({micro(extra.get('quietest_tenth_Jy_beam'))} in the quietest "
            "tenth)"
            + (
                f", beam {extra['beam_arcsec'][0]:.2f}\" x "
                f"{extra['beam_arcsec'][1]:.2f}\""
                if extra.get("beam_arcsec") else ""
            )
            for extra in (entry.get("extra_images") or {}).values()
        )
    )


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

def run(args):
    available = targets_of(args.products)
    args.forced = forced_scans(args.scans, available)
    wanted = args.targets or list(available)
    unknown = [name for name in wanted if name not in available]
    if unknown:
        raise SystemExit(
            f"Not among the calibrated targets: {', '.join(unknown)}. "
            f"Available: {', '.join(available)}."
        )
    workdir = Path(args.workdir)
    if IMAGER.inside(workdir, args.products):
        raise SystemExit(
            f"{workdir} is inside the pipeline's products ({args.products}). "
            "Give --workdir a folder of its own."
        )

    if args.dry_run:
        for target in wanted:
            folder = workdir / target
            print(f"# {target}")
            print(shlex.join(scans_command(
                args, target, available[target], folder)))
            print(shlex.join(field_command(
                args, target, "<kept scans>", folder)))
            print(shlex.join(peel_command(
                args, target, "<calibrated>", "<final image name>",
                "<final image>", folder)))
            print(shlex.join(target_command(args, target, "<peeled>", folder)))
            if not getattr(args, "no_primary_beam", False):
                print(shlex.join(beam_command(args, "<final image name>")))
            for robust, taper in extra_images(args):
                print(shlex.join(extra_image_command(
                    args, target, "<calibrated>", "<final images>", robust,
                    taper)))
        return None

    workdir.mkdir(parents=True, exist_ok=True)
    # A run of an earlier version held the whole folder; one that is still
    # going would not see the claims below.
    earlier = IMAGER.held_by(workdir / IMAGER.LOCK_NAME)
    if earlier is not None:
        raise SystemExit(
            f"{workdir} is held as a whole by another run (process "
            f"{earlier.get('pid', '?')} on {earlier.get('host', '?')}, since "
            f"{earlier.get('since', '?')}). Wait for it to end, or stop it; "
            f"if it is gone, remove {workdir / IMAGER.LOCK_NAME}."
        )
    report_path = workdir / REPORT_NAME
    report = read_json(report_path) or {"targets": {}}
    report["products"] = str(args.products)
    if not isinstance(report.get("targets"), dict):
        report["targets"] = {}
    for name in available:
        recorded = read_json(workdir / name / TARGET_REPORT)
        if recorded is not None:
            report["targets"][name] = recorded

    with contextlib.ExitStack() as held:
        # Each target is held by the run that works on it, for the whole
        # run: several runs, one per target as SLURM jobs, can share the
        # folder, and none takes a target another is on.
        mine = []
        for target in wanted:
            (workdir / target / "logs").mkdir(parents=True, exist_ok=True)
            try:
                held.enter_context(IMAGER.claimed(workdir / target))
            except SystemExit as busy:
                print(f"\n{target}: left to the run that has it. {busy}",
                      flush=True)
                continue
            mine.append(target)

        def save():
            for name in mine:
                if name in report["targets"]:
                    write_json(workdir / name / TARGET_REPORT,
                               report["targets"][name])
            # What the other runs have recorded of theirs, as it now stands.
            for name in available:
                if name in mine:
                    continue
                recorded = read_json(workdir / name / TARGET_REPORT)
                if recorded is not None:
                    report["targets"][name] = recorded
            report["running"] = [
                {"target": name, **entry["running"]}
                for name, entry in report["targets"].items()
                if entry.get("running")
            ]
            report["updated"] = now()
            report["reading"] = [
                summary_line(name, report["targets"][name])
                for name in available if name in report["targets"]
            ]
            write_json(report_path, report)

        def step(target, stage, command):
            entry = report["targets"][target]
            print(f"\n##### {target}: {stage} ({now()}) #####", flush=True)
            entry["running"] = {"stage": stage, "since": now()}
            save()
            IMAGER.run_command(
                command, workdir / target / "logs" / f"{stage}.log")
            entry["stages"][stage] = now()
            entry["running"] = None
            save()

        def stopped(target, stage, failure, own=False):
            """Record a target's failure; say whether the run goes on.

            ``own`` is a failure that says nothing of the other targets.
            """

            entry = report["targets"][target]
            entry["failed"] = {
                "stage": stage, "error": str(failure), "at": now(),
            }
            entry["running"] = None
            print(
                f"\n##### {target}: STOPPED at '{stage}': {failure} #####",
                flush=True,
            )
            save()
            if own or any(report["targets"].get(name, {}).get("finished")
                          for name in available):
                return True
            # Nothing has gone through yet: what stopped this target
            # would as likely stop the next, hours in. Once one target
            # has finished, a failure is that target's own.
            print(
                "\nNo target has finished yet, so the others are not "
                "started: the fault is as likely the chain's as this "
                "target's.",
                flush=True,
            )
            return False

        # Every target's scans first: minutes each, and what they show is
        # worth knowing before the first image, hours long, is started.
        going = True
        for target in mine:
            folder = workdir / target
            entry = report["targets"].setdefault(target, {"stages": {}})
            entry.setdefault("stages", {})
            entry.pop("failed", None)
            entry["running"] = None
            entry["input"] = str(available[target])
            if not available[target].is_dir():
                stopped(
                    target, "scans",
                    f"its calibrated data {available[target]} are missing",
                    own=True,
                )
                continue
            try:
                scans = scans_result(folder, args, target, available[target])
                if scans is None:
                    step(target, "scans", scans_command(
                        args, target, available[target], folder))
                    scans = scans_result(
                        folder, args, target, available[target])
                if scans is None:
                    raise RuntimeError("no scans were chosen")
                entry["scans"] = scans
            except (IMAGER.CommandFailed, RuntimeError) as failure:
                going = stopped(target, "scans", failure)
                if not going:
                    break
        save()
        print("\n##### The scans #####", flush=True)
        for target in mine:
            chosen = report["targets"].get(target, {}).get("scans")
            if chosen:
                print("\n".join(
                    f"{target}: {line}" if line.startswith("scan ") else line
                    for line in chosen["reading"]), flush=True)

        for target in mine if going else []:
            folder = workdir / target
            entry = report["targets"][target]
            if entry.get("failed"):
                continue

            try:
                stage = "field"
                source = entry["scans"]["measurement_set"]
                started_on = recipe_input(folder / "field", target)
                if started_on is not None and started_on != str(source):
                    raise RuntimeError(
                        f"{folder / 'field'} holds a run on {started_on}, "
                        f"not on {source}: the scans to use have changed. "
                        "Give --workdir another folder."
                    )
                field = recipe_result(folder / "field", target)
                if field is None:
                    step(target, stage, field_command(
                        args, target, source, folder))
                    field = recipe_result(folder / "field", target)
                if field is None:
                    raise RuntimeError("the field's recipe left no image")

                stage = "peel"
                peel = read_json(folder / "peel" / "peel_report.json") or {}
                if "chosen" not in peel:
                    step(target, stage, peel_command(
                        args, target, field["measurement_set"], field["name"],
                        field["image"], folder))
                    peel = read_json(folder / "peel" / "peel_report.json") or {}
                chosen = peel.get("chosen")
                if chosen is None:
                    raise RuntimeError("peeling left no record")
                entry["sources_out"] = len(chosen)
                entry["sources"] = chosen

                stage = "target"
                result = field
                if chosen:
                    result = recipe_result(folder / TARGET_FOLDER, target)
                    if result is None:
                        step(target, stage, target_command(
                            args, target, peel["measurement_set"], folder))
                        result = recipe_result(folder / TARGET_FOLDER, target)
                    if result is None:
                        raise RuntimeError("the target's recipe left no image")
                entry["result"] = result

                stage = "primary_beam"
                if getattr(args, "no_primary_beam", False):
                    # A band without a measured beam: the image as the
                    # telescope saw the sky, and the report says so.
                    entry["corrected"] = None
                else:
                    corrected = corrected_image(args, result["name"])
                    if not corrected.is_file():
                        step(target, stage, beam_command(args, result["name"]))
                    if not corrected.is_file():
                        raise RuntimeError("no corrected image was written")
                    entry["corrected"] = str(corrected)

                # Images of other weightings, beside the robust-0 one. The
                # report lists those this run was asked for.
                entry["extra_images"] = {}
                for robust, taper in extra_images(args):
                    label = robust_label(robust, taper)
                    folder_label = extra_folder_label(args, robust, taper)
                    stage = folder_label
                    images = (
                        Path(result["image"]).parents[2] / f"final_{folder_label}"
                    )
                    name = images / target / f"wsclean_{label}_I"
                    image = Path(str(name) + "-MFS-image.fits")
                    if not image.is_file():
                        step(target, stage, extra_image_command(
                            args, target, result["measurement_set"], images,
                            robust, taper))
                        load_recipe().tidy_images(images / target, final=True)
                    if not image.is_file():
                        raise RuntimeError(
                            f"no image of {folder_label} was written")
                    stage = f"{folder_label}_primary_beam"
                    extra = None
                    if not getattr(args, "no_primary_beam", False):
                        extra = corrected_image(args, name)
                        if not extra.is_file():
                            step(target, stage, beam_command(args, name))
                        if not extra.is_file():
                            raise RuntimeError(
                                f"no corrected image of {folder_label} was "
                                "written")
                    entry["extra_images"][folder_label] = (
                        extra_result(name, robust, extra, taper))
                    if taper:
                        cell, imsize = taper_grid(args, taper)
                        entry["extra_images"][folder_label].update({
                            "cell_arcsec": cell, "imsize": imsize,
                            "field_arcmin": round(cell * imsize / 60.0, 1),
                        })
                entry["finished"] = now()
            except (IMAGER.CommandFailed, RuntimeError) as failure:
                # This target stops here; the others do not wait for it.
                if not stopped(target, stage, failure):
                    break
            entry["running"] = None
            save()

        save()
    print("\n" + "\n".join(report["reading"]))
    print(f"\nReport: {report_path}")
    failed = [
        name for name in mine
        if report["targets"].get(name, {}).get("failed")
    ]
    if failed:
        raise SystemExit(
            f"{len(failed)} of {len(mine)} targets stopped: "
            f"{', '.join(failed)}. The same command goes on with them."
        )
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--products", default="uGMRT_Products",
                        help="the pipeline's output folder; read, never "
                        "written to")
    parser.add_argument("--workdir", default="targets_run")
    parser.add_argument("--targets", nargs="+",
                        help="only these targets, in this order (default: "
                        "all the pipeline calibrated)")
    parser.add_argument("--average-channels", type=int, default=1,
                        metavar="N",
                        help="average every N channels of the calibrated "
                        "data before imaging (default 1: none); "
                        "--channel-range still counts the calibrated data's "
                        "channels")
    parser.add_argument("--channel-range", type=int, nargs=2,
                        metavar=("FIRST", "END"),
                        help="channels that hold data, END exclusive "
                        "(needed for uGMRT; for MeerKAT the band's "
                        "guaranteed channels by default)")
    parser.add_argument("--telescope", choices=TELESCOPES,
                        help="default: the telescope the products record")
    parser.add_argument("--scans", nargs="+", metavar="TARGET=SCANS",
                        help="the scans to use for a target, whatever the "
                        "rules say, as in G0212-0335=15 or G0216-0604=5,18")
    parser.add_argument("--least-share", type=float, default=0.10,
                        help="a scan with less than this part of its "
                        "target's weight is left out")
    parser.add_argument("--most-excess", type=float, default=1.4,
                        help="a scan more than this many times as noisy as "
                        "the target's main scan, for the same amount of "
                        "data, is left out")
    parser.add_argument("--imsize", type=int,
                        help="pixels of the image of the whole field "
                        "(default 5120 for uGMRT, the band's own for "
                        "MeerKAT: 6144 at L)")
    parser.add_argument("--inner-imsize", type=int, default=3072,
                        help="pixels of the models once the bright off-axis "
                        "sources are out")
    parser.add_argument("--cell", type=float,
                        help="arcsec (default 1 for uGMRT, the band's own "
                        "for MeerKAT: 1.5 at L)")
    parser.add_argument("--channels-out", type=int, default=12,
                        help="frequency sub-bands WSClean deconvolves in")
    parser.add_argument("--core-fraction", type=float, default=0.13,
                        help="the cores-only models go down to this "
                        "fraction of the field's brightest pixel")
    parser.add_argument("--time-phase-solint", default="1min",
                        help="interval of the phases within a scan; 'none' "
                        "leaves those rounds out, for a field too faint "
                        "for them")
    parser.add_argument("--peel-above", type=float, default=10.0,
                        metavar="MJY")
    parser.add_argument("--keep-radius", type=float, default=4.0,
                        metavar="ARCMIN")
    parser.add_argument("--most", type=int, default=8,
                        help="at most this many sources peeled per target")
    parser.add_argument("--phase-only", action="store_true",
                        help="solve no amplitudes in the field and target "
                        "stages, only phases (and delays): for fields too "
                        "faint for amplitude self-calibration")
    parser.add_argument("--phase-only-after-peel", action="store_true",
                        help="solve no amplitudes once the bright sources "
                        "are peeled: the field stage as usual, the target "
                        "stage in phase only; for a field whose bright "
                        "sources carry its amplitudes and are then taken "
                        "out")
    parser.add_argument("--extra-robust", type=float, nargs="+",
                        metavar="R",
                        help="also image each target's calibrated data with "
                        "these Briggs robust values, beside the robust-0 "
                        "image: lower noise for a broader beam (0.5 is "
                        "about 35 minutes more per target in Band 5)")
    parser.add_argument("--extra-taper", type=float, nargs="+",
                        metavar="ARCSEC",
                        help="also image each target's calibrated data "
                        "tapered to beams of these sizes, at robust 0 on a "
                        "grid of five pixels to the beam: for extended "
                        "emission of low surface brightness")
    parser.add_argument("--taper-field-scale", type=float, default=1.0,
                        metavar="FACTOR",
                        help="make the tapered images this many times as "
                        "wide as the main image, so that the sources in the "
                        "primary beam's first sidelobe are cleaned too "
                        "(default 1); they go to folders of their own, as "
                        "final_robust_0_taper15_field3")
    parser.add_argument("--band",
                        help="the band of the primary beam (default B4 for "
                        "uGMRT; for MeerKAT the band the products record)")
    parser.add_argument("--refant",
                        help="reference antennas of every solve, in order "
                        "of preference (default for MeerKAT: the one the "
                        "calibration used, then MeerKAT's preferred ones; "
                        "for uGMRT the recipe's own)")
    parser.add_argument("--no-primary-beam", action="store_true",
                        help="leave the images as the telescope saw the "
                        "sky: for a band without a measured beam model "
                        "(MeerKAT S), where it is required")
    parser.add_argument("--pbmin", type=float, default=0.3)
    parser.add_argument("--threads", type=int,
                        default=int(os.environ.get("SLURM_CPUS_PER_TASK")
                                    or 6),
                        help="threads for WSClean (default: the CPUs a SLURM "
                        "job was given, else 6)")
    parser.add_argument("--memory-gb", type=float, default=20.0,
                        help="memory WSClean may use, in GB (default 20)")
    parser.add_argument("--parallel-deconvolution", type=int,
                        metavar="PIXELS",
                        help="clean the images in sub-images of this size at "
                        "once (off unless given)")
    parser.add_argument(
        "--wsclean",
        default="singularity exec /software/astro/containers/"
        "oxkat-latest.sif wsclean",
    )
    parser.add_argument("--casa", help="the command that runs a Python with "
                        "casatasks; the scripts' own choice if not given")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the steps and stop")
    args = parser.parse_args(argv)
    if args.average_channels < 1:
        parser.error("--average-channels must be 1 or more")
    too_natural = [
        robust for robust in args.extra_robust or []
        if robust > MOST_UNTAPERED_ROBUST
    ]
    if too_natural:
        parser.error(
            f"--extra-robust {too_natural[0]:g}: untapered, a robust above "
            f"{MOST_UNTAPERED_ROBUST:g} weights the central square so far "
            "above the arms that the clean diverges (all nine Band 5 targets, "
            "2026-10-26); give --extra-taper for a larger beam instead")
    if any(taper <= 0 for taper in args.extra_taper or []):
        parser.error("--extra-taper takes beam sizes above 0 arcsec")
    if args.taper_field_scale < 1:
        parser.error("--taper-field-scale must be 1 or more")
    resolve_telescope(args, parser.error)
    print(
        f"{args.telescope} {args.band}: {args.imsize} pixels of "
        f"{float(args.cell):g}\", channels {args.channel_range[0]} to "
        f"{args.channel_range[1] - 1}"
        + (f", reference antennas {args.refant}" if args.refant else ""),
        flush=True,
    )
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
