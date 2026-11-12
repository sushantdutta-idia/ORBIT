#!/usr/bin/env python3
"""Take the bright off-axis sources out of a field, one after another.

For each source, in the order given:

1. ``experiment_extract_source.py`` leaves it alone in a copy of the data;
2. ``selfcal_recipe.py`` calibrates towards it, with a small image;
3. ``experiment_extract_source.py --peel`` subtracts it with that
   calibration and gives the field back without it.

The next source is taken from what the last one left. On G0212-0335 one
gain per antenna for the whole field left holes of 1 to 3% of the peak
beside every bright source away from the pointing centre; towards such a
source the antennas each have their own gain, more so the higher the
frequency (findings 123 to 125), and a source subtracted with its own
gains takes its holes with it (finding 126).

This script needs neither WSClean nor CASA; it runs the three scripts
above, which do. It writes nothing to the data it is given. A run that is
interrupted continues with the first source not yet taken out when the
same command is given again.

Each source is ``X Y CORES BINS``: its pixel in the model images, the
threshold of its cores-only model in mJy/beam (about an eighth of its
peak), and the channels per gain of its bandpass round, wider for a
fainter source. ``BINS`` of 0 is for a source too faint for delays or a
bandpass at all: only a phase and an amplitude per antenna are solved
towards it. Whatever a solve towards one source cannot solve is left as it
is, never flagged: the data are the whole field's.

Example::

    python3 peel_sources.py --target G0212-0335 \\
        --ms wsclean_test/G0212-0335_scan15_ph_r2_r3_r4_r5_bpa.ms \\
        --models wsclean_test_scan15_bpa_ch12/G0212-0335/wsclean_robust_0_I \\
        --channel-range 51 768 \\
        --source 3162 4227 6 30 --source 1641 2103 4 60
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
import shlex
import shutil
import time
from pathlib import Path


REPORT_NAME = "peel_report.json"
# The disc taken out with a source grows with it: a source as bright as
# this keeps the disc it is given (--radius), a brighter one a disc larger
# by the square root of how much brighter, up to MOST_DISC times as large.
# On G0213-0251, G0216-0604 and G0218-0245 sources of 131 to 257 mJy/beam
# left rings of -0.4 to -1.6 mJy/beam 96 to 180 pixels out, just beyond a
# disc of 90: the field's model held that source's errors there, and they
# were given back as the rest of the field (finding 136).
DISC_FOR_MJY = 30.0
MOST_DISC = 3.0


def load_imager():
    spec = importlib.util.spec_from_file_location(
        "peel_sources_imager", Path(SCRIPT_DIR) / "experiment_wsclean_image.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


IMAGER = load_imager()


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())


def script(name):
    return [sys.executable, "-u", str(Path(SCRIPT_DIR) / name)]


def parallel_options(args):
    """Sub-image size for parallel cleaning, when one is given."""

    if getattr(args, "parallel_deconvolution", None) is None:
        return []
    return ["--parallel-deconvolution", str(int(args.parallel_deconvolution))]


def shared(args):
    """What every step is told alike: channels, tools, the products."""

    options = ["--outdir", str(args.outdir)]
    if args.channel_range:
        options += [
            "--channel-range",
            str(int(args.channel_range[0])), str(int(args.channel_range[1])),
        ]
    options += ["--wsclean", args.wsclean]
    if args.casa:
        options += ["--casa", args.casa]
    return options


def radius_of(args, source):
    """The disc taken out with a source, in pixels."""

    return int(source.get("radius") or args.radius)


def without_options(args, gone):
    """The sources already out, each with its own disc."""

    options = []
    for other in gone:
        options += ["--also-without", str(other["x"]), str(other["y"]),
                    str(radius_of(args, other))]
    return options


def extract_command(args, ms, name, source, gone):
    return script("experiment_extract_source.py") + [
        "--ms", str(ms), "--models", str(args.models),
        "--pixel", str(source["x"]), str(source["y"]),
        "--radius", str(radius_of(args, source)),
        "--name", name, "--workdir", str(args.workdir), "--overwrite",
        "--threads", str(int(args.threads)),
        "--memory-gb", f"{float(args.memory_gb):g}",
        *shared(args),
    ] + without_options(args, gone)


def recipe_command(args, name, source):
    """Self-calibration towards the one source, in a small image."""

    return script("selfcal_recipe.py") + [
        "--target", args.target,
        "--ms", str(Path(args.workdir) / name / f"{name}.ms"),
        "--workdir", str(Path(args.workdir) / f"recipe_{name}"),
        "--imsize", str(int(args.imsize)), "--cell", f"{float(args.cell):g}",
        "--channels-out", str(int(args.channels_out)),
        "--core-threshold", f"{float(source['cores']):g}",
        "--threads", str(int(args.threads)),
        "--memory-gb", f"{float(args.memory_gb):g}",
        *parallel_options(args),
        "--image-first",
        # The field is calibrated already: against a cores-only model a
        # source with anything around its cores has its phases moved, and
        # the image of every source so solved was no better for those
        # rounds, that of a 300 mJy double a third worse (finding 133).
        "--whole-field-only",
        # The data are the field's: a solve towards one source flags none.
        "--keep-unsolved",
        *(["--tidy"] if args.tidy else []),
        *(
            ["--bandpass-channels", str(int(source["bins"]))]
            if source["bins"] else ["--phase-amplitude-only"]
        ),
        # The recipe's own default names uGMRT antennas.
        *(["--refant", args.refant] if getattr(args, "refant", None) else []),
        *shared(args),
    ]


def peel_command(args, name, source, gone):
    return script("experiment_extract_source.py") + [
        "--peel", str(Path(args.workdir) / f"recipe_{name}" / args.target),
        "--models", str(args.models),
        "--pixel", str(source["x"]), str(source["y"]),
        "--radius", str(radius_of(args, source)),
        "--name", name, "--workdir", str(args.workdir), "--overwrite",
        "--threads", str(int(args.threads)),
        "--memory-gb", f"{float(args.memory_gb):g}",
        *shared(args),
    ] + without_options(args, gone)


def image_command(args, ms):
    """The field as it is left, imaged as the field was before."""

    command = script("experiment_wsclean_image.py") + [
        "--outdir", str(args.outdir), "--target", args.target,
        "--ms", str(ms),
        "--workdir", str(Path(args.workdir) / "image"),
        "--channels-out", str(int(args.channels_out)),
        "--imsize", str(int(args.field_imsize)),
        "--cell", f"{float(args.cell):g}",
        "--scales", "0", "5", "15", "30", "--no-local-rms",
        "--auto-threshold", "1.0", "--nmiter", "20",
        "--threads", str(int(args.threads)),
        "--memory-gb", f"{float(args.memory_gb):g}",
        *parallel_options(args),
        "--wsclean", args.wsclean, "--overwrite",
    ]
    if args.channel_range:
        command += [
            "--channel-range",
            str(int(args.channel_range[0])), str(int(args.channel_range[1])),
        ]
    return command


def disc_radius(peak_mjy, offset_pixels, radius, keep_pixels):
    """The disc taken out with a source of this peak, in pixels.

    ``radius`` for a source of ``DISC_FOR_MJY`` or fainter, larger by the
    square root of how much brighter, at most ``MOST_DISC`` times it, and
    never reaching closer to the pointing centre than ``keep_pixels``,
    where the target is: a disc takes out all it holds.
    """

    grown = float(radius) * max(1.0, float(peak_mjy) / DISC_FOR_MJY) ** 0.5
    grown = min(grown, MOST_DISC * float(radius))
    grown = min(grown, float(offset_pixels) - float(keep_pixels))
    return int(round(max(float(radius), grown)))


def choose_sources(plane, cell_arcsec, above_mjy, keep_arcmin, most,
                   radius):
    """The bright sources away from the pointing centre, brightest first.

    ``plane`` is an image as the telescope saw it, indexed [x, y]. A source
    is taken if its peak is at least ``above_mjy`` and it lies more than
    ``keep_arcmin`` from the centre, where the target is and stays. Each
    comes with the numbers its calibration needs: the threshold of its
    cores-only model, an eighth of its peak; the channels per gain of its
    bandpass round, by how bright it is; and the disc taken out with it
    (``disc_radius``). Below 20 mJy/beam there is no bandpass round and no
    delay: on G0212-0335 those gave an 11 mJy source a delay of 18 ns and
    amplitudes that ran by a half. A source inside the disc of a brighter
    one goes out with it and is not taken on its own.
    """

    nx, ny = plane.shape
    keep_pixels = float(keep_arcmin) * 60.0 / float(cell_arcsec)
    chosen = []
    for x, y in IMAGER.brightest_sources(
        plane, 6 * int(most) + 12, separation=2 * int(radius)
    ):
        peak = float(plane[x, y]) * 1e3
        if peak < float(above_mjy):
            break
        offset = ((x - nx // 2) ** 2 + (y - ny // 2) ** 2) ** 0.5
        if offset <= keep_pixels:
            continue
        if any(
            ((x - other["x"]) ** 2 + (y - other["y"]) ** 2) ** 0.5
            <= other["radius"]
            for other in chosen
        ):
            continue
        chosen.append({
            "x": int(x), "y": int(y),
            "cores": round(peak / 8.0, 2),
            "bins": 30 if peak >= 40.0 else 60 if peak >= 20.0 else 0,
            "radius": disc_radius(peak, offset, radius, keep_pixels),
            "peak_mJy": round(peak, 2),
            "offset_arcmin": round(offset * float(cell_arcsec) / 60.0, 1),
        })
        if len(chosen) >= int(most):
            break
    return chosen


def read_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def source_findings(args, name):
    """From the recipe's report: the source before and after."""

    report = read_json(
        Path(args.workdir) / f"recipe_{name}" / args.target
        / "recipe_report.json"
    ) or {}

    def of(image):
        image = image or {}
        return {
            "peak_Jy_beam": image.get("peak_Jy_beam"),
            "most_negative_Jy_beam": image.get("most_negative_Jy_beam"),
            "image_rms_Jy_beam": image.get("image_rms_Jy_beam"),
        }

    return {
        "before": of(report.get("before")),
        "after": of(report.get("final")),
        "rounds": report.get("reading") or [],
    }


def run(args):
    source_ms = Path(str(args.ms).rstrip("/"))
    if not source_ms.is_dir():
        raise SystemExit(f"No Measurement Set at {source_ms}.")
    workdir = Path(args.workdir)
    if IMAGER.inside(workdir, args.outdir):
        raise SystemExit(
            f"{workdir} is inside the pipeline's products ({args.outdir}). "
            "Give --workdir a folder of its own."
        )
    sources = [
        {"x": int(x), "y": int(y), "cores": float(cores), "bins": int(bins)}
        for x, y, cores, bins in args.source or []
    ]
    if args.auto_sources:
        if sources:
            raise SystemExit("Give --source or --auto-sources, not both.")
        image = Path(args.auto_sources)
        if not image.is_file():
            raise SystemExit(f"No image at {image}.")
        sources = choose_sources(
            IMAGER.fits_plane(image), args.cell, args.peel_above,
            args.keep_radius, args.most, args.radius,
        )
        print(
            f"{len(sources)} sources of at least {args.peel_above:g} "
            f"mJy/beam more than {args.keep_radius:g} arcmin from the "
            "centre: " + (", ".join(
                f"{source['peak_mJy']:g} mJy at {source['x']}, {source['y']} "
                f"({source['offset_arcmin']:g} arcmin, a disc of "
                f"{radius_of(args, source)} pixels)" for source in sources
            ) or "none") + ".",
            flush=True,
        )
    elif not sources:
        raise SystemExit("Give --source or --auto-sources.")
    already = [
        {"x": int(x), "y": int(y)} for x, y in args.already_out or []
    ]
    # Numbered on from the sources that were taken out before this run.
    names = [
        f"source{len(already) + index}"
        for index in range(1, len(sources) + 1)
    ]

    if args.dry_run:
        current, gone = source_ms, list(already)
        for name, source in zip(names, sources):
            for command in (
                extract_command(args, current, name, source, gone),
                recipe_command(args, name, source),
                peel_command(args, name, source, gone),
            ):
                print(shlex.join(command))
            gone = gone + [source]
            current = workdir / name / f"{name}_peeled.ms"
        if not args.no_image:
            print(shlex.join(image_command(args, current)))
        return None

    workdir.mkdir(parents=True, exist_ok=True)
    with IMAGER.claimed(workdir):
        report_path = workdir / REPORT_NAME
        report = read_json(report_path) or {
            "input": str(source_ms), "sources": {},
        }
        if report.get("input") != str(source_ms):
            raise SystemExit(
                f"{workdir} holds a run on {report.get('input')}, not "
                f"{source_ms}. Give --workdir another folder."
            )
        (workdir / "logs").mkdir(exist_ok=True)

        def save():
            report["updated"] = now()
            scratch = report_path.with_suffix(".tmp")
            with open(scratch, "w") as handle:
                json.dump(report, handle, indent=2)
            scratch.replace(report_path)

        def step(name, what, command):
            print(f"\n=== {name}: {what} ({now()}) ===", flush=True)
            report["running"] = {"step": f"{name}: {what}", "since": now()}
            save()
            try:
                IMAGER.run_command(command, workdir / "logs" / f"{name}_{what}.log")
            except IMAGER.CommandFailed as failure:
                report["running"] = None
                report["failed"] = {
                    "step": f"{name}: {what}", "error": str(failure),
                }
                save()
                raise SystemExit(
                    f"\nStep '{name}' ({what}) failed: {failure}.\nThe "
                    "sources before it are kept; the same command continues "
                    "with this one."
                ) from None
            report.pop("failed", None)
            report["running"] = None

        current, gone = source_ms, list(already)
        for name, source in zip(names, sources):
            peeled = workdir / name / f"{name}_peeled.ms"
            entry = report["sources"].setdefault(name, dict(source))
            # With --tidy only the last source's field is kept, so a source
            # is out if its field is there or a later source is out.
            later = any(
                report["sources"].get(other, {}).get("peeled")
                for other in names[names.index(name) + 1:]
            )
            if entry.get("peeled") and (peeled.is_dir() or later):
                print(f"{name} is already out.", flush=True)
            else:
                step(name, "extract",
                     extract_command(args, current, name, source, gone))
                step(name, "recipe", recipe_command(args, name, source))
                step(name, "peel", peel_command(args, name, source, gone))
                if not peeled.is_dir():
                    raise SystemExit(
                        f"\n'{name}' was not taken out: {peeled} is missing."
                    )
                entry.update(source_findings(args, name))
                entry["peeled"] = str(peeled)
                entry["finished"] = now()
                save()
                if args.tidy:
                    # The source alone, the copies its calibration worked
                    # on, and the field as the source before it left it.
                    stale = [workdir / name / f"{name}.ms"]
                    stale += sorted(
                        (workdir / f"recipe_{name}" / args.target)
                        .glob("round*.ms")
                    )
                    if current != source_ms:
                        stale.append(current)
                    for path in stale:
                        shutil.rmtree(path, ignore_errors=True)
            gone = gone + [source]
            if peeled.is_dir():
                current = peeled

        report["measurement_set"] = str(current)
        report["chosen"] = sources
        if not args.no_image and report.get("imaged") != str(current):
            step("field", "image", image_command(args, current))
            report["imaged"] = str(current)
        save()
    print(f"\nThe field without {len(sources)} sources: {current}\n"
          f"Report: {report_path}")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", required=True)
    parser.add_argument("--ms", required=True,
                        help="the calibrated data of the field; copied, "
                        "never written to")
    parser.add_argument("--models", required=True,
                        help="the name WSClean was given for the image of "
                        "these data, as in <name>-0000-model.fits")
    parser.add_argument("--auto-sources", metavar="IMAGE",
                        help="instead of --source, take the sources from "
                        "this image of the field (the MFS image of the run "
                        "whose models are given)")
    parser.add_argument("--peel-above", type=float, default=10.0,
                        metavar="MJY",
                        help="with --auto-sources: sources at least this "
                        "bright, as the telescope saw them (default 10 "
                        "mJy/beam)")
    parser.add_argument("--keep-radius", type=float, default=4.0,
                        metavar="ARCMIN",
                        help="with --auto-sources: nothing within this of "
                        "the pointing centre is taken out (default 4)")
    parser.add_argument("--most", type=int, default=8,
                        help="with --auto-sources: at most this many "
                        "sources (default 8)")
    parser.add_argument("--source", nargs=4, action="append",
                        metavar=("X", "Y", "CORES", "BINS"),
                        help="a source to take out, brightest first: its "
                        "pixel, the threshold of its cores-only model in "
                        "mJy/beam, and the channels per gain of its "
                        "bandpass round; may be given several times")
    parser.add_argument("--already-out", type=int, nargs=2, action="append",
                        metavar=("X", "Y"),
                        help="a source that --ms no longer holds, taken out "
                        "in an earlier run; may be given several times")
    parser.add_argument("--workdir", default="peel_field")
    parser.add_argument("--outdir", default="uGMRT_Products",
                        help="the pipeline's products, which are not touched")
    parser.add_argument("--radius", type=int, default=90,
                        help="pixels about each source (default 90)")
    parser.add_argument("--imsize", type=int, default=1024,
                        help="pixels of the image about one source")
    parser.add_argument("--field-imsize", type=int, default=5120,
                        help="pixels of the image of the field at the end")
    parser.add_argument("--cell", type=float, default=1.0, help="arcsec")
    parser.add_argument("--channels-out", type=int, default=12)
    parser.add_argument("--channel-range", type=int, nargs=2,
                        metavar=("FIRST", "END"))
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--memory-gb", type=float, default=20.0,
                        help="memory WSClean may use, in GB (default 20)")
    parser.add_argument("--parallel-deconvolution", type=int,
                        metavar="PIXELS",
                        help="clean the images in sub-images of this size "
                        "at once")
    parser.add_argument(
        "--wsclean",
        default="singularity exec /software/astro/containers/"
        "oxkat-latest.sif wsclean",
    )
    parser.add_argument("--casa", help="the command that runs a Python with "
                        "casatasks; the scripts' own choice if not given")
    parser.add_argument("--tidy", action="store_true",
                        help="once a source is out, remove its own data, "
                        "the copies its calibration worked on, and the "
                        "field as the source before it left it")
    parser.add_argument("--no-image", action="store_true",
                        help="do not image the field at the end")
    parser.add_argument("--refant",
                        help="reference antennas of the solves towards each "
                        "source, in order of preference (default: "
                        "selfcal_recipe.py's, which are uGMRT antennas)")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the steps and stop")
    run(parser.parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
