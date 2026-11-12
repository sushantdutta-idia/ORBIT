#!/usr/bin/env python3
"""Decide which scans of a target go into its image.

A target is observed in scans hours apart, and they are not alike.
G0212-0335 has three: 9 minutes that flagging left of its first scan, a
whole one of 22 minutes, and 78 seconds of its last. Imaged alone, the
whole scan comes to 25 uJy/beam. With the other two the same recipe stops
at 40 (finding 130): the short one cannot be calibrated, half of its
solutions fail on signal-to-noise and the rest fit noise, and the first
one is 1.7 times as noisy as its length allows, which no gain per antenna
repairs. What they add, a floor 7% lower, is not worth what they cost.

So each scan is measured before any image is made, from the data alone.
The difference of the two parallel hands, (RR - LL) / 2, holds little of
the sky: the receivers' noise, whatever else was added to the data, and
what the two hands' beams see differently of the sources away from the
centre. Its scatter per visibility says how noisy a scan is for the amount
of data it has, and its scatter and weights together say what noise the
scan would give an image, and so what share of the target's weight it
holds.

The figures compare the scans of one target, which see the same sky. They
do not compare targets: on the nine fields of this observation the scatter
of whole scans runs from 58 to 171 mJy per visibility, alike in both scans
of a field, which receivers do not do and bright sources do (finding 130).
In such a field the noise given here is above what the receivers set.

Two rules, and a scan is left out if either holds:

- it holds less than ``--least-share`` of the weight: too little for
  solutions of its own, and the image loses little without it;
- it is more than ``--most-excess`` times as noisy, for the same amount of
  data, as the scan that holds the most weight: it carries more than
  receiver noise.

The scan that holds the most weight is always kept. ``--scans`` names the
scans to use and overrules both; what the rules would have done is still
written down.

If every scan is kept, the Measurement Set to go on with is the one given.
Otherwise CASA splits the kept scans into a new one in the work folder.
The Measurement Set given is read, never written to.

Example::

    python3 select_scans.py --target G0212-0335 \\
        --ms uGMRT_Products/calibration/calibrated/G0212-0335.ms \\
        --workdir targets_run/G0212-0335/scans --channel-range 51 768
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
import time
from pathlib import Path


REPORT_NAME = "scans_report.json"
# Stokes codes of the parallel hands: RR and LL, or XX and YY.
PARALLEL_HANDS = ((5, 8), (9, 12))
# Visibilities read at a time, per correlation.
CHUNK_SAMPLES = 4_000_000


def load_sibling(filename, name):
    spec = importlib.util.spec_from_file_location(
        name, Path(SCRIPT_DIR) / filename
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())


def read_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Under CASA
# ---------------------------------------------------------------------------

def parallel_hands(corr_types):
    """Where the two parallel hands are among the correlations."""

    codes = [int(code) for code in corr_types]
    for first, second in PARALLEL_HANDS:
        if first in codes and second in codes:
            return codes.index(first), codes.index(second)
    raise SystemExit(
        f"The correlations {codes} do not hold both parallel hands; the "
        "noise of a scan is measured from their difference."
    )


def data_column(names):
    """The column an imager reads: corrected data if there are any."""

    return "CORRECTED_DATA" if "CORRECTED_DATA" in names else "DATA"


def casa_measure(args):
    """Each scan's amount of data and the scatter of (RR - LL) / 2."""

    import numpy as np
    from casatools import table

    ms = str(args.ms)
    first, end = int(args.channel_range[0]), int(args.channel_range[1])
    polarization = table()
    polarization.open(ms + "/POLARIZATION")
    one, other = parallel_hands(polarization.getcell("CORR_TYPE", 0))
    polarization.close()

    main = table()
    main.open(ms)
    names = main.colnames()
    column = data_column(names)
    if len(set(int(v) for v in main.getcol("DATA_DESC_ID"))) > 1:
        raise SystemExit(
            f"{ms} holds more than one spectral window; --channel-range "
            "names channels of one."
        )
    spectrum = (
        "WEIGHT_SPECTRUM" in names and main.nrows() > 0
        and main.iscelldefined("WEIGHT_SPECTRUM", 0)
    )
    shape = main.getcell("FLAG", 0).shape
    correlations, channels = int(shape[0]), int(shape[1])
    if not 0 <= first < end <= channels:
        raise SystemExit(
            f"--channel-range {first} {end} does not fit {channels} channels."
        )
    blc, trc = [0, first], [correlations - 1, end - 1]
    chunk = max(1, CHUNK_SAMPLES // (end - first))

    scans = []
    for number in sorted(set(int(v) for v in main.getcol("SCAN_NUMBER"))):
        part = main.query(f"SCAN_NUMBER=={number} && ANTENNA1!=ANTENNA2")
        rows = part.nrows()
        times = part.getcol("TIME") if rows else np.zeros(0)
        first_antenna = part.getcol("ANTENNA1") if rows else np.zeros(0, int)
        second_antenna = part.getcol("ANTENNA2") if rows else np.zeros(0, int)
        interval = part.getcol("INTERVAL") if rows else np.zeros(0)
        row_flag = part.getcol("FLAG_ROW") if rows else np.zeros(0, bool)
        sums = {"samples": 0, "sum_w": 0.0, "sum_w2": 0.0,
                "sum_w2v2": 0.0, "sum_v2": 0.0}
        flagged, cells = 0, 0
        medians, with_data, antennas = [], set(), set()
        for start in range(0, rows, chunk):
            count = min(chunk, rows - start)
            data = part.getcolslice(column, blc, trc, [1, 1], start, count)
            flag = part.getcolslice("FLAG", blc, trc, [1, 1], start, count)
            flag = flag | row_flag[None, None, start:start + count]
            if spectrum:
                weight = part.getcolslice(
                    "WEIGHT_SPECTRUM", blc, trc, [1, 1], start, count)
            else:
                weight = np.broadcast_to(
                    part.getcol("WEIGHT", start, count)[:, None, :],
                    data.shape,
                )
            both = np.stack([flag[one], flag[other]])
            flagged += int(both.sum())
            cells += int(both.size)
            mean = 0.5 * (weight[one] + weight[other]).astype(float)
            use = ~flag[one] & ~flag[other] & (mean > 0)
            if not use.any():
                continue
            difference = 0.5 * (data[one] - data[other])
            power = (difference.real.astype(float) ** 2
                     + difference.imag.astype(float) ** 2)[use]
            mean = mean[use]
            sums["samples"] += int(use.sum())
            sums["sum_w"] += float(mean.sum())
            sums["sum_w2"] += float((mean ** 2).sum())
            sums["sum_w2v2"] += float((mean ** 2 * power).sum())
            sums["sum_v2"] += float(power.sum())
            medians.append(float(np.median(np.abs(difference.real[use]))))
            rows_used = use.any(axis=0)
            with_data.update(
                np.round(times[start:start + count][rows_used], 3).tolist())
            antennas.update(
                first_antenna[start:start + count][rows_used].tolist())
            antennas.update(
                second_antenna[start:start + count][rows_used].tolist())
        part.close()
        step = float(np.median(interval)) if rows else 0.0
        scans.append({
            "scan": number,
            "rows": int(rows),
            "start_mjd_s": float(times.min()) if rows else None,
            "end_mjd_s": float(times.max()) if rows else None,
            "integrations": len(with_data),
            "minutes": len(with_data) * step / 60.0,
            "antennas": len(antennas),
            "flagged": flagged / cells if cells else 1.0,
            # 1.4826 times the median distance from zero is the scatter of
            # a normal distribution; rare bright interference barely moves
            # it, and it moves the plain scatter a great deal.
            "robust_jy": (
                1.4826 * float(np.median(medians)) if medians else None),
            **sums,
        })
        print(
            f"scan {number}: {scans[-1]['minutes']:.1f} min, "
            f"{sums['samples']} visibilities", flush=True,
        )
    main.close()

    with open(args.out, "w") as handle:
        json.dump({
            "column": column,
            "weights": "WEIGHT_SPECTRUM" if spectrum else "WEIGHT",
            "scans": scans,
        }, handle, indent=2)
    return 0


def casa_split(args):
    """The kept scans, as a Measurement Set of this run's own."""

    from casatasks import split
    from casatools import table

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    main = table()
    main.open(str(args.ms))
    column = data_column(main.colnames())
    main.close()
    wanted = ",".join(str(int(scan)) for scan in args.scans)
    average = int(getattr(args, "average_channels", 1) or 1)
    print(
        f"Splitting scans {wanted} of {args.ms} into {out}"
        + (f", averaging every {average} channels" if average > 1 else "")
        + " ...",
        flush=True,
    )
    options = {}
    if average > 1:
        # split averages with the weights the channels carry.
        options["width"] = average
    split(
        vis=str(args.ms), outputvis=str(out), scan=wanted,
        datacolumn="corrected" if column == "CORRECTED_DATA" else "data",
        **options,
    )
    if not out.is_dir():
        raise SystemExit("No Measurement Set was written.")
    print(f"Written: {out}")
    return 0


# ---------------------------------------------------------------------------
# What the measurements say
# ---------------------------------------------------------------------------

def floor(scans):
    """The noise an image of these scans has from the data alone (Jy/beam).

    For natural weights and the weights the data carry, which is how an
    imager adds the scans up: not the best that could be made of them.
    """

    weight = sum(scan["sum_w"] for scan in scans)
    if weight <= 0:
        return None
    return math.sqrt(sum(scan["sum_w2v2"] for scan in scans) / 2.0) / weight


def derived(scans):
    """Add to each scan its scatter, its noise in an image, and its share.

    ``rms_jy`` is the scatter of one visibility; ``noise_jy`` what the
    scan alone would give an image; ``share`` its part of the target's
    weight if every scan were weighted by its noise; ``excess`` how many
    times as noisy it is, for the same amount of data, as the scan with
    the largest share; ``outliers`` the plain scatter over the robust one.
    """

    for scan in scans:
        usable = (
            scan.get("samples", 0) > 0 and scan["sum_w"] > 0
            and scan["sum_w2"] > 0 and scan["sum_w2v2"] > 0
        )
        scan["rms_jy"] = (
            math.sqrt(scan["sum_w2v2"] / (2.0 * scan["sum_w2"]))
            if usable else None
        )
        scan["noise_jy"] = floor([scan]) if usable else None
        plain = (
            math.sqrt(scan["sum_v2"] / (2.0 * scan["samples"]))
            if usable else None
        )
        scan["outliers"] = (
            plain / scan["robust_jy"]
            if plain and scan.get("robust_jy") else None
        )
    measured = [scan for scan in scans if scan["noise_jy"]]
    total = sum(scan["noise_jy"] ** -2 for scan in measured)
    for scan in scans:
        scan["share"] = (
            scan["noise_jy"] ** -2 / total if scan["noise_jy"] else 0.0
        )
    reference = max(measured, key=lambda scan: scan["share"], default=None)
    for scan in scans:
        scan["excess"] = (
            scan["rms_jy"] / reference["rms_jy"] if scan["rms_jy"] else None
        )
    return reference["scan"] if reference else None


def choose(scans, least_share, most_excess, forced=None):
    """Mark each scan kept or left out, with the reason; return those kept.

    ``derived`` must have been applied. ``rule`` is what the two rules say
    of a scan and ``kept`` what is done, which differ only for scans named
    with ``--scans``.
    """

    numbers = [scan["scan"] for scan in scans]
    if forced is not None:
        unknown = [number for number in forced if number not in numbers]
        if unknown:
            raise SystemExit(
                f"No scan {', '.join(str(n) for n in unknown)} in these "
                f"data; there are: {', '.join(str(n) for n in numbers)}."
            )
    most = max(scans, key=lambda scan: scan["share"], default=None)
    for scan in scans:
        reason = None
        if not scan["noise_jy"]:
            reason = "it holds no data"
        elif scan is most:
            pass
        elif scan["share"] < least_share:
            reason = (
                f"it holds {100 * scan['share']:.0f}% of the weight, too "
                "little for solutions of its own"
            )
        elif scan["excess"] > most_excess:
            reason = (
                f"it is {scan['excess']:.2f} times as noisy as scan "
                f"{most['scan']} for the same amount of data"
            )
        scan["rule"] = reason
        if forced is None:
            scan["kept"], scan["why"] = reason is None, reason
        elif scan["scan"] in forced:
            scan["kept"], scan["why"] = True, None
            if not scan["noise_jy"]:
                raise SystemExit(f"Scan {scan['scan']} holds no data.")
        else:
            scan["kept"], scan["why"] = False, "it was not named with --scans"
    return [scan["scan"] for scan in scans if scan["kept"]]


def micro(value):
    return "?" if value is None else f"{value * 1e6:.1f}"


def reading(target, scans, forced=False):
    """The scans in words: one line each and what leaving out costs."""

    lines = []
    main = max(scans, key=lambda scan: scan["share"], default=None)
    usual = main["outliers"] if main else None
    for scan in scans:
        if not scan["noise_jy"]:
            lines.append(f"scan {scan['scan']}: no data; left out.")
            continue
        line = (
            f"scan {scan['scan']}: {scan['minutes']:.1f} min, "
            f"{scan['antennas']} antennas, {100 * scan['flagged']:.0f}% "
            f"flagged; alone {micro(scan['noise_jy'])} uJy/beam, "
            f"{100 * scan['share']:.0f}% of the weight, "
            f"{scan['excess']:.2f} times the noise of the main scan for "
            "the same amount of data"
        )
        if usual and scan["outliers"] and scan["outliers"] > 1.5 * usual:
            # Channels and antennas differ in noise, so the plain scatter
            # of a clean scan is already 1.2 to 1.6 times its bulk's; what
            # tells is a scan far above the main one.
            line += (
                f" (scatter {scan['outliers']:.1f} times what its bulk "
                f"has, against {usual:.1f} in the main scan: interference "
                "left in)"
            )
        if scan["kept"]:
            line += "; kept"
            if scan["rule"]:
                line += f", named with --scans though {scan['rule']}"
        else:
            line += f"; left out: {scan['why']}"
            if forced:
                # What the rules say of a scan that was named out is the
                # test of the rules.
                line += (
                    f", and {scan['rule']}" if scan["rule"]
                    else ", though the rules would keep it"
                )
        lines.append(line + ".")
    kept = [scan for scan in scans if scan["kept"]]
    everything = floor([scan for scan in scans if scan["noise_jy"]])
    left = floor(kept) if kept else None
    if not kept:
        lines.append(f"{target}: no scan is kept.")
    elif len(kept) == len(scans):
        lines.append(
            f"{target}: every scan is kept; the data allow "
            f"{micro(everything)} uJy/beam with natural weights."
        )
    else:
        lines.append(
            f"{target}: kept "
            + ", ".join(f"scan {scan['scan']}" for scan in kept)
            + f"; the data then allow {micro(left)} uJy/beam with natural "
            f"weights, against {micro(everything)} with every scan."
        )
    return lines


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

def casa_command(args, part, out, scans=()):
    command = shlex.split(args.casa) + [
        str(Path(SCRIPT_DIR) / "select_scans.py"),
        "--casa-part", part, "--ms", str(args.ms), "--out", str(out),
        "--channel-range",
        str(int(args.channel_range[0])), str(int(args.channel_range[1])),
    ]
    if scans:
        command += ["--scans", *(str(int(scan)) for scan in scans)]
    average = int(getattr(args, "average_channels", 1) or 1)
    if part == "split" and average > 1:
        command += ["--average-channels", str(average)]
    return command


def find_casa(args, recipe):
    """The command that runs CASA: the one given, or the first that answers."""

    candidates = [args.casa] if args.casa else list(recipe.CASA_CANDIDATES)
    failures = []
    for casa in candidates:
        words = shlex.split(casa)
        if not words or shutil.which(words[0]) is None:
            failures.append(
                f"'{casa}': '{words[0] if words else ''}' is not on PATH"
            )
            continue
        ok, detail = recipe.casa_answers(casa, args.ms)
        if ok:
            return casa
        failures.append(f"'{casa}': {detail}")
    raise SystemExit(
        "CASA: no command ran a Python that imports casatasks and opens "
        f"{args.ms}. Pass --casa with one that does. Tried:\n    "
        + "\n    ".join(failures)
    )


def run(args):
    recipe = load_sibling("selfcal_recipe.py", "select_scans_recipe")
    imager = recipe.IMAGER
    source = Path(str(args.ms).rstrip("/"))
    if not source.is_dir():
        raise SystemExit(f"No Measurement Set at {source}.")
    workdir = Path(args.workdir)
    if imager.inside(workdir, args.outdir):
        raise SystemExit(
            f"{workdir} is inside the pipeline's products ({args.outdir}). "
            "Give --workdir a folder of its own."
        )
    args.ms = source
    workdir.mkdir(parents=True, exist_ok=True)
    with imager.claimed(workdir):
        report_path = workdir / REPORT_NAME
        report = read_json(report_path) or {}
        if report.get("input") != str(source):
            report = {"target": args.target, "input": str(source)}
        report["channel_range"] = [int(v) for v in args.channel_range]

        def save():
            report["updated"] = now()
            scratch = report_path.with_suffix(".tmp")
            with open(scratch, "w") as handle:
                json.dump(report, handle, indent=2)
            scratch.replace(report_path)

        (workdir / "logs").mkdir(exist_ok=True)
        if "measured" not in report:
            args.casa = find_casa(args, recipe)
            measured = workdir / "measured.json"
            imager.run_command(
                casa_command(args, "measure", measured),
                workdir / "logs" / "measure.log",
            )
            found = read_json(measured)
            if not found or not found.get("scans"):
                raise SystemExit(f"No scans were measured in {source}.")
            report["measured"] = found
            save()

        scans = [dict(scan) for scan in report["measured"]["scans"]]
        report["main_scan"] = derived(scans)
        forced = [int(scan) for scan in args.scans] if args.scans else None
        kept = choose(scans, args.least_share, args.most_excess, forced)
        if not kept:
            raise SystemExit(f"{args.target}: no scan holds data.")
        average = int(getattr(args, "average_channels", 1) or 1)
        report["rule"] = {
            "least_share": float(args.least_share),
            "most_excess": float(args.most_excess),
            "scans": forced,
        }
        if average > 1:
            report["rule"]["average_channels"] = average
        report["scans"] = scans
        report["kept"] = kept
        report["left_out"] = [
            scan["scan"] for scan in scans if not scan["kept"]
        ]
        measured = [scan for scan in scans if scan["noise_jy"]]
        report["noise_every_scan_jy"] = floor(measured)
        report["noise_kept_jy"] = floor(
            [scan for scan in scans if scan["kept"]])
        report["reading"] = reading(args.target, scans, forced is not None)

        if not report["left_out"] and average == 1:
            report["measurement_set"] = str(source)
        else:
            # Scans left out, or channels averaged: a copy of this run's own.
            out = workdir / f"{args.target}.ms"
            wanted = (
                kept if average == 1
                else {"scans": kept, "average_channels": average}
            )
            if report.get("split") != wanted or not out.is_dir():
                report.pop("split", None)
                report.pop("measurement_set", None)
                save()
                args.casa = args.casa or find_casa(args, recipe)
                imager.run_command(
                    casa_command(args, "split", out, kept),
                    workdir / "logs" / "split.log",
                )
                if not out.is_dir():
                    raise SystemExit(f"The split left no {out}.")
                report["split"] = wanted
            report["measurement_set"] = str(out)
        save()

    print("\n".join(report["reading"]))
    print(f"Data to go on with: {report['measurement_set']}")
    print(f"Report: {report_path}")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", help="the target's name")
    parser.add_argument("--ms", required=True,
                        help="its calibrated data; read, never written to")
    parser.add_argument("--workdir",
                        help="where the report and, if scans are left out, "
                        "the Measurement Set of the kept ones go")
    parser.add_argument("--channel-range", type=int, nargs=2, required=True,
                        metavar=("FIRST", "END"),
                        help="channels that hold data, END exclusive")
    parser.add_argument("--least-share", type=float, default=0.10,
                        help="a scan with less than this part of the "
                        "target's weight is left out (default 0.10)")
    parser.add_argument("--most-excess", type=float, default=1.4,
                        help="a scan more than this many times as noisy as "
                        "the main scan, for the same amount of data, is "
                        "left out (default 1.4)")
    parser.add_argument("--scans", type=int, nargs="+",
                        help="use these scans, whatever the rules say")
    parser.add_argument("--average-channels", type=int, default=1,
                        metavar="N",
                        help="average every N channels in the copy the "
                        "imaging goes on with (default 1: none); the copy "
                        "is then made even if every scan is kept")
    parser.add_argument("--outdir", default="uGMRT_Products",
                        help="the pipeline's output folder, which the work "
                        "folder may not be in")
    parser.add_argument("--casa", help="the command that runs a Python with "
                        "casatasks; the first of selfcal_recipe.py's that "
                        "answers if not given")
    parser.add_argument("--casa-part", choices=("measure", "split"),
                        help=argparse.SUPPRESS)
    parser.add_argument("--out", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.casa_part == "measure":
        return casa_measure(args)
    if args.casa_part == "split":
        return casa_split(args)
    if not args.target or not args.workdir:
        parser.error("--target and --workdir are required")
    if args.average_channels < 1:
        parser.error("--average-channels must be 1 or more")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
