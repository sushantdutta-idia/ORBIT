#!/usr/bin/env python3
"""Solve one phase per antenna against a cores-only model, and apply it.

The second step of a test, after ``experiment_wsclean_image.py --core-model``
has written a model of the bright cores into MODEL_DATA.

Self-calibration cannot tell "the antennas carry fixed phase offsets" from
"every source in the field has the same lopsided shape": both give the same
visibilities. The model decides. A model cleaned deep contains the lopsided
shape, and gains solved against it keep the offsets; a model of compact
positive cores says the shape is not sky, and the solve then finds the
offsets. Because they are fixed, one solution per antenna for the whole scan
is enough, and it has the signal-to-noise of the whole scan.

The phases found are printed: if the explanation holds they are several
degrees apart, and if they are all near zero there were no offsets to find.
The corrected data are written to a new Measurement Set beside the input,
which the third step images.

With a shorter ``--solint`` it shows how each antenna's phase moves through
the scan, and with ``--calmode ap`` it solves amplitudes as well.

``--delays`` solves a delay per antenna first. One phase for the whole band
is right at the band's centre only: an antenna whose signal arrives a
fraction of a nanosecond early or late has a phase that runs across the
band, and after the single phase is removed it is left with errors of
opposite sign at the two ends. In the image that is a lopsided pattern
around every source which differs between the low and the high sub-bands.

``--bandpass N`` solves instead a gain for every N channels, one set for
each scan. A delay and a phase describe a phase that is a straight line
across the band. The ionosphere's part of an antenna's phase goes as one
over the frequency, and from 570 to 850 MHz that is not a straight line:
what a delay took for a slope of 228 degrees leaves a bow of some twenty
degrees in the middle of the band. The solve shows each antenna's phase
against frequency, and by default corrects the phases only. With
``--bandpass-amplitude`` each antenna's amplitude across the band is
corrected too, as its ratio to what all antennas share in that bin.

Run it with the Python that has CASA::

    python3 experiment_phase_offsets.py --ms wsclean_test/G0212-0335_scan15.ms
    python3 experiment_phase_offsets.py --ms copy_ph.ms --solint 1min
"""

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# This folder holds the pipeline's own inspect.py; see diagnose_noise_floor.
if sys.path and os.path.abspath(sys.path[0]) == SCRIPT_DIR:
    sys.path.pop(0)

import argparse
import contextlib
import json
import math
import shutil
import socket
import time
from pathlib import Path


def phase_table(antennas, values, flags, names):
    """Per-antenna phases in degrees from a gain table's columns.

    ``values`` and ``flags`` are indexed [polarization, channel, row] as the
    table tool returns them and ``antennas`` gives the antenna of each row.
    A flagged solution is reported as None. Phases are of the first
    solution interval of each antenna, which is the only one when the
    solve covers the scan.
    """

    rows = {}
    for row, antenna in enumerate(antennas):
        antenna = int(antenna)
        if antenna in rows:
            continue
        phases = []
        for pol in range(len(values)):
            value = complex(values[pol][0][row])
            if bool(flags[pol][0][row]) or value == 0:
                phases.append(None)
            else:
                phases.append(math.degrees(math.atan2(value.imag, value.real)))
        rows[antenna] = phases
    return [
        {
            "antenna": names[antenna] if antenna < len(names) else str(antenna),
            "id": antenna,
            "phase_deg": rows[antenna],
        }
        for antenna in sorted(rows)
    ]


def gain_series(antennas, times, values, flags, names):
    """Per antenna, what its solutions do through the observation.

    For each antenna and polarization: the mean phase, how far the phase
    moves about that mean (rms, degrees), the mean amplitude and its
    scatter, and how many solutions failed. A fixed offset shows as a mean;
    an antenna whose phase wanders shows as a large rms, and one correction
    per scan cannot follow it.
    """

    grouped = {}
    for row, antenna in enumerate(antennas):
        grouped.setdefault(int(antenna), []).append(row)
    series = []
    for antenna in sorted(grouped):
        rows = sorted(grouped[antenna], key=lambda row: float(times[row]))
        entry = {
            "antenna": names[antenna] if antenna < len(names) else str(antenna),
            "id": antenna,
            "solutions": len(rows),
            "polarizations": [],
        }
        for pol in range(len(values)):
            good = [
                complex(values[pol][0][row]) for row in rows
                if not bool(flags[pol][0][row])
                and complex(values[pol][0][row]) != 0
            ]
            if not good:
                entry["polarizations"].append({"failed": len(rows)})
                continue
            total = sum(value / abs(value) for value in good)
            mean = math.atan2(total.imag, total.real)
            deviations = [
                math.atan2(
                    math.sin(math.atan2(v.imag, v.real) - mean),
                    math.cos(math.atan2(v.imag, v.real) - mean),
                )
                for v in good
            ]
            amplitudes = [abs(v) for v in good]
            mean_amplitude = sum(amplitudes) / len(amplitudes)
            entry["polarizations"].append({
                "failed": len(rows) - len(good),
                # The largest single solution, which with a solution per
                # scan is not the mean: +80 degrees in one scan and -80 in
                # the next have a mean of nothing.
                "largest_phase_deg": max(
                    abs(math.degrees(math.atan2(v.imag, v.real)))
                    for v in good
                ),
                "largest_amplitude_step": max(
                    abs(a - 1.0) for a in amplitudes
                ),
                "mean_phase_deg": math.degrees(mean),
                "phase_rms_deg": math.degrees(math.sqrt(
                    sum(d * d for d in deviations) / len(deviations)
                )),
                "mean_amplitude": mean_amplitude,
                "amplitude_rms": math.sqrt(
                    sum((a - mean_amplitude) ** 2 for a in amplitudes)
                    / len(amplitudes)
                ),
            })
        series.append(entry)
    return series


def series_reading(series, amplitude):
    """Which antennas are offset, which wander, in sentences."""

    def worst(key):
        values = [
            pol[key] for entry in series for pol in entry["polarizations"]
            if key in pol
        ]
        return max((abs(v) for v in values), default=None)

    solved = [
        entry for entry in series
        if any("mean_phase_deg" in pol for pol in entry["polarizations"])
    ]
    if not solved:
        return ["No antenna has a solution: nothing was applied."]

    def named(key, limit):
        return [
            entry["antenna"] for entry in solved
            if any(abs(pol.get(key, 0.0)) > limit
                   for pol in entry["polarizations"])
        ]

    wandering = named("phase_rms_deg", 10.0)
    offset = named("mean_phase_deg", 5.0)
    lines = [
        f"{len(solved)} of {len(series)} antennas solved, "
        f"{series[0]['solutions']} solutions each. Mean phases reach "
        f"{worst('mean_phase_deg'):.1f} deg; "
        + (f"offset by more than 5 deg: {', '.join(offset)}."
           if offset else "none is more than 5 deg from the reference.")
    ]
    if series[0]["solutions"] > 1:
        lines.append(
            f"Phases move about their means by up to "
            f"{worst('phase_rms_deg'):.1f} deg rms; "
            + (f"more than 10 deg: {', '.join(wandering)}. One correction "
               "per scan cannot follow those antennas."
               if wandering else "no antenna wanders by more than 10 deg.")
        )
    if amplitude:
        low = [
            entry["antenna"] for entry in solved
            if any(abs(pol.get("mean_amplitude", 1.0) - 1.0) > 0.10
                   for pol in entry["polarizations"])
        ]
        lines.append(
            "Amplitudes more than 10% from the median: "
            + (", ".join(low) if low else "none") + "."
        )
    return lines


def delay_table(antennas, values, flags, names):
    """Per-antenna delays in nanoseconds from a K table's columns."""

    rows = {}
    for row, antenna in enumerate(antennas):
        antenna = int(antenna)
        if antenna in rows:
            continue
        rows[antenna] = [
            None if bool(flags[pol][0][row]) else float(values[pol][0][row])
            for pol in range(len(values))
        ]
    return [
        {
            "antenna": names[antenna] if antenna < len(names) else str(antenna),
            "id": antenna,
            "delay_ns": rows[antenna],
        }
        for antenna in sorted(rows)
    ]


def delay_reading(delays, bandwidth_mhz):
    """What the delays amount to across the band, in sentences.

    A delay of t nanoseconds turns the phase by 0.36 t degrees per MHz, so
    across a band of B MHz the phase runs through 0.36 t B degrees, half of
    it either side of the centre.
    """

    solved = [
        (entry["antenna"], delay)
        for entry in delays for delay in entry["delay_ns"][:1]
        if delay is not None
    ]
    if not solved:
        return ["No antenna has a delay solution."]
    swing = {name: 0.36 * delay * bandwidth_mhz for name, delay in solved}
    worst = max(swing, key=lambda name: abs(swing[name]))
    large = sorted(
        (name for name in swing if abs(swing[name]) > 20.0),
        key=lambda name: -abs(swing[name]),
    )
    lines = [
        f"{len(solved)} antennas have a delay. The largest is "
        f"{dict(solved)[worst]:+.2f} ns on {worst}: its phase runs through "
        f"{abs(swing[worst]):.0f} deg across {bandwidth_mhz:g} MHz."
    ]
    if large:
        lines.append(
            "Phase running through more than 20 deg across the band: "
            + ", ".join(f"{name} ({swing[name]:+.0f})" for name in large)
            + ". One phase per antenna leaves half of that at each end."
        )
    else:
        lines.append(
            "No antenna's phase runs through more than 20 deg across the "
            "band: delays are not what is left."
        )
    return lines


def phase_summary(table):
    """How large the offsets are: rms and extremes over the solved antennas."""

    values = [
        phase for entry in table for phase in entry["phase_deg"]
        if phase is not None
    ]
    solved = sum(
        1 for entry in table
        if any(phase is not None for phase in entry["phase_deg"])
    )
    if not values:
        return {"solved_antennas": 0, "antennas": len(table)}
    mean = sum(values) / len(values)
    rms = math.sqrt(sum((v - mean) ** 2 for v in values) / len(values))
    return {
        "antennas": len(table),
        "solved_antennas": solved,
        "rms_deg": rms,
        "largest_deg": max(abs(v) for v in values),
        "beyond_5_deg": sum(
            1 for entry in table
            if any(p is not None and abs(p) > 5.0 for p in entry["phase_deg"])
        ),
    }


def reading(summary):
    if not summary.get("solved_antennas"):
        return ["No antenna has a solution: nothing was applied."]
    lines = [
        f"{summary['solved_antennas']} of {summary['antennas']} antennas "
        f"solved. Phases have an rms of {summary['rms_deg']:.1f} deg and "
        f"reach {summary['largest_deg']:.1f} deg; "
        f"{summary['beyond_5_deg']} antennas are more than 5 deg from the "
        "reference."
    ]
    if summary["rms_deg"] >= 3.0:
        lines.append(
            "Those are fixed offsets that the earlier self-calibration left "
            "in the data. Imaging the corrected Measurement Set shows "
            "whether the pattern around the sources goes with them."
        )
    else:
        lines.append(
            "That is small: fixed phase offsets of this size do not explain "
            "a pattern at a tenth of each source's peak."
        )
    return lines


def wild_delays(antennas, values, flags, names, limit):
    """Delays too large to be a residual of calibrated data.

    The delays this script solves are what an earlier calibration left: on
    G0212-0335 the largest was 2.3 ns. A delay is found from the baseline
    to the reference antenna alone, and where the model has little flux on
    that baseline the search settles on whatever ripple the data have:
    twice it gave an antenna of the central square 17.9 ns, which turns
    its phase five times round across the band, and the amplitudes solved
    afterwards went to 0.27 trying to follow.
    """

    found = []
    for row, antenna in enumerate(antennas):
        antenna = int(antenna)
        large = [
            float(values[pol][0][row]) for pol in range(len(values))
            if not bool(flags[pol][0][row])
            and abs(float(values[pol][0][row])) > limit
        ]
        if large:
            found.append({
                "antenna": names[antenna] if antenna < len(names)
                else str(antenna),
                "delay_ns": large,
            })
    return found


def process_alive(pid):
    try:
        os.kill(int(pid), 0)
    except (ProcessLookupError, ValueError, TypeError):
        return False
    except PermissionError:
        return True
    return True


@contextlib.contextmanager
def claimed(path):
    """Hold one solve's products; refuse them while another run has them.

    Two solves given the same Measurement Set and suffix write the same
    calibration table. On 2026-10-09 a command started twice, fourteen
    seconds apart, had each remove the table the other was writing, and
    both stopped. The claim is a file beside the products, as in
    ``experiment_wsclean_image.claimed``: one left by a run that was killed
    is taken over when its process is gone from this machine.
    """

    path = Path(path)
    mine = {
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "since": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
    }
    for attempt in range(2):
        try:
            handle = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            try:
                other = json.loads(path.read_text())
            except (OSError, ValueError):
                other = {}
            gone = (
                other.get("host") == mine["host"]
                and not process_alive(other.get("pid"))
            )
            if gone and attempt == 0:
                path.unlink(missing_ok=True)
                continue
            raise SystemExit(
                f"Another run of this solve is going (process "
                f"{other.get('pid', '?')} on {other.get('host', '?')}, since "
                f"{other.get('since', '?')}). Two would remove each other's "
                "calibration table. Wait for it, or give --suffix another "
                f"name; if that run is gone, remove {path}."
            ) from None
        with os.fdopen(handle, "w") as lock:
            json.dump(mine, lock)
        break
    try:
        yield path
    finally:
        path.unlink(missing_ok=True)


def output_path(ms, suffix):
    ms = Path(str(ms).rstrip("/"))
    return ms.with_name(ms.stem + suffix + ms.suffix)


def run(args):
    from casatasks import applycal, gaincal, split
    from casatools import table

    ms = Path(str(args.ms).rstrip("/"))
    if not ms.is_dir():
        raise SystemExit(f"No Measurement Set at {ms}.")
    out_ms = output_path(ms, args.suffix)
    caltable = ms.with_name(ms.stem + args.suffix + ".phase.cal")
    delaytable = ms.with_name(ms.stem + args.suffix + ".delay.cal")
    report_path = ms.with_name(ms.stem + args.suffix + ".json")

    tool = table()
    tool.open(str(ms))
    try:
        columns = set(tool.colnames())
    finally:
        tool.close()
    if "MODEL_DATA" not in columns:
        raise SystemExit(
            f"{ms} has no MODEL_DATA column. Run experiment_wsclean_image.py "
            "--core-model on it first, from the terminal that has WSClean."
        )

    existing = [
        path for path in (out_ms, caltable, delaytable) if path.exists()
    ]
    if existing and not args.overwrite:
        raise SystemExit(
            f"{existing[0]} already exists. Pass --overwrite to redo the "
            "solve, or --suffix for a separate result."
        )
    for path in existing:
        shutil.rmtree(path)

    def antenna_names():
        tool.open(str(ms / "ANTENNA"))
        try:
            return [str(name) for name in tool.getcol("NAME")]
        finally:
            tool.close()

    delays = None
    rejected_delays = []
    prior = []
    if args.delays:
        print(
            f"Solving a delay per antenna against MODEL_DATA of {ms} ...",
            flush=True,
        )
        gaincal(
            vis=str(ms),
            caltable=str(delaytable),
            gaintype="K",
            solint="inf",
            combine="",
            refant=args.refant,
            minsnr=float(args.minsnr),
        )
        if not delaytable.is_dir():
            raise SystemExit(
                "gaincal wrote no delay table: no delay could be solved."
            )
        names_k = antenna_names()
        limit = float(args.max_delay_ns)
        tool.open(str(delaytable), nomodify=False)
        try:
            antennas_k = tool.getcol("ANTENNA1")
            values_k = tool.getcol("FPARAM")
            flags_k = tool.getcol("FLAG")
            rejected_delays = wild_delays(
                antennas_k, values_k, flags_k, names_k, limit
            )
            if rejected_delays:
                # Not believed, so not applied: the antenna keeps the delay
                # it has and still gets its phase.
                values_k = values_k * (abs(values_k) <= limit)
                tool.putcol("FPARAM", values_k)
            delays = delay_table(antennas_k, values_k, flags_k, names_k)
        finally:
            tool.close()
        if rejected_delays:
            print(
                f"Delays beyond {limit:g} ns are not believed and are not "
                "applied: " + ", ".join(
                    f"{entry['antenna']} ("
                    + ", ".join(f"{d:+.2f}" for d in entry["delay_ns"])
                    + " ns)"
                    for entry in rejected_delays
                ) + "."
            )
        print("\nDelay per antenna and polarization (ns):")
        for entry in delays:
            print(f"  {entry['antenna']:<5} " + "  ".join(
                "flagged" if delay is None else f"{delay:+7.3f}"
                for delay in entry["delay_ns"]
            ))
        for line in delay_reading(delays, args.bandwidth_mhz):
            print(line)
        print()
        prior = [str(delaytable)]

    amplitude = args.calmode == "ap"
    print(
        f"Solving {'amplitude and phase' if amplitude else 'phase'} per "
        f"antenna ({args.gaintype}, solint {args.solint}) against "
        f"MODEL_DATA of {ms} ...",
        flush=True,
    )
    solve = dict(
        vis=str(ms),
        caltable=str(caltable),
        calmode=args.calmode,
        solint=args.solint,
        combine="",
        gaintype=args.gaintype,
        refant=args.refant,
        refantmode="flex",
        minsnr=float(args.minsnr),
        minblperant=4,
    )
    if amplitude:
        # The flux scale stays the calibrators': amplitudes are only
        # redistributed between antennas.
        solve.update(solnorm=True, normtype="median")
    if args.uvrange:
        # A cores-only model lacks the extended emission the short
        # baselines see; an amplitude solved on them would be solved to
        # hide it. The delay solve keeps every baseline: it works on those
        # to the reference antenna, which for the central square are short.
        solve.update(uvrange=args.uvrange)
    if prior:
        # The gains are what is left once the delays are taken out.
        solve.update(gaintable=prior, interp=["nearest"])
    gaincal(**solve)
    if not caltable.is_dir():
        raise SystemExit("gaincal wrote no table: no solution was possible.")

    names = antenna_names()
    rejected_amplitudes = []
    rejected_scans = []
    trend = None
    positions = None
    medians = None
    if amplitude and getattr(args, "detrend_radius", None):
        tool.open(str(ms / "ANTENNA"))
        try:
            positions = tool.getcol("POSITION").T
        finally:
            tool.close()
    tool.open(str(caltable), nomodify=False)
    try:
        values = tool.getcol("CPARAM")
        flags = tool.getcol("FLAG")
        antennas = tool.getcol("ANTENNA1")
        times = tool.getcol("TIME")
        if positions is not None:
            trend = radial_trend(antennas, values, flags, positions, names)
            if trend and trend["ratio"] and (
                abs(trend["ratio"] - 1.0) > float(args.detrend_radius)
            ):
                values = detrend_amplitudes(antennas, values, trend)
                tool.putcol("CPARAM", values)
                trend["divided_out"] = True
                print(
                    "Amplitudes run with distance from the array's middle "
                    f"(inner half {trend['inner_median']:.2f}, outer half "
                    f"{trend['outer_median']:.2f}): that is the model's, so "
                    "the run is divided out and only each antenna's own "
                    "departure from it is applied."
                )
            elif trend:
                trend["divided_out"] = False
        if amplitude:
            # CASA's normalisation is undone by the line divided out, and
            # is set again here for every amplitude solve.
            values, medians = normalise_amplitudes(values, flags)
            if any(abs(median - 1.0) > 1e-4 for median in medians):
                tool.putcol("CPARAM", values)
            # An amplitude this far from one is not applied: either the
            # solve went wrong or the antenna's data are too far out to be
            # mended by a factor.
            low, high = (float(v) for v in args.amplitude_limits)
            size = abs(values)
            wild = (flags == 0) & ((size < low) | (size > high))
            if wild.any():
                scans = tool.getcol("SCAN_NUMBER")
                flags = (flags != 0) | wild
                tool.putcol("FLAG", flags)

                def name_of(antenna):
                    return (names[int(antenna)] if int(antenna) < len(names)
                            else str(int(antenna)))

                rejected_amplitudes = sorted({
                    name_of(antenna)
                    for row, antenna in enumerate(antennas)
                    if wild[:, :, row].any()
                })
                rejected_scans = sorted({
                    (name_of(antenna), int(scans[row]))
                    for row, antenna in enumerate(antennas)
                    if wild[:, :, row].any()
                })
                print(
                    f"Amplitudes outside {low:g} to {high:g} are not "
                    "believed and are not applied: "
                    + ", ".join(rejected_amplitudes) + "."
                )
    finally:
        tool.close()

    phases = phase_table(antennas, values, flags, names)
    summary = phase_summary(phases)
    series = gain_series(antennas, times, values, flags, names)
    several = bool(series) and series[0]["solutions"] > 1
    per_scan = str(args.solint).strip().lower() in ("inf", "infinite")
    report = {
        "measurement_set": str(ms),
        "caltable": str(caltable),
        "delays": delays,
        "rejected_delays": rejected_delays,
        "rejected_amplitudes": rejected_amplitudes,
        "radial_trend": (
            {key: value for key, value in trend.items() if key != "line"}
            if trend else None
        ),
        "amplitude_medians": medians,
        "delay_reading": (
            delay_reading(delays, args.bandwidth_mhz) if delays else None
        ),
        "solve": {
            "calmode": args.calmode, "uvrange": args.uvrange,
            "solint": args.solint, "gaintype": args.gaintype,
            "refant": args.refant, "minsnr": float(args.minsnr),
        },
        "phases": phases,
        "summary": summary,
        "series": series,
        "reading": (
            series_reading(series, amplitude) if several or amplitude
            else reading(summary)
        ),
    }

    print(
        "\nPer antenna and polarization: mean phase (deg), phase rms "
        "about it (deg)" + (", mean amplitude" if amplitude else "") + ":"
    )
    for entry in series:
        parts = []
        for pol in entry["polarizations"]:
            if "mean_phase_deg" not in pol:
                parts.append("      failed      ")
                continue
            part = f"{pol['mean_phase_deg']:+7.1f} {pol['phase_rms_deg']:5.1f}"
            if amplitude:
                part += f" {pol['mean_amplitude']:5.2f}"
            parts.append(part)
        print(f"  {entry['antenna']:<5} " + "   ".join(parts))

    if summary.get("solved_antennas"):
        print("\nApplying and writing the corrected data ...", flush=True)
        applycal(
            vis=str(ms),
            gaintable=prior + [str(caltable)],
            # One solution for a scan is that scan's own, also where there
            # are several scans: an offset fixed within one need not be the
            # next one's, and a line drawn between them would be neither.
            # Several within a scan are joined linearly, as the pipeline
            # joins its own.
            interp=["nearest"] * len(prior)
            + ["linear" if several and not per_scan else "nearest"],
            applymode="calonly",
            # An amplitude correction scales an antenna's noise with its
            # signal, so its weights have to follow. Left as they were, an
            # antenna turned up by 1/0.7 keeps the weight of a quiet one:
            # on G0212-0335 that raised the noise of the data from 23.6 to
            # 26.5 uJy/beam.
            calwt=[False] * len(prior) + [amplitude],
        )
        split(vis=str(ms), outputvis=str(out_ms), datacolumn="corrected")
        report["corrected_measurement_set"] = str(out_ms)
        # Data whose amplitude was not believed are not left in at that
        # amplitude. On G0216-0604 one antenna of the south arm, S04, was
        # beyond the limit in one of two scans; left as it was, it put a
        # hole of 9% of the peak beside every source (finding 134). The
        # antenna is flagged for that scan, unless the solve is towards
        # one source and the data are the whole field's.
        flagged = []
        if rejected_scans and per_scan and not args.keep_unsolved:
            from casatasks import flagdata

            for name, scan in rejected_scans:
                if scan < 0:
                    continue
                flagdata(vis=str(out_ms), mode="manual", antenna=name,
                         scan=str(scan), flagbackup=False)
                flagged.append({"antenna": name, "scan": scan})
            if flagged:
                print(
                    "Flagged in the corrected data: "
                    + ", ".join(f"{entry['antenna']} in scan {entry['scan']}"
                                for entry in flagged) + "."
                )
        report["flagged_rejected"] = flagged
        report["reading"].append(
            f"Corrected data: {out_ms}. Image it with "
            "experiment_wsclean_image.py --ms, from the terminal that has "
            "WSClean."
        )

    with open(report_path, "w") as handle:
        json.dump(report, handle, indent=2)
    print("\n" + "\n".join(report["reading"]))
    print(f"\nReport written: {report_path}")
    return report


# ---------------------------------------------------------------------------
# Gains against frequency
# ---------------------------------------------------------------------------

def bandpass_table(antennas, values, flags, names):
    """Per antenna and polarization, the phase and amplitude of each bin.

    ``values`` and ``flags`` are indexed [polarization, bin, row]; the first
    row of an antenna is taken, which is the only one for a single scan. A
    flagged bin is None.
    """

    rows = {}
    for row, antenna in enumerate(antennas):
        rows.setdefault(int(antenna), row)
    table = []
    for antenna in sorted(rows):
        row = rows[antenna]
        entry = {
            "antenna": names[antenna] if antenna < len(names) else str(antenna),
            "id": antenna,
            "polarizations": [],
        }
        for pol in range(len(values)):
            phases, amplitudes = [], []
            for chan in range(len(values[pol])):
                value = complex(values[pol][chan][row])
                if bool(flags[pol][chan][row]) or value == 0:
                    phases.append(None)
                    amplitudes.append(None)
                else:
                    phases.append(
                        math.degrees(math.atan2(value.imag, value.real))
                    )
                    amplitudes.append(abs(value))
            entry["polarizations"].append(
                {"phase_deg": phases, "amplitude": amplitudes}
            )
        table.append(entry)
    return table


def bandpass_statistics(table):
    """How far each antenna's gain moves across the band.

    For each antenna and polarization: the rms of the phase about its mean,
    its span, what is left of it (rms) once the best straight line through
    the bins is taken out, and the rms of the amplitude as a fraction of
    its mean. A phase that is a straight line is a delay that an earlier
    step could have taken; what is left after the line is what no delay
    can take.
    """

    statistics = []
    for entry in table:
        row = {"antenna": entry["antenna"], "polarizations": []}
        for pol in entry["polarizations"]:
            good = [
                (index, phase, amplitude)
                for index, (phase, amplitude) in enumerate(
                    zip(pol["phase_deg"], pol["amplitude"])
                )
                if phase is not None
            ]
            if len(good) < 3:
                row["polarizations"].append({"bins": len(good)})
                continue
            total = sum(
                complex(math.cos(math.radians(p)), math.sin(math.radians(p)))
                for _, p, _ in good
            )
            mean = math.degrees(math.atan2(total.imag, total.real))
            x = [index for index, _, _ in good]
            y = [((p - mean + 180.0) % 360.0) - 180.0 for _, p, _ in good]
            count = len(good)
            x_mean, y_mean = sum(x) / count, sum(y) / count
            spread = sum((a - x_mean) ** 2 for a in x)
            slope = (
                sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y)) / spread
                if spread else 0.0
            )
            bowed = [
                b - y_mean - slope * (a - x_mean) for a, b in zip(x, y)
            ]
            amplitudes = [a for _, _, a in good]
            amplitude_mean = sum(amplitudes) / count
            row["polarizations"].append({
                "bins": count,
                "phase_mean_deg": mean + y_mean,
                "phase_rms_deg": math.sqrt(
                    sum((b - y_mean) ** 2 for b in y) / count
                ),
                "phase_span_deg": max(y) - min(y),
                "slope_deg_across_band": slope * (max(x) - min(x)),
                "phase_rms_after_line_deg": math.sqrt(
                    sum(b * b for b in bowed) / count
                ),
                "amplitude_mean": amplitude_mean,
                "amplitude_rms": math.sqrt(
                    sum((a - amplitude_mean) ** 2 for a in amplitudes) / count
                ) / amplitude_mean if amplitude_mean else None,
            })
        statistics.append(row)
    return statistics


def common_amplitude(table):
    """The amplitude all antennas share in each bin: their median.

    Returned per polarization as a list over the bins, None where no
    antenna has a solution. An antenna's own bandpass is its amplitude
    over this. What every antenna shows alike is not an antenna: it is the
    model. WSClean's model is one image for each of its sub-bands, a
    staircase in frequency, and a sky that falls smoothly across a
    sub-band lies above the step at its low end and below it at the high
    end, in every antenna's solution.
    """

    if not table:
        return []
    common = []
    for pol in range(len(table[0]["polarizations"])):
        bins = len(table[0]["polarizations"][pol]["amplitude"])
        row = []
        for chan in range(bins):
            values = sorted(
                entry["polarizations"][pol]["amplitude"][chan]
                for entry in table
                if entry["polarizations"][pol]["amplitude"][chan] is not None
            )
            if not values:
                row.append(None)
                continue
            middle = len(values) // 2
            row.append(
                values[middle] if len(values) % 2
                else 0.5 * (values[middle - 1] + values[middle])
            )
        common.append(row)
    return common


def common_reading(common):
    """How large the shared ripple is, in a sentence."""

    values = [v for row in common for v in row if v is not None]
    if len(values) < 2:
        return []
    mean = sum(values) / len(values)
    rms = math.sqrt(sum((v - mean) ** 2 for v in values) / len(values))
    return [
        f"The amplitude common to all antennas moves by {100 * rms:.1f}% rms "
        f"across the band ({min(values):.3f} to {max(values):.3f}). That is "
        "the model's steps in frequency, not the antennas, and it is never "
        "applied."
    ]


def bandpass_reading(statistics, channels):
    """Which antennas' phases bend across the band, in sentences."""

    solved = [
        (row["antenna"], pol)
        for row in statistics for pol in row["polarizations"]
        if "phase_rms_deg" in pol
    ]
    if not solved:
        return ["No antenna has a solution: nothing was applied."]
    antennas = sorted({name for name, _ in solved})

    def worst(key):
        return max(solved, key=lambda item: abs(item[1][key]))

    def named(key, limit):
        return sorted({
            name for name, pol in solved if abs(pol[key]) > limit
        })

    name, pol = worst("phase_rms_deg")
    lines = [
        f"{len(antennas)} antennas solved, a gain for every {channels} "
        f"channels. Across the band the phases move by up to "
        f"{pol['phase_rms_deg']:.1f} deg rms ({name}), "
        f"{pol['phase_span_deg']:.0f} deg from lowest to highest."
    ]
    moving = named("phase_rms_deg", 3.0)
    lines.append(
        "More than 3 deg rms: " + (", ".join(moving) if moving else "none")
        + "."
    )
    name, pol = worst("phase_rms_after_line_deg")
    bowed = named("phase_rms_after_line_deg", 3.0)
    lines.append(
        f"With the best straight line taken out, which a delay could have "
        f"been, up to {pol['phase_rms_after_line_deg']:.1f} deg rms is left "
        f"({name}); more than 3 deg: "
        + (", ".join(bowed) if bowed else "none") + "."
    )
    name, pol = worst("amplitude_rms")
    uneven = named("amplitude_rms", 0.05)
    lines.append(
        f"Amplitudes vary across the band by up to "
        f"{100 * pol['amplitude_rms']:.1f}% rms ({name}); more than 5%: "
        + (", ".join(uneven) if uneven else "none") + "."
    )
    return lines


def run_bandpass(args):
    """Solve a gain per channel bin against MODEL_DATA, apply it, split."""

    from casatasks import applycal, bandpass, split
    from casatools import table

    ms = Path(str(args.ms).rstrip("/"))
    if not ms.is_dir():
        raise SystemExit(f"No Measurement Set at {ms}.")
    out_ms = output_path(ms, args.suffix)
    caltable = ms.with_name(ms.stem + args.suffix + ".bandpass.cal")
    report_path = ms.with_name(ms.stem + args.suffix + ".json")
    prior = [str(Path(path)) for path in (args.prior or [])]
    for path in prior:
        if not Path(path).is_dir():
            raise SystemExit(f"No calibration table at {path}.")

    tool = table()
    tool.open(str(ms))
    try:
        columns = set(tool.colnames())
    finally:
        tool.close()
    if "MODEL_DATA" not in columns:
        raise SystemExit(
            f"{ms} has no MODEL_DATA column. Write a model into it first "
            "with experiment_wsclean_image.py."
        )
    existing = [path for path in (out_ms, caltable) if path.exists()]
    if existing and not args.overwrite:
        raise SystemExit(
            f"{existing[0]} already exists. Pass --overwrite to redo the "
            "solve, or --suffix for a separate result."
        )
    for path in existing:
        shutil.rmtree(path)

    channels = int(args.bandpass)
    # Flagging what could not be solved is for a solve on the whole field,
    # where an antenna without a solution has no usable data there. Towards
    # one faint source a bin without a solution only says the source was
    # too faint in it, and the data, which are the field's, are kept: five
    # sources solved one after another had otherwise flagged the top two
    # sub-bands of G0212-0335 away between them.
    flag_unsolved = bool(args.bandpass_amplitude) and not args.keep_unsolved
    print(
        f"Solving a gain per antenna for every {channels} channels against "
        f"MODEL_DATA of {ms} ...",
        flush=True,
    )
    solve = dict(
        vis=str(ms),
        caltable=str(caltable),
        solint=f"inf,{channels}ch",
        combine="",
        bandtype="B",
        refant=args.refant,
        minsnr=float(args.minsnr),
        minblperant=4,
        solnorm=False,
        fillgaps=0,
    )
    if prior:
        # What an earlier solve found is taken out first, as it will be
        # when both are applied.
        solve.update(gaintable=prior, interp=["linear"] * len(prior))
    bandpass(**solve)
    if not caltable.is_dir():
        raise SystemExit("bandpass wrote no table: no solution was possible.")

    tool.open(str(ms / "ANTENNA"))
    try:
        names = [str(name) for name in tool.getcol("NAME")]
        positions = (
            tool.getcol("POSITION").T
            if args.bandpass_amplitude and getattr(args, "detrend_radius", None)
            else None
        )
    finally:
        tool.close()
    rejected_bins = {}
    trend = None
    medians = None
    tool.open(str(caltable), nomodify=False)
    try:
        values = tool.getcol("CPARAM")
        flags = tool.getcol("FLAG")
        antennas = tool.getcol("ANTENNA1")
        gains = bandpass_table(antennas, values, flags, names)
        common = common_amplitude(gains)
        shared = values.real * 0 + 1
        for pol, row in enumerate(common):
            for chan, value in enumerate(row):
                if value:
                    shared[pol, chan, :] = value
        # A bin whose amplitude is far from what the antennas share is a
        # solve that went wrong: it is flagged in the table, so that it is
        # neither applied nor counted.
        low, high = (float(v) for v in args.bin_limits)
        ratio = abs(values) / shared
        wild = (flags == 0) & ((ratio < low) | (ratio > high))
        if wild.any():
            flags = (flags != 0) | wild
            tool.putcol("FLAG", flags)
            for row, antenna in enumerate(antennas):
                count = int(wild[:, :, row].sum())
                if count:
                    name = (
                        names[int(antenna)] if int(antenna) < len(names)
                        else str(int(antenna))
                    )
                    rejected_bins[name] = rejected_bins.get(name, 0) + count
            gains = bandpass_table(antennas, values, flags, names)
        if args.bandpass_amplitude:
            # Each antenna's amplitude over what all antennas share: its
            # own bandpass, without the model's steps in frequency.
            own = values / shared
            if positions is not None:
                # What still runs with distance from the array's middle is
                # the model's, as in the gain solves (finding 148).
                trend = radial_trend(antennas, own, flags, positions, names)
                if trend and trend["ratio"] and (
                    abs(trend["ratio"] - 1.0) > float(args.detrend_radius)
                ):
                    own = detrend_amplitudes(antennas, own, trend)
                    own, medians = normalise_amplitudes(own, flags)
                    trend["divided_out"] = True
                elif trend:
                    trend["divided_out"] = False
            tool.putcol("CPARAM", own)
        else:
            # The phases are what a delay could not follow; every gain is
            # set to unit amplitude.
            size = abs(values)
            tool.putcol("CPARAM", values / (size + (size == 0)))
    finally:
        tool.close()

    statistics = bandpass_statistics(gains)
    solved = any(
        "phase_rms_deg" in pol
        for row in statistics for pol in row["polarizations"]
    )
    report = {
        "measurement_set": str(ms),
        "caltable": str(caltable),
        "prior": prior,
        "solve": {
            "bandpass_channels": channels, "refant": args.refant,
            "minsnr": float(args.minsnr),
            "amplitudes_applied": bool(args.bandpass_amplitude),
        },
        "bandpass": gains,
        "rejected_bins": rejected_bins,
        "radial_trend": (
            {key: value for key, value in trend.items() if key != "line"}
            if trend else None
        ),
        "amplitude_medians": medians,
        "common_amplitude": common,
        "statistics": statistics,
        "reading": bandpass_reading(statistics, channels)
        + common_reading(common),
    }

    print(
        "\nPer antenna and polarization: phase rms across the band (deg), "
        "the same with a straight line taken out, amplitude rms (%):"
    )
    for row in statistics:
        parts = []
        for pol in row["polarizations"]:
            if "phase_rms_deg" not in pol:
                parts.append("      failed      ")
                continue
            parts.append(
                f"{pol['phase_rms_deg']:5.1f} "
                f"{pol['phase_rms_after_line_deg']:5.1f} "
                f"{100 * pol['amplitude_rms']:5.1f}"
            )
        print(f"  {row['antenna']:<5} " + "   ".join(parts))

    if solved:
        print("\nApplying and writing the corrected data ...", flush=True)
        applycal(
            vis=str(ms),
            gaintable=prior + [str(caltable)],
            # In time a scan's solution is that scan's; in frequency the
            # bins are joined. Where amplitudes are corrected, channels of
            # an antenna that a whole scan could not solve are flagged
            # rather than given a neighbour's gain: an antenna whose gain
            # doubles towards a band edge is not to be guessed beyond it.
            interp=["linear"] * len(prior) + [
                "nearest,linearflag" if flag_unsolved else "nearest,linear"
            ],
            applymode="calflag" if flag_unsolved else "calonly",
            # Weights follow every amplitude that is applied.
            calwt=[True] * len(prior) + [bool(args.bandpass_amplitude)],
        )
        split(vis=str(ms), outputvis=str(out_ms), datacolumn="corrected")
        report["corrected_measurement_set"] = str(out_ms)
        report["reading"].append(
            f"Corrected data ("
            + ("phases and amplitudes" if args.bandpass_amplitude
               else "phases only")
            + f" across the band): {out_ms}."
        )

    with open(report_path, "w") as handle:
        json.dump(report, handle, indent=2)
    print("\n" + "\n".join(report["reading"]))
    print(f"\nReport written: {report_path}")
    return report


def radial_trend(antennas, values, flags, positions, names=None):
    """How the amplitudes of a solve run with distance from the array's middle.

    A cores-only model of a resolved source holds the wrong share of its
    flux on short baselines against long ones, and an amplitude solve then
    turns the antennas up or down by where they stand: on G0212-0335 in Band
    5 (2" beam) the central square came out at 1.03-1.38 and the far arm
    antennas at 0.33-0.6, the target's core rose by 84% and a bowl of
    negatives and 40% of its flux were lost. A faulty antenna stands alone;
    the model's error follows the distance. Each antenna's median amplitude
    is set against the logarithm of its distance from the median position
    of the antennas solved, and a Theil-Sen line drawn through them.
    """

    import numpy as np

    antennas = np.asarray(antennas)
    size = np.abs(np.asarray(values))
    good = np.asarray(flags) == 0
    positions = np.asarray(positions, dtype=float)
    rows = []
    for antenna in sorted(set(int(value) for value in antennas)):
        select = antennas == antenna
        amplitudes = size[:, :, select][good[:, :, select]]
        if amplitudes.size and antenna < len(positions):
            rows.append([antenna, float(np.median(amplitudes))])
    if len(rows) < 6:
        return None
    ids = [row[0] for row in rows]
    centre = np.median(positions[ids], axis=0)
    distance = np.array([
        max(float(np.linalg.norm(positions[antenna] - centre)), 100.0)
        for antenna in ids
    ])
    amplitude = np.array([row[1] for row in rows])
    x = np.log10(distance)
    slopes = [
        (amplitude[j] - amplitude[i]) / (x[j] - x[i])
        for i in range(len(x)) for j in range(i + 1, len(x))
        if x[j] != x[i]
    ]
    slope = float(np.median(slopes)) if slopes else 0.0
    intercept = float(np.median(amplitude - slope * x))
    line = np.clip(intercept + slope * x, 0.05, None)
    inner = amplitude[distance <= np.median(distance)]
    outer = amplitude[distance > np.median(distance)]
    inner_median = float(np.median(inner))
    outer_median = float(np.median(outer)) if outer.size else inner_median
    return {
        "inner_median": round(inner_median, 3),
        "outer_median": round(outer_median, 3),
        "ratio": round(inner_median / outer_median, 3) if outer_median else None,
        "slope_per_decade": round(slope, 3),
        "line": {
            int(antenna): float(value) for antenna, value in zip(ids, line)
        },
        "antennas": [
            {
                "antenna": names[antenna] if names and antenna < len(names)
                else str(antenna),
                "distance_km": round(float(d) / 1e3, 3),
                "amplitude": round(float(a), 3),
                "line": round(float(t), 3),
            }
            for antenna, d, a, t in zip(ids, distance, amplitude, line)
        ],
    }


def detrend_amplitudes(antennas, values, trend):
    """The gains with the line of ``radial_trend`` divided out.

    What stays is each antenna's own departure from the line; the phases
    are untouched. The level is not kept: ``normalise_amplitudes`` sets it
    afterwards.
    """

    import numpy as np

    values = np.array(values, dtype=complex)
    line = trend["line"]
    for row, antenna in enumerate(np.asarray(antennas)):
        factor = line.get(int(antenna))
        if factor:
            values[:, :, row] /= factor
    return values


def normalise_amplitudes(values, flags):
    """The gains scaled so that each hand's median amplitude is one.

    The flux scale is the calibrators'; a gain solve may only move
    amplitude between antennas. Scaling each antenna to the median of the
    line drawn by ``radial_trend`` did not keep it: the antennas crowd at
    the array's middle, so the median of the line fell below the median
    antenna, at 0.85 on G0212-0335 in Band 5, and every visibility was
    turned up by 1.6 (2026-10-23). The median is taken per hand over the
    unflagged solutions of every antenna, time and channel. Returns the
    gains and the median each hand had.
    """

    import numpy as np

    values = np.array(values, dtype=complex)
    good = np.asarray(flags) == 0
    medians = []
    for pol in range(values.shape[0]):
        size = np.abs(values[pol])[good[pol]]
        size = size[np.isfinite(size) & (size > 0)]
        median = float(np.median(size)) if size.size else 1.0
        values[pol] /= median
        medians.append(round(median, 4))
    return values, medians


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--ms", required=True,
                        help="a copy of the data with the cores-only model "
                        "in MODEL_DATA")
    parser.add_argument("--refant", default="C00,W03,E05,E02,E06",
                        help="reference antennas in order of preference")
    parser.add_argument("--solint", default="inf",
                        help="solution interval (default inf: one per scan)")
    parser.add_argument("--uvrange", default="",
                        help="baselines used for the gain solve, in CASA's "
                        "syntax, for example '>3klambda' (default: all)")
    parser.add_argument("--delays", action="store_true",
                        help="solve a delay per antenna first, and the "
                        "gains on top of it")
    parser.add_argument("--bandwidth-mhz", type=float, default=280.0,
                        help="usable bandwidth, for stating what a delay "
                        "does across it (default 280)")
    parser.add_argument("--calmode", default="p", choices=("p", "ap"),
                        help="p: phase only; ap: amplitude and phase, "
                        "normalised so the flux scale is kept")
    parser.add_argument("--gaintype", default="G", choices=("G", "T"),
                        help="G: a phase per polarization; T: one for both")
    parser.add_argument("--minsnr", type=float, default=3.0)
    parser.add_argument("--suffix", default="_ph",
                        help="appended to the name of the corrected "
                        "Measurement Set (default _ph)")
    parser.add_argument("--bandpass", type=int, default=0, metavar="N",
                        help="instead of one gain for the band, solve one "
                        "for every N channels and correct the phases with "
                        "it")
    parser.add_argument("--bandpass-amplitude", action="store_true",
                        help="with --bandpass, correct each antenna's "
                        "amplitude across the band as well, relative to "
                        "what all antennas share, and flag what could not "
                        "be solved")
    parser.add_argument("--max-delay-ns", type=float, default=5.0,
                        help="a delay larger than this is not believed and "
                        "not applied (default 5)")
    parser.add_argument("--amplitude-limits", type=float, nargs=2,
                        default=[0.25, 4.0], metavar=("LOW", "HIGH"),
                        help="an amplitude outside these is not believed "
                        "and not applied, and with one solution per scan "
                        "the antenna is flagged for that scan (default "
                        "0.25 4: on this observation real ones reach 2.5 "
                        "and more)")
    parser.add_argument("--bin-limits", type=float, nargs=2,
                        default=[0.4, 2.5], metavar=("LOW", "HIGH"),
                        help="with --bandpass-amplitude, a channel bin "
                        "whose amplitude is outside these, relative to "
                        "what the antennas share, is not believed and not "
                        "applied (default 0.4 2.5)")
    parser.add_argument("--detrend-radius", type=float, metavar="FRACTION",
                        help="with --calmode ap: when the amplitudes of the "
                        "inner and outer halves of the array differ by more "
                        "than this (e.g. 0.1), divide out their run with "
                        "distance from the array's middle, which is the "
                        "model's and not the antennas'")
    parser.add_argument("--keep-unsolved", action="store_true",
                        help="flag nothing: with --bandpass-amplitude, not "
                        "the channels of an antenna that could not be "
                        "solved; with --calmode ap, not an antenna whose "
                        "amplitude was not believed")
    parser.add_argument("--prior", nargs="*", default=[],
                        metavar="CALTABLE",
                        help="with --bandpass, tables of earlier solves on "
                        "this Measurement Set to apply first and together "
                        "with the result")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    if not args.bandpass and (args.prior or args.bandpass_amplitude):
        parser.error("--prior and --bandpass-amplitude go with --bandpass")
    ms = Path(str(args.ms).rstrip("/"))
    if not ms.is_dir():
        raise SystemExit(f"No Measurement Set at {ms}.")
    with claimed(ms.with_name(ms.stem + args.suffix + ".lock")):
        if args.bandpass:
            run_bandpass(args)
        else:
            run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
