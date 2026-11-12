#!/usr/bin/env python3
"""How far each antenna's delay moves during an observation.

The pipeline solves one delay per antenna on the bandpass calibrator, at
one time, and applies it to the whole observation. Where an antenna's
delay drifts, later scans keep a phase slope across the band that one
phase per scan cannot follow, and averaging the band then cancels the
signal. This solves a delay on every scan of the calibrators named, with
the pipeline's delay and bandpass tables applied, and prints what is left:
for each antenna and polarisation, the delay (ns) at every scan and how
far it moves.

Run it with CASA's Python, from the pipeline folder, for example:

    singularity exec IMAGE.sif python3 diagnose_delay_drift.py \\
        --ms Band5_Products/ingest/Band5.ms \\
        --caldir Band5_Products/calibration --fields 0323+055,3C48

The Measurement Set is only read. The delay table it solves goes to
``--out`` (default ``delay_drift/<ms name>.K.cal`` in the current folder),
and a report to the same name with ``.json``.
"""

import argparse
import json
import math
import os
import shutil
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
# Run as `python3 diagnose_delay_drift.py` the pipeline folder comes first
# on sys.path, where its own inspect.py stands in for the standard module.
if sys.path and sys.path[0] and Path(sys.path[0]).resolve() == SCRIPT_DIR:
    sys.path.pop(0)


def excluded_antennas(caldir):
    """The antennas the pipeline excluded from its calibration."""

    path = Path(caldir) / "delay_quarantined_antennas.json"
    try:
        record = json.loads(path.read_text())
    except (OSError, ValueError):
        return []
    return [str(item["name"]) for item in record.get("antennas", [])]


def drift_table(names, rows):
    """Delays per antenna, polarisation and scan, and how far they move.

    ``rows`` holds (antenna id, scan, field, [delay per polarisation],
    [flag per polarisation]). Returns the scans in order and, per antenna,
    the delays and the spread (largest minus smallest) per polarisation.
    """

    scans = sorted({row[1] for row in rows})
    antennas = {}
    for antenna, scan, field, delays, flags in rows:
        entry = antennas.setdefault(antenna, {
            "antenna": names[antenna], "id": int(antenna), "by_scan": {},
        })
        entry["by_scan"][int(scan)] = {
            "field": field,
            "delay_ns": [
                None if flag or not math.isfinite(value) else round(float(value), 3)
                for value, flag in zip(delays, flags)
            ],
        }
    for entry in antennas.values():
        spreads = []
        for pol in range(2):
            values = [
                item["delay_ns"][pol] for item in entry["by_scan"].values()
                if len(item["delay_ns"]) > pol and item["delay_ns"][pol] is not None
            ]
            spreads.append(round(max(values) - min(values), 3) if len(values) > 1 else None)
        entry["spread_ns"] = spreads
    ordered = [antennas[key] for key in sorted(antennas)]
    return scans, ordered


def format_table(scans, antennas, fields):
    lines = []
    header = "ant    " + " ".join(f"{scan:>13}" for scan in scans) + "   spread R/L"
    lines.append(header)
    lines.append("field  " + " ".join(f"{fields.get(scan, '')[:13]:>13}" for scan in scans))
    for entry in antennas:
        cells = []
        for scan in scans:
            item = entry["by_scan"].get(scan)
            if not item:
                cells.append(" " * 13)
                continue
            text = "/".join(
                "  -  " if value is None else f"{value:+5.1f}"
                for value in item["delay_ns"][:2]
            )
            cells.append(f"{text:>13}")
        spread = "/".join(
            " - " if value is None else f"{value:.1f}" for value in entry["spread_ns"]
        )
        lines.append(f"{entry['antenna']:6} " + " ".join(cells) + f"   {spread}")
    return "\n".join(lines)


def run(args):
    from casatasks import gaincal
    from casatools import table

    ms = Path(args.ms).resolve()
    caldir = Path(args.caldir).resolve()
    out = Path(args.out or Path("delay_drift") / f"{ms.stem}.K.cal").resolve()
    for forbidden in (ms,):
        try:
            out.relative_to(forbidden)
        except ValueError:
            continue
        raise SystemExit(f"Refusing to write inside the Measurement Set: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        shutil.rmtree(out)

    excluded = args.exclude if args.exclude is not None else excluded_antennas(caldir)
    antenna = "!" + ",".join(excluded) if excluded else ""
    tables = [str(caldir / "delay.cal"), str(caldir / "bandpass.cal")]
    for path in tables:
        if not Path(path).is_dir():
            raise SystemExit(f"No calibration table {path}.")
    print(f"Solving a delay per scan on {args.fields} of {ms}"
          + (f", without {', '.join(excluded)}" if excluded else "") + " ...",
          flush=True)
    gaincal(
        vis=str(ms), caltable=str(out), field=args.fields,
        antenna=antenna, gaintype="K", solint="inf", combine="",
        refant=args.refant, minsnr=float(args.minsnr),
        gaintable=tables, gainfield=[args.bandpass_field] * 2,
        interp=["linear", "linear,linear"], parang=False,
    )
    if not out.is_dir():
        raise SystemExit("gaincal wrote no delay table.")

    tool = table()
    tool.open(str(ms / "ANTENNA"))
    names = [str(name) for name in tool.getcol("NAME")]
    tool.close()
    tool.open(str(ms / "FIELD"))
    field_names = [str(name) for name in tool.getcol("NAME")]
    tool.close()
    tool.open(str(out))
    try:
        antennas = tool.getcol("ANTENNA1")
        scans = tool.getcol("SCAN_NUMBER")
        fields = tool.getcol("FIELD_ID")
        values = tool.getcol("FPARAM")
        flags = tool.getcol("FLAG")
    finally:
        tool.close()
    rows = []
    for index in range(len(antennas)):
        rows.append((
            int(antennas[index]), int(scans[index]),
            field_names[int(fields[index])],
            [float(values[pol, 0, index]) for pol in range(values.shape[0])],
            [bool(flags[pol, 0, index]) for pol in range(flags.shape[0])],
        ))
    scan_list, table_rows = drift_table(names, rows)
    scan_fields = {row[1]: row[2] for row in rows}
    print(format_table(scan_list, table_rows, scan_fields))
    report = {
        "measurement_set": str(ms),
        "applied": tables,
        "fields": args.fields,
        "excluded": excluded,
        "refant": args.refant,
        "scans": scan_list,
        "scan_fields": {str(key): value for key, value in scan_fields.items()},
        "antennas": table_rows,
    }
    report_path = out.with_suffix(".json")
    report_path.write_text(json.dumps(report, indent=2))
    print(f"Report: {report_path}")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--ms", required=True, help="the calibrated observation; read only")
    parser.add_argument("--caldir", required=True,
                        help="the pipeline's calibration folder (delay.cal, bandpass.cal)")
    parser.add_argument("--fields", required=True,
                        help="calibrators to solve on, e.g. 0323+055,3C48")
    parser.add_argument("--bandpass-field", default="3C48",
                        help="the field the delay and bandpass tables were solved on")
    parser.add_argument("--refant", default="C00")
    parser.add_argument("--minsnr", type=float, default=5.0)
    parser.add_argument("--exclude", nargs="*",
                        help="antennas to leave out (default: those the pipeline excluded)")
    parser.add_argument("--out", help="the delay table to write (default delay_drift/<ms>.K.cal)")
    run(parser.parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
