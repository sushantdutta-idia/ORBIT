#!/usr/bin/env python3
"""How noisy each range of baseline length is, against what its weights say.

The tapered Band 5 images (2026-10-28) were far noisier than the data's
weights promise: at 15" even their Stokes V was 85-200 µJy/beam, four to
seven times the natural-weight figure, so the short baselines that a taper
favours carry more noise than the long ones. Whether their weights know
it decides the remedy: weights that do not are to be recomputed, weights
that do leave the data to be flagged or the shortest baselines left out.

For one calibrated Measurement Set this takes the difference of
neighbouring channels, which removes the sky and keeps the noise, in
Stokes I and V, scaled by the noise the weights predict for it, and gives
for each range of baseline length (kilo-wavelengths) the scatter of that
scaled difference: one where the weights are right. It also gives the
scatter itself in Jy and the share of the unflagged data in each range.

Run it with CASA's Python, from the pipeline folder:

    singularity exec IMAGE.sif python3 diagnose_baseline_noise.py \\
        --ms Band5_targets_0323/G0222-0321/target_whole_field/G0222-0321/round3.ms

The Measurement Set is only read. The report goes to ``--out`` (default
``baseline_noise/<ms name>.json``).
"""

import argparse
import json
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
# Run as `python3 diagnose_baseline_noise.py` the pipeline folder comes
# first on sys.path, where its own inspect.py stands in for the standard one.
if sys.path and sys.path[0] and Path(sys.path[0]).resolve() == SCRIPT_DIR:
    sys.path.pop(0)

LIGHT_SPEED = 299792458.0
DEFAULT_EDGES = (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 200.0)


def robust_sigma(values):
    """The scatter of real values from their median absolute deviation."""

    import numpy as np

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 10:
        return None
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def differences(data, flags, weights):
    """Channel differences of Stokes I and V, raw and scaled by their noise.

    ``data`` and ``flags`` are (corr, chan, row) with the two parallel
    hands first and last; ``weights`` the same shape, or (corr, row). A
    weight is taken as one over the variance of a visibility's real part.
    Returns, for I and V, the real parts of the differences in Jy and the
    same divided by the scatter the weights predict for them.
    """

    import numpy as np

    data = np.asarray(data)
    flags = np.asarray(flags, dtype=bool)
    weights = np.asarray(weights, dtype=float)
    if weights.ndim == 2:
        weights = np.broadcast_to(weights[:, None, :], data.shape)
    hands = [0, data.shape[0] - 1]
    first, last = data[hands[0]], data[hands[1]]
    good = ~(flags[hands[0]] | flags[hands[1]])
    good &= (weights[hands[0]] > 0) & (weights[hands[1]] > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        variance = 1.0 / weights[hands[0]] + 1.0 / weights[hands[1]]
    out = {}
    for name, value in (("I", 0.5 * (first + last)), ("V", 0.5 * (first - last))):
        step = value[1:] - value[:-1]
        both = good[1:] & good[:-1]
        # Each of I and V carries a quarter of the two hands' variance;
        # a difference of two channels doubles it.
        predicted = np.sqrt(0.5 * (variance[1:] + variance[:-1]) * 0.5)
        raw = step.real[both]
        out[name] = (raw, raw / predicted[both])
    return out


def run(args):
    import numpy as np
    from casatools import table

    ms = Path(args.ms).resolve()
    out = Path(args.out or Path("baseline_noise") / f"{ms.stem}.json").resolve()
    try:
        out.relative_to(ms)
    except ValueError:
        pass
    else:
        raise SystemExit(f"Refusing to write inside the Measurement Set: {out}")

    tool = table()
    tool.open(str(ms / "SPECTRAL_WINDOW"))
    freqs = np.asarray(tool.getcol("CHAN_FREQ"))[:, 0]
    tool.close()
    edges = list(args.edges)
    bins = len(edges) - 1
    raw = {name: [[] for _ in range(bins)] for name in ("I", "V")}
    scaled = {name: [[] for _ in range(bins)] for name in ("I", "V")}
    counts = np.zeros(bins)
    tool.open(str(ms))
    try:
        columns = tool.colnames()
        column = args.column if args.column in columns else "DATA"
        weight_column = (
            "WEIGHT_SPECTRUM" if "WEIGHT_SPECTRUM" in columns
            and tool.iscelldefined("WEIGHT_SPECTRUM", 0) else "WEIGHT"
        )
        rows = tool.nrows()
        reference = float(np.median(freqs))
        rng = np.random.default_rng(1)
        for start in range(0, rows, args.chunk):
            count = min(args.chunk, rows - start)
            a1 = tool.getcol("ANTENNA1", start, count)
            a2 = tool.getcol("ANTENNA2", start, count)
            cross = a1 != a2
            if not cross.any():
                continue
            uvw = tool.getcol("UVW", start, count)
            length = np.hypot(uvw[0], uvw[1]) * reference / LIGHT_SPEED / 1e3
            data = tool.getcol(column, start, count)
            flags = tool.getcol("FLAG", start, count)
            weights = tool.getcol(weight_column, start, count)
            for b in range(bins):
                pick = cross & (length >= edges[b]) & (length < edges[b + 1])
                if not pick.any():
                    continue
                w = weights[..., pick] if weights.ndim == 3 else weights[:, pick]
                found = differences(data[..., pick], flags[..., pick], w)
                counts[b] += int((~flags[..., pick]).sum())
                for name, (r, s) in found.items():
                    if r.size > args.sample:
                        keep = rng.choice(r.size, args.sample, replace=False)
                        r, s = r[keep], s[keep]
                    raw[name][b].append(r)
                    scaled[name][b].append(s)
    finally:
        tool.close()

    total = counts.sum() or 1.0
    table_rows = []
    for b in range(bins):
        entry = {
            "klambda": [edges[b], edges[b + 1]],
            "share_of_data": round(float(counts[b] / total), 4),
        }
        for name in ("I", "V"):
            r = np.concatenate(raw[name][b]) if raw[name][b] else np.array([])
            s = np.concatenate(scaled[name][b]) if scaled[name][b] else np.array([])
            entry[f"{name}_jy"] = robust_sigma(r)
            entry[f"{name}_over_weights"] = robust_sigma(s)
        table_rows.append(entry)

    print(f"{ms}: channel differences, {column}, weights from {weight_column}")
    print(f"{'klambda':>12} {'data':>6} {'I/weights':>10} {'V/weights':>10} "
          f"{'I Jy':>10} {'V Jy':>10}")
    for entry in table_rows:
        if not entry["share_of_data"]:
            continue
        low, high = entry["klambda"]

        def show(value, digits=2):
            return f"{value:.{digits}f}" if value is not None else "-"

        print(f"{low:5.1f}-{high:<6.1f} {100 * entry['share_of_data']:5.1f}% "
              f"{show(entry['I_over_weights']):>10} "
              f"{show(entry['V_over_weights']):>10} "
              f"{show(entry['I_jy'], 4):>10} {show(entry['V_jy'], 4):>10}")
    report = {"measurement_set": str(ms), "column": column,
              "weights": weight_column, "bins": table_rows}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"Report: {out}")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--ms", required=True, help="calibrated data; read only")
    parser.add_argument("--column", default="DATA",
                        help="data column (default DATA: a split copy)")
    parser.add_argument("--edges", type=float, nargs="+", default=DEFAULT_EDGES,
                        metavar="KLAMBDA",
                        help="edges of the ranges of baseline length")
    parser.add_argument("--chunk", type=int, default=3000)
    parser.add_argument("--sample", type=int, default=20000,
                        help="values kept per range and chunk")
    parser.add_argument("--out", help="the report")
    run(parser.parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
