#!/usr/bin/env python3
"""The flux density of each calibrator scan, as the calibrated data give it.

In Band 5 the fields' compact sources came out at about 0.8 of FIRST near
the pointing centres, and the phase calibrator 0323+055 was bootstrapped at
2.16 Jy against 2.9-3.0 Jy from NVSS, FIRST and TGSS (2026-10-25). This
reads the calibrated data of the calibrators scan by scan and gives the
vector average of Stokes I on the shorter baselines, for the band and for
a few sub-bands: 3C48 to be set against the Perley-Butler 2017 model,
0323+055 to show whether its flux holds from scan to scan. It also gives,
per scan, the median of the uncalibrated autocorrelations, which follow
the system temperature unless the correlator normalises them away.

Run it with CASA's Python, from the pipeline folder, for example:

    singularity exec IMAGE.sif python3 diagnose_flux_scale.py \\
        --ms Band5_Products/ingest/Band5.ms --fields 3C48,0323+055

``--uv-bins`` adds each calibrator's Stokes I against baseline length, all
scans together, sub-band by sub-band: a point source calibrated against a
point model is flat; one resolved on the longer baselines falls, and its
antenna gains then carry part of its shape into every target. ``--model``
gives a calibrator without a standard model its stated one
(``0323+055=3.08@1250,-0.76``) to compare against.

``--time-profile SCAN ...`` reads the uncorrected data of those scans, with
every flag ignored, and gives each antenna's amplitude minute by minute
against its own level in the second half of the scan. On Band 5 the
amplitude solves of the first ten minutes of the 3C48 scan failed on 22 of
26 antennas while the phase solves held (2026-10-25 run); this shows where
the data themselves settle.

The Measurement Set is only read. The report goes to ``--out`` (default
``flux_scale/<ms name>.json`` in the current folder).
"""

import argparse
import json
import math
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
# Run as `python3 diagnose_flux_scale.py` the pipeline folder comes first
# on sys.path, where its own inspect.py stands in for the standard module.
if sys.path and sys.path[0] and Path(sys.path[0]).resolve() == SCRIPT_DIR:
    sys.path.pop(0)

LIGHT_SPEED = 299792458.0
# Perley & Butler (2017), log10 S[Jy] against log10 nu[GHz].
PERLEY_BUTLER_2017 = {
    "3C48": (1.3253, -0.7553, -0.1914, 0.0498),
    "3C147": (1.4516, -0.6961, -0.2007, 0.0640, -0.0464, 0.0289),
    "3C286": (1.2481, -0.4507, -0.1798, 0.0357),
    "3C138": (1.0088, -0.4981, -0.1550, -0.0100, 0.0220),
}


def model_flux(name, freq_hz, stated=None):
    """The flux density of a calibrator: stated, else Perley-Butler 2017.

    ``stated`` maps a field to (flux in Jy, reference MHz, spectral terms),
    as ``--model`` gives it. None when neither knows the source.
    """

    if stated and name in stated:
        flux, reffreq_mhz, terms = stated[name]
        x = math.log10(freq_hz / 1e6 / reffreq_mhz)
        return flux * 10 ** sum(t * x ** (k + 1) for k, t in enumerate(terms))
    coefficients = PERLEY_BUTLER_2017.get(str(name).upper().replace(" ", ""))
    if coefficients is None:
        return None
    x = math.log10(freq_hz / 1e9)
    return 10 ** sum(c * x ** k for k, c in enumerate(coefficients))


def parse_models(values):
    """``FIELD=FLUX@MHZ[,SPIX[,CURVATURE]]`` entries into model_flux's form."""

    stated = {}
    for value in values or []:
        try:
            field, rest = value.split("=", 1)
            flux, rest = rest.split("@", 1)
            parts = [float(part) for part in rest.split(",")]
            stated[field.strip()] = (float(flux), parts[0], parts[1:])
        except (ValueError, IndexError):
            raise SystemExit(
                f"--model {value!r} is not FIELD=FLUX@MHZ[,SPIX]") from None
    return stated


def uv_profile_sums(data, flags, uvw, freqs, edges_klambda, bands):
    """Vector sums of one chunk's rows by baseline length and sub-band.

    Lengths are in kilolambda at each sub-band's own mean frequency, so a
    bin holds the same angular scale across the band. Returns sums and
    counts shaped (bin, sub-band, hand).
    """

    import numpy as np

    data = np.asarray(data)
    good = ~np.asarray(flags, dtype=bool)
    hands = [0, data.shape[0] - 1] if data.shape[0] > 1 else [0]
    data, good = data[hands], good[hands]
    metres = np.hypot(uvw[0], uvw[1])
    edges = np.linspace(0, len(freqs), bands + 1).astype(int)
    nbins = len(edges_klambda) - 1
    sums = np.zeros((nbins, bands, data.shape[0]), dtype=complex)
    counts = np.zeros((nbins, bands, data.shape[0]))
    for band in range(bands):
        part = slice(edges[band], edges[band + 1])
        klambda = metres * float(np.mean(freqs[part])) / LIGHT_SPEED / 1e3
        which = np.digitize(klambda, edges_klambda) - 1
        for k in range(nbins):
            rows = which == k
            if not rows.any():
                continue
            for corr in range(data.shape[0]):
                mask = good[corr, part][:, rows]
                sums[k, band, corr] += data[corr, part][:, rows][mask].sum()
                counts[k, band, corr] += mask.sum()
    return sums, counts


def time_profile(times, antenna1, antenna2, levels, names, step_s):
    """Each antenna's amplitude per time step against its later level.

    ``levels`` is one robust amplitude per row (the median over channels and
    hands). An antenna's level in a step is the median over its baselines;
    it is divided by its median over the second half of the scan. Returns
    the steps' start offsets (s) and, per antenna, its ratios and the first
    step from which every later one stays within 10% of 1.
    """

    import numpy as np

    times = np.asarray(times, dtype=float)
    antenna1 = np.asarray(antenna1)
    antenna2 = np.asarray(antenna2)
    levels = np.asarray(levels, dtype=float)
    start = float(times.min())
    step = np.floor((times - start) / step_s).astype(int)
    nsteps = int(step.max()) + 1
    late = times >= start + 0.5 * (float(times.max()) - start)
    antennas = {}
    for antenna in sorted(set(antenna1.tolist()) | set(antenna2.tolist())):
        own = ((antenna1 == antenna) | (antenna2 == antenna)) & np.isfinite(
            levels)
        reference = levels[own & late]
        if not reference.size or not np.median(reference) > 0:
            continue
        reference = float(np.median(reference))
        ratios = []
        for k in range(nsteps):
            values = levels[own & (step == k)]
            ratios.append(round(float(np.median(values)) / reference, 3)
                          if values.size else None)
        settled = None
        for k in range(nsteps):
            rest = [r for r in ratios[k:] if r is not None]
            if rest and all(abs(r - 1.0) <= 0.10 for r in rest):
                settled = k
                break
        # Steady only over less than the second half is not steady.
        if settled is not None and settled > nsteps // 2:
            settled = None
        name = names[antenna] if antenna < len(names) else str(antenna)
        antennas[name] = {"ratios": ratios, "settled_step": settled}
    return [round(k * step_s, 1) for k in range(nsteps)], antennas


def stokes_i(sums, counts):
    """Stokes I from the vector sums of the two hands (RR, LL or XX, YY)."""

    means = [s / n if n else None for s, n in zip(sums, counts)]
    present = [m for m in means if m is not None]
    if not present:
        return None, means
    return abs(sum(present) / len(present)), means


def scan_flux(data, flags, uvw, freqs, max_klambda, bands):
    """Vector averages of one scan's rows on baselines within the cut.

    ``data`` is (corr, chan, row), ``uvw`` (3, row), ``freqs`` per channel.
    Returns the sums and counts per sub-band and hand, to be added over
    chunks; the band is the sum of the sub-bands.
    """

    import numpy as np

    data = np.asarray(data)
    good = ~np.asarray(flags, dtype=bool)
    # The parallel hands: RR and LL of four correlations, or both of two.
    hands = [0, data.shape[0] - 1] if data.shape[0] > 1 else [0]
    data, good = data[hands], good[hands]
    reference = float(np.median(freqs))
    length = np.hypot(uvw[0], uvw[1]) * reference / LIGHT_SPEED / 1e3
    rows = length <= max_klambda
    good = good & rows[None, None, :]
    edges = np.linspace(0, len(freqs), bands + 1).astype(int)
    sums = np.zeros((bands, data.shape[0]), dtype=complex)
    counts = np.zeros((bands, data.shape[0]))
    for band in range(bands):
        part = slice(edges[band], edges[band + 1])
        for corr in range(data.shape[0]):
            mask = good[corr, part, :]
            sums[band, corr] = data[corr, part, :][mask].sum()
            counts[band, corr] = mask.sum()
    return sums, counts


def uv_profile_report(name, sums, counts, edges, band_freqs, stated,
                      skipped):
    """One calibrator's Stokes I per baseline-length bin and sub-band."""

    bins = []
    for k in range(len(edges) - 1):
        bands = []
        for band, freq in enumerate(band_freqs):
            flux, _ = stokes_i(sums[k, band], counts[k, band])
            model = model_flux(name, freq, stated)
            bands.append({
                "mhz": round(freq / 1e6, 2),
                "flux_jy": None if flux is None else round(flux, 4),
                "ratio": (round(flux / model, 4)
                          if flux is not None and model else None),
            })
        total, _ = stokes_i(sums[k].sum(axis=0), counts[k].sum(axis=0))
        bins.append({
            "klambda": [edges[k], edges[k + 1]],
            "visibilities": int(counts[k].sum()),
            "flux_jy": None if total is None else round(total, 4),
            "sub_bands": bands,
        })
    return {"skipped_scans": list(skipped or []), "bins": bins}


def time_profiles(tool, ms, scans, step_s, chunk, field_names):
    """Antenna amplitudes through each listed scan, every flag ignored.

    Per row the level is the median amplitude over the middle half of the
    band and the parallel hands of the uncorrected data, so a few channels
    of interference do not move it.
    """

    import numpy as np
    from casatools import table

    antenna_table = table()
    antenna_table.open(str(ms / "ANTENNA"))
    names = [str(name) for name in antenna_table.getcol("NAME")]
    antenna_table.close()
    profiles = {}
    for scan in scans:
        cross = tool.query(f"SCAN_NUMBER=={scan} && ANTENNA1!=ANTENNA2")
        rows = cross.nrows()
        if not rows:
            cross.close()
            raise SystemExit(f"Scan {scan} has no cross-correlations in {ms}.")
        times, first, second, levels, fields = [], [], [], [], set()
        for start in range(0, rows, chunk):
            count = min(chunk, rows - start)
            amplitude = np.abs(np.asarray(cross.getcol("DATA", start, count)))
            nchan = amplitude.shape[1]
            hands = [0, amplitude.shape[0] - 1] if amplitude.shape[0] > 1 else [0]
            middle = amplitude[hands, nchan // 4: max(nchan // 4 + 1,
                                                      3 * nchan // 4), :]
            levels.append(np.median(middle.reshape(-1, middle.shape[-1]),
                                    axis=0))
            times.append(np.asarray(cross.getcol("TIME", start, count)))
            first.append(np.asarray(cross.getcol("ANTENNA1", start, count)))
            second.append(np.asarray(cross.getcol("ANTENNA2", start, count)))
            fields.update(int(field) for field in
                          cross.getcol("FIELD_ID", start, count))
        cross.close()
        times = np.concatenate(times)
        steps, antennas = time_profile(
            times, np.concatenate(first), np.concatenate(second),
            np.concatenate(levels), names, step_s)
        array_median = []
        for k in range(len(steps)):
            values = [entry["ratios"][k] for entry in antennas.values()
                      if entry["ratios"][k] is not None]
            array_median.append(round(float(np.median(values)), 3)
                                if values else None)
        profiles[str(scan)] = {
            "field": ", ".join(field_names[field] for field in sorted(fields)),
            "time_start_s": float(times.min()),
            "step_s": step_s,
            "step_start_s": steps,
            "array_median": array_median,
            "antennas": antennas,
        }
    return profiles


def run(args):
    import numpy as np
    from casatools import table

    ms = Path(args.ms).resolve()
    out = Path(args.out or Path("flux_scale") / f"{ms.stem}.json").resolve()
    try:
        out.relative_to(ms)
    except ValueError:
        pass
    else:
        raise SystemExit(f"Refusing to write inside the Measurement Set: {out}")

    tool = table()
    tool.open(str(ms / "FIELD"))
    field_names = [str(name) for name in tool.getcol("NAME")]
    tool.close()
    tool.open(str(ms / "SPECTRAL_WINDOW"))
    freqs = np.asarray(tool.getcol("CHAN_FREQ"))[:, 0]
    tool.close()
    wanted = [name.strip() for name in args.fields.split(",") if name.strip()]
    missing = [name for name in wanted if name not in field_names]
    if missing:
        raise SystemExit(
            f"Not in {ms}: {', '.join(missing)}. Fields: {', '.join(field_names)}."
        )

    edges = np.linspace(0, len(freqs), args.bands + 1).astype(int)
    band_freqs = [float(np.mean(freqs[edges[i]:edges[i + 1]]))
                  for i in range(args.bands)]
    stated = parse_models(args.model)
    uv_edges = sorted(float(edge) for edge in (args.uv_bins or []))
    if args.uv_bins and len(uv_edges) < 2:
        raise SystemExit("--uv-bins needs at least two edges (kilolambda).")
    tool.open(str(ms))
    columns = tool.colnames()
    column = "CORRECTED_DATA" if "CORRECTED_DATA" in columns else None
    if column is None and not args.time_profile:
        tool.close()
        raise SystemExit(f"{ms} has no CORRECTED_DATA: it is not calibrated.")
    if column is None:
        # The time profile reads the uncorrected data, so an MS fresh from
        # ingest serves it; the flux measurements need the calibration.
        print(f"{ms} has no CORRECTED_DATA: only the time profile is made.")
        wanted = []
    report = {"measurement_set": str(ms), "column": column,
              "max_klambda": args.max_klambda, "band_mhz": [
                  round(f / 1e6, 2) for f in band_freqs], "fields": {}}
    if uv_edges:
        report["uv_profiles"] = {}
    try:
        for name in wanted:
            field_id = field_names.index(name)
            selected = tool.query(f"FIELD_ID=={field_id}")
            scans = sorted(set(int(s) for s in selected.getcol("SCAN_NUMBER")))
            selected.close()
            entries = []
            if uv_edges:
                uv_sums = np.zeros((len(uv_edges) - 1, args.bands, 2),
                                   dtype=complex)
                uv_counts = np.zeros(uv_sums.shape)
            for scan in scans:
                profiled = bool(uv_edges) and scan not in (
                    args.profile_skip_scans or [])
                cross = tool.query(
                    f"FIELD_ID=={field_id} && SCAN_NUMBER=={scan} "
                    "&& ANTENNA1!=ANTENNA2")
                sums = np.zeros((args.bands, 2), dtype=complex)
                counts = np.zeros((args.bands, 2))
                times = cross.getcol("TIME")
                rows = cross.nrows()
                for start in range(0, rows, args.chunk):
                    count = min(args.chunk, rows - start)
                    data = cross.getcol(column, start, count)
                    flags = cross.getcol("FLAG", start, count)
                    uvw = cross.getcol("UVW", start, count)
                    s, n = scan_flux(data, flags, uvw, freqs,
                                     args.max_klambda, args.bands)
                    sums[:, :s.shape[1]] += s
                    counts[:, :n.shape[1]] += n
                    if profiled:
                        s, n = uv_profile_sums(data, flags, uvw, freqs,
                                               uv_edges, args.bands)
                        uv_sums[:, :, :s.shape[2]] += s
                        uv_counts[:, :, :n.shape[2]] += n
                cross.close()
                autos = tool.query(
                    f"FIELD_ID=={field_id} && SCAN_NUMBER=={scan} "
                    "&& ANTENNA1==ANTENNA2")
                auto_level = None
                if autos.nrows():
                    raw = np.abs(autos.getcol("DATA"))
                    raw_flags = autos.getcol("FLAG").astype(bool)
                    values = raw[~raw_flags]
                    auto_level = float(np.median(values)) if values.size else None
                autos.close()
                bands = []
                for band in range(args.bands):
                    flux, _ = stokes_i(sums[band], counts[band])
                    model = model_flux(name, band_freqs[band], stated)
                    bands.append({
                        "mhz": round(band_freqs[band] / 1e6, 2),
                        "flux_jy": None if flux is None else round(flux, 4),
                        "model_jy": None if model is None else round(model, 4),
                        "ratio": (round(flux / model, 4)
                                  if flux is not None and model else None),
                    })
                total, _ = stokes_i(sums.sum(axis=0), counts.sum(axis=0))
                entries.append({
                    "scan": scan,
                    "time_start_s": float(np.min(times)) if len(times) else None,
                    "visibilities": int(counts.sum()),
                    "flux_jy": None if total is None else round(total, 4),
                    "sub_bands": bands,
                    "autocorrelation_median": auto_level,
                })
            report["fields"][name] = entries
            if uv_edges:
                report["uv_profiles"][name] = uv_profile_report(
                    name, uv_sums, uv_counts, uv_edges, band_freqs, stated,
                    args.profile_skip_scans)
        if args.time_profile:
            report["time_profiles"] = time_profiles(
                tool, ms, args.time_profile, args.time_step_min * 60.0,
                args.chunk, field_names)
    finally:
        tool.close()

    print(f"Vector average of Stokes I within {args.max_klambda:g} klambda, "
          f"{column}:")
    for name, entries in report["fields"].items():
        for entry in entries:
            parts = [f"{band['flux_jy']}" + (
                f" ({band['ratio']:.3f} of the model)"
                if band["ratio"] else "")
                for band in entry["sub_bands"]]
            print(f"  {name:10} scan {entry['scan']:>3}: band "
                  f"{entry['flux_jy']} Jy; sub-bands "
                  f"{', '.join(parts)}; autocorrelations "
                  f"{entry['autocorrelation_median']}")
    for name, profile in report.get("uv_profiles", {}).items():
        print(f"{name} against baseline length (all scans"
              + (f" but {profile['skipped_scans']}"
                 if profile["skipped_scans"] else "") + "):")
        for entry in profile["bins"]:
            ratios = [band["ratio"] for band in entry["sub_bands"]]
            print(f"  {entry['klambda'][0]:g}-{entry['klambda'][1]:g} "
                  f"klambda, {entry['visibilities']} visibilities: "
                  f"{entry['flux_jy']} Jy; of the model by sub-band "
                  + ", ".join("-" if r is None else f"{r:.3f}"
                              for r in ratios))
    for scan, profile in report.get("time_profiles", {}).items():
        late = {name: entry["settled_step"]
                for name, entry in profile["antennas"].items()
                if entry["settled_step"] != 0}
        print(f"Scan {scan} ({profile['field']}), amplitude per "
              f"{args.time_step_min:g} min against the second half: "
              + ("every antenna within 10% from the start" if not late else
                 "within 10% only from " + ", ".join(
                     f"{name} " + ("never" if k is None else
                                   f"{profile['step_start_s'][k] / 60:g} min")
                     for name, k in sorted(
                         late.items(),
                         key=lambda item: -(1e9 if item[1] is None
                                            else item[1])))))
        print("  array median by step: " + ", ".join(
            "-" if r is None else f"{r:.2f}"
            for r in profile["array_median"]))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"Report: {out}")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--ms", required=True,
                        help="the calibrated observation; read only")
    parser.add_argument("--fields", default="3C48,0323+055",
                        help="calibrators to measure (default 3C48,0323+055)")
    parser.add_argument("--max-klambda", type=float, default=5.0,
                        help="baselines up to this length at the band's "
                        "middle (default 5: the central square, on which "
                        "both calibrators are points)")
    parser.add_argument("--bands", type=int, default=4,
                        help="sub-bands to report (default 4)")
    parser.add_argument("--chunk", type=int, default=4000,
                        help="rows read at a time (default 4000)")
    parser.add_argument("--model", action="append", metavar="FIELD=FLUX@MHZ",
                        help="a stated model for a calibrator without a "
                        "standard one, e.g. 0323+055=3.08@1250,-0.76 "
                        "(repeatable)")
    parser.add_argument("--uv-bins", nargs="+", type=float,
                        metavar="KLAMBDA",
                        help="edges of baseline-length bins in kilolambda; "
                        "adds each calibrator's Stokes I per bin, all scans "
                        "together")
    parser.add_argument("--profile-skip-scans", nargs="+", type=int,
                        metavar="SCAN",
                        help="scans left out of the --uv-bins profile")
    parser.add_argument("--time-profile", nargs="+", type=int,
                        metavar="SCAN",
                        help="scans whose uncorrected amplitudes to follow "
                        "antenna by antenna, flags ignored")
    parser.add_argument("--time-step-min", type=float, default=1.0,
                        help="time step of --time-profile in minutes "
                        "(default 1)")
    parser.add_argument("--out", help="the report (default "
                        "flux_scale/<ms>.json)")
    run(parser.parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
