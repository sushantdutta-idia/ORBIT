#!/usr/bin/env python3
"""Headless, CASA-friendly QA plotting helpers for the uGMRT pipeline."""

import json
import math
import warnings
from pathlib import Path


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of producing the old plots.
PIPELINE_INTERFACE_VERSION = 3


def _plot_backend():
    """Return pyplot or an explanatory error string without requiring a GUI."""

    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        return plt, None
    except Exception as exc:
        return None, str(exc)


def _numpy():
    try:
        import numpy as np

        return np
    except Exception as exc:
        raise RuntimeError(f"NumPy is unavailable: {exc}") from exc


def write_status(path, status):
    """Write a JSON status file that explains available and missing plots."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as handle:
        json.dump(status, handle, indent=2, default=str)


def plot_spw_coverage(spws, outfile, title="Spectral-window frequency coverage"):
    """Plot frequency extents supplied by an inspection-report ``spws`` map."""

    plt, error = _plot_backend()
    if plt is None:
        raise RuntimeError(f"Matplotlib is unavailable: {error}")

    rows = []
    for key, value in (spws or {}).items():
        if not isinstance(value, dict):
            continue
        try:
            spw_id = int(value.get("id", key))
            low = float(value["fmin_MHz"])
            high = float(value["fmax_MHz"])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(low) and math.isfinite(high):
            rows.append((spw_id, min(low, high), max(low, high)))

    if not rows:
        return None

    rows.sort()
    figure, axis = plt.subplots(figsize=(9, max(3.5, 0.55 * len(rows) + 1)))
    for index, (spw_id, low, high) in enumerate(rows):
        axis.barh(index, high - low, left=low, height=0.55, color="tab:blue")
        axis.text(high, index, f"  SPW {spw_id}", va="center", fontsize=8)

    axis.set_yticks([])
    axis.set_xlabel("Frequency (MHz)")
    axis.set_title(title)
    axis.grid(axis="x", alpha=0.25)

    outfile = Path(outfile)
    figure.savefig(str(outfile), dpi=160, bbox_inches="tight")
    plt.close(figure)
    return outfile


def _read_caltable(caltable):
    """Read what calibration QA needs from a CASA calibration table.

    Besides the solutions themselves: their flags (flagged solutions are
    stored as 1+0j, and drawn as real ones they put false lines at amplitude
    1 and phase 0), the field of each solution, and the antenna, field and
    frequency subtables, so plots can say what they show.
    """

    try:
        from casatools import table
    except Exception as exc:
        raise RuntimeError(f"CASA table tool is unavailable: {exc}") from exc

    def subtable_column(name, column):
        sub = table()
        try:
            sub.open(str(Path(caltable) / name))
            return sub.getcol(column)
        except Exception:
            return None
        finally:
            try:
                sub.close()
            except Exception:
                pass

    tool = table()
    opened = False
    try:
        tool.open(str(caltable))
        opened = True
        columns = set(tool.colnames())
        parameter_column = (
            "CPARAM" if "CPARAM" in columns
            else "FPARAM" if "FPARAM" in columns
            else None
        )
        if parameter_column is None:
            raise RuntimeError(
                "Neither CPARAM nor FPARAM is present in the calibration table."
            )

        result = {
            "parameters": tool.getcol(parameter_column),
            "parameter_column": parameter_column,
        }
        for key, column in (
            ("time", "TIME"),
            ("antenna", "ANTENNA1"),
            ("field", "FIELD_ID"),
            ("spw", "SPECTRAL_WINDOW_ID"),
            ("flags", "FLAG"),
        ):
            result[key] = tool.getcol(column) if column in columns else None
    finally:
        if opened:
            tool.close()

    result["antenna_names"] = subtable_column("ANTENNA", "NAME")
    result["field_names"] = subtable_column("FIELD", "NAME")
    frequencies = subtable_column("SPECTRAL_WINDOW", "CHAN_FREQ")
    result["channel_frequencies"] = frequencies
    return result


def _masked_solutions(parameters, flags):
    """Solutions as a (correlation, channel, row) array, flagged set to NaN."""

    np = _numpy()
    values = np.asarray(parameters)
    if values.ndim == 1:
        values = values[np.newaxis, np.newaxis, :]
    elif values.ndim == 2:
        values = values[:, np.newaxis, :]
    values = values.astype(complex if np.iscomplexobj(values) else float)
    if flags is not None:
        flags = np.asarray(flags, dtype=bool).reshape(values.shape)
        values = np.where(flags, np.nan, values)
    return values


def _time_series(parameters, flags=None):
    """One gain per correlation per row, flagged solutions as NaN."""

    return _masked_solutions(parameters, flags)[:, 0, :]


def _bandpass_series(parameters, quantity, flags=None):
    """Per-correlation median and 16th/84th percentiles across solutions.

    Returns (correlation, statistic, channel) with statistic 0 = median,
    1 = 16th and 2 = 84th percentile, or for phase the circular mean and
    the circular spread; flagged channels are NaN, so they plot as gaps.
    """

    np = _numpy()
    values = _masked_solutions(parameters, flags)
    if values.shape[1] < 2:
        return None

    if quantity == "amplitude":
        amplitudes = np.abs(values)
        # Fully flagged channels are all-NaN on purpose: they plot as gaps.
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.stack([
                np.nanmedian(amplitudes, axis=-1),
                np.nanpercentile(amplitudes, 16.0, axis=-1),
                np.nanpercentile(amplitudes, 84.0, axis=-1),
            ], axis=1)

    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        unit = values / np.abs(values)
        mean = np.nanmean(unit, axis=-1)
    centre = np.degrees(np.angle(mean))
    # Circular standard deviation, from the mean resultant length.
    length = np.clip(np.abs(mean), 1e-12, 1.0)
    spread = np.degrees(np.sqrt(-2.0 * np.log(length)))
    return np.stack([centre, centre - spread, centre + spread], axis=1)


def plot_caltable_time_series(caltable, outfile, quantity, title):
    """Gain amplitude or phase against time, by calibrator field.

    Flagged solutions are left out, both polarizations are shown, and each
    calibrator has its own colour, so a flux calibrator at the start and
    end is not mistaken for a jump in the phase calibrator's gains.
    """

    plt, error = _plot_backend()
    if plt is None:
        raise RuntimeError(f"Matplotlib is unavailable: {error}")

    data = _read_caltable(caltable)
    np = _numpy()
    values = _time_series(data["parameters"], data.get("flags"))
    if values.size == 0:
        return None

    if quantity == "amplitude":
        values = np.abs(values)
        ylabel = "Gain amplitude"
    elif quantity == "phase":
        values = np.degrees(np.angle(values))
        ylabel = "Gain phase relative to the reference antenna (deg)"
    else:
        raise ValueError(f"Unsupported gain quantity: {quantity}")

    times = data.get("time")
    if times is None:
        x_values = np.arange(values.shape[-1], dtype=float)
        xlabel = "Solution row"
    else:
        x_values = np.asarray(times, dtype=float).reshape(-1)
        x_values = (x_values - np.nanmin(x_values)) / 60.0
        xlabel = "Time since first solution (min)"

    fields = data.get("field")
    fields = (
        np.zeros(values.shape[-1], dtype=int) if fields is None
        else np.asarray(fields, dtype=int).reshape(-1)
    )
    names = list(data.get("field_names") or [])

    figure, axis = plt.subplots(figsize=(9, 5))
    colours = plt.get_cmap("tab10")
    markers = ["o", "^", "s", "v"]
    plotted = False
    for index, field in enumerate(sorted(set(fields.tolist()))):
        rows = fields == field
        label = names[field] if 0 <= field < len(names) else f"field {field}"
        for correlation in range(values.shape[0]):
            y = values[correlation, rows]
            x = x_values[rows]
            keep = np.isfinite(y) & np.isfinite(x)
            if not np.any(keep):
                continue
            axis.scatter(
                x[keep], y[keep], s=4, alpha=0.5,
                color=colours(index % 10),
                marker=markers[correlation % len(markers)],
                edgecolors="none",
                label=f"{label}, pol {correlation}",
            )
            plotted = True
    if not plotted:
        plt.close(figure)
        return None

    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.set_title(title + " (flagged solutions omitted)")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=7, markerscale=2, loc="best")

    outfile = Path(outfile)
    figure.savefig(str(outfile), dpi=160, bbox_inches="tight")
    plt.close(figure)
    return outfile


def plot_caltable_bandpass(caltable, outfile, quantity, title):
    """A bandpass's median across antennas and solutions, with its spread.

    Flagged channels are gaps, both polarizations are drawn, and the x axis
    is frequency when the table records it.
    """

    plt, error = _plot_backend()
    if plt is None:
        raise RuntimeError(f"Matplotlib is unavailable: {error}")

    data = _read_caltable(caltable)
    np = _numpy()
    series = _bandpass_series(data["parameters"], quantity, data.get("flags"))
    if series is None:
        return None

    nchan = series.shape[-1]
    x_values = np.arange(nchan, dtype=float)
    xlabel = "Channel"
    frequencies = data.get("channel_frequencies")
    if frequencies is not None:
        frequencies = np.asarray(frequencies, dtype=float)
        if frequencies.ndim == 2 and frequencies.shape[0] == nchan:
            x_values = frequencies[:, 0] / 1e6
            xlabel = "Frequency (MHz)"
        elif frequencies.ndim == 1 and frequencies.size == nchan:
            x_values = frequencies / 1e6
            xlabel = "Frequency (MHz)"

    figure, axis = plt.subplots(figsize=(9, 5))
    colours = ["tab:green", "tab:purple", "tab:orange", "tab:red"]
    plotted = False
    for correlation in range(series.shape[0]):
        centre, low, high = series[correlation]
        if not np.any(np.isfinite(centre)):
            continue
        colour = colours[correlation % len(colours)]
        axis.plot(x_values, centre, linewidth=1.0, color=colour,
                  label=f"pol {correlation}")
        axis.fill_between(x_values, low, high, color=colour, alpha=0.18,
                          linewidth=0)
        plotted = True
    if not plotted:
        plt.close(figure)
        return None

    axis.set_xlabel(xlabel)
    axis.set_ylabel(
        "Bandpass amplitude" if quantity == "amplitude"
        else "Bandpass phase (deg)"
    )
    axis.set_title(
        title + (
            " (band: 16th-84th percentile across antennas)"
            if quantity == "amplitude"
            else " (band: circular spread across antennas)"
        ),
        fontsize=10,
    )
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)

    outfile = Path(outfile)
    figure.savefig(str(outfile), dpi=160, bbox_inches="tight")
    plt.close(figure)
    return outfile


def phase_transfer_stability(times, antennas, fields, phases_deg, flags=None,
                             reference=None):
    """How far each antenna's phase moves between calibrator visits.

    The target's phases are interpolated between scan-averaged solutions
    on the phase calibrator, so the change between consecutive solutions
    bounds how wrong that interpolation can be. Differences are wrapped to
    +/-180 degrees; the reference antenna, fixed at zero, is left out.
    Returns, per field, the median gap and the median and 90th-percentile
    absolute phase change, plus the antennas that change most.
    """

    np = _numpy()
    times = np.asarray(times, dtype=float).reshape(-1)
    antennas = np.asarray(antennas, dtype=int).reshape(-1)
    fields = np.asarray(fields, dtype=int).reshape(-1)
    phases = np.asarray(phases_deg, dtype=float)
    if phases.ndim == 1:
        phases = phases[np.newaxis, :]
    if flags is not None:
        flags = np.asarray(flags, dtype=bool).reshape(phases.shape)
        phases = np.where(flags, np.nan, phases)

    results = {}
    for field in sorted(set(fields.tolist())):
        changes = []
        gaps = []
        per_antenna = {}
        for antenna in sorted(set(antennas[fields == field].tolist())):
            if reference is not None and antenna == reference:
                continue
            rows = np.where((fields == field) & (antennas == antenna))[0]
            if rows.size < 2:
                continue
            rows = rows[np.argsort(times[rows])]
            for correlation in range(phases.shape[0]):
                series = phases[correlation, rows]
                keep = np.isfinite(series)
                if np.count_nonzero(keep) < 2:
                    continue
                values = series[keep]
                if np.all(np.abs(values) < 1e-6):
                    continue  # a reference antenna, fixed at zero
                stamps = times[rows][keep]
                delta = (np.diff(values) + 180.0) % 360.0 - 180.0
                changes.extend(np.abs(delta).tolist())
                gaps.extend((np.diff(stamps) / 60.0).tolist())
                per_antenna.setdefault(antenna, []).extend(
                    np.abs(delta).tolist()
                )
        if not changes:
            continue
        worst = sorted(
            ((float(np.median(v)), a) for a, v in per_antenna.items()),
            reverse=True,
        )[:5]
        results[field] = {
            "transitions": len(changes),
            "median_gap_min": float(np.median(gaps)),
            "median_abs_change_deg": float(np.median(changes)),
            "p90_abs_change_deg": float(np.percentile(changes, 90.0)),
            "worst_antennas": [
                {"antenna": int(a), "median_abs_change_deg": m}
                for m, a in worst
            ],
        }
    return results


def antenna_gain_outliers(amplitudes, antennas, flags=None, low=0.5,
                          high=2.0):
    """Antenna/polarization pairs whose typical gain is far from the array's.

    ``amplitudes`` is (correlation, row). A gain several times below the
    array's median usually means an antenna or receiver chain that is not
    working as the others are; its data may still be usable with statistical
    weights, but it deserves a look.
    """

    np = _numpy()
    values = np.asarray(amplitudes, dtype=float)
    if values.ndim == 1:
        values = values[np.newaxis, :]
    antennas = np.asarray(antennas, dtype=int).reshape(-1)
    if flags is not None:
        values = np.where(
            np.asarray(flags, dtype=bool).reshape(values.shape), np.nan, values
        )

    medians = {}
    for antenna in sorted(set(antennas.tolist())):
        for correlation in range(values.shape[0]):
            series = values[correlation, antennas == antenna]
            series = series[np.isfinite(series)]
            if series.size:
                medians[(antenna, correlation)] = float(np.median(series))
    if not medians:
        return {"array_median": None, "outliers": []}

    array_median = float(np.median(list(medians.values())))
    outliers = [
        {
            "antenna": antenna,
            "correlation": correlation,
            "median_gain": value,
            "ratio_to_array": value / array_median,
        }
        for (antenna, correlation), value in sorted(medians.items())
        if array_median > 0 and not (
            low * array_median <= value <= high * array_median
        )
    ]
    return {"array_median": array_median, "outliers": outliers}


def _name_list(values):
    """Names as a list, whether CASA returned a list, an array or nothing.

    ``values or []`` raises on a numpy array of more than one name, which is
    what the table tool returns: the gain-outlier summary was computed and
    then lost to that on every run with more than one antenna.
    """

    if values is None:
        return []
    return [str(value) for value in list(values)]


def caltable_stability_summary(phase_scan_table, flux_table):
    """Read two calibration tables and summarise what they say."""

    np = _numpy()
    summary = {}
    try:
        data = _read_caltable(phase_scan_table)
        values = _time_series(data["parameters"], None)
        summary["phase_transfer"] = phase_transfer_stability(
            data["time"], data["antenna"], data["field"],
            np.degrees(np.angle(values)),
            None if data.get("flags") is None
            else np.asarray(data["flags"])[:, 0, :],
        )
        names = _name_list(data.get("field_names"))
        summary["field_names"] = names
    except Exception as exc:
        summary["phase_transfer_error"] = str(exc)
    try:
        data = _read_caltable(flux_table)
        values = np.abs(_time_series(data["parameters"], None))
        flags = (
            None if data.get("flags") is None
            else np.asarray(data["flags"])[:, 0, :]
        )
        result = antenna_gain_outliers(values, data["antenna"], flags)
        names = _name_list(data.get("antenna_names"))
        for item in result["outliers"]:
            index = item["antenna"]
            item["name"] = names[index] if 0 <= index < len(names) else None
        summary["antenna_gains"] = result
    except Exception as exc:
        summary["antenna_gains_error"] = str(exc)
    return summary


def radio_display_scaling(pixels, np=None):
    """Noise and display limits suited to a radio image.

    A percentile stretch cannot show a radio image: well over 99% of pixels
    are noise, so the 1st to 99.5th percentiles span about -2.3 to +2.6
    sigma, every source saturates, and a thermal-noise-limited image looks
    like static. The noise is instead measured robustly (median absolute
    deviation) near the centre of the field, where a primary-beam correction
    has barely amplified it, and the display runs from -3 sigma
    to the peak, linear within 3 sigma and logarithmic above, so the noise
    floor, the artefacts and the sources can all be seen at once.
    """

    np = np or _numpy()
    pixels = np.asarray(pixels, dtype=float)
    finite = pixels[np.isfinite(pixels)]
    if finite.size == 0:
        return None

    # The middle quarter of each axis: close enough to the pointing centre
    # that a primary-beam correction has barely raised the noise, and still
    # tens of thousands of pixels in any real image.
    rows, columns = pixels.shape
    centre = pixels[
        rows // 2 - rows // 8: rows // 2 + rows // 8,
        columns // 2 - columns // 8: columns // 2 + columns // 8,
    ]
    centre = centre[np.isfinite(centre)]
    sample = centre if centre.size >= 100 else finite
    median = float(np.median(sample))
    sigma = 1.4826 * float(np.median(np.abs(sample - median)))
    peak = float(np.max(finite))
    if not np.isfinite(sigma) or sigma <= 0:
        low, high = float(np.min(finite)), float(np.max(finite))
        return {"sigma": None, "peak": peak, "vmin": low, "vmax": high}

    return {
        "sigma": sigma,
        "peak": peak,
        "dynamic_range": peak / sigma,
        "vmin": -3.0 * sigma,
        "vmax": max(peak, 10.0 * sigma),
        "linthresh": 3.0 * sigma,
    }


def plot_casa_image(imagename, outfile, title):
    """Render a decimated CASA image plane without loading a huge full image."""

    plt, error = _plot_backend()
    if plt is None:
        raise RuntimeError(f"Matplotlib is unavailable: {error}")

    try:
        from casatools import image
    except Exception as exc:
        raise RuntimeError(f"CASA image tool is unavailable: {exc}") from exc

    np = _numpy()
    tool = image()
    opened = False
    try:
        tool.open(str(imagename))
        opened = True
        shape = [int(value) for value in tool.shape()]
        if len(shape) < 2:
            return None

        blc = [0] * len(shape)
        trc = [max(0, value - 1) for value in shape]
        increment = [1] * len(shape)
        for axis in range(2):
            increment[axis] = max(1, int(math.ceil(shape[axis] / 2048.0)))
        for axis in range(2, len(shape)):
            trc[axis] = 0

        pixels = np.asarray(
            tool.getchunk(blc=blc, trc=trc, inc=increment),
            dtype=float,
        )
    finally:
        if opened:
            tool.close()

    pixels = np.squeeze(pixels)
    if pixels.ndim != 2:
        return None

    scaling = radio_display_scaling(pixels, np)
    if scaling is None:
        return None

    norm = None
    if scaling.get("linthresh"):
        from matplotlib.colors import SymLogNorm
        norm = SymLogNorm(
            linthresh=scaling["linthresh"],
            linscale=1.0,
            vmin=scaling["vmin"],
            vmax=scaling["vmax"],
            base=10,
        )

    figure, axis = plt.subplots(figsize=(7, 6))
    rendered = axis.imshow(
        pixels.T,
        origin="lower",
        cmap="inferno",
        norm=norm,
        vmin=None if norm else scaling["vmin"],
        vmax=None if norm else scaling["vmax"],
        interpolation="nearest",
    )
    axis.set_xlabel("Pixel")
    axis.set_ylabel("Pixel")
    if scaling.get("sigma"):
        title = (
            f"{title}\ncentral noise {scaling['sigma'] * 1e6:.1f} uJy/beam, "
            f"peak {scaling['peak'] * 1e3:.2f} mJy/beam, "
            f"peak/noise {scaling['dynamic_range']:.0f}"
        )
    axis.set_title(title, fontsize=10)
    colorbar = figure.colorbar(rendered, ax=axis)
    colorbar.set_label("Jy/beam (linear within 3 sigma, logarithmic above)")

    outfile = Path(outfile)
    figure.savefig(str(outfile), dpi=160, bbox_inches="tight")
    plt.close(figure)
    return outfile


def plot_dd_directions(directions, outfile, title="DD-calibration directions"):
    """Plot selected direction-dependent calibration centres on the sky."""

    plt, error = _plot_backend()
    if plt is None:
        raise RuntimeError(f"Matplotlib is unavailable: {error}")

    rows = []
    for direction in directions or []:
        try:
            if isinstance(direction, dict):
                name = direction.get("name", "direction")
                ra = float(direction["ra_deg"])
                dec = float(direction["dec_deg"])
                flux = float(direction.get("flux_jy", 0.0))
            else:
                name = direction.name
                ra = float(direction.ra_deg)
                dec = float(direction.dec_deg)
                flux = float(direction.flux_jy)
        except (AttributeError, KeyError, TypeError, ValueError):
            continue

        if math.isfinite(ra) and math.isfinite(dec) and math.isfinite(flux):
            rows.append((str(name), ra, dec, max(0.0, flux)))

    if not rows:
        return None

    np = _numpy()
    names, ra, dec, flux = zip(*rows)
    flux = np.asarray(flux, dtype=float)
    sizes = 40.0 + 260.0 * flux / max(float(np.nanmax(flux)), 1e-6)

    figure, axis = plt.subplots(figsize=(7, 6))
    scatter = axis.scatter(ra, dec, s=sizes, c=flux, cmap="viridis", alpha=0.8)
    for name, x_value, y_value in zip(names, ra, dec):
        axis.annotate(name, (x_value, y_value), xytext=(4, 4),
                      textcoords="offset points", fontsize=8)
    axis.set_xlabel("Right ascension (deg)")
    axis.set_ylabel("Declination (deg)")
    axis.set_title(title)
    axis.invert_xaxis()
    axis.grid(alpha=0.25)
    colorbar = figure.colorbar(scatter, ax=axis)
    colorbar.set_label("Direction flux density (Jy)")

    outfile = Path(outfile)
    figure.savefig(str(outfile), dpi=160, bbox_inches="tight")
    plt.close(figure)
    return outfile
