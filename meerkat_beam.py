#!/usr/bin/env python3
"""MeerKAT primary-beam correction for FITS images written by WSClean.

WSClean's own beam correction (EveryBeam) is validated here only for LOFAR, so
MeerKAT images are corrected afterwards with the measured katbeam model, the
same model the CASA backend applies. The two backends therefore share a beam
model and an image-plane convention: for the same pixel they evaluate the
same katbeam function at the same sky offsets.

What differs, deliberately, is how the frequency dependence is handled. A
MeerKAT band spans a factor of two, and a primary beam's width scales as
1/frequency, so a single-frequency beam is wrong across most of an MFS image.

* Each sub-band image is divided by the beam at its own effective frequency:
  the centre of the part of the sub-band that survived flagging, both the
  science-range edges and the known-RFI masks.
* The MFS image is divided by the band-averaged beam, which is exact for a
  flat-spectrum source. WSClean makes the MFS image as the mean of the
  sub-band images weighted by their imaging weights (WSCIMGWG), so such a
  source is attenuated by the same weighted mean of their beams, and those
  are the weights used. Under Briggs weighting they rise steeply with
  frequency: on MeerKAT L-band 1708779077 the 1444 MHz sub-band carried 2.8
  times the weight of the 1016 MHz one, and an equal-weight average beam left
  the corrected MFS image 12% low 35 arcmin out. Spectral index is not
  corrected, and the product says so.

The module needs numpy, astropy and katbeam, and no CASA.
"""

import json
import math
import os
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of failing inside a stage.
PIPELINE_INTERFACE_VERSION = 1

METHOD = "katbeam_wsclean_mfs"
PROVIDER = "katbeam"
# Same bar as the CASA path: a beam whose peak is not unity would rescale the
# published flux-density scale.
PEAK_TOLERANCE = (0.95, 1.05)
SAMPLES_PER_BAND = 9
# The header keyword in which WSClean records an image's imaging weight, the
# weight that image carries in the MFS image.
IMAGE_WEIGHT_KEYWORD = "WSCIMGWG"


class BeamCorrectionError(RuntimeError):
    """A MeerKAT beam correction could not be made safely."""


# ------------------------------------------------------------------
# katbeam
# ------------------------------------------------------------------

def load_katbeam():
    """Return ``(JimBeam, version)``; refuse a package with no version."""

    try:
        import katbeam
        beam_class = katbeam.JimBeam
    except Exception as exc:
        raise BeamCorrectionError(
            "MeerKAT primary-beam correction requires the katbeam package in "
            "the Python that runs the WSClean imaging stage (--python). "
            "Install it there (pip install katbeam), or image with "
            "--imager casa."
        ) from exc

    version = getattr(katbeam, "__version__", None)
    if not version:
        try:
            from importlib import metadata
            version = metadata.version("katbeam")
        except Exception:
            version = None
    if not version or str(version).strip().lower() in {"", "unknown", "none"}:
        raise BeamCorrectionError(
            "The installed katbeam package has no discoverable version. "
            "Recording an unversioned beam model would make the published "
            "flux scale unreproducible; install a packaged katbeam release."
        )
    return beam_class, str(version)


# ------------------------------------------------------------------
# FITS geometry
# ------------------------------------------------------------------

def _astropy_fits():
    try:
        from astropy.io import fits
    except Exception as exc:
        raise BeamCorrectionError(
            "MeerKAT primary-beam correction needs astropy in the Python "
            "that runs the WSClean imaging stage (--python)."
        ) from exc
    return fits


def _angle_to_degrees(value, unit):
    unit = str(unit or "deg").strip().lower()
    factors = {"deg": 1.0, "degree": 1.0, "rad": 180.0 / math.pi,
               "arcmin": 1.0 / 60.0, "arcsec": 1.0 / 3600.0}
    if unit not in factors:
        raise BeamCorrectionError(f"Unsupported FITS angle unit {unit!r}.")
    return float(value) * factors[unit]


def _frequency_to_hz(value, unit):
    unit = str(unit or "Hz").strip().lower()
    factors = {"hz": 1.0, "khz": 1e3, "mhz": 1e6, "ghz": 1e9}
    if unit not in factors:
        raise BeamCorrectionError(f"Unsupported FITS frequency unit {unit!r}.")
    return float(value) * factors[unit]


def image_geometry(header):
    """Sky geometry and frequency of a single-plane WSClean image header."""

    naxis = int(header.get("NAXIS", 0))
    if naxis < 2:
        raise BeamCorrectionError("FITS image has fewer than two axes.")

    ctype1 = str(header.get("CTYPE1", "")).upper()
    ctype2 = str(header.get("CTYPE2", "")).upper()
    if not (ctype1.startswith("RA") and ctype2.startswith("DEC")):
        raise BeamCorrectionError(
            "Axes 1 and 2 of the FITS image are not RA and Dec "
            f"({ctype1!r}, {ctype2!r}); only WSClean's layout is supported."
        )
    projection = ctype1.rsplit("-", 1)[-1]
    if projection not in {"SIN", "TAN"}:
        raise BeamCorrectionError(
            f"Unsupported FITS projection {projection!r}; expected SIN or TAN."
        )

    frequency_hz = None
    bandwidth_hz = None
    for axis in range(3, naxis + 1):
        length = int(header.get(f"NAXIS{axis}", 1))
        if length != 1:
            raise BeamCorrectionError(
                f"FITS axis {axis} has {length} planes; each WSClean image "
                "must be a single plane."
            )
        if str(header.get(f"CTYPE{axis}", "")).upper().startswith("FREQ"):
            unit = header.get(f"CUNIT{axis}", "Hz")
            frequency_hz = _frequency_to_hz(header[f"CRVAL{axis}"], unit)
            delta = header.get(f"CDELT{axis}")
            if delta not in (None, 0):
                bandwidth_hz = abs(_frequency_to_hz(delta, unit))
    if frequency_hz is None or frequency_hz <= 0:
        raise BeamCorrectionError(
            "The FITS image records no frequency axis, so its primary beam "
            "cannot be evaluated."
        )

    return {
        "nx": int(header["NAXIS1"]),
        "ny": int(header["NAXIS2"]),
        "crpix_x": float(header["CRPIX1"]),
        "crpix_y": float(header["CRPIX2"]),
        "cdelt_x_deg": _angle_to_degrees(
            header["CDELT1"], header.get("CUNIT1", "deg")
        ),
        "cdelt_y_deg": _angle_to_degrees(
            header["CDELT2"], header.get("CUNIT2", "deg")
        ),
        "frequency_hz": frequency_hz,
        "bandwidth_hz": bandwidth_hz,
        "projection": projection,
    }


# ------------------------------------------------------------------
# The beam on an image grid
# ------------------------------------------------------------------

def beam_plane(beam_object, geometry, frequency_mhz, rows_per_chunk=256):
    """The katbeam Stokes-I beam on an image grid, in bounded row chunks.

    Offsets are measured from the reference pixel with the image's own pixel
    increments, signs included: exactly the grid the CASA backend evaluates,
    so a pixel gets the same beam value from either backend.
    """

    nx, ny = geometry["nx"], geometry["ny"]
    x_deg = (
        np.arange(nx, dtype=float) - (geometry["crpix_x"] - 1.0)
    ) * geometry["cdelt_x_deg"]
    result = np.empty((ny, nx), dtype=np.float32)

    for start in range(0, ny, int(rows_per_chunk)):
        stop = min(ny, start + int(rows_per_chunk))
        y_deg = (
            np.arange(start, stop, dtype=float) - (geometry["crpix_y"] - 1.0)
        ) * geometry["cdelt_y_deg"]
        x_grid, y_grid = np.meshgrid(x_deg, y_deg)
        chunk = np.asarray(
            beam_object.I(x_grid, y_grid, float(frequency_mhz)), dtype=float
        )
        if chunk.shape != x_grid.shape:
            raise BeamCorrectionError(
                "katbeam returned an invalid primary-beam plane."
            )
        if not np.isfinite(chunk).all():
            raise BeamCorrectionError(
                f"katbeam returned {int(np.count_nonzero(~np.isfinite(chunk)))} "
                "non-finite primary-beam pixels; refusing to divide an image "
                "by a beam containing NaN or Inf."
            )
        result[start:stop] = chunk
    return result


def check_peak(beam, label="primary beam"):
    peak = float(np.max(beam))
    if not PEAK_TOLERANCE[0] <= peak <= PEAK_TOLERANCE[1]:
        raise BeamCorrectionError(
            f"The {label} peaks at {peak:.4f} rather than unity; a "
            "non-normalized beam would rescale the published flux-density "
            "scale."
        )
    return peak


# ------------------------------------------------------------------
# Frequencies
# ------------------------------------------------------------------

def usable_intervals(low_hz, high_hz, masked_hz=None):
    """``[low_hz, high_hz]`` less the masked ranges, as ``(low, high)`` pairs."""

    low_hz, high_hz = float(low_hz), float(high_hz)
    pieces = [(low_hz, high_hz)] if high_hz > low_hz else []
    for mask in masked_hz or ():
        mask_low, mask_high = sorted((float(mask[0]), float(mask[1])))
        remaining = []
        for start, stop in pieces:
            if mask_high <= start or mask_low >= stop:
                remaining.append((start, stop))
                continue
            if mask_low > start:
                remaining.append((start, mask_low))
            if mask_high < stop:
                remaining.append((mask_high, stop))
        pieces = remaining
    return [(start, stop) for start, stop in pieces if stop > start]


def effective_band(center_hz, bandwidth_hz, limits_hz, masked_hz=None):
    """The part of a nominal band that still holds data, and its centre.

    WSClean divides the whole selected band into equal output channels and
    labels each by its nominal centre. The flagging stage has removed
    everything outside the receiver's science range and inside the
    telescope's known-RFI masks (``masked_hz``, ``(low, high)`` pairs), so a
    sub-band holds data only where neither applies: the centre of that
    remainder is its effective frequency, and its share of the nominal width
    is its coverage. On MeerKAT L-band, with eight sub-bands, the
    1163-1299 MHz mask empties the 1177-1284 MHz sub-band and the
    1524-1630 MHz mask leaves 1498-1524 MHz of the 1498-1605 MHz one.

    Returns ``(effective_center_hz, covered_fraction)``.
    """

    if not bandwidth_hz or bandwidth_hz <= 0:
        return float(center_hz), 1.0
    low = center_hz - 0.5 * bandwidth_hz
    high = center_hz + 0.5 * bandwidth_hz
    if limits_hz:
        low = max(low, float(limits_hz[0]))
        high = min(high, float(limits_hz[1]))
    pieces = usable_intervals(low, high, masked_hz)
    covered = sum(stop - start for start, stop in pieces)
    if covered <= 0:
        return float(center_hz), 0.0
    centre = sum(
        0.5 * (start + stop) * (stop - start) for start, stop in pieces
    ) / covered
    return centre, covered / float(bandwidth_hz)


def band_samples_mhz(center_hz, bandwidth_hz, limits_hz, count=SAMPLES_PER_BAND,
                     masked_hz=None):
    """Evenly spaced frequencies across the usable part of one wide band.

    With masks, the samples are spaced evenly through the unflagged
    frequencies only, so a masked range contributes nothing.
    """

    if not bandwidth_hz or bandwidth_hz <= 0:
        return [center_hz / 1e6]
    low = center_hz - 0.5 * bandwidth_hz
    high = center_hz + 0.5 * bandwidth_hz
    if limits_hz:
        low = max(low, float(limits_hz[0]))
        high = min(high, float(limits_hz[1]))
    if high <= low:
        return [center_hz / 1e6]
    if not masked_hz:
        return [
            float(value) / 1e6 for value in np.linspace(low, high, int(count))
        ]
    pieces = usable_intervals(low, high, masked_hz)
    total = sum(stop - start for start, stop in pieces)
    if total <= 0:
        return [center_hz / 1e6]
    samples = []
    for index in range(int(count)):
        position = (index + 0.5) / int(count) * total
        for start, stop in pieces:
            if position <= stop - start:
                samples.append((start + position) / 1e6)
                break
            position -= stop - start
    return samples


def wsclean_image_weight(header):
    """The imaging weight WSClean recorded for an image, or None."""

    value = header.get(IMAGE_WEIGHT_KEYWORD)
    if value is None or isinstance(value, bool):
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value) or value < 0:
        return None
    return value


# ------------------------------------------------------------------
# What a single-frequency correction would have done
# ------------------------------------------------------------------

def flat_spectrum_bias(beam_object, reference_mhz, frequencies_mhz, weights,
                       levels, limit_deg=4.0, step_deg=0.002):
    """Flux error a reference-frequency beam leaves on a flat-spectrum source.

    At the radius where the reference-frequency beam falls to each level, the
    MFS image of a flat-spectrum source is attenuated by the weighted mean of
    the beams across the band, not by the reference beam, so dividing by the
    latter is off by ``mean(B) / B_ref - 1``. Reported along the x axis from
    the pointing centre, for the levels given.
    """

    radii = np.arange(0.0, float(limit_deg), float(step_deg))
    # Two-dimensional, like the grids the beam is evaluated on elsewhere.
    row_x = radii[np.newaxis, :]
    row_y = np.zeros_like(row_x)

    def beam_along_x(freq):
        return np.asarray(
            beam_object.I(row_x, row_y, float(freq)), dtype=float
        ).reshape(-1)

    reference = beam_along_x(reference_mhz)
    stack = np.asarray([beam_along_x(freq) for freq in frequencies_mhz])
    if not (np.isfinite(reference).all() and np.isfinite(stack).all()):
        raise BeamCorrectionError(
            "katbeam returned non-finite values along the radial profile."
        )
    weights = np.asarray(weights, dtype=float)
    mean = (stack * weights[:, None]).sum(axis=0) / weights.sum()

    result = {}
    for level in levels:
        below = np.where(reference <= float(level))[0]
        if below.size == 0 or below[0] == 0:
            continue
        index = int(below[0])
        result[f"{float(level):g}"] = {
            "radius_deg": float(radii[index]),
            "reference_beam": float(reference[index]),
            "flux_bias_of_single_frequency_correction": float(
                mean[index] / reference[index] - 1.0
            ),
        }
    return result


# ------------------------------------------------------------------
# FITS products
# ------------------------------------------------------------------

def _write_fits_atomically(path, data, header, comment_lines,
                           dimensionless=False):
    fits = _astropy_fits()
    path = Path(path)
    temporary = path.with_name(f"{path.name}.partial.{os.getpid()}")
    new_header = header.copy()
    if dimensionless:
        # A primary beam has no flux unit and no restoring beam of its own.
        for key in ("BUNIT", "BMAJ", "BMIN", "BPA"):
            new_header.remove(key, ignore_missing=True)
    for line in comment_lines:
        new_header["HISTORY"] = line
    try:
        fits.PrimaryHDU(data=data, header=new_header).writeto(
            str(temporary), overwrite=True
        )
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def correct_one(pb_path, pbcor_path, beam, geometry, *, pbmin, model,
                version, frequency_label, source_data, source_header):
    """Write the beam and the beam-corrected version of one image."""

    shape = source_data.shape
    plane = np.asarray(source_data, dtype=np.float32).reshape(
        geometry["ny"], geometry["nx"]
    )

    valid = beam >= float(pbmin)
    corrected = np.full(plane.shape, np.nan, dtype=np.float32)
    np.divide(plane, beam, out=corrected, where=valid)
    corrected[~np.isfinite(plane)] = np.nan

    note = (
        f"Divided by the katbeam {model} Stokes-I beam ({frequency_label}), "
        f"pixels with beam < {float(pbmin):g} blanked; katbeam {version}."
    )
    _write_fits_atomically(
        pb_path, beam.reshape(shape), source_header,
        [f"katbeam {model} primary beam, {frequency_label}."],
        dimensionless=True,
    )
    _write_fits_atomically(
        pbcor_path, corrected.reshape(shape), source_header, [note],
    )
    return {
        "pixels": int(corrected.size),
        "pixels_inside_pbmin": int(np.count_nonzero(valid)),
        "beam_minimum_kept": (
            float(beam[valid].min()) if valid.any() else None
        ),
    }


def sibling(path, suffix):
    """``name-image.fits`` -> ``name-image<suffix>.fits``."""

    path = Path(path)
    return path.with_name(path.name[:-len(".fits")] + suffix + ".fits")


def pb_name(path):
    """The beam image that goes with a restored image."""

    name = Path(path).name
    if name.endswith("-image.fits"):
        return Path(path).with_name(name[:-len("-image.fits")] + "-katbeam-pb.fits")
    return sibling(path, "-katbeam-pb")


def pbcor_name(path):
    """WSClean's own convention for a primary-beam-corrected image."""

    return sibling(path, "-pb")


def remove_stale_products(path):
    """Remove the beam and corrected image an earlier correction made for ``path``.

    A sub-band that is now skipped gets no new products, so any an earlier
    correction left would sit beside a report that does not describe them.
    Only the two names this module writes for ``path`` are touched.
    """

    removed = []
    for candidate in (pb_name(path), pbcor_name(path)):
        if candidate.is_file():
            resolved = str(candidate.resolve())
            candidate.unlink()
            removed.append(resolved)
    return removed


def correct_products(mfs_image, subband_images, *, model, pbmin,
                     science_limits_hz, identity, mode="band_average",
                     masked_ranges_hz=None, rows_per_chunk=256):
    """Beam-correct a WSClean MFS image and its sub-band images.

    ``masked_ranges_hz`` are the frequency ranges the flagging stage removed,
    as ``(low, high)`` pairs in Hz; each sub-band is corrected at the centre
    of what survived them and the science range. ``identity`` is a callable
    returning a content identity for a file (path, size, SHA-256), used so
    that the provenance written here can be checked against the files later.
    Returns the record the imaging report keeps.
    """

    if mode not in {"band_average", "reference_frequency"}:
        raise BeamCorrectionError(
            f"Unknown MFS beam mode {mode!r}; use band_average or "
            "reference_frequency."
        )
    beam_class, version = load_katbeam()
    fits = _astropy_fits()
    beam_object = beam_class(model)
    limits = tuple(science_limits_hz) if science_limits_hz else None
    masks = sorted(
        sorted((float(low), float(high)))
        for low, high in (masked_ranges_hz or ())
    )

    def load(path):
        with fits.open(str(path), memmap=False) as handle:
            header = handle[0].header.copy()
            data = np.array(handle[0].data, dtype=np.float32)
        return data, header, image_geometry(header)

    paths = sorted(subband_images)
    # WSClean's MFS image is the mean of its sub-band images weighted by the
    # imaging weight each header records, so the MFS image of a flat-spectrum
    # source is attenuated by the same weighted mean of their beams. Those
    # weights are used when every sub-band has one; otherwise each sub-band
    # counts by its unflagged share of the band.
    image_weights = []
    for path in paths:
        with fits.open(str(path), memmap=False) as handle:
            image_weights.append(wsclean_image_weight(handle[0].header))
    use_image_weights = bool(paths) and all(
        weight is not None for weight in image_weights
    )
    weighting = (
        "wsclean_image_weight" if use_image_weights else "covered_fraction"
    )

    records = []
    beam_sum = None
    weight_sum = 0.0
    effective = []

    # Sub-bands first: each at its own effective frequency, and together they
    # give the band-averaged beam the MFS image is divided by.
    for path, image_weight in zip(paths, image_weights):
        data, header, geometry = load(path)
        center_hz = geometry["frequency_hz"]
        eff_hz, fraction = effective_band(
            center_hz, geometry["bandwidth_hz"], limits, masks
        )
        base = {
            "source_image": str(Path(path).resolve()),
            "nominal_frequency_mhz": center_hz / 1e6,
            "frequency_mhz": eff_hz / 1e6,
            "covered_fraction": fraction,
            "image_weight": image_weight,
        }
        reason = None
        if fraction <= 0.0:
            inside = effective_band(
                center_hz, geometry["bandwidth_hz"], limits
            )[1] > 0.0
            reason = (
                "every frequency in it is masked by the flagging stage"
                if inside else "outside the science frequency range"
            )
        elif use_image_weights and image_weight <= 0.0:
            reason = (
                "WSClean gave it no imaging weight: none of its data are "
                "unflagged"
            )
        elif not np.isfinite(data).any():
            reason = "the image holds no finite pixels"
        if reason:
            records.append(dict(
                base, skipped=reason,
                removed_stale_products=remove_stale_products(path),
            ))
            continue

        weight = image_weight if use_image_weights else fraction
        beam = beam_plane(beam_object, geometry, eff_hz / 1e6, rows_per_chunk)
        check_peak(beam, f"{eff_hz / 1e6:.1f} MHz primary beam")
        pb_path, pbcor_path = pb_name(path), pbcor_name(path)
        stats = correct_one(
            pb_path, pbcor_path, beam, geometry, pbmin=pbmin,
            model=model, version=version,
            frequency_label=(
                f"{eff_hz / 1e6:.2f} MHz, the centre of the unflagged part "
                f"of this {center_hz / 1e6:.2f} MHz sub-band"
            ),
            source_data=data, source_header=header,
        )
        if beam_sum is None:
            beam_sum = np.zeros(beam.shape, dtype=np.float64)
        beam_sum += beam.astype(np.float64) * weight
        weight_sum += weight
        effective.append((eff_hz / 1e6, weight))
        records.append(dict(
            stats,
            **base,
            pb_image=str(pb_path.resolve()),
            pbcor_image=str(pbcor_path.resolve()),
            mfs_weight=weight,
        ))

    data, header, geometry = load(mfs_image)
    reference_mhz = geometry["frequency_hz"] / 1e6

    # The frequencies the data actually span: the sub-bands' effective
    # frequencies when there are sub-bands, else samples across the usable
    # band. They describe the data, so they are the same whichever way the
    # MFS beam is built, and they are what a single-frequency beam's error
    # is measured against.
    if effective:
        frequencies = [freq for freq, _ in effective]
        weights = [weight for _, weight in effective]
    else:
        frequencies = band_samples_mhz(
            geometry["frequency_hz"], geometry["bandwidth_hz"], limits,
            masked_hz=masks,
        )
        weights = [1.0] * len(frequencies)
        weighting = "uniform_band_samples"
    # The frequency the MFS image describes. WSClean labels it with the
    # centre of the whole band, but the weights put it nearer the top: on
    # MeerKAT L-band 1708779077 it was 1400 MHz, not 1284 MHz.
    mfs_effective_mhz = float(
        sum(freq * weight for freq, weight in zip(frequencies, weights)) /
        float(sum(weights))
    )

    if mode == "reference_frequency":
        beam = beam_plane(beam_object, geometry, reference_mhz, rows_per_chunk)
        label = f"{reference_mhz:.2f} MHz"
    elif beam_sum is not None and weight_sum > 0:
        beam = (beam_sum / weight_sum).astype(np.float32)
        label = (
            f"band average of {len(frequencies)} sub-bands, "
            f"{min(frequencies):.1f}-{max(frequencies):.1f} MHz, weighted by "
            + (
                "WSClean's imaging weights"
                if use_image_weights else "their unflagged coverage"
            )
            + f"; effective frequency {mfs_effective_mhz:.1f} MHz"
        )
    else:
        accumulated = np.zeros((geometry["ny"], geometry["nx"]),
                               dtype=np.float64)
        for frequency in frequencies:
            accumulated += beam_plane(
                beam_object, geometry, frequency, rows_per_chunk
            )
        beam = (accumulated / len(frequencies)).astype(np.float32)
        label = (
            f"band average of {len(frequencies)} frequencies, "
            f"{min(frequencies):.1f}-{max(frequencies):.1f} MHz"
        )
    used_mode = mode

    check_peak(beam, "MFS primary beam")
    mfs_pb, mfs_pbcor = pb_name(mfs_image), pbcor_name(mfs_image)
    stats = correct_one(
        mfs_pb, mfs_pbcor, beam, geometry, pbmin=pbmin,
        model=model, version=version, frequency_label=label,
        source_data=data, source_header=header,
    )

    # Informational. The image does not depend on it, so a problem here is
    # recorded and must never stop the correction that has just succeeded.
    try:
        bias = flat_spectrum_bias(
            beam_object, reference_mhz, frequencies, weights,
            levels=sorted({float(pbmin), 0.5}, reverse=True),
        )
    except Exception as exc:
        bias = {"error": f"{type(exc).__name__}: {exc}"}

    record = {
        # 2: the MFS beam follows WSClean's imaging weights and each sub-band
        # is evaluated where its unflagged data lie (finding 165).
        "schema_version": 2,
        "provider": PROVIDER,
        "method": METHOD,
        "model": model,
        "katbeam_version": version,
        "pbmin": float(pbmin),
        "mfs_beam_mode": used_mode,
        "spectral_terms_corrected": False,
        "frequency_mhz": reference_mhz,
        "effective_frequencies_mhz": [float(value) for value in frequencies],
        "frequency_weights": [float(value) for value in weights],
        "mfs_weighting": weighting,
        "mfs_effective_frequency_mhz": mfs_effective_mhz,
        "masked_ranges_hz": [[low, high] for low, high in masks],
        "science_range_hz": list(limits) if limits else None,
        "source_image": str(Path(mfs_image).resolve()),
        "pb_image": str(mfs_pb.resolve()),
        "pbcor_image": str(mfs_pbcor.resolve()),
        "geometry": {key: geometry[key] for key in (
            "nx", "ny", "crpix_x", "crpix_y", "cdelt_x_deg", "cdelt_y_deg",
            "projection",
        )},
        "statistics": stats,
        "single_frequency_flux_bias": bias,
        "subbands": records,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    # Identities last, so they describe the files as they are on disk.
    for key in ("source_image", "pb_image", "pbcor_image"):
        record[key.replace("_image", "")] = identity(record[key])
    for item in records:
        if "pbcor_image" in item:
            for key in ("source_image", "pb_image", "pbcor_image"):
                item[key.replace("_image", "")] = identity(item[key])
    return record


def write_provenance(path, record):
    path = Path(path)
    temporary = path.with_name(f"{path.name}.partial.{os.getpid()}")
    try:
        temporary.write_text(json.dumps(record, indent=2, sort_keys=True))
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


_SUBBAND_SUFFIX = re.compile(r"^-(\d{4})-image\.fits$")


def classify_outputs(paths, prefix_name):
    """Split WSClean's restored images into the MFS image and sub-bands.

    Only the part of each name after ``prefix_name`` is examined, so a target
    whose own name contains "model", "beam" or a channel-like number cannot
    be mistaken for a product type.
    """

    mfs = None
    plain = None
    subbands = []
    for path in paths:
        name = Path(path).name
        if not name.startswith(prefix_name):
            continue
        suffix = name[len(prefix_name):]
        if suffix == "-MFS-image.fits":
            mfs = Path(path)
        elif suffix == "-image.fits":
            plain = Path(path)
        elif _SUBBAND_SUFFIX.match(suffix):
            subbands.append(Path(path))
    return (mfs or plain), sorted(subbands)
