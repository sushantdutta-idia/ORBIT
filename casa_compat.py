#!/usr/bin/env python3
"""Small compatibility layer for CASA 6 task keyword drift."""

import functools
import importlib
import inspect
import math
import re


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of failing inside a stage.
PIPELINE_INTERFACE_VERSION = 2


def _task_keyword_set(task):

    try:
        signature = inspect.signature(task)
    except (TypeError, ValueError):
        return None

    keywords = set()

    for parameter in signature.parameters.values():

        if parameter.kind == parameter.VAR_KEYWORD:
            return None

        if parameter.kind in {
            parameter.POSITIONAL_OR_KEYWORD,
            parameter.KEYWORD_ONLY,
        }:
            keywords.add(parameter.name)

    return keywords


def _warn_dropped(task_name, names):

    if not names:
        return

    print(
        "CASA compatibility: "
        f"{task_name} does not support "
        f"{', '.join(sorted(names))}; dropping."
    )


def _extract_unsupported_keywords(message, current_keys):

    patterns = [
        r"unexpected keyword argument ['\"]([^'\"]+)['\"]",
        r"unexpected keyword ['\"]([^'\"]+)['\"]",
        r"unknown keyword ['\"]([^'\"]+)['\"]",
        r"unknown keyword\s+([A-Za-z_][A-Za-z0-9_]*)",
        r"unknown parameter ['\"]([^'\"]+)['\"]",
        r"unknown parameter\s+([A-Za-z_][A-Za-z0-9_]*)",
        r"unrecognized parameter ['\"]([^'\"]+)['\"]",
        r"unrecognized parameter\s+([A-Za-z_][A-Za-z0-9_]*)",
        r"parameter ['\"]([^'\"]+)['\"] is not",
        r"parameter ([A-Za-z_][A-Za-z0-9_]*) is not",
    ]

    matches = set()

    for pattern in patterns:
        matches.update(re.findall(pattern, message, flags=re.I))

    if matches:

        normalized_matches = {
            match.lower()
            for match in matches
        }

        return {
            key for key in current_keys
            if key.lower() in normalized_matches
        }

    # No pattern matched. The previous fallback returned every supplied
    # keyword whose name appeared anywhere in the message whenever the text
    # contained a token such as "invalid parameter". A genuine CASA failure
    # like "Invalid parameter: field 'XYZ' does not exist; check vis, field,
    # spw" therefore looked like a signature mismatch on vis/field/spw, and
    # the real diagnosis was replaced by a misleading "does not support
    # required science keyword" error. Attribute nothing on that basis.
    return set()


def compat_task(task, name=None, optional_keywords=()):
    """Wrap a CASA task without silently weakening science parameters.

    CASA releases occasionally remove cosmetic keywords.  Callers may name
    those explicitly in *optional_keywords*.  Every other supplied keyword is
    treated as part of the scientific recipe and an unsupported keyword is a
    hard error rather than being dropped silently.
    """

    task_name = name or getattr(task, "__name__", "CASA task")
    supported_keywords = _task_keyword_set(task)
    optional_keywords = set(optional_keywords or ())

    def reject_required(names, cause=None):
        required = sorted(set(names).difference(optional_keywords))
        if required:
            error = RuntimeError(
                f"CASA task {task_name!r} does not support required science "
                f"keyword(s): {', '.join(required)}. Refusing to run a "
                "weakened recipe; use a supported CASA 6 release or update "
                "the telescope recipe explicitly."
            )
            raise error from cause

    @functools.wraps(task)
    def wrapper(*args, **kwargs):

        call_kwargs = dict(kwargs)

        if supported_keywords is not None:
            dropped = [
                key for key in list(call_kwargs)
                if key not in supported_keywords
            ]

            reject_required(dropped)

            for key in dropped:
                call_kwargs.pop(key, None)

            _warn_dropped(task_name, dropped)

        while True:

            try:
                return task(*args, **call_kwargs)

            except Exception as exc:
                # Only a Python signature error can be attributed to an
                # unsupported keyword with confidence. Anything else is a
                # real task failure and must propagate with its own message:
                # re-invoking a CASA task after a partial failure could apply
                # flags or write a calibration table twice.
                if not isinstance(exc, TypeError):
                    raise

                bad_keys = _extract_unsupported_keywords(
                    str(exc),
                    call_kwargs,
                )

                if not bad_keys:
                    raise

                reject_required(bad_keys, cause=exc)

                for key in bad_keys:
                    call_kwargs.pop(key, None)

                _warn_dropped(task_name, bad_keys)

    wrapper.casa_task_available = True
    return wrapper


def _missing_task(name):

    def wrapper(*args, **kwargs):
        raise RuntimeError(
            f"CASA task '{name}' is not available in this CASA 6 "
            "installation."
        )

    wrapper.__name__ = name
    wrapper.casa_task_available = False
    return wrapper


def import_casa_task(name, required=True, optional_keywords=()):

    try:
        casatasks = importlib.import_module("casatasks")
        task = getattr(casatasks, name)
    except Exception as exc:
        if required:
            raise ImportError(
                f"Required CASA task '{name}' is not available."
            ) from exc

        return _missing_task(name)

    return compat_task(
        task,
        name=name,
        optional_keywords=optional_keywords,
    )


def task_available(task):

    return bool(
        getattr(task, "casa_task_available", False)
    )


def spectral_window_ids(msmd):
    """Return MS spectral-window IDs across CASA 6 releases.

    ``msmetadata.nspw()`` is the CASA 6 public API. Some local wrappers
    expose a non-standard ``spws()`` helper, so retain it only as a fallback.
    Spectral-window table rows are zero-based and contiguous by definition.
    """

    nspw = getattr(msmd, "nspw", None)

    if callable(nspw):
        try:
            count = int(nspw())
        except Exception as exc:
            raise RuntimeError(
                "Could not determine the number of spectral windows "
                "with msmetadata.nspw()."
            ) from exc

        if count < 0:
            raise RuntimeError(
                "msmetadata.nspw() returned an invalid negative value."
            )

        return list(range(count))

    spws = getattr(msmd, "spws", None)

    if callable(spws):
        return [int(spw) for spw in spws()]

    raise RuntimeError(
        "This CASA msmetadata tool exposes neither nspw() nor spws()."
    )


# tclean's stop codes, as CASA 6 documents them in its return dictionary.
TCLEAN_STOP_REASONS = {
    0: "did not run a minor cycle",
    1: "reached the iteration limit (niter)",
    2: "reached the threshold",
    3: "was stopped by force",
    4: "saw no change in the peak residual across two major cycles",
    5: "saw the peak residual rise by more than 3x in one major cycle",
    6: "saw the peak residual rise by more than 3x from its minimum",
    7: "found the mask empty",
    8: "reached the nsigma threshold",
    9: "reached the major-cycle limit (nmajor)",
}

# Stop codes after which the deconvolution can be said to have finished:
# the residual fell to the requested depth, or there was nothing masked left
# to clean. Every other code means CLEAN stopped with work outstanding.
TCLEAN_CONVERGED_STOPCODES = {2, 7, 8}


def summarize_tclean(result):
    """The parts of tclean's return value that say why it stopped.

    tclean returns a dictionary when it runs minor cycles. Older releases, a
    PSF-only call and a test stub return nothing; those give None rather than
    a guess.
    """

    if not isinstance(result, dict):
        return None

    def number(key, cast):
        try:
            value = result.get(key)
            return None if value is None else cast(value)
        except (TypeError, ValueError):
            return None

    stopcode = number("stopcode", int)
    summary = {
        "stopcode": stopcode,
        "stop_reason": TCLEAN_STOP_REASONS.get(
            stopcode,
            None if stopcode is None else f"stop code {stopcode}",
        ),
        "converged": (
            None if stopcode is None
            else stopcode in TCLEAN_CONVERGED_STOPCODES
        ),
        "iterdone": number("iterdone", int),
        "niter": number("niter", int),
        "nmajordone": number("nmajordone", int),
        "max_psf_sidelobe": number("maxpsfsidelobe", float),
    }
    description = result.get("stopDescription")
    if description:
        summary["casa_stop_description"] = str(description)
    return summary


def peak_and_noise(imstat, imagename, **selection):
    """An image's true extremes and its robust noise, from two imstat calls.

    The clipping algorithms (chauvenet, fit-half, ...) exclude outliers
    before computing *every* statistic they return, the extremes and their
    positions included. That is what makes their rms a good noise estimate,
    and it is also why their ``max`` can never be a source: CASA reports a
    60-sigma pixel as a 4-sigma maximum somewhere else. So the extremes come
    from the classic algorithm, which sees every pixel, and the noise from
    Chauvenet clipping. ``selection`` (box, region, mask, ...) applies to
    both, so they describe the same pixels.
    """

    everything = imstat(imagename=str(imagename), **selection)
    clipped = imstat(
        imagename=str(imagename),
        algorithm="chauvenet",
        maxiter=5,
        **selection,
    )

    def first(stats, key):
        try:
            return float(stats[key][0])
        except (KeyError, IndexError, TypeError, ValueError):
            return None

    def position(stats, key):
        try:
            return [int(value) for value in list(stats[key])[:2]]
        except (KeyError, TypeError, ValueError):
            return None

    rms = first(clipped, "rms")
    maximum = first(everything, "max")
    minimum = first(everything, "min")
    extreme = None
    extreme_position = None
    if maximum is not None or minimum is not None:
        if minimum is not None and (
            maximum is None or abs(minimum) > abs(maximum)
        ):
            extreme, extreme_position = abs(minimum), position(
                everything, "minpos"
            )
        else:
            extreme, extreme_position = abs(maximum), position(
                everything, "maxpos"
            )

    return {
        "max": maximum,
        "min": minimum,
        "max_position": position(everything, "maxpos"),
        "min_position": position(everything, "minpos"),
        "peak_abs": extreme,
        "peak_abs_position": extreme_position,
        "rms": rms,
        "dynamic_range": (
            maximum / rms if maximum is not None and rms and rms > 0 else None
        ),
        "peak_abs_snr": (
            extreme / rms if extreme is not None and rms and rms > 0 else None
        ),
        "sum": first(everything, "sum"),
        "sumsq": first(everything, "sumsq"),
        "npts": first(everything, "npts"),
        # What noise alone would reach somewhere in these pixels: the floor
        # below which a "significant" peak is not significant.
        "noise_only_peak_snr": noise_peak_snr(first(everything, "npts")),
    }


def bootstrap_pbmask(configured, pblimit, default=0.2):
    """The primary-beam level the bootstrap clean may reach.

    The bootstrap exists to take bright sources and their sidelobes out
    before auto-multithresh measures the noise. A bright source in the outer
    primary beam does the most damage, so by default the bootstrap reaches as
    far as the deconvolution itself does, which is set by ``pblimit``. An
    explicit ``clean_bootstrap_pbmask`` still wins.
    """

    if configured is not None:
        return float(configured)
    try:
        level = abs(float(pblimit))
    except (TypeError, ValueError):
        return float(default)
    return level if 0.0 < level < 1.0 else float(default)


def adaptive_sidelobethreshold(configured, max_psf_sidelobe,
                                max_fraction=0.6, minimum=1.0):
    """The auto-multithresh sidelobe threshold this PSF can actually use.

    auto-multithresh only masks emission above ``sidelobethreshold *
    max_psf_sidelobe * peak_residual``. For a full-synthesis PSF with
    sidelobes near 0.1 the usual 2.5 puts that at a quarter of the peak. A
    short observation has sidelobes of 0.4 or more, and the same 2.5 puts it
    above the peak itself: no pixel can ever qualify, CASA reports "zero
    mask", and nothing is cleaned however bright the field is.

    The threshold is therefore capped so that the criterion sits at no more
    than ``max_fraction`` of the residual peak. A PSF with low sidelobes is
    unaffected. Returns the threshold to use and whether it was lowered.
    """

    try:
        value = float(configured)
    except (TypeError, ValueError):
        return configured, False
    try:
        sidelobe = float(max_psf_sidelobe)
        fraction = float(max_fraction)
    except (TypeError, ValueError):
        return value, False
    if (
        not math.isfinite(sidelobe) or sidelobe <= 0.0 or
        not math.isfinite(fraction) or fraction <= 0.0
    ):
        return value, False

    capped = max(float(minimum), fraction / sidelobe)
    if capped < value:
        return capped, True
    return value, False


def is_zero_mask_stop(summary):
    """Whether a tclean pass ended because its mask held no pixels."""

    if not isinstance(summary, dict):
        return False
    if summary.get("stopcode") == 7:
        return True
    return "zero mask" in str(summary.get("casa_stop_description", "")).lower()


def noise_peak_snr(npixels, pixels_per_beam=10.0, false_alarm=0.01):
    """The |peak|/sigma that noise alone exceeds anywhere with probability
    ``false_alarm``.

    An image of N independent beams of Gaussian noise has its largest
    absolute value above t with probability about 2 N Q(t), where Q is the
    Gaussian upper tail. Any "significant residual" or "strong feature" test
    on the true image maximum must sit above that level, or it fires on
    noise: a 2048-pixel image holds about 400,000 beams at 10 pixels each,
    and noise alone reaches 5 sigma in most of them. ``pixels_per_beam`` is
    the beam area in pixels; 10 corresponds to about 3 pixels across it.
    """

    try:
        beams = float(npixels) / max(1.0, float(pixels_per_beam))
    except (TypeError, ValueError):
        return None
    if beams < 1:
        return None
    target = float(false_alarm) / (2.0 * beams)

    def tail(value):
        return 0.5 * math.erfc(value / math.sqrt(2.0))

    low, high = 0.0, 40.0
    for _ in range(100):
        middle = 0.5 * (low + high)
        if tail(middle) > target:
            low = middle
        else:
            high = middle
    return high
