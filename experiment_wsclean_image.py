#!/usr/bin/env python3
"""Image one self-calibrated uGMRT target with WSClean, as a test.

The pipeline images uGMRT with CASA's tclean (w-projection, two Taylor
terms). This makes the same field from the same self-calibrated data with
WSClean and measures the result the way ``diagnose_noise_floor.py`` measures
the tclean image, so the two compare directly. What differs is the imager:

* the band is deconvolved in several sub-bands jointly, with a polynomial
  through them, so the model can follow a source far off axis whose apparent
  spectrum is the primary beam's and not a power law two Taylor terms fit;
* the mask and the stopping depth follow a local noise estimate;
* the w-term is handled by WGridder.

Two WSClean runs are made. A dirty Stokes V image first: with the same
weighting it is the noise of the data as WSClean weights them, and it takes
a few minutes. Then the Stokes I image. The report gives both and their
ratio. No primary-beam correction is applied: like the pipeline's
self-calibration images this is a flat-noise image, which is what a noise
comparison needs.

The Measurement Set is not changed: WSClean is told not to write a model
column, and its scratch files go to the work directory.

``--time-split N`` asks a different question of the same data. It makes N
dirty images from N consecutive stretches of the observation and measures
the brightest sources in each. A point source has one flux; if every source
is brighter in one stretch than in another, the amplitude scale moved while
the data were taken, and an image of the whole shows it as the same pattern
around every source, on axis or off. It takes minutes.

``--core-model MJY`` prepares a self-calibration step instead of an image.
It cleans only what is brighter than MJY, with positive components, and
lets WSClean write that model into the Measurement Set's MODEL_DATA. It is
for the case where every source carries the same lopsided pattern: a model
cleaned deep contains the pattern and gains solved against it keep it, while
a model of the cores alone says sources are compact and lets the solve
remove it. This is the one mode that writes to the Measurement Set, so it
refuses one that lies inside the pipeline's products; give it a copy.

Once the phases are right, amplitudes want the whole sky as their model, not
its cores: ``--no-negative --write-model`` on the ordinary image run cleans
everything with positive components and leaves that model in MODEL_DATA. The
image it makes is also the honest one, since whatever negative the data hold
beside a source stays in the residual and is not hidden in the model.

``--uv-split`` does the same by baseline length (edges in kilo-wavelengths).
A point source also has one flux on every baseline, and no source has more
on long baselines than on short ones. If every source is brighter on the
long baselines, the antennas of the central square and of the arms are on
different amplitude scales. If the sources are fainter there, by amounts
that differ from one to the next, they are resolved: that is their
structure, and it says the field needs a multi-scale clean, not that the
gains are wrong.

This script is what is run, with any Python 3; it starts WSClean itself.
Where WSClean lives in a container, give the whole command that reaches it
as ``--wsclean``. CASA is used only for what can also be given by hand: the
grid (``--imsize``, ``--cell``), the channels that hold data
(``--channel-range``) and the tclean image to compare with. Where CASA is not
importable those are taken from the command line and the comparison is left
out. The images are measured with numpy alone; without numpy they are made
and ``--measure-only`` measures them later from a Python that has it.

Examples::

    python3 experiment_wsclean_image.py --outdir uGMRT_Products \\
        --target G0212-0335
    python3 experiment_wsclean_image.py --outdir uGMRT_Products \\
        --target G0212-0335 \\
        --wsclean "singularity exec /software/astro/containers/oxkat-latest.sif wsclean"
"""

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# This folder holds the pipeline's own inspect.py. Run from here, anything
# that imports the standard library's inspect (numpy and dataclasses do)
# would get that file instead.
if sys.path and os.path.abspath(sys.path[0]) == SCRIPT_DIR:
    sys.path.pop(0)

import argparse
import contextlib
import importlib.util
import json
import re
import shlex
import shutil
import socket
import subprocess
import time
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def load_diagnostics():
    spec = importlib.util.spec_from_file_location(
        "experiment_wsclean_diagnostics",
        SCRIPT_DIR / "diagnose_noise_floor.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DIAGNOSTICS = load_diagnostics()


# ---------------------------------------------------------------------------
# Channels
# ---------------------------------------------------------------------------

def unflagged_channel_range(per_channel, spw=0, fully_flagged=0.99):
    """First channel with data and one past the last, from a flag summary.

    ``per_channel`` is the ``spw:channel`` part of flagdata's summary. A band
    whose edge is flagged throughout (uGMRT Band 4 is recorded to 950 MHz and
    is usable to 850) would otherwise give WSClean sub-bands with nothing in
    them. Returns None when no channel has data or the summary is unusable.
    """

    channels = []
    for key, counts in (per_channel or {}).items():
        try:
            key_spw, channel = (int(part) for part in str(key).split(":"))
            total = float(counts["total"])
            flagged = float(counts["flagged"])
        except (KeyError, TypeError, ValueError):
            continue
        if key_spw == int(spw) and total > 0 and flagged / total < fully_flagged:
            channels.append(channel)
    if not channels:
        return None
    return min(channels), max(channels) + 1


def measured_channel_range(ms):
    """The channel range that holds data, read from the Measurement Set."""

    from casatasks import flagdata

    summary = flagdata(
        vis=str(ms),
        mode="summary",
        spwchan=True,
        action="calculate",
    )
    if not isinstance(summary, dict):
        return None
    return unflagged_channel_range(summary.get("spw:channel"))


# ---------------------------------------------------------------------------
# WSClean
# ---------------------------------------------------------------------------

def wsclean_capabilities(prefix):
    """The version of a WSClean and the options its help lists.

    Option names have changed between releases (``-use-wgridder`` became
    ``-gridder wgridder``), and a rejected option stops WSClean only after
    it has been launched. Asking first costs two seconds.
    """

    def output(arguments):
        try:
            done = subprocess.run(
                list(prefix) + arguments,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=180,
            )
        except (OSError, subprocess.SubprocessError):
            return ""
        return done.stdout or ""

    found = re.search(r"WSClean version (\d+)\.(\d+)", output(["--version"]))
    options = set(re.findall(r"(?m)^(-[a-z][a-z0-9-]+)", output([])))
    return {
        "version": (int(found.group(1)), int(found.group(2))) if found else None,
        # A help text that lists almost nothing was not a help text.
        "options": options if len(options) >= 20 else set(),
    }


# The faintest a cores-only model goes, in units of the image's noise.
CORE_MODEL_MIN_SIGMA = 5.0


def gridder_arguments(capabilities):
    """How this WSClean is told to use WGridder."""

    options = (capabilities or {}).get("options") or set()
    version = (capabilities or {}).get("version")
    if "-gridder" in options:
        return ["-gridder", "wgridder"]
    if "-use-wgridder" in options:
        return ["-use-wgridder"]
    if version is not None and version < (3, 2):
        return ["-use-wgridder"]
    return ["-gridder", "wgridder"]


def parallel_deconvolution_size(args, capabilities):
    """Sub-image size for parallel deconvolution, or 0 to leave it off.

    WSClean 3.1.0 aborts with "double free or corruption" as multi-scale
    deconvolution starts when the image is split for parallel cleaning
    (uGMRT Band 4, 5120 pixels, six joined sub-bands, local rms). Releases
    before 3.2 therefore clean the image whole unless a size is asked for.
    """

    if args.parallel_deconvolution is not None:
        return max(0, int(args.parallel_deconvolution))
    version = (capabilities or {}).get("version")
    if version is not None and version < (3, 2):
        return 0
    return 2048


def supported(capabilities, option):
    """Whether an option may be passed: yes unless the help rules it out."""

    options = (capabilities or {}).get("options") or set()
    return not options or option in options


FEATURES = ("multiscale", "local_rms", "spectral_fit")


def feature_ladder(args):
    """Deconvolution set-ups to try in turn, fullest first.

    WSClean 3.1.0 in the container on ilifu aborts with a corrupted heap as
    deconvolution starts, with multi-scale, the local-rms mask and the
    spectral fit all on, whether or not the image is split for parallel
    cleaning. Which of them it cannot do is found by running: each rung
    drops one, and the report says which rung made the image.
    """

    rung = {
        "multiscale": bool(args.scales),
        "local_rms": not args.no_local_rms,
        "spectral_fit": int(args.spectral_terms) > 0,
    }
    ladder = [dict(rung)]
    for feature in FEATURES:
        if rung[feature]:
            rung = dict(rung, **{feature: False})
            ladder.append(rung)
    return ladder


def describe(features):
    names = {
        "multiscale": "multi-scale",
        "local_rms": "local-rms mask",
        "spectral_fit": "spectral fit",
    }
    if features.get("scales"):
        names["multiscale"] = f"multi-scale to {max(features['scales'])} px"
    kept = [names[key] for key in FEATURES if features[key]]
    return ", ".join(kept) if kept else "plain joined-channel clean"


# Fewest scales a multi-scale clean is left with after one diverged.
FEWEST_SCALES = 3


def calmer(features, scales):
    """The deconvolution to try after one that diverged, or None.

    G0213-0251's robust-0.5 image (Band 4, one 22-minute scan, 2026-10-29)
    diverged from its fourth major iteration: the 30-pixel scale took
    components of tens of Jy, the residual's noise went from 59 uJy/beam to
    19 mJy/beam, and the image came out 2200 times the noise of the data.
    Its robust-0 image of the same data, with the same scales, did not. The
    largest scale goes first; with fewer than three scales left,
    multi-scale goes.
    """

    if not features["multiscale"]:
        return None
    kept = sorted(features.get("scales") or scales)[:-1]
    gentler = {key: value for key, value in features.items()
               if key != "scales"}
    if len(kept) < FEWEST_SCALES:
        gentler["multiscale"] = False
    else:
        gentler["scales"] = kept
    return gentler


def correlations_of(outdir):
    """The correlations the pipeline's inspection found in the data, or None."""

    try:
        with open(Path(outdir) / "inspection" / "inspection_report.json") as handle:
            report = json.load(handle)
    except (OSError, ValueError, TypeError):
        return None
    value = report.get("correlations") if isinstance(report, dict) else None
    if isinstance(value, dict):
        value = [
            item
            for items in (value.get("setups") or {}).values()
            for item in (items or [])
        ]
    if not isinstance(value, (list, tuple)):
        return None
    return {str(item).upper() for item in value} or None


def noise_stokes(args):
    """The Stokes parameter of the dirty image that gives the data's noise.

    Stokes V holds almost none of the sky, so its dirty image is the noise
    of the data under the image's weighting. WSClean forms it from the
    circular hands or from all four linear correlations. MeerKAT's data hold
    XX and YY only, from which it forms I and Q alone, and Q = (XX - YY) / 2
    holds what the two hands' beams see differently of every bright source
    off axis: on 1708779077 that is J2341-3506, 0.8 Jy 35 arcmin out. For
    such data no noise image is made unless ``--noise-stokes`` asks for one.
    """

    choice = str(getattr(args, "noise_stokes", None) or "auto").strip().upper()
    if choice == "NONE":
        return None
    if choice != "AUTO":
        return choice
    found = correlations_of(getattr(args, "outdir", None) or "")
    if not found:
        return "V"
    if found & {"RR", "LL", "RL", "LR"} or {"XY", "YX"} <= found:
        return "V"
    return None


def image_label(robust, taper=0):
    """As the images are named: robust_0, robust_0p5, robust_0_taper15."""

    label = f"robust_{float(robust):g}"
    if taper and float(taper) > 0:
        label += f"_taper{float(taper):g}"
    return label.replace("-", "m").replace(".", "p")


def wsclean_command(args, ms, name, imsize, cell, channel_range, stokes,
                    temp_dir, capabilities=None, features=None):
    """The WSClean command for one polarization."""

    if features is None:
        features = feature_ladder(args)[0]

    writes_model = bool(args.write_model) and stokes.upper() == "I"
    command = shlex.split(args.wsclean) + [
        "-name", str(name),
        "-data-column", args.data_column,
    ]
    if not writes_model:
        # The Measurement Set is a pipeline product: no MODEL_DATA from
        # here, unless a copy was given and its model asked for.
        command.append("-no-update-model-required")
    command += [
        "-temp-dir", str(temp_dir),
        "-j", str(int(args.threads)),
        "-abs-mem", f"{float(args.memory_gb):g}",
        *gridder_arguments(capabilities),
        "-size", str(int(imsize)), str(int(imsize)),
        "-scale", f"{cell:.4f}asec",
        "-weight", "briggs", f"{float(args.robust):g}",
    ]
    taper = float(getattr(args, "taper", 0) or 0)
    if taper > 0:
        # Extended emission of low surface brightness needs a larger
        # beam: the long baselines are weighted down by a Gaussian.
        command += ["-taper-gaussian", f"{taper:g}asec"]
    command += [
        "-pol", stokes,
        "-join-channels",
        "-channels-out", str(int(args.channels_out)),
    ]
    if channel_range:
        command += [
            "-channel-range", str(int(channel_range[0])),
            str(int(channel_range[1])),
        ]
    if stokes.upper() != "I":
        # The noise of the data under this weighting; nothing to clean.
        return command + ["-niter", "0", str(ms)]

    command += [
        "-niter", str(int(args.niter)),
        "-nmiter", str(int(args.nmiter)),
        "-mgain", f"{float(args.mgain):g}",
        "-auto-mask", f"{float(args.auto_mask):g}",
        "-auto-threshold", f"{float(args.auto_threshold):g}",
    ]
    if features["spectral_fit"]:
        command += ["-fit-spectral-pol", str(int(args.spectral_terms))]
    if args.no_negative:
        # For a model to calibrate against: sky is positive, and the
        # negative lobes that calibration errors leave are not sky.
        command.append("-no-negative")
    if features["local_rms"] and supported(capabilities, "-local-rms"):
        command.append("-local-rms")
    parallel = parallel_deconvolution_size(args, capabilities)
    if parallel and supported(capabilities, "-parallel-deconvolution"):
        command += ["-parallel-deconvolution", str(parallel)]
    if features["multiscale"]:
        scales = features.get("scales") or args.scales
        command += [
            "-multiscale",
            "-multiscale-scales", ",".join(str(int(s)) for s in scales),
        ]
    return command + [str(ms)]


class CommandFailed(RuntimeError):
    """A command that did not exit cleanly; ``status`` is its exit status.

    A negative status is the signal that killed it: -6 is an abort, which
    is how a C++ program ends on a corrupted heap, and -11 a segmentation
    fault. Those are faults of the program, not of what it was asked.
    """

    def __init__(self, message, status):
        super().__init__(message)
        self.status = status

    @property
    def crashed(self):
        # A wrapper such as singularity may report the signal as 128 + n.
        return self.status < 0 or self.status in (134, 139)


class Diverged(CommandFailed):
    """A clean stopped because it was diverging; its image is no image."""

    def __init__(self, message):
        super().__init__(message, 0)


NOISE_ESTIMATE = re.compile(
    r"Estimated standard deviation of background noise:\s*"
    r"([0-9.]+(?:[eE][-+]?[0-9]+)?)\s*(nJy|mJy|µJy|μJy|uJy|Jy)"
)
JANSKY = {"nJy": 1e-9, "mJy": 1e-3, "µJy": 1e-6, "μJy": 1e-6,
          "uJy": 1e-6, "Jy": 1.0}
# How far the residual's noise may climb above its lowest before the clean
# is taken to be diverging. In a clean that converges it does not rise: in
# 254 other cleans of more than one major iteration in the Band 4 target run
# of 2026-10-29 (models of the self-calibration rounds, peel images, final
# images) it never went above its lowest. In G0213-0251's robust-0.5 image
# it went 58.6, 62.9, 109, 258, 769 uJy/beam and on to 19 mJy/beam.
DIVERGED_RISE = 2.0


class DivergenceWatch:
    """Reads WSClean's output as it comes; says when the clean diverges.

    WSClean gives the residual's noise at the start of each major
    iteration. A rise to twice its lowest is a clean adding error, not
    taking out sky, and it does not come back: WSClean's own reset of a
    sub-image that diverged did not stop it. Stopping there saves the hour
    such a run goes on for and lets a gentler clean have the time.
    """

    def __init__(self, rise=DIVERGED_RISE):
        self.rise = rise
        self.noise = []

    def __call__(self, line):
        found = NOISE_ESTIMATE.search(line)
        if found is None:
            return None
        value = float(found.group(1)) * JANSKY[found.group(2)]
        self.noise.append(value)
        lowest = min(self.noise)
        if value > self.rise * lowest:
            return (
                "the residual's noise rose from "
                f"{lowest * 1e6:.1f} to {value * 1e6:.1f} uJy/beam between "
                "major iterations: the clean was diverging"
            )
        return None


def stop(process, grace=60):
    """End a running command, by force if it does not end when asked."""

    process.terminate()
    try:
        process.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        process.kill()


def run_command(command, log_path, watch=None):
    """Run a command, showing its output and keeping it in a log file.

    ``watch`` is given each line of output; a reason it returns stops the
    command, and ``Diverged`` is raised with it.
    """

    reason = None
    with open(log_path, "w") as log:
        log.write(" ".join(shlex.quote(part) for part in command) + "\n\n")
        log.flush()
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
            log.flush()
            reason = watch(line) if watch is not None else None
            if reason:
                stop(process)
                log.write(f"\nStopped: {reason}.\n")
                sys.stdout.write(f"\nStopped: {reason}.\n")
                sys.stdout.flush()
                break
        process.stdout.close()
        code = process.wait()
    if reason:
        raise Diverged(f"{command[0]} was stopped: {reason}; see {log_path}")
    if code != 0:
        raise CommandFailed(
            f"{command[0]} exited with status {code}; see {log_path}", code
        )


def product(name, kind):
    """A WSClean output image: the joined (MFS) one, however it is named."""

    name = Path(name)
    for pattern in (
        f"{name.name}-MFS-{kind}.fits",
        f"{name.name}-MFS-?-{kind}.fits",
        f"{name.name}-{kind}.fits",
    ):
        found = sorted(name.parent.glob(pattern))
        if found:
            return found[0]
    return None


def fits_plane(path):
    """First plane of a FITS image as a 2-D array indexed [x, y].

    WSClean writes plain single-HDU images, which numpy reads directly:
    80-character header cards in 2880-byte blocks, then big-endian pixels
    with the first axis varying fastest. No astropy or CASA is needed.
    """

    import numpy as np

    header = {}
    with open(path, "rb") as handle:
        finished = False
        while not finished:
            block = handle.read(2880)
            if len(block) < 2880:
                raise ValueError(f"{path} is not a complete FITS file")
            for start in range(0, 2880, 80):
                card = block[start:start + 80].decode("ascii", "replace")
                key = card[:8].strip()
                if key == "END":
                    finished = True
                    break
                if card[8:10] == "= ":
                    header[key] = card[10:].split("/")[0].strip()
        bitpix = int(header["BITPIX"])
        shape = [
            int(header[f"NAXIS{axis}"])
            for axis in range(1, int(header["NAXIS"]) + 1)
        ]
        if len(shape) < 2:
            raise ValueError(f"{path} holds no image")
        dtype = {
            -32: ">f4", -64: ">f8", 16: ">i2", 32: ">i4", 8: "u1",
        }[bitpix]
        count = shape[0] * shape[1]
        data = np.frombuffer(
            handle.read(count * abs(bitpix) // 8), dtype=dtype, count=count
        )
    plane = data.reshape(shape[1], shape[0]).astype(float)
    plane = plane * float(header.get("BSCALE", 1.0)) + float(
        header.get("BZERO", 0.0)
    )
    # FITS rows are declination; the CASA tools index [x, y].
    return plane.T


def read_image(path):
    """A CASA image (a directory) or a FITS file as a 2-D array."""

    if Path(path).is_dir():
        return DIAGNOSTICS.read_plane(path)[0]
    return fits_plane(path)


def image_statistics(path, grid_cells):
    plane = read_image(path)
    grid = DIAGNOSTICS.rms_grid(plane, grid_cells)
    return {
        "image": str(path),
        "robust_rms_Jy_beam": DIAGNOSTICS.robust_rms(plane),
        "max_Jy_beam": float(plane.max()),
        "min_Jy_beam": float(plane.min()),
        "grid": DIAGNOSTICS.summarize_grid(grid, plane.shape),
        "grid_rms_Jy_beam": grid,
    }


def data_noise_of(report):
    """The noise image's measurement, under its name or the earlier one."""

    report = report or {}
    return report.get("data_noise") or report.get("stokes_v") or {}


def reading(report):
    """The outcome in sentences."""

    def micro(value):
        return f"{value * 1e6:.1f}"

    lines = []
    measured = data_noise_of(report)
    noise = measured.get("robust_rms_Jy_beam")
    result = report.get("stokes_i") or {}
    image = result.get("robust_rms_Jy_beam")
    quiet = (result.get("grid") or {}).get("quietest_tenth_rms")
    if noise:
        lines.append(
            f"WSClean, robust {report['robust']:g}: the data's noise "
            f"(dirty Stokes {measured.get('stokes', 'V')}) is "
            f"{micro(noise)} uJy/beam."
        )
    elif "noise_stokes" in report and report["noise_stokes"] is None:
        lines.append(
            "No noise image: the data hold no Stokes V (XX and YY only), and "
            "Stokes Q carries the bright sources' off-axis leakage."
        )
    if image:
        text = f"The Stokes I image is at {micro(image)} uJy/beam"
        if quiet:
            text += f", {micro(quiet)} in its quietest tenth"
        if noise:
            text += f": {image / noise:.1f} times the noise"
            if quiet:
                text += f", {quiet / noise:.1f} times in the quietest tenth"
        lines.append(text + ".")
    asymmetry = report.get("source_asymmetry")
    if asymmetry:
        lines.append(
            f"Around the {asymmetry['sources']} brightest sources the "
            "lopsided part of the image reaches "
            f"{100 * asymmetry['antisymmetric_peak']:.1f}% of the peak "
            f"(rms {100 * asymmetry['antisymmetric_rms']:.1f}%), and the "
            "deepest negative beside them is "
            f"{100 * asymmetry['deepest_negative_median']:.1f}% of the peak."
        )
    attempts = report.get("stokes_i_attempts") or []
    diverged = [a for a in attempts
                if str(a.get("outcome", "")).startswith("diverged")]
    crashed = [a for a in attempts if a.get("outcome") != "completed"
               and a not in diverged]
    if image and crashed:
        lines.append(
            f"WSClean aborted {len(crashed)} time(s) before this; the image "
            f"was made with: {attempts[-1]['deconvolution']}."
        )
    if image and diverged:
        lines.append(
            f"The clean diverged {len(diverged)} time(s) and was stopped "
            f"({', '.join(a['deconvolution'] for a in diverged)}); the image "
            f"was made with: {attempts[-1]['deconvolution']}."
        )
    reference = report.get("reference") or {}
    before = reference.get("robust_rms_Jy_beam")
    if before and image:
        before_quiet = (reference.get("grid") or {}).get("quietest_tenth_rms")
        text = (
            f"The tclean image it is compared with is at {micro(before)} "
            "uJy/beam"
        )
        if before_quiet:
            text += f" ({micro(before_quiet)} in its quietest tenth)"
        lines.append(
            text + ". The two have different beams unless the weighting is "
            "the same, so compare each with its own noise."
        )
    return lines


def brightest_sources(plane, count, separation=40):
    """Pixels of the ``count`` brightest peaks, at least ``separation`` apart."""

    import numpy as np

    work = np.array(plane, dtype=float, copy=True)
    work[~np.isfinite(work)] = -np.inf
    found = []
    for _ in range(int(count)):
        x, y = np.unravel_index(int(np.argmax(work)), work.shape)
        if not np.isfinite(work[x, y]):
            break
        found.append((int(x), int(y)))
        work[
            max(0, x - separation):x + separation + 1,
            max(0, y - separation):y + separation + 1,
        ] = -np.inf
    return found


def peak_near(plane, pixel, radius=2):
    x, y = pixel
    box = plane[
        max(0, x - radius):x + radius + 1,
        max(0, y - radius):y + radius + 1,
    ]
    return float(box.max())


def time_split_table(planes, count=8, minimum_ratio=20.0, locate=None,
                     typical_noise=False):
    """Flux of the brightest sources in each interval, against their mean.

    Sources are found in the mean of the interval images (or in ``locate``)
    and measured at the same pixels in each. Only those well above the
    noise of one interval are used: a dirty image of a fraction of the data
    is noisy, and a faint source would say more about that than about the
    gains. ``typical_noise`` takes the median of the images' noise where
    one of them (the shortest baselines) is far noisier than the rest.
    """

    import numpy as np

    mean = sum(planes) / float(len(planes))
    levels = [DIAGNOSTICS.robust_rms(plane) or 0.0 for plane in planes]
    noise = float(np.median(levels)) if typical_noise else max(levels)
    sources = []
    for pixel in brightest_sources(mean if locate is None else locate, count):
        fluxes = [peak_near(plane, pixel) for plane in planes]
        average = float(np.mean(fluxes))
        if noise > 0 and average < minimum_ratio * noise:
            continue
        sources.append({
            "pixel": list(pixel),
            "mean_Jy_beam": average,
            "flux_Jy_beam": fluxes,
            "ratio": [flux / average for flux in fluxes],
        })
    intervals = []
    for index in range(len(planes)):
        ratios = [source["ratio"][index] for source in sources]
        intervals.append({
            "interval": index,
            "sources": len(ratios),
            "median_ratio": float(np.median(ratios)) if ratios else None,
            "scatter": float(np.std(ratios)) if len(ratios) > 1 else None,
        })
    return {
        "interval_noise_Jy_beam": noise,
        "sources": sources,
        "intervals": intervals,
    }


def time_split_reading(table):
    """What the interval fluxes say, in sentences."""

    medians = [
        entry["median_ratio"] for entry in table["intervals"]
        if entry["median_ratio"] is not None
    ]
    if not medians:
        return [
            "No source is bright enough against the noise of one interval "
            "to compare the intervals."
        ]
    lines = [
        "Flux of the " + str(table["intervals"][0]["sources"])
        + " brightest sources in each interval, relative to their mean: "
        + ", ".join(f"{value:.3f}" for value in medians) + "."
    ]
    scatters = [
        entry["scatter"] for entry in table["intervals"]
        if entry["scatter"] is not None
    ]
    spread = max(medians) - min(medians)
    typical = max(scatters) if scatters else 0.0
    if spread > max(0.05, 3 * typical / max(1, len(table["sources"])) ** 0.5):
        lines.append(
            f"The sources change together by {100 * spread:.0f}% between "
            "intervals: the amplitude scale moves within the observation, "
            "which one gain per scan cannot follow."
        )
    else:
        lines.append(
            f"The intervals agree to {100 * spread:.1f}% (source-to-source "
            f"scatter up to {100 * typical:.1f}%): the amplitude scale is "
            "steady through the observation."
        )
    return lines


def time_split_command(args, ms, name, imsize, cell, channel_range,
                       temp_dir, capabilities):
    """Dirty Stokes I images of consecutive stretches of the data."""

    command = shlex.split(args.wsclean) + [
        "-name", str(name),
        "-data-column", args.data_column,
        "-no-update-model-required",
        "-temp-dir", str(temp_dir),
        "-j", str(int(args.threads)),
        "-abs-mem", f"{float(args.memory_gb):g}",
        *gridder_arguments(capabilities),
        "-size", str(int(imsize)), str(int(imsize)),
        "-scale", f"{cell:.4f}asec",
        "-weight", "briggs", f"{float(args.robust):g}",
        "-pol", "I",
        "-intervals-out", str(int(args.time_split)),
        "-niter", "0",
    ]
    if channel_range:
        command += [
            "-channel-range", str(int(channel_range[0])),
            str(int(channel_range[1])),
        ]
    return command + [str(ms)]


def run_time_split(args, ms, workdir, label, imsize, cell, channel_range,
                   temp_dir, capabilities):
    """Image the intervals and compare the sources between them."""

    name = workdir / f"timesplit_{label}"
    outfile = workdir / f"timesplit_{label}.json"
    for path in workdir.glob(name.name + "-*"):
        path.unlink()
    command = time_split_command(
        args, ms, name, imsize, cell, channel_range, temp_dir, capabilities
    )
    report = {
        "target": args.target,
        "measurement_set": str(ms),
        "robust": float(args.robust),
        "interval_count": int(args.time_split),
        "command": " ".join(shlex.quote(part) for part in command),
    }
    print(
        f"\nWSClean: {args.time_split} dirty Stokes I images, one per "
        "stretch of the observation ...",
        flush=True,
    )
    run_command(command, workdir / f"timesplit_{label}.log")
    shutil.rmtree(temp_dir, ignore_errors=True)

    images = []
    for index in range(int(args.time_split)):
        interval = Path(f"{name}-t{index:04d}")
        image = product(interval, "image") or product(interval, "dirty")
        if image is None:
            raise RuntimeError(
                f"WSClean wrote no image for interval {index}."
            )
        images.append(image)
    report["images"] = [str(image) for image in images]
    report.update(time_split_table([fits_plane(image) for image in images]))
    report["reading"] = time_split_reading(report)
    with open(outfile, "w") as handle:
        json.dump(DIAGNOSTICS._plain(report), handle, indent=2,
                  allow_nan=False)
    print("\n" + "\n".join(report["reading"]))
    print(f"\nReport written: {outfile}")
    return report


def source_asymmetry(plane, count=8, radius=36, margin=60,
                     minimum_ratio=25.0):
    """How lopsided the brightest sources are, in fractions of their peaks.

    Cut-outs around the brightest sources are divided by their peaks and
    averaged. A point-spread function is unchanged by a half-turn about its
    peak, so whatever in the average changes sign under one is not the
    beam: it is phase error common to the field (or, for one source alone,
    its own shape, which is why several are averaged). The deepest negative
    beside each source is given as well, since that is what the eye finds.
    """

    import numpy as np

    noise = DIAGNOSTICS.robust_rms(plane) or 0.0
    size_x, size_y = plane.shape
    width = 2 * radius + 1
    stack = np.zeros((width, width))
    used, negatives = [], []
    for x, y in brightest_sources(plane, 4 * count, separation=70):
        peak = float(plane[x, y])
        if noise > 0 and peak < minimum_ratio * noise:
            break
        if not (margin <= x < size_x - margin and margin <= y < size_y - margin):
            continue
        stack += plane[x - radius:x + radius + 1,
                       y - radius:y + radius + 1] / peak
        negatives.append(float(
            plane[x - margin:x + margin + 1, y - margin:y + margin + 1].min()
            / peak
        ))
        used.append([int(x), int(y)])
        if len(used) >= count:
            break
    if not used:
        return None
    stack /= len(used)
    odd = 0.5 * (stack - stack[::-1, ::-1])
    xx, yy = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    distance = np.hypot(xx, yy)
    ring = (distance >= 6) & (distance <= radius)
    index = np.unravel_index(
        int(np.argmax(np.where(ring, np.abs(odd), 0.0))), odd.shape
    )
    return {
        "sources": len(used),
        "pixels": used,
        "antisymmetric_peak": float(abs(odd[index])),
        "antisymmetric_peak_offset": [
            int(index[0] - radius), int(index[1] - radius)
        ],
        "antisymmetric_rms": float(np.sqrt(np.mean(odd[ring] ** 2))),
        "deepest_negative_median": float(np.median(negatives)),
    }


def core_model_command(args, ms, name, imsize, cell, channel_range,
                       temp_dir, capabilities):
    """Clean the cores of the bright sources and write the model column."""

    command = shlex.split(args.wsclean) + [
        "-name", str(name),
        "-data-column", args.data_column,
        # No -no-update-model-required: the model column is the product.
        "-temp-dir", str(temp_dir),
        "-j", str(int(args.threads)),
        "-abs-mem", f"{float(args.memory_gb):g}",
        *gridder_arguments(capabilities),
        "-size", str(int(imsize)), str(int(imsize)),
        "-scale", f"{cell:.4f}asec",
        "-weight", "briggs", f"{float(args.robust):g}",
        "-pol", "I",
        "-join-channels",
        "-channels-out", str(int(args.channels_out)),
    ]
    if channel_range:
        command += [
            "-channel-range", str(int(channel_range[0])),
            str(int(channel_range[1])),
        ]
    if args.core_fraction:
        # To a fraction of the field's own brightest pixel, which need not
        # be known beforehand: WSClean ends its first major iteration when
        # the peak has fallen to (1 - mgain) of what it began at, and with
        # one major iteration allowed it ends there altogether.
        command += [
            "-niter", "50000",
            "-nmiter", "1",
            "-mgain", f"{1.0 - float(args.core_fraction):g}",
            # Never into the noise: in a faint field a fraction of the
            # brightest pixel can sit at 4 sigma (G0231-0533 in Band 5),
            # where WSClean cleans noise for hours and the cores model fills
            # with noise peaks. WSClean stops at whichever comes first.
            "-auto-threshold", f"{CORE_MODEL_MIN_SIGMA:g}",
            "-no-negative",
        ]
    else:
        command += [
            "-niter", "50000",
            "-nmiter", "8",
            "-mgain", f"{float(args.mgain):g}",
            # Jy; nothing fainter than this goes into the model, and nothing
            # negative: the sources are to be compact and positive.
            "-threshold", f"{float(args.core_model) * 1e-3:g}",
            "-no-negative",
        ]
    if int(args.spectral_terms) > 0:
        command += ["-fit-spectral-pol", str(int(args.spectral_terms))]
    return command + [str(ms)]


def inside(path, parent):
    """Whether ``path`` lies in ``parent``, both as they resolve on disk."""

    try:
        Path(path).resolve().relative_to(Path(parent).resolve())
    except ValueError:
        return False
    return True


LOCK_NAME = "running.lock"


def process_alive(pid):
    try:
        os.kill(int(pid), 0)
    except (ProcessLookupError, ValueError, TypeError):
        return False
    except PermissionError:
        return True
    return True


def slurm_job_gone(job):
    """Whether a SLURM job has ended, or None when that cannot be asked.

    A job that has ended is no longer listed by ``squeue``; one that is
    queued, running or being requeued is.
    """

    if not job or shutil.which("squeue") is None:
        return None
    try:
        done = subprocess.run(
            ["squeue", "-h", "-j", str(job), "-o", "%T"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if done.returncode != 0:
        # squeue refuses a job id it no longer knows.
        return "Invalid job id" in (done.stderr or "")
    return not done.stdout.strip()


def claim_of_this_run():
    """What a lock says of the run that holds it."""

    mine = {
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "since": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
    }
    if os.environ.get("SLURM_JOB_ID"):
        mine["slurm_job"] = os.environ["SLURM_JOB_ID"]
    return mine


def held_by(path):
    """The run that holds a lock, or None if there is none or it is gone.

    A claim is gone when its process is gone from this machine, or when it
    was made by a SLURM job that has ended: a job on another node cannot
    be seen any other way, and a killed job leaves its claim behind.
    """

    path = Path(path)
    try:
        other = json.loads(path.read_text())
    except FileNotFoundError:
        return None
    except (OSError, ValueError):
        other = {}
    if other.get("host") == socket.gethostname() and other.get("pid"):
        if not process_alive(other.get("pid")):
            return None
    elif slurm_job_gone(other.get("slurm_job")):
        return None
    return other


@contextlib.contextmanager
def claimed(workdir):
    """Hold a work folder for one run; refuse it while another run has it.

    Two runs given the same folder write the same image names and the same
    temporary files, and each removes the other's as it tidies up: on
    2026-10-08 a command started twice lost both runs that way, minutes in.
    A claim left by a run that was killed is taken over when its process is
    gone from this machine, or when the SLURM job that made it has ended;
    from another machine otherwise that cannot be seen, so the message says
    which file to remove.
    """

    path = Path(workdir) / LOCK_NAME
    mine = claim_of_this_run()
    for attempt in range(2):
        try:
            handle = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            other = held_by(path)
            if other is None and attempt == 0:
                path.unlink(missing_ok=True)
                continue
            other = other or {}
            raise SystemExit(
                f"{workdir} is in use by another run (process "
                f"{other.get('pid', '?')} on {other.get('host', '?')}"
                + (f", SLURM job {other['slurm_job']}"
                   if other.get("slurm_job") else "")
                + f", since {other.get('since', '?')}). Two runs in one "
                "folder remove each other's files. Wait for it, or give "
                "--workdir another folder; if that run is gone, remove "
                f"{path}."
            ) from None
        with os.fdopen(handle, "w") as lock:
            json.dump(mine, lock)
        break
    try:
        yield path
    finally:
        path.unlink(missing_ok=True)


def run_core_model(args, ms, workdir, label, imsize, cell, channel_range,
                   temp_dir, capabilities):
    """Make the cores-only model and leave it in MODEL_DATA."""

    if inside(ms, args.outdir):
        raise SystemExit(
            f"{ms} is inside the pipeline's products ({args.outdir}). "
            "--core-model writes a model column: give it a copy of the "
            "data with --ms."
        )
    name = workdir / f"coremodel_{label}"
    outfile = workdir / f"coremodel_{label}.json"
    for path in workdir.glob(name.name + "-*"):
        path.unlink()
    command = core_model_command(
        args, ms, name, imsize, cell, channel_range, temp_dir, capabilities
    )
    report = {
        "target": args.target,
        "measurement_set": str(ms),
        "threshold_Jy_beam": (
            None if args.core_fraction else float(args.core_model) * 1e-3
        ),
        "fraction_of_peak": (
            float(args.core_fraction) if args.core_fraction else None
        ),
        "command": " ".join(shlex.quote(part) for part in command),
    }
    print(
        "\nWSClean: model of everything above "
        + (f"{100 * float(args.core_fraction):g}% of the brightest pixel"
           if args.core_fraction else f"{args.core_model:g} mJy/beam")
        + f", positive components only, written to MODEL_DATA of {ms} ...",
        flush=True,
    )
    run_command(command, workdir / f"coremodel_{label}.log")
    shutil.rmtree(temp_dir, ignore_errors=True)

    model = product(name, "model")
    if model is None:
        raise RuntimeError("WSClean wrote no model image.")
    lopsided = None
    try:
        import numpy as np

        # The dirty image is the data before any cleaning: its beam is
        # unchanged by a half-turn, so the lopsided part measured here is
        # the state of the calibration and nothing else. It is the quick
        # reading between rounds, in place of a full image.
        dirty = product(name, "dirty")
        if dirty is not None:
            lopsided = source_asymmetry(fits_plane(dirty))
            report["dirty_asymmetry"] = lopsided
        plane = fits_plane(model)
        report["model"] = {
            "image": str(model),
            "total_flux_Jy": float(plane.sum()),
            "pixels_with_flux": int(np.count_nonzero(plane)),
            "most_negative_Jy": float(plane.min()),
        }
        lines = [
            f"The model holds {report['model']['total_flux_Jy'] * 1e3:.0f} "
            f"mJy in {report['model']['pixels_with_flux']} pixels, none "
            "negative." if report["model"]["most_negative_Jy"] >= 0 else
            "WARNING: the model has negative pixels; -no-negative was not "
            "honoured."
        ]
    except ImportError:
        lines = ["numpy is not available, so the model was not measured."]
    if lopsided:
        lines.insert(0, (
            f"Before any cleaning, around the {lopsided['sources']} "
            "brightest sources the lopsided part of the image reaches "
            f"{100 * lopsided['antisymmetric_peak']:.1f}% of the peak "
            f"(rms {100 * lopsided['antisymmetric_rms']:.2f}%)."
        ))
    lines.append(
        "MODEL_DATA is written. Next, from the Python that has CASA: "
        f"python3 experiment_phase_offsets.py --ms {ms}"
    )
    report["reading"] = lines
    with open(outfile, "w") as handle:
        json.dump(DIAGNOSTICS._plain(report), handle, indent=2,
                  allow_nan=False)
    print("\n" + "\n".join(lines))
    print(f"\nReport written: {outfile}")
    return report


def uv_split_command(args, ms, name, imsize, cell, channel_range,
                     temp_dir, capabilities, low, high):
    """A dirty Stokes I image from one range of baseline length."""

    command = shlex.split(args.wsclean) + [
        "-name", str(name),
        "-data-column", args.data_column,
        "-no-update-model-required",
        "-temp-dir", str(temp_dir),
        "-j", str(int(args.threads)),
        "-abs-mem", f"{float(args.memory_gb):g}",
        *gridder_arguments(capabilities),
        "-size", str(int(imsize)), str(int(imsize)),
        "-scale", f"{cell:.4f}asec",
        # Every baseline in the range counts alike: the question is what
        # flux the range sees, not what beam it makes.
        "-weight", "natural",
        "-pol", "I",
        "-minuv-l", f"{1000.0 * float(low):g}",
        "-maxuv-l", f"{1000.0 * float(high):g}",
        "-niter", "0",
    ]
    if channel_range:
        command += [
            "-channel-range", str(int(channel_range[0])),
            str(int(channel_range[1])),
        ]
    return command + [str(ms)]


def uv_split_reading(table, edges):
    """What the fluxes by baseline length say, in sentences."""

    ranges = [
        f"{low:g}-{high:g}" for low, high in zip(edges[:-1], edges[1:])
    ]
    medians = [entry["median_ratio"] for entry in table["intervals"]]
    if any(value is None for value in medians):
        return [
            "No source is bright enough against the noise of these images "
            "to compare the baseline ranges."
        ]
    lines = [
        "Flux of the " + str(table["intervals"][0]["sources"])
        + " brightest sources by baseline length, relative to their mean: "
        + ", ".join(
            f"{name} klambda {value:.3f}"
            for name, value in zip(ranges, medians)
        ) + "."
    ]
    # The shortest range blends neighbouring sources; judge on the rest.
    first = 1 if len(medians) > 2 else 0
    change = medians[-1] - medians[first]
    last = [source["ratio"][-1] / source["ratio"][first]
            for source in table["sources"] if source["ratio"][first] > 0]
    if change > 0.05:
        lines.append(
            f"Beyond the shortest range the sources are "
            f"{100 * change:.0f}% brighter on the longest baselines than on "
            "the shorter ones. No source is: the amplitude scale depends on "
            "baseline length, which is what antennas on different scales "
            "give."
        )
    elif change < -0.05:
        lines.append(
            "The sources are fainter on the longest baselines than on the "
            "shorter ones, by factors from "
            f"{min(last):.2f} to {max(last):.2f}: they are resolved. That "
            "is their structure, not the amplitude scale, which this test "
            "cannot then judge; a field of resolved sources needs a "
            "multi-scale clean."
        )
    else:
        lines.append(
            f"Beyond the shortest range the baseline ranges agree to "
            f"{100 * abs(change):.1f}%: the sources are compact and the "
            "amplitude scale does not depend on baseline length."
        )
    return lines


def run_uv_split(args, ms, workdir, label, imsize, cell, channel_range,
                 temp_dir, capabilities):
    """Image each range of baseline length and compare the sources."""

    edges = [float(edge) for edge in args.uv_split]
    if len(edges) < 3 or any(b <= a for a, b in zip(edges[:-1], edges[1:])):
        raise SystemExit(
            "--uv-split needs at least three increasing edges in "
            "kilo-wavelengths, for example --uv-split 0.6 3 8 20 35 70."
        )
    outfile = workdir / f"uvsplit_{label}.json"
    report = {
        "target": args.target,
        "measurement_set": str(ms),
        "edges_klambda": edges,
        "commands": [],
    }
    planes = []
    for index, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
        name = workdir / f"uvsplit_{label}_b{index}"
        for path in workdir.glob(name.name + "-*"):
            path.unlink()
        command = uv_split_command(
            args, ms, name, imsize, cell, channel_range, temp_dir,
            capabilities, low, high,
        )
        report["commands"].append(
            " ".join(shlex.quote(part) for part in command)
        )
        print(
            f"\nWSClean: dirty Stokes I from baselines of {low:g} to "
            f"{high:g} klambda ...",
            flush=True,
        )
        run_command(command, workdir / f"uvsplit_{label}_b{index}.log")
        image = product(name, "image") or product(name, "dirty")
        if image is None:
            raise RuntimeError(
                f"WSClean wrote no image for {low:g}-{high:g} klambda."
            )
        planes.append(fits_plane(image))
    shutil.rmtree(temp_dir, ignore_errors=True)

    # Sources are located where they are sharpest: the longest baselines.
    report.update(time_split_table(
        planes, minimum_ratio=8.0, locate=planes[-1], typical_noise=True,
    ))
    report["reading"] = uv_split_reading(report, edges)
    with open(outfile, "w") as handle:
        json.dump(DIAGNOSTICS._plain(report), handle, indent=2,
                  allow_nan=False)
    print("\n" + "\n".join(report["reading"]))
    print(f"\nReport written: {outfile}")
    return report


def reference_grid(reference_image):
    """Size and cell of the tclean image, or None where CASA is absent."""

    try:
        plane, cell = DIAGNOSTICS.read_plane(reference_image)
    except ImportError:
        return None
    return plane.shape[0], cell


def measure(report, names, grid_cells):
    """Measure the WSClean images that exist; say what could not be done."""

    notes = []
    try:
        import numpy  # noqa: F401
    except ImportError:
        notes.append(
            "numpy is not available to this Python, so the images were not "
            "measured. Run the same command with --measure-only from a "
            "Python that has numpy (CASA's does)."
        )
        report["notes"] = notes
        return report

    noise_name = names.get("noise") or names.get("V")
    noise_image = (
        product(noise_name, "image") or product(noise_name, "dirty")
        if noise_name is not None else None
    )
    if noise_image is not None:
        plane = fits_plane(noise_image)
        entry = {
            "image": str(noise_image),
            "stokes": report.get("noise_stokes") or "V",
            "robust_rms_Jy_beam": DIAGNOSTICS.robust_rms(plane),
            "max_abs_Jy_beam": float(abs(plane).max()),
        }
        report["data_noise"] = entry
        if entry["stokes"] == "V":
            # The name earlier reports, and their readers, use.
            report["stokes_v"] = entry
        del plane

    image = product(names["I"], "image")
    if image is not None:
        report["stokes_i"] = image_statistics(image, grid_cells)
        report["source_asymmetry"] = source_asymmetry(fits_plane(image))
        residual = product(names["I"], "residual")
        if residual is not None:
            report["stokes_i_residual"] = image_statistics(
                residual, grid_cells
            )
            report["stokes_i_residual"].pop("grid_rms_Jy_beam", None)

    reference_image = report.get("reference_image")
    if reference_image and Path(reference_image).is_dir():
        try:
            report["reference"] = image_statistics(reference_image, grid_cells)
        except ImportError:
            notes.append(
                "CASA is not importable here, so the tclean image was not "
                "measured for comparison; --measure-only from CASA's Python "
                "adds it."
            )
    report["notes"] = notes
    report["reading"] = reading(report)
    return report


def run(args):
    # Whatever branch claims the work folder, the claim ends with the run.
    with contextlib.ExitStack() as claims:
        return make(args, claims)


def make(args, claims):
    ms, reference_image = DIAGNOSTICS.resolve_products(
        args.outdir, args.target
    )
    if args.ms:
        ms = Path(args.ms)
    if args.reference:
        reference_image = Path(args.reference)
    if ms is None or not Path(ms).is_dir():
        raise SystemExit(
            f"No self-calibrated Measurement Set found for {args.target} "
            f"under {args.outdir}. Pass --ms."
        )

    label = image_label(args.robust, getattr(args, "taper", 0))
    workdir = Path(args.workdir) / args.target
    temp_dir = workdir / "tmp"
    outfile = workdir / f"wsclean_{label}.json"
    noise = noise_stokes(args)
    if getattr(args, "measure_only", False):
        # Measured as it was made: an earlier report names its noise image.
        try:
            with open(outfile) as handle:
                noise = json.load(handle).get("noise_stokes", "V")
        except (OSError, ValueError, AttributeError):
            pass
    names = {"I": workdir / f"wsclean_{label}_I"}
    if noise:
        names["noise"] = workdir / f"wsclean_{label}_{noise}"

    def save(report):
        with open(outfile, "w") as handle:
            json.dump(DIAGNOSTICS._plain(report), handle, indent=2,
                      allow_nan=False)

    def finish(report):
        measure(report, names, args.grid)
        save(report)
        for line in (report.get("reading") or []) + report["notes"]:
            print(line)
        print(f"\nReport written: {outfile}")
        return report

    if args.measure_only:
        try:
            with open(outfile) as handle:
                report = json.load(handle)
        except (OSError, ValueError) as exc:
            raise SystemExit(
                f"No report of a WSClean run at {outfile} ({exc}). Run "
                "without --measure-only first."
            ) from exc
        return finish(report)

    prefix = shlex.split(args.wsclean)
    if not prefix:
        raise SystemExit("--wsclean is empty.")
    if not args.dry_run and shutil.which(prefix[0]) is None:
        raise SystemExit(
            f"'{prefix[0]}' is not on PATH. Pass --wsclean with the command "
            "that runs WSClean here, in quotes if it has several words, for "
            "example --wsclean \"singularity exec /path/to/container.sif "
            "wsclean\"."
        )

    imsize, cell = args.imsize, args.cell
    have_reference = (
        reference_image is not None and Path(reference_image).is_dir()
    )
    if not (imsize and cell):
        grid = reference_grid(reference_image) if have_reference else None
        if grid is None:
            raise SystemExit(
                "The image grid could not be read from the tclean image"
                + (" because CASA is not importable in this Python"
                   if have_reference else "")
                + ". Pass --imsize and --cell (arcsec), for example "
                "--imsize 5120 --cell 1.0."
            )
        imsize, cell = imsize or grid[0], cell or grid[1]

    if args.channel_range:
        channel_range = tuple(args.channel_range)
    elif args.all_channels:
        channel_range = None
    else:
        print("Finding the channels that hold data ...", flush=True)
        try:
            channel_range = measured_channel_range(ms)
        except ImportError:
            raise SystemExit(
                "The channels that hold data are measured with CASA, which "
                "is not importable in this Python. Pass --channel-range "
                "FIRST END (END exclusive), or --all-channels."
            ) from None
        if channel_range is None:
            print(
                "WARNING: the flag summary gave no channel range; imaging "
                "every channel.",
                flush=True,
            )
        else:
            print(
                f"  channels {channel_range[0]} to {channel_range[1] - 1}",
                flush=True,
            )

    capabilities = None
    if not args.dry_run:
        capabilities = wsclean_capabilities(prefix)
        version = capabilities["version"]
        print(
            "WSClean "
            + (f"{version[0]}.{version[1]}" if version else "(version unknown)")
            + f": {' '.join(gridder_arguments(capabilities))}",
            flush=True,
        )
        if version is None and not capabilities["options"]:
            raise SystemExit(
                f"'{args.wsclean} --version' did not answer as WSClean. "
                "Check the command given as --wsclean."
            )

    if args.write_model and not args.dry_run and inside(ms, args.outdir):
        raise SystemExit(
            f"{ms} is inside the pipeline's products ({args.outdir}). "
            "--write-model writes a model column: give it a copy of the "
            "data with --ms."
        )

    if args.core_fraction and not 0.0 < float(args.core_fraction) < 1.0:
        raise SystemExit("--core-model-fraction is between 0 and 1.")
    if args.core_model or args.core_fraction:
        if args.dry_run:
            print(" ".join(shlex.quote(part) for part in core_model_command(
                args, ms, workdir / f"coremodel_{label}", imsize, cell,
                channel_range, temp_dir, capabilities,
            )))
            return None
        if inside(ms, args.outdir):
            raise SystemExit(
                f"{ms} is inside the pipeline's products ({args.outdir}). "
                "--core-model writes a model column: give it a copy of the "
                "data with --ms."
            )
        workdir.mkdir(parents=True, exist_ok=True)
        claims.enter_context(claimed(workdir))
        temp_dir.mkdir(parents=True, exist_ok=True)
        return run_core_model(
            args, ms, workdir, label, imsize, cell, channel_range, temp_dir,
            capabilities,
        )

    if args.uv_split:
        if args.dry_run:
            raise SystemExit("--dry-run does not cover --uv-split.")
        workdir.mkdir(parents=True, exist_ok=True)
        claims.enter_context(claimed(workdir))
        temp_dir.mkdir(parents=True, exist_ok=True)
        return run_uv_split(
            args, ms, workdir, label, imsize, cell, channel_range, temp_dir,
            capabilities,
        )

    if args.time_split:
        if args.dry_run:
            print(" ".join(shlex.quote(part) for part in time_split_command(
                args, ms, workdir / f"timesplit_{label}", imsize, cell,
                channel_range, temp_dir, capabilities,
            )))
            return None
        workdir.mkdir(parents=True, exist_ok=True)
        claims.enter_context(claimed(workdir))
        temp_dir.mkdir(parents=True, exist_ok=True)
        return run_time_split(
            args, ms, workdir, label, imsize, cell, channel_range, temp_dir,
            capabilities,
        )

    order = ([noise] if noise else []) + ["I"]
    commands = {
        stokes: wsclean_command(
            args, ms, names["I" if stokes == "I" else "noise"], imsize, cell,
            channel_range, stokes, temp_dir, capabilities,
        )
        for stokes in order
    }
    if args.dry_run:
        for stokes in order:
            print(" ".join(shlex.quote(part) for part in commands[stokes]))
        return None

    workdir.mkdir(parents=True, exist_ok=True)
    claims.enter_context(claimed(workdir))
    existing = [
        path for name in names.values()
        for path in workdir.glob(name.name + "-*")
    ]
    if existing and not args.overwrite:
        raise SystemExit(
            f"{workdir} already holds WSClean images for robust "
            f"{args.robust:g}. Pass --overwrite to remake them, "
            "--measure-only to measure them, or --workdir for a separate "
            "run."
        )
    for path in existing:
        path.unlink()
    temp_dir.mkdir(parents=True, exist_ok=True)

    report = {
        "target": args.target,
        "measurement_set": str(ms),
        "robust": float(args.robust),
        "imsize": int(imsize),
        "cell_arcsec": float(cell),
        "channel_range": list(channel_range) if channel_range else None,
        "wsclean_version": (capabilities or {}).get("version"),
        "noise_stokes": noise,
        "commands": {
            stokes: " ".join(shlex.quote(part) for part in command)
            for stokes, command in commands.items()
        },
        "reference_image": str(reference_image) if have_reference else None,
    }
    save(report)

    if noise:
        print(f"\nWSClean: dirty Stokes {noise} (the noise of the data) ...",
              flush=True)
        run_command(commands[noise], workdir / f"wsclean_{label}_{noise}.log")
        if (product(names["noise"], "image")
                or product(names["noise"], "dirty")) is None:
            raise RuntimeError(f"WSClean wrote no Stokes {noise} image.")
        measure(report, names, args.grid)
        save(report)
        measured = data_noise_of(report).get("robust_rms_Jy_beam")
        if measured:
            print(f"  robust rms {measured * 1e6:.1f} uJy/beam", flush=True)
    else:
        print(
            "\nNo noise image: the data hold XX and YY only, with no Stokes "
            "V, and Q carries the bright sources' off-axis leakage.",
            flush=True,
        )

    report["stokes_i_attempts"] = []
    ladder = feature_ladder(args)
    index = -1
    while index + 1 < len(ladder):
        index += 1
        features = ladder[index]
        command = wsclean_command(
            args, ms, names["I"], imsize, cell, channel_range, "I",
            temp_dir, capabilities, features,
        )
        attempt = {
            "deconvolution": describe(features),
            "command": " ".join(shlex.quote(part) for part in command),
        }
        report["stokes_i_attempts"].append(attempt)
        report["commands"]["I"] = attempt["command"]
        print(f"\nWSClean: Stokes I ({attempt['deconvolution']}) ...",
              flush=True)
        log_path = workdir / f"wsclean_{label}_I.log"
        try:
            run_command(command, log_path, watch=DivergenceWatch())
        except Diverged as failure:
            attempt["outcome"] = f"diverged: {failure}"
            save(report)
            gentler = calmer(features, args.scales)
            if gentler is None:
                raise
            # Keep the log of what diverged; clear what it left.
            log_path.replace(
                workdir / f"wsclean_{label}_I.diverged{index + 1}.log")
            for path in workdir.glob(names["I"].name + "-*"):
                path.unlink()
            ladder = ladder[:index + 1] + [gentler] + [
                rung for rung in ladder[index + 1:] if rung != gentler
            ]
            print(
                f"\nWSClean's clean diverged. Trying again with: "
                f"{describe(gentler)}.",
                flush=True,
            )
            continue
        except CommandFailed as failure:
            attempt["outcome"] = f"failed: {failure}"
            save(report)
            if not failure.crashed or index + 1 == len(ladder):
                raise
            # Keep the log of the crash; clear what it left of the image.
            log_path.replace(workdir / f"wsclean_{label}_I.crash{index + 1}.log")
            for path in workdir.glob(names["I"].name + "-*"):
                path.unlink()
            print(
                f"\nWSClean aborted (status {failure.status}). Trying again "
                f"with: {describe(ladder[index + 1])}.",
                flush=True,
            )
            continue
        attempt["outcome"] = "completed"
        break

    if product(names["I"], "image") is None:
        raise RuntimeError("WSClean wrote no Stokes I image.")
    shutil.rmtree(temp_dir, ignore_errors=True)
    print()
    return finish(report)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--outdir", default="uGMRT_Products")
    parser.add_argument("--target", required=True)
    parser.add_argument("--ms", help="Measurement Set, if not the target's "
                        "final self-calibrated one")
    parser.add_argument("--reference", help="CASA image to compare with, if "
                        "not the target's final self-calibration image")
    parser.add_argument("--workdir", default="wsclean_test")
    parser.add_argument("--wsclean", default="wsclean",
                        help="the command that runs WSClean; quote it when "
                        "it has several words, as with a container")
    parser.add_argument("--data-column", default="DATA",
                        help="DATA holds the self-calibrated visibilities "
                        "in the pipeline's self-calibration products")
    parser.add_argument("--robust", type=float, default=0.0)
    parser.add_argument("--taper", type=float, default=0.0, metavar="ARCSEC",
                        help="Gaussian taper of the visibility weights, as "
                        "the beam it gives in arcsec (default 0: none); "
                        "named in the images, as robust_0_taper15")
    parser.add_argument("--channels-out", type=int, default=6,
                        help="sub-bands deconvolved jointly (default 6)")
    parser.add_argument("--spectral-terms", type=int, default=3,
                        help="terms of the polynomial fitted through the "
                        "sub-bands; 0 for none (default 3)")
    parser.add_argument("--channel-range", type=int, nargs=2,
                        metavar=("FIRST", "END"),
                        help="channels to image, END exclusive; measured "
                        "from the flags if not given")
    parser.add_argument("--all-channels", action="store_true",
                        help="image every channel without measuring")
    parser.add_argument("--niter", type=int, default=500000)
    parser.add_argument("--nmiter", type=int, default=15,
                        help="most major iterations (default 15)")
    parser.add_argument("--mgain", type=float, default=0.8)
    parser.add_argument("--auto-mask", type=float, default=4.0)
    parser.add_argument("--auto-threshold", type=float, default=1.5)
    parser.add_argument("--scales", type=int, nargs="*", default=[0, 5, 15],
                        help="multi-scale sizes in pixels; none for a "
                        "point-source clean (default 0 5 15)")
    parser.add_argument("--no-local-rms", action="store_true",
                        help="mask and stop on a single noise value for "
                        "the whole image")
    parser.add_argument("--parallel-deconvolution", type=int,
                        help="sub-image size for parallel deconvolution, 0 "
                        "for none (default 2048, and none before WSClean "
                        "3.2, where it crashes)")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory-gb", type=float, default=20.0,
                        help="WSClean sizes itself from the node's memory "
                        "unless told (default 20)")
    parser.add_argument("--imsize", type=int)
    parser.add_argument("--cell", type=float, help="arcsec")
    parser.add_argument("--grid", type=int, default=10)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the WSClean commands and stop")
    parser.add_argument("--time-split", type=int, default=0, metavar="N",
                        help="instead of the image, make N dirty images of "
                        "consecutive stretches of the data and compare the "
                        "brightest sources between them")
    parser.add_argument("--no-negative", action="store_true",
                        help="clean with positive components only")
    parser.add_argument("--write-model", action="store_true",
                        help="let WSClean write the Stokes I model to "
                        "MODEL_DATA of the Measurement Set (a copy, not a "
                        "pipeline product), for a calibration step to use")
    parser.add_argument("--core-model", type=float, default=0.0,
                        metavar="MJY",
                        help="instead of the image, clean only what is "
                        "brighter than MJY mJy/beam with positive components "
                        "and write that model to MODEL_DATA of the "
                        "Measurement Set (a copy, not a pipeline product)")
    parser.add_argument("--core-model-fraction", dest="core_fraction",
                        type=float, default=0.0, metavar="FRACTION",
                        help="as --core-model, with the threshold a "
                        "fraction of the field's own brightest pixel "
                        "instead of a number of mJy/beam")
    parser.add_argument("--uv-split", type=float, nargs="+", metavar="KL",
                        help="instead of the image, make a dirty image "
                        "from each range of baseline length between these "
                        "edges (kilo-wavelengths) and compare the brightest "
                        "sources between them")
    parser.add_argument("--measure-only", action="store_true",
                        help="measure the images of an earlier run and "
                        "rewrite its report; runs no WSClean")
    parser.add_argument("--noise-stokes", default="auto",
                        choices=("auto", "V", "Q", "none"),
                        help="the dirty image that gives the data's noise "
                        "(default auto: Stokes V where the correlations the "
                        "pipeline found can form it, none for XX and YY "
                        "alone)")
    run(parser.parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
