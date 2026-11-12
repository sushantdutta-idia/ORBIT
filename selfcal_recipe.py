#!/usr/bin/env python3
"""Run the self-calibration recipe on one Measurement Set, start to finish.

The recipe is the sequence that was worked out by hand on uGMRT Band 4,
G0212-0335 (SCIENCE_AUDIT_PATCH_NOTES.md, findings 118 and 119). Each round
is a model made by WSClean and a solve made by CASA against it:

1. a model of the bright cores only, positive, and one phase per antenna
   for each scan; again until the phases found are small, since the first
   model is made from data that are still distorted;
2. the same model, and a delay per antenna with a phase on top;
3. a model of the whole field and amplitude and phase per antenna and
   polarization on every baseline, the weights following the amplitudes;
   again until the amplitudes found are small;
4. a model of the whole field again and a gain for every few channels,
   which corrects each antenna's own run of phase and amplitude across
   the band;
5. the image.

There is no amplitude solve against the cores alone. A model without the
extended emission gives the antennas of the central square, whose short
baselines see it, amplitudes above the others', and a usual clean made
afterwards takes that tilt for the shape of the sources and cannot give
it back (finding 125). What differs from antenna to antenna is found
against the whole field; how the array's short and long baselines stand to
each other is left as the calibrators set it.

Positive components are for the cores-only models, where they are what
breaks the degeneracy below. The whole-field models and the image are
usual cleans: a model that may not go negative ends above the data by 15
to 40% around every resolved source (finding 122), and an amplitude solved
against it scales the short baselines up.

A cores-only model first is the point of the recipe. Self-calibration
cannot tell an offset antenna from a sky in which every source has the same
lopsided shape; only a model that insists on compact, positive sources
can, and a model cleaned deep does not insist.

WSClean and CASA live in different containers, so this script, which needs
neither, runs each step with the command given for it. It works on a copy
of the Measurement Set and never writes to the one it is given. A run that
is interrupted continues from the last finished step when the same command
is given again.

Example::

    python3 selfcal_recipe.py --target G0212-0335 \\
        --ms wsclean_test/G0212-0335_scan15.ms \\
        --imsize 5120 --cell 1.0 --channel-range 51 768
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
import subprocess
import time
from pathlib import Path


CONTAINERS = "/software/astro/containers"
DEFAULT_WSCLEAN = f"singularity exec {CONTAINERS}/oxkat-latest.sif wsclean"
# Tried in turn when --casa is not given; the first that answers is used.
CASA_CANDIDATES = tuple(
    f"singularity exec {CONTAINERS}/{name} python3"
    for name in (
        "casa-modular-py3.12-latest.sif",
        "casa-modular-latest.sif",
        "casa-modular-py3.10-latest.sif",
        "casa-6.5.0-modular.sif",
    )
)
REPORT_NAME = "recipe_report.json"
AMPLITUDE_CORES_DEPTH = 0.5
# When the amplitudes of a solve differ between the inner and outer halves
# of the array by more than this, their run with distance from the array's
# middle is divided out (findings 147, 148). It is the model's: a cores-only
# model, and a usual clean too, can hold the wrong share of a resolved
# source's flux on short baselines against long ones. Every amplitude solve
# and the bandpass round are corrected, so that small runs do not add up
# from round to round.
AMPLITUDE_DETREND = 0.05

# The rounds, in order. ``most`` and ``until`` repeat a round until the
# corrections it finds are small: a phase in degrees, an amplitude as its
# distance from one. The numbers are those of G0212-0335, where the phases
# found went 82, 16, 7 degrees and the amplitudes against the whole field
# were all within 7% of one.
DEFAULT_PLAN = (
    {"id": "phase", "model": {"core": "phase"},
     "solve": {"calmode": "p"},
     "most": 4, "until": {"phase_deg": 10.0}},
    {"id": "delay", "model": {"core": "phase"},
     "solve": {"calmode": "p", "delays": True}},
    {"id": "amplitude", "model": {"full": True},
     "solve": {"calmode": "ap"},
     "most": 3, "until": {"phase_deg": 5.0, "amplitude": 0.10}},
    # Each antenna's own bandpass, last: on G0212-0335 one antenna was at
    # twice its amplitude over a third of the band (finding 123).
    {"id": "bandpass", "model": {"full": True},
     "solve": {"bandpass": True}},
)


def load_sibling(filename, name):
    spec = importlib.util.spec_from_file_location(
        name, Path(SCRIPT_DIR) / filename
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


IMAGER = load_sibling("experiment_wsclean_image.py", "selfcal_recipe_imager")


def robust_label(robust):
    return f"robust_{robust:g}".replace("-", "m").replace(".", "p")


def rounds_of(stage):
    """The names of a stage's rounds: one, or numbered when it may repeat."""

    most = int(stage.get("most", 1))
    if most == 1:
        return [stage["id"]]
    return [f"{stage['id']}{turn}" for turn in range(1, most + 1)]


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def image_options(args, imsize=None):
    """What every WSClean step shares: grid, channels, resources."""

    options = [
        "--outdir", str(args.outdir),
        "--target", args.target,
        "--imsize", str(int(imsize or args.imsize)),
        "--cell", f"{float(args.cell):g}",
        "--robust", f"{float(args.robust):g}",
        "--channels-out", str(int(args.channels_out)),
        "--threads", str(int(args.threads)),
        "--memory-gb", f"{float(args.memory_gb):g}",
        "--wsclean", args.wsclean,
    ]
    if args.parallel_deconvolution is not None:
        options += [
            "--parallel-deconvolution", str(int(args.parallel_deconvolution)),
        ]
    if args.channel_range:
        options += [
            "--channel-range",
            str(int(args.channel_range[0])), str(int(args.channel_range[1])),
        ]
    else:
        options.append("--all-channels")
    return options


def deep_image_options(args):
    """The deconvolution that gave 27 to 30 uJy/beam on G0212-0335.

    Multi-scale, with the mask on one noise value for the image: with a
    local noise map WSClean 3.1.0 aborts, and on a dirty image that map
    freezes the mask after a few hundred components.
    """

    return [
        "--scales", *[str(int(scale)) for scale in args.scales],
        "--no-local-rms",
        "--auto-threshold", f"{float(args.auto_threshold):g}",
        "--nmiter", str(int(args.nmiter)),
        "--overwrite",
    ]


def model_command(args, ms, workdir, model):
    """The command that leaves a model of the sky in MODEL_DATA of ``ms``."""

    command = [
        sys.executable, str(Path(SCRIPT_DIR) / "experiment_wsclean_image.py"),
        "--ms", str(ms), "--workdir", str(workdir), *image_options(args),
    ]
    if "core" in model:
        # The amplitudes' cores go twice as deep as the phases': what a
        # source has around its peak is in them, and the baselines that an
        # amplitude is solved on see it. By hand: 6 and 3 mJy/beam.
        depth = AMPLITUDE_CORES_DEPTH if model["core"] == "amplitude" else 1.0
        if args.core_fraction:
            return command + [
                "--core-model-fraction",
                f"{float(args.core_fraction) * depth:g}",
            ]
        return command + [
            "--core-model", f"{float(args.core_threshold) * depth:g}",
        ]
    return command + deep_image_options(args) + ["--write-model"]


def image_command(args, ms, workdir, imsize=None):
    """The image of the calibrated data; nothing is written to them."""

    return [
        sys.executable, str(Path(SCRIPT_DIR) / "experiment_wsclean_image.py"),
        "--ms", str(ms), "--workdir", str(workdir),
        *image_options(args, imsize), *deep_image_options(args),
    ]


def solve_command(args, ms, solve):
    """The CASA step: solve against MODEL_DATA, apply, write a new copy."""

    command = shlex.split(args.casa) + [
        str(Path(SCRIPT_DIR) / "experiment_phase_offsets.py"),
        "--ms", str(ms),
        "--suffix", "_solved",
        "--overwrite",
        "--refant", args.refant,
    ]
    detrend = float(getattr(args, "amplitude_detrend", AMPLITUDE_DETREND) or 0)
    if solve.get("bandpass"):
        return command + [
            "--bandpass", str(int(args.bandpass_channels)),
            "--bandpass-amplitude",
        ] + (["--detrend-radius", f"{detrend:g}"] if detrend > 0 else []) + (
            ["--keep-unsolved"] if args.keep_unsolved else []
        )
    command += [
        "--calmode", solve.get("calmode", "p"),
        "--solint", solve.get("solint", "inf"),
        "--gaintype", solve.get("gaintype", "G"),
    ]
    if solve.get("delays"):
        command.append("--delays")
    if solve.get("uvrange"):
        command += ["--uvrange", solve["uvrange"]]
    if solve.get("calmode") == "ap" and detrend > 0:
        command += ["--detrend-radius", f"{detrend:g}"]
    if args.keep_unsolved:
        # A solve towards one source, or on what is left of a field: an
        # amplitude that is not believed flags none of the field's data.
        command.append("--keep-unsolved")
    return command


def solved_products(ms):
    """Where the solve step leaves its copy of the data and its report."""

    ms = Path(ms)
    return (
        ms.with_name(ms.stem + "_solved" + ms.suffix),
        ms.with_name(ms.stem + "_solved.json"),
    )


# ---------------------------------------------------------------------------
# Checks before anything long is started
# ---------------------------------------------------------------------------

CASA_CHECK = (
    "import sys, casatasks, casatools; "
    "t = casatools.table(); t.open(sys.argv[1]); n = t.nrows(); t.close(); "
    "print('CASA_OK', "
    "getattr(casatasks, 'version_string', lambda: 'unknown')(), n)"
)


def casa_answers(casa, ms):
    """Whether a command runs a Python that has CASA and can open ``ms``.

    Opening the Measurement Set, read-only, also shows that its folder is
    visible from inside a container, which is otherwise found out only
    when the first solve starts.
    """

    try:
        done = subprocess.run(
            shlex.split(casa) + ["-c", CASA_CHECK, str(ms)],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            timeout=900,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return False, str(exc)
    for line in (done.stdout or "").splitlines():
        if line.startswith("CASA_OK"):
            return True, line.split()[1]
    return False, (done.stdout or "").strip()[-600:]


def preflight(args, ms):
    """Find out now, not in the third hour, that a tool does not answer.

    Returns the problems found and what was learnt; sets ``args.casa`` to
    the command that answered when none was given.
    """

    problems, found = [], {}
    words = shlex.split(args.wsclean)
    if not words or shutil.which(words[0]) is None:
        problems.append(
            f"WSClean: '{words[0] if words else ''}' is not on PATH. This "
            "script is started from the terminal that has singularity."
        )
    else:
        capabilities = IMAGER.wsclean_capabilities(words)
        if capabilities["version"] is None and not capabilities["options"]:
            problems.append(
                f"WSClean: '{args.wsclean} --version' did not answer as "
                "WSClean. Check --wsclean."
            )
        else:
            found["wsclean_version"] = capabilities["version"]

    candidates = [args.casa] if args.casa else list(CASA_CANDIDATES)
    failures = []
    for casa in candidates:
        words = shlex.split(casa)
        if not words or shutil.which(words[0]) is None:
            failures.append(
                f"'{casa}': '{words[0] if words else ''}' is not on PATH"
            )
            continue
        ok, detail = casa_answers(casa, ms)
        if ok:
            args.casa = casa
            found["casa"] = casa
            found["casa_version"] = detail
            break
        failures.append(f"'{casa}': {detail}")
    else:
        problems.append(
            "CASA: no command ran a Python that imports casatasks and opens "
            f"{ms}. Pass --casa with one that does. Tried:\n    "
            + "\n    ".join(failures)
        )
    return problems, found


# ---------------------------------------------------------------------------
# What each step found
# ---------------------------------------------------------------------------

def read_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def model_findings(report):
    """From a cores-only model: how lopsided the data were, and its flux."""

    report = report or {}
    asymmetry = report.get("dirty_asymmetry") or {}
    return {
        "lopsided_peak": asymmetry.get("antisymmetric_peak"),
        "lopsided_rms": asymmetry.get("antisymmetric_rms"),
        "model_flux_Jy": (report.get("model") or {}).get("total_flux_Jy"),
    }


def solve_findings(report):
    """From a solve: how large the corrections were, and on which antenna."""

    report = report or {}
    if "statistics" in report:
        return bandpass_findings(report)
    phases, steps, amplitudes = [], [], []
    for entry in report.get("series") or []:
        for pol in entry.get("polarizations") or []:
            if "mean_phase_deg" not in pol:
                continue
            # The largest single solution where the solve gives it: with a
            # solution per scan, or several within one, a mean hides what
            # one of them corrected.
            phases.append((
                pol.get("largest_phase_deg", abs(pol["mean_phase_deg"])),
                entry["antenna"],
            ))
            amplitude = pol.get("mean_amplitude", 1.0)
            amplitudes.append(amplitude)
            steps.append((
                pol.get("largest_amplitude_step", abs(amplitude - 1.0)),
                entry["antenna"],
            ))
    delays = [
        (abs(delay), entry["antenna"])
        for entry in report.get("delays") or []
        for delay in entry["delay_ns"] if delay is not None
    ]
    found = {
        "calmode": (report.get("solve") or {}).get("calmode"),
        "antennas_solved": (report.get("summary") or {}).get(
            "solved_antennas"),
        "applied": bool(report.get("corrected_measurement_set")),
    }
    rejected = [
        f"delay of {entry['antenna']}"
        for entry in report.get("rejected_delays") or []
    ] + [
        f"amplitude of {name}"
        for name in report.get("rejected_amplitudes") or []
    ]
    if rejected:
        found["rejected"] = rejected
    trend = report.get("radial_trend") or {}
    if trend.get("divided_out"):
        found["radial_trend"] = [trend["inner_median"], trend["outer_median"]]
    found.update(rescaled(report))
    moving = [
        (pol.get("amplitude_rms", 0.0), pol.get("phase_rms_deg", 0.0),
         entry["antenna"])
        for entry in report.get("series") or []
        if entry.get("solutions", 1) > 1
        for pol in entry.get("polarizations") or []
        if "mean_phase_deg" in pol
    ]
    if moving:
        found["time_amplitude"], _, found["time_amplitude_antenna"] = max(
            moving
        )
        found["time_phase_deg"] = max(item[1] for item in moving)
    if phases:
        found["phase_deg"], found["phase_antenna"] = max(phases)
    if delays:
        found["delay_ns"], found["delay_antenna"] = max(delays)
    if found["calmode"] == "ap" and steps:
        found["amplitude"], found["amplitude_antenna"] = max(steps)
        found["amplitude_range"] = [min(amplitudes), max(amplitudes)]
    return found


def rescaled(report):
    """The median amplitudes a solve was brought back from, when far from one."""

    medians = report.get("amplitude_medians") or []
    if any(abs(float(median) - 1.0) > 0.02 for median in medians):
        return {"rescaled": [float(median) for median in medians]}
    return {}


def bandpass_findings(report):
    """From a solve per channel bin: how far the gains run across the band."""

    phases, amplitudes = [], []
    for row in report.get("statistics") or []:
        for pol in row.get("polarizations") or []:
            if "phase_rms_deg" not in pol:
                continue
            phases.append((pol["phase_rms_deg"], row["antenna"]))
            if pol.get("amplitude_rms") is not None:
                amplitudes.append((pol["amplitude_rms"], row["antenna"]))
    found = {
        "calmode": "bandpass",
        "antennas_solved": len({name for _, name in phases}),
        "applied": bool(report.get("corrected_measurement_set")),
    }
    if report.get("rejected_bins"):
        found["rejected"] = [
            f"{count} bins of {name}"
            for name, count in sorted(report["rejected_bins"].items())
        ]
    trend = report.get("radial_trend") or {}
    if trend.get("divided_out"):
        found["radial_trend"] = [trend["inner_median"], trend["outer_median"]]
    found.update(rescaled(report))
    if phases:
        found["band_phase_deg"], found["band_phase_antenna"] = max(phases)
    if amplitudes:
        found["band_amplitude"], found["band_amplitude_antenna"] = max(
            amplitudes
        )
    return found


def image_findings(report):
    """From an image: its noise, and what is left beside the sources."""

    report = report or {}
    image = report.get("stokes_i") or {}
    asymmetry = report.get("source_asymmetry") or {}
    return {
        "image": image.get("image"),
        # Stokes V where the data have it; none for XX and YY alone.
        "data_noise_Jy_beam": IMAGER.data_noise_of(report).get(
            "robust_rms_Jy_beam"),
        "image_rms_Jy_beam": image.get("robust_rms_Jy_beam"),
        "quietest_tenth_Jy_beam": (image.get("grid") or {}).get(
            "quietest_tenth_rms"),
        "peak_Jy_beam": image.get("max_Jy_beam"),
        "most_negative_Jy_beam": image.get("min_Jy_beam"),
        "negative_beside_sources": asymmetry.get("deepest_negative_median"),
        "lopsided_peak": asymmetry.get("antisymmetric_peak"),
    }


def tidy_images(folder, final=False):
    """Remove from a finished step's images what nothing will read again.

    A step at 5120 pixels and twelve sub-bands writes 6 GB: five images for
    each sub-band and for the band. Of a round that is done, the band's
    image, model and residual are kept, to look at. Of the final image the
    sub-bands' images and models are kept as well: the models are what a
    source is peeled with, and the images carry the frequency and weight
    the primary-beam correction needs.
    """

    folder = Path(folder)
    removed = 0
    for path in folder.glob("*.fits"):
        name = path.name
        sub_band = "-MFS-" not in name
        keep = (
            name.endswith(("-MFS-image.fits", "-MFS-model.fits",
                           "-MFS-residual.fits"))
            or (final and sub_band and "_I-" in name
                and name.endswith(("-image.fits", "-model.fits")))
        )
        if not keep:
            path.unlink()
            removed += 1
    return removed


def settled(solve, until):
    """Whether a solve found corrections small enough to stop repeating.

    A quantity the solve did not report counts as not settled: a round is
    not ended on a number that is missing.
    """

    if not until:
        return False
    return all(
        solve.get(key) is not None and solve[key] <= limit
        for key, limit in until.items()
    )


def percent(value):
    return "?" if value is None else f"{100 * value:.1f}%"


def micro(value):
    return "?" if value is None else f"{value * 1e6:.1f}"


def image_line(name, image):
    line = (
        f"{name}: image {micro(image['image_rms_Jy_beam'])} uJy/beam "
        f"({micro(image['quietest_tenth_Jy_beam'])} in its quietest tenth)"
    )
    if image.get("data_noise_Jy_beam") is not None:
        line += f" against {micro(image['data_noise_Jy_beam'])} for the data"
    if image.get("peak_Jy_beam") is not None:
        line += f"; peak {image['peak_Jy_beam'] * 1e3:.1f} mJy/beam"
    return line + (
        "; negative beside the sources "
        f"{percent(image['negative_beside_sources'])} of their peak."
    )


def reading(report):
    """The run so far in a few lines: each round's correction, the images."""

    lines = []
    for entry in report["rounds"]:
        model, solve = entry.get("model") or {}, entry.get("solve")
        if (entry.get("image") or {}).get("image_rms_Jy_beam"):
            lines.append(image_line(entry["id"] + " model", entry["image"]))
        if solve is None:
            continue
        parts = []
        if model.get("lopsided_peak") is not None:
            parts.append(
                "lopsided part of the dirty image "
                f"{percent(model['lopsided_peak'])} of the peaks"
            )
        if "delay_ns" in solve:
            parts.append(
                f"largest delay {solve['delay_ns']:.2f} ns "
                f"({solve['delay_antenna']})"
            )
        if "phase_deg" in solve:
            parts.append(
                f"largest phase {solve['phase_deg']:.0f} deg "
                f"({solve['phase_antenna']})"
            )
        if "amplitude_range" in solve:
            low, high = solve["amplitude_range"]
            parts.append(f"amplitudes {low:.2f} to {high:.2f}")
        if "time_amplitude" in solve:
            parts.append(
                f"within the scan, amplitudes move by up to "
                f"{100 * solve['time_amplitude']:.1f}% rms "
                f"({solve['time_amplitude_antenna']}) and phases by up to "
                f"{solve['time_phase_deg']:.1f} deg"
            )
        if "band_phase_deg" in solve:
            parts.append(
                f"across the band, phases up to {solve['band_phase_deg']:.1f} "
                f"deg rms ({solve['band_phase_antenna']})"
            )
        if "band_amplitude" in solve:
            parts.append(
                f"amplitudes up to {100 * solve['band_amplitude']:.1f}% rms "
                f"({solve['band_amplitude_antenna']})"
            )
        if solve.get("rejected"):
            parts.append(
                "not believed and not applied: " + ", ".join(solve["rejected"])
            )
        if solve.get("radial_trend"):
            inner, outer = solve["radial_trend"]
            parts.append(
                f"amplitudes ran with distance from the array's middle "
                f"(inner half {inner:.2f}, outer half {outer:.2f}), the "
                "model's doing: divided out"
            )
        if solve.get("rescaled"):
            parts.append(
                "median amplitude "
                + "/".join(f"{value:.2f}" for value in solve["rescaled"])
                + " set back to 1, so the flux scale stays the calibrators'"
            )
        if not solve.get("applied"):
            parts.append("NOTHING APPLIED: no antenna had a solution")
        lines.append(f"{entry['id']}: " + "; ".join(parts) + ".")
    final = report.get("final") or {}
    if final.get("image_rms_Jy_beam"):
        lines.append(image_line("final", final))
    before = report.get("before") or {}
    if before.get("image_rms_Jy_beam"):
        lines.insert(0, image_line("before", before))
    return lines


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

def print_plan(args, plan, source, root):
    """What a run would do, without doing it."""

    args.casa = args.casa or CASA_CANDIDATES[0]
    ms = root / "round0.ms"
    print(f"copy {source} -> {ms}")
    if args.image_first:
        print(shlex.join(image_command(args, ms, root / "images" / "before")))
    for stage in plan:
        names = rounds_of(stage)
        if len(names) > 1:
            limits = ", ".join(
                f"{key} <= {limit:g}"
                for key, limit in (stage.get("until") or {}).items()
            )
            print(f"# {stage['id']}: up to {len(names)} times, until {limits}")
        images = root / "images" / names[0]
        print(shlex.join(model_command(args, ms, images, stage["model"])))
        print(shlex.join(solve_command(args, ms, stage["solve"])))
    print(shlex.join(image_command(
        args, ms, root / "images" / "final", args.final_imsize
    )))


def run(args, plan=DEFAULT_PLAN):
    source = Path(str(args.ms).rstrip("/"))
    if not source.is_dir():
        raise SystemExit(f"No Measurement Set at {source}.")
    root = Path(args.workdir) / args.target
    if IMAGER.inside(root, args.outdir):
        raise SystemExit(
            f"{root} is inside the pipeline's products ({args.outdir}). "
            "Give --workdir a folder of its own."
        )
    if IMAGER.inside(source, root):
        raise SystemExit(f"{source} is inside the work folder {root}.")
    if args.dry_run:
        print_plan(args, plan, source, root)
        return None

    restart = args.restart and root.exists()
    if restart and not (root / REPORT_NAME).is_file():
        raise SystemExit(
            f"{root} holds no {REPORT_NAME}: it is not an earlier run of "
            "this script, and --restart will not remove it."
        )

    root.mkdir(parents=True, exist_ok=True)
    # A second copy of this run would remove the first one's data from
    # under it, so the folder is held before anything in it is touched.
    with IMAGER.claimed(root) as claim:
        if restart:
            for path in root.iterdir():
                if path == claim:
                    continue
                if path.is_dir() and not path.is_symlink():
                    shutil.rmtree(path)
                else:
                    path.unlink()
        return follow(args, plan, source, root)


def follow(args, plan, source, root):
    """Carry the plan out in a work folder this run holds."""

    report_path = root / REPORT_NAME
    label = robust_label(args.robust)
    report = read_json(report_path) or {
        "target": args.target, "input": str(source),
        "done": [], "times": {}, "rounds": [],
    }
    if report.get("input") != str(source):
        raise SystemExit(
            f"{root} holds a run on {report.get('input')}, not {source}. "
            "Give --workdir another folder, or --restart to discard it."
        )

    def save():
        report["reading"] = reading(report)
        report["updated"] = now()
        root.mkdir(parents=True, exist_ok=True)
        scratch = report_path.with_suffix(".tmp")
        with open(scratch, "w") as handle:
            json.dump(report, handle, indent=2)
        scratch.replace(report_path)

    def done(key):
        return key in report["done"]

    def mark(key):
        report["done"].append(key)
        report["times"][key] = now()
        report["running"] = None
        save()

    def record(name):
        for entry in report["rounds"]:
            if entry["id"] == name:
                return entry
        report["rounds"].append({"id": name})
        return report["rounds"][-1]

    def step(name, what, command):
        print(f"\n=== {name}: {what} ({now()}) ===", flush=True)
        report["running"] = {"step": f"{name}: {what}", "since": now()}
        save()
        log = root / "logs" / f"{name}_{what.split()[0]}.log"
        try:
            IMAGER.run_command(command, log)
        except IMAGER.CommandFailed as failure:
            report["running"] = None
            report["failed"] = {"step": f"{name}: {what}", "error": str(failure)}
            save()
            raise SystemExit(
                f"\nStep '{name}' ({what}) failed: {failure}.\nEverything "
                "before it is kept; the same command continues from this "
                "step."
            ) from None
        report.pop("failed", None)

    problems, found = preflight(args, source)
    if problems:
        raise SystemExit("Cannot start:\n  " + "\n  ".join(problems))
    version = found.get("wsclean_version")
    print(
        "WSClean "
        + (f"{version[0]}.{version[1]}" if version else "(version unknown)")
        + f"; CASA {found['casa_version']} ({found['casa']}).",
        flush=True,
    )
    report["versions"] = found
    (root / "logs").mkdir(parents=True, exist_ok=True)

    count = 0
    current = root / "round0.ms"
    if not done("copy"):
        if current.exists():
            shutil.rmtree(current)
        print(f"\nCopying {source} to {current} ...", flush=True)
        shutil.copytree(source, current)
        mark("copy")

    if args.image_first and not done("before:image"):
        # What the data give as they arrive, to set the result against.
        step("before", f"image of {current.name}",
             image_command(args, current, root / "images" / "before"))
        report["before"] = image_findings(read_json(
            root / "images" / "before" / args.target / f"wsclean_{label}.json"
        ))
        mark("before:image")

    for stage in plan:
        for name in rounds_of(stage):
            entry = record(name)
            images = root / "images" / name
            following = root / f"round{count + 1}.ms"

            if not done(f"{name}:model"):
                step(name, f"model for {current.name}",
                     model_command(args, current, images, stage["model"]))
                if "core" in stage["model"]:
                    entry["model"] = model_findings(read_json(
                        images / args.target / f"coremodel_{label}.json"
                    ))
                else:
                    entry["image"] = image_findings(read_json(
                        images / args.target / f"wsclean_{label}.json"
                    ))
                mark(f"{name}:model")

            if not done(f"{name}:solve"):
                step(name, f"solve on {current.name}",
                     solve_command(args, current, stage["solve"]))
                solved, solve_report = solved_products(current)
                entry["solve"] = solve_findings(read_json(solve_report))
                if not solved.is_dir():
                    save()
                    raise SystemExit(
                        f"\nThe solve of '{name}' applied nothing: no antenna "
                        f"had a solution. See {root / 'logs'}."
                    )
                if following.exists():
                    shutil.rmtree(following)
                solved.rename(following)
                mark(f"{name}:solve")

            if not args.keep_intermediate and current.exists():
                # This round's input has served: with its model and its
                # corrected columns it is three times the size of the data.
                shutil.rmtree(current)
            if args.tidy:
                tidy_images(images / args.target)
            count += 1
            current = following
            if settled(entry.get("solve") or {}, stage.get("until")):
                break

    report["calibrated_measurement_set"] = str(current)
    if not done("final:image"):
        step("final", f"image of {current.name}",
             image_command(args, current, root / "images" / "final",
                           args.final_imsize))
        report["final"] = image_findings(read_json(
            root / "images" / "final" / args.target / f"wsclean_{label}.json"
        ))
        mark("final:image")
    if args.tidy:
        tidy_images(root / "images" / "final" / args.target, final=True)
        tidy_images(root / "images" / "before" / args.target)

    save()
    print("\n" + "\n".join(report["reading"]))
    print(f"\nCalibrated data: {current}")
    if (report.get("final") or {}).get("image"):
        print(f"Image: {report['final']['image']}")
    print(f"Report: {report_path}")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", required=True)
    parser.add_argument("--ms", required=True,
                        help="the Measurement Set to calibrate; it is copied, "
                        "never written to")
    parser.add_argument("--workdir", default="selfcal_recipe")
    parser.add_argument("--outdir", default="uGMRT_Products",
                        help="the pipeline's products, which are not touched")
    parser.add_argument("--imsize", type=int, required=True)
    parser.add_argument("--final-imsize", type=int,
                        help="pixels of the final image, if larger than "
                        "the models need to be: the solves want the sources "
                        "that matter, the image the whole field")
    parser.add_argument("--cell", type=float, required=True, help="arcsec")
    parser.add_argument("--channel-range", type=int, nargs=2,
                        metavar=("FIRST", "END"),
                        help="channels that hold data, END exclusive; every "
                        "channel if not given")
    parser.add_argument("--robust", type=float, default=0.0)
    parser.add_argument("--core-threshold", type=float, default=6.0,
                        metavar="MJY",
                        help="the cores-only model for phases and delays "
                        "holds what is brighter than this (default 6 "
                        "mJy/beam, a quarter of the brightest peak of "
                        "G0212-0335 before correction)")
    parser.add_argument("--core-fraction", type=float, default=0.0,
                        metavar="FRACTION",
                        help="make the cores-only models down to this "
                        "fraction of the field's own brightest pixel "
                        "instead of --core-threshold, for a field whose "
                        "brightness is not known beforehand (0.13 is what "
                        "6 mJy/beam was for G0212-0335)")
    parser.add_argument("--time-phase-solint", metavar="INTERVAL",
                        help="after the phase rounds per scan, solve "
                        "phases at this interval within each scan, for "
                        "example 1min: for data that have had no "
                        "self-calibration yet")
    parser.add_argument("--amplitude-cores", action="store_true",
                        help="after the delays, solve amplitude and phase "
                        "against a cores-only model on the longer "
                        "baselines, until what is found is small: for data "
                        "whose amplitudes have had no self-calibration, "
                        "where a usual clean takes an antenna's error for "
                        "the shape of the sources (finding 131)")
    parser.add_argument("--amplitude-detrend", type=float,
                        default=AMPLITUDE_DETREND, metavar="FRACTION",
                        help="divide out the amplitudes' run with distance "
                        "from the array's middle when the inner and outer "
                        "halves differ by more than this (default "
                        f"{AMPLITUDE_DETREND:g}; 0 leaves them as solved)")
    parser.add_argument("--amplitude-cores-uvrange", default=">3klambda",
                        metavar="UVRANGE",
                        help="the baselines of that solve (default "
                        ">3klambda): the cores-only model lacks the "
                        "extended emission that shorter ones see")
    parser.add_argument("--time-solint", metavar="INTERVAL",
                        help="add a last round that solves amplitude and "
                        "phase at this interval within each scan, for "
                        "example 2min")
    parser.add_argument("--channels-out", type=int, default=12,
                        help="sub-bands the models and images are made in "
                        "(default 12: with six the model's steps in "
                        "frequency cost 2%% of the noise)")
    parser.add_argument("--bandpass-channels", type=int, default=30,
                        help="channels per gain in the bandpass round "
                        "(default 30)")
    parser.add_argument("--keep-unsolved", action="store_true",
                        help="in the bandpass round, leave what could not "
                        "be solved unflagged: for a solve towards one "
                        "source, whose data are the whole field's")
    parser.add_argument("--whole-field-only", action="store_true",
                        help="no rounds against cores-only models: only "
                        "the amplitude-and-phase rounds and the bandpass "
                        "round, against images of the whole field; for "
                        "data that are self-calibrated already, where a "
                        "cores-only model lacks most of what the short "
                        "baselines see and a phase solved against it is "
                        "wrong for them (finding 132)")
    parser.add_argument("--phase-amplitude-only", action="store_true",
                        help="solve only a phase and an amplitude per "
                        "antenna and scan: no delays, no bandpass; for a "
                        "source too faint for more")
    parser.add_argument("--phase-only", action="store_true",
                        help="no amplitude solved at all: the rounds "
                        "against the whole field solve phases only, and "
                        "the cores-only amplitude rounds and the bandpass "
                        "round are left out; for a field too faint for "
                        "amplitudes")
    parser.add_argument("--image-first", action="store_true",
                        help="image the data before any round as well, to "
                        "compare the result with")
    parser.add_argument("--scales", type=int, nargs="+",
                        default=[0, 5, 15, 30])
    parser.add_argument("--auto-threshold", type=float, default=1.0)
    parser.add_argument("--nmiter", type=int, default=20)
    parser.add_argument("--refant", default="C00,W03,E05,E02,E06")
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--memory-gb", type=float, default=20.0)
    parser.add_argument("--parallel-deconvolution", type=int,
                        metavar="PIXELS",
                        help="clean the images in sub-images of this size at "
                        "once; nothing given leaves it to "
                        "experiment_wsclean_image.py (off before WSClean "
                        "3.2)")
    parser.add_argument("--wsclean", default=DEFAULT_WSCLEAN,
                        help="the command that runs WSClean")
    parser.add_argument("--casa",
                        help="the command that runs a Python with casatasks; "
                        "the CASA containers on ilifu are tried if not given")
    parser.add_argument("--tidy", action="store_true",
                        help="remove the sub-band images of each finished "
                        "round, and of the final image all but what "
                        "peeling and the primary-beam correction read")
    parser.add_argument("--keep-intermediate", action="store_true",
                        help="keep every round's copy of the data")
    parser.add_argument("--restart", action="store_true",
                        help="discard an earlier run in this work folder")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the steps and stop")
    args = parser.parse_args(argv)
    plan = DEFAULT_PLAN
    if args.whole_field_only:
        for name in ("time_phase_solint", "amplitude_cores"):
            if getattr(args, name):
                parser.error(
                    f"--{name.replace('_', '-')} is a round against a "
                    "cores-only model, which --whole-field-only leaves out"
                )
        # Cores-only models are what first breaks the degeneracy between
        # a lopsided source and an offset one, on data whose phases are
        # tens of degrees out. On data calibrated already they cost more
        # than they give: the brightest sources are what holds them, and
        # with those out the model has a tenth of the field's flux.
        plan = tuple(stage for stage in plan if stage["model"].get("full"))
    if args.phase_amplitude_only:
        # A delay or a gain per channel bin has too many unknowns for a
        # faint source: on an 11 mJy source the delay round gave one antenna
        # 18 ns and the bandpass round amplitudes that ran by 53%.
        plan = tuple(
            stage for stage in plan if stage["id"] in ("phase", "amplitude")
        )
    if args.time_phase_solint:
        # Data that have had no self-calibration yet carry phases that
        # drift within a scan, which one phase per scan cannot follow.
        # After the rounds per scan, the cores-only model again and a
        # phase at this interval, until what it finds is small.
        index = [stage["id"] for stage in plan].index("phase") + 1
        plan = plan[:index] + ({
            "id": "phase_time", "model": {"core": "phase"},
            "solve": {"calmode": "p", "solint": args.time_phase_solint},
            "most": 3, "until": {"phase_deg": 5.0},
        },) + plan[index:]
    if args.amplitude_cores:
        # A usual clean reproduces whatever the data are: an antenna that
        # is too high or too low looks to it like negative components
        # beside every source, and a solve against that model returns one.
        # Only a model that cannot hold those, positive cores, shows the
        # amplitudes. It goes before the whole field's rounds.
        index = [stage["id"] for stage in plan].index("amplitude")
        plan = plan[:index] + ({
            "id": "amplitude_cores", "model": {"core": "amplitude"},
            "solve": {"calmode": "ap",
                      "uvrange": args.amplitude_cores_uvrange},
            "most": 3, "until": {"phase_deg": 5.0, "amplitude": 0.10},
        },) + plan[index:]
    if args.time_solint:
        plan += ({
            "id": "time", "model": {"full": True},
            "solve": {"calmode": "ap", "solint": args.time_solint},
        },)
    if args.phase_only:
        plan = phase_only(plan)
    run(args, plan)
    return 0


def phase_only(plan):
    """The plan without amplitudes: whole-field rounds solve phases only.

    In the faintest Band 5 fields the amplitude rounds lowered the compact
    sources' peaks round by round while the Stokes V noise fell with them:
    G0231-0533 (one 12-minute scan, brightest source 4 mJy/beam) lost 14%
    of its brightest compact source's peak between the last phase-only
    round and its final image, G0216-0604's target stage 5 to 8%
    (2026-10-24).
    """

    renamed = {"amplitude": "phase_field", "time": "phase_field_time"}
    kept = []
    for stage in plan:
        solve = stage["solve"]
        if solve.get("bandpass"):
            continue
        if solve.get("calmode") == "ap":
            if "core" in stage["model"]:
                continue
            stage = dict(
                stage, id=renamed.get(stage["id"], stage["id"] + "_phase"),
                solve=dict(solve, calmode="p"),
                until={"phase_deg": 5.0},
            )
        kept.append(stage)
    return tuple(kept)


if __name__ == "__main__":
    raise SystemExit(main())
