#!/usr/bin/env python3
"""Explain where calibrator data were lost: by scan, by field, by channel.

CASA's gain solvers report "Insufficient unflagged antennas to proceed with
this solve" and "Found no unflagged data" wherever flags leave too little data.
Those messages say *that* data are missing, not *why*. This tool reads the
products the pipeline already wrote and answers the why:

* for each scan, how much was flagged on arrival and after each flagging step,
  and which step removed most of it;
* for each run of fully flagged channels, whether it is a deliberate mask
  (receiver edge or out-of-band) or a data-driven RFI flag;
* optionally, with CASA's ``casatools`` importable, how many usable
  calibration solutions each scan has in the current and round-1 tables.

It is read-only and never modifies a Measurement Set or a report.

Examples::

    python3 diagnose_flag_losses.py --outdir Products
    python3 diagnose_flag_losses.py --outdir Products --time 2020/11/22/14:05:38
"""

import argparse
import json
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
MJD_EPOCH = datetime(1858, 11, 17, tzinfo=timezone.utc)
LOSS_THRESHOLD = 0.90
FULLY_FLAGGED = 0.99


def read_json(path, default=None):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError, TypeError):
        return default


def configured_outdir():
    """The pipeline's configured output directory (config.yaml ``outdir``)."""

    try:
        import importlib.util

        spec = importlib.util.spec_from_file_location(
            "diagnose_flag_losses_config",
            SCRIPT_DIR / "config.py",
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return str(module.PIPELINE_RUN_CONFIG.get("outdir") or "Products")
    except Exception:
        return "Products"


def resolve_outdir(value, script_dir=None):
    """Find the output directory the way the pipeline's batch job uses it.

    A relative ``outdir`` is interpreted by the job from the pipeline
    directory, so look there as well as in the current directory.
    """

    path = Path(value).expanduser()
    if path.is_absolute():
        return path, [path]

    bases = [Path.cwd(), Path(script_dir) if script_dir else SCRIPT_DIR]
    searched = []
    for base in bases:
        candidate = (base / path).resolve()
        if candidate in searched:
            continue
        searched.append(candidate)
        if (candidate / "inspection" / "inspection_report.json").is_file():
            return candidate, searched
    return searched[0], searched


def mjd_seconds_to_utc(seconds):
    moment = MJD_EPOCH + timedelta(seconds=float(seconds))
    return moment.strftime("%Y/%m/%d/%H:%M:%S.") + f"{moment.microsecond // 100000}"


def casa_time_to_mjd_seconds(text):
    """Parse CASA log times such as ``2020/11/22/14:05:38.4``."""

    match = re.fullmatch(
        r"\s*(\d{4})/(\d{2})/(\d{2})/(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)\s*",
        str(text),
    )
    if not match:
        raise ValueError(f"Unrecognised CASA time: {text!r}")
    year, month, day, hour, minute = (int(value) for value in match.groups()[:5])
    second = float(match.group(6))
    moment = datetime(year, month, day, hour, minute, tzinfo=timezone.utc)
    return (moment - MJD_EPOCH).total_seconds() + second


def flagged_fraction(entry):
    try:
        total = float(entry["total"])
        flagged = float(entry["flagged"])
    except (KeyError, TypeError, ValueError):
        return None
    if total <= 0:
        return None
    return flagged / total


def group_fractions(summary, category):
    groups = (summary or {}).get(category, {})
    if not isinstance(groups, dict):
        return {}
    result = {}
    for key, entry in groups.items():
        fraction = flagged_fraction(entry)
        if fraction is not None:
            result[str(key)] = fraction
    return result


def flagging_stages(outdir):
    """Return (label, summary) for each saved point in the flagging history."""

    root = Path(outdir) / "flagging"
    stages = []

    for label, name in (
        ("on arrival", "preflag_summary.json"),
        ("after fixed masks", "iterative_baseline_summary.json"),
    ):
        summary = read_json(root / name)
        if isinstance(summary, dict):
            stages.append((label, summary))

    iterations = sorted(
        root.glob("iter_*_summary.json"),
        key=lambda path: int(re.search(r"iter_(\d+)_", path.name).group(1)),
    )
    for path in iterations:
        summary = read_json(path)
        if isinstance(summary, dict):
            number = re.search(r"iter_(\d+)_", path.name).group(1)
            stages.append((f"after auto-flag {number}", summary))

    summary = read_json(root / "postflag_summary.json")
    if isinstance(summary, dict):
        stages.append(("after extension", summary))

    return stages


def postcal_stages(outdir):
    """Per-group fractions around the corrected-data flagging pass, if saved."""

    stages = []
    for path in sorted((Path(outdir) / "calibration").glob(
        "*_postcal_flagging_report.json"
    )):
        report = read_json(path, {}) or {}
        breakdowns = report.get("group_flag_fractions") or {}
        if "after" in breakdowns:
            stages.append(("after post-cal flagging", breakdowns["after"]))
    return stages


def attribute_loss(fractions):
    """Name the step responsible for most of a group's flagged data.

    ``fractions`` is an ordered list of (label, fraction or None).
    """

    known = [(label, value) for label, value in fractions if value is not None]
    if not known:
        return "no data recorded"

    first_label, first_value = known[0]
    best_label, best_increase = first_label, first_value
    previous = first_value
    for label, value in known[1:]:
        increase = value - previous
        if increase > best_increase:
            best_label, best_increase = label, increase
        previous = value

    if best_label == first_label:
        return "already flagged on arrival"
    return best_label


def parse_spw_selection(selection):
    """Expand a CASA ``spw:chan~chan;chan`` selector into channel sets."""

    channels = {}
    for part in str(selection or "").split(","):
        part = part.strip()
        if ":" not in part:
            continue
        spw, ranges = part.split(":", 1)
        target = channels.setdefault(spw.strip(), set())
        for item in ranges.split(";"):
            item = item.strip()
            if not item:
                continue
            if "~" in item:
                low, high = (int(value) for value in item.split("~", 1))
                target.update(range(low, high + 1))
            else:
                target.add(int(item))
    return channels


def deliberate_channel_masks(outdir):
    """Channels flagged on purpose: receiver edges and out-of-band spectrum."""

    report = read_json(Path(outdir) / "flagging" / "edge_flagging_report.json", {})
    report = report or {}
    selections = [report.get("spw_string", "")]
    for key in ("guaranteed_science_range_mask", "percentage_edge_mask"):
        nested = report.get(key) or {}
        selections.append(nested.get("spw_selection", ""))

    masks = {}
    for selection in selections:
        for spw, channels in parse_spw_selection(selection).items():
            masks.setdefault(spw, set()).update(channels)
    return masks


def channel_fractions(outdir):
    """Final per-channel flagged fraction, keyed by SPW then channel."""

    statistics = read_json(Path(outdir) / "flagging" / "spw_statistics.json", {})
    by_channel = (statistics or {}).get("by_spw_channel") or {}
    entries = by_channel.get("spw:channel")
    if not isinstance(entries, dict):
        entries = next(
            (
                value
                for key, value in by_channel.items()
                if "channel" in str(key) and isinstance(value, dict)
            ),
            {},
        )

    result = {}
    for key, entry in entries.items():
        if ":" not in str(key):
            continue
        spw, channel = str(key).split(":", 1)
        fraction = flagged_fraction(entry)
        if fraction is None:
            continue
        try:
            result.setdefault(spw, {})[int(channel)] = fraction
        except ValueError:
            continue
    return result


def channel_runs(fractions, threshold=FULLY_FLAGGED):
    """Contiguous channel runs at or above ``threshold``."""

    runs = []
    start = previous = None
    for channel in sorted(fractions):
        if fractions[channel] >= threshold:
            if start is None:
                start = channel
            elif channel != previous + 1:
                runs.append((start, previous))
                start = channel
            previous = channel
        elif start is not None:
            runs.append((start, previous))
            start = previous = None
    if start is not None:
        runs.append((start, previous))
    return runs


def classify_run(run, deliberate, hanning_guard=1):
    """Is a flagged run a deliberate mask, an RFI flag, or both?

    Hanning smoothing in calibration flags each neighbour of a flagged
    channel, so a deliberate mask may appear ``hanning_guard`` channels wider.
    """

    low, high = run
    widened = set()
    for channel in deliberate:
        widened.update(range(channel - hanning_guard, channel + hanning_guard + 1))
    channels = range(low, high + 1)
    covered = sum(1 for channel in channels if channel in widened)
    if covered == len(channels):
        return "deliberate mask"
    if covered == 0:
        return "data-driven flags (RFI / auto-flagging)"
    return (
        f"mixed: {covered}/{len(channels)} channels deliberate, "
        "rest data-driven"
    )


def channel_frequency_mhz(inspection, spw, channel):
    info = ((inspection or {}).get("spws") or {}).get(str(spw)) or {}
    frequencies = info.get("channel_frequencies_MHz") or []
    if 0 <= channel < len(frequencies):
        return float(frequencies[channel])
    return None


def field_names(outdir, inspection):
    names = {}
    for item in read_json(Path(outdir) / "inspection" / "fields.json", []) or []:
        if isinstance(item, dict) and "id" in item and "name" in item:
            names[int(item["id"])] = str(item["name"])

    roles = {}
    for key, label in (
        ("flux_calibrators", "flux cal"),
        ("bandpass_calibrators", "bandpass cal"),
        ("phase_calibrators", "phase cal"),
        ("target_fields", "target"),
    ):
        for name in (inspection or {}).get(key, []) or []:
            roles.setdefault(str(name), []).append(label)
    return names, roles


def caltable_scan_coverage(caltable, scan_spans, table_factory):
    """Usable solutions per scan in one calibration table (needs casatools)."""

    tool = table_factory()
    opened = False
    try:
        tool.open(str(caltable))
        opened = True
        times = list(tool.getcol("TIME"))
        flags = tool.getcol("FLAG")
    finally:
        if opened:
            tool.close()

    shape = getattr(flags, "shape", ())
    rows = shape[-1] if shape else len(times)
    coverage = {}
    for row in range(rows):
        # A solution row is usable if any of its elements is unflagged.
        if shape:
            usable = not bool(flags[..., row].all())
        else:
            usable = not bool(flags[row])
        for scan, span in scan_spans.items():
            if span["start"] - 1.0 <= float(times[row]) <= span["end"] + 1.0:
                entry = coverage.setdefault(scan, {"rows": 0, "usable_rows": 0})
                entry["rows"] += 1
                entry["usable_rows"] += int(usable)
                break
    return coverage


def percent(value):
    return "   -  " if value is None else f"{100.0 * value:5.1f}%"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--outdir",
        default=None,
        help=(
            "Pipeline output directory. Defaults to outdir in config.yaml, "
            "looked up in the current directory and next to this script."
        ),
    )
    parser.add_argument(
        "--time",
        help="CASA log time to locate, e.g. 2020/11/22/14:05:38.4",
    )
    parser.add_argument(
        "--all-scans",
        action="store_true",
        help="List every scan, not only those that lost most of their data.",
    )
    args = parser.parse_args(argv)
    outdir, searched = resolve_outdir(args.outdir or configured_outdir())

    inspection = read_json(outdir / "inspection" / "inspection_report.json", {})
    if not inspection:
        raise SystemExit(
            "No pipeline output found (inspection/inspection_report.json). "
            "Looked in: " + ", ".join(str(item) for item in searched) +
            ". Run this where the pipeline wrote its products, or pass "
            "--outdir with that directory's path."
        )
    print(f"Output directory: {outdir}\n")

    names, roles = field_names(outdir, inspection)
    stages = flagging_stages(outdir) + postcal_stages(outdir)
    labels = [label for label, _ in stages]

    def describe_field(field_id):
        name = names.get(field_id, f"field {field_id}")
        role = ", ".join(roles.get(name, [])) or "unclassified"
        return f"{field_id} {name} ({role})"

    print("Flagged fraction by field at each step")
    print("=" * 78)
    field_ids = sorted(names)
    for field_id in field_ids:
        name = names[field_id]
        values = [
            (label, group_fractions(summary, "field").get(name))
            for label, summary in stages
        ]
        print(f"\n{describe_field(field_id)}")
        for label, value in values:
            print(f"    {label:26} {percent(value)}")

    spans = ((inspection.get("time_coverage") or {}).get("scan_spans")) or {}
    target_time = None
    if args.time:
        target_time = casa_time_to_mjd_seconds(args.time)

    print("\n\nScans")
    print("=" * 78)
    if not spans:
        print("No scan timing was recorded by inspection; rerun inspection.")
    for scan in sorted(spans, key=int):
        span = spans[scan]
        values = [
            (label, group_fractions(summary, "scan").get(str(scan)))
            for label, summary in stages
        ]
        final = next(
            (value for _, value in reversed(values) if value is not None),
            None,
        )
        contains = (
            target_time is not None and
            span["start"] - 30.0 <= target_time <= span["end"] + 30.0
        )
        lost = final is not None and final >= LOSS_THRESHOLD
        if not (args.all_scans or lost or contains):
            continue

        fields = ", ".join(describe_field(int(item)) for item in span["fields"])
        marker = "  <-- contains --time" if contains else ""
        print(
            f"\nscan {scan}: {mjd_seconds_to_utc(span['start'])} -> "
            f"{mjd_seconds_to_utc(span['end'])}{marker}"
        )
        print(f"    field: {fields}")
        for label, value in values:
            print(f"    {label:26} {percent(value)}")
        print(f"    most data removed: {attribute_loss(values)}")

    print("\n\nFully flagged channel runs after flagging (>= 99%)")
    print("=" * 78)
    masks = deliberate_channel_masks(outdir)
    per_channel = channel_fractions(outdir)
    if not per_channel:
        print("No per-channel statistics were saved (flagging/spw_statistics.json).")
    for spw in sorted(per_channel, key=lambda value: int(value)):
        runs = channel_runs(per_channel[spw])
        print(f"\nSPW {spw}: {len(per_channel[spw])} channels, {len(runs)} run(s)")
        for low, high in runs:
            f_low = channel_frequency_mhz(inspection, spw, low)
            f_high = channel_frequency_mhz(inspection, spw, high)
            frequency = (
                f"{min(f_low, f_high):8.3f}-{max(f_low, f_high):8.3f} MHz"
                if f_low is not None and f_high is not None
                else ""
            )
            kind = classify_run((low, high), masks.get(spw, set()))
            print(f"    chan {low:>5}-{high:<5} {frequency}  {kind}")

    try:
        from casatools import table as table_factory
    except Exception:
        table_factory = None

    if table_factory is not None and spans:
        print("\n\nPhase-solution coverage by scan (casatools)")
        print("=" * 78)
        caldir = outdir / "calibration"
        for label, path in (
            ("current", caldir / "phase.cal"),
            ("round 1", caldir / "round1" / "phase.cal"),
        ):
            if not path.is_dir():
                continue
            try:
                coverage = caltable_scan_coverage(path, spans, table_factory)
            except Exception as exc:
                print(f"{label}: could not read {path}: {exc}")
                continue
            empty = sorted(
                (scan for scan, item in coverage.items() if item["usable_rows"] == 0),
                key=int,
            )
            print(
                f"{label:8} {path}: solutions in {len(coverage)} scans; "
                f"{len(empty)} scan(s) with no usable solution: "
                f"{', '.join(empty) or 'none'}"
            )
    else:
        print(
            "\n(casatools unavailable: run inside CASA to also compare "
            "current and round-1 phase solutions per scan.)"
        )


if __name__ == "__main__":
    main()
