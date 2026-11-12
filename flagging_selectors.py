#!/usr/bin/env python3
"""CASA-free construction of channel-accurate frequency flag selectors."""

from __future__ import annotations

import math
import os
import re
from typing import Mapping


_FREQUENCY_VALUE = re.compile(
    r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"\s*(GHz|MHz|kHz|Hz)?\s*$",
    re.IGNORECASE,
)


def _spw_sort_key(value):
    text = str(value)
    try:
        return (0, int(text))
    except ValueError:
        return (1, text)


def _channel_frequencies_mhz(spw, info):
    if not isinstance(info, Mapping):
        raise ValueError(f"SPW {spw} metadata must be a mapping.")

    values = info.get("channel_frequencies_MHz")
    if not isinstance(values, (list, tuple)) or not values:
        raise ValueError(
            f"SPW {spw} has no channel_frequencies_MHz metadata. Rerun "
            "the inspection stage before applying frequency hard masks."
        )

    try:
        frequencies = [float(value) for value in values]
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"SPW {spw} contains a non-numeric channel frequency."
        ) from exc

    if not all(math.isfinite(value) for value in frequencies):
        raise ValueError(f"SPW {spw} contains a non-finite channel frequency.")

    try:
        expected = int(info.get("nchan", len(frequencies)))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"SPW {spw} has an invalid nchan value.") from exc
    if expected != len(frequencies):
        raise ValueError(
            f"SPW {spw} reports nchan={expected}, but inspection recorded "
            f"{len(frequencies)} channel frequencies."
        )

    return frequencies


def _compress_indices(indices):
    ordered = sorted(set(int(index) for index in indices))
    if not ordered:
        return []

    ranges = []
    start = previous = ordered[0]
    for index in ordered[1:]:
        if index == previous + 1:
            previous = index
            continue
        ranges.append((start, previous))
        start = previous = index
    ranges.append((start, previous))
    return ranges


def _format_ranges(ranges):
    return ";".join(
        str(start) if start == end else f"{start}~{end}"
        for start, end in ranges
    )


def _spw_selector(spw, indices):
    ranges = _compress_indices(indices)
    if not ranges:
        return None, []
    formatted = _format_ranges(ranges)
    return f"{spw}:{formatted}", [list(item) for item in ranges]


def build_science_range_channel_masks(spws, low_mhz, high_mhz):
    """Select only channel centres outside a guaranteed science range."""

    try:
        low = float(low_mhz)
        high = float(high_mhz)
    except (TypeError, ValueError) as exc:
        raise ValueError("Science-range limits must be numeric MHz values.") from exc
    if not (math.isfinite(low) and math.isfinite(high) and low < high):
        raise ValueError("Science-range limits must be finite and increasing.")
    if not isinstance(spws, Mapping) or not spws:
        raise ValueError("Inspection report contains no spectral-window metadata.")

    selections = []
    reports = []
    for spw in sorted(spws, key=_spw_sort_key):
        frequencies = _channel_frequencies_mhz(spw, spws[spw])
        outside = [
            index
            for index, frequency in enumerate(frequencies)
            if frequency < low or frequency > high
        ]
        selector, channel_ranges = _spw_selector(spw, outside)
        if selector:
            selections.append(selector)

        reports.append({
            "spw": str(spw),
            "coverage_MHz": [min(frequencies), max(frequencies)],
            "total_channels": len(frequencies),
            "masked_channels": len(outside),
            "channel_ranges": channel_ranges,
            "selection": selector,
            "status": (
                "outside_science_range_masked"
                if selector else
                "within_science_range"
            ),
        })

    return {
        "scheme": "outside_guaranteed_science_range",
        "science_range_MHz": [low, high],
        "spw_selection": ",".join(selections),
        "spws": reports,
    }


def _parse_frequency_value(value):
    match = _FREQUENCY_VALUE.fullmatch(str(value))
    if not match:
        raise ValueError(f"Invalid frequency value '{value}'.")
    number = float(match.group(1))
    unit = (match.group(2) or "").lower()
    return number, unit


def _to_mhz(number, unit):
    factors = {
        "ghz": 1e3,
        "mhz": 1.0,
        "khz": 1e-3,
        "hz": 1e-6,
    }
    return number * factors[unit]


def parse_frequency_range_mhz(value):
    """Parse a CASA-style two-endpoint frequency range into MHz."""

    text = str(value).strip()
    if text.startswith("*:"):
        text = text[2:]
    parts = text.split("~")
    if len(parts) != 2:
        raise ValueError(
            f"Frequency mask '{value}' must contain exactly one '~'."
        )

    left_number, left_unit = _parse_frequency_value(parts[0])
    right_number, right_unit = _parse_frequency_value(parts[1])
    if not left_unit and not right_unit:
        raise ValueError(
            f"Frequency mask '{value}' must include Hz, kHz, MHz, or GHz."
        )
    left_unit = left_unit or right_unit
    right_unit = right_unit or left_unit
    low = _to_mhz(left_number, left_unit)
    high = _to_mhz(right_number, right_unit)
    if not (math.isfinite(low) and math.isfinite(high)):
        raise ValueError(f"Frequency mask '{value}' is not finite.")
    if low > high:
        low, high = high, low
    return low, high


def build_intersected_static_masks(spws, frequency_masks):
    """Intersect static masks with channel centres in each actual SPW."""

    if not isinstance(spws, Mapping) or not spws:
        raise ValueError("Inspection report contains no spectral-window metadata.")

    parsed_masks = [
        (str(mask), *parse_frequency_range_mhz(mask))
        for mask in (frequency_masks or [])
    ]
    selected_by_spw = {str(spw): set() for spw in spws}
    applied = []
    skipped = []

    for spw in sorted(spws, key=_spw_sort_key):
        spw_text = str(spw)
        frequencies = _channel_frequencies_mhz(spw, spws[spw])
        coverage_low = min(frequencies)
        coverage_high = max(frequencies)

        for configured, mask_low, mask_high in parsed_masks:
            overlap_low = max(coverage_low, mask_low)
            overlap_high = min(coverage_high, mask_high)
            base = {
                "spw": spw_text,
                "configured_range": configured,
                "configured_range_MHz": [mask_low, mask_high],
                "spw_coverage_MHz": [coverage_low, coverage_high],
            }
            if overlap_low > overlap_high:
                skipped.append({
                    **base,
                    "reason": "no_spw_overlap",
                })
                continue

            indices = [
                index
                for index, frequency in enumerate(frequencies)
                if mask_low <= frequency <= mask_high
            ]
            selector, channel_ranges = _spw_selector(spw_text, indices)
            if not selector:
                skipped.append({
                    **base,
                    "intersection_MHz": [overlap_low, overlap_high],
                    "reason": "no_channel_centres_in_overlap",
                })
                continue

            selected_by_spw[spw_text].update(indices)
            applied.append({
                **base,
                "intersection_MHz": [overlap_low, overlap_high],
                "channel_count": len(indices),
                "channel_ranges": channel_ranges,
                "selection": selector,
            })

    selections = []
    combined = []
    for spw in sorted(selected_by_spw, key=_spw_sort_key):
        selector, channel_ranges = _spw_selector(
            spw,
            selected_by_spw[spw],
        )
        if not selector:
            continue
        selections.append(selector)
        combined.append({
            "spw": spw,
            "channel_count": len(selected_by_spw[spw]),
            "channel_ranges": channel_ranges,
            "selection": selector,
        })

    return {
        "scheme": "spw_intersected_static_frequency_masks",
        "configured_frequency_ranges": [item[0] for item in parsed_masks],
        "spw_selection": ",".join(selections),
        "applied_intersections": applied,
        "skipped_non_overlaps": skipped,
        "combined_spw_selections": combined,
        "status": "applied" if selections else "no_overlapping_channels",
    }


# ---------------------------------------------------------------------------
# flagdata time-chunk sizing
# ---------------------------------------------------------------------------
#
# CASA flagdata holds one time chunk per field and spectral window in memory.
# With ntime='scan' that chunk is an entire scan: every baseline, channel and
# correlation at once. Its size scales with scan length, channel count and
# array size together, so on high-resolution uGMRT GWB or MeerKAT 32k data one
# scan can exceed the node's memory limit. The kernel or Slurm then kills CASA
# with SIGKILL and nothing is written.

# Visibility cube (complex float, 8 bytes) + flag cube (1 byte) + the handler's
# copy of the flags used to detect modification (1 byte), rounded up for the
# auto-flagger's working arrays. This is an engineering estimate; calibrate it
# against `sacct --format=MaxRSS` for the data in hand.
FLAGDATA_BYTES_PER_ELEMENT = 12.0
# Headroom for CASA itself, Python, flag-version writes and estimate error.
FLAGDATA_MEMORY_FRACTION = 0.5
# TFCrop fits a trend in time inside each chunk; fewer integrations than this
# do not constrain that fit meaningfully.
FLAGDATA_MINIMUM_INTEGRATIONS = 4

# cgroup v1 reports "no limit" as a huge page-aligned integer.
_UNLIMITED_CGROUP_BYTES = 1 << 60


def _finite_positive(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and number > 0 else None


def _read_text(path):
    try:
        with open(path) as handle:
            return handle.read()
    except OSError:
        return None


def _cgroup_limits(read_text):
    """Collect every numeric cgroup memory limit on this process's path.

    A Slurm job's limit is normally set on an ancestor of the task cgroup, so
    every level from the process's own cgroup up to the root is checked.
    """

    limits = []
    membership = read_text("/proc/self/cgroup") or ""

    for line in membership.splitlines():
        parts = line.strip().split(":", 2)
        if len(parts) != 3:
            continue
        _hierarchy, controllers, relative = parts

        if controllers == "":
            root, filename = "/sys/fs/cgroup", "memory.max"
        elif "memory" in controllers.split(","):
            root, filename = "/sys/fs/cgroup/memory", "memory.limit_in_bytes"
        else:
            continue

        segments = [item for item in relative.split("/") if item]
        for depth in range(len(segments), -1, -1):
            directory = "/".join([root] + segments[:depth])
            text = read_text(f"{directory}/{filename}")
            if text is None:
                continue
            value = text.strip()
            if value in {"", "max"}:
                continue
            number = _finite_positive(value)
            if number is not None and number < _UNLIMITED_CGROUP_BYTES:
                limits.append(int(number))

    return min(limits) if limits else None


def detect_memory_limit_bytes(
    environ=None,
    read_text=None,
    physical_bytes=None,
    configured_gb=None,
):
    """Return the tightest memory limit that applies to this process.

    Every available source is consulted and the smallest wins: an explicit
    configuration value, the Slurm allocation, the cgroup limit the kernel
    actually enforces, and physical memory. Any one of them can kill the job.
    """

    environ = os.environ if environ is None else environ
    read_text = _read_text if read_text is None else read_text
    sources = {}

    configured = _finite_positive(configured_gb)
    if configured is not None:
        sources["configured_memory_limit_gb"] = int(configured * 1024 ** 3)

    per_node_mb = _finite_positive(environ.get("SLURM_MEM_PER_NODE"))
    if per_node_mb is not None:
        sources["slurm_mem_per_node"] = int(per_node_mb * 1024 ** 2)

    per_cpu_mb = _finite_positive(environ.get("SLURM_MEM_PER_CPU"))
    cpus = _finite_positive(environ.get("SLURM_CPUS_ON_NODE"))
    if per_cpu_mb is not None and cpus is not None:
        sources["slurm_mem_per_cpu"] = int(per_cpu_mb * cpus * 1024 ** 2)

    cgroup = _cgroup_limits(read_text)
    if cgroup is not None:
        sources["cgroup"] = cgroup

    if physical_bytes is None:
        try:
            physical_bytes = (
                os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
            )
        except (AttributeError, OSError, ValueError):
            physical_bytes = None
    physical = _finite_positive(physical_bytes)
    if physical is not None:
        sources["physical_memory"] = int(physical)

    return {
        "limit_bytes": min(sources.values()) if sources else None,
        "limiting_source": (
            min(sources, key=sources.get) if sources else None
        ),
        "sources_bytes": sources,
    }


def summarize_scan_chunks(ms_path, msmd_factory, table_factory):
    """Describe the time chunks flagdata will load, one per scan and SPW.

    CASA tools are injected so this module stays importable without CASA.
    Rows per integration assume every cross-correlation plus autocorrelation
    is present; that is an upper bound, which is the safe direction here.
    """

    ncorr = 1
    table_tool = table_factory()
    opened = False
    try:
        table_tool.open(f"{ms_path}/POLARIZATION")
        opened = True
        values = table_tool.getcol("NUM_CORR")
        flat = getattr(values, "flat", values)
        ncorr = max(int(value) for value in flat)
    finally:
        if opened:
            table_tool.close()

    msmd = msmd_factory()
    msmd.open(str(ms_path))
    groups = []
    try:
        for scan in msmd.scannumbers():
            scan = int(scan)
            times = sorted({float(value) for value in msmd.timesforscan(scan)})
            if not times:
                continue

            steps = sorted(
                later - earlier
                for earlier, later in zip(times, times[1:])
                if later > earlier
            )
            antennas = len(list(msmd.antennasforscan(scan)))
            rows_per_integration = antennas * (antennas + 1) // 2

            for spw in msmd.spwsforscan(scan):
                spw = int(spw)
                integration = steps[len(steps) // 2] if steps else None
                if integration is None:
                    try:
                        integration = float(
                            msmd.exposuretime(scan=scan, spwid=spw)["value"]
                        )
                    except Exception:
                        integration = 1.0

                groups.append({
                    "scan": scan,
                    "spw": spw,
                    "rows_per_integration": rows_per_integration,
                    "integration_s": float(integration),
                    "n_integrations": len(times),
                    "nchan": int(msmd.nchan(spw)),
                    "ncorr": ncorr,
                })
    finally:
        msmd.close()

    return groups


def plan_flagdata_time_chunk(
    groups,
    limit_bytes,
    *,
    memory_fraction=FLAGDATA_MEMORY_FRACTION,
    bytes_per_element=FLAGDATA_BYTES_PER_ELEMENT,
    minimum_integrations=FLAGDATA_MINIMUM_INTEGRATIONS,
):
    """Choose ntime so no flagdata chunk exceeds the memory budget.

    Whole-scan chunks are kept wherever they fit, so validated behaviour is
    unchanged for data that were already safe. Only when the largest scan
    would not fit is a time interval chosen, and only as long as the most
    demanding scan/SPW allows. A chunk never crosses a scan boundary either
    way (combinescans=False).

    ``ntime`` is not purely an execution setting: TFCrop fits a time trend
    within each chunk and flag extension's ``growtime`` is evaluated over the
    chunk. A bounded value is therefore always reported, never silent.
    """

    groups = [dict(group) for group in groups or ()]
    fraction = _finite_positive(memory_fraction)
    element_bytes = _finite_positive(bytes_per_element)
    if fraction is None or fraction > 1.0:
        raise ValueError("memory_fraction must lie in (0, 1].")
    if element_bytes is None:
        raise ValueError("bytes_per_element must be positive.")

    result = {
        "ntime": "scan",
        "memory_limit_bytes": limit_bytes,
        "memory_fraction": fraction,
        "bytes_per_element": element_bytes,
        "minimum_integrations": int(minimum_integrations),
    }

    if not groups:
        result.update(status="not_measured", reason="no scan geometry")
        return result

    for group in groups:
        group["integration_bytes"] = (
            group["rows_per_integration"] *
            group["nchan"] *
            group["ncorr"] *
            element_bytes
        )
        group["scan_chunk_bytes"] = (
            group["integration_bytes"] * max(1, group["n_integrations"])
        )

    largest = max(groups, key=lambda item: item["scan_chunk_bytes"])
    result["largest_scan_chunk"] = {
        key: largest[key]
        for key in (
            "scan", "spw", "rows_per_integration", "integration_s",
            "n_integrations", "nchan", "ncorr", "scan_chunk_bytes",
        )
    }

    budget = _finite_positive(limit_bytes)
    if budget is None:
        result.update(
            status="no_memory_limit_detected",
            reason="whole-scan chunks retained; no limit was available",
        )
        return result

    budget *= fraction
    result["budget_bytes"] = budget

    if largest["scan_chunk_bytes"] <= budget:
        result.update(
            status="scan_fits",
            reason="the largest whole-scan chunk fits the memory budget",
        )
        return result

    constrained = [
        group for group in groups if group["scan_chunk_bytes"] > budget
    ]
    chunk_s = min(
        budget * group["integration_s"] / group["integration_bytes"]
        for group in constrained
    )
    chunk_s = math.floor(chunk_s)

    shortfall = [
        group
        for group in constrained
        if chunk_s < minimum_integrations * group["integration_s"]
    ]
    if chunk_s < 1 or shortfall:
        worst = max(constrained, key=lambda item: item["integration_bytes"])
        required = (
            worst["integration_bytes"] * minimum_integrations / fraction
        )
        result.update(
            status="insufficient_memory",
            required_memory_bytes=required,
            reason=(
                "even a chunk of "
                f"{minimum_integrations} integrations exceeds the budget"
            ),
        )
        return result

    result.update(
        ntime=f"{int(chunk_s)}s",
        status="bounded",
        chunk_s=int(chunk_s),
        reason=(
            "the largest whole-scan chunk exceeds the memory budget, so "
            "flagdata time chunks are bounded"
        ),
    )
    return result

