"""Memory evidence for pipeline stages, independent of CASA.

A stage ended by SIGKILL leaves nothing behind that says why. On a batch node
three different things end a CASA process that way, and each needs a
different remedy:

- the job exceeded its own memory limit. The kernel enforces that limit
  through the job's cgroup and counts an out-of-memory kill there;
- memory ran out outside the job. A job's limit caps what the job may use; it
  does not guarantee the node still has that much free, for example when
  other processes on the node use more than their own share;
- the job was cancelled, reached its time limit, or the node failed.

CASA cannot tell these apart either. casacore reads the job limit only from a
cgroup v1 hierarchy and otherwise reports the node's available memory, so
tclean's "Memory available" warning is often a node figure, not a job one.

This module records the job's cgroup accounting and the node's available
memory before a stage and while it runs, then states which cause the counters
support when the stage is killed. It reads only /proc and /sys/fs/cgroup.
"""

import os
import sys
import time

try:
    import resource
except ImportError:  # pragma: no cover - not available on every platform
    resource = None


GIB = 1024 ** 3

# cgroup v1 reports "no limit" as a huge page-aligned integer.
_UNLIMITED_BYTES = 1 << 60


def read_text(path):
    try:
        with open(path) as handle:
            return handle.read()
    except OSError:
        return None


def _integer(value):
    try:
        number = int(str(value).strip())
    except (TypeError, ValueError):
        return None
    return number if number >= 0 else None


def _limit(value):
    number = _integer(value)
    if number is None or number >= _UNLIMITED_BYTES:
        return None
    return number


def format_bytes(value):
    if value is None:
        return "unknown"
    return f"{value / GIB:.1f} GiB"


def parse_meminfo(text):
    """Return the node's total and available memory from /proc/meminfo."""

    values = {}
    for line in (text or "").splitlines():
        name, _, rest = line.partition(":")
        fields = rest.split()
        if not fields or not fields[0].isdigit():
            continue
        factor = 1024 if len(fields) > 1 and fields[1].lower() == "kb" else 1
        values[name.strip()] = int(fields[0]) * factor

    total = values.get("MemTotal")
    if total is None:
        return None

    available = values.get("MemAvailable")
    estimated = False
    if available is None and "MemFree" in values:
        # Kernels before 3.14 have no MemAvailable.
        available = (
            values["MemFree"] +
            values.get("Cached", 0) +
            values.get("Buffers", 0)
        )
        estimated = True

    return {
        "total_bytes": total,
        "available_bytes": available,
        "available_is_estimate": estimated,
    }


def parse_keyed_counts(text):
    """Parse cgroup files made of ``name value`` lines."""

    values = {}
    for line in (text or "").splitlines():
        fields = line.split()
        if len(fields) == 2:
            number = _integer(fields[1])
            if number is not None:
                values[fields[0]] = number
    return values


def _read_level(version, directory, read):
    def get(name):
        return read(f"{directory}/{name}")

    if version == 2:
        limit_text = get("memory.max")
        usage = _integer(get("memory.current"))
        events = parse_keyed_counts(get("memory.events"))
        if limit_text is None and usage is None and not events:
            return None
        stat = parse_keyed_counts(get("memory.stat"))
        return {
            "version": 2,
            "limit_bytes": _limit(limit_text),
            "usage_bytes": usage,
            "peak_bytes": _integer(get("memory.peak")),
            "oom_kill": events.get("oom_kill"),
            "limit_hits": events.get("max"),
            "anon_bytes": stat.get("anon"),
            "file_bytes": stat.get("file"),
            "shmem_bytes": stat.get("shmem"),
        }

    limit_text = get("memory.limit_in_bytes")
    usage = _integer(get("memory.usage_in_bytes"))
    if limit_text is None and usage is None:
        return None
    stat = parse_keyed_counts(get("memory.stat"))
    return {
        "version": 1,
        "limit_bytes": _limit(limit_text),
        "usage_bytes": usage,
        "peak_bytes": _integer(get("memory.max_usage_in_bytes")),
        # Present on kernels 4.13 and later.
        "oom_kill": parse_keyed_counts(get("memory.oom_control")).get(
            "oom_kill"
        ),
        "limit_hits": _integer(get("memory.failcnt")),
        "anon_bytes": stat.get("total_rss", stat.get("rss")),
        "file_bytes": stat.get("total_cache", stat.get("cache")),
        "shmem_bytes": stat.get("total_shmem", stat.get("shmem")),
    }


def cgroup_levels(read=read_text):
    """Memory accounting for every cgroup from this process up to the root.

    Returned deepest first. A Slurm job's limit is usually set on an ancestor
    of the cgroup the process runs in, so every level is read.
    """

    levels = []
    membership = read("/proc/self/cgroup") or ""

    for line in membership.splitlines():
        parts = line.strip().split(":", 2)
        if len(parts) != 3:
            continue
        _hierarchy, controllers, relative = parts
        if controllers == "":
            version, root = 2, "/sys/fs/cgroup"
        elif "memory" in controllers.split(","):
            version, root = 1, "/sys/fs/cgroup/memory"
        else:
            continue

        segments = [item for item in relative.split("/") if item]
        for depth in range(len(segments), -1, -1):
            path = "/" + "/".join(segments[:depth])
            directory = root if path == "/" else root + path
            level = _read_level(version, directory, read)
            if level is not None:
                level["path"] = path
                levels.append(level)

    return levels


def job_limit_level(levels):
    """The level whose limit binds this process (smallest; deepest on ties)."""

    limited = [level for level in levels if level.get("limit_bytes")]
    if not limited:
        return None
    return min(limited, key=lambda level: level["limit_bytes"])


def _within(level, container):
    if level.get("version") != container.get("version"):
        return False
    base = container["path"].rstrip("/")
    return level["path"] == container["path"] or level["path"].startswith(
        base + "/"
    )


def job_levels(snapshot):
    """Levels from this process's cgroup up to the one carrying the limit.

    Levels above the job limit are excluded: under cgroup v2 their event
    counters also include other jobs on the node.
    """

    container = snapshot.get("job_limit")
    if not container:
        return []
    return [
        level
        for level in snapshot.get("cgroup_levels", [])
        if _within(level, container)
    ]


def children_max_rss_bytes():
    """Peak resident memory of the largest child process that has exited."""

    if resource is None:
        return None
    try:
        value = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    except (OSError, ValueError):
        return None
    # Linux reports KiB; macOS reports bytes.
    return int(value) if sys.platform == "darwin" else int(value) * 1024


def snapshot(read=read_text, environ=None, clock=time.time, child_rss=None):
    environ = os.environ if environ is None else environ
    levels = cgroup_levels(read)
    limit_level = job_limit_level(levels)
    slurm_mb = _integer(environ.get("SLURM_MEM_PER_NODE"))
    per_cpu_mb = _integer(environ.get("SLURM_MEM_PER_CPU"))
    cpus_on_node = _integer(environ.get("SLURM_CPUS_ON_NODE"))

    return {
        "time_unix": clock(),
        "node": parse_meminfo(read("/proc/meminfo")),
        "cgroup_levels": levels,
        "job_limit": dict(limit_level) if limit_level else None,
        "slurm": {
            "job_id": environ.get("SLURM_JOB_ID"),
            "node_name": environ.get("SLURMD_NODENAME"),
            "mem_per_node_bytes": (
                slurm_mb * 1024 ** 2 if slurm_mb else None
            ),
            "mem_per_cpu_bytes": (
                per_cpu_mb * 1024 ** 2 if per_cpu_mb else None
            ),
            "cpus_on_node": cpus_on_node,
        },
        "children_max_rss_bytes": (
            children_max_rss_bytes() if child_rss is None else child_rss
        ),
    }


def _node_available(value):
    node = (value or {}).get("node") or {}
    return node.get("available_bytes")


def summary_line(value):
    limit = (value or {}).get("job_limit")
    node = (value or {}).get("node") or {}
    parts = []
    if limit:
        parts.append(
            f"job limit {format_bytes(limit['limit_bytes'])} "
            f"(cgroup v{limit['version']}), job usage "
            f"{format_bytes(limit.get('usage_bytes'))}"
        )
    else:
        parts.append("no job memory limit visible in the cgroup hierarchy")
    if node:
        parts.append(
            f"node available {format_bytes(node.get('available_bytes'))} "
            f"of {format_bytes(node.get('total_bytes'))}"
        )
    return "Memory: " + "; ".join(parts) + "."


def low_node_memory_warning(value, fraction=0.5):
    """Warn when the node cannot supply most of the job's remaining allowance."""

    limit = (value or {}).get("job_limit")
    available = _node_available(value)
    if not limit or available is None:
        return None

    headroom = limit["limit_bytes"] - (limit.get("usage_bytes") or 0)
    if headroom <= 0 or available >= fraction * headroom:
        return None

    return (
        f"WARNING: the node has {format_bytes(available)} available, less "
        f"than {fraction:.0%} of this job's remaining "
        f"{format_bytes(headroom)} memory allowance. A job's memory limit "
        "caps what it may use; it does not keep that memory free on the "
        "node, so this stage can be killed well below the job limit."
    )


class MemoryWatch:
    """Track the extremes of memory use while a stage runs."""

    def __init__(self, start):
        self.samples = 0
        self.lowest_node_available_bytes = _node_available(start)
        self.lowest_node_available_elapsed_s = 0.0
        limit = (start or {}).get("job_limit") or {}
        self.highest_job_usage_bytes = limit.get("usage_bytes")
        self.highest_job_usage_elapsed_s = 0.0

    def observe(self, value, elapsed_s):
        self.samples += 1
        available = _node_available(value)
        if available is not None and (
            self.lowest_node_available_bytes is None or
            available < self.lowest_node_available_bytes
        ):
            self.lowest_node_available_bytes = available
            self.lowest_node_available_elapsed_s = float(elapsed_s)

        usage = ((value or {}).get("job_limit") or {}).get("usage_bytes")
        if usage is not None and (
            self.highest_job_usage_bytes is None or
            usage > self.highest_job_usage_bytes
        ):
            self.highest_job_usage_bytes = usage
            self.highest_job_usage_elapsed_s = float(elapsed_s)

    def as_dict(self):
        return {
            "samples": self.samples,
            "lowest_node_available_bytes": self.lowest_node_available_bytes,
            "lowest_node_available_elapsed_s": (
                self.lowest_node_available_elapsed_s
            ),
            "highest_job_usage_bytes": self.highest_job_usage_bytes,
            "highest_job_usage_elapsed_s": self.highest_job_usage_elapsed_s,
        }


def _counter_delta(before, after, key):
    """Largest increase of a counter across the job's cgroup levels.

    The largest, not the sum: cgroup v2 counts one kill at every level.
    """

    earlier = {
        (level["version"], level["path"]): level.get(key)
        for level in job_levels(before)
    }
    deltas = []
    readable = False
    for level in job_levels(after):
        current = level.get(key)
        if current is None:
            continue
        readable = True
        previous = earlier.get((level["version"], level["path"]))
        deltas.append(current - previous if previous is not None else current)
    return (max(deltas) if deltas else None), readable


def summary_at_exit(before, after, watch=None):
    """One line of memory context for a stage that died on any signal."""

    limit = (after or {}).get("job_limit") or (before or {}).get("job_limit")
    watch_values = watch.as_dict() if watch is not None else {}
    parts = []
    if limit:
        parts.append(
            f"job limit {format_bytes(limit.get('limit_bytes'))}, "
            f"peak {format_bytes(limit.get('peak_bytes'))}, "
            f"in use {format_bytes(limit.get('usage_bytes'))}"
        )
    highest = watch_values.get("highest_job_usage_bytes")
    if highest is not None:
        parts.append(f"highest sampled {format_bytes(highest)}")
    available = _node_available(after) or _node_available(before)
    if available is not None:
        parts.append(f"node available {format_bytes(available)}")
    kills, _ = _counter_delta(before, after, "oom_kill")
    if kills is not None:
        parts.append(
            f"{kills} out-of-memory kill(s) counted against the job"
            if kills else "no out-of-memory kill counted against the job"
        )
    rss_after = (after or {}).get("children_max_rss_bytes")
    rss_before = (before or {}).get("children_max_rss_bytes")
    if rss_after is not None and (rss_before is None or rss_after > rss_before):
        parts.append(f"largest process {format_bytes(rss_after)} resident")
    return "; ".join(parts)


def attribute_kill(before, after, watch=None):
    """State what the memory counters say about a SIGKILLed stage."""

    watch_values = watch.as_dict() if watch is not None else {}
    limit = (after or {}).get("job_limit") or (before or {}).get("job_limit")
    kills, kills_readable = _counter_delta(before, after, "oom_kill")
    hits, _ = _counter_delta(before, after, "limit_hits")

    rss_before = (before or {}).get("children_max_rss_bytes")
    rss_after = (after or {}).get("children_max_rss_bytes")
    stage_peak_rss = (
        rss_after
        if rss_after is not None and (rss_before is None or rss_after > rss_before)
        else None
    )

    node = (before or {}).get("node") or {}
    evidence = {
        "job_limit_bytes": (limit or {}).get("limit_bytes"),
        "job_limit_cgroup": (limit or {}).get("path"),
        "cgroup_version": (limit or {}).get("version"),
        "job_usage_before_bytes": (
            ((before or {}).get("job_limit") or {}).get("usage_bytes")
        ),
        "job_usage_after_bytes": (
            ((after or {}).get("job_limit") or {}).get("usage_bytes")
        ),
        "job_peak_bytes": (
            ((after or {}).get("job_limit") or {}).get("peak_bytes")
        ),
        "oom_kills_during_stage": kills,
        "limit_hits_during_stage": hits,
        "node_total_bytes": node.get("total_bytes"),
        "node_available_before_bytes": _node_available(before),
        "node_available_after_bytes": _node_available(after),
        "stage_peak_rss_bytes": stage_peak_rss,
        **watch_values,
    }

    lowest = watch_values.get("lowest_node_available_bytes")
    node_sentence = ""
    if evidence["node_available_before_bytes"] is not None:
        node_sentence = (
            " The node had "
            f"{format_bytes(evidence['node_available_before_bytes'])} "
            "available when the stage started"
        )
        if lowest is not None and lowest < evidence[
            "node_available_before_bytes"
        ]:
            node_sentence += (
                f" and {format_bytes(lowest)} at the lowest sample"
            )
        node_sentence += "."

    rss_sentence = (
        " The killed stage's largest process reached "
        f"{format_bytes(stage_peak_rss)} resident memory."
        if stage_peak_rss is not None else ""
    )

    if kills is not None and kills > 0:
        cause = "job_memory_limit"
        explanation = (
            f"The job's memory cgroup recorded {kills} out-of-memory kill(s) "
            "during this stage, so the stage exceeded the job's "
            f"{format_bytes(evidence['job_limit_bytes'])} memory limit "
            f"(cgroup {evidence['job_limit_cgroup']})"
        )
        if evidence["job_peak_bytes"] is not None:
            explanation += (
                f" (job peak {format_bytes(evidence['job_peak_bytes'])})"
            )
        explanation += "." + rss_sentence
    elif limit and kills_readable:
        cause = "outside_job_limit"
        explanation = (
            "No out-of-memory kill was recorded against the job's "
            f"{format_bytes(evidence['job_limit_bytes'])} memory limit, so "
            "the stage did not exceed it. The SIGKILL came from outside the "
            "job: memory exhaustion elsewhere on the node, cancellation or a "
            "time limit, or a node failure." + node_sentence + rss_sentence
        )
    elif limit:
        cause = "undetermined"
        explanation = (
            "The job has a "
            f"{format_bytes(evidence['job_limit_bytes'])} cgroup memory limit, "
            "but this kernel does not report out-of-memory kill counts, so "
            "whether the stage exceeded it cannot be established here."
        )
        if hits:
            explanation += (
                f" The job reached its limit {hits} time(s) during the stage."
            )
        explanation += node_sentence + rss_sentence
    else:
        cause = "no_job_limit_visible"
        explanation = (
            "No cgroup memory limit is visible to the pipeline, so this job's "
            "memory is not enforced separately from the rest of the node and "
            "a kill here is most likely node-wide memory exhaustion, "
            "cancellation, or a node failure." + node_sentence + rss_sentence
        )

    return {"cause": cause, "explanation": explanation, "evidence": evidence}


def effective_limit_bytes(value):
    """The most memory this process may use: its cgroup limit, else the node."""

    limit = (value or {}).get("job_limit")
    if limit and limit.get("limit_bytes"):
        return limit["limit_bytes"]
    return ((value or {}).get("node") or {}).get("total_bytes")


def requested_slurm_bytes(value):
    slurm = (value or {}).get("slurm") or {}
    if slurm.get("mem_per_node_bytes"):
        return slurm["mem_per_node_bytes"]
    if slurm.get("mem_per_cpu_bytes") and slurm.get("cpus_on_node"):
        return slurm["mem_per_cpu_bytes"] * slurm["cpus_on_node"]
    return None


def allocation_warning(value, fraction=0.5):
    """Explain a memory limit far below what the run was meant to have."""

    limit = (value or {}).get("job_limit")
    if not limit:
        return None
    slurm = (value or {}).get("slurm") or {}
    amount = format_bytes(limit["limit_bytes"])

    if not slurm.get("job_id"):
        return (
            "This process is not running inside a Slurm job, and its cgroup "
            f"({limit['path']}) limits it to {amount}. Login nodes and "
            "interactive sessions commonly impose limits like this; submit "
            "the pipeline with sbatch."
        )

    requested = requested_slurm_bytes(value)
    if requested and limit["limit_bytes"] < fraction * requested:
        return (
            f"This process may use only {amount} (cgroup {limit['path']}), "
            f"although Slurm job {slurm['job_id']} requested "
            f"{format_bytes(requested)} on this node. It is running in a job "
            "step or session with a smaller share of the memory, for example "
            "under srun or in an interactive shell. Run pipeline.py directly "
            "from the batch script."
        )
    return None


# ------------------------------------------------------------------------------
# CASA w-projection imaging memory
# ------------------------------------------------------------------------------

# From CASA's WPConvFunc::findConvFunction: the memory figure it plans from is
# never taken below 4000 MB, and at most 100 w-planes are computed per chunk.
CASA_WPROJECT_MINIMUM_PLANNING_MB = 4000
CASA_WPROJECT_MAX_PLANES_PER_CHUNK = 100


def composite_even_at_least(value):
    """Smallest even 2^a 3^b 5^c not below ``value``, as CASA pads grids."""

    value = max(2, int(value))
    best = None
    two = 2
    while two < 4 * value:
        three = two
        while three < 4 * value:
            number = three
            while number < 4 * value:
                if number >= value and (best is None or number < best):
                    best = number
                number *= 5
            three *= 3
        two *= 2
    return best


def _convolution_peak(size, planes, threads, per_chunk):
    """Peak bytes of one convolution-function build at a given chunk size.

    Shared: an FFTPACK save array and the first screen (8 bytes x size^2
    each). Per chunk: a complex cube of (size/2 - 1)^2 per plane. Per thread:
    a complex screen and an FFTPACK work array (16 bytes x size^2). Chunks run
    last first, and each sets the OpenMP thread count to at most its own
    plane count for the rest of the process.
    """

    plane = size * size
    quarter = (size // 2 - 1) ** 2
    chunks = planes // per_chunk + (1 if planes % per_chunk else 0)
    sizes = [per_chunk] * chunks
    if planes % per_chunk:
        sizes[-1] = planes % per_chunk

    active = max(1, int(threads))
    peak = 0
    for index in range(chunks - 1, -1, -1):
        active = min(active, sizes[index])
        peak = max(
            peak,
            16 * plane + sizes[index] * 8 * quarter + active * 16 * plane,
        )
    return peak


def _planes_per_chunk(size, budget_bytes):
    """How many w-planes CASA computes at once for a given memory figure."""

    megabytes = max(
        CASA_WPROJECT_MINIMUM_PLANNING_MB,
        int(budget_bytes) // (1024 * 1024),
    )
    per_chunk = int((megabytes / 8.0 * 1024 * 1024) / (size * size))
    return max(1, min(CASA_WPROJECT_MAX_PLANES_PER_CHUNK, per_chunk))


def wproject_imaging_memory(imsize, nterms, wprojplanes, threads, padding=1.2,
                            budget_bytes=None):
    """Memory range for one serial tclean with gridder='wproject'.

    CASA chooses how many w-planes to compute at once from its own view of
    free memory, which is often the node's rather than the job's (see the
    module docstring). ``minimum_bytes`` is the least any such choice can use,
    so an allocation below it cannot succeed. ``largest_bytes`` is the most
    the convolution functions can take under any choice, plus the same
    image products and grid. Neither counts CASA's other working memory.

    ``budget_bytes`` is the figure CASA will actually plan from, when it is
    known: the run's casarc sets one, so what CASA will do is a number rather
    than a range, and ``planned_bytes`` reports it.
    """

    if isinstance(imsize, (list, tuple)):
        imsize = max(int(item) for item in imsize)
    imsize = int(imsize)
    nterms = max(1, int(nterms))
    wprojplanes = int(wprojplanes)

    size = composite_even_at_least(int(imsize * padding))
    plane = size * size

    if wprojplanes > 1:
        smallest_chunk = int(
            (CASA_WPROJECT_MINIMUM_PLANNING_MB / 8.0 * 1024 * 1024) / plane
        )
        smallest_chunk = max(
            1,
            min(CASA_WPROJECT_MAX_PLANES_PER_CHUNK, smallest_chunk),
        )
        peaks = [
            _convolution_peak(size, wprojplanes, threads, chunk)
            for chunk in range(
                smallest_chunk,
                CASA_WPROJECT_MAX_PLANES_PER_CHUNK + 1,
            )
        ]
        convolution_minimum = min(peaks)
        convolution_largest = max(peaks)
    else:
        convolution_minimum = convolution_largest = 0

    # Residual and model per Taylor term, and 2*nterms-1 PSF terms, as 4-byte
    # images; one padded complex grid for the major cycle.
    images = (4 * nterms - 1) * 4 * imsize * imsize
    grid = 8 * plane

    requirement = {
        "imsize": imsize,
        "nterms": nterms,
        "wprojplanes": wprojplanes,
        "threads": max(1, int(threads)),
        "convolution_size": size,
        "convolution_minimum_bytes": convolution_minimum,
        "convolution_largest_bytes": convolution_largest,
        "image_and_grid_bytes": images + grid,
        "minimum_bytes": convolution_minimum + images + grid,
        "largest_bytes": convolution_largest + images + grid,
    }

    if budget_bytes:
        if wprojplanes > 1:
            per_chunk = _planes_per_chunk(size, budget_bytes)
            convolution = _convolution_peak(
                size, wprojplanes, threads, per_chunk
            )
        else:
            per_chunk = 0
            convolution = 0
        requirement["planning_budget_bytes"] = int(budget_bytes)
        requirement["planned_planes_per_chunk"] = per_chunk
        requirement["planned_bytes"] = convolution + images + grid

    return requirement
