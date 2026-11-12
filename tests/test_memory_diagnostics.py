import contextlib
import importlib.util
import io
import json
import os
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, SOURCE_DIR / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MEMORY = load("test_radio_continuum_memory_diagnostics", "memory_diagnostics.py")
PIPELINE = load("test_radio_continuum_pipeline_memory", "pipeline.py")

GIB = 1024 ** 3
JOB = "/system.slice/slurmstepd.scope/job_42"
TASK = JOB + "/step_batch/user/task_0"


def meminfo(total_gib, available_gib):
    return (
        f"MemTotal:       {total_gib * 1024 * 1024} kB\n"
        f"MemFree:        {available_gib * 1024 * 1024 // 2} kB\n"
        f"MemAvailable:   {available_gib * 1024 * 1024} kB\n"
    )


def events(oom_kill, limit_hits=0):
    return f"low 0\nhigh 0\nmax {limit_hits}\noom 0\noom_kill {oom_kill}\n"


def cgroup_v2(job_kills=0, task_kills=0, neighbour_kills=0, available_gib=3,
              job_usage_gib=2):
    root = "/sys/fs/cgroup"
    files = {
        "/proc/self/cgroup": f"0::{TASK}\n",
        "/proc/meminfo": meminfo(503, available_gib),
        f"{root}{TASK}/memory.max": "max\n",
        f"{root}{TASK}/memory.current": str(job_usage_gib * GIB),
        f"{root}{TASK}/memory.events": events(task_kills),
        f"{root}{JOB}/step_batch/memory.max": "max\n",
        f"{root}{JOB}/step_batch/memory.current": str(job_usage_gib * GIB),
        f"{root}{JOB}/step_batch/memory.events": events(job_kills),
        f"{root}{JOB}/memory.max": str(232 * GIB),
        f"{root}{JOB}/memory.current": str(job_usage_gib * GIB),
        f"{root}{JOB}/memory.peak": str(40 * GIB),
        f"{root}{JOB}/memory.events": events(job_kills),
        # Above the job: counts every job on the node.
        f"{root}/system.slice/slurmstepd.scope/memory.max": "max\n",
        f"{root}/system.slice/slurmstepd.scope/memory.current": str(400 * GIB),
        f"{root}/system.slice/slurmstepd.scope/memory.events": events(
            neighbour_kills
        ),
    }
    return files


def snapshot(files, rss=None):
    return MEMORY.snapshot(read=files.get, environ={}, clock=lambda: 0.0,
                           child_rss=rss)


class SnapshotTests(unittest.TestCase):

    def test_cgroup_v2_job_limit_is_found_on_an_ancestor(self):
        value = snapshot(cgroup_v2())
        self.assertEqual(value["job_limit"]["path"], JOB)
        self.assertEqual(value["job_limit"]["limit_bytes"], 232 * GIB)
        self.assertEqual(value["node"]["available_bytes"], 3 * GIB)
        self.assertIn("job limit 232.0 GiB (cgroup v2)",
                      MEMORY.summary_line(value))

    def test_cgroup_v1_job_limit_ignores_the_unlimited_sentinel(self):
        base = "/sys/fs/cgroup/memory/slurm/uid_1/job_42"
        files = {
            "/proc/self/cgroup": (
                "12:memory:/slurm/uid_1/job_42/step_batch/task_0\n"
                "4:cpuset:/slurm/uid_1/job_42/step_batch/task_0\n"
            ),
            f"{base}/step_batch/task_0/memory.limit_in_bytes":
                "9223372036854771712\n",
            f"{base}/step_batch/task_0/memory.usage_in_bytes": str(GIB),
            f"{base}/memory.limit_in_bytes": str(232 * GIB),
            f"{base}/memory.usage_in_bytes": str(GIB),
        }
        value = snapshot(files)
        self.assertEqual(value["job_limit"]["version"], 1)
        self.assertEqual(value["job_limit"]["limit_bytes"], 232 * GIB)

    def test_meminfo_without_mem_available_is_estimated(self):
        parsed = MEMORY.parse_meminfo(
            "MemTotal: 1024 kB\nMemFree: 100 kB\nCached: 50 kB\n"
        )
        self.assertEqual(parsed["available_bytes"], 150 * 1024)
        self.assertTrue(parsed["available_is_estimate"])


class AttributionTests(unittest.TestCase):

    def test_kill_counted_in_the_job_is_the_job_limit(self):
        before = snapshot(cgroup_v2())
        after = snapshot(cgroup_v2(job_kills=1, task_kills=1))
        result = MEMORY.attribute_kill(before, after)
        self.assertEqual(result["cause"], "job_memory_limit")
        # v2 counts one kill at every level; it is still one kill.
        self.assertEqual(result["evidence"]["oom_kills_during_stage"], 1)
        self.assertIn("exceeded the job's 232.0 GiB memory limit",
                      result["explanation"])

    def test_other_jobs_kills_on_the_node_are_not_this_jobs(self):
        before = snapshot(cgroup_v2(neighbour_kills=5))
        watch = MEMORY.MemoryWatch(before)
        watch.observe(snapshot(cgroup_v2(available_gib=1)), 120.0)
        after = snapshot(cgroup_v2(neighbour_kills=6, available_gib=2))

        result = MEMORY.attribute_kill(before, after, watch)

        self.assertEqual(result["cause"], "outside_job_limit")
        self.assertIn("did not exceed it", result["explanation"])
        self.assertIn("3.0 GiB available when the stage started",
                      result["explanation"])
        self.assertIn("1.0 GiB at the lowest sample", result["explanation"])

    def test_old_v1_kernel_without_kill_counts_is_undetermined(self):
        base = "/sys/fs/cgroup/memory/slurm/job_42"

        def files(failcnt):
            return {
                "/proc/self/cgroup": "9:memory:/slurm/job_42\n",
                f"{base}/memory.limit_in_bytes": str(232 * GIB),
                f"{base}/memory.usage_in_bytes": str(200 * GIB),
                f"{base}/memory.failcnt": str(failcnt),
                f"{base}/memory.oom_control": "oom_kill_disable 0\nunder_oom 0\n",
            }

        result = MEMORY.attribute_kill(snapshot(files(0)), snapshot(files(7)))
        self.assertEqual(result["cause"], "undetermined")
        self.assertIn("reached its limit 7 time(s)", result["explanation"])

    def test_no_cgroup_limit(self):
        files = {"/proc/meminfo": meminfo(64, 1)}
        result = MEMORY.attribute_kill(snapshot(files), snapshot(files))
        self.assertEqual(result["cause"], "no_job_limit_visible")

    def test_stage_peak_rss_is_reported_only_when_this_stage_set_it(self):
        before = snapshot(cgroup_v2(), rss=30 * GIB)
        unchanged = MEMORY.attribute_kill(
            before, snapshot(cgroup_v2(), rss=30 * GIB)
        )
        self.assertIsNone(unchanged["evidence"]["stage_peak_rss_bytes"])

        grew = MEMORY.attribute_kill(
            before, snapshot(cgroup_v2(), rss=220 * GIB)
        )
        self.assertEqual(grew["evidence"]["stage_peak_rss_bytes"], 220 * GIB)
        self.assertIn("reached 220.0 GiB resident memory", grew["explanation"])


class WarningTests(unittest.TestCase):

    def test_node_short_of_the_job_allowance_is_warned(self):
        warning = MEMORY.low_node_memory_warning(snapshot(cgroup_v2()))
        self.assertIn("node has 3.0 GiB available", warning)

    def test_ample_node_memory_is_quiet(self):
        self.assertIsNone(
            MEMORY.low_node_memory_warning(
                snapshot(cgroup_v2(available_gib=400))
            )
        )


class ImagingMemoryTests(unittest.TestCase):

    def test_padded_grid_matches_casa_composite_sizes(self):
        self.assertEqual(MEMORY.composite_even_at_least(int(8192 * 1.2)), 10000)
        self.assertEqual(MEMORY.composite_even_at_least(int(2048 * 1.2)), 2500)

    def test_selfcal_b4_image_cannot_fit_in_three_gib(self):
        needed = MEMORY.wproject_imaging_memory(8192, 2, 128, 64)
        self.assertEqual(needed["convolution_size"], 10000)
        self.assertGreater(needed["minimum_bytes"], 7.5 * GIB)
        self.assertLess(needed["minimum_bytes"], 8.5 * GIB)
        self.assertGreater(needed["largest_bytes"], 100 * GIB)

        single = MEMORY.wproject_imaging_memory(8192, 2, 128, 1)
        self.assertLess(single["largest_bytes"], needed["largest_bytes"])
        self.assertGreater(single["minimum_bytes"], 3 * GIB)

    def test_without_w_planes_only_images_and_grid_count(self):
        needed = MEMORY.wproject_imaging_memory(2048, 1, 1, 8)
        self.assertEqual(needed["convolution_minimum_bytes"], 0)
        self.assertEqual(needed["minimum_bytes"], needed["largest_bytes"])


class AllocationWarningTests(unittest.TestCase):

    def environment(self, **values):
        return dict({"SLURM_JOB_ID": "42", "SLURM_MEM_PER_NODE": str(232 * 1024)},
                    **values)

    def limited(self, gib, environ):
        files = cgroup_v2()
        files[f"/sys/fs/cgroup{JOB}/memory.max"] = str(gib * GIB)
        return MEMORY.snapshot(read=files.get, environ=environ,
                               clock=lambda: 0.0, child_rss=None)

    def test_limit_far_below_the_slurm_request_is_explained(self):
        warning = MEMORY.allocation_warning(self.limited(3, self.environment()))
        self.assertIn("may use only 3.0 GiB", warning)
        self.assertIn("requested 232.0 GiB", warning)
        self.assertIn("srun", warning)

    def test_limit_outside_slurm_points_to_sbatch(self):
        warning = MEMORY.allocation_warning(self.limited(3, {}))
        self.assertIn("not running inside a Slurm job", warning)

    def test_matching_allocation_is_quiet(self):
        self.assertIsNone(
            MEMORY.allocation_warning(self.limited(232, self.environment()))
        )


class FakeKilledProcess:

    def __init__(self, *args, **kwargs):
        pass

    def wait(self, timeout=None):
        return -9

    def poll(self):
        return -9


class StageRunnerTests(unittest.TestCase):

    def setUp(self):
        # These tests describe a job too small for a large image. The size
        # is stated here instead of being taken from the band default, which
        # is a science setting and free to change.
        patcher = mock.patch.dict(
            PIPELINE.CONFIG.SELFCAL_IMAGING_CONFIG["B4"],
            {"imsize": 8192},
        )
        patcher.start()
        self.addCleanup(patcher.stop)

    def selfcal_runner(self, root):
        target = root / "target.ms"
        target.mkdir()
        (target / "table.dat").touch()
        outdir = root / "Products"
        outdir.mkdir()
        args = PIPELINE.build_parser().parse_args([
            "--outdir", str(outdir),
            "--telescope", "ugmrt",
            "--receiver-band", "B4",
            "--stages", "selfcal",
            "--selfcal-ms", str(target),
            "--dry-run",
        ])
        runner = PIPELINE.PipelineRunner(args)
        runner.args.dry_run = False
        return runner, outdir

    def test_selfcal_that_cannot_fit_is_refused_before_casa_starts(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, outdir = self.selfcal_runner(Path(temporary))
            files = cgroup_v2()
            files[f"/sys/fs/cgroup{JOB}/memory.max"] = str(3 * GIB)
            small = MEMORY.snapshot(
                read=files.get,
                environ={"SLURM_JOB_ID": "42",
                         "SLURM_MEM_PER_NODE": str(232 * 1024)},
                clock=lambda: 0.0,
                child_rss=None,
            )
            launched = mock.Mock()
            with mock.patch.object(
                runner, "command_for_stage", return_value=["casa", "-c", "x"]
            ), mock.patch.object(
                PIPELINE.subprocess, "Popen", launched
            ), mock.patch.object(
                PIPELINE.MEMORY_DIAGNOSTICS, "snapshot", return_value=small
            ):
                with self.assertRaises(RuntimeError) as raised:
                    runner.run_stage("selfcal")

                launched.assert_not_called()
                message = str(raised.exception)
                self.assertIn("cannot run in the 3.0 GiB", message)
                self.assertIn("10000-pixel square grid", message)
                self.assertIn("requested 232.0 GiB", message)
                self.assertNotIn(
                    "selfcal",
                    runner.load_resume_state().get("stages", {}),
                )

                # A dry run reports the same problem without stopping.
                runner.args.dry_run = True
                self.assertEqual(runner.run_stage("selfcal"), 0)

    def test_a_crash_also_carries_memory_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, outdir = self.selfcal_runner(Path(temporary))

            class Segfault(FakeKilledProcess):
                def wait(self, timeout=None):
                    return -11

            files = cgroup_v2(available_gib=180, job_usage_gib=27)
            files[f"/sys/fs/cgroup{JOB}/memory.max"] = str(28 * GIB)
            files[f"/sys/fs/cgroup{JOB}/memory.peak"] = str(28 * GIB)
            tight = MEMORY.snapshot(read=files.get, environ={},
                                    clock=lambda: 0.0, child_rss=None)

            with mock.patch.object(
                runner, "command_for_stage", return_value=["casa", "-c", "x"]
            ), mock.patch.object(
                PIPELINE.subprocess, "Popen", Segfault
            ), mock.patch.object(
                PIPELINE.MEMORY_DIAGNOSTICS, "snapshot", return_value=tight
            ):
                with self.assertRaises(RuntimeError) as raised:
                    runner.run_stage("selfcal")

            message = str(raised.exception)
            self.assertIn("SIGSEGV", message)
            self.assertIn("Memory at the time:", message)
            self.assertIn("job limit 28.0 GiB", message)
            self.assertIn("no out-of-memory kill counted", message)
            written = list(
                (outdir / "diagnostics").glob("selfcal_sigkill_*.json")
            )
            self.assertEqual(len(written), 1)
            self.assertEqual(json.loads(written[0].read_text())["signal"],
                             "SIGSEGV")

    def test_an_interrupted_stage_is_not_a_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, outdir = self.selfcal_runner(Path(temporary))

            class Interrupted(FakeKilledProcess):
                def wait(self, timeout=None):
                    return PIPELINE.STAGE_INTERRUPTED_EXIT_CODE

            with mock.patch.object(
                runner, "command_for_stage", return_value=["casa", "-c", "x"]
            ), mock.patch.object(
                PIPELINE.subprocess, "Popen", Interrupted
            ), mock.patch.object(
                PIPELINE.MEMORY_DIAGNOSTICS,
                "snapshot",
                return_value=snapshot(cgroup_v2(available_gib=400)),
            ):
                with self.assertRaises(PIPELINE.StageInterrupted) as raised:
                    runner.run_stage("selfcal")

            self.assertIn("interrupted from outside", str(raised.exception))
            entry = runner.load_resume_state()["stages"]["selfcal"]
            self.assertEqual(entry["status"], "interrupted")

    def test_sigkill_message_and_evidence_name_the_cause(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, outdir = self.selfcal_runner(Path(temporary))

            before = snapshot(cgroup_v2())
            after = snapshot(cgroup_v2(available_gib=2))
            with mock.patch.object(
                runner, "command_for_stage", return_value=["casa", "-c", "x"]
            ), mock.patch.object(
                PIPELINE.subprocess, "Popen", FakeKilledProcess
            ), mock.patch.object(
                PIPELINE.MEMORY_DIAGNOSTICS,
                "snapshot",
                side_effect=[before, after],
            ):
                with self.assertRaises(RuntimeError) as raised:
                    runner.run_stage("selfcal")

            message = str(raised.exception)
            self.assertIn("terminated by SIGKILL", message)
            self.assertIn("No out-of-memory kill was recorded", message)
            self.assertIn("Memory evidence:", message)

            written = list((outdir / "diagnostics").glob("selfcal_sigkill_*.json"))
            self.assertEqual(len(written), 1)
            record = json.loads(written[0].read_text())
            self.assertEqual(record["attribution"]["cause"], "outside_job_limit")

            entry = runner.load_resume_state()["stages"]["selfcal"]
            self.assertEqual(entry["memory_kill"]["cause"], "outside_job_limit")


if __name__ == "__main__":
    unittest.main()


class CasaResourceBudgetTests(unittest.TestCase):
    """CASA is told the job's limit, not the node's free memory."""

    def runner(self, temporary):
        outdir = Path(temporary) / "Products"
        outdir.mkdir()
        args = PIPELINE.build_parser().parse_args([
            "--outdir", str(outdir),
            "--telescope", "ugmrt",
            "--receiver-band", "B4",
            "--stages", "selfcal",
            "--selfcal-ms", str(Path(temporary) / "target.ms"),
            "--dry-run",
        ])
        return PIPELINE.PipelineRunner(args), outdir

    def test_the_budget_is_a_fraction_of_the_job_limit(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, outdir = self.runner(temporary)
            files = cgroup_v2(available_gib=176)
            files[f"/sys/fs/cgroup{JOB}/memory.max"] = str(28 * GIB)
            tight = MEMORY.snapshot(read=files.get, environ={},
                                    clock=lambda: 0.0, child_rss=None)

            with mock.patch.object(
                PIPELINE.MEMORY_DIAGNOSTICS, "snapshot", return_value=tight
            ):
                path = runner.casa_resource_rc("selfcal", 4)

            self.assertIsNotNone(path)
            written = Path(path).read_text()
            # Half of 28 GiB, in megabytes, not the node's 176 GiB.
            self.assertIn("system.resources.memory: 14336", written)
            self.assertIn("system.resources.cores: 4", written)

    def test_the_run_file_comes_before_the_users_own(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, _ = self.runner(temporary)
            order = runner.casa_resource_file_list(
                "/run/selfcal_resources.rc",
                {"HOME": "/home/observer",
                 "CASAPATH": "/opt/casa linux site host"},
            ).split(":")

            self.assertEqual(order[0], "/run/selfcal_resources.rc")
            self.assertEqual(order[1], "/home/observer/.casarc")
            # casacore's own order, so nothing the user set is lost.
            self.assertEqual(order[3], "/home/observer/.aipsrc")
            self.assertEqual(
                order[4:],
                [
                    "/opt/casa/.aipsrc",
                    "/opt/casa/linux/site/host/aipsrc",
                    "/opt/casa/linux/site/aipsrc",
                    "/opt/casa/linux/aipsrc",
                ],
            )

    def test_an_existing_setting_is_kept_behind_the_run_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, _ = self.runner(temporary)
            self.assertEqual(
                runner.casa_resource_file_list(
                    "/run/rc", {"CASARCFILES": "/custom/rc"}
                ),
                "/run/rc:/custom/rc",
            )

    def test_no_limit_or_a_disabled_fraction_leaves_casa_alone(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, _ = self.runner(temporary)
            without_limit = MEMORY.snapshot(
                read={"/proc/meminfo": meminfo(64, 60)}.get,
                environ={}, clock=lambda: 0.0, child_rss=None,
            )
            without_limit["node"] = None
            with mock.patch.object(
                PIPELINE.MEMORY_DIAGNOSTICS,
                "snapshot",
                return_value=without_limit,
            ):
                self.assertIsNone(runner.casa_resource_rc("selfcal", 4))

            with mock.patch.dict(PIPELINE.RUN_CONFIG,
                                 {"casa_memory_fraction": 0}):
                self.assertIsNone(runner.casa_resource_rc("selfcal", 4))


class CasaLogCaptureTests(StageRunnerTests):
    """A crash inside CASA is reported with CASA's own account of it."""

    def crash(self, runner, returncode=-11):
        class Crashed(FakeKilledProcess):
            def wait(self, timeout=None):
                return returncode

            def poll(self):
                return returncode

        with mock.patch.object(
            runner, "command_for_stage", return_value=["casa", "-c", "x"]
        ), mock.patch.object(
            PIPELINE.subprocess, "Popen", Crashed
        ), mock.patch.object(
            PIPELINE.MEMORY_DIAGNOSTICS,
            "snapshot",
            return_value=snapshot(cgroup_v2()),
        ):
            with self.assertRaises(RuntimeError) as raised:
                runner.run_stage("selfcal")
        return str(raised.exception)

    def written_record(self, outdir):
        written = list((outdir / "diagnostics").glob("selfcal_sigkill_*.json"))
        self.assertEqual(len(written), 1)
        return json.loads(written[0].read_text())

    def test_the_log_of_the_crashing_run_is_kept_and_shown(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, outdir = self.selfcal_runner(Path(temporary))
            log = outdir / "casa-20260924-120000.log"
            log.write_text(
                "".join(f"line {number}\n" for number in range(300))
                + "Memory used in gridding function = 14000 MB "
                  "from maximum 14336 MB\n"
            )

            message = self.crash(runner)

            record = self.written_record(outdir)["casa_log"]
            self.assertEqual(Path(record["path"]).resolve(), log.resolve())
            # Bounded, and it is the *end* of the log that is kept.
            self.assertEqual(len(record["lines"]), 200)
            self.assertIn(
                "Memory used in gridding function = 14000 MB "
                "from maximum 14336 MB",
                record["lines"][-1],
            )
            self.assertIn(log.name, message)

    def test_a_log_from_an_earlier_run_is_not_offered_as_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, outdir = self.selfcal_runner(Path(temporary))
            stale = outdir / "casa-20200101-000000.log"
            stale.write_text("from a run that finished last year\n")
            old = time.time() - 86400
            os.utime(stale, (old, old))

            self.crash(runner)

            self.assertIsNone(self.written_record(outdir)["casa_log"])

    def test_the_newest_log_wins_and_a_kill_carries_it_too(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, outdir = self.selfcal_runner(Path(temporary))
            for name, mtime in (
                ("casa-old.log", time.time() - 30),
                ("casa-new.log", time.time()),
            ):
                path = outdir / name
                path.write_text(f"{name}\n")
                os.utime(path, (mtime, mtime))

            self.crash(runner, returncode=-9)

            record = self.written_record(outdir)["casa_log"]
            self.assertEqual(
                Path(record["path"]).resolve(),
                (outdir / "casa-new.log").resolve(),
            )

    def test_an_unreadable_log_never_hides_the_real_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, outdir = self.selfcal_runner(Path(temporary))
            with mock.patch.object(
                PIPELINE.PipelineRunner,
                "casa_log_tail",
                side_effect=OSError("gone"),
            ):
                with self.assertRaises(OSError):
                    self.crash(runner)

            # The capture is best-effort: a log that cannot be read leaves
            # the stage failure itself intact.
            runner.casa_log_tail = lambda *a, **k: None
            self.assertIn("SIGSEGV", self.crash(runner))


class PlannedImagingFootprintTests(StageRunnerTests):
    """The pre-flight figure is what CASA will do, not a range."""

    def preflight(self, limit_gib, runner=None, root=None):
        files = cgroup_v2(available_gib=180, job_usage_gib=2)
        files[f"/sys/fs/cgroup{JOB}/memory.max"] = str(int(limit_gib * GIB))
        snapshot_value = MEMORY.snapshot(
            read=files.get, environ={}, clock=lambda: 0.0, child_rss=None,
        )
        printed = io.StringIO()
        with contextlib.redirect_stdout(printed):
            runner.check_stage_memory("selfcal", snapshot_value)
        return printed.getvalue()

    def test_the_plan_follows_the_budget_casa_is_given(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, _ = self.selfcal_runner(Path(temporary))
            requirement = runner.stage_memory_requirement(
                "selfcal",
                MEMORY.snapshot(
                    read=(lambda path, files=cgroup_v2(): files.get(path)),
                    environ={}, clock=lambda: 0.0, child_rss=None,
                ),
            )
        # The same 8192-pixel image plans differently at different budgets,
        # and the budget is half the limit.
        wide = MEMORY.wproject_imaging_memory(
            8192, 2, 128, 8, budget_bytes=176 * GIB
        )
        narrow = MEMORY.wproject_imaging_memory(
            8192, 2, 128, 8, budget_bytes=14 * GIB
        )
        self.assertEqual(wide["planned_planes_per_chunk"], 100)
        self.assertEqual(narrow["planned_planes_per_chunk"], 18)
        self.assertLess(narrow["planned_bytes"], wide["planned_bytes"] / 3)
        self.assertIn("planned_bytes", requirement)

    def test_a_job_that_fits_is_reported_without_a_warning(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, _ = self.selfcal_runner(Path(temporary))
            printed = self.preflight(28, runner)

        self.assertIn("w-planes at a time", printed)
        self.assertNotIn("WARNING: this process may use", printed)

    def test_a_job_too_small_names_both_levers(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, _ = self.selfcal_runner(Path(temporary))
            printed = self.preflight(20, runner)

        self.assertIn("selfcal_imsize", printed)
        self.assertIn("--threads", printed)
        self.assertIn("convolution grid", printed)

    def test_leaving_casas_view_alone_falls_back_to_the_range(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner, _ = self.selfcal_runner(Path(temporary))
            with mock.patch.dict(PIPELINE.RUN_CONFIG,
                                 {"casa_memory_fraction": 0}):
                printed = self.preflight(28, runner)

        self.assertNotIn("w-planes at a time", printed)
        self.assertIn("if CASA plans its w-projection work from a", printed)


class SelfcalImagingOverrideTests(unittest.TestCase):
    """Self-cal imaging can be sized to fit the job that has to run it."""

    def config(self, **yaml_settings):
        module = load("test_radio_continuum_config_selfcal_size", "config.py")
        module.CONFIG_YAML.clear()
        module.CONFIG_YAML.update(yaml_settings)
        return module

    def test_imaging_size_knobs_reach_every_band(self):
        module = self.config(
            selfcal_imsize="4096",
            selfcal_cell="1.0arcsec",
            selfcal_wprojplanes=64,
        )
        module._apply_config_yaml_overrides()

        for band, settings in module.SELFCAL_IMAGING_CONFIG.items():
            with self.subTest(band=band):
                self.assertEqual(settings["imsize"], 4096)
                self.assertEqual(settings["cell"], "1.0arcsec")
                self.assertEqual(settings["wprojplanes"], 64)

    def test_the_published_image_keeps_its_own_size(self):
        module = self.config(selfcal_imsize=4096)
        published = dict(module.FINAL_IMAGING_CONFIG["B4"])
        module._apply_config_yaml_overrides()

        self.assertEqual(module.FINAL_IMAGING_CONFIG["B4"], published)

    def test_untouched_settings_keep_their_band_defaults(self):
        module = self.config(selfcal_imsize=4096)
        cells = {band: settings["cell"]
                 for band, settings in module.SELFCAL_IMAGING_CONFIG.items()}
        module._apply_config_yaml_overrides()

        self.assertEqual(
            {band: settings["cell"]
             for band, settings in module.SELFCAL_IMAGING_CONFIG.items()},
            cells,
        )
        self.assertGreater(len(set(cells.values())), 1)

    def test_the_specific_key_beats_the_flat_one_and_is_recorded(self):
        module = self.config(imsize=2048, selfcal_imsize=4096)
        module._apply_config_yaml_overrides()

        self.assertEqual(module.SELFCAL_IMAGING_CONFIG["B4"]["imsize"], 4096)
        self.assertEqual(module.FINAL_IMAGING_CONFIG["B4"]["imsize"], 2048)

        # The displaced band default is reported the same way any other
        # flattened imaging key is, so self-cal says what it imaged.
        recorded = [
            item for item in module.FLAT_IMAGING_OVERRIDES
            if item["table"] == "SELFCAL_IMAGING_CONFIG" and
            item["band"] == "B4" and item["config_yaml_value"] == 4096
        ]
        self.assertEqual(len(recorded), 1)
        self.assertEqual(recorded[0]["key"], "imsize")

    def test_a_meaningless_size_is_refused_at_load(self):
        for settings in ({"selfcal_imsize": 0},
                         {"selfcal_wprojplanes": "many"}):
            with self.subTest(settings=settings):
                module = self.config(**settings)
                with self.assertRaises(ValueError):
                    module._apply_config_yaml_overrides()


if __name__ == "__main__":
    unittest.main()
