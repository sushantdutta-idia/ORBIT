import contextlib
import importlib.util
import io
import json
import os
import socket
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_pipeline_module():
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_pipeline_run_modes",
        SOURCE_DIR / "pipeline.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PIPELINE = load_pipeline_module()


def completed(returncode=0, stdout="", stderr=""):
    return subprocess.CompletedProcess([], returncode, stdout, stderr)


def shell_with_module(name_that_loads):
    """Run the pipeline's module script in bash with a stand-in `module`."""

    def runner(command, **kwargs):
        definition = (
            "module() { [ \"$1\" = load ] && [ \"$2\" = "
            f"{name_that_loads} ] || return 1; "
            "export PATH=/fake/casa/bin:$PATH; export CASA_FROM_MODULE=yes; }; "
        )
        return subprocess.run(["bash", "-c", definition + command[2]], **kwargs)

    return runner


class ModuleEnvironmentTests(unittest.TestCase):
    def test_the_module_script_loads_and_reports_the_environment(self):
        environment = PIPELINE.module_environment(
            "casa/6.7.0", runner=shell_with_module("casa/6.7.0")
        )
        self.assertEqual(environment["CASA_FROM_MODULE"], "yes")
        self.assertTrue(environment["PATH"].startswith("/fake/casa/bin:"))

    def test_a_module_that_does_not_load_is_named(self):
        with self.assertRaisesRegex(RuntimeError, "module avail casa"):
            PIPELINE.module_environment(
                "casa/9.9", runner=shell_with_module("casa/6.7.0")
            )

    def test_a_failing_shell_is_reported(self):
        def runner(command, **kwargs):
            return subprocess.CompletedProcess(command, 2, b"", b"bash: bad")

        with self.assertRaisesRegex(RuntimeError, "status 2.*bad"):
            PIPELINE.module_environment("casa/6.7.0", runner=runner)

    def test_apply_keeps_old_path_entries_and_skips_shell_internals(self):
        environment = {"PATH": "/opt/venv/bin:/usr/bin", "HOME": "/home/x"}
        changed = PIPELINE.apply_module_environment(
            {
                "PATH": "/casa/bin:/usr/bin",
                "HOME": "/home/x",
                "CASAPATH": "/casa",
                "SHLVL": "2",
                "BASH_FUNC_module%%": "() { :; }",
            },
            environment,
        )
        self.assertEqual(environment["PATH"], "/casa/bin:/usr/bin:/opt/venv/bin")
        self.assertEqual(environment["CASAPATH"], "/casa")
        self.assertNotIn("SHLVL", environment)
        self.assertNotIn("BASH_FUNC_module%%", environment)
        self.assertEqual(sorted(changed), ["CASAPATH", "PATH"])


class RunnerModuleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.casa_dir = self.root / "casa" / "bin"
        self.casa_dir.mkdir(parents=True)
        casa = self.casa_dir / "casa"
        casa.write_text("#!/bin/sh\nexit 0\n")
        casa.chmod(0o755)
        self.measurement_set = self.root / "raw.ms"
        self.measurement_set.mkdir()
        (self.measurement_set / "table.dat").touch()

    def tearDown(self):
        self.temporary.cleanup()

    def runner(self, *extra):
        args = PIPELINE.build_parser().parse_args([
            str(self.measurement_set),
            "--outdir", str(self.root / "Products"),
            "--mpi", "off",
            "--telescope", "ugmrt",
            "--stages", "inspection",
            "--casa", "casa",
            *extra,
        ])
        return PIPELINE.PipelineRunner(args)

    def test_loading_the_module_makes_casa_findable(self):
        runner = self.runner("--casa-module", "casa/6.7.0")
        captured = dict(os.environ, PATH=f"{self.casa_dir}:/usr/bin")
        with mock.patch.dict(os.environ, {"PATH": "/usr/bin:/bin"}, clear=False):
            with mock.patch.object(PIPELINE, "module_environment", return_value=captured) as load:
                with contextlib.redirect_stdout(io.StringIO()) as output:
                    runner.load_casa_module()
                    runner.load_casa_module()
                self.assertEqual(
                    PIPELINE.shutil.which("casa"), str(self.casa_dir / "casa")
                )
        load.assert_called_once_with("casa/6.7.0")
        self.assertIn(str(self.casa_dir / "casa"), output.getvalue())

    def test_no_module_means_nothing_is_loaded(self):
        runner = self.runner()
        with mock.patch.object(PIPELINE, "module_environment") as load:
            runner.load_casa_module()
        load.assert_not_called()

    def test_a_module_failure_stops_a_run_but_only_warns_in_a_dry_run(self):
        failure = RuntimeError("'module load casa/9.9' failed here.")
        with mock.patch.object(PIPELINE, "module_environment", side_effect=failure):
            with self.assertRaisesRegex(RuntimeError, "casa/9.9"):
                self.runner("--casa-module", "casa/9.9").load_casa_module()
            with contextlib.redirect_stdout(io.StringIO()) as output:
                self.runner("--casa-module", "casa/9.9", "--dry-run").load_casa_module()
        self.assertIn("WARNING", output.getvalue())

    def test_missing_casa_message_points_at_casa_module(self):
        runner = self.runner()
        runner.args.casa = str(self.root / "nowhere" / "casa")
        with self.assertRaisesRegex(FileNotFoundError, "casa_module"):
            runner.preflight(["inspection"])


class LockTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.lock = self.root / ".pipeline.lock"
        self.lock.write_text("")

    def tearDown(self):
        self.temporary.cleanup()

    def owner(self, record, **values):
        return PIPELINE.live_lock_owner(
            json.dumps(record), self.lock, host="jupyter-013", **values
        )

    def test_a_run_on_another_machine_holds_the_lock_while_it_refreshes_it(self):
        record = {"pid": 1, "slurm_job_id": "14300001", "host": "compute-104"}
        owner = self.owner(record, job_active=lambda job: True)
        self.assertIn("compute-104 (SLURM job 14300001)", owner)
        # squeue cannot be asked: the fresh lock still decides.
        self.assertIn("compute-104", self.owner(record, job_active=lambda job: None))
        # The job has left the queue: whatever the lock's time, it is free.
        self.assertIsNone(self.owner(record, job_active=lambda job: False))
        # A listed job whose run stopped refreshing the lock holds nothing.
        old = time.time() - PIPELINE.LOCK_STALE_SECONDS - 60
        os.utime(self.lock, (old, old))
        self.assertIsNone(self.owner(record, job_active=lambda job: True))

    def test_a_jupyter_session_job_does_not_hold_a_dead_runs_lock(self):
        # On ilifu the Jupyter server is a SLURM job that lives for days. A
        # run started in it that has died must not keep the folder.
        finished = subprocess.Popen([sys.executable, "-c", "pass"])
        finished.wait()
        record = {"pid": finished.pid, "host": "jupyter-013", "slurm_job_id": "999"}
        self.assertIsNone(self.owner(record, job_active=lambda job: True))

    def test_a_live_process_on_this_host_holds_the_lock(self):
        self.assertIn(
            "on jupyter-013",
            self.owner({"pid": os.getppid(), "host": "jupyter-013"}),
        )
        finished = subprocess.Popen([sys.executable, "-c", "pass"])
        finished.wait()
        self.assertIsNone(self.owner({"pid": finished.pid, "host": "jupyter-013"}))

    def test_another_host_holds_it_only_while_it_refreshes_the_lock(self):
        record = {"pid": 4242, "host": "compute-240"}
        self.assertIn("compute-240", self.owner(record))
        old = time.time() - PIPELINE.LOCK_STALE_SECONDS - 60
        os.utime(self.lock, (old, old))
        self.assertIsNone(self.owner(record))

    def test_empty_old_or_broken_records_hold_nothing(self):
        self.assertIsNone(PIPELINE.live_lock_owner("", self.lock))
        self.assertIsNone(PIPELINE.live_lock_owner("not json", self.lock))
        # Records written before hosts were recorded are left to flock.
        self.assertIsNone(self.owner({"pid": os.getppid()}))

    def test_run_lock_records_the_host_and_is_cleared_on_release(self):
        args = PIPELINE.build_parser().parse_args([
            str(self.root / "raw.ms"), "--outdir", str(self.root),
            "--mpi", "off", "--telescope", "ugmrt", "--stages", "inspection",
        ])
        runner = PIPELINE.PipelineRunner(args)
        runner.acquire_run_lock()
        try:
            record = json.loads(self.lock.read_text())
            self.assertEqual(record["host"], socket.gethostname())
            self.assertEqual(record["pid"], os.getpid())
            self.assertIsNotNone(runner._lock_heartbeat)
        finally:
            runner.release_run_lock()
        self.assertEqual(self.lock.read_text(), "")
        self.assertIsNone(runner._lock_heartbeat)

    def test_run_lock_refuses_a_directory_a_slurm_job_holds(self):
        self.lock.write_text(json.dumps(
            {"pid": 1, "host": "compute-104", "slurm_job_id": "14300001"}
        ))
        args = PIPELINE.build_parser().parse_args([
            str(self.root / "raw.ms"), "--outdir", str(self.root),
            "--mpi", "off", "--telescope", "ugmrt", "--stages", "inspection",
        ])
        runner = PIPELINE.PipelineRunner(args)
        with mock.patch.object(PIPELINE, "slurm_job_active", return_value=True):
            with self.assertRaisesRegex(RuntimeError, "SLURM job 14300001"):
                runner.acquire_run_lock()
        self.assertIsNone(runner.run_lock_handle)


class SlurmRequestTests(unittest.TestCase):
    def test_memory_sizes(self):
        self.assertEqual(PIPELINE.memory_mb("120GB"), 122880)
        self.assertEqual(PIPELINE.memory_mb("120G"), 122880)
        self.assertEqual(PIPELINE.memory_mb("237568"), 237568)
        self.assertEqual(PIPELINE.memory_mb("500M"), 500)
        self.assertEqual(PIPELINE.memory_mb("1T"), 1048576)
        self.assertEqual(PIPELINE.memory_mb("232000+"), 232000)
        self.assertIsNone(PIPELINE.memory_mb("lots"))

    def sinfo(self, text, returncode=0):
        calls = []

        def runner(command, **kwargs):
            calls.append(command)
            return completed(returncode, text)

        return runner, calls

    def test_a_job_larger_than_any_node_is_explained(self):
        runner, calls = self.sinfo("32 237568\n32 237568\n")
        problem = PIPELINE.slurm_request_problem(
            {"partition": "Main", "ntasks": 32, "cpus_per_task": 16, "mem": "120GB"},
            runner=runner,
        )
        self.assertIn("512 CPUs", problem)
        self.assertIn("at most 32 CPUs", problem)
        self.assertEqual(calls[0][-2:], ["-p", "Main"])

    def test_memory_beyond_the_nodes_is_explained(self):
        runner, _ = self.sinfo("32 237568\n")
        problem = PIPELINE.slurm_request_problem(
            {"partition": "Main", "ntasks": 1, "cpus_per_task": 16, "mem": "500GB"},
            runner=runner,
        )
        self.assertIn("232 GB", problem)

    def test_a_job_that_fits_or_an_unknown_answer_passes(self):
        runner, _ = self.sinfo("32 237568\n")
        self.assertIsNone(PIPELINE.slurm_request_problem(
            {"partition": "Main", "ntasks": 1, "cpus_per_task": 16, "mem": "120GB"},
            runner=runner,
        ))
        failing, _ = self.sinfo("", returncode=1)
        self.assertIsNone(PIPELINE.slurm_request_problem(
            {"partition": "Main", "ntasks": 64, "cpus_per_task": 16, "mem": "120GB"},
            runner=failing,
        ))

    def test_launch_options_are_removed_for_the_run_itself(self):
        self.assertEqual(
            PIPELINE.strip_launch_options([
                "--to-stage", "calibration", "--background",
                "--write-sbatch", "job.sbatch", "--submit", "--mpi", "off",
            ]),
            ["--to-stage", "calibration", "--mpi", "off"],
        )
        self.assertEqual(
            PIPELINE.strip_launch_options(["--write-sbatch", "--dry-run"]),
            ["--dry-run"],
        )

    def test_background_and_submit_cannot_be_combined(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                PIPELINE.build_parser().parse_args(["--background", "--submit"])


class LaunchTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.outdir = self.root / "Products"
        self.environment = mock.patch.dict(os.environ, {}, clear=False)
        self.environment.start()
        os.environ.pop("SLURM_JOB_ID", None)
        self.processes = []

    def tearDown(self):
        for process in self.processes:
            if process.poll() is None:
                process.kill()
            process.wait()
        self.environment.stop()
        self.temporary.cleanup()

    def args(self, *extra):
        return PIPELINE.build_parser().parse_args([
            "--outdir", str(self.outdir), "--mpi", "off", *extra,
        ])

    def fake_popen(self, code):
        self.started = []

        def popen(command, **kwargs):
            self.started.append((command, kwargs))
            process = subprocess.Popen(
                [sys.executable, "-c", code],
                stdin=kwargs["stdin"],
                stdout=kwargs["stdout"],
                stderr=kwargs["stderr"],
                start_new_session=True,
            )
            self.processes.append(process)
            return process

        return popen

    def test_background_run_is_detached_and_reported(self):
        popen = self.fake_popen("import time; print('working'); time.sleep(30)")
        with contextlib.redirect_stdout(io.StringIO()) as output:
            status = PIPELINE.start_in_background(
                self.args("--to-stage", "calibration", "--background"),
                ["--to-stage", "calibration", "--background"],
                popen=popen,
                wait_seconds=1,
            )
        command, options = self.started[0]
        self.assertEqual(status, 0)
        self.assertEqual(command[1:3], ["-u", str(SOURCE_DIR / "pipeline.py")])
        self.assertEqual(command[3:], ["--to-stage", "calibration"])
        self.assertTrue(options["start_new_session"])
        self.assertEqual(options["stdin"], subprocess.DEVNULL)
        self.assertIsNone(self.processes[0].poll())
        text = output.getvalue()
        self.assertIn(f"as process {self.processes[0].pid}", text)
        self.assertIn(str(self.outdir / "pipeline_run.log"), text)
        self.assertIn("started on", (self.outdir / "pipeline_run.log").read_text())

    def test_a_run_that_stops_at_once_shows_its_log(self):
        popen = self.fake_popen("import sys; print('ERROR: no CASA here'); sys.exit(3)")
        with contextlib.redirect_stdout(io.StringIO()) as output:
            status = PIPELINE.start_in_background(
                self.args(), [], popen=popen, wait_seconds=10
            )
        self.assertEqual(status, 3)
        self.assertIn("no CASA here", output.getvalue())

    def test_background_from_a_jupyter_session_that_is_a_slurm_job(self):
        os.environ["SLURM_JOB_ID"] = "14290001"
        popen = self.fake_popen("import time; time.sleep(30)")
        with contextlib.redirect_stdout(io.StringIO()) as output:
            status = PIPELINE.start_in_background(
                self.args(), [], popen=popen, wait_seconds=1
            )
        self.assertEqual(status, 0)
        self.assertIn("itself SLURM job 14290001", output.getvalue())

    def test_background_refuses_beside_a_live_run(self):
        self.outdir.mkdir()
        (self.outdir / ".pipeline.lock").write_text(json.dumps(
            {"pid": 1, "host": "compute-1", "slurm_job_id": "77"}
        ))
        with self.assertRaisesRegex(RuntimeError, "SLURM job 77"):
            PIPELINE.start_in_background(
                self.args(), [], job_active=lambda job: True
            )
        with self.assertRaisesRegex(RuntimeError, "dry-run"):
            PIPELINE.start_in_background(self.args("--dry-run"), [])

    def slurm(self, responses):
        self.calls = []

        def runner(command, **kwargs):
            self.calls.append(command)
            return responses[Path(command[0]).name](command)

        return runner

    def submit(self, *extra, responses, which="/usr/bin/sbatch"):
        job_file = self.root / "job.sbatch"
        settings = {
            "partition": "Main", "ntasks": 1, "cpus_per_task": 16,
            "mem": "120GB", "log_dir": ".", "job_name": "radio_continuum",
            "output": "radio_continuum-%j.out", "error": "radio_continuum-%j.err",
        }
        with mock.patch.dict(PIPELINE.CONFIG.ILIFU_SLURM_CONFIG, settings):
            with contextlib.redirect_stdout(io.StringIO()) as output:
                status = PIPELINE.submit_to_slurm(
                    self.args("--write-sbatch", str(job_file), *extra),
                    runner=self.slurm(responses),
                    which=lambda name: which,
                )
        return status, output.getvalue(), job_file

    def test_submit_writes_checks_and_submits(self):
        status, text, job_file = self.submit(
            "--casa-module", "casa/6.7.0", "--to-stage", "calibration",
            responses={
                "sinfo": lambda command: completed(0, "32 237568\n"),
                "sbatch": lambda command: completed(0, "14300042;cluster\n"),
            },
        )
        self.assertEqual(status, 0)
        self.assertIn("Submitted SLURM job 14300042", text)
        self.assertIn("radio_continuum-14300042.out", text)
        self.assertEqual(self.calls[-1], ["/usr/bin/sbatch", "--parsable", str(job_file)])
        content = job_file.read_text()
        self.assertIn("--casa-module casa/6.7.0", content)
        self.assertIn("--to-stage calibration", content)
        self.assertNotIn("--submit", content)

    def test_submit_stops_before_sbatch_when_no_node_fits(self):
        with self.assertRaisesRegex(RuntimeError, "at most 8 CPUs"):
            self.submit(responses={
                "sinfo": lambda command: completed(0, "8 122880\n"),
                "sbatch": lambda command: self.fail("sbatch must not run"),
            })

    def test_submit_from_a_jupyter_session_that_is_a_slurm_job(self):
        os.environ["SLURM_JOB_ID"] = "14290001"
        status, text, _ = self.submit(responses={
            "sinfo": lambda command: completed(0, "32 237568\n"),
            "sbatch": lambda command: completed(0, "14300043\n"),
        })
        self.assertEqual(status, 0)
        self.assertIn("Submitted SLURM job 14300043", text)

    def test_submit_reports_a_refusal_and_a_missing_sbatch(self):
        status, text, _ = self.submit(responses={
            "sinfo": lambda command: completed(1),
            "sbatch": lambda command: completed(
                1, "", "sbatch: error: Requested node configuration is not available"
            ),
        })
        self.assertEqual(status, 1)
        self.assertIn("Requested node configuration", text)
        self.assertIn("sinfo -s", text)
        with self.assertRaisesRegex(RuntimeError, "--background"):
            self.submit(responses={}, which=None)



ILIFU_IMAGE = (
    "/software/astro/containers/casa-modular/"
    "casa-modular-v6.7.0-py3.10-2025-12-01.sif"
)


def shell_with_ilifu_casa_module():
    """bash with a `module` that does what ilifu's casa/6.7.0.lua does."""

    def runner(command, **kwargs):
        definition = (
            'module() { [ "$2" = casa/6.7.0 ] || return 1; '
            f'casapy() {{ singularity exec {ILIFU_IMAGE} xvfb-run -a ipython "$@"; }}; '
            "}; "
        )
        return subprocess.run(["bash", "-c", definition + command[2]], **kwargs)

    return runner


class ModularCasaTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.measurement_set = self.root / "raw.ms"
        self.measurement_set.mkdir()
        (self.measurement_set / "table.dat").touch()

    def tearDown(self):
        self.temporary.cleanup()

    def runner(self, *extra):
        args = PIPELINE.build_parser().parse_args([
            str(self.measurement_set),
            "--outdir", str(self.root / "Products"),
            "--telescope", "ugmrt",
            "--stages", "inspection",
            *extra,
        ])
        return PIPELINE.PipelineRunner(args)

    def test_the_container_behind_ilifus_casapy_is_found(self):
        self.assertEqual(
            PIPELINE.module_casa_container(
                "casa/6.7.0", runner=shell_with_ilifu_casa_module()
            ),
            f"singularity exec {ILIFU_IMAGE} python3",
        )
        self.assertIsNone(PIPELINE.module_casa_container(
            "casa/9.9", runner=shell_with_ilifu_casa_module()
        ))

    def test_stages_run_with_the_modular_python(self):
        runner = self.runner(
            "--casa-python", f"singularity exec {ILIFU_IMAGE} python3",
            "--mpi", "auto",
        )
        command = runner.command_for_stage("inspection")
        self.assertEqual(
            command[:5],
            ["singularity", "exec", ILIFU_IMAGE, "python3",
             str(SOURCE_DIR / "inspect.py")],
        )
        self.assertNotIn("--nologger", command)
        self.assertFalse(runner.mpi_enabled())
        with self.assertRaisesRegex(RuntimeError, "mpicasa"):
            runner.casa_command(SOURCE_DIR / "imaging.py", mpi=True)

    def test_mpi_on_is_refused_with_a_modular_casa(self):
        runner = self.runner("--casa-python", "python3", "--mpi", "on")
        with mock.patch.object(PIPELINE.platform, "system", return_value="Linux"):
            with self.assertRaisesRegex(RuntimeError, "--mpi on needs mpicasa"):
                runner.preflight(["inspection"])

    def test_preflight_starts_the_modular_python_once(self):
        runner = self.runner("--casa-python", "fake-python", "--mpi", "off")
        with mock.patch.object(
            PIPELINE, "casa_python_version", return_value="6.7.0.31"
        ) as probe:
            with contextlib.redirect_stdout(io.StringIO()) as output:
                runner.preflight(["inspection"])
        probe.assert_called_once_with(["fake-python"])
        self.assertIn("CASA 6.7.0.31 via fake-python", output.getvalue())

    def test_a_module_without_casa_uses_its_container(self):
        runner = self.runner("--casa-module", "casa/6.7.0", "--mpi", "off")
        with mock.patch.dict(os.environ, {"PATH": "/usr/bin:/bin"}):
            with mock.patch.object(
                PIPELINE, "module_environment",
                return_value={"PATH": "/usr/bin:/bin"},
            ), mock.patch.object(
                PIPELINE, "module_casa_container",
                return_value=f"singularity exec {ILIFU_IMAGE} python3",
            ):
                with contextlib.redirect_stdout(io.StringIO()) as output:
                    runner.load_casa_module()
        self.assertEqual(
            runner.args.casa_python, f"singularity exec {ILIFU_IMAGE} python3"
        )
        self.assertIn("no `casa` program", output.getvalue())
        self.assertEqual(
            runner.command_for_stage("inspection")[:4],
            ["singularity", "exec", ILIFU_IMAGE, "python3"],
        )

    def test_casa_python_version(self):
        def answering(command, **kwargs):
            self.assertEqual(command[:2], ["python3", "-c"])
            return completed(0, "some banner\nCASA 6.7.0.31\n")

        self.assertEqual(
            PIPELINE.casa_python_version(["python3"], runner=answering), "6.7.0.31"
        )
        # A Python without CASA is refused, with what it said.
        with self.assertRaisesRegex(RuntimeError, "cannot import casatasks"):
            PIPELINE.casa_python_version([sys.executable, "-I"])

    def test_generated_job_keeps_the_modular_python(self):
        job = self.root / "job.sbatch"
        args = PIPELINE.build_parser().parse_args([
            str(self.measurement_set), "--outdir", str(self.root / "P"),
            "--casa-python", f"singularity exec {ILIFU_IMAGE} python3",
            "--write-sbatch", str(job),
        ])
        with contextlib.redirect_stdout(io.StringIO()):
            PIPELINE.write_sbatch_from_config(args)
        self.assertIn(
            f"--casa-python 'singularity exec {ILIFU_IMAGE} python3'",
            job.read_text(),
        )

    def test_stage_scripts_start_under_plain_python(self):
        # As a modular CASA starts them; inspect.py must not hide the
        # standard inspect module from the script.
        for script in ["fits_ingest.py", "lta_convert.py"]:
            with self.subTest(script=script):
                result = subprocess.run(
                    [sys.executable, str(SOURCE_DIR / script), "--help"],
                    capture_output=True, text=True, cwd=SOURCE_DIR,
                )
                self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
