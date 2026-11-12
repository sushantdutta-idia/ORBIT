"""A missing tool stops the run at the start, not after calibration."""

import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_resume_safety import PIPELINE, make_casa_table


SOURCE_DIR = Path(__file__).resolve().parents[1]


class PreflightFixture(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.ms = self.root / "in.ms"
        make_casa_table(self.ms)
        self.wsclean = self.root / "wsclean"
        self.wsclean.write_text("#!/bin/sh\n")

    def tearDown(self):
        self.temporary.cleanup()

    def runner(self, *extra, dry_run=False, wsclean=None):
        argv = [
            str(self.ms), "--outdir", str(self.root / "P"),
            "--no-copy-input-ms", "--wsclean",
            str(wsclean if wsclean is not None else self.wsclean), *extra,
        ]
        if dry_run:
            argv.append("--dry-run")
        return PIPELINE.PipelineRunner(PIPELINE.build_parser().parse_args(argv))


class RequiredFilesTests(PreflightFixture):

    def required(self, *extra, stages=("imaging",)):
        runner = self.runner(*extra, "--stages", ",".join(stages),
                             dry_run=True)
        return set(runner.required_source_versions(list(stages)))

    def test_ugmrt_never_demands_the_wsclean_files(self):
        required = self.required("--telescope", "ugmrt",
                                 "--receiver-band", "B4")
        self.assertNotIn("wsclean_imaging.py", required)
        self.assertNotIn("meerkat_beam.py", required)

    def test_lofar_needs_the_stage_but_not_the_meerkat_beam(self):
        required = self.required("--telescope", "lofar",
                                 "--receiver-band", "HBA",
                                 "--data-state", "calibrated")
        self.assertIn("wsclean_imaging.py", required)
        self.assertNotIn("meerkat_beam.py", required)

    def test_meerkat_needs_both(self):
        required = self.required("--telescope", "meerkat",
                                 "--receiver-band", "L")
        self.assertLessEqual({"wsclean_imaging.py", "meerkat_beam.py"},
                             required)

    def test_an_unknown_telescope_could_be_either(self):
        required = self.required()
        self.assertLessEqual({"wsclean_imaging.py", "meerkat_beam.py"},
                             required)

    def test_casa_or_no_imaging_needs_neither(self):
        for extra, stages in (
            (("--telescope", "meerkat", "--receiver-band", "L",
              "--imager", "casa"), ("imaging",)),
            (("--telescope", "meerkat", "--receiver-band", "L"),
             ("calibration",)),
        ):
            with self.subTest(extra=extra, stages=stages):
                required = self.required(*extra, stages=stages)
                self.assertNotIn("wsclean_imaging.py", required)
                self.assertNotIn("meerkat_beam.py", required)

    def test_the_shipped_files_satisfy_every_run(self):
        for extra in (
            ("--telescope", "meerkat", "--receiver-band", "L"),
            ("--telescope", "lofar", "--receiver-band", "HBA",
             "--data-state", "calibrated"),
            ("--telescope", "ugmrt", "--receiver-band", "B4"),
            (),
        ):
            with self.subTest(extra=extra):
                runner = self.runner(*extra, "--stages", "imaging",
                                     dry_run=True)
                self.assertEqual(
                    PIPELINE.outdated_source_files(
                        SOURCE_DIR, runner.required_source_versions(
                            ["imaging"])),
                    [])

    def test_a_missing_beam_module_is_named_for_meerkat_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in PIPELINE.SOURCE_INTERFACE_VERSIONS:
                (root / name).write_text(
                    f"PIPELINE_INTERFACE_VERSION = "
                    f"{PIPELINE.SOURCE_INTERFACE_VERSIONS[name]}\n")
            (root / "wsclean_imaging.py").write_text(
                "PIPELINE_INTERFACE_VERSION = 1\n")
            meerkat = self.runner("--telescope", "meerkat",
                                  "--receiver-band", "L",
                                  "--stages", "imaging", dry_run=True)
            problems = PIPELINE.outdated_source_files(
                root, meerkat.required_source_versions(["imaging"]))
            self.assertEqual(problems, ["meerkat_beam.py is missing"])
            lofar = self.runner("--telescope", "lofar", "--receiver-band",
                                "HBA", "--data-state", "calibrated",
                                "--stages", "imaging", dry_run=True)
            self.assertEqual(PIPELINE.outdated_source_files(
                root, lofar.required_source_versions(["imaging"])), [])


class ToolCheckTests(PreflightFixture):

    MEERKAT = ("--telescope", "meerkat", "--receiver-band", "L",
               "--stages", "imaging")

    def completed(self, returncode=0, stderr=""):
        return subprocess.CompletedProcess([], returncode, "", stderr)

    def test_a_dry_run_and_the_casa_backend_check_nothing(self):
        for runner in (
            self.runner(*self.MEERKAT, dry_run=True,
                        wsclean=self.root / "absent"),
            self.runner(*self.MEERKAT, "--imager", "casa",
                        wsclean=self.root / "absent"),
            self.runner("--telescope", "ugmrt", "--receiver-band", "B4",
                        "--stages", "imaging", wsclean=self.root / "absent"),
        ):
            with mock.patch.object(PIPELINE.subprocess, "run") as run:
                runner.check_wsclean_tools(["imaging"])
                run.assert_not_called()

    def test_stages_without_imaging_check_nothing(self):
        runner = self.runner("--telescope", "meerkat", "--receiver-band", "L",
                             "--stages", "calibration",
                             wsclean=self.root / "absent")
        runner.check_wsclean_tools(["calibration"])

    def test_a_missing_wsclean_stops_a_meerkat_run_and_offers_casa(self):
        runner = self.runner(*self.MEERKAT, wsclean=self.root / "absent")
        with self.assertRaises(FileNotFoundError) as raised:
            runner.check_wsclean_tools(["imaging"])
        message = str(raised.exception)
        self.assertIn("meerkat imaging uses WSClean by default", message)
        self.assertIn("--imager casa", message)

    def test_a_missing_wsclean_does_not_offer_casa_for_lofar(self):
        runner = self.runner(
            "--telescope", "lofar", "--receiver-band", "HBA",
            "--data-state", "calibrated", "--stages", "imaging",
            "--lofar-beam-mode", "full", wsclean=self.root / "absent")
        with self.assertRaises(FileNotFoundError) as raised:
            runner.check_wsclean_tools(["imaging"])
        self.assertNotIn("--imager casa", str(raised.exception))

    def test_katbeam_and_astropy_are_checked_in_the_stage_python(self):
        runner = self.runner(*self.MEERKAT)
        with mock.patch.object(
            PIPELINE.subprocess, "run", return_value=self.completed()
        ) as run:
            runner.check_wsclean_tools(["imaging"])
        command = run.call_args[0][0]
        self.assertEqual(command[0], runner.args.python)
        self.assertIn("katbeam", command[-1])
        self.assertIn("astropy", command[-1])

    def test_a_missing_katbeam_stops_the_run_at_the_start(self):
        runner = self.runner(*self.MEERKAT)
        with mock.patch.object(
            PIPELINE.subprocess, "run",
            return_value=self.completed(
                1, "Traceback...\nModuleNotFoundError: No module named "
                   "'katbeam'"),
        ):
            with self.assertRaises(RuntimeError) as raised:
                runner.check_wsclean_tools(["imaging"])
        message = str(raised.exception)
        self.assertIn("No module named 'katbeam'", message)
        self.assertIn("--python", message)
        self.assertIn("--imager casa", message)

    def test_an_unrunnable_python_is_reported(self):
        runner = self.runner(*self.MEERKAT)
        with mock.patch.object(PIPELINE.subprocess, "run",
                               side_effect=OSError("no such file")):
            with self.assertRaisesRegex(RuntimeError, "Could not run"):
                runner.check_wsclean_tools(["imaging"])

    def test_lofar_does_not_need_katbeam(self):
        runner = self.runner(
            "--telescope", "lofar", "--receiver-band", "HBA",
            "--data-state", "calibrated", "--stages", "imaging",
            "--lofar-beam-mode", "full")
        with mock.patch.object(PIPELINE.subprocess, "run") as run:
            runner.check_wsclean_tools(["imaging"])
            run.assert_not_called()

    def test_preflight_runs_the_check(self):
        runner = self.runner(*self.MEERKAT, wsclean=self.root / "absent")
        with self.assertRaises(FileNotFoundError):
            runner.preflight(["imaging"])


if __name__ == "__main__":
    unittest.main()
