import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_flagging_selectors import load_flagging_module_without_casa


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_pipeline():
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_pipeline_recipe_refresh",
        SOURCE_DIR / "pipeline.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PIPELINE = load_pipeline()
FLAGGING = load_flagging_module_without_casa()


def make_ms(path):
    path.mkdir(parents=True)
    (path / "table.dat").write_bytes(b"stable")
    for name in (
        "ANTENNA",
        "DATA_DESCRIPTION",
        "OBSERVATION",
        "POLARIZATION",
        "SPECTRAL_WINDOW",
    ):
        (path / name).mkdir()
        (path / name / "table.dat").write_bytes(b"stable")
    return path


def refresh_flagging_provenance(outdir, ms, fingerprint):
    """Re-enter a completed flagging stage with the real flagging.py code."""

    flagger = FLAGGING.CASAFlagger.__new__(FLAGGING.CASAFlagger)
    flagger.root = outdir / "flagging"
    flagger.donefile = flagger.root / "flagging.done"
    flagger.ms = str(ms)
    flagger.telescope = "ugmrt"
    flagger.band = "B4"
    flagger.provenance_fingerprint = fingerprint
    with mock.patch.object(flagger, "regenerate_qa_plots"):
        flagger.run()


class CalibrationRecipeAcrossFlaggingRefreshTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        self.ms = make_ms(root / "input.ms")
        self.outdir = root / "Products"
        (self.outdir / "inspection").mkdir(parents=True)
        (self.outdir / "flagging").mkdir()
        (self.outdir / "inspection" / "inspection_report.json").write_text(
            json.dumps({"telescope": "ugmrt", "band": "B4"})
        )
        (self.outdir / "flagging" / "flagging.done").touch()

        # Written exactly as flagging.py writes its post-flag summary.
        writer = FLAGGING.CASAFlagger.__new__(FLAGGING.CASAFlagger)
        writer.save_json(
            {
                "flagged": 1234.0,
                "total": 5000.0,
                "field": {"G0231-0533": {"flagged": 1000.0, "total": 4000.0}},
                "telescope": "ugmrt",
                "band": "B4",
                "measurement_set": str(self.ms.resolve()),
                "provenance_fingerprint": "a" * 64,
            },
            self.post_flag_path,
        )

    def tearDown(self):
        self.temporary.cleanup()

    @property
    def post_flag_path(self):
        return self.outdir / "flagging" / "postflag_summary.json"

    def runner(self):
        args = PIPELINE.build_parser().parse_args([
            str(self.ms),
            "--outdir", str(self.outdir),
            "--telescope", "ugmrt",
            "--receiver-band", "B4",
            "--stages", "calibration",
            "--flux-calibrator", "3C147",
            "--phase-calibrator", "0240-231",
            "--dry-run",
            "--no-copy-input-ms",
        ])
        return PIPELINE.PipelineRunner(args)

    def failed_attempt_state(self, recipe):
        return {
            "schema_version": 1,
            "created_utc": PIPELINE.utc_now(),
            "stages": {
                "calibration": {
                    "status": "failed",
                    "attempt_stage_fingerprint": "0" * 64,
                    "attempt_recipe_fingerprint": recipe,
                }
            },
        }

    def test_failed_calibration_resumes_after_a_flagging_refresh(self):
        recorded = self.runner().stage_recipe_fingerprint("calibration")

        refresh_flagging_provenance(self.outdir, self.ms, "b" * 64)
        refresh_flagging_provenance(self.outdir, self.ms, "c" * 64)
        refreshed = json.loads(self.post_flag_path.read_text())
        self.assertEqual(len(refreshed["provenance_migration_history"]), 2)

        runner = self.runner()
        self.assertNotEqual(
            runner.computed_stage_recipe_fingerprint("calibration"),
            recorded,
        )
        runner.resume_state = self.failed_attempt_state(recorded)

        self.assertFalse(runner.should_skip_stage("calibration"))

        entry = runner.resume_state["stages"]["calibration"]
        self.assertEqual(entry["attempt_recipe_fingerprint"], recorded)
        self.assertIn(
            "provenance refreshes",
            entry["implementation_resume_history"][-1]["reason"],
        )
        # Calibration checkpoints bound to the recorded recipe stay valid.
        command = runner.command_for_stage("calibration")
        self.assertEqual(
            command[command.index("--recipe-fingerprint") + 1],
            recorded,
        )

    def test_science_change_alongside_a_refresh_still_fails_closed(self):
        recorded = self.runner().stage_recipe_fingerprint("calibration")

        refresh_flagging_provenance(self.outdir, self.ms, "b" * 64)
        report = json.loads(self.post_flag_path.read_text())
        report["flagged"] = 2345.0
        self.post_flag_path.write_text(json.dumps(report, indent=4))

        runner = self.runner()
        runner.resume_state = self.failed_attempt_state(recorded)
        with self.assertRaisesRegex(
            RuntimeError,
            "science-recipe fingerprint changed",
        ):
            runner.should_skip_stage("calibration")

    def test_changed_calibration_setting_after_a_refresh_fails_closed(self):
        recorded = self.runner().stage_recipe_fingerprint("calibration")
        refresh_flagging_provenance(self.outdir, self.ms, "b" * 64)

        runner = self.runner()
        runner.resume_state = self.failed_attempt_state(recorded)
        changed = dict(PIPELINE.CONFIG.CALIBRATION_TASK_CONFIG)
        changed["test_only_setting"] = 1
        with mock.patch.object(
            PIPELINE.CONFIG,
            "CALIBRATION_TASK_CONFIG",
            changed,
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "science-recipe fingerprint changed",
            ):
                runner.should_skip_stage("calibration")


class AdoptedReportProvenanceTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        self.ms = make_ms(root / "input.ms")
        self.outdir = root / "Products"
        (self.outdir / "flagging").mkdir(parents=True)
        (self.outdir / "flagging" / "flagging.done").touch()
        args = PIPELINE.build_parser().parse_args([
            str(self.ms),
            "--outdir", str(self.outdir),
            "--telescope", "ugmrt",
            "--receiver-band", "B4",
            "--stages", "flagging",
            "--dry-run",
            "--no-copy-input-ms",
        ])
        self.runner = PIPELINE.PipelineRunner(args)
        self.current = self.runner.stage_fingerprint("flagging")
        (self.outdir / "flagging" / "postflag_summary.json").write_text(
            json.dumps({
                "telescope": "ugmrt",
                "band": "B4",
                "measurement_set": str(self.ms.resolve()),
                "provenance_fingerprint": "a" * 64,
            })
        )

    def tearDown(self):
        self.temporary.cleanup()

    def state(self, product, history):
        self.runner.resume_state = {
            "schema_version": 1,
            "stages": {
                "flagging": {
                    "status": "completed",
                    "product_stage_fingerprint": product,
                    "implementation_adoption_history": history,
                }
            },
        }

    def test_adopted_report_is_current_on_the_next_run(self):
        self.state(self.current, [
            {
                "previous_stage_fingerprint": "a" * 64,
                "current_stage_fingerprint": "b" * 64,
            },
            {
                "previous_stage_fingerprint": "b" * 64,
                "current_stage_fingerprint": self.current,
            },
        ])
        valid, reason = self.runner.validate_stage_completion("flagging")
        self.assertTrue(valid, reason)
        self.assertTrue(self.runner.should_skip_stage("flagging"))

    def test_unadopted_report_is_still_stale(self):
        self.state(self.current, [
            {
                "previous_stage_fingerprint": "d" * 64,
                "current_stage_fingerprint": self.current,
            },
        ])
        valid, reason = self.runner.validate_stage_completion("flagging")
        self.assertFalse(valid)
        self.assertIn("provenance fingerprint", reason)

    def test_adoption_chain_must_reach_the_recorded_product(self):
        self.state("e" * 64, [
            {
                "previous_stage_fingerprint": "a" * 64,
                "current_stage_fingerprint": self.current,
            },
        ])
        self.assertFalse(
            self.runner.report_provenance_is_current("flagging", "a" * 64)
        )


class SelfcalRecipeTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        self.target_ms = make_ms(root / "targets" / "G0231-0533.ms")
        self.outdir = root / "Products"
        (self.outdir / "calibration").mkdir(parents=True)
        self.calibration_report = (
            self.outdir / "calibration" / "calibration_report.json"
        )
        self.write_calibration_report("a" * 64, solint="20s")

    def tearDown(self):
        self.temporary.cleanup()

    def write_calibration_report(self, provenance, solint):
        self.calibration_report.write_text(json.dumps({
            "provenance_fingerprint": provenance,
            "telescope": "ugmrt",
            "band": "B4",
            "calibration_parameters": {"phase_solint": solint},
        }))

    def runner(self):
        args = PIPELINE.build_parser().parse_args([
            "--outdir", str(self.outdir),
            "--telescope", "ugmrt",
            "--receiver-band", "B4",
            "--stages", "selfcal",
            "--selfcal-ms", str(self.target_ms),
            "--target", "G0231-0533",
            "--dry-run",
        ])
        return PIPELINE.PipelineRunner(args)

    def test_recipe_is_not_the_implementation_fingerprint(self):
        runner = self.runner()
        self.assertNotEqual(
            runner.stage_recipe_fingerprint("selfcal"),
            runner.stage_fingerprint("selfcal"),
        )

    def test_recipe_ignores_upstream_implementation_provenance(self):
        before = self.runner().stage_recipe_fingerprint("selfcal")
        self.write_calibration_report("b" * 64, solint="20s")
        self.assertEqual(self.runner().stage_recipe_fingerprint("selfcal"), before)

    def test_recipe_tracks_upstream_science_and_selfcal_settings(self):
        before = self.runner().stage_recipe_fingerprint("selfcal")

        self.write_calibration_report("a" * 64, solint="60s")
        self.assertNotEqual(
            self.runner().stage_recipe_fingerprint("selfcal"),
            before,
        )
        self.write_calibration_report("a" * 64, solint="20s")

        changed = dict(PIPELINE.CONFIG.SELFCAL_TASK_CONFIG)
        changed["amp_start_dynamic_range"] = 30.0
        with mock.patch.object(PIPELINE.CONFIG, "SELFCAL_TASK_CONFIG", changed):
            self.assertNotEqual(
                self.runner().stage_recipe_fingerprint("selfcal"),
                before,
            )

        other_target = self.runner()
        other_target.args.target = "another_target"
        self.assertNotEqual(
            other_target.stage_recipe_fingerprint("selfcal"),
            before,
        )

    def test_killed_selfcal_resumes_after_a_code_only_update(self):
        runner = self.runner()
        recipe = runner.stage_recipe_fingerprint("selfcal")
        runner.resume_state = {
            "schema_version": 1,
            "stages": {
                "selfcal": {
                    "status": "failed",
                    "attempt_stage_fingerprint": "0" * 64,
                    "attempt_recipe_fingerprint": recipe,
                }
            },
        }

        self.assertFalse(runner.should_skip_stage("selfcal"))
        entry = runner.resume_state["stages"]["selfcal"]
        self.assertEqual(len(entry["implementation_resume_history"]), 1)

    def failed_attempt(self, runner, recipe):
        runner.resume_state = {
            "schema_version": 1,
            "stages": {
                "selfcal": {
                    "status": "failed",
                    "attempt_stage_fingerprint": "0" * 64,
                    "attempt_recipe_fingerprint": recipe,
                }
            },
        }
        partial = self.outdir / "selfcal" / "images" / "initial.image.tt0"
        partial.mkdir(parents=True)
        return partial

    def test_attempt_under_another_recipe_is_set_aside_and_restarted(self):
        runner = self.runner()
        runner.args.dry_run = False
        self.failed_attempt(runner, "0" * 64)

        self.assertFalse(runner.should_skip_stage("selfcal"))

        archived = self.outdir / "selfcal.superseded.1"
        self.assertFalse((self.outdir / "selfcal").exists())
        self.assertTrue((archived / "images" / "initial.image.tt0").is_dir())

        entry = runner.resume_state["stages"]["selfcal"]
        self.assertEqual(len(entry["superseded_attempts"]), 1)
        self.assertIn("predates", entry["superseded_attempts"][0]["reason"])
        self.assertEqual(
            Path(entry["superseded_attempts"][0]["archived_directory"]),
            archived.resolve(),
        )
        self.assertEqual(
            entry["attempt_recipe_fingerprint"],
            runner.stage_recipe_fingerprint("selfcal"),
        )
        self.assertNotIn("implementation_resume_history", entry)

    def test_changed_selfcal_setting_restarts_and_keeps_earlier_archives(self):
        (self.outdir / "selfcal.superseded.1").mkdir()
        runner = self.runner()
        recipe = runner.stage_recipe_fingerprint("selfcal")
        self.failed_attempt(runner, recipe)

        changed = dict(PIPELINE.CONFIG.SELFCAL_TASK_CONFIG)
        changed["phase_improvement_stop"] = 0.05
        with mock.patch.object(PIPELINE.CONFIG, "SELFCAL_TASK_CONFIG", changed):
            fresh = self.runner()
            fresh.args.dry_run = False
            fresh.resume_state = runner.resume_state
            self.assertFalse(fresh.should_skip_stage("selfcal"))

        self.assertTrue((self.outdir / "selfcal.superseded.2" / "images").is_dir())
        entry = fresh.resume_state["stages"]["selfcal"]
        self.assertIn("changed", entry["superseded_attempts"][0]["reason"])

    def test_a_failed_imaging_attempt_is_set_aside_too(self):
        runner = self.runner()
        runner.args.dry_run = False
        partial = self.outdir / "imaging" / "images"
        partial.mkdir(parents=True)
        runner.resume_state = {
            "schema_version": 1,
            "stages": {
                "imaging": {
                    "status": "failed",
                    "attempt_stage_fingerprint": "0" * 64,
                    "attempt_recipe_fingerprint": "0" * 64,
                    "started_utc": "2000-01-01T00:00:00+00:00",
                }
            },
        }

        self.assertFalse(runner.should_skip_stage("imaging"))

        self.assertFalse((self.outdir / "imaging").exists())
        self.assertTrue((self.outdir / "imaging.superseded.1" / "images").is_dir())
        entry = runner.resume_state["stages"]["imaging"]
        self.assertEqual(len(entry["superseded_attempts"]), 1)

    def test_dry_run_reports_but_does_not_move_products(self):
        runner = self.runner()
        partial = self.failed_attempt(runner, "0" * 64)
        runner.args.dry_run = True

        self.assertFalse(runner.should_skip_stage("selfcal"))
        self.assertTrue(partial.is_dir())
        self.assertFalse((self.outdir / "selfcal.superseded.1").exists())

    def test_run_specific_calibration_fields_do_not_change_the_recipe(self):
        def write(**extra):
            report = {
                "provenance_fingerprint": "a" * 64,
                "telescope": "ugmrt",
                "band": "B4",
                "calibration_parameters": {"phase_solint": "20s"},
                "calibration_tables": {"phase": "/cal/phase.cal"},
                "split_target_identities": {"targets": {"T": {"path": "/t.ms"}}},
                "flux_model_provenance": {
                    "identity": {
                        "mode": "builtin-setjy",
                        "flux_standard": "Perley-Butler 2017",
                        "pipeline_recipe_fingerprint": "r" * 64,
                    },
                },
            }
            report.update(extra)
            self.calibration_report.write_text(json.dumps(report))

        write()
        before = self.runner().stage_recipe_fingerprint("selfcal")

        write(
            recipe_fingerprint="s" * 64,
            measurement_set_stat_sha256="m" * 64,
            checkpoints={"setjy": True},
            qa_plot_status="/elsewhere/qa_plot_status.json",
            flux_model_provenance={
                "identity": {
                    "mode": "builtin-setjy",
                    "flux_standard": "Perley-Butler 2017",
                    "pipeline_recipe_fingerprint": "t" * 64,
                },
                "identity_migration_history": [{"migrated_utc": "now"}],
            },
        )
        self.assertEqual(self.runner().stage_recipe_fingerprint("selfcal"), before)

        write(calibration_tables={"phase": "/cal/phase.round2.cal"})
        self.assertNotEqual(
            self.runner().stage_recipe_fingerprint("selfcal"),
            before,
        )

        write(split_target_identities={"targets": {"T": {"path": "/new.ms"}}})
        self.assertNotEqual(
            self.runner().stage_recipe_fingerprint("selfcal"),
            before,
        )


if __name__ == "__main__":
    unittest.main()


class QuarantinePolicyIsNotRecipeTests(CalibrationRecipeAcrossFlaggingRefreshTests):
    """Raising the antenna-quarantine limit must not orphan a calibration."""

    def with_task(self, **settings):
        changed = dict(PIPELINE.CONFIG.CALIBRATION_TASK_CONFIG)
        changed.update(settings)
        return mock.patch.object(
            PIPELINE.CONFIG, "CALIBRATION_TASK_CONFIG", changed
        )

    def test_a_failed_attempt_resumes_after_the_limit_is_raised(self):
        # The observer's case: recorded with the shipped limit of 1, then
        # C03 and S03 both dead, and the limit raised to 2.
        recorded = self.runner().stage_recipe_fingerprint("calibration")

        with self.with_task(max_auto_quarantine_delay_antennas=2):
            runner = self.runner()
            self.assertNotEqual(
                runner.computed_stage_recipe_fingerprint("calibration"),
                recorded,
            )
            runner.resume_state = self.failed_attempt_state(recorded)

            self.assertFalse(runner.should_skip_stage("calibration"))

            # The recorded value is pinned, so the calibration checkpoints
            # bound to it (the flux model) stay valid.
            self.assertEqual(
                runner.stage_recipe_fingerprint("calibration"), recorded
            )
            self.assertEqual(
                runner.recipe_failure_policy("calibration")[
                    "max_auto_quarantine_delay_antennas"
                ],
                1,
            )

    def test_the_same_holds_across_a_flagging_refresh(self):
        recorded = self.runner().stage_recipe_fingerprint("calibration")
        refresh_flagging_provenance(self.outdir, self.ms, "b" * 64)

        with self.with_task(max_auto_quarantine_delay_antennas=3,
                            max_auto_quarantine_delay_fraction=0.15):
            runner = self.runner()
            runner.resume_state = self.failed_attempt_state(recorded)
            self.assertFalse(runner.should_skip_stage("calibration"))

    def test_a_recorded_non_default_policy_is_reproduced(self):
        with self.with_task(max_auto_quarantine_delay_antennas=2):
            recorded = self.runner().stage_recipe_fingerprint("calibration")
        state = self.failed_attempt_state(recorded)
        state["stages"]["calibration"]["attempt_failure_policy"] = dict(
            PIPELINE.CALIBRATION_FAILURE_POLICY_DEFAULTS,
            max_auto_quarantine_delay_antennas=2,
        )

        with self.with_task(max_auto_quarantine_delay_antennas=3):
            runner = self.runner()
            runner.resume_state = state
            self.assertFalse(runner.should_skip_stage("calibration"))

    def test_a_science_setting_changed_with_the_limit_still_fails_closed(self):
        recorded = self.runner().stage_recipe_fingerprint("calibration")

        with self.with_task(max_auto_quarantine_delay_antennas=2,
                            gain_minsnr=7.0):
            runner = self.runner()
            runner.resume_state = self.failed_attempt_state(recorded)
            with self.assertRaisesRegex(
                RuntimeError, "science-recipe fingerprint changed"
            ):
                runner.should_skip_stage("calibration")

    def test_new_attempts_record_the_policy_behind_their_fingerprint(self):
        runner = self.runner()
        runner.resume_state = {"schema_version": 1, "stages": {}}
        runner.save_resume_state = lambda: None
        with self.with_task(max_auto_quarantine_delay_antennas=2):
            runner._stage_recipe_fingerprints.clear()
            runner.update_stage_state("calibration", "running")
        entry = runner.resume_state["stages"]["calibration"]
        self.assertEqual(
            entry["attempt_failure_policy"]["max_auto_quarantine_delay_antennas"],
            2,
        )
