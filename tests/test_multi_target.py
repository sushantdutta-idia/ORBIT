import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

from test_selfcal_resume_safety import load_selfcal_module


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_radio_continuum_pipeline_multi_target",
    SOURCE_DIR / "pipeline.py",
)
PIPELINE = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PIPELINE)
SELFCAL = load_selfcal_module()

TARGETS = ("G0231-0533", "3C48A", "NGC 1234")


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


class MultiTargetScopeTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name).resolve()
        self.ms = make_ms(root / "input.ms")
        self.outdir = root / "Products"
        (self.outdir / "inspection").mkdir(parents=True)
        (self.outdir / "calibration").mkdir()
        (self.outdir / "inspection" / "inspection_report.json").write_text(
            json.dumps({
                "telescope": "ugmrt",
                "band": "B4",
                "fields": {"targets": list(TARGETS)},
            })
        )
        split = {}
        for name in TARGETS:
            target_ms = make_ms(
                self.outdir / "calibration" / "targets" / f"{name}.ms"
            )
            split[name] = str(target_ms)
        (self.outdir / "calibration" / "split_targets.json").write_text(
            json.dumps(split)
        )

    def tearDown(self):
        self.temporary.cleanup()

    def runner(self, *extra):
        args = PIPELINE.build_parser().parse_args([
            str(self.ms),
            "--outdir", str(self.outdir),
            "--telescope", "ugmrt",
            "--receiver-band", "B4",
            "--stages", "selfcal,imaging",
            "--dry-run",
            "--no-copy-input-ms",
            *extra,
        ])
        return PIPELINE.PipelineRunner(args)

    def scopes(self, runner):
        return [
            (target, str(root))
            for target, root in runner.science_target_scopes(
                ["selfcal", "imaging"]
            )
        ]

    def test_every_calibrated_target_is_imaged_by_default(self):
        scopes = self.scopes(self.runner())

        self.assertEqual([target for target, _ in scopes], list(TARGETS))
        for target, root in scopes:
            self.assertEqual(
                Path(root),
                self.outdir / "targets" /
                PIPELINE.target_directory_name(target),
            )
        # Directory names are filesystem safe without colliding.
        self.assertEqual(
            PIPELINE.target_directory_name("NGC 1234"),
            "NGC_1234",
        )

    def test_products_of_an_earlier_single_target_run_stay_in_place(self):
        reports = self.outdir / "selfcal" / "reports"
        reports.mkdir(parents=True)
        (reports / "selfcal_report.json").write_text(
            json.dumps({"target": "3C48A"})
        )

        scopes = dict(self.scopes(self.runner()))

        self.assertEqual(scopes["3C48A"], str(self.outdir))
        self.assertEqual(
            scopes["G0231-0533"],
            str(self.outdir / "targets" / "G0231-0533"),
        )

    def test_a_named_target_or_subset_restricts_the_run(self):
        self.assertEqual(
            [target for target, _ in self.scopes(
                self.runner("--target", "3C48A")
            )],
            ["3C48A"],
        )
        self.assertEqual(
            [target for target, _ in self.scopes(
                self.runner("--targets", "NGC 1234,G0231-0533")
            )],
            ["G0231-0533", "NGC 1234"],
        )

    def test_an_unknown_target_fails_closed(self):
        with self.assertRaisesRegex(RuntimeError, "not among"):
            self.scopes(self.runner("--targets", "3C48A,NOT_OBSERVED"))

    def test_other_selection_modes_are_kept(self):
        explicit = self.runner("--target-selection", "explicit")
        with self.assertRaisesRegex(RuntimeError, "explicit"):
            self.scopes(explicit)

        (self.outdir / "flagging").mkdir()
        (self.outdir / "flagging" / "field_statistics.json").write_text(
            json.dumps({"by_field": {
                "G0231-0533": {"total": 100, "flagged": 90},
                "3C48A": {"total": 100, "flagged": 10},
                "NGC 1234": {"total": 100, "flagged": 20},
            }})
        )
        chosen = self.scopes(self.runner("--target-selection", "auto"))
        self.assertEqual([target for target, _ in chosen], ["3C48A"])

    def test_an_explicit_downstream_ms_is_one_dataset(self):
        runner = self.runner("--selfcal-ms", str(self.ms))
        self.assertEqual(self.scopes(runner), [(None, str(self.outdir))])


class TargetScopePathTests(MultiTargetScopeTests):

    def test_stage_commands_products_and_state_are_per_target(self):
        runner = self.runner()
        seen = {}
        for target, root in runner.science_target_scopes(["selfcal"]):
            with runner.target_scope(target, root):
                command = [
                    str(part) for part in runner.command_for_stage("selfcal")
                ]
                seen[target] = command
                self.assertEqual(
                    command[command.index("--outdir") + 1],
                    str(root),
                )
                # Inspection and calibration stay at the run's root.
                self.assertEqual(
                    command[command.index("--upstream-outdir") + 1],
                    str(self.outdir),
                )
                self.assertEqual(
                    command[command.index("--target") + 1],
                    target,
                )
                measurement_set = command[command.index("--outdir") - 1]
                self.assertTrue(
                    measurement_set.endswith(f"{target}.ms"),
                    measurement_set,
                )
                self.assertEqual(
                    runner.report_path("selfcal"),
                    Path(root) / "selfcal" / "reports" / "selfcal_report.json",
                )
                self.assertEqual(
                    runner.report_path("calibration"),
                    self.outdir / "calibration" / "calibration_report.json",
                )
                self.assertEqual(
                    runner.stage_state_key("selfcal"),
                    f"selfcal@{target}",
                )

        self.assertEqual(len(seen), len(TARGETS))
        # Leaving the scope restores the run-level view.
        self.assertEqual(runner.stage_outdir, self.outdir)
        self.assertEqual(runner.stage_state_key("selfcal"), "selfcal")

    def test_source_detection_is_scoped_too(self):
        runner = self.runner()
        target, root = runner.science_target_scopes(["source_detection"])[0]
        with runner.target_scope(target, root):
            command = [
                str(part)
                for part in runner.command_for_stage("source_detection")
            ]
        self.assertEqual(command[command.index("--outdir") + 1], str(root))

    def test_each_target_keeps_its_own_resume_record(self):
        runner = self.runner()
        for target, root in runner.science_target_scopes(["selfcal"]):
            with runner.target_scope(target, root):
                runner.update_stage_state("selfcal", "failed", returncode=1)

        recorded = sorted(runner.load_resume_state()["stages"])
        self.assertEqual(
            recorded,
            sorted(f"selfcal@{target}" for target in TARGETS),
        )


class PerTargetIsolationTests(MultiTargetScopeTests):

    def test_validators_look_inside_the_target_directory(self):
        runner = self.runner()
        target, root = runner.science_target_scopes(["source_detection"])[1]
        with runner.target_scope(target, root):
            valid, reason = runner.validate_stage_completion(
                "source_detection"
            )
        self.assertFalse(valid)
        # The marker is looked for where that target's stage writes it.
        self.assertIn(str(Path(root) / "source_detection"), reason)

    def test_one_failing_target_does_not_stop_the_others(self):
        runner = self.runner()
        attempted = []

        def fake_stage(stage, stages, unavailable):
            attempted.append((runner.target_scope_name, stage))
            if runner.target_scope_name == "3C48A":
                raise RuntimeError("selfcal exploded")

        runner.run_selected_stage = fake_stage

        with self.assertRaisesRegex(RuntimeError, "1 of 3 targets") as raised:
            runner.run_per_target_stages(
                ["selfcal", "imaging"],
                ["selfcal", "imaging"],
                set(),
            )

        self.assertIn("3C48A (selfcal exploded)", str(raised.exception))
        self.assertEqual(
            [target for target, _ in attempted],
            ["G0231-0533", "G0231-0533", "3C48A", "NGC 1234", "NGC 1234"],
        )

    def test_an_interruption_stops_the_whole_run(self):
        runner = self.runner()
        attempted = []

        def fake_stage(stage, stages, unavailable):
            attempted.append(runner.target_scope_name)
            raise PIPELINE.StageInterrupted("cancelled")

        runner.run_selected_stage = fake_stage

        # A failed target is skipped past; a cancelled job is not.
        with self.assertRaises(PIPELINE.StageInterrupted):
            runner.run_per_target_stages(["selfcal"], ["selfcal"], set())

        self.assertEqual(attempted, ["G0231-0533"])

    def test_a_single_target_failure_is_raised_directly(self):
        runner = self.runner("--target", "3C48A")

        def fake_stage(stage, stages, unavailable):
            raise RuntimeError("selfcal exploded")

        runner.run_selected_stage = fake_stage

        with self.assertRaisesRegex(RuntimeError, "selfcal exploded"):
            runner.run_per_target_stages(["selfcal"], ["selfcal"], set())


class StageScriptUpstreamTests(unittest.TestCase):

    def test_selfcal_reads_shared_reports_from_the_run_root(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            (root / "inspection").mkdir()
            (root / "calibration").mkdir()
            (root / "inspection" / "inspection_report.json").write_text(
                json.dumps({"telescope": "ugmrt"})
            )
            (root / "calibration" / "calibration_report.json").write_text(
                json.dumps({"reference_antenna": "C03"})
            )

            manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
            manager.outdir = root / "targets" / "NGC_1234"
            manager.upstream_outdir = root
            manager.load_inspection_report()
            manager.load_calibration_report()

            self.assertEqual(manager.inspection["telescope"], "ugmrt")
            self.assertEqual(
                manager.calibration["reference_antenna"],
                "C03",
            )


if __name__ == "__main__":
    unittest.main()
