import importlib.util
import json
import os
import tempfile
import unittest
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_pipeline():
    path = SOURCE_DIR / "pipeline.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_pipeline_dependencies",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PIPELINE = load_pipeline()


class PipelineDependencySafetyTests(unittest.TestCase):

    def test_raw_chain_has_fail_closed_dependencies(self):
        stages = [
            "inspection",
            "flagging",
            "calibration",
            "selfcal",
            "imaging",
            "source_detection",
        ]

        self.assertEqual(
            PIPELINE.PipelineRunner.selected_stage_dependencies(
                "calibration", stages
            ),
            {"inspection", "flagging"},
        )
        self.assertEqual(
            PIPELINE.PipelineRunner.selected_stage_dependencies(
                "imaging", stages
            ),
            {"selfcal"},
        )
        self.assertEqual(
            PIPELINE.PipelineRunner.selected_stage_dependencies(
                "source_detection", stages
            ),
            {"imaging"},
        )

    def test_multiple_targets_fail_immediately_before_downstream_stages(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").touch()
            outdir = root / "Products"
            inspection = outdir / "inspection"
            inspection.mkdir(parents=True)
            (inspection / "inspection_report.json").write_text(json.dumps({
                "fields": {
                    "targets": ["target_a", "target_b"],
                },
            }))
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--target-selection", "explicit",
                "--dry-run",
                "--no-copy-input-ms",
            ])
            runner = PIPELINE.PipelineRunner(args)

            with self.assertRaisesRegex(
                RuntimeError,
                "target_a, target_b.*--target TARGET_NAME",
            ):
                runner.validate_downstream_target_selection([
                    "inspection", "flagging", "calibration", "selfcal"
                ])

            runner.args.target = "target_b"
            runner.validate_downstream_target_selection([
                "inspection", "flagging", "calibration", "selfcal"
            ])

            runner.args.target = None
            runner.validate_downstream_target_selection([
                "inspection", "flagging", "calibration"
            ])

    def test_multiple_targets_are_ranked_by_usable_postflag_data(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").touch()
            outdir = root / "Products"
            inspection = outdir / "inspection"
            flagging = outdir / "flagging"
            inspection.mkdir(parents=True)
            flagging.mkdir(parents=True)
            (inspection / "inspection_report.json").write_text(json.dumps({
                "fields": {
                    "targets": ["target_a", "target_b", "target_c"],
                },
            }))
            (flagging / "field_statistics.json").write_text(json.dumps({
                "by_field": {
                    "target_a": {"total": 100, "flagged": 20},
                    "target_b": {"total": 140, "flagged": 25},
                    "target_c": {"total": 100, "flagged": 90},
                },
            }))
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                # Ranked selection of a single target is now opt-in; the
                # default images every calibrated target.
                "--target-selection", "auto",
                "--no-copy-input-ms",
            ])
            runner = PIPELINE.PipelineRunner(args)

            selected = runner.validate_downstream_target_selection([
                "inspection", "flagging", "calibration", "selfcal"
            ])

            self.assertEqual(selected, "target_b")
            self.assertEqual(runner.args.target, "target_b")
            report = json.loads(
                (outdir / "target_selection.json").read_text()
            )
            self.assertEqual(report["selected_target"], "target_b")
            self.assertEqual(
                [item["target"] for item in report["ranking"]],
                ["target_b", "target_a", "target_c"],
            )
            self.assertFalse(report["ranking"][-1]["eligible"])
            self.assertEqual(
                report["policy"],
                "max_usable_postflag_visibilities",
            )

    def test_auto_target_selection_is_deterministic_and_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").touch()
            outdir = root / "Products"
            inspection = outdir / "inspection"
            flagging = outdir / "flagging"
            inspection.mkdir(parents=True)
            flagging.mkdir(parents=True)
            (inspection / "inspection_report.json").write_text(json.dumps({
                "fields": {"targets": ["target_b", "target_a"]},
            }))
            statistics_path = flagging / "field_statistics.json"
            statistics_path.write_text(json.dumps({
                "by_field": {
                    "target_a": {"total": 100, "flagged": 10},
                    "target_b": {"total": 100, "flagged": 10},
                },
            }))
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--target-selection", "auto",
                "--no-copy-input-ms",
            ])
            runner = PIPELINE.PipelineRunner(args)

            runner.validate_downstream_target_selection(["selfcal"])
            self.assertEqual(runner.args.target, "target_a")

            runner.args.target = None
            runner._target_selection_origin = None
            statistics_path.write_text(json.dumps({
                "by_field": {
                    "target_a": {"total": 100, "flagged": 10},
                },
            }))
            with self.assertRaisesRegex(
                RuntimeError,
                "target_b: missing summary",
            ):
                runner.validate_downstream_target_selection(["selfcal"])

    def test_explicit_downstream_only_run_has_no_phantom_raw_dependency(self):
        stages = ["imaging", "source_detection"]

        self.assertEqual(
            PIPELINE.PipelineRunner.selected_stage_dependencies(
                "imaging", stages
            ),
            set(),
        )
        self.assertEqual(
            PIPELINE.PipelineRunner.selected_stage_dependencies(
                "source_detection", stages
            ),
            {"imaging"},
        )

    def test_source_ms_identity_detects_same_path_metadata_replacement(self):
        with tempfile.TemporaryDirectory() as temporary:
            ms = Path(temporary) / "input.ms"
            ms.mkdir()
            table_dat = ms / "table.dat"
            table_dat.write_bytes(b"original")

            before = PIPELINE.PipelineRunner.measurement_set_identity(ms)
            original_stat = table_dat.stat()
            table_dat.write_bytes(b"replaced")
            os.utime(
                table_dat,
                ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns),
            )
            after = PIPELINE.PipelineRunner.measurement_set_identity(ms)

            self.assertNotEqual(before["stat_sha256"], after["stat_sha256"])

    def test_completed_calibration_can_repair_legacy_working_provenance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "GWB.ms"
            source.mkdir()
            (source / "table.dat").write_bytes(b"immutable source")
            outdir = root / "Products"
            working = outdir / "working_ms" / "GWB.ms"
            working.mkdir(parents=True)
            (working / "table.dat").write_bytes(b"calibrated working MS")
            args = PIPELINE.build_parser().parse_args([
                str(source),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "flagging,calibration",
            ])
            runner = PIPELINE.PipelineRunner(args)
            source_identity = runner.measurement_set_identity(source)
            (outdir / "pipeline_run_manifest.json").write_text(json.dumps({
                "schema_version": 3,
                "created_utc": PIPELINE.utc_now(),
                "input_kind": "ms",
                "source_ms": str(source.resolve()),
                "source_ms_identity": source_identity,
                "working_ms": str(working.resolve()),
            }))
            inspection = outdir / "inspection"
            inspection.mkdir()
            (inspection / "inspection_report.json").write_text(json.dumps({
                "measurement_set": str(working.resolve()),
            }))
            calibration = outdir / "calibration"
            calibration.mkdir()
            (calibration / "calibration.done").touch()
            (calibration / "calibration_report.json").write_text(json.dumps({
                "measurement_set_used": str(working.resolve()),
                "measurement_set_stat_sha256": (
                    PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(working)
                ),
            }))
            legacy_provenance = runner.working_ms_manifest_path(working)
            legacy_provenance.write_text(json.dumps({
                "source_ms": str(source.resolve()),
                "working_ms": str(working.resolve()),
                "source_identity": {"schema_version": 1, "legacy": True},
                "copy_method": "legacy-copy",
            }))

            runner.prepare_working_copy(["flagging", "calibration"])

            provenance = json.loads(
                legacy_provenance.read_text()
            )
            self.assertEqual(provenance["source_identity"], source_identity)
            self.assertEqual(
                provenance["adoption"]["method"],
                "completed-calibration-evidence-v1",
            )
            archive = Path(
                provenance["adoption"]["previous_provenance_archive"]
            )
            self.assertTrue(archive.is_file())
            self.assertIn("legacy", archive.name)
            self.assertEqual(runner.args.ms, str(working.resolve()))

    def test_legacy_working_provenance_repair_rejects_changed_ms(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "GWB.ms"
            source.mkdir()
            (source / "table.dat").write_bytes(b"immutable source")
            outdir = root / "Products"
            working = outdir / "working_ms" / "GWB.ms"
            working.mkdir(parents=True)
            (working / "table.dat").write_bytes(b"calibrated working MS")
            args = PIPELINE.build_parser().parse_args([
                str(source),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "flagging,calibration",
            ])
            runner = PIPELINE.PipelineRunner(args)
            source_identity = runner.measurement_set_identity(source)
            (outdir / "pipeline_run_manifest.json").write_text(json.dumps({
                "schema_version": 3,
                "input_kind": "ms",
                "source_ms": str(source.resolve()),
                "source_ms_identity": source_identity,
                "working_ms": str(working.resolve()),
            }))
            inspection = outdir / "inspection"
            inspection.mkdir()
            (inspection / "inspection_report.json").write_text(json.dumps({
                "measurement_set": str(working.resolve()),
            }))
            calibration = outdir / "calibration"
            calibration.mkdir()
            (calibration / "calibration.done").touch()
            recorded_identity = (
                PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(working)
            )
            (calibration / "calibration_report.json").write_text(json.dumps({
                "measurement_set_used": str(working.resolve()),
                "measurement_set_stat_sha256": recorded_identity,
            }))
            (working / "unexpected-replacement").write_bytes(b"changed")

            with self.assertRaisesRegex(
                RuntimeError,
                "live working-MS contents do not match",
            ):
                runner.prepare_working_copy(["flagging", "calibration"])

    def test_source_ms_identity_ignores_transient_casa_locks(self):
        with tempfile.TemporaryDirectory() as temporary:
            ms = Path(temporary) / "input.ms"
            ms.mkdir()
            (ms / "table.dat").write_bytes(b"stable metadata")

            before = PIPELINE.PipelineRunner.measurement_set_identity(ms)
            (ms / "table.lock").write_bytes(b"transient reader lock")
            after = PIPELINE.PipelineRunner.measurement_set_identity(ms)

            self.assertEqual(before, after)
            self.assertEqual(before["schema_version"], 2)

    def test_ms_instance_identity_survives_main_edits_but_rejects_swap(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "working.ms"
            ms.mkdir()
            (ms / "table.dat").write_bytes(b"main before flagging")
            for name in [
                "ANTENNA",
                "DATA_DESCRIPTION",
                "OBSERVATION",
                "POLARIZATION",
                "SPECTRAL_WINDOW",
            ]:
                subtable = ms / name
                subtable.mkdir()
                (subtable / "table.dat").write_bytes(
                    ("stable-" + name).encode("ascii")
                )

            before = PIPELINE.PipelineRunner.measurement_set_instance_identity(
                ms
            )
            (ms / "table.dat").write_bytes(b"main after calibration")
            (ms / "table.f0").write_bytes(b"new corrected-data column")
            after_science_write = (
                PIPELINE.PipelineRunner.measurement_set_instance_identity(ms)
            )
            self.assertEqual(before, after_science_write)

            replacement = root / "replacement.ms"
            replacement.mkdir()
            (replacement / "table.dat").write_bytes(b"replacement main")
            for name in [
                "ANTENNA",
                "DATA_DESCRIPTION",
                "OBSERVATION",
                "POLARIZATION",
                "SPECTRAL_WINDOW",
            ]:
                subtable = replacement / name
                subtable.mkdir()
                (subtable / "table.dat").write_bytes(
                    ("stable-" + name).encode("ascii")
                )
            old = root / "old.ms"
            ms.rename(old)
            replacement.rename(ms)

            after_swap = (
                PIPELINE.PipelineRunner.measurement_set_instance_identity(ms)
            )
            self.assertNotEqual(before, after_swap)

    def test_flagging_fingerprint_binds_managed_ms_instance_not_main_writes(self):
        def make_structured_ms(path, main_value):
            path.mkdir(parents=True)
            (path / "table.dat").write_bytes(main_value)
            for name in [
                "ANTENNA",
                "DATA_DESCRIPTION",
                "OBSERVATION",
                "POLARIZATION",
                "SPECTRAL_WINDOW",
            ]:
                subtable = path / name
                subtable.mkdir()
                (subtable / "table.dat").write_bytes(name.encode("ascii"))

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.ms"
            working = root / "Products" / "working.ms"
            make_structured_ms(source, b"immutable source")
            make_structured_ms(working, b"pre-flag working data")
            args = PIPELINE.build_parser().parse_args([
                str(source),
                "--outdir", str(root / "Products"),
                "--telescope", "ugmrt",
                "--receiver-band", "B3",
                "--stages", "flagging",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)
            runner.source_ms = source.resolve()
            runner.working_ms = working.resolve()
            runner.args.ms = str(working.resolve())
            before = runner.stage_fingerprint("flagging")

            (working / "table.dat").write_bytes(b"post-flag working data")
            (working / "table.f0").write_bytes(b"new FLAG column")
            runner._stage_fingerprints.clear()
            after_flagging = runner.stage_fingerprint("flagging")
            self.assertEqual(before, after_flagging)

            replacement = root / "replacement.ms"
            make_structured_ms(replacement, b"replacement")
            working.rename(root / "old-working.ms")
            replacement.rename(working)
            runner._stage_fingerprints.clear()
            after_swap = runner.stage_fingerprint("flagging")
            self.assertNotEqual(before, after_swap)

    def test_flagging_resume_is_bound_to_input_and_recipe(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").touch()
            outdir = root / "Products"
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B3",
                "--stages", "flagging",
                "--no-copy-input-ms",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)
            flagdir = outdir / "flagging"
            flagdir.mkdir(parents=True)
            (flagdir / "flagging.done").touch()
            report = {
                "telescope": "ugmrt",
                "band": "B3",
                "measurement_set": str(ms.resolve()),
                "provenance_fingerprint": runner.stage_fingerprint("flagging"),
            }
            runner.report_path("flagging_post").write_text(json.dumps(report))

            self.assertTrue(runner.validate_stage_completion("flagging")[0])

            report["provenance_fingerprint"] = "0" * 64
            runner.report_path("flagging_post").write_text(json.dumps(report))
            valid, reason = runner.validate_stage_completion("flagging")
            self.assertFalse(valid)
            self.assertIn("provenance", reason)

    def test_mutable_raw_stages_reject_in_place_measurement_set(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").touch()
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(root / "Products"),
                "--telescope", "ugmrt",
                "--receiver-band", "B3",
                "--stages", "flagging",
                "--no-copy-input-ms",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)

            with self.assertRaisesRegex(RuntimeError, "managed working copy"):
                runner.preflight(["flagging"])

    def test_nonmutating_calibrated_lofar_imaging_allows_no_copy(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").touch()
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(root / "Products"),
                "--telescope", "lofar",
                "--receiver-band", "HBA",
                "--data-state", "calibrated",
                "--imager", "wsclean",
                "--lofar-beam-mode", "full",
                "--stages", "imaging",
                "--no-copy-input-ms",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)

            runner.preflight(["imaging"])

    def test_lofar_beam_none_rejects_pbcor_catalogue_default(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").touch()
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(root / "Products"),
                "--telescope", "lofar",
                "--receiver-band", "HBA",
                "--data-state", "calibrated",
                "--lofar-beam-mode", "none",
                "--stages", "imaging,source_detection",
                "--no-copy-input-ms",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)

            with self.assertRaisesRegex(RuntimeError, "quick-look image"):
                runner.preflight(["imaging", "source_detection"])

    def test_calibration_fingerprint_binds_external_model_contents(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").touch()
            model = root / "calibrator.model"
            model.write_bytes(b"first absolute flux model")
            argv = [
                str(ms),
                "--outdir", str(root / "Products"),
                "--telescope", "meerkat",
                "--receiver-band", "L",
                "--stages", "calibration",
                "--flux-model", str(model),
                "--dry-run",
                "--no-copy-input-ms",
            ]
            first = PIPELINE.PipelineRunner(
                PIPELINE.build_parser().parse_args(argv)
            )
            first_fingerprint = first.stage_fingerprint("calibration")
            first_identity = first.run_identity()["calibrator_model"]

            model.write_bytes(b"replacement absolute flux model")
            second = PIPELINE.PipelineRunner(
                PIPELINE.build_parser().parse_args(argv)
            )

            self.assertNotEqual(
                first_fingerprint,
                second.stage_fingerprint("calibration"),
            )
            self.assertNotEqual(
                first_identity["flux_model_sha256"],
                second.run_identity()["calibrator_model"]["flux_model_sha256"],
            )

    def test_calibration_recipe_fingerprint_ignores_code_only_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").write_bytes(b"stable")
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(root / "Products"),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "calibration",
                "--flux-calibrator", "3C147",
                "--bandpass-calibrator", "3C147",
                "--phase-calibrator", "0323+055",
                "--dry-run",
                "--no-copy-input-ms",
            ])
            baseline = PIPELINE.PipelineRunner(args)
            full_before = baseline.stage_fingerprint("calibration")
            recipe_before = baseline.stage_recipe_fingerprint("calibration")
            script = baseline.root / "calibration.py"
            original_read_bytes = Path.read_bytes

            def changed_bytes(path):
                value = original_read_bytes(path)
                return (
                    value + b"\n# simulated implementation-only fix\n"
                    if path == script else value
                )

            changed = PIPELINE.PipelineRunner(args)
            with mock.patch.object(Path, "read_bytes", changed_bytes):
                full_after = changed.stage_fingerprint("calibration")
                recipe_after = changed.stage_recipe_fingerprint("calibration")

            self.assertNotEqual(full_before, full_after)
            self.assertEqual(recipe_before, recipe_after)

    def test_imaging_recipe_fingerprint_ignores_code_only_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            final_ms = root / "final_selfcal.ms"
            final_ms.mkdir()
            (final_ms / "table.dat").write_bytes(b"stable")
            args = PIPELINE.build_parser().parse_args([
                "--outdir", str(root / "Products"),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "imaging",
                "--imaging-ms", str(final_ms),
                "--target", "science_target",
                "--dry-run",
            ])
            baseline = PIPELINE.PipelineRunner(args)
            full_before = baseline.stage_fingerprint("imaging")
            recipe_before = baseline.stage_recipe_fingerprint("imaging")
            script = baseline.root / "imaging.py"
            original_read_bytes = Path.read_bytes

            def changed_bytes(path):
                value = original_read_bytes(path)
                return (
                    value + b"\n# simulated PB implementation fix\n"
                    if path == script else value
                )

            changed = PIPELINE.PipelineRunner(args)
            with mock.patch.object(Path, "read_bytes", changed_bytes):
                full_after = changed.stage_fingerprint("imaging")
                recipe_after = changed.stage_recipe_fingerprint("imaging")

            self.assertNotEqual(full_before, full_after)
            self.assertEqual(recipe_before, recipe_after)

    def test_legacy_failed_calibration_resumes_when_command_recipe_unchanged(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").write_bytes(b"stable")
            outdir = root / "Products"
            inspection = outdir / "inspection"
            flagging = outdir / "flagging"
            inspection.mkdir(parents=True)
            flagging.mkdir()
            (inspection / "inspection_report.json").write_text(json.dumps({
                "measurement_set": str(ms.resolve()),
                "telescope": "ugmrt",
                "band": "B4",
            }))
            (flagging / "post_flagging_report.json").write_text(json.dumps({
                "measurement_set": str(ms.resolve()),
                "telescope": "ugmrt",
                "band": "B4",
            }))
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "calibration",
                "--flux-calibrator", "3C147",
                "--bandpass-calibrator", "3C147",
                "--secondary-flux-calibrator", "3C48",
                "--phase-calibrator", "0323+055",
                "--dry-run",
                "--no-copy-input-ms",
            ])
            runner = PIPELINE.PipelineRunner(args)
            current_command = runner.command_for_stage("calibration")
            legacy_command = []
            skip = False
            for part in current_command:
                if skip:
                    skip = False
                    continue
                if part == "--recipe-fingerprint":
                    skip = True
                    continue
                legacy_command.append(part)
            provenance_index = legacy_command.index(
                "--provenance-fingerprint"
            ) + 1
            legacy_command[provenance_index] = "0" * 64
            started = datetime.now(timezone.utc) + timedelta(minutes=5)
            runner.resume_state = {
                "schema_version": 1,
                "created_utc": PIPELINE.utc_now(),
                "stages": {
                    "calibration": {
                        "status": "failed",
                        "started_utc": started.isoformat(),
                        "attempt_stage_fingerprint": "0" * 64,
                        "command": legacy_command,
                    }
                },
            }

            self.assertFalse(runner.should_skip_stage("calibration"))
            entry = runner.resume_state["stages"]["calibration"]
            self.assertEqual(
                entry["attempt_recipe_fingerprint"],
                runner.stage_recipe_fingerprint("calibration"),
            )
            self.assertEqual(len(entry["implementation_resume_history"]), 1)

    def test_legacy_failed_imaging_resumes_after_code_only_pb_fix(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            final_ms = root / "final_selfcal.ms"
            final_ms.mkdir()
            (final_ms / "table.dat").write_bytes(b"stable")
            outdir = root / "Products"
            report_path = outdir / "selfcal" / "reports"
            report_path.mkdir(parents=True)
            final_identity = PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(
                final_ms
            )
            (report_path / "selfcal_report.json").write_text(json.dumps({
                "target": "science_target",
                "final_ms": str(final_ms.resolve()),
                "final_ms_stat_sha256": final_identity,
            }))
            args = PIPELINE.build_parser().parse_args([
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "imaging",
                "--imaging-ms", str(final_ms),
                "--target", "science_target",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)
            current_command = runner.command_for_stage("imaging")
            legacy_command = list(current_command)
            provenance_index = legacy_command.index(
                "--provenance-fingerprint"
            ) + 1
            legacy_command[provenance_index] = "0" * 64
            started = datetime.now(timezone.utc) + timedelta(minutes=5)
            runner.resume_state = {
                "schema_version": 1,
                "created_utc": PIPELINE.utc_now(),
                "stages": {
                    "imaging": {
                        "status": "failed",
                        "started_utc": started.isoformat(),
                        "attempt_stage_fingerprint": "0" * 64,
                        "attempt_recipe_fingerprint": "0" * 64,
                        "command": legacy_command,
                    }
                },
            }

            self.assertFalse(runner.should_skip_stage("imaging"))
            entry = runner.resume_state["stages"]["imaging"]
            self.assertEqual(
                entry["attempt_recipe_fingerprint"],
                runner.stage_recipe_fingerprint("imaging"),
            )
            self.assertEqual(len(entry["implementation_resume_history"]), 1)

    def test_failed_imaging_resumes_when_its_recipe_matches(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            final_ms = root / "final_selfcal.ms"
            final_ms.mkdir()
            (final_ms / "table.dat").write_bytes(b"stable")
            outdir = root / "Products"
            report_path = outdir / "selfcal" / "reports"
            report_path.mkdir(parents=True)
            (report_path / "selfcal_report.json").write_text(json.dumps({
                "target": "science_target",
                "final_ms": str(final_ms.resolve()),
            }))
            args = PIPELINE.build_parser().parse_args([
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "imaging",
                "--imaging-ms", str(final_ms),
                "--target", "science_target",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)
            recipe = runner.stage_recipe_fingerprint("imaging")
            runner.resume_state = {
                "schema_version": 1,
                "created_utc": PIPELINE.utc_now(),
                "stages": {
                    "imaging": {
                        "status": "failed",
                        # A newer configuration file would fail the legacy
                        # audit; the recorded recipe settles it instead.
                        "started_utc": "2000-01-01T00:00:00+00:00",
                        "attempt_stage_fingerprint": "0" * 64,
                        "attempt_recipe_fingerprint": recipe,
                    }
                },
            }

            self.assertFalse(runner.should_skip_stage("imaging"))
            entry = runner.resume_state["stages"]["imaging"]
            self.assertIn(
                "science-recipe fingerprint matches",
                entry["implementation_resume_history"][-1]["reason"],
            )

    def completed_imaging_with_another_recipe(self, root, *telescope):
        final_ms = root / "final_selfcal.ms"
        final_ms.mkdir()
        (final_ms / "table.dat").write_bytes(b"stable")
        outdir = root / "Products"
        report_path = outdir / "selfcal" / "reports"
        report_path.mkdir(parents=True)
        (report_path / "selfcal_report.json").write_text(json.dumps({
            "target": "science_target",
            "final_ms": str(final_ms.resolve()),
        }))
        args = PIPELINE.build_parser().parse_args([
            "--outdir", str(outdir), *telescope,
            "--stages", "imaging",
            "--imaging-ms", str(final_ms),
            "--target", "science_target",
            "--dry-run",
        ])
        runner = PIPELINE.PipelineRunner(args)
        runner.resume_state = {
            "schema_version": 1,
            "created_utc": PIPELINE.utc_now(),
            "stages": {
                "imaging": {
                    "status": "completed",
                    "product_stage_fingerprint": "0" * 64,
                    "product_recipe_fingerprint": "1" * 64,
                }
            },
        }
        return runner

    def test_a_wsclean_image_with_a_new_recipe_reruns_in_place(self):
        # Finding 165: a MeerKAT beam-correction fix changes the imaging
        # recipe. The WSClean stage's own checkpoint is bound to WSClean's
        # command and data, so the stage reruns where it is.
        with tempfile.TemporaryDirectory() as temporary:
            runner = self.completed_imaging_with_another_recipe(
                Path(temporary), "--telescope", "meerkat",
                "--receiver-band", "L")
            self.assertEqual(runner.imaging_backend(), "wsclean")
            with mock.patch("builtins.print"):
                self.assertFalse(runner.should_skip_stage("imaging"))
            entry = runner.resume_state["stages"]["imaging"]
            self.assertEqual(entry["status"], "stale")
            self.assertIn("recipe", entry["validation_error"])

    def test_a_casa_image_with_a_new_recipe_still_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            runner = self.completed_imaging_with_another_recipe(
                Path(temporary), "--telescope", "ugmrt",
                "--receiver-band", "B4")
            self.assertEqual(runner.imaging_backend(), "casa")
            with self.assertRaisesRegex(RuntimeError, "Refusing to rerun"):
                runner.should_skip_stage("imaging")

    def test_failed_calibration_recipe_change_still_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").write_bytes(b"stable")
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(root / "Products"),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "calibration",
                "--dry-run",
                "--no-copy-input-ms",
            ])
            runner = PIPELINE.PipelineRunner(args)
            runner.resume_state = {
                "schema_version": 1,
                "created_utc": PIPELINE.utc_now(),
                "stages": {
                    "calibration": {
                        "status": "failed",
                        "attempt_stage_fingerprint": "0" * 64,
                        "attempt_recipe_fingerprint": "1" * 64,
                    }
                },
            }

            with self.assertRaisesRegex(
                RuntimeError,
                "science-recipe fingerprint changed",
            ):
                runner.should_skip_stage("calibration")

    def test_flagging_fingerprint_binds_selector_helper(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").touch()
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(root / "Products"),
                "--telescope", "meerkat",
                "--receiver-band", "L",
                "--stages", "flagging",
                "--dry-run",
                "--no-copy-input-ms",
            ])
            runner = PIPELINE.PipelineRunner(args)
            original = runner.stage_fingerprint("flagging")
            helper = runner.root / "flagging_selectors.py"
            original_read_bytes = Path.read_bytes

            def changed_bytes(path):
                value = original_read_bytes(path)
                if path == helper:
                    return value + b"\n# simulated selector change\n"
                return value

            fresh = PIPELINE.PipelineRunner(args)
            with mock.patch.object(Path, "read_bytes", changed_bytes):
                changed = fresh.stage_fingerprint("flagging")

            self.assertNotEqual(original, changed)

    def test_inputless_downstream_resume_keeps_manifest_source_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            outdir = Path(temporary) / "Products"
            outdir.mkdir()
            args = PIPELINE.build_parser().parse_args([
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--from-stage", "imaging",
                "--resume",
            ])
            runner = PIPELINE.PipelineRunner(args)
            manifest = runner.run_identity()
            manifest.update({
                "input_data": "/archive/original/GWB.ms",
                "input_kind": "ms",
                "source_ms": "/archive/original/GWB.ms",
                "source_ms_identity": {
                    "schema_version": 2,
                    "stat_sha256": "recorded-immutable-source",
                },
                "working_ms": str(
                    outdir / "working_ms" / "GWB.ms"
                ),
            })
            (outdir / "pipeline_run_manifest.json").write_text(
                json.dumps(manifest)
            )

            runner.validate_outdir_support(["imaging"])

            preserved = json.loads(
                (outdir / "pipeline_run_manifest.json").read_text()
            )
            self.assertEqual(
                preserved["source_ms_identity"],
                manifest["source_ms_identity"],
            )

    def test_downstream_resume_ignores_ms_supplied_only_by_config(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            configured_ms = root / "configured-but-unused.ms"
            configured_ms.mkdir()
            (configured_ms / "table.dat").write_bytes(b"configured source")
            outdir = root / "Products"
            outdir.mkdir()
            with mock.patch.dict(
                PIPELINE.RUN_CONFIG,
                {"ms": str(configured_ms)},
            ):
                args = PIPELINE.build_parser().parse_args([
                    "--outdir", str(outdir),
                    "--telescope", "ugmrt",
                    "--receiver-band", "B4",
                    "--from-stage", "imaging",
                    "--resume",
                ])
                runner = PIPELINE.PipelineRunner(args)

            manifest = runner.run_identity()
            manifest["source_ms_identity"] = {
                "schema_version": 2,
                "stat_sha256": "identity-from-original-run",
            }
            (outdir / "pipeline_run_manifest.json").write_text(
                json.dumps(manifest)
            )

            stages = runner.selected_stages()
            runner.discard_unused_configured_primary_input(stages)
            self.assertIsNone(runner.args.ms)
            self.assertIsNone(runner.input_data)
            self.assertIsNone(runner.input_kind)
            runner.validate_outdir_support(stages)

    def test_inputless_resume_restores_managed_ms_for_any_raw_stage(self):
        for resume_stage in ["inspection", "flagging", "calibration"]:
            with self.subTest(stage=resume_stage):
                with tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    outdir = root / "Products"
                    working = outdir / "working_ms" / "GWB.ms"
                    working.mkdir(parents=True)
                    (working / "table.dat").write_bytes(b"managed working")
                    source = root / "offline-original" / "GWB.ms"
                    recorded_identity = {
                        "schema_version": 2,
                        "root_device": 10,
                        "root_inode": 20,
                        "entry_count": 1,
                        "stat_sha256": "recorded-source-identity",
                    }
                    args = PIPELINE.build_parser().parse_args([
                        "--outdir", str(outdir),
                        "--telescope", "ugmrt",
                        "--receiver-band", "B4",
                        "--resume-from", resume_stage,
                        "--resume",
                    ])
                    runner = PIPELINE.PipelineRunner(args)
                    manifest = runner.run_identity()
                    manifest.update({
                        "input_data": str(source.resolve()),
                        "input_kind": "ms",
                        "source_ms": str(source.resolve()),
                        "source_ms_identity": recorded_identity,
                        "working_ms": str(working.resolve()),
                    })
                    (outdir / "pipeline_run_manifest.json").write_text(
                        json.dumps(manifest)
                    )
                    provenance = runner.working_ms_manifest_path(working)
                    provenance.write_text(json.dumps({
                        "schema_version": 2,
                        "source_ms": str(source.resolve()),
                        "working_ms": str(working.resolve()),
                        "source_identity": recorded_identity,
                        "copy_method": "recorded managed copy",
                    }))

                    stages = runner.selected_stages()
                    runner.discard_unused_configured_primary_input(stages)
                    runner.rehydrate_resume_run_context(stages)
                    runner.prepare_working_copy(stages)

                    self.assertEqual(Path(runner.args.ms), working.resolve())
                    self.assertEqual(runner.source_ms, source.resolve())
                    self.assertEqual(
                        runner.source_ms_identity_for_provenance(),
                        recorded_identity,
                    )
                    self.assertTrue(runner._prepared_working_from_manifest)

    def test_explicit_resume_from_alias_matches_from_stage(self):
        parser = PIPELINE.build_parser()
        self.assertEqual(
            parser.parse_args(["--resume-from", "selfcal"]).from_stage,
            "selfcal",
        )
        self.assertEqual(
            parser.parse_args(["--from-stage", "selfcal"]).from_stage,
            "selfcal",
        )

    def test_explicit_source_still_requires_exact_manifest_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "GWB.ms"
            source.mkdir()
            (source / "table.dat").write_bytes(b"current source")
            outdir = root / "Products"
            outdir.mkdir()
            args = PIPELINE.build_parser().parse_args([
                str(source),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "imaging",
                "--no-copy-input-ms",
            ])
            runner = PIPELINE.PipelineRunner(args)
            manifest = runner.run_identity()
            manifest["source_ms_identity"] = {
                "schema_version": 2,
                "stat_sha256": "different-source",
            }
            (outdir / "pipeline_run_manifest.json").write_text(
                json.dumps(manifest)
            )

            with self.assertRaisesRegex(
                RuntimeError,
                "explicitly supplied source Measurement Set",
            ):
                runner.validate_outdir_support(["imaging"])

    def test_no_input_resume_status_restores_managed_ms_context_first(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.ms"
            working = root / "Products" / "working_ms" / "source.ms"
            for ms in [source, working]:
                ms.mkdir(parents=True)
                (ms / "table.dat").write_bytes(b"stable")
            outdir = root / "Products"
            manifest = {
                "schema_version": 3,
                "input_data": str(source),
                "input_kind": "ms",
                "source_ms": str(source),
                "source_ms_identity": (
                    PIPELINE.PipelineRunner.measurement_set_identity(source)
                ),
                "working_ms": str(working),
                "instrument_request": {
                    "telescope": "ugmrt",
                    "receiver_band": "B3",
                    "data_state": "raw",
                    "imager": "auto",
                    "lofar_beam_mode": None,
                },
                "execution_profile": "science",
                "field_selection": {},
                "downstream_selection": {},
                "calibrator_model": {},
            }
            (outdir / "pipeline_run_manifest.json").write_text(
                json.dumps(manifest)
            )
            (outdir / "pipeline_state.json").write_text(
                json.dumps({"stages": {"inspection": {}, "flagging": {}}})
            )
            args = PIPELINE.build_parser().parse_args([
                "--outdir", str(outdir), "--resume-status"
            ])
            runner = PIPELINE.PipelineRunner(args)
            captured = []
            runner.print_resume_status = lambda stages: captured.extend(stages)

            runner.run()

            self.assertEqual(runner.input_data, source.resolve())
            self.assertEqual(runner.input_kind, "ms")
            self.assertEqual(runner.source_ms, source.resolve())
            self.assertEqual(runner.working_ms, working.resolve())
            self.assertEqual(Path(runner.args.ms), working.resolve())
            self.assertEqual(captured, ["inspection", "flagging"])

    def test_no_input_fits_status_selects_ingest_after_rehydration(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "observation.uvfits"
            source.write_bytes(b"visibility fits placeholder")
            outdir = root / "Products"
            imported = outdir / "ingest" / "observation.ms"
            imported.mkdir(parents=True)
            (imported / "table.dat").touch()
            manifest = {
                "schema_version": 3,
                "input_data": str(source),
                "input_kind": "fits",
                "fits_format": "uvfits",
                "imported_ms": str(imported),
                "gmrt_flagfiles": [],
                "instrument_request": {
                    "telescope": "ugmrt",
                    "receiver_band": "B3",
                    "data_state": "raw",
                    "imager": "auto",
                    "lofar_beam_mode": None,
                },
                "execution_profile": "science",
                "field_selection": {},
                "downstream_selection": {},
                "calibrator_model": {},
            }
            (outdir / "pipeline_run_manifest.json").write_text(
                json.dumps(manifest)
            )
            (outdir / "pipeline_state.json").write_text(
                json.dumps({"stages": {"ingest": {}, "inspection": {}}})
            )
            args = PIPELINE.build_parser().parse_args([
                "--outdir", str(outdir), "--resume-status"
            ])
            runner = PIPELINE.PipelineRunner(args)
            captured = []
            runner.print_resume_status = lambda stages: captured.extend(stages)

            runner.run()

            self.assertEqual(runner.input_kind, "fits")
            self.assertEqual(runner.input_data, source.resolve())
            self.assertEqual(runner.imported_ms, imported.resolve())
            self.assertEqual(captured, ["ingest", "inspection"])

    def test_resume_restores_omitted_science_context_but_not_cli_override(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            outdir = root / "Products"
            working = outdir / "working_ms" / "GWB.ms"
            working.mkdir(parents=True)
            (working / "table.dat").write_bytes(b"managed working")
            source = root / "offline" / "GWB.ms"
            args = PIPELINE.build_parser().parse_args([
                "--outdir", str(outdir),
                "--resume-from", "calibration",
                "--resume",
                "--bandpass-calibrator", "CLI_BANDPASS",
            ])
            runner = PIPELINE.PipelineRunner(args)
            manifest = runner.run_identity()
            manifest.update({
                "input_data": str(source.resolve()),
                "input_kind": "ms",
                "source_ms": str(source.resolve()),
                "source_ms_identity": {
                    "schema_version": 2,
                    "stat_sha256": "recorded-source",
                },
                "working_ms": str(working.resolve()),
                "instrument_request": {
                    "telescope": "ugmrt",
                    "receiver_band": "B4",
                    "data_state": "raw",
                    "imager": "casa",
                    "lofar_beam_mode": None,
                },
                "field_selection": {
                    "flux_calibrator": "3C147,3C48",
                    "bandpass_calibrator": "3C147",
                    "secondary_flux_calibrator": "3C48",
                    "phase_calibrator": "0323+055",
                    "target_phase_map": None,
                    "target_fields": None,
                    "reference_antenna": "C05",
                },
            })
            (outdir / "pipeline_run_manifest.json").write_text(
                json.dumps(manifest)
            )

            runner.rehydrate_resume_run_context(runner.selected_stages())

            self.assertEqual(runner.args.receiver_band, "B4")
            self.assertEqual(runner.args.flux_calibrator, "3C147,3C48")
            self.assertEqual(runner.args.phase_calibrator, "0323+055")
            self.assertEqual(runner.args.reference_antenna, "C05")
            self.assertEqual(
                runner.args.bandpass_calibrator,
                "CLI_BANDPASS",
            )
            self.assertEqual(Path(runner.args.ms), working.resolve())

    def test_inputless_resume_status_allows_offline_archival_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            outdir = root / "Products"
            working = outdir / "working_ms" / "source.ms"
            working.mkdir(parents=True)
            (working / "table.dat").write_bytes(b"managed working")
            offline_source = root / "unmounted" / "source.ms"
            manifest = {
                "schema_version": 3,
                "input_data": str(offline_source),
                "input_kind": "ms",
                "source_ms": str(offline_source),
                "source_ms_identity": {
                    "schema_version": 2,
                    "stat_sha256": "recorded-source",
                },
                "working_ms": str(working),
                "instrument_request": {
                    "telescope": "ugmrt",
                    "receiver_band": "B4",
                    "data_state": "raw",
                    "imager": "casa",
                    "lofar_beam_mode": None,
                },
                "execution_profile": "science",
                "field_selection": {},
                "downstream_selection": {},
                "calibrator_model": {},
            }
            (outdir / "pipeline_run_manifest.json").write_text(
                json.dumps(manifest)
            )
            (outdir / "pipeline_state.json").write_text(
                json.dumps({"stages": {"inspection": {}}})
            )
            args = PIPELINE.build_parser().parse_args([
                "--outdir", str(outdir), "--resume-status"
            ])
            runner = PIPELINE.PipelineRunner(args)
            captured = []
            runner.print_resume_status = lambda stages: captured.extend(stages)

            runner.run()

            self.assertEqual(runner.input_data, offline_source.resolve())
            self.assertEqual(Path(runner.args.ms), working.resolve())
            self.assertEqual(captured, ["inspection"])

    @staticmethod
    def _write_valid_inspection_products(outdir, ms):
        inspection = outdir / "inspection"
        inspection.mkdir(parents=True)
        (inspection / "inspection_report.json").write_text(json.dumps({
            "measurement_set": str(ms.resolve()),
            "telescope": "ugmrt",
            "band": "B4",
        }))
        for filename in [
            "fields.json",
            "antennas.json",
            "qa_plot_status.json",
        ]:
            (inspection / filename).write_text("{}")

    def test_completed_inspection_is_adopted_after_code_only_upgrade(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").write_bytes(b"stable")
            outdir = root / "Products"
            self._write_valid_inspection_products(outdir, ms)
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "inspection",
                "--resume",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)
            runner.resume_state = {
                "schema_version": 1,
                "created_utc": PIPELINE.utc_now(),
                "stages": {
                    "inspection": {
                        "status": "completed",
                        "product_stage_fingerprint": "0" * 64,
                        "product_recipe_fingerprint": "0" * 64,
                        "command": runner.command_for_stage("inspection"),
                    },
                },
            }

            self.assertTrue(runner.should_skip_stage("inspection"))
            entry = runner.resume_state["stages"]["inspection"]
            self.assertEqual(
                entry["product_stage_fingerprint"],
                runner.stage_fingerprint("inspection"),
            )
            self.assertEqual(len(entry["implementation_adoption_history"]), 1)

    def test_stale_inspection_reports_are_overwritten_after_code_upgrade(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").write_bytes(b"stable")
            outdir = root / "Products"
            self._write_valid_inspection_products(outdir, ms)
            (outdir / "inspection" / "fields.json").unlink()
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "inspection",
                "--resume",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)
            runner.resume_state = {
                "schema_version": 1,
                "created_utc": PIPELINE.utc_now(),
                "stages": {
                    "inspection": {
                        "status": "completed",
                        "product_stage_fingerprint": "0" * 64,
                        "command": runner.command_for_stage("inspection"),
                    },
                },
            }

            self.assertFalse(runner.should_skip_stage("inspection"))
            self.assertEqual(
                runner.resume_state["stages"]["inspection"]["status"],
                "stale",
            )

    def test_completed_stage_upgrade_does_not_hide_changed_science_command(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").write_bytes(b"stable")
            outdir = root / "Products"
            self._write_valid_inspection_products(outdir, ms)
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "inspection",
                "--resume",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)
            old_command = runner.command_for_stage("inspection") + [
                "--changed-science-option"
            ]
            runner.resume_state = {
                "schema_version": 1,
                "created_utc": PIPELINE.utc_now(),
                "stages": {
                    "inspection": {
                        "status": "completed",
                        "product_stage_fingerprint": "0" * 64,
                        "command": old_command,
                    },
                },
            }

            with self.assertRaisesRegex(
                RuntimeError,
                "recorded science command changed",
            ):
                runner.should_skip_stage("inspection")

    def test_completed_flagging_is_adopted_across_code_only_upgrade(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ms = root / "input.ms"
            ms.mkdir()
            (ms / "table.dat").write_bytes(b"stable")
            outdir = root / "Products"
            flagging = outdir / "flagging"
            flagging.mkdir(parents=True)
            (flagging / "flagging.done").touch()
            (flagging / "postflag_summary.json").write_text(json.dumps({
                "measurement_set": str(ms.resolve()),
                "telescope": "ugmrt",
                "band": "B4",
                "provenance_fingerprint": "0" * 64,
            }))
            args = PIPELINE.build_parser().parse_args([
                str(ms),
                "--outdir", str(outdir),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "flagging",
                "--resume",
                "--dry-run",
                "--no-copy-input-ms",
            ])
            runner = PIPELINE.PipelineRunner(args)
            runner.resume_state = {
                "schema_version": 1,
                "created_utc": PIPELINE.utc_now(),
                "stages": {
                    "flagging": {
                        "status": "completed",
                        "product_stage_fingerprint": "0" * 64,
                        "command": runner.command_for_stage("flagging"),
                    },
                },
            }

            self.assertTrue(runner.should_skip_stage("flagging"))
            entry = runner.resume_state["stages"]["flagging"]
            self.assertEqual(
                entry["product_stage_fingerprint"],
                runner.stage_fingerprint("flagging"),
            )
            self.assertEqual(len(entry["implementation_adoption_history"]), 1)


if __name__ == "__main__":
    unittest.main()
