import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_pipeline_module():
    path = SOURCE_DIR / "pipeline.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_pipeline_cli",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PIPELINE = load_pipeline_module()


class PipelineTelescopeCLITests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.outdir = self.root / "Products"
        self.ms = self.root / "input.ms"

    def tearDown(self):
        self.temporary.cleanup()

    def parse(self, *extra):
        argv = [
            str(self.ms),
            "--outdir",
            str(self.outdir),
            "--dry-run",
            "--mpi",
            "off",
            "--no-copy-input-ms",
            *extra,
        ]
        return PIPELINE.build_parser().parse_args(argv)

    def write_inspection(self, telescope, band):
        path = self.outdir / "inspection" / "inspection_report.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(
                {
                    "measurement_set": str(self.ms),
                    "telescope": {"name": telescope},
                    "band": {"name": band},
                }
            )
        )

    def test_parser_accepts_supported_telescopes_and_rejects_unknown(self):
        for telescope in ["auto", "ugmrt", "meerkat", "lofar"]:
            with self.subTest(telescope=telescope):
                args = self.parse("--telescope", telescope)
                self.assertEqual(args.telescope, telescope)

        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as raised:
                self.parse("--telescope", "vla")
        self.assertEqual(raised.exception.code, 2)

    def test_inspection_command_forwards_instrument_context(self):
        args = self.parse(
            "--telescope",
            "meerkat",
            "--receiver-band",
            "L",
            "--data-state",
            "raw",
            "--stages",
            "inspection",
        )
        command = PIPELINE.PipelineRunner(args).command_for_stage(
            "inspection"
        )

        self.assertIn(str(SOURCE_DIR / "inspect.py"), command)
        self.assertEqual(command[command.index("--telescope") + 1], "meerkat")
        self.assertEqual(command[command.index("--receiver-band") + 1], "L")
        self.assertEqual(command[command.index("--data-state") + 1], "raw")

    def test_calibration_command_forwards_model_and_field_overrides(self):
        args = self.parse(
            "--telescope",
            "meerkat",
            "--stages",
            "calibration",
            "--flux-calibrator",
            "J1939-6342,J0408-6545",
            "--bandpass-calibrator",
            "J1939-6342",
            "--phase-calibrator",
            "J1234-5678,J2345-6789",
            "--target-phase-map",
            "target_a:J1234-5678,target_b:J2345-6789",
            "--target-fields",
            "target_a,target_b",
            "--flux-standard",
            "Perley-Butler 2017",
            "--flux-model",
            "/models/meerkat.im",
        )
        command = PIPELINE.PipelineRunner(args).command_for_stage(
            "calibration"
        )

        expected = {
            "--telescope": "meerkat",
            "--flux-calibrator": "J1939-6342,J0408-6545",
            "--bandpass-calibrator": "J1939-6342",
            "--phase-calibrator": "J1234-5678,J2345-6789",
            "--target-phase-map": (
                "target_a:J1234-5678,target_b:J2345-6789"
            ),
            "--targets": "target_a,target_b",
            "--flux-standard": "Perley-Butler 2017",
            "--flux-model": "/models/meerkat.im",
        }
        for option, value in expected.items():
            with self.subTest(option=option):
                self.assertEqual(command[command.index(option) + 1], value)

    def test_downstream_target_is_not_turned_into_calibration_targets(self):
        args = self.parse(
            "--telescope",
            "meerkat",
            "--stages",
            "calibration",
            "--target",
            "downstream_product_target",
        )
        command = PIPELINE.PipelineRunner(args).command_for_stage(
            "calibration"
        )

        self.assertNotIn("--targets", command)
        self.assertNotIn("downstream_product_target", command)

        # The same option remains available to genuinely downstream stages.
        source_command = PIPELINE.PipelineRunner(args).command_for_stage(
            "source_detection"
        )
        self.assertEqual(
            source_command[source_command.index("--target") + 1],
            "downstream_product_target",
        )

    def test_lofar_auto_imager_routes_to_wsclean(self):
        calibrated_ms = self.root / "lofar-calibrated.ms"
        args = self.parse(
            "--telescope",
            "lofar",
            "--receiver-band",
            "HBA",
            "--data-state",
            "calibrated",
            "--stages",
            "imaging",
            "--imaging-ms",
            str(calibrated_ms),
            "--target",
            "field_0",
            "--imager",
            "auto",
            "--lofar-beam-mode",
            "differential",
        )
        runner = PIPELINE.PipelineRunner(args)
        command = runner.command_for_stage("imaging")

        self.assertEqual(runner.imaging_backend(), "wsclean")
        self.assertEqual(command[0], sys.executable)
        self.assertEqual(command[1], str(SOURCE_DIR / "wsclean_imaging.py"))
        self.assertEqual(command[command.index("--telescope") + 1], "lofar")
        self.assertEqual(command[command.index("--receiver-band") + 1], "HBA")
        self.assertEqual(
            command[command.index("--beam-mode") + 1],
            "differential",
        )
        self.assertEqual(
            command[command.index("--datacolumn") + 1],
            "CORRECTED_DATA",
        )

    def test_auto_telescope_uses_inspection_report_for_backend_and_band(self):
        self.write_inspection("lofar", "LOFAR_HBA")
        calibrated_ms = self.root / "lofar-calibrated.ms"
        args = self.parse(
            "--telescope",
            "auto",
            "--data-state",
            "calibrated",
            "--stages",
            "imaging",
            "--imaging-ms",
            str(calibrated_ms),
            "--lofar-beam-mode",
            "full",
        )
        runner = PIPELINE.PipelineRunner(args)

        self.assertEqual(runner.resolved_telescope(), "lofar")
        self.assertEqual(runner.resolved_receiver_band(), "LOFAR_HBA")
        self.assertEqual(runner.imaging_backend(), "wsclean")
        self.assertIn(
            str(SOURCE_DIR / "wsclean_imaging.py"),
            runner.command_for_stage("imaging"),
        )

    def test_lofar_imaging_requires_an_explicit_beam_mode(self):
        calibrated_ms = self.root / "lofar-calibrated.ms"
        args = self.parse(
            "--telescope",
            "lofar",
            "--receiver-band",
            "HBA",
            "--data-state",
            "calibrated",
            "--stages",
            "imaging",
            "--imaging-ms",
            str(calibrated_ms),
            "--imager",
            "auto",
        )

        self.assertIsNone(args.lofar_beam_mode)
        with self.assertRaisesRegex(RuntimeError, "explicit beam state"):
            PIPELINE.PipelineRunner(args).command_for_stage("imaging")

    def meerkat_args(self, *extra):
        return self.parse(
            "--telescope",
            "meerkat",
            "--receiver-band",
            "L",
            "--stages",
            "imaging",
            "--imaging-ms",
            str(self.root / "meerkat.ms"),
            *extra,
        )

    def test_meerkat_auto_imager_routes_to_wsclean(self):
        args = self.meerkat_args("--imager", "auto")
        runner = PIPELINE.PipelineRunner(args)
        command = runner.command_for_stage("imaging")

        self.assertEqual(runner.imaging_backend(), "wsclean")
        self.assertEqual(command[0], sys.executable)
        self.assertEqual(command[1], str(SOURCE_DIR / "wsclean_imaging.py"))
        self.assertEqual(command[command.index("--telescope") + 1], "meerkat")
        # MeerKAT has no LOFAR beam state to state, and none is passed.
        self.assertNotIn("--beam-mode", command)
        # Its data column is chosen like the CASA backend chooses it.
        self.assertEqual(command[command.index("--datacolumn") + 1], "auto")

    def test_meerkat_does_not_need_a_lofar_beam_mode(self):
        args = self.meerkat_args()
        self.assertIsNone(args.lofar_beam_mode)
        PIPELINE.PipelineRunner(args).command_for_stage("imaging")

    def test_meerkat_ignores_a_stray_lofar_beam_mode(self):
        args = self.meerkat_args("--lofar-beam-mode", "full")
        command = PIPELINE.PipelineRunner(args).command_for_stage("imaging")
        self.assertNotIn("--beam-mode", command)

    def test_an_explicit_data_column_is_honoured_for_meerkat(self):
        args = self.meerkat_args("--wsclean-data-column", "DATA")
        command = PIPELINE.PipelineRunner(args).command_for_stage("imaging")
        self.assertEqual(command[command.index("--datacolumn") + 1], "DATA")

    def test_meerkat_can_still_be_imaged_with_casa(self):
        args = self.meerkat_args("--imager", "casa")
        runner = PIPELINE.PipelineRunner(args)
        command = runner.command_for_stage("imaging")

        self.assertEqual(runner.imaging_backend(), "casa")
        self.assertEqual(command[0], args.casa)
        self.assertIn(str(SOURCE_DIR / "imaging.py"), command)
        self.assertIn("--no-subbands", command)

    def test_ugmrt_stays_on_tclean_whatever_the_imager_setting(self):
        for imager in ("auto", "casa"):
            with self.subTest(imager=imager):
                args = self.parse(
                    "--telescope", "ugmrt", "--receiver-band", "B4",
                    "--stages", "imaging",
                    "--imaging-ms", str(self.root / "ugmrt.ms"),
                    "--imager", imager,
                )
                runner = PIPELINE.PipelineRunner(args)
                command = runner.command_for_stage("imaging")
                self.assertEqual(runner.imaging_backend(), "casa")
                self.assertIn(str(SOURCE_DIR / "imaging.py"), command)
                self.assertNotIn(
                    str(SOURCE_DIR / "wsclean_imaging.py"), command)

        # Asking for WSClean on uGMRT is refused before anything runs.
        args = self.parse(
            "--telescope", "ugmrt", "--receiver-band", "B4",
            "--stages", "imaging",
            "--imaging-ms", str(self.root / "ugmrt.ms"),
            "--imager", "wsclean",
        )
        with self.assertRaisesRegex(ValueError, "CASA backend"):
            PIPELINE.PipelineRunner(args).preflight(["imaging"])

    def ugmrt_imaging_command(self, *extra):
        imaging_ms = self.root / "ugmrt.ms"
        imaging_ms.mkdir(exist_ok=True)
        (imaging_ms / "table.dat").write_text("table")
        args = self.parse(
            "--telescope",
            "ugmrt",
            "--receiver-band",
            "B3",
            "--stages",
            "imaging",
            "--imaging-ms",
            str(imaging_ms),
            *extra,
        )
        return PIPELINE.PipelineRunner(args).command_for_stage("imaging")

    def test_subband_imaging_is_off_unless_requested(self):
        default = self.ugmrt_imaging_command()
        self.assertIn("--no-subbands", default)
        self.assertNotIn("--make-subbands", default)
        self.assertFalse(
            PIPELINE.IMAGING_TASK_CONFIG.get("make_subbands", False)
        )

        count_only = self.ugmrt_imaging_command("--subbands", "4")
        self.assertIn("--no-subbands", count_only)

        requested = self.ugmrt_imaging_command(
            "--make-subbands",
            "--subbands",
            "4",
        )
        self.assertIn("--make-subbands", requested)
        self.assertNotIn("--no-subbands", requested)

    def test_preflight_rejects_unsupported_telescope_stage_combinations(self):
        cases = [
            (
                [
                    "--telescope",
                    "lofar",
                    "--data-state",
                    "raw",
                    "--stages",
                    "calibration",
                ],
                ["calibration"],
                "Raw LOFAR calibration",
            ),
            (
                [
                    "--telescope",
                    "lofar",
                    "--data-state",
                    "calibrated",
                    "--stages",
                    "selfcal",
                ],
                ["selfcal"],
                "self-calibration stage",
            ),
            (
                [
                    "--telescope",
                    "meerkat",
                    "--data-state",
                    "raw",
                    "--stages",
                    "ddcal",
                    "--include-ddcal",
                ],
                ["ddcal"],
                "validated only for uGMRT",
            ),
        ]

        for extra, stages, message in cases:
            with self.subTest(stages=stages, message=message):
                runner = PIPELINE.PipelineRunner(self.parse(*extra))
                with self.assertRaisesRegex(ValueError, message):
                    runner.preflight(stages)

    def test_run_identity_records_canonical_instrument_context(self):
        args = self.parse(
            "--telescope",
            "meerkat",
            "--receiver-band",
            "L",
            "--data-state",
            "raw",
            "--imager",
            "auto",
            "--stages",
            "inspection",
            "--threads",
            "1",
        )
        identity = PIPELINE.PipelineRunner(args).run_identity()

        self.assertEqual(identity["schema_version"], 3)
        self.assertEqual(
            identity["instrument_request"]["telescope"],
            "meerkat",
        )
        self.assertEqual(
            identity["instrument_context"],
            {
                "telescope": "meerkat",
                "receiver_band": "MEERKAT_L",
                "imaging_backend": "wsclean",
            },
        )
        self.assertEqual(
            identity["instrument_request"]["profile_schema_version"],
            PIPELINE.TELESCOPE_PROFILE_MODULE.PROFILE_SCHEMA_VERSION,
        )

    def test_generated_sbatch_embeds_telescope_backend_options(self):
        outfile = self.root / "pipeline.sbatch"
        args = self.parse(
            "--telescope",
            "lofar",
            "--receiver-band",
            "HBA",
            "--data-state",
            "calibrated",
            "--imager",
            "wsclean",
            "--wsclean",
            "/opt/wsclean/bin/wsclean",
            "--lofar-beam-mode",
            "full",
            "--stages",
            "imaging",
            "--write-sbatch",
            str(outfile),
        )

        PIPELINE.write_sbatch_from_config(args)
        content = outfile.read_text()

        for expected in [
            "--telescope lofar",
            "--receiver-band HBA",
            "--data-state calibrated",
            "--imager wsclean",
            "--wsclean /opt/wsclean/bin/wsclean",
            "--lofar-beam-mode full",
        ]:
            with self.subTest(expected=expected):
                self.assertIn(expected, content)


if __name__ == "__main__":
    unittest.main()
