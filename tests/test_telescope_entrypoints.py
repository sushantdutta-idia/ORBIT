import contextlib
import importlib.util
import io
import sys
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]

if str(SOURCE_DIR) not in sys.path:
    sys.path.insert(0, str(SOURCE_DIR))

import telescope_entrypoint as ENTRYPOINT


def load_launcher(name):
    path = SOURCE_DIR / f"run_{name}.py"
    spec = importlib.util.spec_from_file_location(
        f"test_{name}_launcher",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TelescopeEntrypointTests(unittest.TestCase):

    def test_each_launcher_selects_only_its_named_telescope(self):
        for telescope in ("ugmrt", "meerkat", "lofar"):
            with self.subTest(telescope=telescope):
                launcher = load_launcher(telescope)

                with mock.patch.object(
                    launcher,
                    "run_telescope_cli",
                    return_value=0,
                ) as delegated:
                    self.assertEqual(launcher.main(["example.ms"]), 0)

                delegated.assert_called_once_with(
                    telescope,
                    ["example.ms"],
                )

    def test_shared_config_defaults_remain_available(self):
        args = ENTRYPOINT.parse_for_telescope("meerkat", [])

        self.assertEqual(args.telescope, "meerkat")
        self.assertEqual(
            args.outdir,
            ENTRYPOINT.PIPELINE.RUN_CONFIG["outdir"],
        )
        self.assertEqual(
            args.casa,
            ENTRYPOINT.PIPELINE.RUN_CONFIG["casa"],
        )

    def test_cli_values_still_override_shared_config(self):
        args = ENTRYPOINT.parse_for_telescope(
            "meerkat",
            [
                "observation.ms",
                "--receiver-band",
                "L",
                "--outdir",
                "MeerKAT_products",
            ],
        )

        self.assertEqual(args.ms, "observation.ms")
        self.assertEqual(args.receiver_band, "L")
        self.assertEqual(args.outdir, "MeerKAT_products")

    def test_lofar_launcher_defaults_to_calibrated_wsclean_path(self):
        args = ENTRYPOINT.parse_for_telescope(
            "lofar",
            ["target.ms", "--lofar-beam-mode", "full"],
        )

        self.assertEqual(args.telescope, "lofar")
        self.assertEqual(args.data_state, "calibrated")
        self.assertEqual(args.imager, "wsclean")
        self.assertEqual(args.lofar_beam_mode, "full")

    def test_explicit_lofar_values_are_not_silently_rewritten(self):
        args = ENTRYPOINT.parse_for_telescope(
            "lofar",
            [
                "target.ms",
                "--data-state=raw",
                "--imager",
                "auto",
            ],
        )

        self.assertEqual(args.data_state, "raw")
        self.assertEqual(args.imager, "auto")

    def test_the_ugmrt_launcher_is_pinned_to_tclean(self):
        args = ENTRYPOINT.parse_for_telescope("ugmrt", ["target.ms"])
        self.assertEqual(args.imager, "casa")

    def test_a_shared_config_cannot_send_ugmrt_to_wsclean(self):
        # The same config.yaml often serves several launchers; an
        # `imager: wsclean` meant for MeerKAT must not reach a uGMRT run.
        with mock.patch.dict(ENTRYPOINT.PIPELINE.RUN_CONFIG,
                             {"imager": "wsclean"}):
            ugmrt = ENTRYPOINT.parse_for_telescope("ugmrt", ["target.ms"])
            meerkat = ENTRYPOINT.parse_for_telescope("meerkat", ["target.ms"])
        self.assertEqual(ugmrt.imager, "casa")
        # MeerKAT's launcher leaves the configured choice alone.
        self.assertEqual(meerkat.imager, "wsclean")

    def test_asking_ugmrt_for_wsclean_on_the_command_line_still_fails(self):
        args = ENTRYPOINT.parse_for_telescope(
            "ugmrt", ["target.ms", "--imager", "wsclean",
                      "--receiver-band", "B4", "--stages", "imaging"])
        self.assertEqual(args.imager, "wsclean")
        with self.assertRaisesRegex(ValueError, "CASA backend"):
            ENTRYPOINT.PIPELINE.PipelineRunner(args).preflight(["imaging"])

    def test_meerkat_launcher_defaults_to_the_profile_choice(self):
        args = ENTRYPOINT.parse_for_telescope("meerkat", ["target.ms"])
        self.assertEqual(args.imager, "auto")
        runner = ENTRYPOINT.PIPELINE.PipelineRunner(args)
        self.assertEqual(runner.imaging_backend(), "wsclean")

    def test_launcher_rejects_a_conflicting_telescope(self):
        stderr = io.StringIO()

        with contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as raised:
                ENTRYPOINT.parse_for_telescope(
                    "meerkat",
                    ["target.ms", "--telescope", "lofar"],
                )

        self.assertEqual(raised.exception.code, 2)
        self.assertIn("fixed to --telescope meerkat", stderr.getvalue())

    def test_launcher_delegates_execution_to_pipeline_runner(self):
        with mock.patch.object(
            ENTRYPOINT.PIPELINE,
            "PipelineRunner",
        ) as runner_class:
            result = ENTRYPOINT.run_telescope_cli(
                "ugmrt",
                ["target.ms", "--dry-run"],
            )

        self.assertEqual(result, 0)
        runner_class.assert_called_once()
        runner_class.return_value.run.assert_called_once_with()
        parsed = runner_class.call_args.args[0]
        self.assertEqual(parsed.telescope, "ugmrt")
        self.assertTrue(parsed.dry_run)


if __name__ == "__main__":
    unittest.main()
