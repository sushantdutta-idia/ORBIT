import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_selfcal_resume_safety import load_selfcal_module


SELFCAL = load_selfcal_module()


class FakeColumnTable:

    def __init__(self, columns):
        self.columns = columns

    def open(self, path):
        return None

    def close(self):
        return None

    def colnames(self):
        return list(self.columns)


class FakeFlagdata:
    """Stand-in recording flagdata calls and reporting flagged fractions."""

    def __init__(self, fractions):
        self.fractions = list(fractions)
        self.calls = []

    def __call__(self, **kwargs):
        self.calls.append(kwargs)
        if kwargs.get("mode") == "summary":
            flagged = self.fractions.pop(0)
            return {"total": 1000.0, "flagged": flagged * 1000.0}
        return None


class ResidualFlaggingTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        self.imagedir = root / "images"
        self.reportdir = root / "reports"
        self.checkdir = root / "checkpoints"
        for directory in (self.imagedir, self.reportdir, self.checkdir):
            directory.mkdir()
        model = self.imagedir / "phase_1.model.tt0"
        model.mkdir()
        (model / "table.dat").touch()

        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.imagedir = self.imagedir
        manager.reportdir = self.reportdir
        manager.checkdir = self.checkdir
        manager.current_ms = str(root / "phase1.ms")
        manager.phase_history = [{"loop": 1}]
        manager.task_config = {
            "residual_flagging": True,
            "residual_flag_timedevscale": 5.0,
            "residual_flag_freqdevscale": 5.0,
            "residual_flag_winsize": 5,
            "residual_flag_max_new_fraction": 0.10,
            "residual_flag_extend": True,
            "residual_flag_growtime": 80.0,
            "residual_flag_growfreq": 80.0,
        }
        manager.ensure_model_column = mock.Mock()
        # A deconvolved field by default; individual tests override this.
        manager.deconvolution_diagnostics = lambda base: {
            "model_nonzero": True
        }
        self.manager = manager

    def tearDown(self):
        self.temporary.cleanup()

    def run_pass(self, fractions, columns=("DATA", "MODEL_DATA")):
        flagdata = FakeFlagdata(fractions)
        flagmanager = mock.Mock()
        table = FakeColumnTable(dict.fromkeys(columns, None))
        with mock.patch.object(SELFCAL, "flagdata", flagdata), \
                mock.patch.object(SELFCAL, "flagmanager", flagmanager), \
                mock.patch.object(SELFCAL, "tb", table, create=True):
            self.manager.flag_residual_rfi()
        return flagdata, flagmanager

    def record(self):
        return json.loads(
            (self.reportdir / "residual_flagging.json").read_text()
        )

    def test_residuals_are_flagged_then_extended_within_each_scan(self):
        flagdata, flagmanager = self.run_pass([0.40, 0.44])

        modes = [call["mode"] for call in flagdata.calls]
        self.assertEqual(modes, ["summary", "rflag", "extend", "summary"])

        rflag = flagdata.calls[1]
        # The model is subtracted, so real emission is not what is detected.
        self.assertEqual(rflag["datacolumn"], "residual_data")
        self.assertFalse(rflag["combinescans"])
        self.assertFalse(rflag["extendflags"])
        self.assertEqual(rflag["timedevscale"], 5.0)

        extend = flagdata.calls[2]
        self.assertFalse(extend["combinescans"])
        self.assertEqual(extend["growtime"], 80.0)

        record = self.record()
        self.assertEqual(record["status"], "completed")
        self.assertAlmostEqual(record["new_fraction"], 0.04)
        self.assertTrue((self.checkdir / "residual_flagging.done").exists())
        self.manager.ensure_model_column.assert_called_once()
        # The flag state is saved before anything is flagged.
        self.assertEqual(flagmanager.call_args_list[0].kwargs["mode"], "save")

    def model(self, loop):
        model = self.imagedir / f"phase_{loop}.model.tt0"
        model.mkdir(exist_ok=True)
        (model / "table.dat").touch()
        return self.imagedir / f"phase_{loop}"

    def test_the_model_subtracted_is_the_last_accepted_loops(self):
        # Three loops accepted and a fourth rejected, as on the uGMRT target
        # that showed this: the data kept are loop 3's, and loop 4's image,
        # made from the discarded candidate, is still on disk.
        for loop in (1, 2, 3, 4):
            self.model(loop)
        self.manager.phase_history = [
            {"loop": 1, "accepted": True},
            {"loop": 2, "accepted": True},
            {"loop": 3, "accepted": True},
            {"loop": 4, "accepted": False},
        ]

        self.assertEqual(
            self.manager.latest_image_base(), self.imagedir / "phase_3"
        )
        self.run_pass([0.40, 0.41])
        self.manager.ensure_model_column.assert_called_once_with(
            self.manager.current_ms, self.imagedir / "phase_3"
        )
        self.assertEqual(
            self.record()["parameters"]["image"],
            str(self.imagedir / "phase_3"),
        )

    def test_with_every_loop_rejected_the_initial_model_is_used(self):
        initial = self.imagedir / "initial.model.tt0"
        initial.mkdir()
        self.manager.phase_history = [{"loop": 1, "accepted": False}]

        self.assertEqual(
            self.manager.latest_image_base(), self.imagedir / "initial"
        )

    def test_corrected_column_is_used_when_present(self):
        flagdata, _ = self.run_pass(
            [0.40, 0.41],
            columns=("DATA", "CORRECTED_DATA", "MODEL_DATA"),
        )
        self.assertEqual(flagdata.calls[1]["datacolumn"], "residual")

    def test_an_over_flagging_pass_is_undone_in_full(self):
        flagdata, flagmanager = self.run_pass([0.40, 0.65])

        restore = flagmanager.call_args_list[-1].kwargs
        self.assertEqual(restore["mode"], "restore")
        self.assertEqual(
            restore["versionname"],
            flagmanager.call_args_list[0].kwargs["versionname"],
        )

        record = self.record()
        self.assertEqual(record["status"], "reverted")
        self.assertIn("25.0%", record["reason"])

    def test_it_is_skipped_when_disabled_or_without_a_model(self):
        self.manager.task_config["residual_flagging"] = False
        flagdata, _ = self.run_pass([])
        self.assertEqual(flagdata.calls, [])

        self.manager.task_config["residual_flagging"] = True
        for product in self.imagedir.iterdir():
            product.rename(product.with_suffix(".unused"))
        flagdata, _ = self.run_pass([])
        self.assertEqual(flagdata.calls, [])

    def test_a_measurement_set_without_a_model_is_left_alone(self):
        flagdata, _ = self.run_pass([], columns=("DATA",))
        self.assertEqual(
            [call["mode"] for call in flagdata.calls],
            [],
        )

    def test_completed_pass_is_not_repeated_on_resume(self):
        (self.checkdir / "residual_flagging.done").touch()
        flagdata, _ = self.run_pass([])
        self.assertEqual(flagdata.calls, [])

    def test_an_empty_clean_model_is_skipped_not_failed(self):
        self.manager.deconvolution_diagnostics = lambda base: {
            "model_nonzero": False,
            "residual_peak_snr": 2.0,
        }

        flagdata, _ = self.run_pass([])

        self.assertEqual(flagdata.calls, [])
        record = self.record()
        self.assertEqual(record["status"], "skipped")
        self.assertIn("empty clean model", record["reason"])

    def test_an_unexpected_failure_is_recorded_not_raised(self):
        self.manager.ensure_model_column = mock.Mock(
            side_effect=RuntimeError("tclean did not produce MODEL_DATA")
        )
        self.manager.deconvolution_diagnostics = lambda base: {
            "model_nonzero": True
        }

        # The stage keeps going: the images do not depend on this pass.
        self.run_pass([0.4])

        record = self.record()
        self.assertEqual(record["status"], "failed")
        self.assertIn("tclean did not produce MODEL_DATA", record["reason"])

    def test_an_empty_model_is_reported_at_the_end_of_the_stage(self):
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.target = "G0231-0533"
        manager.clean_diagnostics = {
            "final": {
                "model_nonzero": False,
                "residual_rms_Jy_beam": 1.0e-4,
                "residual_peak_abs_Jy_beam": 5.0e-4,
                "residual_peak_snr": 5.0,
            }
        }

        self.assertTrue(manager.sky_model_is_empty())
        self.assertTrue(manager.warn_if_no_sky_model())

        manager.clean_diagnostics["final"]["model_nonzero"] = True
        self.assertFalse(manager.sky_model_is_empty())
        self.assertFalse(manager.warn_if_no_sky_model())

    def test_it_runs_before_the_final_image(self):
        order = []
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        for name in (
            "print_summary",
            "initial_imaging",
            "estimate_snr",
            "determine_solints",
            "run_phase_selfcal",
            "run_amp_phase_selfcal",
            "flag_residual_rfi",
            "final_imaging",
            "generate_qa_plots",
            "write_report",
        ):
            setattr(
                manager,
                name,
                (lambda step: lambda *a, **k: order.append(step))(name),
            )
        manager.mark_checkpoint = lambda name: order.append(f"mark:{name}")
        manager.telescope = "ugmrt"

        manager.run()

        self.assertLess(
            order.index("flag_residual_rfi"),
            order.index("final_imaging"),
        )


if __name__ == "__main__":
    unittest.main()
