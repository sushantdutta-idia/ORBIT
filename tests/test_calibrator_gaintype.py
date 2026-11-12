"""The calibrator gain solves fall back from T to G when the hands disagree.

The case is the Band 5 run of 2026-10-21: the phase solve with one phase for
both hands (T) kept 9-28% of the solutions of C11, S02, S03 and the reference
C00, which had kept all of their per-hand (G) pre-bandpass solutions.
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_calibration_table_qa import CALIBRATION


T_FAILURE = (
    "ANTENNA1=10: flagged fraction 0.914530 exceeds the configured group "
    "maximum 0.800000"
)


class CalibratorGaintypeTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        manager = CALIBRATION.CalibrationManager.__new__(
            CALIBRATION.CalibrationManager
        )
        manager.outdir = root
        manager.caldir = root / "calibration"
        manager.caldir.mkdir()
        manager.input_ms = str(root / "input.ms")
        manager.ms = manager.input_ms
        manager.flagging = {"bad_antennas": ["C03", "C05"]}
        manager.task_config = {
            "phase_gaintype": "T",
            "amp_gaintype": "T",
            "refantmode": "flex",
            "minblperant": 4,
            "gain_solmode": "L1R",
            "gain_minsnr": 3.0,
            "parang": False,
            "target_phase_transfer_solint": "inf",
        }
        manager.band_config = {"combine_delay": "", "phase_solint": "10s"}
        manager.get_phase_calibrators = lambda: ["0323+055"]
        manager.get_flux_calibrators = lambda: ["3C48"]
        manager.get_bandpass_calibrator = lambda: "3C48"
        manager.get_reference_antenna = lambda: "C00"
        manager.caltable_checkpoint_ready = lambda *args: False
        self.checkpoints = []
        manager.mark_checkpoint = self.checkpoints.append
        manager.require_caltable = lambda filename: None
        self.manager = manager

    def tearDown(self):
        self.temporary.cleanup()

    def solve(self, method, results):
        """Run a solve; results are successive table statuses."""

        manager = self.manager
        calls = []
        results = list(results)

        def gaincal(**kwargs):
            calls.append(kwargs)
            table = Path(kwargs["caltable"])
            table.mkdir(exist_ok=True)
            (table / "table.dat").touch()

        manager.caltable_status = lambda filename: results.pop(0)
        with mock.patch.object(CALIBRATION, "gaincal", gaincal), \
                mock.patch("builtins.print"):
            getattr(manager, method)()
        return calls

    def test_a_good_t_solve_is_kept(self):
        calls = self.solve("run_phasecal", [(True, "")])
        self.assertEqual([call["gaintype"] for call in calls], ["T"])
        self.assertEqual(calls[0]["calmode"], "p")
        self.assertEqual(calls[0]["field"], "0323+055,3C48")
        self.assertEqual(calls[0]["antenna"], "!C03,C05")
        self.assertEqual(self.checkpoints, ["phase"])
        self.assertFalse(self.manager.gaintype_switch_path().exists())
        self.assertEqual(self.manager.calibrator_gaintype("amp"), "T")

    def test_incoherent_hands_switch_the_run_to_g(self):
        calls = self.solve("run_phasecal", [(False, T_FAILURE), (True, "")])

        self.assertEqual([call["gaintype"] for call in calls], ["T", "G"])
        self.assertTrue(
            (self.manager.caldir / "phase.cal.incomplete.1").is_dir()
        )
        record = json.loads(self.manager.gaintype_switch_path().read_text())
        self.assertEqual(record["gaintype"], "G")
        self.assertEqual(record["table"], "phase.cal")
        self.assertIn("ANTENNA1=10", record["reason"])
        self.assertEqual(record["configured_gaintype"], {"phase": "T", "amp": "T"})
        self.assertEqual(self.checkpoints, ["phase"])
        # Every later calibrator gain solve of the run follows.
        self.assertEqual(self.manager.calibrator_gaintype("phase"), "G")
        self.assertEqual(self.manager.calibrator_gaintype("amp"), "G")
        transfer = self.solve("run_phase_transfer_cal", [])
        self.assertEqual(transfer[0]["gaintype"], "G")

    def test_when_g_fails_too_both_reasons_are_given(self):
        with self.assertRaises(RuntimeError) as raised:
            self.solve(
                "run_phasecal",
                [(False, T_FAILURE), (False, "ANTENNA1=10: group has no usable solution")],
            )
        message = str(raised.exception)
        self.assertIn("T: ANTENNA1=10: flagged fraction", message)
        self.assertIn("G: ANTENNA1=10: group has no usable solution", message)
        self.assertFalse(self.manager.gaintype_switch_path().exists())
        self.assertEqual(self.checkpoints, [])

    def test_a_resumed_run_solves_with_g_at_once(self):
        self.manager.record_gaintype_switch = (
            CALIBRATION.CalibrationManager.record_gaintype_switch.__get__(
                self.manager
            )
        )
        with mock.patch("builtins.print"):
            self.manager.record_gaintype_switch("phase.cal", T_FAILURE)
        calls = self.solve("run_phasecal", [(True, "")])
        self.assertEqual([call["gaintype"] for call in calls], ["G"])

    def test_a_switch_recorded_for_another_measurement_set_is_ignored(self):
        self.manager.gaintype_switch_path().write_text(json.dumps({
            "gaintype": "G",
            "measurement_set": "/elsewhere/other.ms",
        }))
        self.assertEqual(self.manager.calibrator_gaintype("phase"), "T")

    def test_a_g_configuration_is_not_tried_twice(self):
        self.manager.task_config["phase_gaintype"] = "G"

        def require(filename):
            raise RuntimeError(f"{filename} unusable")

        self.manager.require_caltable = require
        with self.assertRaisesRegex(RuntimeError, "phase.cal unusable"):
            self.solve("run_phasecal", [(False, T_FAILURE), (True, "")])
        self.assertFalse(self.manager.gaintype_switch_path().exists())


if __name__ == "__main__":
    unittest.main()
