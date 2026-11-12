"""A dead antenna in the pre-bandpass solve is excluded, not fatal."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_calibration_table_qa import CALIBRATION


class PrebpQuarantineTests(unittest.TestCase):

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
        manager.flagging = {"bad_antennas": []}
        manager.delay_quarantine = {
            "schema_version": 1,
            "measurement_set": str(Path(manager.input_ms).resolve()),
            "antennas": [],
        }
        manager.task_config = {
            "auto_quarantine_failed_delay_antennas": True,
            "max_auto_quarantine_delay_antennas": 1,
            "max_auto_quarantine_delay_fraction": 0.10,
            "min_antennas_after_delay_quarantine": 8,
            "prebp_gaintype": "G",
            "gain_minsnr": 3.0,
            "minblperant": 4,
            "gain_solmode": "",
            "parang": False,
            "delay_minsnr": 3.0,
        }
        manager.band_config = {"combine_delay": "", "delay_solint": "inf"}
        (root / "inspection").mkdir()
        (root / "inspection" / "antennas.json").write_text(json.dumps([
            {"id": value, "name": f"C{value:02d}"} for value in range(30)
        ]))
        manager.active_measurement_set_antenna_ids = lambda: set(range(30))
        manager.get_reference_antenna = lambda: "C00"
        manager.get_bandpass_calibrator = lambda: "3C48"
        manager.caltable_checkpoint_ready = lambda *args: False
        manager.mark_checkpoint = lambda *args: None
        self.manager = manager

    def tearDown(self):
        self.temporary.cleanup()

    def run_prebp(self, failed_ids, prebp_results):
        """prebp_results: successive caltable_status results for prebp.cal."""

        manager = self.manager
        calls = []
        results = list(prebp_results)

        def gaincal(**kwargs):
            calls.append(kwargs)
            table = Path(kwargs["caltable"])
            table.mkdir(exist_ok=True)
            (table / "table.dat").touch()

        def status(filename):
            if filename == "prebp.cal":
                return results.pop(0)
            return True, ""

        manager.caltable_status = status
        manager.require_caltable = lambda filename: None
        manager.failed_delay_antenna_ids = lambda filename: list(failed_ids)
        with mock.patch.object(CALIBRATION, "gaincal", gaincal), \
                mock.patch("builtins.print"):
            manager.run_prebp()
        return calls

    def test_a_dead_antenna_is_quarantined_and_both_solves_redone(self):
        calls = self.run_prebp(
            [3],
            [(False, "ANTENNA1=3: group has no usable solution"), (True, "")],
        )

        tables = [Path(call["caltable"]).name for call in calls]
        # The failed solve, then delay and pre-bandpass without C03.
        self.assertEqual(tables, ["prebp.cal", "delay.cal", "prebp.cal"])
        self.assertEqual(calls[0]["antenna"], "")
        self.assertEqual(calls[1]["antenna"], "!C03")
        self.assertEqual(calls[2]["antenna"], "!C03")

        caldir = self.manager.caldir
        self.assertTrue((caldir / "prebp.cal.incomplete.1").is_dir())
        record = json.loads(self.manager.delay_quarantine_path().read_text())
        item = record["antennas"][0]
        self.assertEqual(item["name"], "C03")
        self.assertEqual(item["reason"], "no usable pre-bandpass phase solution")
        self.assertTrue(item["source_table"].endswith("prebp.cal"))

    def test_more_than_one_dead_antenna_still_stops_the_run(self):
        with self.assertRaises(RuntimeError) as raised:
            self.run_prebp(
                [3, 7],
                [(False, "ANTENNA1=3: group has no usable solution")],
            )
        message = str(raised.exception)
        # It names the antennas and the setting that would let it proceed.
        self.assertIn("C03, C07", message)
        self.assertIn("calibration_max_auto_quarantine_delay_antennas", message)
        self.assertFalse(self.manager.delay_quarantine_path().exists())

    def test_a_failure_that_is_not_one_antenna_is_not_quarantined(self):
        with self.assertRaisesRegex(RuntimeError, "SPECTRAL_WINDOW_ID"):
            self.run_prebp(
                [],
                [(False, "SPECTRAL_WINDOW_ID=0: group has no usable solution")],
            )

    def test_the_array_limit_counts_earlier_quarantines(self):
        # Three already excluded out of 30; a fourth would exceed 10%.
        self.manager.delay_quarantine["antennas"] = [
            {"id": value, "name": f"C{value:02d}"} for value in (5, 6, 7)
        ]
        self.assertEqual(
            self.manager.quarantine_unsolvable_antennas(
                [3], "prebp.cal", "no usable pre-bandpass phase solution"
            ),
            [],
        )

    def test_the_reference_antenna_is_never_quarantined(self):
        self.manager.get_reference_antenna = lambda: "C03"
        with self.assertRaisesRegex(RuntimeError, "reference-antenna"):
            self.manager.quarantine_unsolvable_antennas(
                [3], "prebp.cal", "no usable pre-bandpass phase solution"
            )


if __name__ == "__main__":
    unittest.main()
