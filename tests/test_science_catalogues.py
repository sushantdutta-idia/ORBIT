import unittest
import importlib.util
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


class ScienceCatalogueTests(unittest.TestCase):

    def test_solution_application_is_fail_closed_by_default(self):
        path = SOURCE_DIR / "config.py"
        spec = importlib.util.spec_from_file_location(
            "test_radio_continuum_science_config",
            path,
        )
        config = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(config)

        self.assertEqual(
            config.CALIBRATION_TASK_CONFIG["applymode"],
            "calflagstrict",
        )
        # Self-calibration is the exception, and deliberately so: a phase
        # loop whose increment cannot be solved for some data leaves them
        # with the calibration they already carry instead of flagging them.
        self.assertEqual(
            config.SELFCAL_TASK_CONFIG["applymode"],
            "calonly",
        )
        # Total-intensity calibration loses nothing by re-referencing, and
        # 'strict' flagged whole science scans next to a calibrator scan in
        # which the reference antenna was flagged.
        self.assertEqual(
            config.CALIBRATION_TASK_CONFIG["refantmode"],
            "flex",
        )
        self.assertEqual(config.FLAGGING_CONFIG["quack_interval_s"], 0.0)
        self.assertEqual(config.FLAGGING_CONFIG["max_iterations"], 1)
        self.assertFalse(config.FLAGGING_CONFIG["precal_rflag"])
        self.assertEqual(
            config.TELESCOPE_FLAGGING_CONFIG["meerkat"]["quack_interval_s"],
            0.0,
        )

    def test_3c295_aliases_are_correct_and_oq208_is_not_a_flux_standard(self):
        entries = {
            line.strip().upper()
            for line in (SOURCE_DIR / "flux-cals.list").read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }

        self.assertIn("3C295", entries)
        self.assertIn("1411+522", entries)
        self.assertIn("J1411+5212", entries)
        self.assertNotIn("1407+282", entries)
        self.assertNotIn("J1407+2827", entries)


if __name__ == "__main__":
    unittest.main()
