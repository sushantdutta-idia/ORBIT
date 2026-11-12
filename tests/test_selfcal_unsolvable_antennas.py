"""Data that self-calibration cannot correct is not imaged with data it has.

On the uGMRT target that prompted this, five antennas kept 6 to 7% of their
data. They got no self-calibration solutions, so their visibilities, about 3%
of the total, went into the image with the calibrator-transfer phases that
phase self-calibration had improved on by 70% for everything else.
"""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_selfcal_resume_safety import load_selfcal_module


SELFCAL = load_selfcal_module()

CENSUS = {
    "status": "measured",
    "antennas": {
        "C00": {"id": 0, "unflagged_samples": 15.4e6, "solvable": True},
        "C01": {"id": 1, "unflagged_samples": 15.3e6, "solvable": True},
        "C02": {"id": 2, "unflagged_samples": 0.97e6, "solvable": False},
        "C03": {"id": 3, "unflagged_samples": 0.0, "solvable": False},
        "W06": {"id": 29, "unflagged_samples": 1.18e6, "solvable": False},
    },
}


class RowTable:

    def open(self, path):
        return None

    def close(self):
        return None

    def colnames(self):
        return ["DATA", "MODEL_DATA", "CORRECTED_DATA"]

    def nrows(self):
        return 100


class UnsolvableAntennaTests(unittest.TestCase):

    def manager(self, census=CENSUS, **task_config):
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.task_config = dict(task_config)
        manager.antenna_data = census
        return manager

    def flag(self, manager):
        calls = []
        with mock.patch.object(
            SELFCAL, "flagdata", lambda **kwargs: calls.append(kwargs)
        ), mock.patch("builtins.print"):
            record = manager.flag_unsolvable_antenna_data("candidate.ms")
        return calls, record

    def test_the_remnant_of_an_unsolvable_antenna_is_flagged(self):
        calls, record = self.flag(self.manager())

        self.assertEqual(len(calls), 1)
        call = calls[0]
        self.assertEqual(call["vis"], "candidate.ms")
        self.assertEqual(call["mode"], "manual")
        # C03 has no data left and needs nothing.
        self.assertEqual(call["antenna"], "C02,W06")
        self.assertEqual(record["antennas"], ["C02", "W06"])
        self.assertAlmostEqual(
            record["fraction_of_antenna_samples"],
            (0.97e6 + 1.18e6) / (15.4e6 + 15.3e6 + 0.97e6 + 1.18e6),
        )

    def test_an_array_of_working_and_dead_antennas_is_left_alone(self):
        census = {
            "status": "measured",
            "antennas": {
                "C00": {"unflagged_samples": 15.4e6, "solvable": True},
                "C03": {"unflagged_samples": 0.0, "solvable": False},
            },
        }
        calls, record = self.flag(self.manager(census))
        self.assertEqual(calls, [])
        self.assertIsNone(record)

    def test_it_can_be_switched_off(self):
        calls, _ = self.flag(self.manager(flag_unsolvable_antennas=False))
        self.assertEqual(calls, [])

    def test_without_a_census_nothing_is_flagged(self):
        calls, _ = self.flag(self.manager({"status": "unavailable"}))
        self.assertEqual(calls, [])

    def test_it_happens_in_the_candidate_before_the_solutions_are_applied(self):
        with tempfile.TemporaryDirectory() as temporary:
            msdir = Path(temporary)
            manager = self.manager(
                applymode="calonly", calwt=False, parang=False
            )
            manager.msdir = msdir
            manager.current_ms = str(msdir / "target.ms")
            manager.checkpoint_exists = lambda name: False
            manager.caltable_status = lambda caltable: (True, "")
            order = []

            def split(**kwargs):
                order.append(("split", kwargs["outputvis"]))
                Path(kwargs["outputvis"]).mkdir()

            with mock.patch.object(SELFCAL, "split", split), \
                    mock.patch.object(
                        SELFCAL, "flagdata",
                        lambda **kw: order.append(("flagdata", kw["vis"])),
                    ), mock.patch.object(
                        SELFCAL, "applycal",
                        lambda **kw: order.append(("applycal", kw["vis"])),
                    ), mock.patch.object(
                        SELFCAL, "tb", RowTable(), create=True
                    ), mock.patch("builtins.print"):
                manager.apply_phase_selfcal(1, "phase_1.cal")

            candidate = str(msdir / "phase1.candidate.ms")
            self.assertEqual(
                order,
                [("split", candidate), ("flagdata", candidate),
                 ("applycal", candidate)],
            )
            # The accepted Measurement Set is never the one flagged.
            self.assertNotIn(("flagdata", manager.msdir / "target.ms"), order)


if __name__ == "__main__":
    unittest.main()
