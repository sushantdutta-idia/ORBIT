"""What a self-calibration stage reported must survive its being resumed.

A resumed stage restores accepted loops from their checkpoints. It images and
flags nothing again for them, so anything the first run only held in memory
is absent from the report the resumed run writes.
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_selfcal_caltable_safety import SELFCAL


def manager(reportdir):
    value = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
    value.reportdir = reportdir
    value.task_config = {}
    return value


def clean_summary(result):
    return {"stop_reason": "reached its threshold", "nmajordone": result}


class CleanPassRecordTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.reportdir = Path(self.temporary.name)
        patches = [
            mock.patch.object(SELFCAL, "tclean", lambda **kwargs: kwargs["n"]),
            mock.patch.object(
                SELFCAL.CASA_COMPAT, "summarize_tclean", clean_summary
            ),
        ]
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)

    def tearDown(self):
        self.temporary.cleanup()

    def saved(self):
        return json.loads((self.reportdir / "clean_passes.json").read_text())

    def test_each_pass_is_written_as_it_finishes(self):
        first = manager(self.reportdir)

        first.tclean_pass({"n": 9}, "phase_1", "bootstrap")
        first.tclean_pass({"n": 11}, "phase_1", "auto-multithresh")

        self.assertEqual(
            [item["name"] for item in self.saved()["phase_1"]],
            ["bootstrap", "auto-multithresh"],
        )
        self.assertEqual(self.saved()["phase_1"][1]["nmajordone"], 11)

    def test_a_resumed_stage_reports_the_passes_of_restored_loops(self):
        first = manager(self.reportdir)
        first.tclean_pass({"n": 9}, "phase_1", "bootstrap")

        resumed = manager(self.reportdir)
        resumed.tclean_pass({"n": 8}, "amp_phase", "bootstrap")

        self.assertEqual(
            sorted(resumed.recorded_clean_passes()), ["amp_phase", "phase_1"]
        )
        self.assertEqual(sorted(self.saved()), ["amp_phase", "phase_1"])

    def test_a_stage_that_images_nothing_still_reports_them(self):
        manager(self.reportdir).tclean_pass({"n": 9}, "phase_1", "bootstrap")

        resumed = manager(self.reportdir)

        self.assertEqual(
            resumed.recorded_clean_passes()["phase_1"][0]["nmajordone"], 9
        )

    def test_an_image_made_again_replaces_its_earlier_account(self):
        first = manager(self.reportdir)
        first.tclean_pass({"n": 9}, "amp_phase", "bootstrap")
        first.tclean_pass({"n": 20}, "amp_phase", "auto-multithresh")

        resumed = manager(self.reportdir)
        with mock.patch.object(
            resumed, "bootstrap_clean", side_effect=RuntimeError("stop")
        ):
            with self.assertRaises(RuntimeError):
                resumed.run_tclean_with_recovery({"n": 1}, "amp_phase")
        resumed.tclean_pass({"n": 7}, "amp_phase", "bootstrap")

        self.assertEqual(
            [item["nmajordone"] for item in self.saved()["amp_phase"]], [7]
        )

    def test_an_unwritable_report_folder_does_not_stop_the_clean(self):
        value = manager(self.reportdir / "missing")

        summary = value.tclean_pass({"n": 9}, "phase_1", "bootstrap")

        self.assertEqual(summary["nmajordone"], 9)
        self.assertEqual(len(value.clean_passes["phase_1"]), 1)

    def test_an_unreadable_record_starts_empty(self):
        (self.reportdir / "clean_passes.json").write_text("[not a record")

        self.assertEqual(manager(self.reportdir).recorded_clean_passes(), {})

    def test_without_a_report_folder_the_record_stays_in_memory(self):
        value = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)

        value.tclean_pass({"n": 9}, "phase_1", "bootstrap")

        self.assertEqual(len(value.clean_passes["phase_1"]), 1)


class FlaggingRecordTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.reportdir = self.root / "reports"
        self.checkdir = self.root / "checkpoints"
        self.reportdir.mkdir()
        self.checkdir.mkdir()

    def tearDown(self):
        self.temporary.cleanup()

    def test_residual_flagging_already_done_is_still_reported(self):
        record = {"status": "applied", "new_flagged_fraction": 0.0195}
        (self.reportdir / "residual_flagging.json").write_text(
            json.dumps(record)
        )
        (self.checkdir / "residual_flagging.done").touch()
        value = manager(self.reportdir)
        value.checkdir = self.checkdir

        value.attempt_residual_flagging()

        self.assertEqual(value.residual_flagging, record)

    def test_unsolvable_antenna_flagging_is_written_down(self):
        value = manager(self.reportdir)
        value.antenna_data = {
            "status": "measured",
            "antennas": {
                "C00": {"unflagged_samples": 900.0, "solvable": True},
                "W06": {"unflagged_samples": 100.0, "solvable": False},
            },
        }

        with mock.patch.object(SELFCAL, "flagdata", lambda **kwargs: None):
            value.flag_unsolvable_antenna_data("candidate.ms")

        resumed = manager(self.reportdir)
        saved = resumed.saved_record("unsolvable_antenna_flagging.json")
        self.assertEqual(saved["antennas"], ["W06"])
        self.assertAlmostEqual(saved["fraction_of_antenna_samples"], 0.1)


if __name__ == "__main__":
    unittest.main()
