"""An antenna with no signal at the bandpass calibrator is excluded, not fatal.

The case is the Band 5 run of 2026-10-20: C03 and C05 failed SNR 3 in every
channel of the bandpass solve (C05 kept 1% of its solutions), after C12 and
W06 had been excluded at the delay and pre-bandpass solves. On the rerun of
2026-11-03, with the 3C48 scan cut to its 9.7 settled minutes (finding 161),
noise alone kept 10.0% and 9.5% of their channels: past the old fixed limit
of 5% usable, yet failing the per-antenna check (more than 80% flagged), so
the calibration stopped. An antenna is now excluded whenever it alone fails
that check.
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_calibration_table_qa import CALIBRATION


def antenna_group(value, usable, total=8192):
    flagged = total - usable
    group = {
        "value": value,
        "rows": 1,
        "total_solutions": total,
        "flagged_solutions": flagged,
        "usable_solutions": usable,
        "flagged_fraction": flagged / total,
        "usable_fraction": usable / total,
        "status": "passed",
        "reason": "",
    }
    if usable == 0:
        group.update(status="failed", reason="group has no usable solution")
    elif flagged / total > 0.8:
        group.update(
            status="failed",
            reason=(
                f"flagged fraction {flagged / total:.6f} exceeds the "
                "configured group maximum 0.800000"
            ),
        )
    return group


class BandpassQuarantineTests(unittest.TestCase):

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
        manager.flagging = {"bad_antennas": ["C11", "C29"]}
        manager.delay_quarantine = {
            "schema_version": 1,
            "measurement_set": str(Path(manager.input_ms).resolve()),
            "antennas": [
                {"id": 11, "name": "C11"},
                {"id": 29, "name": "C29"},
            ],
        }
        manager.task_config = {
            "auto_quarantine_failed_delay_antennas": True,
            "max_auto_quarantine_delay_antennas": 2,
            "max_auto_quarantine_delay_fraction": 0.15,
            "min_antennas_after_delay_quarantine": 8,
            "prebp_gaintype": "G",
            "gain_minsnr": 3.0,
            "minblperant": 4,
            "gain_solmode": "",
            "parang": False,
            "delay_minsnr": 3.0,
            "bandpass_fillgaps": 8,
            "bandpass_solint": "inf",
            "bandpass_solnorm": True,
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
        self.checkpoints = []
        manager.mark_checkpoint = self.checkpoints.append
        self.manager = manager

    def tearDown(self):
        self.temporary.cleanup()

    def write_bandpass_qa(self, usable_by_antenna):
        groups = [
            antenna_group(value, usable_by_antenna.get(value, 7704))
            for value in range(30)
            if value not in (11, 29)
        ]
        qa = {
            "nonfinite_usable_solutions": 0,
            "group_metrics": {
                "ANTENNA1": {
                    "status": (
                        "failed"
                        if any(g["status"] == "failed" for g in groups)
                        else "passed"
                    ),
                    "groups": groups,
                },
                "SPECTRAL_WINDOW_ID": {"status": "passed", "groups": []},
            },
        }
        path = self.manager.caltable_qa_path("bandpass.cal")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(qa))

    def run_bandpass(self, bandpass_results):
        manager = self.manager
        calls = []
        results = list(bandpass_results)

        def solve(**kwargs):
            calls.append(kwargs)
            table = Path(kwargs["caltable"])
            table.mkdir(exist_ok=True)
            (table / "table.dat").touch()

        def status(filename):
            if filename == "bandpass.cal":
                return results.pop(0)
            return True, ""

        manager.caltable_status = status
        manager.require_caltable = lambda filename: None
        with mock.patch.object(CALIBRATION, "gaincal", solve), \
                mock.patch.object(CALIBRATION, "bandpass", solve), \
                mock.patch("builtins.print"):
            manager.run_bandpass()
        return calls

    def test_dead_and_nearly_dead_antennas_are_found(self):
        # C03 has nothing, C05 1%: both count for the bandpass.
        self.write_bandpass_qa({3: 0, 5: 82})
        self.assertEqual(
            self.manager.failed_delay_antenna_ids(
                "bandpass.cal",
                max_usable_fraction=(
                    self.manager.bandpass_quarantine_max_usable_fraction()
                ),
            ),
            [3, 5],
        )
        # The delay and phase solves still require nothing at all.
        self.assertEqual(self.manager.failed_delay_antenna_ids("bandpass.cal"), [])

    def test_the_limit_is_where_the_per_antenna_check_fails(self):
        self.assertAlmostEqual(
            self.manager.bandpass_quarantine_max_usable_fraction(), 0.20)
        self.manager.task_config["max_caltable_group_flagged_fraction"] = 0.9
        self.assertAlmostEqual(
            self.manager.bandpass_quarantine_max_usable_fraction(), 0.10)

    def test_noise_passing_in_a_tenth_of_the_channels_is_still_dead(self):
        # The rerun of 2026-11-03: 10.0% and 9.5% usable on 9.7 minutes.
        self.write_bandpass_qa({3: 819, 5: 778})
        self.assertEqual(
            self.manager.failed_delay_antenna_ids(
                "bandpass.cal",
                max_usable_fraction=(
                    self.manager.bandpass_quarantine_max_usable_fraction()
                ),
            ),
            [3, 5],
        )

    def test_an_antenna_that_passes_the_check_is_kept(self):
        # 25% usable passes the per-antenna check, so nothing failed there and
        # nothing is excluded, whatever else failed.
        self.write_bandpass_qa({3: 0, 5: 2048})
        self.assertEqual(
            self.manager.failed_delay_antenna_ids(
                "bandpass.cal",
                max_usable_fraction=(
                    self.manager.bandpass_quarantine_max_usable_fraction()
                ),
            ),
            [3],
        )

    def test_they_are_excluded_and_all_three_solves_redone(self):
        self.write_bandpass_qa({3: 0, 5: 82})
        calls = self.run_bandpass([
            (False, "ANTENNA1=3: group has no usable solution"),
            (True, ""),
        ])

        tables = [Path(call["caltable"]).name for call in calls]
        self.assertEqual(
            tables, ["bandpass.cal", "delay.cal", "prebp.cal", "bandpass.cal"]
        )
        self.assertEqual(calls[0]["antenna"], "!C11,C29")
        for call in calls[1:]:
            self.assertEqual(call["antenna"], "!C11,C29,C03,C05")
        self.assertEqual(self.checkpoints, ["bandpass"])

        caldir = self.manager.caldir
        self.assertTrue((caldir / "bandpass.cal.incomplete.1").is_dir())
        record = json.loads(self.manager.delay_quarantine_path().read_text())
        names = [item["name"] for item in record["antennas"]]
        self.assertEqual(names, ["C03", "C05", "C11", "C29"])
        new = [item for item in record["antennas"] if item["name"] in ("C03", "C05")]
        self.assertTrue(all(
            item["reason"].startswith("no usable bandpass solution")
            for item in new
        ))

    def test_a_good_bandpass_is_solved_once(self):
        calls = self.run_bandpass([(True, "")])
        self.assertEqual([Path(c["caltable"]).name for c in calls], ["bandpass.cal"])
        self.assertEqual(self.checkpoints, ["bandpass"])

    def test_the_limit_that_stops_it_is_named_with_what_would_pass(self):
        # The default 10% of the array: 4 of 30 is 13%.
        self.manager.task_config["max_auto_quarantine_delay_fraction"] = 0.10
        self.write_bandpass_qa({3: 0, 5: 82})
        with self.assertRaises(RuntimeError) as raised:
            self.run_bandpass([
                (False, "ANTENNA1=3: group has no usable solution"),
            ])
        message = str(raised.exception)
        self.assertIn("C03, C05", message)
        self.assertIn("with the 2 excluded before (C11, C29)", message)
        self.assertIn("4 of 30 antennas (13%)", message)
        self.assertIn("calibration_max_auto_quarantine_delay_fraction (0.1)", message)
        self.assertIn("0.14 would allow it", message)
        self.assertFalse(self.manager.delay_quarantine_path().exists())

    def test_too_many_at_once_is_named(self):
        self.manager.task_config["max_auto_quarantine_delay_antennas"] = 1
        self.write_bandpass_qa({3: 0, 5: 82})
        with self.assertRaisesRegex(
            RuntimeError,
            r"2 at once is more than "
            r"calibration_max_auto_quarantine_delay_antennas \(1\)",
        ):
            self.run_bandpass([
                (False, "ANTENNA1=3: group has no usable solution"),
            ])


if __name__ == "__main__":
    unittest.main()
