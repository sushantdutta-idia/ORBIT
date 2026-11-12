"""One failed calibrator solution must not cost a target the scans around it.

applycal flags a visibility when either of the two solutions bracketing it is
flagged. In the uGMRT Band 4 run that prompted these tests, the phase
calibrator was visited ten times; five antennas with one to three failed
solutions kept 6 to 7% of their data on a target, and the target was imaged
from 26% of its visibilities.
"""

import importlib.util
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from test_postcal_flagging import CALIBRATION
from test_calibration_weights import calibrator


SOURCE_DIR = Path(__file__).resolve().parents[1]
HOUR = 3600.0
# Ten visits to the phase calibrator, 54 minutes apart.
VISITS = np.arange(10) * 54.0 * 60.0

bridge = CALIBRATION.CalibrationManager.bridge_solution_series


def phasor(degrees, amplitude=1.0):
    return amplitude * np.exp(1j * np.deg2rad(degrees))


class SeriesTests(unittest.TestCase):

    def series(self, failed, gains=None):
        gains = (
            np.array([phasor(10.0 * index) for index in range(VISITS.size)])
            if gains is None else np.array(gains, dtype=complex)
        )
        flags = np.zeros(VISITS.size, dtype=bool)
        flags[list(failed)] = True
        # CASA leaves a flagged solution at 1+0j; it must never be used.
        gains[flags] = 1.0 + 0.0j
        return gains, flags

    def test_a_failure_between_two_good_solutions_is_interpolated(self):
        gains, flags = self.series([4])
        new_gains, new_flags, fills = bridge(VISITS, gains, flags, 2 * HOUR)

        self.assertFalse(new_flags.any())
        self.assertAlmostEqual(np.rad2deg(np.angle(new_gains[4])), 40.0)
        self.assertAlmostEqual(abs(new_gains[4]), 1.0)
        self.assertEqual(len(fills), 1)
        self.assertEqual(fills[0]["kind"], "interpolated")
        self.assertAlmostEqual(fills[0]["span_s"], 108.0 * 60.0)
        # The inputs are not modified in place.
        self.assertTrue(flags[4])

    def test_amplitude_is_interpolated_linearly(self):
        gains = [phasor(0.0, 2.0)] * VISITS.size
        gains[3] = phasor(0.0, 1.0)
        gains[5] = phasor(0.0, 3.0)
        gains, flags = self.series([4], gains)
        new_gains, _, _ = bridge(VISITS, gains, flags, 2 * HOUR)
        self.assertAlmostEqual(abs(new_gains[4]), 2.0)

    def test_phase_takes_the_shorter_arc(self):
        gains = [phasor(0.0)] * VISITS.size
        gains[3] = phasor(170.0)
        gains[5] = phasor(-170.0)
        gains, flags = self.series([4], gains)
        new_gains, _, _ = bridge(VISITS, gains, flags, 2 * HOUR)
        # Half-way from 170 to -170 degrees is 180, not 0.
        self.assertAlmostEqual(abs(np.rad2deg(np.angle(new_gains[4]))), 180.0)

    def test_two_failures_in_a_row_are_too_far_apart_by_default(self):
        gains, flags = self.series([4, 5])
        new_gains, new_flags, fills = bridge(VISITS, gains, flags, 2 * HOUR)
        # Their good neighbours are 162 minutes apart.
        self.assertEqual(fills, [])
        self.assertTrue(new_flags[4] and new_flags[5])

    def test_a_wider_limit_bridges_them_from_the_good_solutions_only(self):
        gains, flags = self.series([4, 5])
        new_gains, new_flags, fills = bridge(VISITS, gains, flags, 3 * HOUR)
        self.assertFalse(new_flags.any())
        # Both come from solutions 3 and 6, not from each other.
        self.assertAlmostEqual(np.rad2deg(np.angle(new_gains[4])), 40.0)
        self.assertAlmostEqual(np.rad2deg(np.angle(new_gains[5])), 50.0)
        self.assertEqual(len(fills), 2)

    def test_a_failed_first_visit_takes_the_next_one(self):
        # The case of a target observed between the first two visits.
        gains, flags = self.series([0])
        new_gains, new_flags, fills = bridge(VISITS, gains, flags, 2 * HOUR)
        self.assertFalse(new_flags[0])
        self.assertAlmostEqual(new_gains[0], phasor(10.0))
        self.assertEqual(fills[0]["kind"], "nearest")

    def test_an_edge_is_not_extended_beyond_half_the_limit(self):
        gains, flags = self.series([0, 1])
        _, new_flags, fills = bridge(VISITS, gains, flags, 2 * HOUR)
        # Visit 1 is 54 minutes from a good one and is filled; visit 0 is
        # 108 minutes from it and is not.
        self.assertFalse(new_flags[1])
        self.assertTrue(new_flags[0])
        self.assertEqual(len(fills), 1)

    def test_a_failing_antenna_is_left_flagged(self):
        gains, flags = self.series([0, 2, 4, 6, 8, 9])
        _, new_flags, fills = bridge(VISITS, gains, flags, 2 * HOUR)
        self.assertEqual(fills, [])
        self.assertEqual(int(new_flags.sum()), 6)

    def test_nothing_to_do_is_nothing_done(self):
        gains, flags = self.series([])
        self.assertEqual(bridge(VISITS, gains, flags, 2 * HOUR)[2], [])
        gains, flags = self.series(range(VISITS.size))
        self.assertEqual(bridge(VISITS, gains, flags, 2 * HOUR)[2], [])

    def test_a_non_finite_solution_counts_as_failed(self):
        gains, flags = self.series([])
        gains[4] = complex(float("nan"), 0.0)
        new_gains, _, fills = bridge(VISITS, gains, flags, 2 * HOUR)
        self.assertEqual(len(fills), 1)
        self.assertTrue(np.isfinite(new_gains[4]))


class GainTable:
    """A T-Jones gain table as the table tool presents it."""

    def __init__(self, antennas=3, fail=None):
        # Field 1 is the phase calibrator (ten visits); field 0 the flux
        # calibrator, seen once.
        rows = []
        for antenna in range(antennas):
            rows.append((0.0 - 600.0, antenna, 0))
            rows.extend((time, antenna, 1) for time in VISITS)
        self.columns = {
            "TIME": np.array([row[0] for row in rows]),
            "ANTENNA1": np.array([row[1] for row in rows]),
            "FIELD_ID": np.array([row[2] for row in rows]),
            "SPECTRAL_WINDOW_ID": np.zeros(len(rows), dtype=int),
            "CPARAM": np.ones((1, 1, len(rows)), dtype=complex) * phasor(30.0),
            "FLAG": np.zeros((1, 1, len(rows)), dtype=bool),
        }
        self.fail = fail
        self.opened = []
        self.written = {}

    def flag(self, antenna, field, visit):
        row = np.flatnonzero(
            (self.columns["ANTENNA1"] == antenna) &
            (self.columns["FIELD_ID"] == field)
        )[visit]
        self.columns["FLAG"][0, 0, row] = True
        self.columns["CPARAM"][0, 0, row] = 1.0
        return int(row)

    def open(self, path, nomodify=True):
        if self.fail:
            raise RuntimeError(self.fail)
        self.opened.append((path, nomodify))

    def close(self):
        return None

    def getcol(self, name):
        return self.columns[name].copy()

    def putcol(self, name, value):
        self.written[name] = np.array(value)


class TableTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.source = self.root / "phase_scan.cal"
        self.source.mkdir()
        (self.source / "table.dat").write_text("solved")
        self.destination = self.root / "phase_scan.target.cal"
        self.manager = calibrator(self.root, target_transfer_bridge=True)
        self.manager.write_atomic_json = lambda path, data: Path(
            path
        ).write_text(json.dumps(data))

    def tearDown(self):
        self.temporary.cleanup()

    def bridge(self, gain_table):
        with mock.patch.object(CALIBRATION, "table", lambda: gain_table), \
                mock.patch("builtins.print") as printed:
            record = self.manager.bridge_transfer_table(
                self.source, self.destination
            )
        return record, " ".join(
            str(argument) for call in printed.call_args_list
            for argument in call.args
        )

    def test_the_copy_is_bridged_and_the_solved_table_is_untouched(self):
        gain_table = GainTable()
        row = gain_table.flag(antenna=1, field=1, visit=4)
        record, _ = self.bridge(gain_table)

        # The copy, not the solved table, is opened for writing.
        self.assertEqual(gain_table.opened, [(str(self.destination), False)])
        self.assertEqual((self.source / "table.dat").read_text(), "solved")
        self.assertTrue(self.destination.is_dir())

        self.assertFalse(gain_table.written["FLAG"].any())
        self.assertAlmostEqual(
            gain_table.written["CPARAM"][0, 0, row], phasor(30.0)
        )
        self.assertEqual(record["bridged"], 1)
        self.assertEqual(record["flagged_before"], 1)
        self.assertEqual(record["antennas"], [1])
        fill = record["fills"][0]
        self.assertEqual((fill["antenna_id"], fill["field_id"], fill["row"]),
                         (1, 1, row))
        self.assertAlmostEqual(fill["span_minutes"], 108.0)

    def test_fields_are_bridged_separately(self):
        # The flux calibrator's single solution has no neighbour in its own
        # field; the phase calibrator's solutions are not borrowed for it.
        gain_table = GainTable()
        row = gain_table.flag(antenna=2, field=0, visit=0)
        record, _ = self.bridge(gain_table)

        self.assertEqual(record["bridged"], 0)
        self.assertNotIn("FLAG", gain_table.written)
        self.assertTrue(gain_table.columns["FLAG"][0, 0, row])

    def test_a_table_that_cannot_be_read_falls_back_to_the_solved_one(self):
        record, printed = self.bridge(GainTable(fail="table is locked"))

        self.assertIsNone(record)
        self.assertFalse(self.destination.exists())
        self.assertIn("Targets receive the table as solved", printed)

    def test_targets_get_the_bridged_tables_and_a_record_of_them(self):
        flux = self.root / "flux.cal"
        flux.mkdir()
        gain_table = GainTable()
        gain_table.flag(antenna=0, field=1, visit=0)
        with mock.patch.object(CALIBRATION, "table", lambda: gain_table), \
                mock.patch("builtins.print"):
            tables = self.manager.target_transfer_tables(self.source, flux)

        self.assertEqual(tables["phase"], str(self.destination))
        self.assertEqual(tables["flux"], str(self.root / "flux.target.cal"))
        report = json.loads(
            (self.root / "target_transfer_bridge.json").read_text()
        )
        self.assertEqual(report["tables"]["phase"]["status"], "bridged")
        self.assertEqual(report["tables"]["phase"]["bridged"], 1)
        self.assertEqual(report["max_gap_minutes"], 120.0)

    def test_switched_off_the_targets_get_the_tables_as_solved(self):
        self.manager.task_config["target_transfer_bridge"] = False
        with mock.patch.object(
            CALIBRATION, "table",
            lambda: self.fail("no table should be opened"),
        ):
            tables = self.manager.target_transfer_tables(
                self.source, self.root / "flux.cal"
            )
        self.assertEqual(tables["phase"], str(self.source))
        self.assertEqual(tables["flux"], str(self.root / "flux.cal"))
        self.assertFalse((self.root / "target_transfer_bridge.json").exists())


class ApplicationTests(unittest.TestCase):

    def test_calibrators_keep_the_solved_tables_and_targets_do_not(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manager = calibrator(
                root, applymode="calflagstrict", parang=False, calwt=True,
            )
            manager.ms = str(root / "working.ms")
            manager.checkpoint_exists = lambda name: False
            manager.mark_checkpoint = lambda name: None
            manager.get_phase_calibrators = lambda: ["0323+055"]
            manager.get_flux_calibrators = lambda: ["3C48"]
            manager.get_bandpass_calibrator = lambda: "3C48"
            manager.get_targets = lambda: ["G0212-0335"]
            manager.get_target_phase_map = lambda: {"G0212-0335": "0323+055"}
            manager.target_transfer_tables = lambda phase, flux: {
                "phase": str(root / "phase_scan.target.cal"),
                "flux": str(root / "flux.target.cal"),
            }
            calls = []
            with mock.patch.object(
                CALIBRATION, "applycal", lambda **kwargs: calls.append(kwargs)
            ), mock.patch("builtins.print"):
                manager.run_applycal()

        by_field = {call["field"]: call for call in calls}
        self.assertEqual(set(by_field), {"3C48", "0323+055", "G0212-0335"})
        for name in ("3C48", "0323+055"):
            tables = [Path(item).name for item in by_field[name]["gaintable"]]
            self.assertEqual(
                tables, ["delay.cal", "bandpass.cal", "phase.cal", "flux.cal"]
            )
        target = by_field["G0212-0335"]
        self.assertEqual(
            [Path(item).name for item in target["gaintable"]],
            ["delay.cal", "bandpass.cal", "phase_scan.target.cal",
             "flux.target.cal"],
        )
        self.assertEqual(
            target["gainfield"], ["3C48", "3C48", "0323+055", "0323+055"]
        )
        # Whatever the final solutions cannot calibrate is still flagged.
        self.assertEqual(target["applymode"], "calflagstrict")


class ConfigurationTests(unittest.TestCase):

    def test_ugmrt_bridges_and_the_other_telescopes_are_unchanged(self):
        spec = importlib.util.spec_from_file_location(
            "test_radio_continuum_config_bridge", SOURCE_DIR / "config.py"
        )
        config = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(config)
        ugmrt = config.TELESCOPE_CALIBRATION_TASK_CONFIG["ugmrt"]
        self.assertTrue(ugmrt["target_transfer_bridge"])
        self.assertEqual(ugmrt["target_transfer_bridge_max_gap_minutes"], 120.0)
        for telescope in ("meerkat", "lofar"):
            self.assertNotIn(
                "target_transfer_bridge",
                config.TELESCOPE_CALIBRATION_TASK_CONFIG[telescope],
            )


if __name__ == "__main__":
    unittest.main()
