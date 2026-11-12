"""Flags written by provisional solutions must not outlive them.

With corrected-data flagging enabled, calibration is solved twice. The first
applycal exists to give the interference detectors something to read, but it
also flags every visibility whose provisional solution failed. In the uGMRT
Band 4 run that prompted these tests that took the flagged fraction from
32.7% to 49.2% before a final solution had been attempted, and nothing ever
gave those visibilities back.
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from test_postcal_flagging import CALIBRATION


def flags(*columns):
    """A (ncorr=1, nchan, nrow) flag array from per-row channel lists."""

    return np.array([list(column) for column in columns], dtype=bool).T[None]


class FakeTables:
    """Enough of casatools.table to hold an MS and its flag versions."""

    def __init__(self):
        self.columns = {}
        self.writes = []
        self.variable_shape = set()
        self.fail_on_open = None

    def add(self, path, flag, flag_row=None):
        flag = np.array(flag, dtype=bool)
        store = {"FLAG": flag}
        if flag_row is not False:
            store["FLAG_ROW"] = (
                np.zeros(flag.shape[-1], dtype=bool)
                if flag_row is None else np.array(flag_row, dtype=bool)
            )
        self.columns[str(path)] = store

    def tool(self):
        owner = self

        class Tool:

            def open(self, path, nomodify=True):
                if owner.fail_on_open and str(path).endswith(
                    owner.fail_on_open
                ):
                    raise RuntimeError(f"cannot open {path}")
                self.path = str(path)
                self.writable = not nomodify
                self.data = owner.columns[self.path]

            def close(self):
                return None

            def flush(self):
                return None

            def nrows(self):
                return self.data["FLAG"].shape[-1]

            def colnames(self):
                return list(self.data)

            def getcol(self, name, startrow=0, nrow=-1):
                if name == "FLAG" and self.path in owner.variable_shape:
                    raise RuntimeError("Table array conformance error")
                stop = None if nrow < 0 else startrow + nrow
                return self.data[name][..., startrow:stop].copy()

            def getcell(self, name, row):
                return self.data[name][..., row].copy()

            def putcol(self, name, value, startrow=0, nrow=-1):
                assert self.writable, "write to a table opened read-only"
                stop = None if nrow < 0 else startrow + nrow
                self.data[name][..., startrow:stop] = value
                owner.writes.append((self.path, name, startrow))

            def putcell(self, name, row, value):
                assert self.writable, "write to a table opened read-only"
                self.data[name][..., row] = value
                owner.writes.append((self.path, name, row))

        return Tool()


class ProvisionalFlagReleaseTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.manager = object.__new__(CALIBRATION.CalibrationManager)
        self.manager.ms = str(self.root / "hann.ms")
        self.manager.caldir = self.root / "calibration"
        self.manager.caldir.mkdir()
        self.manager.telescope = "ugmrt"
        self.manager.band = "B4"
        self.manager.task_config = {"postcal_flagging": True}
        self.tables = FakeTables()
        self.flagmanager_calls = []

        # Four rows of four channels.
        #   A: what the raw-data flagging left.
        #   B: A plus what the provisional applycal flagged (rows 1 and 3).
        #   C: B plus what tfcrop/rflag found (one channel of row 0, one of
        #      row 2).
        self.a = flags([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0])
        self.b = flags([1, 0, 0, 0], [1, 1, 1, 1], [0, 0, 0, 1], [1, 1, 1, 1])
        self.c = flags([1, 0, 1, 0], [1, 1, 1, 1], [0, 1, 0, 1], [1, 1, 1, 1])
        # What round 2 should start from: A plus the detectors' two finds.
        self.expected = flags(
            [1, 0, 1, 0], [0, 0, 0, 0], [0, 1, 0, 1], [0, 0, 0, 0]
        )

    def tearDown(self):
        self.temporary.cleanup()

    # ---- helpers ----------------------------------------------------------

    def version_path(self, name):
        return Path(self.manager.ms + ".flagversions") / f"flags.{name}"

    def install(self, row_a=None, row_b=None, row_c=None, saved=True):
        for name, flag, row in (("A", self.a, row_a), ("B", self.b, row_b)):
            path = self.version_path(name)
            path.mkdir(parents=True, exist_ok=True)
            self.tables.add(path, flag, row)
        self.tables.add(self.manager.ms, self.c, row_c)
        if saved:
            (self.manager.caldir / "round1_flag_state.json").write_text(
                json.dumps({"before_applycal_version": "A"})
            )
        (self.manager.caldir / "ugmrt_postcal_flagging_state.json").write_text(
            json.dumps({"flag_version": "B"})
        )

    def release(self):
        def fake_flagmanager(**kwargs):
            self.flagmanager_calls.append(kwargs)

        with mock.patch.object(
            CALIBRATION, "table", self.tables.tool
        ), mock.patch.object(
            CALIBRATION, "flagmanager", fake_flagmanager
        ), mock.patch("builtins.print"):
            return self.manager.release_provisional_solution_flags()

    def ms_flags(self):
        return self.tables.columns[self.manager.ms]["FLAG"]

    def report(self):
        return json.loads(
            (self.manager.caldir / "round1_flag_release_report.json").read_text()
        )

    # ---- the release itself -----------------------------------------------

    def test_only_the_detectors_flags_survive_into_round_two(self):
        self.install()

        report = self.release()

        np.testing.assert_array_equal(self.ms_flags(), self.expected)
        self.assertEqual(report["status"], "completed")
        # 12 of 16 were flagged, 4 remain: the two full rows are back.
        self.assertEqual(report["flagged_before_release"], 12)
        self.assertEqual(report["flagged_after_release"], 4)
        self.assertEqual(report["released"], 8)
        self.assertTrue(
            (self.manager.caldir / "ugmrt_round1_flags_released.done").exists()
        )

    def test_the_saved_versions_are_never_written_to(self):
        self.install()
        before_a, before_b = self.a.copy(), self.b.copy()

        self.release()

        written = {path for path, _, _ in self.tables.writes}
        self.assertEqual(written, {self.manager.ms})
        np.testing.assert_array_equal(
            self.tables.columns[str(self.version_path("A"))]["FLAG"], before_a
        )
        np.testing.assert_array_equal(
            self.tables.columns[str(self.version_path("B"))]["FLAG"], before_b
        )

    def test_the_current_flags_are_saved_before_anything_is_changed(self):
        self.install()

        report = self.release()

        self.assertEqual(self.flagmanager_calls[0]["mode"], "save")
        self.assertEqual(
            self.flagmanager_calls[0]["versionname"], report["safety_version"]
        )
        self.assertEqual(len(self.flagmanager_calls), 1)

    def test_row_flags_are_released_the_same_way(self):
        self.install(
            row_a=[0, 0, 0, 0], row_b=[0, 1, 0, 1], row_c=[1, 1, 0, 1]
        )

        self.release()

        # Row 0 was flagged by a detector and stays; rows 1 and 3 only by
        # the provisional applycal and come back.
        np.testing.assert_array_equal(
            self.tables.columns[self.manager.ms]["FLAG_ROW"],
            np.array([1, 0, 0, 0], dtype=bool),
        )

    def test_blocks_of_rows_give_the_same_answer_as_one_read(self):
        self.install()

        with mock.patch.object(CALIBRATION, "table", self.tables.tool):
            counts = self.manager.combine_released_flags(
                self.version_path("A"), self.version_path("B"),
                chunk_values=4,           # one row per block
            )

        np.testing.assert_array_equal(self.ms_flags(), self.expected)
        self.assertEqual(counts["rows"], 4)
        self.assertEqual(counts["rows_handled_singly"], 0)

    def test_rows_of_differing_shape_are_done_one_at_a_time(self):
        self.install()
        self.tables.variable_shape.add(self.manager.ms)

        report = self.release()

        np.testing.assert_array_equal(self.ms_flags(), self.expected)
        self.assertEqual(report["rows_handled_singly"], 4)

    def test_repeating_an_interrupted_release_changes_nothing(self):
        self.install()
        self.release()
        (self.manager.caldir / "ugmrt_round1_flags_released.done").unlink()

        self.release()

        np.testing.assert_array_equal(self.ms_flags(), self.expected)

    def test_flag_versions_without_row_flags_are_still_usable(self):
        self.install(row_a=False, row_b=False)

        report = self.release()

        self.assertEqual(report["status"], "completed")
        np.testing.assert_array_equal(self.ms_flags(), self.expected)

    # ---- it must never make things worse ----------------------------------

    def test_a_failure_restores_the_flags_and_lets_round_two_proceed(self):
        self.install()
        self.tables.fail_on_open = "flags.B"

        result = self.release()

        self.assertIsNone(result)
        self.assertEqual(
            [call["mode"] for call in self.flagmanager_calls],
            ["save", "restore"],
        )
        self.assertEqual(
            self.flagmanager_calls[1]["versionname"],
            self.flagmanager_calls[0]["versionname"],
        )
        self.assertEqual(self.report()["status"], "failed_restored")
        self.assertFalse(
            (self.manager.caldir / "ugmrt_round1_flags_released.done").exists()
        )

    def test_mismatched_versions_are_refused_before_any_write(self):
        self.install()
        path = str(self.version_path("B"))
        self.tables.columns[path]["FLAG"] = self.b[..., :3]
        self.tables.columns[path]["FLAG_ROW"] = np.zeros(3, dtype=bool)

        self.release()

        self.assertEqual(self.tables.writes, [])
        self.assertEqual(self.report()["status"], "failed_restored")
        self.assertIn("row count", self.report()["reason"])

    def test_a_run_started_before_this_step_existed_is_left_alone(self):
        self.install(saved=False)

        result = self.release()

        self.assertIsNone(result)
        self.assertEqual(self.tables.writes, [])
        self.assertEqual(self.flagmanager_calls, [])
        self.assertEqual(self.report()["status"], "skipped")

    def test_flags_are_not_changed_once_round_two_has_begun(self):
        self.install()
        (self.manager.caldir / "ugmrt_round2_prepared.done").touch()

        self.assertIsNone(self.release())
        self.assertEqual(self.tables.writes, [])

    def test_it_is_skipped_without_a_second_round(self):
        self.install()
        self.manager.task_config = {"postcal_flagging": False}

        self.assertIsNone(self.release())
        self.assertEqual(self.tables.writes, [])

    def test_it_can_be_switched_off(self):
        self.install()
        self.manager.task_config["release_provisional_flags"] = False

        self.assertIsNone(self.release())
        self.assertEqual(self.tables.writes, [])


class SaveBeforeRoundOneTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.manager = object.__new__(CALIBRATION.CalibrationManager)
        self.manager.ms = str(self.root / "hann.ms")
        self.manager.caldir = self.root / "calibration"
        self.manager.caldir.mkdir()
        self.manager.telescope = "ugmrt"
        self.manager.task_config = {"postcal_flagging": True}
        self.calls = []

    def tearDown(self):
        self.temporary.cleanup()

    def save(self):
        with mock.patch.object(
            CALIBRATION, "flagmanager",
            lambda **kwargs: self.calls.append(kwargs),
        ), mock.patch("builtins.print"):
            return self.manager.save_flags_before_round1_applycal()

    def test_the_state_before_the_provisional_applycal_is_saved(self):
        version = self.save()

        self.assertEqual(self.calls[0]["mode"], "save")
        self.assertEqual(self.calls[0]["versionname"], version)
        state = json.loads(
            (self.manager.caldir / "round1_flag_state.json").read_text()
        )
        self.assertEqual(state["before_applycal_version"], version)

    def test_an_interrupted_applycal_keeps_the_version_it_first_saved(self):
        first = self.save()
        (Path(self.manager.ms + ".flagversions") / f"flags.{first}").mkdir(
            parents=True
        )

        # The applycal was killed partway and is being redone: the flags
        # now include some it wrote, so they must not be saved again.
        second = self.save()

        self.assertEqual(second, first)
        self.assertEqual(len(self.calls), 1)

    def test_nothing_is_saved_once_the_first_applycal_has_run(self):
        for marker in ("applycal.done", "ugmrt_postcal_flagging.done",
                       "ugmrt_round2_prepared.done"):
            (self.manager.caldir / marker).touch()
            self.assertIsNone(self.save())
            (self.manager.caldir / marker).unlink()
        self.assertEqual(self.calls, [])

    def test_a_single_round_calibration_saves_nothing(self):
        self.manager.task_config = {"postcal_flagging": False}

        self.assertIsNone(self.save())
        self.assertEqual(self.calls, [])


class WiringTests(unittest.TestCase):

    def test_round_two_releases_before_it_archives_and_solves(self):
        manager = object.__new__(CALIBRATION.CalibrationManager)
        manager.telescope = "ugmrt"
        manager.band = "B4"
        manager.task_config = {"postcal_flagging": True}
        order = []
        for name in (
            "run_postcal_flagging", "release_provisional_solution_flags",
            "prepare_second_round", "run_delaycal", "run_prebp",
            "run_bandpass", "run_phasecal", "run_phase_transfer_cal",
            "run_ampcal", "run_fluxscale", "run_applycal",
        ):
            setattr(manager, name, lambda name=name: order.append(name))
        manager.checkpoint_exists = lambda name: False
        manager.mark_checkpoint = lambda name: None
        manager.write_atomic_json = lambda path, value: None
        manager.caldir = Path("calibration")

        with mock.patch("builtins.print"):
            manager.run_second_round()

        self.assertEqual(
            order[:4],
            ["run_postcal_flagging", "release_provisional_solution_flags",
             "prepare_second_round", "run_delaycal"],
        )
        self.assertEqual(order[-1], "run_applycal")


if __name__ == "__main__":
    unittest.main()
