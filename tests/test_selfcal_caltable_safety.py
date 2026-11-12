import importlib.util
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SOURCE_DIR = Path(__file__).resolve().parents[1]


def _unused_task(*args, **kwargs):
    return None


class _ImportOnlyTable:
    pass


def load_selfcal_module():
    casatasks = types.ModuleType("casatasks")
    for name in [
        "applycal",
        "flagdata",
        "flagmanager",
        "gaincal",
        "imhead",
        "imstat",
        "rmtables",
        "split",
        "tclean",
    ]:
        setattr(casatasks, name, _unused_task)

    casatools = types.ModuleType("casatools")
    casatools.table = _ImportOnlyTable
    module_path = SOURCE_DIR / "selfcal.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_selfcal_caltable_safety",
        module_path,
    )
    module = importlib.util.module_from_spec(spec)
    original_sys_path = list(sys.path)
    try:
        with mock.patch.dict(
            sys.modules,
            {"casatasks": casatasks, "casatools": casatools},
        ):
            spec.loader.exec_module(module)
    finally:
        sys.path[:] = original_sys_path
    return module


SELFCAL = load_selfcal_module()


class FakeTableTool:

    def __init__(self, dataset):
        self.dataset = dataset

    def open(self, path):
        return None

    def close(self):
        return None

    def nrows(self):
        return self.dataset["nrows"]

    def colnames(self):
        return list(self.dataset["columns"])

    def getcol(self, name):
        return self.dataset["columns"][name]


class PathTableTool:
    """A table tool that serves a different table for each path."""

    def __init__(self, datasets, failing=()):
        self.datasets = {
            str(Path(path).resolve()): dataset
            for path, dataset in datasets.items()
        }
        self.failing = {str(Path(path).resolve()) for path in failing}
        self.opened = []
        self.dataset = None

    def open(self, path):
        self.opened.append(str(path))
        key = str(Path(path).resolve())
        if key in self.failing:
            raise RuntimeError(f"cannot open {path}")
        self.dataset = self.datasets[key]

    def close(self):
        self.dataset = None

    def nrows(self):
        return self.dataset["nrows"]

    def colnames(self):
        return list(self.dataset["columns"])

    def getcol(self, name):
        return self.dataset["columns"][name]


class SelfCalCaltableSafetyTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.caltabledir = self.root / "caltables"
        self.caltabledir.mkdir()
        self.table = self.caltabledir / "phase_1.cal"
        self.table.mkdir()
        (self.table / "table.dat").touch()

        self.manager = SELFCAL.SelfCalManager.__new__(
            SELFCAL.SelfCalManager
        )
        self.manager.caltabledir = self.caltabledir
        self.manager.caltable_qa = {}
        self.manager.task_config = {
            "max_caltable_unusable_fraction": 0.50,
            "max_caltable_group_unusable_fraction": 0.80,
        }

    def tearDown(self):
        self.temporary.cleanup()

    @staticmethod
    def columns(flags=None, values=None):
        if values is None:
            values = np.ones((2, 1, 4), dtype=np.complex128)
        if flags is None:
            flags = np.zeros(values.shape, dtype=bool)
        return {
            "CPARAM": values,
            "FLAG": flags,
            "ANTENNA1": np.asarray([0, 0, 1, 1]),
            "SPECTRAL_WINDOW_ID": np.asarray([0, 1, 0, 1]),
            "FIELD_ID": np.asarray([7, 7, 7, 7]),
            "OBSERVATION_ID": np.asarray([0, 0, 0, 0]),
            "TIME": np.asarray([1.0, 2.0, 1.0, 2.0]),
        }

    def set_columns(self, columns):
        SELFCAL.tb = FakeTableTool({
            "nrows": 4,
            "columns": columns,
        })

    def qa(self):
        return json.loads(
            (self.caltabledir / "phase_1.cal.qa.json").read_text()
        )

    def test_healthy_table_passes_and_records_group_metrics(self):
        self.set_columns(self.columns())

        valid, reason = self.manager.caltable_status(self.table)

        self.assertTrue(valid, reason)
        qa = self.qa()
        self.assertEqual(qa["status"], "passed")
        self.assertEqual(qa["usable_solutions"], 8)
        self.assertIn("ANTENNA1", qa["group_metrics"])
        self.assertIn(str(self.table.resolve()), self.manager.caltable_qa)

    def test_global_sparse_solution_table_is_rejected(self):
        flags = np.zeros((2, 1, 4), dtype=bool)
        flags.reshape(-1)[:5] = True
        self.set_columns(self.columns(flags=flags))

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("unusable calibration-solution fraction", reason)

    def test_bad_antenna_group_is_rejected_even_at_global_limit(self):
        flags = np.zeros((2, 1, 4), dtype=bool)
        flags[..., :2] = True
        self.set_columns(self.columns(flags=flags))

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("ANTENNA1=0", reason)

    def test_missing_flag_column_is_rejected(self):
        columns = self.columns()
        del columns["FLAG"]
        self.set_columns(columns)

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("no FLAG column", reason)

    def test_misaligned_flag_shape_is_rejected(self):
        columns = self.columns()
        columns["FLAG"] = np.zeros((1, 1, 4), dtype=bool)
        self.set_columns(columns)

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("shapes differ", reason)

    def test_unflagged_invalid_gain_is_rejected(self):
        values = np.ones((2, 1, 4), dtype=np.complex128)
        values[0, 0, 0] = complex(float("nan"), 0.0)
        self.set_columns(self.columns(values=values))

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("unflagged non-finite", reason)

    # ---- antennas that hold no data must not veto the solution -----------

    @staticmethod
    def array_columns(antennas=10, intervals=7, dead=(), failed=()):
        """A T-type table: one row per antenna and interval.

        ``dead`` antennas return only flagged solutions; ``failed`` is a set
        of (antenna, interval) solutions flagged for want of S/N.
        """

        rows = antennas * intervals
        antenna = np.repeat(np.arange(antennas), intervals)
        time = np.tile(np.arange(intervals, dtype=float), antennas)
        flags = np.zeros((1, 1, rows), dtype=bool)
        for index in range(rows):
            if antenna[index] in dead or (
                (int(antenna[index]), int(time[index])) in failed
            ):
                flags[0, 0, index] = True
        return rows, {
            "CPARAM": np.ones((1, 1, rows), dtype=np.complex128),
            "FLAG": flags,
            "ANTENNA1": antenna,
            "SPECTRAL_WINDOW_ID": np.zeros(rows, dtype=int),
            "FIELD_ID": np.zeros(rows, dtype=int),
            "OBSERVATION_ID": np.zeros(rows, dtype=int),
            "TIME": time,
        }

    def set_array(self, **kwargs):
        rows, columns = self.array_columns(**kwargs)
        SELFCAL.tb = FakeTableTool({"nrows": rows, "columns": columns})

    def test_flagged_antennas_reject_a_table_when_they_are_not_known(self):
        # The run this guards against: three quarantined antennas and five
        # that lost their calibration upstream, eight of thirty in all,
        # returned flagged solutions, and the per-antenna test rejected the
        # healthy solutions of the other twenty-two.
        self.set_array(antennas=30, dead=range(8))

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("ANTENNA1=0", reason)

    def test_antennas_without_data_are_left_out_of_the_coverage_test(self):
        self.set_array(antennas=30, dead=range(8))
        self.manager.unsolvable_antenna_ids = set(range(8))

        valid, reason = self.manager.caltable_status(self.table)

        self.assertTrue(valid, reason)
        qa = self.qa()
        self.assertEqual(qa["excluded_unsolvable_antenna_ids"], list(range(8)))
        self.assertEqual(qa["excluded_solution_rows"], 8 * 7)
        # Coverage is judged on the antennas that could be solved.
        self.assertEqual(qa["total_solutions"], 22 * 7)
        self.assertEqual(qa["unusable_fraction"], 0.0)
        self.assertEqual(
            len(qa["group_metrics"]["ANTENNA1"]), 22
        )

    def test_a_few_unsolved_antennas_are_tolerated_and_recorded(self):
        # Two of twenty antennas with data fail every interval.
        self.set_array(antennas=20, dead=(3, 11))

        valid, reason = self.manager.caltable_status(self.table)

        self.assertTrue(valid, reason)
        self.assertEqual(
            self.qa()["tolerated_failed_groups"], {"ANTENNA1": [3, 11]}
        )

    def test_too_many_unsolved_antennas_still_reject_the_table(self):
        # Six of twenty is beyond the 25% that may go uncorrected.
        self.set_array(antennas=20, dead=range(6))

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("6 of 20 ANTENNA1 groups failed", reason)

    def test_the_tolerance_is_configurable(self):
        self.set_array(antennas=20, dead=(3,))
        self.manager.task_config["max_caltable_failed_antenna_fraction"] = 0.0

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("ANTENNA1=3", reason)

    def test_one_lost_interval_is_tolerated(self):
        failed = {(antenna, 0) for antenna in range(10)}
        self.set_array(antennas=10, intervals=7, failed=failed)

        valid, reason = self.manager.caltable_status(self.table)

        self.assertTrue(valid, reason)
        self.assertEqual(
            self.qa()["tolerated_failed_groups"], {"TIME": [0.0]}
        )

    # ---- a failed interval is weighed by the data it covers --------------

    def set_array_and_ms(self, integrations, ms_fails=False, **kwargs):
        """A table and the Measurement Set it was solved from.

        ``integrations`` maps an integration time to its number of rows.
        """

        rows, columns = self.array_columns(**kwargs)
        ms = self.root / "phase4.ms"
        times = np.concatenate([
            np.full(count, float(time))
            for time, count in integrations.items()
        ])
        SELFCAL.tb = PathTableTool({
            str(self.table): {"nrows": rows, "columns": columns},
            str(ms): {"nrows": len(times), "columns": {"TIME": times}},
        }, failing=[str(ms)] if ms_fails else [])
        self.manager.current_ms = str(ms)

    @staticmethod
    def scans(*lengths, rows=5):
        """Integrations 10 s apart in scans an hour apart.

        Solution ``index`` of an ``array_columns`` table sits at time
        ``index``, so scan ``index`` is centred there, in units where an
        hour is one.
        """

        integrations = {}
        for index, length in enumerate(lengths):
            for step in range(length):
                offset = (step - (length - 1) / 2.0) * (10.0 / 3600.0)
                integrations[index + offset] = rows
        return integrations

    def test_a_short_failed_scan_does_not_cost_the_whole_table(self):
        # uGMRT Band 4, G0212-0335: one amplitude solution per scan over
        # scans of 53, 125 and 8 integrations. No antenna reached the
        # signal-to-noise threshold in the last, which is a third of the
        # intervals and 4% of the observation.
        failed = {(antenna, 2) for antenna in range(27)}
        self.set_array_and_ms(
            self.scans(53, 125, 8), antennas=27, intervals=3, failed=failed
        )

        valid, reason = self.manager.caltable_status(self.table)

        self.assertTrue(valid, reason)
        qa = self.qa()
        self.assertEqual(qa["tolerated_failed_groups"], {"TIME": [2.0]})
        self.assertAlmostEqual(qa["failed_observing_time_share"], 8 / 186)
        self.assertEqual(
            [
                round(group["observing_time_share"] * 186)
                for group in qa["group_metrics"]["TIME"]
            ],
            [53, 125, 8],
        )

    def test_a_long_failed_scan_is_not_hidden_by_short_ones(self):
        # One interval of five is within the 25% by count; it is the scan
        # that holds most of the data.
        failed = {(antenna, 1) for antenna in range(10)}
        self.set_array_and_ms(
            self.scans(10, 120, 10, 10, 10),
            antennas=10, intervals=5, failed=failed,
        )

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("TIME=1.0", reason)
        self.assertIn("1 of 5 TIME groups failed", reason)
        self.assertIn("covering 75.0% of the observing time", reason)

    def test_rows_not_integrations_set_the_share(self):
        # The failed scan has as many integrations as the others and a
        # tenth of their baselines.
        failed = {(antenna, 0) for antenna in range(10)}
        integrations = self.scans(20, 20, 20, rows=50)
        for time in list(integrations)[:20]:
            integrations[time] = 5
        self.set_array_and_ms(
            integrations, antennas=10, intervals=3, failed=failed
        )

        valid, reason = self.manager.caltable_status(self.table)

        self.assertTrue(valid, reason)
        self.assertAlmostEqual(
            self.qa()["failed_observing_time_share"], 100 / 2100
        )

    def test_without_a_measurement_set_the_intervals_are_counted(self):
        failed = {(antenna, 2) for antenna in range(27)}
        self.set_array(antennas=27, intervals=3, failed=failed)

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("1 of 3 TIME groups failed; at most 25%", reason)
        self.assertNotIn("failed_observing_time_share", self.qa())

    def test_an_unreadable_measurement_set_falls_back_to_counting(self):
        failed = {(antenna, 2) for antenna in range(27)}
        self.set_array_and_ms(
            self.scans(53, 125, 8), ms_fails=True,
            antennas=27, intervals=3, failed=failed,
        )

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("1 of 3 TIME groups failed; at most 25%", reason)

    def test_a_measurement_set_from_another_time_is_not_used(self):
        failed = {(antenna, 2) for antenna in range(27)}
        elsewhere = {
            time + 1000.0: rows
            for time, rows in self.scans(53, 125, 8).items()
        }
        self.set_array_and_ms(
            elsewhere, antennas=27, intervals=3, failed=failed
        )

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("1 of 3 TIME groups failed; at most 25%", reason)

    def test_the_antenna_tolerance_still_counts_antennas(self):
        self.set_array_and_ms(
            self.scans(53, 125, 8), antennas=20, intervals=3, dead=range(6)
        )

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("6 of 20 ANTENNA1 groups failed; at most 25%", reason)

    def test_the_integrations_of_a_measurement_set_are_read_once(self):
        self.set_array_and_ms(self.scans(53, 125, 8), antennas=5, intervals=3)

        self.manager.caltable_status(self.table)
        self.manager.caltable_status(self.table)

        self.assertEqual(
            SELFCAL.tb.opened.count(self.manager.current_ms), 1
        )

    def test_a_table_of_only_unsolvable_antennas_is_rejected(self):
        self.set_array(antennas=4, dead=range(4))
        self.manager.unsolvable_antenna_ids = set(range(4))

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("no usable data", reason)

    def test_a_spectral_window_without_solutions_is_never_tolerated(self):
        rows, columns = self.array_columns(antennas=10, intervals=4)
        spw = np.zeros(rows, dtype=int)
        spw[::4] = 1                      # one interval per antenna in SPW 1
        columns["SPECTRAL_WINDOW_ID"] = spw
        columns["FLAG"][0, 0, ::4] = True  # and none of it solved
        SELFCAL.tb = FakeTableTool({"nrows": rows, "columns": columns})

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("SPECTRAL_WINDOW_ID=1", reason)

    def test_required_group_columns_are_fail_closed(self):
        columns = self.columns()
        del columns["SPECTRAL_WINDOW_ID"]
        self.set_columns(columns)

        valid, reason = self.manager.caltable_status(self.table)

        self.assertFalse(valid)
        self.assertIn("SPECTRAL_WINDOW_ID", reason)


class AntennaDataCensusTests(unittest.TestCase):
    """Which antennas have data to solve from is measured, not assumed."""

    NAMES = ["C00", "C01", "C02", "C03", "E02", "S01", "W01", "W06"]

    def manager(self):
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.task_config = {}
        manager.reference_antenna = "C00"
        return manager

    def census(self, counts, manager=None):
        manager = manager or self.manager()
        summary = {
            "antenna": {
                name: {"total": 100.0, "flagged": 100.0 - unflagged}
                for name, unflagged in counts.items()
            }
        }
        table = FakeTableTool({
            "nrows": len(self.NAMES),
            "columns": {"NAME": list(self.NAMES)},
        })
        with mock.patch.object(SELFCAL, "flagdata", lambda **kw: summary), \
                mock.patch.object(SELFCAL, "tb", table), \
                mock.patch("builtins.print"):
            record = manager.measure_antenna_data("target.ms")
        return manager, record

    def test_dead_and_nearly_dead_antennas_are_identified(self):
        # 35% unflagged is healthy for this band; 2% is an antenna that lost
        # its calibration for the main scan; 0% is a quarantined antenna.
        counts = {"C00": 35, "C01": 35, "C02": 2, "C03": 0,
                  "E02": 36, "S01": 34, "W01": 35}       # W06 has no rows
        manager, record = self.census(counts)

        self.assertEqual(record["status"], "measured")
        self.assertEqual(record["unsolvable_names"], ["C02", "C03", "W06"])
        self.assertEqual(manager.unsolvable_antenna_ids, {2, 3, 7})

    def test_an_antenna_with_most_of_its_data_is_kept(self):
        counts = {name: 35 for name in self.NAMES}
        counts["W06"] = 12                        # a third of the median
        manager, record = self.census(counts)

        self.assertEqual(record["unsolvable_names"], [])
        self.assertTrue(record["antennas"]["W06"]["solvable"])

    def test_reference_antennas_are_ranked_by_unflagged_data(self):
        counts = {"C00": 0, "C01": 30, "C02": 2, "C03": 0,
                  "E02": 36, "S01": 34, "W01": 35, "W06": 33}
        manager, _ = self.census(counts)

        refants = manager.phase_reference_antennas("target.ms").split(",")

        # The configured C00 holds no data here and must not be offered.
        self.assertNotIn("C00", refants)
        self.assertEqual(refants[:3], ["E02", "W01", "S01"])
        self.assertNotIn("C02", refants)

    def test_a_usable_configured_reference_stays_first(self):
        counts = {name: 30 for name in self.NAMES}
        counts["E02"] = 36
        manager, _ = self.census(counts)

        refants = manager.phase_reference_antennas("target.ms").split(",")

        self.assertEqual(refants[0], "C00")
        self.assertEqual(refants[1], "E02")

    def test_a_resumed_run_reuses_the_census_it_accepted_tables_under(self):
        counts = {"C00": 35, "C01": 35, "C02": 2, "C03": 0,
                  "E02": 36, "S01": 34, "W01": 35}
        with tempfile.TemporaryDirectory() as temporary:
            first = self.manager()
            first.reportdir = Path(temporary)
            summary = {
                "antenna": {
                    name: {"total": 100.0, "flagged": 100.0 - unflagged}
                    for name, unflagged in counts.items()
                }
            }
            table = FakeTableTool({
                "nrows": len(self.NAMES),
                "columns": {"NAME": list(self.NAMES)},
            })
            with mock.patch.object(SELFCAL, "flagdata", lambda **kw: summary), \
                    mock.patch.object(SELFCAL, "tb", table), \
                    mock.patch("builtins.print"):
                first.ensure_antenna_census("target.ms")
            self.assertTrue((Path(temporary) / "antenna_data.json").is_file())

            # The flags may have changed since (residual flagging); the
            # resumed run must not measure them again.
            resumed = self.manager()
            resumed.reportdir = Path(temporary)

            def must_not_run(**kwargs):
                raise AssertionError("the census was measured again")

            with mock.patch.object(SELFCAL, "flagdata", must_not_run), \
                    mock.patch("builtins.print"):
                record = resumed.ensure_antenna_census("target.ms")

        self.assertEqual(record["unsolvable_names"], ["C02", "C03", "W06"])
        self.assertEqual(resumed.unsolvable_antenna_ids, {2, 3, 7})

    def test_a_failed_summary_excludes_nothing(self):
        manager = self.manager()

        def broken(**kwargs):
            raise RuntimeError("no such table")

        with mock.patch.object(SELFCAL, "flagdata", broken):
            record = manager.measure_antenna_data("target.ms")

        self.assertEqual(record["status"], "unavailable")
        self.assertEqual(manager.unsolvable_antenna_ids, set())


if __name__ == "__main__":
    unittest.main()
