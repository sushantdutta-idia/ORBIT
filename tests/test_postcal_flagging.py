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


def load_calibration_module_without_casa():
    casatasks = types.ModuleType("casatasks")

    def unused_task(**kwargs):
        return None

    for name in [
        "setjy",
        "ft",
        "gaincal",
        "bandpass",
        "fluxscale",
        "applycal",
        "split",
        "flagdata",
        "flagmanager",
        "hanningsmooth",
    ]:
        setattr(casatasks, name, unused_task)

    casatools = types.ModuleType("casatools")

    class UnusedTable:
        pass

    casatools.table = UnusedTable

    path = SOURCE_DIR / "calibration.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_postcal_flagging",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    with mock.patch.dict(
        sys.modules,
        {
            "casatasks": casatasks,
            "casatools": casatools,
        },
    ):
        spec.loader.exec_module(module)
    return module


CALIBRATION = load_calibration_module_without_casa()


class CorrTypes:

    def __init__(self, values):
        self.values = values

    def ravel(self):
        return list(self.values)


class PolarizationTable:

    def __init__(self, values):
        self.values = values
        self.opened = None
        self.closed = False

    def open(self, path):
        self.opened = path

    def getcol(self, name):
        if name != "CORR_TYPE":
            raise AssertionError(name)
        return CorrTypes(self.values)

    def close(self):
        self.closed = True


class PostcalFlaggingTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def manager(self, telescope="meerkat"):
        manager = object.__new__(CALIBRATION.CalibrationManager)
        manager.ms = str(self.root / "input.ms")
        manager.caldir = self.root / "calibration"
        manager.caldir.mkdir()
        manager.telescope = telescope
        manager.band = "MEERKAT_L" if telescope == "meerkat" else "B3"
        manager.task_config = {
            "postcal_flagging": True,
            "postcal_tfcrop_timecutoff": 5.0,
            "postcal_tfcrop_freqcutoff": 4.5,
            "postcal_rflag_timedevscale": 6.0,
            "postcal_rflag_freqdevscale": 6.0,
            "postcal_rflag_winsize": 5,
            "postcal_max_new_fraction": 0.10,
            "postcal_max_total_fraction": 0.70,
            "postcal_max_group_new_fraction": 0.15,
            "postcal_max_group_total_fraction": 0.85,
        }
        manager.get_flux_calibrators = lambda: ["flux"]
        manager.get_flux_calibrator = lambda: "flux"
        manager.get_bandpass_calibrator = lambda: "bandpass"
        manager.get_secondary_flux_calibrator = lambda: None
        manager.get_phase_calibrators = lambda: ["phase"]
        manager.get_phase_calibrator = lambda: "phase"
        manager.get_targets = lambda: ["target"]
        return manager

    def test_ugmrt_defaults_enable_a_conservative_guarded_second_round(self):
        settings = (
            CALIBRATION.CONFIG.TELESCOPE_CALIBRATION_TASK_CONFIG["ugmrt"]
        )

        self.assertTrue(settings["postcal_flagging"])
        self.assertGreaterEqual(settings["postcal_tfcrop_timecutoff"], 4.5)
        self.assertGreaterEqual(settings["postcal_rflag_timedevscale"], 6.0)
        self.assertLessEqual(settings["postcal_max_new_fraction"], 0.10)
        self.assertLessEqual(settings["postcal_max_total_fraction"], 0.70)

    def test_prebandpass_solve_establishes_receptor_coherence_with_g(self):
        manager = self.manager("meerkat")
        manager.task_config.update({
            "prebp_gaintype": "G",
            "gain_minsnr": 3.0,
            "minblperant": 4,
            "gain_solmode": "L1R",
            "parang": False,
        })
        manager.band_config = {"combine_delay": ""}
        manager.get_reference_antenna = lambda: "m000"
        manager.get_antenna_selection = lambda: ""
        manager.caltable_checkpoint_ready = lambda *args: False
        manager.caltable_status = lambda *args: (True, "")
        manager.require_caltable = lambda *args: None
        manager.mark_checkpoint = lambda *args: None
        calls = []

        with mock.patch.object(
            CALIBRATION,
            "gaincal",
            side_effect=lambda **kwargs: calls.append(kwargs),
        ):
            manager.run_prebp()

        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["gaintype"], "G")
        self.assertEqual(calls[0]["calmode"], "p")
        self.assertTrue(calls[0]["gaintable"][0].endswith("delay.cal"))

    def test_parallel_selector_is_absolute_and_excludes_cross_hands(self):
        manager = self.manager()
        tool = PolarizationTable([9, 10, 11, 12])

        with mock.patch.object(CALIBRATION, "table", return_value=tool):
            selector = manager.absolute_parallel_hand_correlation()

        self.assertEqual(selector, "ABS_XX,YY")
        self.assertTrue(tool.opened.endswith("/POLARIZATION"))
        self.assertTrue(tool.closed)

    def test_parallel_selector_fails_closed_on_unknown_correlations(self):
        manager = self.manager()
        tool = PolarizationTable([1, 2, 3, 4])

        with mock.patch.object(CALIBRATION, "table", return_value=tool):
            with self.assertRaisesRegex(RuntimeError, "parallel-hand"):
                manager.absolute_parallel_hand_correlation()

    def test_meerkat_corrected_data_detectors_use_scan_parallel_hands(self):
        manager = self.manager("meerkat")
        manager.absolute_parallel_hand_correlation = lambda: "ABS_XX,YY"
        # Pin the chunk plan so this asserts whole-scan behaviour for data
        # that fit, instead of passing through the no-CASA fallback.
        manager._postcal_flagdata_chunk_plan = {
            "ntime": "scan",
            "status": "scan_fits",
        }
        calls = []
        summaries = [
            {
                "total": 1000,
                "flagged": 100,
                "field": {"target": {"total": 500, "flagged": 50}},
                "spw": {"0": {"total": 1000, "flagged": 100}},
            },
            {
                "total": 1000,
                "flagged": 130,
                "field": {"target": {"total": 500, "flagged": 65}},
                "spw": {"0": {"total": 1000, "flagged": 130}},
            },
        ]

        def fake_flagdata(**kwargs):
            calls.append(kwargs)
            if kwargs["mode"] == "summary":
                return summaries.pop(0)
            return None

        manager_calls = []

        with mock.patch.object(CALIBRATION, "flagdata", fake_flagdata), \
                mock.patch.object(
                    CALIBRATION,
                    "flagmanager",
                    lambda **kwargs: manager_calls.append(kwargs),
                ):
            manager.run_postcal_flagging()

        detectors = {
            call["mode"]: call
            for call in calls
            if call["mode"] in {"tfcrop", "rflag"}
        }
        self.assertEqual(set(detectors), {"tfcrop", "rflag"})
        for call in detectors.values():
            self.assertEqual(call["datacolumn"], "CORRECTED")
            self.assertEqual(call["correlation"], "ABS_XX,YY")
            self.assertEqual(call["ntime"], "scan")
            self.assertFalse(call["combinescans"])
            self.assertFalse(call["extendflags"])

        self.assertEqual([call["mode"] for call in manager_calls], ["save"])
        report = json.loads(
            (manager.caldir / "meerkat_postcal_flagging_report.json")
            .read_text()
        )
        self.assertEqual(report["status"], "completed")
        self.assertEqual(report["correlation"], "ABS_XX,YY")
        self.assertTrue(
            (manager.caldir / "meerkat_postcal_flagging.done").is_file()
        )

    def summaries(self):
        return [
            {
                "total": 1000,
                "flagged": 100,
                "field": {"target": {"total": 500, "flagged": 50}},
                "spw": {"0": {"total": 1000, "flagged": 100}},
            },
            {
                "total": 1000,
                "flagged": 130,
                "field": {"target": {"total": 500, "flagged": 65}},
                "spw": {"0": {"total": 1000, "flagged": 130}},
            },
        ]

    def run_postcal(self, telescope="meerkat", **task_config):
        manager = self.manager(telescope)
        manager.absolute_parallel_hand_correlation = lambda: "ABS_XX,YY"
        manager._postcal_flagdata_chunk_plan = {
            "ntime": "scan",
            "status": "scan_fits",
        }
        manager.task_config.update(task_config)
        calls = []
        summaries = self.summaries()

        def fake_flagdata(**kwargs):
            calls.append(kwargs)
            if kwargs["mode"] == "summary":
                return summaries.pop(0)
            return None

        with mock.patch.object(CALIBRATION, "flagdata", fake_flagdata), \
                mock.patch.object(CALIBRATION, "flagmanager", lambda **k: None):
            manager.run_postcal_flagging()
        return manager, calls

    def test_accepted_flags_are_grown_over_neighbours_within_each_scan(self):
        manager, calls = self.run_postcal()

        modes = [call["mode"] for call in calls]
        # Detection first, then the growth around what it found.
        self.assertEqual(
            modes,
            ["summary", "tfcrop", "rflag", "extend", "summary"],
        )

        extend = calls[modes.index("extend")]
        self.assertEqual(extend["datacolumn"], "CORRECTED")
        self.assertEqual(extend["ntime"], "scan")
        self.assertFalse(extend["combinescans"])
        self.assertEqual(extend["growtime"], 80.0)
        self.assertEqual(extend["growfreq"], 80.0)
        self.assertFalse(extend["growaround"])
        self.assertFalse(extend["flagneartime"])
        self.assertFalse(extend["flagnearfreq"])

        report = json.loads(
            (manager.caldir / "meerkat_postcal_flagging_report.json")
            .read_text()
        )
        recorded = report["parameters"]["extend"]
        self.assertTrue(recorded["applied"])
        self.assertFalse(recorded["combinescans"])
        self.assertEqual(recorded["growtime"], 80.0)

    def test_extension_can_be_turned_off(self):
        manager, calls = self.run_postcal(postcal_extend=False)

        self.assertNotIn("extend", [call["mode"] for call in calls])
        report = json.loads(
            (manager.caldir / "meerkat_postcal_flagging_report.json")
            .read_text()
        )
        self.assertFalse(report["parameters"]["extend"]["applied"])

    def test_bounded_chunk_reaches_detectors_and_report(self):
        manager = self.manager("meerkat")
        manager.absolute_parallel_hand_correlation = lambda: "ABS_XX,YY"
        manager._postcal_flagdata_chunk_plan = {
            "ntime": "304s",
            "status": "bounded",
        }
        summaries = self.summaries()
        calls = []

        def fake_flagdata(**kwargs):
            calls.append(kwargs)
            if kwargs["mode"] == "summary":
                return summaries.pop(0)
            return None

        with mock.patch.object(CALIBRATION, "flagdata", fake_flagdata), \
                mock.patch.object(
                    CALIBRATION, "flagmanager", lambda **kwargs: None,
                ):
            manager.run_postcal_flagging()

        self.assertTrue(calls)
        self.assertTrue(all(call["ntime"] == "304s" for call in calls))
        report = json.loads(
            (manager.caldir / "meerkat_postcal_flagging_report.json")
            .read_text()
        )
        self.assertEqual(report["parameters"]["ntime"], "304s")

    def test_insufficient_memory_stops_before_flags_are_saved(self):
        manager = self.manager("meerkat")
        manager.absolute_parallel_hand_correlation = lambda: "ABS_XX,YY"

        class Metadata:
            def open(self, path):
                pass

            def close(self):
                pass

            def scannumbers(self):
                return [1]

            def timesforscan(self, scan):
                return [8.0 * index for index in range(150)]

            def antennasforscan(self, scan):
                return list(range(64))

            def spwsforscan(self, scan):
                return [0]

            def nchan(self, spw):
                return 32768

        class Polarization:
            def open(self, path):
                pass

            def getcol(self, name):
                return [4]

            def close(self):
                pass

        casatools = types.ModuleType("casatools")
        casatools.msmetadata = Metadata
        limit = {
            "limit_bytes": 1024 ** 3,
            "limiting_source": "slurm_mem_per_node",
            "sources_bytes": {"slurm_mem_per_node": 1024 ** 3},
        }
        calls = []
        manager_calls = []

        with mock.patch.dict(sys.modules, {"casatools": casatools}), \
                mock.patch.object(CALIBRATION, "table", Polarization), \
                mock.patch.object(
                    CALIBRATION.FLAGGING_SELECTORS,
                    "detect_memory_limit_bytes",
                    lambda **kwargs: limit,
                ), \
                mock.patch.object(
                    CALIBRATION, "flagdata",
                    lambda **kwargs: calls.append(kwargs),
                ), \
                mock.patch.object(
                    CALIBRATION, "flagmanager",
                    lambda **kwargs: manager_calls.append(kwargs),
                ):
            with self.assertRaisesRegex(RuntimeError, "No flags were written"):
                manager.run_postcal_flagging()

        self.assertEqual(calls, [])
        self.assertEqual(manager_calls, [])
        self.assertFalse(
            (manager.caldir / "meerkat_postcal_flagging.done").exists()
        )

    def test_postcal_group_guard_restores_flags_and_withholds_marker(self):
        manager = self.manager("ugmrt")
        manager.absolute_parallel_hand_correlation = lambda: "ABS_RR,LL"
        summaries = [
            {
                "total": 1000,
                "flagged": 100,
                "field": {"target": {"total": 100, "flagged": 5}},
            },
            {
                "total": 1000,
                "flagged": 150,
                "field": {"target": {"total": 100, "flagged": 30}},
            },
        ]

        def fake_flagdata(**kwargs):
            if kwargs["mode"] == "summary":
                return summaries.pop(0)
            return None

        manager_calls = []
        with mock.patch.object(CALIBRATION, "flagdata", fake_flagdata), \
                mock.patch.object(
                    CALIBRATION,
                    "flagmanager",
                    lambda **kwargs: manager_calls.append(kwargs),
                ):
            with self.assertRaisesRegex(RuntimeError, "field target"):
                manager.run_postcal_flagging()

        self.assertEqual(
            [call["mode"] for call in manager_calls],
            ["save", "restore"],
        )
        state = json.loads(
            (manager.caldir / "ugmrt_postcal_flagging_state.json")
            .read_text()
        )
        self.assertEqual(state["status"], "failed_restored")
        self.assertFalse(
            (manager.caldir / "ugmrt_postcal_flagging.done").exists()
        )

    def test_ugmrt_second_round_runs_a_complete_fresh_solve(self):
        manager = self.manager("ugmrt")
        operations = []
        manager.checkpoint_exists = lambda name: False
        manager.mark_checkpoint = (
            lambda name: operations.append(f"mark:{name}")
        )
        manager.run_postcal_flagging = lambda: operations.append("postcal")
        manager.prepare_second_round = lambda: operations.append("archive")
        for method in [
            "run_delaycal",
            "run_prebp",
            "run_bandpass",
            "run_phasecal",
            "run_phase_transfer_cal",
            "run_ampcal",
            "run_fluxscale",
            "run_applycal",
        ]:
            setattr(
                manager,
                method,
                lambda method=method: operations.append(method),
            )

        manager.run_second_round()

        self.assertEqual(
            operations,
            [
                "postcal",
                "archive",
                "run_delaycal",
                "run_prebp",
                "run_bandpass",
                "run_phasecal",
                # The scan-averaged transfer table is part of the fresh
                # round-2 solve; reusing the round-1 one would mix rounds.
                "run_phase_transfer_cal",
                "run_ampcal",
                "run_fluxscale",
                "run_applycal",
                "mark:ugmrt_round2_complete",
            ],
        )
        self.assertEqual(manager.calibration_rounds, 2)
        report = json.loads(
            (manager.caldir / "ugmrt_round2_report.json").read_text()
        )
        self.assertEqual(report["status"], "completed")


if __name__ == "__main__":
    unittest.main()
