import importlib.util
import json
import math
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SOURCE_DIR = Path(__file__).resolve().parents[1]


def _no_op_task(*args, **kwargs):
    return None


class _ImportOnlyTable:
    pass


def load_calibration_module():
    casatasks = types.ModuleType("casatasks")
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
        setattr(casatasks, name, _no_op_task)

    casatools = types.ModuleType("casatools")
    casatools.table = _ImportOnlyTable
    module_path = SOURCE_DIR / "calibration.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_calibration_table_qa",
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


CALIBRATION = load_calibration_module()


class FakeTableTool:

    def __init__(self, datasets):
        self.datasets = datasets
        self.current = None

    def open(self, path):
        name = Path(path).name
        self.current = self.datasets.get(name, self.datasets.get("*"))
        if self.current is None:
            raise RuntimeError(f"No mocked table data for {name}")

    def close(self):
        self.current = None

    def nrows(self):
        return self.current.get("nrows", 1)

    def colnames(self):
        return list(self.current.get("columns", {}))

    def getcol(self, name):
        return self.current["columns"][name]


def make_table(path):
    path.mkdir(parents=True, exist_ok=True)
    (path / "table.dat").touch()


class CalibrationTableQATests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.caldir = Path(self.temporary.name)
        self.manager = CALIBRATION.CalibrationManager.__new__(
            CALIBRATION.CalibrationManager
        )
        self.manager.caldir = self.caldir
        self.manager.outdir = self.caldir
        self.manager.task_config = {
            "max_caltable_flagged_fraction": 0.50,
            "max_caltable_group_flagged_fraction": 0.80,
        }
        make_table(self.caldir / "test.cal")

    def tearDown(self):
        self.temporary.cleanup()

    def set_dataset(self, columns, nrows=1):
        CALIBRATION.tb = FakeTableTool({
            "test.cal": {
                "nrows": nrows,
                "columns": columns,
            }
        })

    def qa(self):
        return json.loads(
            (self.caldir / "test.cal.qa.json").read_text()
        )

    def test_valid_cparam_table_writes_metrics_sidecar(self):
        self.set_dataset({
            "CPARAM": [1 + 0j, 0 + 1j, complex(math.nan, 0)],
            "FLAG": [False, False, True],
        })

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertTrue(valid, reason)
        qa = self.qa()
        self.assertEqual(qa["status"], "passed")
        self.assertEqual(qa["parameter_column"], "CPARAM")
        self.assertEqual(qa["total_solutions"], 3)
        self.assertEqual(qa["flagged_solutions"], 1)
        self.assertAlmostEqual(qa["flagged_fraction"], 1 / 3)
        self.assertEqual(qa["nonfinite_usable_solutions"], 0)

    def test_missing_solution_column_is_rejected(self):
        self.set_dataset({"FLAG": [False]})

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertFalse(valid)
        self.assertIn("neither CPARAM nor FPARAM", reason)
        self.assertEqual(self.qa()["status"], "failed")

    def test_missing_flag_column_is_rejected(self):
        self.set_dataset({"CPARAM": [1 + 0j]})

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertFalse(valid)
        self.assertIn("no FLAG column", reason)
        self.assertEqual(self.qa()["status"], "failed")

    def test_all_flagged_table_is_rejected(self):
        self.set_dataset({
            "CPARAM": [1 + 0j, 2 + 0j],
            "FLAG": [True, True],
        })

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertFalse(valid)
        self.assertIn("all calibration solutions are flagged", reason)

    def test_unflagged_nonfinite_solution_is_rejected(self):
        self.set_dataset({
            "CPARAM": [1 + 0j, complex(math.inf, 0)],
            "FLAG": [False, False],
        })

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertFalse(valid)
        self.assertIn("non-finite", reason)
        self.assertEqual(self.qa()["nonfinite_usable_solutions"], 1)

    def test_zero_delay_fparam_table_is_allowed(self):
        self.set_dataset({
            "FPARAM": [0.0, 0.0],
            "FLAG": [False, False],
        })

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertTrue(valid, reason)

    def test_partly_zero_fparam_table_is_allowed(self):
        self.set_dataset({
            "FPARAM": [0.0, 1.25],
            "FLAG": [False, False],
        })

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertTrue(valid, reason)
        self.assertEqual(self.qa()["zero_usable_solutions"], 0)

    def test_any_unflagged_zero_cparam_is_rejected(self):
        self.set_dataset({
            "CPARAM": [1 + 0j, 0 + 0j],
            "FLAG": [False, False],
        })

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertFalse(valid)
        self.assertIn("CPARAM", reason)
        self.assertEqual(self.qa()["zero_usable_solutions"], 1)

    def test_configured_flagged_fraction_limit_is_enforced(self):
        self.set_dataset({
            "CPARAM": [1 + 0j, 1 + 0j, 1 + 0j],
            "FLAG": [False, True, True],
        })

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertFalse(valid)
        self.assertIn("exceeds the configured maximum", reason)
        self.assertAlmostEqual(self.qa()["flagged_fraction"], 2 / 3)

    def test_group_with_no_usable_solution_is_rejected(self):
        self.set_dataset({
            "CPARAM": [[1 + 0j, 1 + 0j], [1 + 0j, 1 + 0j]],
            "FLAG": [[False, True], [False, True]],
            "ANTENNA1": [10, 11],
        }, nrows=2)

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertFalse(valid)
        self.assertIn("ANTENNA1=11", reason)
        self.assertIn("no usable solution", reason)
        antenna_qa = self.qa()["group_metrics"]["ANTENNA1"]
        self.assertEqual(antenna_qa["status"], "failed")
        failed = [
            group for group in antenna_qa["groups"]
            if group["value"] == 11
        ][0]
        self.assertEqual(failed["flagged_fraction"], 1.0)
        self.assertEqual(failed["usable_fraction"], 0.0)
        self.assertEqual(failed["usable_solutions"], 0)

    def test_placeholder_for_antenna_absent_from_main_is_excluded(self):
        self.set_dataset({
            "FPARAM": [[1.0, 0.0], [1.0, 0.0]],
            "FLAG": [[False, True], [False, True]],
            "ANTENNA1": [10, 11],
        }, nrows=2)
        self.manager.ms = self.caldir / "input.ms"
        self.manager.active_measurement_set_antenna_ids = lambda: {10}

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertTrue(valid, reason)
        qa = self.qa()
        self.assertEqual(qa["excluded_nonparticipating_antenna_ids"], [11])
        self.assertEqual(qa["ignored_group_row_indices"], [1])
        antenna_qa = qa["group_metrics"]["ANTENNA1"]
        self.assertEqual(antenna_qa["status"], "passed")
        self.assertEqual(antenna_qa["group_count"], 1)
        self.assertEqual(antenna_qa["groups"][0]["value"], 10)

    def test_placeholder_for_deliberately_excluded_antenna_is_excluded(self):
        self.set_dataset({
            "FPARAM": [[1.0, 0.0], [1.0, 0.0]],
            "FLAG": [[False, True], [False, True]],
            "ANTENNA1": [10, 11],
        }, nrows=2)
        self.manager.ms = self.caldir / "input.ms"
        self.manager.active_measurement_set_antenna_ids = lambda: {10, 11}
        self.manager.flagging = {"bad_antennas": ["C11"]}
        inspection = self.caldir / "inspection"
        inspection.mkdir()
        (inspection / "antennas.json").write_text(json.dumps([
            {"id": 10, "name": "C10"},
            {"id": 11, "name": "C11"},
        ]))

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertTrue(valid, reason)
        qa = self.qa()
        self.assertEqual(qa["excluded_solver_antenna_ids"], [11])
        self.assertEqual(qa["ignored_group_row_indices"], [1])
        self.assertEqual(
            qa["group_metrics"]["ANTENNA1"]["groups"][0]["value"],
            10,
        )

    def test_legacy_name_only_inventory_maps_ids_by_row_order(self):
        self.manager.outdir = self.caldir
        inspection = self.caldir / "inspection"
        inspection.mkdir()
        (inspection / "antennas.json").write_text(json.dumps([
            {"name": "C00"},
            {"name": "C01"},
            {"name": "C02"},
            {"name": "C03"},
        ]))

        self.assertEqual(
            self.manager.inspection_antenna_id_name_map(),
            {0: "C00", 1: "C01", 2: "C02", 3: "C03"},
        )

    def test_group_fraction_limit_catches_locally_bad_antenna(self):
        self.set_dataset({
            "CPARAM": [[1 + 0j, 1 + 0j] for _ in range(10)],
            "FLAG": (
                [[False, True] for _ in range(9)] +
                [[False, False]]
            ),
            "ANTENNA1": [10, 11],
        }, nrows=2)

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertFalse(valid)
        self.assertIn("ANTENNA1=11", reason)
        self.assertIn("0.900000", reason)
        qa = self.qa()
        self.assertAlmostEqual(qa["flagged_fraction"], 0.45)
        self.assertEqual(qa["group_qa_status"], "failed")

    def test_single_unsolved_delay_antenna_is_bounded_and_recorded(self):
        manager = self.manager
        manager.input_ms = str(self.caldir / "input.ms")
        manager.ms = manager.input_ms
        manager.outdir = self.caldir
        manager.caldir = self.caldir / "calibration"
        manager.caldir.mkdir()
        manager.flagging = {"bad_antennas": []}
        manager.delay_quarantine = {
            "schema_version": 1,
            "measurement_set": str(Path(manager.input_ms).resolve()),
            "antennas": [],
        }
        manager.task_config.update({
            "auto_quarantine_failed_delay_antennas": True,
            "max_auto_quarantine_delay_antennas": 1,
            "max_auto_quarantine_delay_fraction": 0.10,
            "min_antennas_after_delay_quarantine": 8,
        })
        (self.caldir / "inspection").mkdir()
        (self.caldir / "inspection" / "antennas.json").write_text(
            json.dumps([
                {"id": value, "name": f"C{value:02d}"}
                for value in range(30)
            ])
        )
        manager.active_measurement_set_antenna_ids = lambda: set(range(30))
        manager.get_reference_antenna = lambda: "C00"

        names = manager.quarantine_failed_delay_antennas([3])

        self.assertEqual(names, ["C03"])
        self.assertEqual(manager.get_antenna_selection(), "!C03")
        record = json.loads(manager.delay_quarantine_path().read_text())
        self.assertEqual(
            [(item["id"], item["name"]) for item in record["antennas"]],
            [(3, "C03")],
        )
        self.assertEqual(record["antennas"][0]["reason"],
                         "no usable delay solution")
        self.assertTrue(
            record["antennas"][0]["source_table"].endswith("delay.cal")
        )

    def test_multiple_unsolved_delay_antennas_are_not_auto_quarantined(self):
        manager = self.manager
        manager.input_ms = str(self.caldir / "input.ms")
        manager.ms = manager.input_ms
        manager.outdir = self.caldir
        manager.flagging = {"bad_antennas": []}
        manager.delay_quarantine = {
            "schema_version": 1,
            "measurement_set": str(Path(manager.input_ms).resolve()),
            "antennas": [],
        }
        manager.task_config.update({
            "auto_quarantine_failed_delay_antennas": True,
            "max_auto_quarantine_delay_antennas": 1,
            "max_auto_quarantine_delay_fraction": 0.10,
            "min_antennas_after_delay_quarantine": 8,
        })
        manager.active_measurement_set_antenna_ids = lambda: set(range(30))

        self.assertEqual(
            manager.quarantine_failed_delay_antennas([3, 4]),
            [],
        )
        self.assertEqual(manager.flagging["bad_antennas"], [])

    def test_bandpass_uses_only_supported_reference_controls(self):
        manager = self.manager
        manager.ms = str(self.caldir / "input.ms")
        manager.band_config = {"combine_delay": "scan"}
        manager.flagging = {"bad_antennas": ["C03"]}
        manager.task_config.update({
            "bandpass_fillgaps": 8,
            "bandpass_solint": "inf",
            "bandpass_solnorm": True,
            "gain_minsnr": 3.0,
            "minblperant": 4,
            "parang": False,
        })
        manager.caltable_checkpoint_ready = lambda *args: False
        manager.get_bandpass_calibrator = lambda: "3C147"
        manager.get_reference_antenna = lambda: "C00"
        manager.require_caltable = lambda filename: None
        # run_bandpass reads the table's status before requiring it.
        manager.caltable_status = lambda filename: (True, "")
        marked = []
        manager.mark_checkpoint = lambda stage: marked.append(stage)
        calls = []

        with mock.patch.object(
            CALIBRATION,
            "bandpass",
            lambda **kwargs: calls.append(kwargs),
        ):
            manager.run_bandpass()

        self.assertEqual(len(calls), 1)
        self.assertNotIn("refantmode", calls[0])
        self.assertEqual(calls[0]["refant"], "C00")
        self.assertEqual(calls[0]["minblperant"], 4)
        self.assertEqual(calls[0]["antenna"], "!C03")
        self.assertEqual(marked, ["bandpass"])

    def test_all_available_identifier_columns_are_reported(self):
        self.set_dataset({
            "CPARAM": [[1 + 0j, 1 + 0j] for _ in range(4)],
            "FLAG": [
                [False, False],
                [True, False],
                [False, True],
                [False, False],
            ],
            "ANTENNA1": [10, 11],
            "SPECTRAL_WINDOW_ID": [2, 2],
            "FIELD_ID": [7, 7],
            "OBSERVATION_ID": [0, 1],
        }, nrows=2)

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertTrue(valid, reason)
        qa = self.qa()
        self.assertEqual(qa["group_qa_status"], "passed")
        self.assertEqual(
            set(qa["group_metrics"]),
            {
                "ANTENNA1",
                "SPECTRAL_WINDOW_ID",
                "FIELD_ID",
                "OBSERVATION_ID",
            },
        )
        self.assertEqual(
            qa["group_metrics"]["ANTENNA1"]["group_count"],
            2,
        )
        self.assertAlmostEqual(
            qa["group_metrics"]["ANTENNA1"]["groups"][0][
                "usable_fraction"
            ],
            0.75,
        )
        for status in qa["group_column_status"].values():
            self.assertEqual(status["status"], "evaluated")

    def test_misaligned_optional_group_column_is_recorded_and_skipped(self):
        self.set_dataset({
            "CPARAM": [[1 + 0j, 1 + 0j]],
            "FLAG": [[False, False]],
            "ANTENNA1": [10],
        }, nrows=2)

        valid, reason = self.manager.caltable_status("test.cal")

        self.assertTrue(valid, reason)
        qa = self.qa()
        self.assertEqual(
            qa["group_column_status"]["ANTENNA1"]["status"],
            "skipped",
        )
        self.assertEqual(qa["group_qa_status"], "not_evaluated")

    def test_multiple_phase_calibrators_require_an_explicit_target_map(self):
        manager = self.manager
        manager.report = {
            "fields": {
                "flux_calibrators": ["flux_a"],
                "bandpass_calibrators": ["flux_a"],
                "phase_calibrators": ["phase_a", "phase_b"],
                "targets": ["target_a", "target_b"],
            }
        }
        manager.field_overrides = {
            "flux_calibrator": ["flux_a"],
            "bandpass_calibrator": ["flux_a"],
            "secondary_flux_calibrator": [],
            "phase_calibrator": ["phase_a", "phase_b"],
            "target_phase_map": None,
            "target_fields": ["target_a", "target_b"],
            "reference_antenna": None,
        }

        with self.assertRaisesRegex(
            RuntimeError,
            "Multiple phase calibrators are available",
        ):
            manager.get_target_phase_map()

    def test_explicit_primary_flux_breaks_vetted_bandpass_tie(self):
        manager = self.manager
        manager.report = {
            "fields": {
                "flux_calibrators": ["flux_a", "flux_b"],
                "bandpass_calibrators": ["flux_a", "flux_b"],
                "phase_calibrators": ["phase"],
                "targets": ["target"],
            }
        }
        manager.field_overrides = {
            "flux_calibrator": ["flux_a"],
            "bandpass_calibrator": [],
            "secondary_flux_calibrator": ["flux_b"],
            "phase_calibrator": ["phase"],
            "target_phase_map": None,
            "target_fields": ["target"],
        }

        self.assertEqual(manager.get_bandpass_calibrator(), "flux_a")

    def test_inspection_strategy_breaks_inferred_bandpass_tie(self):
        manager = self.manager
        manager.report = {
            "fields": {
                "flux_calibrators": ["3C147", "3C48"],
                "bandpass_calibrators": ["3C147", "3C48"],
                "phase_calibrators": ["phase"],
                "targets": ["target"],
            },
            "calibration_strategy": {
                "primary_flux_calibrator": "3C147",
                "bandpass_calibrator": "3C147",
                "secondary_flux_calibrator": "3C48",
                "phase_calibrator": "phase",
                "target_fields": ["target"],
            },
        }
        manager.field_overrides = {
            "flux_calibrator": [],
            "bandpass_calibrator": [],
            "secondary_flux_calibrator": [],
            "phase_calibrator": [],
            "target_phase_map": None,
            "target_fields": [],
        }

        self.assertEqual(manager.get_bandpass_calibrator(), "3C147")

    def band4_report(self, strategy_primary="3C147"):
        return {
            "fields": {
                "flux_calibrators": ["3C147", "3C48"],
                "bandpass_calibrators": ["3C147", "3C48"],
                "phase_calibrators": ["0323+055"],
                "targets": ["target"],
            },
            "calibration_strategy": {
                "primary_flux_calibrator": strategy_primary,
                "bandpass_calibrator": "3C147",
                "secondary_flux_calibrator": "3C48",
                "phase_calibrator": "0323+055",
                "target_fields": ["target"],
            },
        }

    def band4_overrides(self):
        return {
            "flux_calibrator": [],
            "bandpass_calibrator": [],
            "secondary_flux_calibrator": ["0323+055"],
            "phase_calibrator": [],
            "target_phase_map": None,
            "target_fields": [],
        }

    def test_named_secondary_takes_inspection_primary_of_two_standards(self):
        manager = self.manager
        manager.report = self.band4_report()
        manager.field_overrides = self.band4_overrides()

        self.assertEqual(manager.get_flux_calibrators(), ["3C147", "0323+055"])
        self.assertEqual(manager.get_secondary_flux_calibrator(), "0323+055")
        self.assertEqual(manager.get_bandpass_calibrator(), "3C147")

    def test_named_primary_and_bandpass_override_inspection_strategy(self):
        manager = self.manager
        manager.report = self.band4_report()
        manager.field_overrides = dict(
            self.band4_overrides(),
            flux_calibrator=["3C48"],
            bandpass_calibrator=["3C48"],
        )

        self.assertEqual(manager.get_flux_calibrators(), ["3C48", "0323+055"])
        self.assertEqual(manager.get_secondary_flux_calibrator(), "0323+055")
        self.assertEqual(manager.get_bandpass_calibrator(), "3C48")

    def test_named_secondary_with_unvetted_strategy_primary_still_fails(self):
        manager = self.manager
        manager.report = self.band4_report(strategy_primary="not_a_standard")
        manager.field_overrides = self.band4_overrides()

        with self.assertRaisesRegex(
            RuntimeError,
            "Multiple primary flux calibrators",
        ):
            manager.get_flux_calibrators()

    def test_named_secondary_with_one_standard_is_unchanged(self):
        manager = self.manager
        manager.report = {
            "fields": {
                "flux_calibrators": ["3C48"],
                "bandpass_calibrators": ["3C48"],
                "phase_calibrators": ["0323+055"],
                "targets": ["target"],
            },
            "calibration_strategy": {"primary_flux_calibrator": "3C48"},
        }
        manager.field_overrides = self.band4_overrides()

        self.assertEqual(manager.get_flux_calibrators(), ["3C48", "0323+055"])

    def test_unvetted_inspection_bandpass_strategy_is_not_trusted(self):
        manager = self.manager
        manager.report = {
            "fields": {
                "flux_calibrators": ["flux"],
                "bandpass_calibrators": ["bp_a", "bp_b"],
                "phase_calibrators": ["phase"],
                "targets": ["target"],
            },
            "calibration_strategy": {
                "bandpass_calibrator": "not_a_candidate",
            },
        }
        manager.field_overrides = {
            "flux_calibrator": [],
            "bandpass_calibrator": [],
            "secondary_flux_calibrator": [],
            "phase_calibrator": [],
            "target_phase_map": None,
            "target_fields": [],
        }

        with self.assertRaisesRegex(
            RuntimeError,
            "Multiple bandpass calibrators",
        ):
            manager.get_bandpass_calibrator()

    def test_bandpass_tie_still_fails_when_primary_is_not_candidate(self):
        manager = self.manager
        manager.report = {
            "fields": {
                "flux_calibrators": ["flux"],
                "bandpass_calibrators": ["bp_a", "bp_b"],
                "phase_calibrators": ["phase"],
                "targets": ["target"],
            }
        }
        manager.field_overrides = {
            "flux_calibrator": ["flux"],
            "bandpass_calibrator": [],
            "secondary_flux_calibrator": [],
            "phase_calibrator": ["phase"],
            "target_phase_map": None,
            "target_fields": ["target"],
        }

        with self.assertRaisesRegex(
            RuntimeError,
            "Multiple bandpass calibrators",
        ):
            manager.get_bandpass_calibrator()

    def test_explicit_target_phase_map_is_validated_and_target_ordered(self):
        manager = self.manager
        manager.report = {
            "fields": {
                "flux_calibrators": ["flux_a"],
                "bandpass_calibrators": ["flux_a"],
                "phase_calibrators": ["phase_a", "phase_b"],
                "targets": ["target_a", "target_b"],
            }
        }
        manager.field_overrides = {
            "flux_calibrator": ["flux_a"],
            "bandpass_calibrator": ["flux_a"],
            "secondary_flux_calibrator": [],
            "phase_calibrator": ["phase_a", "phase_b"],
            "target_phase_map": {
                "target_b": "phase_b",
                "target_a": "phase_a",
            },
            "target_fields": ["target_a", "target_b"],
            "reference_antenna": None,
        }

        self.assertEqual(
            manager.get_target_phase_map(),
            {
                "target_a": "phase_a",
                "target_b": "phase_b",
            },
        )

    def test_applycal_uses_all_flux_fields_and_explicit_target_mapping(self):
        manager = self.manager
        manager.ms = str(self.caldir / "input.ms")
        make_table(Path(manager.ms))
        manager.field_overrides = {
            "flux_calibrator": ["flux_a", "flux_b"],
            "bandpass_calibrator": ["flux_a"],
            "secondary_flux_calibrator": [],
            "phase_calibrator": ["phase_a", "phase_b"],
            "target_phase_map": {
                "target_a": "phase_a",
                "target_b": "phase_b",
            },
            "target_fields": ["target_a", "target_b"],
            "reference_antenna": None,
        }
        manager.report = {
            "fields": {
                "flux_calibrators": ["flux_a", "flux_b"],
                "bandpass_calibrators": ["flux_a"],
                "phase_calibrators": ["phase_a", "phase_b"],
                "targets": ["target_a", "target_b"],
            }
        }
        manager.task_config.update({
            "applymode": "calflagstrict",
            "calwt": True,
            "parang": False,
        })
        manager.checkpoint_exists = lambda stage: False
        marked = []
        manager.mark_checkpoint = lambda stage: marked.append(stage)
        calls = []

        with mock.patch.object(
            CALIBRATION,
            "applycal",
            lambda **kwargs: calls.append(kwargs),
        ):
            manager.run_applycal()

        self.assertEqual(
            [call["field"] for call in calls],
            [
                "flux_a",
                "flux_b",
                "phase_a",
                "phase_b",
                "target_a",
                "target_b",
            ],
        )
        self.assertEqual(
            calls[-2]["gainfield"][-1],
            "phase_a",
        )
        self.assertEqual(
            calls[-1]["gainfield"][-1],
            "phase_b",
        )
        self.assertEqual(marked, ["applycal"])

    def test_write_report_embeds_all_table_qa_sidecars(self):
        filenames = [
            "delay.cal",
            "prebp.cal",
            "bandpass.cal",
            "phase.cal",
            "phase_scan.cal",
            "amp.cal",
            "flux.cal",
        ]
        for filename in filenames:
            make_table(self.caldir / filename)
        CALIBRATION.tb = FakeTableTool({
            "*": {
                "nrows": 1,
                "columns": {
                    "CPARAM": [1 + 0j, 1 + 0j],
                    "FLAG": [False, False],
                },
            }
        })

        manager = self.manager
        input_ms = self.caldir / "input.ms"
        make_table(input_ms)
        manager.report = {"measurement_set": str(input_ms)}
        manager.ms = str(input_ms)
        manager.telescope = "ugmrt"
        manager.band = "B3"
        manager.flux_model = None
        manager.flux_model_type = "casa-image"
        manager.provenance_fingerprint = "f" * 64
        manager.band_config = {"hanningsmooth": True}
        manager.flagging = {"bad_antennas": []}
        manager.task_config.update({
            "flux_standard": "Perley-Butler 2017",
            "postcal_flagging": False,
            "bandpass_solnorm": True,
            "prebp_gaintype": "G",
            "applymode": "calflag",
            "calwt": True,
            "parang": False,
        })
        manager.get_flux_calibrators = lambda: ["flux", "flux_2"]
        manager.get_flux_calibrator = lambda: "flux"
        manager.get_secondary_flux_calibrator = lambda: "flux_2"
        manager.get_bandpass_calibrator = lambda: "flux"
        manager.get_phase_calibrators = lambda: ["phase_a", "phase_b"]
        manager.get_phase_calibrator = lambda: "phase"
        manager.get_target_phase_map = lambda: {"target": "phase_a"}
        manager.get_targets = lambda: ["target"]
        manager.get_reference_antenna = lambda: "C00"
        manager.checkpoint_exists = lambda stage: True
        manager.write_atomic_json(
            manager.flux_model_provenance_path(),
            {
                "schema_version": 1,
                "status": "completed",
                "identity": {"mode": "builtin-setjy"},
                "model_data_qa": {"flux": {"status": "passed"}},
            },
        )
        manager.write_atomic_json(
            self.caldir / "split_targets.provenance.json",
            {"schema_version": 1, "targets": {}},
        )

        manager.write_report()

        report = json.loads(
            (self.caldir / "calibration_report.json").read_text()
        )
        self.assertEqual(
            set(report["calibration_table_qa"]),
            {
                "delay",
                "prebp",
                "bandpass",
                "phase",
                "phase_scan",
                "amp",
                "flux",
            },
        )
        for table_qa in report["calibration_table_qa"].values():
            self.assertEqual(table_qa["status"], "passed")
            self.assertTrue(Path(table_qa["sidecar"]).is_file())
            self.assertIn("group_metrics", table_qa)
        self.assertEqual(
            report["calibration_parameters"][
                "max_caltable_group_flagged_fraction"
            ],
            0.80,
        )
        self.assertEqual(
            report["measurement_set_stat_sha256"],
            CALIBRATION.IMAGING_SAFETY.path_stat_fingerprint(input_ms),
        )
        self.assertEqual(
            report["split_target_identities"],
            {"schema_version": 1, "targets": {}},
        )
        self.assertEqual(report["provenance_fingerprint"], "f" * 64)
        self.assertEqual(
            report["flux_model_provenance"]["status"],
            "completed",
        )
        self.assertEqual(report["flux_calibrators"], ["flux", "flux_2"])
        self.assertEqual(
            report["phase_calibrators"],
            ["phase_a", "phase_b"],
        )
        self.assertEqual(
            report["target_phase_map"],
            {"target": "phase_a"},
        )


if __name__ == "__main__":
    unittest.main()
