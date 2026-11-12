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


def _unused_task(*args, **kwargs):
    return None


class _ImportTable:
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
        setattr(casatasks, name, _unused_task)

    casatools = types.ModuleType("casatools")
    casatools.table = _ImportTable
    path = SOURCE_DIR / "calibration.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_flux_model_safety",
        path,
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


def setjy_result(i0=5.0, i1=4.5):
    return {
        "0": {
            "fieldName": "primary",
            "0": {
                "fluxd": [i0, 0.0, 0.0, 0.0],
                "numSol": 1,
            },
            "1": {
                "fluxd": [i1, 0.0, 0.0, 0.0],
                "numSol": 1,
            },
        }
    }


class ModelSelection:

    def __init__(self, values, rows=1):
        self.values = values
        self.rows = rows
        self.closed = False

    def nrows(self):
        return self.rows

    def getcol(self, name, startrow=0, nrow=-1):
        # Mirrors casatools: an explicit row range reads only those rows.
        # Here each list entry stands for one row.
        if name != "MODEL_DATA":
            raise AssertionError(name)
        self.reads = getattr(self, "reads", []) + [(startrow, nrow)]
        if nrow < 0:
            return self.values[startrow:]
        return self.values[startrow:startrow + nrow]

    def close(self):
        self.closed = True


class ModelTable:

    def __init__(self, values, rows=1):
        self.values = values
        self.rows = rows
        self.expression = None
        self.closed = False

    def open(self, path):
        self.path = path

    def colnames(self):
        return ["DATA", "MODEL_DATA", "FIELD_ID"]

    def query(self, expression, columns=None):
        self.expression = expression
        self.columns = columns
        return ModelSelection(self.values, self.rows)

    def close(self):
        self.closed = True


class ShapedModelSelection:
    """Generates (ncorr, nchan, nrow) MODEL_DATA blocks on demand."""

    def __init__(self, rows, nchan=64, ncorr=4, nan_row=None, zero=False):
        self.rows = rows
        self.nchan = nchan
        self.ncorr = ncorr
        self.nan_row = nan_row
        self.zero = zero
        self.reads = []

    def nrows(self):
        return self.rows

    def getcol(self, name, startrow=0, nrow=-1):
        import numpy as np

        if nrow < 0:
            nrow = self.rows - startrow
        self.reads.append((startrow, nrow))
        value = 0j if self.zero else 1.5 + 0j
        block = np.full((self.ncorr, self.nchan, nrow), value, dtype=complex)
        if self.nan_row is not None and startrow <= self.nan_row < startrow + nrow:
            block[0, 0, self.nan_row - startrow] = complex(float("nan"), 0.0)
        return block

    def close(self):
        pass


class ShapedModelTable(ModelTable):

    def __init__(self, selection):
        super().__init__(values=[], rows=selection.rows)
        self.selection = selection

    def query(self, expression, columns=None):
        self.expression = expression
        return self.selection


class FluxModelSafetyTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.outdir = self.root / "Products"
        self.caldir = self.outdir / "calibration"
        self.caldir.mkdir(parents=True)
        self.original_ft = CALIBRATION.ft
        self.original_setjy = CALIBRATION.setjy
        self.original_tb = CALIBRATION.tb

    def tearDown(self):
        CALIBRATION.ft = self.original_ft
        CALIBRATION.setjy = self.original_setjy
        CALIBRATION.tb = self.original_tb
        self.temporary.cleanup()

    def model_path(self, name="primary.image"):
        path = self.root / name
        path.mkdir()
        (path / "table.dat").write_bytes(b"model-v1")
        return path

    def manager(self,
                flux_model=None,
                model_type="casa-image",
                bandpass="primary"):
        manager = CALIBRATION.CalibrationManager.__new__(
            CALIBRATION.CalibrationManager
        )
        manager.ms = str(self.root / "input.ms")
        manager.input_ms = manager.ms
        manager.outdir = self.outdir
        manager.caldir = self.caldir
        manager.plotdir = self.caldir / "plots"
        manager.splitdir = self.caldir / "calibrated"
        manager.splitdir.mkdir(exist_ok=True)
        manager.telescope = "ugmrt"
        manager.band = "B3"
        manager.band_config = {"hanningsmooth": False}
        manager.report = {
            "spws": {
                "0": {"id": 0},
                "1": {"id": 1},
            },
            "nspws": 2,
        }
        manager.task_config = {
            "flux_standard": "Perley-Butler 2017",
        }
        manager.flux_model = str(flux_model) if flux_model else None
        manager.flux_model_type = model_type
        manager.flux_standard_override = None
        manager.provenance_fingerprint = "a" * 64
        manager.recipe_fingerprint = "r" * 64
        manager.get_flux_calibrators = lambda: ["primary"]
        manager.get_flux_calibrator = lambda: "primary"
        manager.get_bandpass_calibrator = lambda: bandpass
        manager.validate_selected_model_data = lambda field: {
            "status": "passed",
            "field": field,
            "field_ids": [0],
            "selected_rows": 2,
            "total_values": 8,
            "finite_values": 8,
            "nonfinite_values": 0,
            "nonzero_values": 4,
            "reason": "",
        }
        return manager

    def test_external_image_is_inserted_with_ft_not_manual_setjy(self):
        model = self.model_path()
        manager = self.manager(model)
        ft_calls = []
        setjy_calls = []
        CALIBRATION.ft = lambda **kwargs: ft_calls.append(kwargs)
        CALIBRATION.setjy = lambda **kwargs: setjy_calls.append(kwargs)

        manager.run_setjy()

        self.assertEqual(len(ft_calls), 1)
        self.assertEqual(setjy_calls, [])
        self.assertEqual(ft_calls[0]["model"], str(model.resolve()))
        self.assertNotIn("complist", ft_calls[0])
        self.assertTrue(ft_calls[0]["usescratch"])
        self.assertFalse(ft_calls[0]["incremental"])
        self.assertTrue((self.caldir / "setjy.done").is_file())
        provenance = json.loads(
            (self.caldir / "flux_model_provenance.json").read_text()
        )
        self.assertEqual(provenance["status"], "completed")
        self.assertEqual(provenance["identity"]["mode"], "external-ft")
        self.assertEqual(
            provenance["identity"]["flux_model_sha256"],
            CALIBRATION.sha256_file_or_directory(model),
        )

    def test_component_list_uses_ft_complist_keyword(self):
        model = self.model_path("primary.cl")
        manager = self.manager(model, model_type="component-list")
        calls = []
        CALIBRATION.ft = lambda **kwargs: calls.append(kwargs)

        manager.run_setjy()

        self.assertEqual(calls[0]["complist"], str(model.resolve()))
        self.assertNotIn("model", calls[0])

    def test_builtin_setjy_result_is_captured_and_uses_scratch_model(self):
        manager = self.manager()
        calls = []

        def fake_setjy(**kwargs):
            calls.append(kwargs)
            return setjy_result()

        CALIBRATION.setjy = fake_setjy

        manager.run_setjy()

        self.assertEqual(len(calls), 1)
        self.assertTrue(calls[0]["usescratch"])
        self.assertTrue(calls[0]["scalebychan"])
        provenance = json.loads(
            (self.caldir / "flux_model_provenance.json").read_text()
        )
        flux_qa = provenance["applications"][0]["flux_qa"]
        self.assertEqual(flux_qa["returned_spws"], ["0", "1"])
        self.assertEqual(flux_qa["records"][0]["flux_density_jy"][0], 5.0)

    def test_all_one_jy_setjy_fallback_fails_without_checkpoint(self):
        manager = self.manager()
        CALIBRATION.setjy = lambda **kwargs: setjy_result(1.0, 1.0)

        with self.assertRaisesRegex(RuntimeError, "unsafe fallback"):
            manager.run_setjy()

        self.assertFalse((self.caldir / "setjy.done").exists())
        provenance = json.loads(
            (self.caldir / "flux_model_provenance.json").read_text()
        )
        self.assertEqual(provenance["status"], "failed")

    def test_distinct_bandpass_calibrator_must_have_valid_model(self):
        manager = self.manager(bandpass="bandpass")

        def fake_setjy(**kwargs):
            if kwargs["field"] == "bandpass":
                return setjy_result(1.0, 1.0)
            return setjy_result()

        CALIBRATION.setjy = fake_setjy

        with self.assertRaisesRegex(
            RuntimeError,
            "distinct bandpass calibrator.*modeled primary",
        ):
            manager.run_setjy()

        self.assertFalse((self.caldir / "setjy.done").exists())

    def test_distinct_bandpass_builtin_model_is_validated_and_recorded(self):
        manager = self.manager(bandpass="bandpass")
        calls = []

        def fake_setjy(**kwargs):
            calls.append(kwargs["field"])
            return (
                setjy_result(6.0, 5.5)
                if kwargs["field"] == "bandpass"
                else setjy_result()
            )

        CALIBRATION.setjy = fake_setjy

        manager.run_setjy()

        self.assertEqual(calls, ["primary", "bandpass"])
        provenance = json.loads(
            (self.caldir / "flux_model_provenance.json").read_text()
        )
        self.assertEqual(len(provenance["applications"]), 2)
        self.assertEqual(
            provenance["applications"][1]["role"],
            "distinct-bandpass-calibrator",
        )

    def test_external_model_content_hash_invalidates_checkpoint_reuse(self):
        model = self.model_path()
        manager = self.manager(model)
        calls = []
        CALIBRATION.ft = lambda **kwargs: calls.append(kwargs)

        manager.run_setjy()
        manager.run_setjy()
        self.assertEqual(len(calls), 1)

        # A code-only implementation fingerprint change must not rerun an
        # already validated model application.
        manager.provenance_fingerprint = "b" * 64
        manager.run_setjy()
        self.assertEqual(len(calls), 1)

        # A changed science-recipe fingerprint invalidates the checkpoint.
        manager.recipe_fingerprint = "s" * 64
        manager.run_setjy()
        self.assertEqual(len(calls), 2)

        (model / "table.dat").write_bytes(b"model-v2")
        manager.run_setjy()
        self.assertEqual(len(calls), 3)
        provenance = json.loads(
            (self.caldir / "flux_model_provenance.json").read_text()
        )
        self.assertEqual(
            provenance["identity"]["flux_model_sha256"],
            CALIBRATION.sha256_file_or_directory(model),
        )

    def test_setjy_validator_requires_every_spw_and_rejects_nonfinite(self):
        incomplete = {"0": {"0": {"fluxd": [5.0, 0, 0, 0]}}}
        with self.assertRaisesRegex(RuntimeError, "SPW.*1"):
            CALIBRATION.validate_setjy_flux_result(
                incomplete,
                expected_spws=[0, 1],
            )

        nonfinite = setjy_result(math.nan, 4.0)
        with self.assertRaisesRegex(RuntimeError, "non-finite"):
            CALIBRATION.validate_setjy_flux_result(
                nonfinite,
                expected_spws=[0, 1],
            )

    def test_model_data_validation_requires_finite_nonzero_selected_field(self):
        manager = self.manager()
        fields = self.outdir / "inspection" / "fields.json"
        fields.parent.mkdir(parents=True)
        fields.write_text(json.dumps([{"id": 3, "name": "primary"}]))

        CALIBRATION.tb = ModelTable([0 + 0j, 2 + 0j], rows=2)
        qa = CALIBRATION.CalibrationManager.validate_selected_model_data(
            manager,
            "primary",
        )
        self.assertEqual(qa["status"], "passed")
        self.assertEqual(qa["field_ids"], [3])
        self.assertIn("3", CALIBRATION.tb.expression)

        CALIBRATION.tb = ModelTable([0 + 0j, 0 + 0j], rows=2)
        with self.assertRaisesRegex(RuntimeError, "entirely zero"):
            CALIBRATION.CalibrationManager.validate_selected_model_data(
                manager,
                "primary",
            )


    def shaped_validation(self, selection, max_values, block_values):
        manager = self.manager()
        fields = self.outdir / "inspection" / "fields.json"
        fields.parent.mkdir(parents=True, exist_ok=True)
        fields.write_text(json.dumps([{"id": 3, "name": "primary"}]))
        CALIBRATION.tb = ShapedModelTable(selection)
        with mock.patch.object(
            CALIBRATION, "MODEL_DATA_CHECK_MAX_VALUES", max_values,
        ), mock.patch.object(
            CALIBRATION, "MODEL_DATA_CHECK_BLOCK_VALUES", block_values,
        ):
            return CALIBRATION.CalibrationManager.validate_selected_model_data(
                manager,
                "primary",
            )

    def test_small_field_is_checked_completely(self):
        selection = ShapedModelSelection(rows=50)
        qa = self.shaped_validation(selection, 1 << 20, 1 << 16)
        self.assertEqual(qa["coverage"], "complete")
        self.assertEqual(qa["rows_checked"], 50)
        self.assertEqual(qa["total_values"], 50 * 64 * 4)
        self.assertEqual(qa["status"], "passed")

    def test_large_field_is_sampled_with_bounded_reads(self):
        values_per_row = 64 * 4
        block_values = 64 * values_per_row      # 64 rows per read
        max_values = 8 * block_values          # ~8 blocks in total
        selection = ShapedModelSelection(rows=100000)

        qa = self.shaped_validation(selection, max_values, block_values)

        self.assertEqual(qa["coverage"], "sampled")
        self.assertEqual(qa["values_per_row"], values_per_row)
        self.assertLessEqual(qa["total_values"], max_values + block_values)
        block_reads = selection.reads[1:]      # first read is the probe
        for start, count in block_reads:
            self.assertLessEqual(count * values_per_row, block_values)
        starts = [start for start, _ in block_reads]
        self.assertIn(0, starts)
        self.assertIn(100000 - 64, starts)

    def test_bad_value_in_last_row_is_found_when_sampled(self):
        values_per_row = 64 * 4
        block_values = 64 * values_per_row
        selection = ShapedModelSelection(rows=100000, nan_row=99999)
        with self.assertRaisesRegex(RuntimeError, "non-finite"):
            self.shaped_validation(selection, 8 * block_values, block_values)

    def test_zero_model_is_rejected_when_sampled(self):
        values_per_row = 64 * 4
        block_values = 64 * values_per_row
        selection = ShapedModelSelection(rows=100000, zero=True)
        with self.assertRaisesRegex(RuntimeError, "entirely zero"):
            self.shaped_validation(selection, 8 * block_values, block_values)

if __name__ == "__main__":
    unittest.main()
