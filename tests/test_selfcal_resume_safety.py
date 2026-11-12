import importlib.util
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
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_selfcal_resume_safety",
        SOURCE_DIR / "selfcal.py",
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


class Stop(Exception):
    """Raised by a patched step to end a run at a known point."""


class FakeMainTable:
    """MAIN-table stand-in with CASA's (npol, nchan, nrow) column shape."""

    def __init__(self, columns):
        self.columns = columns
        self.reads = []

    def open(self, path):
        return None

    def close(self):
        return None

    def nrows(self):
        return next(iter(self.columns.values())).shape[-1]

    def colnames(self):
        return list(self.columns)

    def getcol(self, name, startrow=0, nrow=-1):
        column = self.columns[name]
        stop = column.shape[-1] if nrow < 0 else startrow + nrow
        self.reads.append((name, startrow, stop - startrow))
        return column[..., startrow:stop]


def manager_with_task_config(**task_config):
    manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
    manager.task_config = dict(task_config)
    return manager


class ModelColumnStatusTests(unittest.TestCase):

    def status(self, columns, **task_config):
        table = FakeMainTable(columns)
        manager = manager_with_task_config(**task_config)
        with mock.patch.object(SELFCAL, "tb", table, create=True):
            return manager.model_column_status("target.ms"), table

    @staticmethod
    def model(nrow, written_rows=None):
        values = np.ones((2, 4, nrow), dtype=np.complex64)
        if written_rows is not None:
            values[..., written_rows:] = 0
        return values

    def test_complete_model_passes(self):
        (ok, reason), _ = self.status({"MODEL_DATA": self.model(100)})
        self.assertTrue(ok, reason)

    def test_absent_model_column_is_reported(self):
        columns = {"DATA": self.model(10)}
        (ok, reason), _ = self.status(columns)
        self.assertFalse(ok)
        self.assertEqual(reason, "MODEL_DATA is absent")

    def test_partially_written_model_fails(self):
        # An interrupted write left rows 60 onward at zero. A check that stops
        # at the first non-zero value would accept this model.
        with mock.patch.object(SELFCAL, "MODEL_CHECK_BLOCK_VALUES", 80):
            (ok, reason), _ = self.status(
                {"MODEL_DATA": self.model(100, written_rows=60)}
            )
        self.assertFalse(ok)
        self.assertIn("rows 60-69", reason)
        self.assertIn("incompletely written", reason)

    def test_sampled_check_reaches_the_end_of_a_large_table(self):
        columns = {"MODEL_DATA": self.model(1000, written_rows=900)}
        with mock.patch.object(SELFCAL, "MODEL_CHECK_BLOCK_VALUES", 80):
            (ok, reason), table = self.status(columns, model_check_rows=40)
        self.assertFalse(ok)
        self.assertIn("rows 990-999", reason)

        rows_read = sum(n for _, start, n in table.reads[1:])
        self.assertLessEqual(rows_read, 40)

    def test_value_budget_bounds_the_rows_read(self):
        columns = {"MODEL_DATA": self.model(1000)}
        with mock.patch.object(SELFCAL, "MODEL_CHECK_MAX_VALUES", 160), \
                mock.patch.object(SELFCAL, "MODEL_CHECK_BLOCK_VALUES", 80):
            (ok, reason), table = self.status(columns)
        self.assertTrue(ok, reason)
        rows_read = sum(n for _, start, n in table.reads[1:])
        self.assertLessEqual(rows_read * 8, 160)
        self.assertEqual(table.reads[-1][1] + table.reads[-1][2], 1000)

    def test_flagged_zero_rows_are_not_a_truncation(self):
        model = self.model(100, written_rows=60)
        flags = np.zeros(model.shape, dtype=bool)
        flags[..., 60:] = True
        with mock.patch.object(SELFCAL, "MODEL_CHECK_BLOCK_VALUES", 80):
            (ok, reason), _ = self.status({"MODEL_DATA": model, "FLAG": flags})
        self.assertTrue(ok, reason)

    def test_unflagged_non_finite_model_fails(self):
        model = self.model(20)
        model[0, 1, 7] = np.nan
        (ok, reason), _ = self.status({"MODEL_DATA": model})
        self.assertFalse(ok)
        self.assertIn("non-finite", reason)

    def test_all_zero_model_fails(self):
        (ok, reason), _ = self.status({"MODEL_DATA": self.model(20, 0)})
        self.assertFalse(ok)
        self.assertIn("zero", reason)


class InterruptedProductRecoveryTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        self.msdir = root / "ms"
        self.imagedir = root / "images"
        self.caltabledir = root / "caltables"
        for directory in (self.msdir, self.imagedir, self.caltabledir):
            directory.mkdir()

        self.manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        self.manager.msdir = self.msdir
        self.manager.imagedir = self.imagedir
        self.manager.caltabledir = self.caltabledir
        self.manager.checkpoint_exists = lambda name: False

    def tearDown(self):
        self.temporary.cleanup()

    def make(self, directory, *names):
        for name in names:
            path = directory / name
            path.mkdir()
            (path / "table.dat").touch()

    def test_phase_loop_recovery_preserves_everything_it_moves(self):
        self.make(self.msdir, "phase2.candidate.ms", "phase2.ms", "phase1.ms")
        self.make(
            self.imagedir,
            "phase_2.image.tt0",
            "phase_2.model.tt0",
            "phase_1.image.tt0",
            "phase_20.image.tt0",
        )

        self.manager.recover_interrupted_phase_loop(2)

        for name in ("phase2.candidate.ms", "phase2.ms"):
            self.assertFalse((self.msdir / name).exists())
            self.assertTrue((self.msdir / f"{name}.incomplete.1").is_dir())
        for name in ("phase_2.image.tt0", "phase_2.model.tt0"):
            self.assertFalse((self.imagedir / name).exists())
            self.assertTrue((self.imagedir / f"{name}.incomplete.1").is_dir())

        # The accepted loop and an unrelated loop number stay in place.
        self.assertTrue((self.msdir / "phase1.ms").is_dir())
        self.assertTrue((self.imagedir / "phase_1.image.tt0").is_dir())
        self.assertTrue((self.imagedir / "phase_20.image.tt0").is_dir())
        self.assertEqual(len(self.manager.archived_incomplete_products), 4)

    def test_repeated_interruptions_never_overwrite_an_archive(self):
        self.make(self.imagedir, "final.image.tt0")
        self.manager.archive_incomplete_images(self.imagedir / "final", "x")
        self.make(self.imagedir, "final.image.tt0")
        self.manager.archive_incomplete_images(self.imagedir / "final", "x")

        self.assertTrue((self.imagedir / "final.image.tt0.incomplete.1").is_dir())
        self.assertTrue((self.imagedir / "final.image.tt0.incomplete.2").is_dir())
        self.assertFalse(
            (self.imagedir / "final.image.tt0.incomplete.1.incomplete.1").exists()
        )

    def test_interrupted_phase_loop_is_redone_instead_of_stopping(self):
        self.make(self.msdir, "phase1.candidate.ms")
        self.make(self.imagedir, "phase_1.image.tt0")
        manager = self.manager
        manager.selfcal_strategy = {
            "phase_selfcal": True,
            "phase_solints": ["20s"],
        }
        manager.image_stats = {"dynamic_range": 50.0}
        manager.current_ms = str(self.msdir / "target.ms")
        manager.phase_rejected = lambda loop: False
        manager.ensure_model_column = mock.Mock(side_effect=Stop)

        with self.assertRaises(Stop):
            manager.run_phase_selfcal()

        self.assertTrue(
            (self.msdir / "phase1.candidate.ms.incomplete.1").is_dir()
        )
        self.assertTrue(
            (self.imagedir / "phase_1.image.tt0.incomplete.1").is_dir()
        )
        manager.ensure_model_column.assert_called_once()

    def test_interrupted_amp_selfcal_is_redone_instead_of_stopping(self):
        # Killed while the candidate was being imaged, or during the split
        # that commits it: all three are left behind.
        self.make(self.msdir, "final_selfcal.ms", "amp_phase.candidate.ms")
        self.make(self.imagedir, "amp_phase.image.tt0")
        manager = self.manager
        manager.task_config = {"enable_amp_selfcal": True}
        manager.band_config = {"amp_solint": "inf"}
        manager.selfcal_strategy = {"amp_selfcal": True, "amp_solint": "inf"}
        manager.checkdir = self.msdir.parent
        manager.phase_tables = [str(self.caltabledir / "phase_1.cal")]
        manager.current_ms = str(self.msdir / "phase1.ms")
        manager.ensure_model_column = mock.Mock(side_effect=Stop)

        with self.assertRaises(Stop):
            manager.run_amp_phase_selfcal()

        for name in ("final_selfcal.ms", "amp_phase.candidate.ms"):
            self.assertFalse((self.msdir / name).exists())
            self.assertTrue((self.msdir / f"{name}.incomplete.1").is_dir())
        self.assertTrue(
            (self.imagedir / "amp_phase.image.tt0.incomplete.1").is_dir()
        )

    def test_partial_images_are_moved_before_tclean_reuses_them(self):
        for method, base in (
            ("initial_imaging", "initial"),
            ("final_imaging", "final"),
        ):
            with self.subTest(method=method):
                archive = mock.Mock(side_effect=Stop)
                self.manager.archive_incomplete_images = archive
                with self.assertRaises(Stop):
                    getattr(self.manager, method)()
                self.assertEqual(
                    Path(archive.call_args[0][0]),
                    self.imagedir / base,
                )

    def test_checkpointed_images_are_left_alone(self):
        self.make(self.imagedir, "final.image.tt0")
        self.manager.checkpoint_exists = lambda name: name == "final_image"
        self.manager.clean_diagnostics = {}
        self.manager.deconvolution_diagnostics = lambda base: {}
        self.manager.archive_incomplete_images = mock.Mock()

        self.manager.final_imaging()

        self.manager.archive_incomplete_images.assert_not_called()
        self.assertTrue((self.imagedir / "final.image.tt0").is_dir())


if __name__ == "__main__":
    unittest.main()
