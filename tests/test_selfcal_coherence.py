import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from test_selfcal_resume_safety import load_selfcal_module


SELFCAL = load_selfcal_module()


class FakeVisibilityTable:
    """MAIN-table stand-in with CASA's (npol, nchan, nrow) column shape."""

    def __init__(self, data, flags=None, antenna1=None, antenna2=None,
                 corrected=False):
        self.columns = {"CORRECTED_DATA" if corrected else "DATA": data}
        if flags is not None:
            self.columns["FLAG"] = flags
        nrow = data.shape[-1]
        self.columns["ANTENNA1"] = (
            np.zeros(nrow, dtype=int) if antenna1 is None else antenna1
        )
        self.columns["ANTENNA2"] = (
            np.ones(nrow, dtype=int) if antenna2 is None else antenna2
        )
        self.reads = []

    def open(self, path):
        return None

    def close(self):
        return None

    def nrows(self):
        return self.columns["ANTENNA1"].shape[-1]

    def colnames(self):
        return list(self.columns)

    def getcol(self, name, startrow=0, nrow=-1):
        column = self.columns[name]
        stop = column.shape[-1] if nrow < 0 else startrow + nrow
        self.reads.append((name, startrow, stop - startrow))
        return column[..., startrow:stop]


def manager(**task_config):
    value = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
    value.task_config = dict(task_config)
    value.current_ms = "target.ms"
    return value


class CoherenceTests(unittest.TestCase):

    @staticmethod
    def visibilities(phases, amplitude=1.0, shape=(2, 4)):
        rows = len(phases)
        values = np.empty(shape + (rows,), dtype=np.complex128)
        for index, phase in enumerate(phases):
            values[..., index] = amplitude * np.exp(1j * phase)
        return values

    def measure(self, table, **task_config):
        with mock.patch.object(SELFCAL, "tb", table, create=True):
            return manager(**task_config).measure_coherence("target.ms")

    @staticmethod
    def noisy_field(flux, sigma, rows, seed=7):
        """Visibilities of a point source at the phase centre, plus noise."""

        rng = np.random.default_rng(seed)
        shape = (2, 4, rows)
        noise = (
            rng.normal(scale=sigma, size=shape) +
            1j * rng.normal(scale=sigma, size=shape)
        )
        return flux + noise

    def test_a_real_field_is_coherent_even_though_noise_dominates(self):
        # 4.5 mJy of flux under 0.27 Jy of per-visibility noise, as in a
        # uGMRT target field: the vector/scalar ratio is ~1%, and the
        # vector average is still many sigma from zero.
        result = self.measure(
            FakeVisibilityTable(self.noisy_field(4.5e-3, 0.274, 20000))
        )

        self.assertEqual(result["status"], "measured")
        self.assertLess(result["vector_to_scalar_ratio"], 0.05)
        self.assertGreater(result["significance"], 5.0)
        self.assertAlmostEqual(
            result["component_sigma"], 0.274, delta=0.02
        )

    def test_random_phases_leave_no_significant_coherent_flux(self):
        rng = np.random.default_rng(2)
        phases = rng.uniform(0.0, 2.0 * np.pi, size=20000)
        result = self.measure(FakeVisibilityTable(self.visibilities(phases)))

        self.assertLess(result["significance"], 5.0)

    def test_the_uncertainty_falls_as_the_sample_count_grows(self):
        small = self.measure(
            FakeVisibilityTable(self.noisy_field(4.5e-3, 0.274, 500))
        )
        large = self.measure(
            FakeVisibilityTable(self.noisy_field(4.5e-3, 0.274, 20000))
        )
        self.assertGreater(small["vector_mean_error"],
                           large["vector_mean_error"])

    def test_flagged_samples_and_autocorrelations_are_excluded(self):
        phases = np.zeros(10)
        data = self.visibilities(phases)
        data[..., 5:] = 50.0  # flagged junk
        flags = np.zeros(data.shape, dtype=bool)
        flags[..., 5:] = True
        antenna1 = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        antenna2 = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        table = FakeVisibilityTable(data, flags, antenna1, antenna2)

        result = self.measure(table)

        self.assertAlmostEqual(result["scalar_mean_amplitude"], 1.0, places=6)
        self.assertEqual(result["autocorrelations_skipped"], 5)

    def test_reads_are_bounded_and_spread_over_the_table(self):
        phases = np.zeros(5000)
        table = FakeVisibilityTable(self.visibilities(phases))

        with mock.patch.object(SELFCAL, "tb", table, create=True):
            manager().measure_coherence(
                "target.ms",
                rows_wanted=200,
                block_rows=50,
            )

        data_reads = [read for read in table.reads if read[0] == "DATA"]
        self.assertLessEqual(sum(read[2] for read in data_reads), 200)
        self.assertGreater(data_reads[-1][1], 0)

    def test_corrected_column_is_preferred(self):
        table = FakeVisibilityTable(
            self.visibilities(np.zeros(20)),
            corrected=True,
        )
        self.assertEqual(self.measure(table)["column"], "CORRECTED_DATA")

    def test_a_failure_is_reported_not_raised(self):
        class Broken:
            def open(self, path):
                raise RuntimeError("table is locked")

            def close(self):
                return None

        with mock.patch.object(SELFCAL, "tb", Broken(), create=True):
            value = manager()
            self.assertIsNone(value.measure_coherence("target.ms"))
            self.assertIsNone(value.check_calibrated_coherence())

    def test_an_insignificant_result_is_reported_as_inconclusive(self):
        rng = np.random.default_rng(3)
        table = FakeVisibilityTable(
            self.visibilities(rng.uniform(0.0, 2.0 * np.pi, size=20000))
        )
        value = manager(coherence_min_significance=5.0)
        with mock.patch.object(SELFCAL, "tb", table, create=True):
            measurement = value.check_calibrated_coherence()

        self.assertEqual(measurement["status"], "measured")
        self.assertIs(value.coherence, measurement)
        self.assertLess(measurement["significance"], 5.0)


if __name__ == "__main__":
    unittest.main()


class EmptyModelRecoveryTests(unittest.TestCase):
    """The relaxed-mask retry is what recovers a field CLEAN could not mask."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.imagename = str(Path(self.temporary.name) / "initial")
        self.make_image(".pb.tt0")

    def tearDown(self):
        self.temporary.cleanup()

    def make_image(self, suffix):
        path = Path(self.imagename + suffix)
        path.mkdir(parents=True, exist_ok=True)
        (path / "table.dat").touch()
        return path

    def manager(self, diagnostics, **task_config):
        value = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        value.task_config = dict(
            {
                "automask_retry_enabled": True,
                # These cases are about what happens after auto-multithresh
                # opens no mask; the bootstrap that usually prevents that is
                # exercised separately.
                "clean_bootstrap_enabled": False,
            },
            **task_config,
        )
        self.sequence = list(diagnostics)
        value.deconvolution_diagnostics = lambda base: dict(
            self.sequence.pop(0)
        )
        value.automask_kwargs = lambda relaxed=False: {
            "usemask": "auto-multithresh",
            "sidelobethreshold": 1.5 if relaxed else 2.5,
        }
        return value

    def test_a_masked_out_field_is_retried_at_five_sigma(self):
        empty = {
            "model_nonzero": False,
            "residual_peak_snr": 6.0,
            "residual_rms_Jy_beam": 1e-4,
            "residual_peak_abs_Jy_beam": 6e-4,
        }
        recovered = dict(empty, model_nonzero=True)
        manager = self.manager([empty, recovered])
        calls = []

        with mock.patch.object(
            SELFCAL, "tclean", lambda **kwargs: calls.append(kwargs)
        ):
            diagnostics = manager.run_tclean_with_recovery(
                {"imagename": self.imagename, "vis": "target.ms"},
                "initial",
            )

        self.assertEqual(len(calls), 2)
        self.assertEqual(calls[1]["sidelobethreshold"], 1.5)
        self.assertTrue(diagnostics["recovery_attempted"])
        self.assertTrue(diagnostics["recovery_succeeded"])

    def test_a_field_no_mask_can_hold_is_cleaned_inside_the_beam(self):
        # A snapshot PSF whose sidelobes sit near the peak: auto-multithresh
        # opens no mask, so neither the first pass nor the relaxed retry
        # produces a model.
        empty = {
            "model_nonzero": False,
            "residual_peak_snr": 12.0,
            "residual_rms_Jy_beam": 2e-4,
            "residual_peak_abs_Jy_beam": 2.4e-3,
        }
        recovered = dict(empty, model_nonzero=True)
        manager = self.manager(
            [empty, empty, recovered],
            automask_fallback_nsigma=6.0,
            automask_fallback_niter=5000,
            automask_fallback_pbmask=0.2,
        )
        calls = []

        with mock.patch.object(
            SELFCAL, "tclean", lambda **kwargs: calls.append(kwargs)
        ):
            diagnostics = manager.run_tclean_with_recovery(
                {"imagename": self.imagename, "vis": "target.ms",
                 "niter": 100000, "nsigma": 3.0},
                "initial",
            )

        self.assertEqual(len(calls), 3)
        fallback = calls[2]
        self.assertEqual(fallback["usemask"], "pb")
        self.assertEqual(fallback["pbmask"], 0.2)
        self.assertEqual(fallback["nsigma"], 6.0)
        # Bounded: it cannot chase sidelobes down into the noise.
        self.assertEqual(fallback["niter"], 5000)
        self.assertTrue(diagnostics["fallback_attempted"])
        self.assertTrue(diagnostics["fallback_succeeded"])

    def test_a_stale_mask_is_removed_before_each_new_pass(self):
        # Restarting tclean over the previous pass's mask mixes two ways of
        # making one, which is where CASA became fragile.
        self.make_image(".mask")
        empty = {
            "model_nonzero": False,
            "residual_peak_snr": 12.0,
            "residual_rms_Jy_beam": 2e-4,
            "residual_peak_abs_Jy_beam": 2.4e-3,
        }
        manager = self.manager([empty, empty, dict(empty, model_nonzero=True)])
        removed = []

        # Every tclean pass writes its own mask.
        def fake_tclean(**kwargs):
            self.make_image(".mask")

        with mock.patch.object(SELFCAL, "tclean", fake_tclean), \
                mock.patch.object(
                    SELFCAL,
                    "rmtables",
                    lambda tablenames: removed.append(tablenames),
                ):
            manager.run_tclean_with_recovery(
                {"imagename": self.imagename, "vis": "target.ms"},
                "initial",
            )

        self.assertEqual(removed, [self.imagename + ".mask"] * 2)

    def test_the_beam_clean_needs_a_primary_beam_image(self):
        # usemask='pb' reads it, and calcpsf=False does not make one.
        Path(self.imagename + ".pb.tt0").rename(
            Path(self.imagename + ".pb.unused")
        )
        empty = {
            "model_nonzero": False,
            "residual_peak_snr": 12.0,
            "residual_rms_Jy_beam": 2e-4,
            "residual_peak_abs_Jy_beam": 2.4e-3,
        }
        manager = self.manager([empty, empty])
        calls = []

        with mock.patch.object(
            SELFCAL, "tclean", lambda **kwargs: calls.append(kwargs)
        ), mock.patch.object(SELFCAL, "rmtables", lambda tablenames: None):
            diagnostics = manager.run_tclean_with_recovery(
                {"imagename": self.imagename, "vis": "target.ms"},
                "initial",
            )

        # The first pass and the relaxed retry only.
        self.assertEqual(len(calls), 2)
        self.assertNotIn("fallback_attempted", diagnostics)

    def test_the_fallback_respects_its_own_threshold_and_switch(self):
        empty = {
            "model_nonzero": False,
            "residual_peak_snr": 3.0,
            "residual_rms_Jy_beam": 2e-4,
            "residual_peak_abs_Jy_beam": 6e-4,
        }
        manager = self.manager([empty])
        calls = []
        with mock.patch.object(
            SELFCAL, "tclean", lambda **kwargs: calls.append(kwargs)
        ):
            manager.run_tclean_with_recovery(
                {"imagename": self.imagename, "vis": "target.ms"}, "initial"
            )
        self.assertEqual(len(calls), 1)

        bright = dict(empty, residual_peak_snr=30.0)
        manager = self.manager(
            [bright, bright],
            automask_fallback_enabled=False,
        )
        calls = []
        with mock.patch.object(
            SELFCAL, "tclean", lambda **kwargs: calls.append(kwargs)
        ):
            diagnostics = manager.run_tclean_with_recovery(
                {"imagename": self.imagename, "vis": "target.ms"}, "initial"
            )
        # The relaxed retry still runs; the beam clean does not.
        self.assertEqual(len(calls), 2)
        self.assertNotIn("fallback_attempted", diagnostics)

    def test_a_field_below_the_retry_threshold_is_left_alone(self):
        empty = {
            "model_nonzero": False,
            "residual_peak_snr": 3.0,
            "residual_rms_Jy_beam": 1e-4,
            "residual_peak_abs_Jy_beam": 3e-4,
        }
        manager = self.manager([empty])
        calls = []

        with mock.patch.object(
            SELFCAL, "tclean", lambda **kwargs: calls.append(kwargs)
        ):
            diagnostics = manager.run_tclean_with_recovery(
                {"imagename": self.imagename, "vis": "target.ms"},
                "initial",
            )

        self.assertEqual(len(calls), 1)
        self.assertFalse(diagnostics["recovery_attempted"])
