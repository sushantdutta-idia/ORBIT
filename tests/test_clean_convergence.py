"""CLEAN is judged on what is really in the image, and runs until done."""

import importlib.util
import math
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from test_clean_bootstrap import IMAGING
from test_selfcal_coherence import SELFCAL


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_radio_continuum_casa_compat_clean", SOURCE_DIR / "casa_compat.py"
)
COMPAT = importlib.util.module_from_spec(spec)
spec.loader.exec_module(COMPAT)


def casa_like_imstat(pixels):
    """imstat as CASA 6 behaves: clipping removes outliers from *every*
    statistic, including max, min and their positions. Checked against
    casatools.image.statistics: a 60-sigma pixel is reported by chauvenet
    as a 4-sigma maximum elsewhere."""

    pixels = np.asarray(pixels, dtype=float)

    def imstat(imagename=None, algorithm="classic", maxiter=-1, box=None,
               **_):
        data = pixels
        offset = (0, 0)
        if box:
            x0, y0, x1, y1 = (int(v) for v in box.split(","))
            data = pixels[x0:x1 + 1, y0:y1 + 1]
            offset = (x0, y0)
        values = data.copy()
        keep = np.ones(values.shape, dtype=bool)
        if algorithm == "chauvenet":
            for _ in range(max(1, maxiter)):
                sample = values[keep]
                centre, spread = sample.mean(), sample.std()
                new = np.abs(values - centre) <= 4.0 * spread
                if np.array_equal(new, keep):
                    break
                keep = new
        kept = np.where(keep, values, np.nan)
        maxpos = np.unravel_index(np.nanargmax(kept), kept.shape)
        minpos = np.unravel_index(np.nanargmin(kept), kept.shape)
        return {
            "max": [np.nanmax(kept)],
            "min": [np.nanmin(kept)],
            "maxpos": [maxpos[0] + offset[0], maxpos[1] + offset[1], 0, 0],
            "minpos": [minpos[0] + offset[0], minpos[1] + offset[1], 0, 0],
            "rms": [math.sqrt(np.nanmean(kept ** 2))],
            "sum": [np.nansum(kept)],
            "sumsq": [np.nansum(kept ** 2)],
        }

    return imstat


def field(peak=6e-3, where=(200, 50), sigma=1e-4, size=256, seed=0):
    image = np.random.default_rng(seed).normal(0, sigma, (size, size))
    image[where] = peak
    return image


class PeakAndNoiseTests(unittest.TestCase):

    def test_clipped_statistics_cannot_see_a_source(self):
        imstat = casa_like_imstat(field())
        clipped = imstat(algorithm="chauvenet", maxiter=5)
        # This is what the pipeline used to call the peak.
        self.assertLess(clipped["max"][0] / clipped["rms"][0], 5.0)

    def test_the_true_peak_and_the_robust_noise_are_combined(self):
        result = COMPAT.peak_and_noise(casa_like_imstat(field()), "x.image")

        self.assertAlmostEqual(result["max"], 6e-3)
        self.assertEqual(result["max_position"], [200, 50])
        self.assertAlmostEqual(result["rms"] / 1e-4, 1.0, delta=0.05)
        self.assertAlmostEqual(result["dynamic_range"], 60.0, delta=3.0)

    def test_a_negative_extreme_is_found_too(self):
        result = COMPAT.peak_and_noise(
            casa_like_imstat(field(peak=-4e-3, where=(10, 20))), "x.image"
        )
        self.assertAlmostEqual(result["peak_abs"], 4e-3)
        self.assertEqual(result["peak_abs_position"], [10, 20])
        self.assertAlmostEqual(result["peak_abs_snr"], 40.0, delta=2.0)


class NoiseFloorTests(unittest.TestCase):
    """A significant peak is one noise alone would not reach."""

    def test_the_floor_grows_with_the_image(self):
        small = COMPAT.noise_peak_snr(512 ** 2)
        large = COMPAT.noise_peak_snr(8192 ** 2)
        self.assertGreater(large, small)
        self.assertAlmostEqual(small, 5.08, delta=0.05)
        self.assertAlmostEqual(large, 6.05, delta=0.05)
        self.assertIsNone(COMPAT.noise_peak_snr(0))

    def test_the_floor_is_what_noise_actually_reaches(self):
        # 200 independent 512x512 noise images at one pixel per beam: the
        # 1% level should be exceeded about 1% of the time, 5 sigma far more.
        rng = np.random.default_rng(11)
        floor = COMPAT.noise_peak_snr(512 ** 2, pixels_per_beam=1.0)
        peaks = [
            np.max(np.abs(rng.standard_normal((512, 512))))
            for _ in range(200)
        ]
        self.assertLessEqual(sum(peak > floor for peak in peaks), 6)
        self.assertGreater(sum(peak > 5.0 for peak in peaks), 20)

    def test_a_pure_noise_field_is_not_refused_as_undeconvolved(self):
        manager = IMAGING.FinalImagingManager.__new__(
            IMAGING.FinalImagingManager
        )
        manager.clean_diagnostics = {"final_continuum": {
            "model_nonzero": False,
            "residual_peak_snr": 5.3,
            "residual_rms_Jy_beam": 1e-4,
            "residual_peak_abs_Jy_beam": 5.3e-4,
            "noise_only_peak_snr": COMPAT.noise_peak_snr(2048 ** 2),
        }}
        manager.check_pixel_sampling = lambda label: None
        manager.describe_residual_peak = lambda label: None
        with mock.patch("builtins.print"):
            manager.verify_deconvolution("final_continuum")

    def test_a_real_residual_source_is_still_refused(self):
        manager = IMAGING.FinalImagingManager.__new__(
            IMAGING.FinalImagingManager
        )
        manager.clean_diagnostics = {"final_continuum": {
            "model_nonzero": False,
            "residual_peak_snr": 45.0,
            "residual_rms_Jy_beam": 1.34e-4,
            "residual_peak_abs_Jy_beam": 6.07e-3,
            "noise_only_peak_snr": COMPAT.noise_peak_snr(2048 ** 2),
        }}
        manager.check_pixel_sampling = lambda label: None
        manager.describe_residual_peak = lambda label: None
        with self.assertRaisesRegex(RuntimeError, "dirty map"):
            manager.verify_deconvolution("final_continuum")

    def test_triggers_use_the_larger_of_config_and_noise(self):
        for module in (IMAGING, SELFCAL):
            with self.subTest(module=module.__name__):
                self.assertEqual(
                    module.significance_floor(5.0,
                                              {"noise_only_peak_snr": 5.6}),
                    5.6,
                )
                self.assertEqual(
                    module.significance_floor(8.0,
                                              {"noise_only_peak_snr": 5.6}),
                    8.0,
                )
                self.assertEqual(module.significance_floor(5.0, {}), 5.0)

    def test_noise_does_not_start_a_retry(self):
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.task_config = {"clean_bootstrap_enabled": False,
                               "automask_retry_enabled": True,
                               "automask_fallback_enabled": True}
        manager.clean_diagnostics = {}
        manager.clear_stale_clean_mask = lambda base: None
        manager.deconvolution_diagnostics = lambda base: {
            "model_nonzero": False,
            "residual_peak_snr": 5.4,
            "residual_rms_Jy_beam": 1e-4,
            "residual_peak_abs_Jy_beam": 5.4e-4,
            "noise_only_peak_snr": 6.05,
        }
        manager.automask_kwargs = lambda relaxed=False: {}
        calls = []
        with mock.patch.object(
            SELFCAL, "tclean", lambda **kw: calls.append(kw)
        ), mock.patch("builtins.print"):
            diagnostics = manager.run_tclean_with_recovery(
                {"imagename": "initial"}, "initial"
            )
        self.assertEqual(len(calls), 1)
        self.assertFalse(diagnostics["recovery_attempted"])


class StopSummaryTests(unittest.TestCase):

    def test_casa_stop_codes_are_named_and_judged(self):
        niter = COMPAT.summarize_tclean({
            "stopcode": 1, "iterdone": 50000, "niter": 50000,
            "nmajordone": 7, "maxpsfsidelobe": 0.18,
            "stopDescription": "iteration limit",
        })
        self.assertEqual(niter["stop_reason"],
                         "reached the iteration limit (niter)")
        self.assertFalse(niter["converged"])
        self.assertEqual(niter["casa_stop_description"], "iteration limit")
        self.assertAlmostEqual(niter["max_psf_sidelobe"], 0.18)

        for code in (2, 7, 8):
            self.assertTrue(COMPAT.summarize_tclean({"stopcode": code})
                            ["converged"])
        for code in (1, 3, 4, 5, 6, 9):
            self.assertFalse(COMPAT.summarize_tclean({"stopcode": code})
                             ["converged"])

    def test_no_return_value_is_not_guessed_at(self):
        self.assertIsNone(COMPAT.summarize_tclean(None))
        self.assertIsNone(COMPAT.summarize_tclean({})["converged"])


class BootstrapRegionTests(unittest.TestCase):

    def test_the_bootstrap_reaches_as_far_as_the_deconvolution(self):
        self.assertEqual(COMPAT.bootstrap_pbmask(None, 0.05), 0.05)
        self.assertEqual(COMPAT.bootstrap_pbmask(None, -0.05), 0.05)
        self.assertEqual(COMPAT.bootstrap_pbmask(0.3, 0.05), 0.3)
        self.assertEqual(COMPAT.bootstrap_pbmask(None, None), 0.2)
        self.assertEqual(COMPAT.bootstrap_pbmask(None, 0), 0.2)

    def test_a_selfcal_bootstrap_uses_the_selfcal_pblimit(self):
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.task_config = {}
        manager.clean_diagnostics = {}
        manager.clear_stale_clean_mask = lambda base: None
        manager.deconvolution_diagnostics = lambda base: {
            "model_nonzero": True
        }
        calls = []
        with mock.patch.object(
            SELFCAL, "tclean", lambda **kw: calls.append(kw)
        ), mock.patch("builtins.print"):
            manager.bootstrap_clean(
                {"imagename": "initial", "pblimit": 0.05}, "initial"
            )
        self.assertEqual(calls[0]["pbmask"], 0.05)


class ConvergenceTests(unittest.TestCase):

    def imaging(self, sequence):
        value = IMAGING.FinalImagingManager.__new__(
            IMAGING.FinalImagingManager
        )
        value.clean_diagnostics = {}
        self.diagnostics = list(sequence)
        value.deconvolution_diagnostics = lambda base: dict(
            self.diagnostics.pop(0)
        )
        return value

    def run_continue(self, returns, sequence, **config):
        value = self.imaging(sequence)
        results = list(returns)
        calls = []

        def fake_tclean(**kwargs):
            calls.append(kwargs)
            return results.pop(0)

        with mock.patch.object(IMAGING, "tclean", fake_tclean), \
                mock.patch.dict(IMAGING.IMAGING_TASK_CONFIG, config), \
                mock.patch("builtins.print"):
            first = COMPAT.summarize_tclean({"stopcode": 1})
            diagnostics = value.continue_to_convergence(
                {"imagename": "final", "niter": 50000, "restart": False,
                 "calcpsf": True, "calcres": True},
                "final_continuum",
                first,
                {"residual_peak_snr": 45.0, "bootstrap": {"b": 1}},
            )
        return calls, diagnostics, value

    def test_an_unfinished_clean_is_carried_on_until_it_converges(self):
        calls, diagnostics, value = self.run_continue(
            returns=[{"stopcode": 1}, {"stopcode": 8}],
            sequence=[{"residual_peak_snr": 12.0},
                      {"residual_peak_snr": 4.1}],
        )
        self.assertEqual(len(calls), 2)
        for call in calls:
            self.assertTrue(call["restart"])
            self.assertFalse(call["calcpsf"])
            self.assertFalse(call["calcres"])
        self.assertEqual(diagnostics["continuations"], 2)
        self.assertTrue(diagnostics["final_stop"]["converged"])
        self.assertEqual(diagnostics["bootstrap"], {"b": 1})
        self.assertEqual(
            [p["name"] for p in value.clean_passes["final_continuum"]],
            ["continuation 1", "continuation 2"],
        )

    def test_continuation_is_bounded(self):
        calls, diagnostics, _ = self.run_continue(
            returns=[{"stopcode": 1}] * 5,
            sequence=[{"residual_peak_snr": 30.0}] * 5,
            clean_max_continuations=2,
        )
        self.assertEqual(len(calls), 2)
        self.assertFalse(diagnostics["final_stop"]["converged"])

    def test_a_clean_that_finished_is_left_alone(self):
        value = self.imaging([])
        with mock.patch.object(IMAGING, "tclean", mock.Mock()) as tclean:
            diagnostics = value.continue_to_convergence(
                {"imagename": "final"}, "final_continuum",
                COMPAT.summarize_tclean({"stopcode": 8}),
                {"residual_peak_snr": 45.0},
            )
        tclean.assert_not_called()
        self.assertEqual(diagnostics["continuations"], 0)

    def test_a_faint_residual_does_not_justify_more_iterations(self):
        value = self.imaging([])
        with mock.patch.object(IMAGING, "tclean", mock.Mock()) as tclean:
            value.continue_to_convergence(
                {"imagename": "final"}, "final_continuum",
                COMPAT.summarize_tclean({"stopcode": 1}),
                {"residual_peak_snr": 4.0},
            )
        tclean.assert_not_called()


class ResidualPeakTests(unittest.TestCase):
    """The largest residual is located and its cause named."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = str(Path(self.temporary.name) / "final")
        size = 256
        yy, xx = np.mgrid[:size, :size]
        radius = np.hypot(xx - 128, yy - 128)
        self.pb = np.exp(-4 * np.log(2) * (radius / 120.0) ** 2)
        self.mask = (self.pb > 0.3).astype(float)
        self.residual = field(peak=6.07e-3, where=(230, 200), sigma=1.34e-4)
        for suffix in (".residual.tt0", ".pb.tt0", ".mask", ".image.tt0"):
            path = Path(self.base + suffix)
            path.mkdir(parents=True)
            (path / "table.dat").touch()

    def tearDown(self):
        self.temporary.cleanup()

    def manager(self, stop, nterms=2):
        value = IMAGING.FinalImagingManager.__new__(
            IMAGING.FinalImagingManager
        )
        value.imagename = self.base
        value.image_tt0 = self.base + ".image.tt0"
        value.config = {"nterms": nterms}
        value.science_reffreq = {"valid_hz": [550e6, 750e6]}
        diagnostics = COMPAT.peak_and_noise(
            casa_like_imstat(self.residual), "residual"
        )
        value.clean_diagnostics = {"final_continuum": {
            "residual_peak_position": diagnostics["peak_abs_position"],
            "residual_peak_abs_Jy_beam": diagnostics["peak_abs"],
            "residual_peak_snr": diagnostics["peak_abs_snr"],
            "final_stop": stop,
        }}
        return value

    def describe(self, value, mask=None):
        images = {
            ".pb.tt0": self.pb,
            ".mask": self.mask if mask is None else mask,
            ".residual.tt0": self.residual,
        }

        def imstat(imagename=None, **kwargs):
            for suffix, pixels in images.items():
                if str(imagename).endswith(suffix):
                    return casa_like_imstat(pixels)(**kwargs)
            raise AssertionError(imagename)

        header = {
            "refpix": [128.0, 128.0, 0.0, 0.0],
            "incr": [-2.0 / 206264.8, 2.0 / 206264.8, 1.0, 1.0],
            "axisunits": ["rad", "rad", "", "Hz"],
        }
        printed = []
        with mock.patch.object(IMAGING, "imstat", imstat), \
                mock.patch.object(IMAGING, "imhead",
                                  lambda **kw: header), \
                mock.patch("builtins.print",
                           lambda *a, **k: printed.append(" ".join(map(str, a)))):
            record = value.describe_residual_peak("final_continuum")
        return record, " ".join(printed)

    def test_a_source_outside_the_mask_is_named_as_never_deconvolved(self):
        record, text = self.describe(
            self.manager(COMPAT.summarize_tclean({"stopcode": 8}))
        )
        self.assertEqual(record["pixel"], [230, 200])
        self.assertEqual(record["cause"], "outside_mask")
        self.assertFalse(record["inside_clean_mask"])
        self.assertAlmostEqual(
            record["offset_arcmin"],
            math.hypot(102, 72) * 2.0 / 60.0,
            places=3,
        )
        self.assertLess(record["casa_primary_beam"], 0.3)
        self.assertIn("never deconvolved", text)
        self.assertIn("6.07 mJy/beam (45 sigma)", text)

    def test_an_unfinished_clean_is_blamed_when_it_is_the_cause(self):
        record, text = self.describe(
            self.manager(COMPAT.summarize_tclean({"stopcode": 1})),
            mask=np.ones_like(self.pb),
        )
        self.assertEqual(record["cause"], "not_converged")
        self.assertIn("iteration limit", text)

    def test_a_modelled_region_gives_the_beam_curvature_it_implies(self):
        record, text = self.describe(
            self.manager(COMPAT.summarize_tclean({"stopcode": 8})),
            mask=np.ones_like(self.pb),
        )
        self.assertEqual(record["cause"], "not_modelled")
        pb = record["casa_primary_beam"]
        half = 0.5 * math.log(750.0 / 550.0)
        self.assertAlmostEqual(
            record["beam_spectral_misfit"],
            (2.0 / 3.0) * abs(2.0 * math.log(pb)) * half ** 2,
        )
        self.assertIn("nterms=3", text)

    def test_nterms_three_is_not_blamed_for_curvature(self):
        record, _ = self.describe(
            self.manager(COMPAT.summarize_tclean({"stopcode": 8}), nterms=3),
            mask=np.ones_like(self.pb),
        )
        self.assertNotIn("beam_spectral_misfit", record)

    def test_an_unknown_stop_is_not_called_convergence(self):
        record, text = self.describe(self.manager({}),
                                     mask=np.ones_like(self.pb))
        self.assertEqual(record["cause"], "unknown")


class PixelSamplingTests(unittest.TestCase):

    def check(self, minor_arcsec, cell_arcsec):
        value = IMAGING.FinalImagingManager.__new__(
            IMAGING.FinalImagingManager
        )
        value.image_tt0 = "final.image.tt0"
        value.clean_diagnostics = {}
        header = {
            "refpix": [1024.0, 1024.0],
            "incr": [-cell_arcsec / 206264.806, cell_arcsec / 206264.806],
            "axisunits": ["rad", "rad"],
            "restoringbeam": {
                "major": {"value": minor_arcsec * 1.3, "unit": "arcsec"},
                "minor": {"value": minor_arcsec, "unit": "arcsec"},
            },
        }
        printed = []
        with mock.patch.object(IMAGING, "imhead", lambda **kw: header), \
                mock.patch("builtins.print",
                           lambda *a, **k: printed.append(" ".join(map(str, a)))):
            pixels = value.check_pixel_sampling("final_continuum")
        return pixels, " ".join(printed)

    def test_coarse_pixels_are_warned_about_with_the_cell_to_use(self):
        pixels, text = self.check(4.2, 2.0)
        self.assertAlmostEqual(pixels, 2.1, places=3)
        self.assertIn("WARNING", text)
        self.assertIn('at most 1.05"', text)

    def test_adequate_sampling_is_only_reported(self):
        pixels, text = self.check(4.2, 1.0)
        self.assertAlmostEqual(pixels, 4.2, places=3)
        self.assertNotIn("WARNING", text)


class SelfcalStartUsesTheTruePeakTests(unittest.TestCase):

    def test_an_image_at_dynamic_range_126_is_self_calibrated(self):
        # The observer's image: 19.3 mJy/beam over 153 uJy/beam. Clipped
        # statistics made this look like dynamic range 4.
        image = field(peak=19.3e-3, where=(128, 128), sigma=1.53e-4, seed=4)
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.initial_image = "initial.image.tt0"
        with mock.patch.object(SELFCAL, "imstat", casa_like_imstat(image)), \
                mock.patch("builtins.print"):
            manager.estimate_snr()
        self.assertAlmostEqual(
            manager.image_stats["dynamic_range"], 126.0, delta=6.0
        )
        self.assertAlmostEqual(manager.image_stats["peak_Jy_beam"], 19.3e-3)


if __name__ == "__main__":
    unittest.main()
