"""Self-calibration starts when its solutions can be trusted, not later."""

import importlib.util
import math
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from test_clean_bootstrap import IMAGING
# The same module instance the shared helpers build managers from, so that
# patching its table tool reaches them.
from test_selfcal_coherence import SELFCAL, FakeVisibilityTable, manager


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_radio_continuum_qa_products_stretch", SOURCE_DIR / "qa_products.py"
)
QA = importlib.util.module_from_spec(spec)
spec.loader.exec_module(QA)


def field(rows=4000, antennas=6, times=200, interval=8.0, sigma=0.3, seed=3):
    """A small MAIN table: every baseline at every timestamp, two scans."""

    rng = np.random.default_rng(seed)
    pairs = [(a, b) for a in range(antennas) for b in range(a + 1, antennas)]
    rows = len(pairs) * times
    shape = (2, 4, rows)
    data = (rng.normal(scale=sigma, size=shape)
            + 1j * rng.normal(scale=sigma, size=shape))
    table = FakeVisibilityTable(
        data,
        flags=np.zeros(shape, dtype=bool),
        antenna1=np.array([a for _ in range(times) for a, _ in pairs]),
        antenna2=np.array([b for _ in range(times) for _, b in pairs]),
    )
    timestamps = np.repeat(np.arange(times) * interval, len(pairs))
    table.columns["TIME"] = timestamps
    table.columns["INTERVAL"] = np.full(rows, interval)
    table.columns["DATA_DESC_ID"] = np.zeros(rows, dtype=int)
    table.columns["SCAN_NUMBER"] = np.where(
        np.arange(rows) < rows // 2, 1, 2
    )
    return table


class FullDatasetMeasurementTests(unittest.TestCase):

    def measure(self, table, rows_wanted):
        with mock.patch.object(SELFCAL, "tb", table, create=True):
            return manager().measure_coherence(
                "target.ms", rows_wanted=rows_wanted, block_rows=150
            )

    def test_thermal_noise_comes_from_all_the_data_not_the_sample(self):
        table = field()
        sampled = self.measure(table, rows_wanted=600)
        whole = self.measure(field(), rows_wanted=10 ** 6)

        self.assertLess(sampled["rows_sampled"], sampled["rows"])
        # The sample's own uncertainty is larger than the floor of the whole
        # dataset, by about the square root of the fraction read.
        fraction = sampled["rows_sampled"] / float(sampled["rows"])
        self.assertAlmostEqual(
            sampled["thermal_noise_Jy_beam"],
            sampled["vector_mean_error"] * math.sqrt(fraction),
            delta=1e-12,
        )
        # And it agrees with reading everything.
        self.assertAlmostEqual(
            sampled["thermal_noise_Jy_beam"] / whole["thermal_noise_Jy_beam"],
            1.0,
            delta=0.05,
        )

    def test_the_shape_a_gain_solve_will_see_is_recorded(self):
        result = self.measure(field(), rows_wanted=10 ** 6)

        self.assertEqual(result["antennas_with_data"], 6)
        self.assertEqual(result["spectral_windows"], 1)
        self.assertEqual(result["scans"], 2)
        self.assertEqual(result["integration_s"], 8.0)
        self.assertEqual(result["integrations"], 200)
        self.assertEqual(result["on_source_s"], 1600.0)
        self.assertEqual(result["correlations"], 2)

    def test_integrations_are_counted_a_block_at_a_time(self):
        table = field()
        with mock.patch.object(SELFCAL, "tb", table, create=True):
            count = SELFCAL.SelfCalManager.count_integrations(
                table.nrows(), block_rows=97
            )
        self.assertEqual(count, 200)
        time_reads = [n for name, _, n in table.reads if name == "TIME"]
        self.assertLessEqual(max(time_reads), 97)


def solvable_manager(amplitude, gaintype="G", **coherence):
    value = manager(phase_gaintype=gaintype)
    value.band_config = {"min_snr_phase": 5.0}
    value.coherence = dict({
        "component_sigma": 0.274,
        "estimated_unflagged_samples": 1.5e10,
        "antennas_with_data": 30,
        "spectral_windows": 1,
        "correlations": 2,
        "scans": 6,
        "integration_s": 2.68,
        "on_source_s": 10800.0,
    }, **coherence)
    value.clean_diagnostics = {
        "initial": {"model_rms_visibility_Jy": amplitude}
    }
    return value


class SolutionSnrTests(unittest.TestCase):

    def test_the_estimate_is_the_textbook_expression(self):
        value = solvable_manager(0.007)
        # One antenna, one hand (G), one 4-minute interval.
        samples = 1.5e10 * (2 / 30.0) * (1 / 2.0) * (240.0 / 10800.0)
        expected = 0.007 * math.sqrt(samples) / 0.274

        self.assertAlmostEqual(value.phase_solution_snr("4min"), expected)
        # Combining the two hands (T) gains sqrt(2).
        self.assertAlmostEqual(
            solvable_manager(0.007, "T").phase_solution_snr("4min"),
            expected * math.sqrt(2),
        )

    def test_solints_are_understood_and_kept_within_a_scan(self):
        value = solvable_manager(0.007)
        self.assertAlmostEqual(
            value.phase_solution_snr("int") / value.phase_solution_snr("4min"),
            math.sqrt(2.68 / 240.0),
        )
        # combine='' never spans scans: an 'inf' or over-long solint is a
        # scan, 10800 s / 6 scans here.
        self.assertAlmostEqual(
            value.phase_solution_snr("inf"),
            value.phase_solution_snr("30min"),
        )
        self.assertAlmostEqual(
            value.phase_solution_snr("5h"), value.phase_solution_snr("inf")
        )
        self.assertIsNone(value.phase_solution_snr("sometimes"))

    def test_missing_ingredients_give_no_estimate(self):
        self.assertIsNone(solvable_manager(None).phase_solution_snr("4min"))
        value = solvable_manager(0.007)
        value.coherence = {"status": "unavailable"}
        self.assertIsNone(value.phase_solution_snr("4min"))


class PhaseStartTests(unittest.TestCase):

    def decide(self, dynamic_range, amplitude, **task_config):
        value = solvable_manager(amplitude)
        value.task_config.update(task_config)
        with mock.patch("builtins.print"):
            started = value.phase_selfcal_is_warranted(
                dynamic_range, ["4min", "2min", "1min", "int"]
            )
        return started, value.phase_start_evidence

    def test_a_field_held_back_by_gain_errors_is_self_calibrated(self):
        # Spokes keep the image at dynamic range 12, but a 7 mJy model gives
        # each antenna's 4-minute solution an S/N near 100.
        started, evidence = self.decide(12.0, 0.007)

        self.assertTrue(started)
        self.assertEqual(evidence["reason"], "solution_snr")
        self.assertGreater(evidence["expected_solution_snr"], 50)

    def test_a_high_dynamic_range_is_still_enough_on_its_own(self):
        started, evidence = self.decide(40.0, None)
        self.assertTrue(started)
        self.assertEqual(evidence["reason"], "dynamic_range")

    def test_a_model_too_weak_to_solve_against_is_not_used(self):
        # 0.1 mJy of model: S/N about 1.2 per solution.
        started, evidence = self.decide(12.0, 1.0e-4)
        self.assertFalse(started)
        self.assertEqual(evidence["reason"], "neither")
        self.assertLess(evidence["expected_solution_snr"], 10)

    def test_the_new_gate_can_be_turned_off(self):
        started, _ = self.decide(
            12.0, 0.007, phase_start_min_solution_snr=0
        )
        self.assertFalse(started)

    def test_the_skip_message_names_both_measurements(self):
        value = solvable_manager(1.0e-4)
        value.image_stats = {
            "dynamic_range": 12.0, "peak_Jy_beam": 1e-3, "rms_Jy_beam": 8e-5
        }
        value.selfcal_strategy = {"phase_selfcal": False}
        with mock.patch("builtins.print"):
            value.phase_selfcal_is_warranted(12.0, ["4min"])
        printed = []
        with mock.patch(
            "builtins.print", lambda *a, **k: printed.append(" ".join(map(str, a)))
        ):
            value.run_phase_selfcal()
        text = " ".join(printed)
        self.assertIn("dynamic range of 12.0", text)
        self.assertIn("S/N 1.2", text)


class RadioStretchTests(unittest.TestCase):

    def test_noise_is_measured_where_the_beam_has_not_amplified_it(self):
        rng = np.random.default_rng(5)
        size = 400
        yy, xx = np.mgrid[:size, :size]
        radius = np.hypot(xx - size / 2, yy - size / 2) / (size / 2)
        beam = np.exp(-4 * np.log(2) * (radius / 1.2) ** 2)
        pbcor = rng.normal(0, 30e-6, (size, size)) / beam

        scaling = QA.radio_display_scaling(pbcor, np)
        self.assertAlmostEqual(scaling["sigma"] / 30e-6, 1.0, delta=0.15)

    def test_sources_and_noise_are_both_inside_the_range(self):
        rng = np.random.default_rng(6)
        image = rng.normal(0, 30e-6, (300, 300))
        image[150, 150] = 3e-3

        # The old stretch: about -2.3 to +2.6 sigma, so the source is white
        # and the display is all noise.
        low, high = np.percentile(image, [1.0, 99.5])
        self.assertLess((high - low) / 30e-6, 5.5)

        scaling = QA.radio_display_scaling(image, np)
        self.assertAlmostEqual(scaling["vmax"], 3e-3)
        self.assertAlmostEqual(scaling["vmin"], -3 * scaling["sigma"])
        self.assertGreater(scaling["dynamic_range"], 50)

    def test_a_blank_image_does_not_break_the_plot(self):
        self.assertIsNone(
            QA.radio_display_scaling(np.full((10, 10), np.nan), np)
        )
        flat = QA.radio_display_scaling(np.zeros((20, 20)), np)
        self.assertIsNone(flat["sigma"])


class NoiseAgainstThermalTests(unittest.TestCase):

    def imaging(self, selfcal):
        value = IMAGING.FinalImagingManager.__new__(IMAGING.FinalImagingManager)
        value.selfcal = selfcal
        value.clean_diagnostics = {"final_continuum": {}}
        return value

    def report(self, selfcal, rms):
        printed = []
        with mock.patch(
            "builtins.print", lambda *a, **k: printed.append(" ".join(map(str, a)))
        ):
            ratio = self.imaging(selfcal).report_noise_against_thermal(
                "final_continuum", rms
            )
        return ratio, " ".join(printed)

    def test_an_error_limited_image_says_so_and_why(self):
        ratio, text = self.report(
            {
                "calibrated_coherence": {"thermal_noise_Jy_beam": 10e-6},
                "selfcal_strategy": {"phase_selfcal": False},
                "phase_start_evidence": {"expected_solution_snr": 3.2},
            },
            150e-6,
        )
        self.assertAlmostEqual(ratio, 15.0)
        self.assertIn("15.0 times", text)
        self.assertIn("limited by errors in the data", text)
        self.assertIn("did not run (expected solution S/N 3.2)", text)

    def test_a_noise_limited_image_is_not_warned_about(self):
        ratio, text = self.report(
            {"calibrated_coherence": {"thermal_noise_Jy_beam": 10e-6}},
            18e-6,
        )
        self.assertAlmostEqual(ratio, 1.8)
        self.assertNotIn("WARNING", text)

    def test_without_a_floor_nothing_is_claimed(self):
        ratio, text = self.report({}, 150e-6)
        self.assertIsNone(ratio)
        self.assertEqual(text, "")


if __name__ == "__main__":
    unittest.main()
