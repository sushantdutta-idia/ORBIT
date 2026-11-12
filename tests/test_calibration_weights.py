"""Calibrated targets carry per-channel weights, and the QA says so."""

import importlib.util
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from test_postcal_flagging import CALIBRATION
from test_selfcal_coherence import SELFCAL, FakeVisibilityTable, manager


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_radio_continuum_qa_products_calibration", SOURCE_DIR / "qa_products.py"
)
QA = importlib.util.module_from_spec(spec)
spec.loader.exec_module(QA)
spec = importlib.util.spec_from_file_location(
    "test_radio_continuum_config_weights", SOURCE_DIR / "config.py"
)
CONFIG = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CONFIG)


def fake_statwt(calls):
    def statwt(**kwargs):
        calls.append(kwargs)
    statwt.casa_task_available = True
    return statwt


class WeightTable:
    """WEIGHT_SPECTRUM and FLAG as CASA shapes them: (corr, chan, row)."""

    def __init__(self, weights, flags=None):
        self.columns = {"WEIGHT_SPECTRUM": weights,
                        "FLAG": np.zeros(weights.shape, dtype=bool)
                        if flags is None else flags}

    def open(self, path):
        return None

    def close(self):
        return None

    def colnames(self):
        return list(self.columns)

    def nrows(self):
        return self.columns["FLAG"].shape[-1]

    def getcol(self, name, startrow=0, nrow=-1):
        return self.columns[name][..., startrow:startrow + nrow]


def calibrator(root, **config):
    value = CALIBRATION.CalibrationManager.__new__(
        CALIBRATION.CalibrationManager
    )
    value.task_config = dict(config)
    value.caldir = Path(root)
    return value


class StatisticalWeightTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.ms = self.root / "TARGET.ms"
        self.ms.mkdir()
        (self.ms / "table.dat").touch()

    def tearDown(self):
        self.temporary.cleanup()

    def test_ugmrt_weights_its_targets_by_default(self):
        ugmrt = CONFIG.TELESCOPE_CALIBRATION_TASK_CONFIG["ugmrt"]
        self.assertTrue(ugmrt["statwt_enabled"])
        self.assertEqual(ugmrt["statwt_chanbin"], 16)
        for telescope in ("meerkat", "lofar"):
            self.assertFalse(
                CONFIG.TELESCOPE_CALIBRATION_TASK_CONFIG[telescope].get(
                    "statwt_enabled", False
                )
            )

    def test_statwt_writes_per_channel_weights_and_a_marker(self):
        value = calibrator(self.root, statwt_enabled=True)
        calls = []
        value.channel_weight_spread = lambda outms: {
            "channels_with_weights": 2800,
            "channel_weight_ratio_p95_p5": 180.0,
            "equal_weight_noise_penalty": 4.2,
        }
        with mock.patch.object(CALIBRATION, "statwt", fake_statwt(calls)), \
                mock.patch("builtins.print"):
            record = value.set_statistical_weights(self.ms, "TARGET")

        self.assertEqual(len(calls), 1)
        call = calls[0]
        self.assertEqual(call["vis"], str(self.ms))
        # The split target's calibrated data are in DATA.
        self.assertEqual(call["datacolumn"], "data")
        # A channel bin, not "spw": that is what creates WEIGHT_SPECTRUM.
        self.assertEqual(call["chanbin"], 16)
        self.assertEqual(call["timebin"], "120s")
        self.assertEqual(call["minsamp"], 16)
        self.assertTrue(call["flagbackup"])

        marker = json.loads(value.statwt_marker(self.ms).read_text())
        self.assertEqual(marker["status"], "applied")
        self.assertEqual(marker["equal_weight_noise_penalty"], 4.2)
        self.assertEqual(record["status"], "applied")

    def test_disabled_does_nothing(self):
        value = calibrator(self.root, statwt_enabled=False)
        calls = []
        with mock.patch.object(CALIBRATION, "statwt", fake_statwt(calls)):
            record = value.set_statistical_weights(self.ms, "TARGET")
        self.assertEqual(calls, [])
        self.assertEqual(record["status"], "disabled")
        self.assertFalse(value.statwt_marker(self.ms).exists())

    def test_a_half_weighted_target_is_split_again_not_reused(self):
        # An interrupted run left the target MS but no statwt marker.
        value = calibrator(self.root, statwt_enabled=True)
        value.splitdir = self.root
        value.ms = str(self.root / "input.ms")
        value.checkpoint_exists = lambda name: False
        value.get_targets = lambda: ["TARGET"]
        value.write_atomic_json = lambda path, data: Path(path).write_text(
            json.dumps(data)
        )
        value.mark_checkpoint = lambda name: None
        value.channel_weight_spread = lambda outms: None
        splits, weights = [], []

        def split(**kwargs):
            splits.append(kwargs)
            out = Path(kwargs["outputvis"])
            out.mkdir()
            (out / "table.dat").touch()

        with mock.patch.object(CALIBRATION, "split", split), \
                mock.patch.object(CALIBRATION, "statwt",
                                  fake_statwt(weights)), \
                mock.patch.object(
                    CALIBRATION.IMAGING_SAFETY,
                    "measurement_set_instance_identity",
                    lambda path: {"id": 1},
                ), mock.patch("builtins.print"):
            value.split_targets()

        self.assertEqual(len(splits), 1)
        self.assertEqual(len(weights), 1)
        self.assertTrue(
            (self.root / "TARGET.ms.incomplete.1").is_dir(),
            "the half-weighted MS must be preserved, not overwritten",
        )
        self.assertTrue(value.statwt_marker(self.ms).is_file())

    def test_a_missing_statwt_is_an_error_not_a_silent_skip(self):
        value = calibrator(self.root, statwt_enabled=True)
        with self.assertRaisesRegex(RuntimeError, "statwt_enabled: false"):
            value.set_statistical_weights(self.ms, "TARGET")

    def test_the_equal_weight_penalty_is_measured_from_the_weights(self):
        # Channel noise rising 20x across the band: weights fall 400x.
        channels = 64
        sigma = np.linspace(1.0, 20.0, channels)
        weights = np.broadcast_to(
            (1.0 / sigma ** 2)[np.newaxis, :, np.newaxis], (2, channels, 500)
        ).copy()
        flags = np.zeros(weights.shape, dtype=bool)
        flags[:, 10, :] = True
        with mock.patch.object(CALIBRATION, "table",
                               lambda: WeightTable(weights, flags)):
            spread = CALIBRATION.CalibrationManager.channel_weight_spread(
                self.ms, rows_wanted=200, block_rows=50
            )

        kept = np.delete(sigma, 10)
        expected = math.sqrt(np.mean(kept ** 2) * np.mean(1.0 / kept ** 2))
        self.assertAlmostEqual(spread["equal_weight_noise_penalty"], expected)
        self.assertEqual(spread["channels_with_weights"], channels - 1)
        low, high = np.percentile(1.0 / kept ** 2, [5.0, 95.0])
        self.assertAlmostEqual(
            spread["channel_weight_ratio_p95_p5"], high / low
        )


class PlotDataTests(unittest.TestCase):
    """What the calibration plots draw is what the tables hold."""

    def test_flagged_solutions_are_not_drawn(self):
        gains = np.full((2, 1, 6), 0.55 + 0j)
        gains[:, :, 2] = 1.0 + 0j  # CASA's value for a flagged solution
        flags = np.zeros(gains.shape, dtype=bool)
        flags[:, :, 2] = True
        series = QA._time_series(gains, flags)

        self.assertEqual(series.shape, (2, 6))
        self.assertTrue(np.all(np.isnan(series[:, 2])))
        self.assertTrue(np.allclose(np.abs(series[:, [0, 1, 3, 4, 5]]), 0.55))

    def test_flagged_channels_are_gaps_in_the_bandpass(self):
        bandpass = np.full((2, 8, 30), 1.2 + 0j)
        flags = np.zeros(bandpass.shape, dtype=bool)
        flags[:, 5:, :] = True
        bandpass[:, 5:, :] = 1.79  # a fill value, not a solution
        series = QA._bandpass_series(bandpass, "amplitude", flags)

        self.assertEqual(series.shape, (2, 3, 8))
        self.assertTrue(np.allclose(series[:, 0, :5], 1.2))
        self.assertTrue(np.all(np.isnan(series[:, 0, 5:])))

    def test_phase_spread_is_circular(self):
        phases = np.deg2rad(np.array([179.0, -179.0, 178.0, -178.0]))
        bandpass = np.exp(1j * phases)[np.newaxis, np.newaxis, :].repeat(
            2, axis=1
        )
        series = QA._bandpass_series(bandpass, "phase")
        self.assertAlmostEqual(abs(series[0, 0, 0]), 180.0, delta=0.5)
        self.assertLess(series[0, 2, 0] - series[0, 0, 0], 3.0)


class PhaseStabilityTests(unittest.TestCase):

    def test_changes_between_visits_are_wrapped_and_measured(self):
        # Two antennas, four visits 55 min apart; antenna 0 is the reference.
        times = np.repeat(np.arange(4) * 3300.0, 2)
        antennas = np.tile([0, 1], 4)
        fields = np.zeros(8, dtype=int)
        phase_1 = [170.0, -170.0, -150.0, -120.0]  # +20, +20, +30 wrapped
        phases = np.zeros((1, 8))
        phases[0, 1::2] = phase_1

        result = QA.phase_transfer_stability(times, antennas, fields, phases)
        field = result[0]
        self.assertEqual(field["transitions"], 3)
        self.assertAlmostEqual(field["median_gap_min"], 55.0)
        self.assertAlmostEqual(field["median_abs_change_deg"], 20.0)
        self.assertEqual(field["worst_antennas"][0]["antenna"], 1)

    def test_flagged_solutions_do_not_count(self):
        times = np.arange(3) * 3300.0
        phases = np.array([[10.0, 90.0, 20.0]])
        flags = np.array([[False, True, False]])
        result = QA.phase_transfer_stability(
            times, np.ones(3, dtype=int), np.zeros(3, dtype=int), phases, flags
        )
        self.assertAlmostEqual(result[0]["median_abs_change_deg"], 10.0)


class StabilitySummaryTests(unittest.TestCase):
    """The summary survives the names CASA actually returns."""

    def tables(self):
        # Four antennas, three visits; antenna 2's gain is a third of the
        # others'. The table tool returns names as numpy arrays.
        antennas = np.tile(np.arange(4), 3)
        times = np.repeat(np.arange(3) * 3300.0, 4)
        gains = np.full((1, 1, 12), 0.55 + 0j)
        gains[0, 0, antennas == 2] = 0.18
        return {
            "parameters": gains,
            "time": times,
            "antenna": antennas,
            "field": np.ones(12, dtype=int),
            "flags": np.zeros((1, 1, 12), dtype=bool),
            "antenna_names": np.array(["C00", "C01", "C02", "C04"]),
            "field_names": np.array(["3C48", "0323+055"]),
        }

    def test_numpy_name_arrays_do_not_lose_the_summary(self):
        with mock.patch.object(QA, "_read_caltable", lambda path: self.tables()):
            summary = QA.caltable_stability_summary("phase_scan.cal", "flux.cal")

        self.assertNotIn("phase_transfer_error", summary)
        self.assertNotIn("antenna_gains_error", summary)
        self.assertEqual(summary["field_names"], ["3C48", "0323+055"])
        outliers = summary["antenna_gains"]["outliers"]
        self.assertEqual([item["name"] for item in outliers], ["C02"])
        self.assertAlmostEqual(outliers[0]["ratio_to_array"], 0.18 / 0.55)

    def test_missing_names_are_tolerated(self):
        tables = self.tables()
        tables["antenna_names"] = None
        tables["field_names"] = None
        with mock.patch.object(QA, "_read_caltable", lambda path: tables):
            summary = QA.caltable_stability_summary("phase_scan.cal", "flux.cal")

        self.assertEqual(summary["field_names"], [])
        self.assertIsNone(summary["antenna_gains"]["outliers"][0]["name"])


class AntennaGainTests(unittest.TestCase):

    def test_a_low_gain_chain_is_named(self):
        antennas = np.repeat(np.arange(10), 5)
        gains = np.full((2, 50), 0.55)
        gains[1, antennas == 7] = 0.18
        result = QA.antenna_gain_outliers(gains, antennas)

        self.assertAlmostEqual(result["array_median"], 0.55)
        self.assertEqual(len(result["outliers"]), 1)
        item = result["outliers"][0]
        self.assertEqual((item["antenna"], item["correlation"]), (7, 1))
        self.assertAlmostEqual(item["ratio_to_array"], 0.18 / 0.55)

    def test_flagged_gains_are_ignored(self):
        antennas = np.repeat(np.arange(4), 3)
        gains = np.full((1, 12), 0.55)
        gains[0, 0] = 1.0
        flags = np.zeros((1, 12), dtype=bool)
        flags[0, 0] = True
        result = QA.antenna_gain_outliers(gains, antennas, flags)
        self.assertEqual(result["outliers"], [])


class WeightedThermalNoiseTests(unittest.TestCase):

    def weighted_field(self, weight=4.0):
        rng = np.random.default_rng(2)
        rows = 400
        shape = (2, 4, rows)
        data = rng.normal(size=shape) + 1j * rng.normal(size=shape)
        table = FakeVisibilityTable(
            data, flags=np.zeros(shape, dtype=bool),
            antenna1=np.zeros(rows, dtype=int),
            antenna2=np.ones(rows, dtype=int),
        )
        table.columns["WEIGHT_SPECTRUM"] = np.full(shape, weight)
        return table

    def test_statistical_weights_set_the_floor(self):
        table = self.weighted_field(weight=4.0)
        value = manager()
        value.target = "TARGET"
        value.calibration = {
            "statistical_weights": {"TARGET": {"status": "applied"}}
        }
        with mock.patch.object(SELFCAL, "tb", table, create=True):
            measurement = value.measure_coherence(
                "target.ms", rows_wanted=10 ** 6, block_rows=100
            )
        value.adopt_statistical_thermal_noise(measurement)

        samples = 2 * 4 * 400
        self.assertAlmostEqual(
            measurement["thermal_noise_Jy_beam"],
            1.0 / math.sqrt(4.0 * samples),
        )
        self.assertEqual(measurement["thermal_noise_source"],
                         "statistical_weights")
        self.assertIsNotNone(
            measurement["thermal_noise_equal_weights_Jy_beam"]
        )

    def test_row_weights_that_are_not_statistical_are_not_trusted(self):
        table = self.weighted_field()
        value = manager()
        value.target = "TARGET"
        value.calibration = {}
        with mock.patch.object(SELFCAL, "tb", table, create=True):
            measurement = value.measure_coherence("target.ms")
        before = measurement["thermal_noise_Jy_beam"]
        self.assertIsNone(value.adopt_statistical_thermal_noise(measurement))
        self.assertEqual(measurement["thermal_noise_Jy_beam"], before)

    def test_a_bright_field_is_warned_about(self):
        value = manager()
        value.target = "TARGET"
        value.calibration = {
            "statistical_weights": {"TARGET": {"status": "applied"}}
        }
        value.coherence = {"component_sigma": 0.274}
        value.clean_diagnostics = {
            "initial": {"model_rms_visibility_Jy": 0.007}
        }
        printed = []
        with mock.patch(
            "builtins.print", lambda *a, **k: printed.append(" ".join(map(str, a)))
        ):
            self.assertAlmostEqual(value.check_weights_against_sky(),
                                   0.007 / 0.274)
            self.assertEqual(printed, [])
            value.clean_diagnostics["initial"]["model_rms_visibility_Jy"] = 0.1
            value.check_weights_against_sky()
        self.assertIn("under-weighted", " ".join(printed))


if __name__ == "__main__":
    unittest.main()
