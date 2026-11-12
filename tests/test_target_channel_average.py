"""A calibrated target is written at the channel width imaging needs.

The uGMRT Band 4 run that prompted this kept the correlator's 4096 channels
of 97.7 kHz in every target. Each self-calibration image then read and
gridded all of them in every major cycle: seven minutes a cycle, about five
hours an image, and more than a day of self-calibration for one 33-minute
target.
"""

import importlib.util
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_postcal_flagging import CALIBRATION
from test_calibration_weights import calibrator, fake_statwt


SOURCE_DIR = Path(__file__).resolve().parents[1]
NATIVE_HZ = 400.0e6 / 4096


def load_config(name):
    spec = importlib.util.spec_from_file_location(name, SOURCE_DIR / "config.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SpectralWindowTable:
    """SPECTRAL_WINDOW as the table tool reads it, one row per window."""

    def __init__(self, windows, fail=False):
        self.windows = windows
        self.fail = fail
        self.opened = None

    def open(self, path):
        if self.fail:
            raise RuntimeError("table does not exist")
        self.opened = path

    def close(self):
        self.opened = None

    def nrows(self):
        return len(self.windows)

    def getcell(self, column, row):
        count, width = self.windows[row]
        if column == "NUM_CHAN":
            return count
        if column == "CHAN_WIDTH":
            return [width] * count
        raise KeyError(column)


def bare_ms(path):
    path.mkdir()
    (path / "table.dat").touch()


def smearing_peak_loss(channel_hz, radius_arcsec, beam_arcsec, nu_hz):
    """Bridle & Schwab: Gaussian beam, rectangular bandpass, radial smear."""
    beta = (channel_hz / nu_hz) * (radius_arcsec / beam_arcsec)
    return 1.0 - 1.0645 / beta * math.erf(0.8326 * beta)


class ChannelAveragePlanTests(unittest.TestCase):

    def plan(self, limit_khz, windows, fail=False):
        value = calibrator("/tmp")
        value.ms = "calibrated.ms"
        value.band_config = (
            {} if limit_khz is None
            else {"target_max_channel_width_khz": limit_khz}
        )
        table = SpectralWindowTable(windows, fail=fail)
        printed = []
        with mock.patch.object(CALIBRATION, "table", lambda: table), \
                mock.patch("builtins.print",
                           lambda *a, **k: printed.append(" ".join(map(str, a)))):
            plan = value.target_channel_average_plan()
        return plan, printed

    def test_band_4_averages_four_channels(self):
        plan, _ = self.plan(400.0, [(4096, NATIVE_HZ)])

        self.assertEqual(plan["width"], 4)
        self.assertAlmostEqual(plan["native_channel_width_khz"], 97.65625)
        self.assertAlmostEqual(plan["output_channel_width_khz"], 390.625)
        self.assertEqual(plan["output_channels"], [1024])

    def test_the_smearing_quoted_for_band_4_holds(self):
        # 5120 pixels of 1 arcsec, a 4.9 x 4.5 arcsec beam at 750 MHz. The
        # ratio of smear to beam does not depend on the channel's frequency,
        # because both scale as 1/frequency.
        beam = math.sqrt(4.88 * 4.51)
        edge = smearing_peak_loss(390.625e3, 2560.0, beam, 750.0e6)
        inner = smearing_peak_loss(390.625e3, 1800.0, beam, 750.0e6)
        self.assertLess(edge, 0.02)
        self.assertLess(inner, 0.01)

    def test_a_band_without_the_setting_is_left_alone(self):
        plan, _ = self.plan(None, [(4096, NATIVE_HZ)])
        self.assertEqual(plan["width"], 1)
        self.assertNotIn("native_channel_width_khz", plan)

    def test_zero_switches_it_off(self):
        plan, _ = self.plan(0.0, [(4096, NATIVE_HZ)])
        self.assertEqual(plan["width"], 1)

    def test_channels_already_wide_enough_are_not_averaged(self):
        plan, _ = self.plan(150.0, [(4096, NATIVE_HZ)])
        self.assertEqual(plan["width"], 1)

    def test_the_factor_divides_every_window(self):
        # 8 would fit the width, but 990 channels do not divide by 8, 7 or 6.
        plan, _ = self.plan(800.0, [(1000, 100.0e3), (990, 100.0e3)])
        self.assertEqual(plan["width"], 5)
        self.assertEqual(plan["output_channels"], [200, 198])

    def test_the_widest_window_sets_the_factor(self):
        plan, _ = self.plan(400.0, [(4096, NATIVE_HZ), (2048, 2 * NATIVE_HZ)])
        self.assertEqual(plan["width"], 2)

    def test_unreadable_widths_keep_native_channels_and_say_so(self):
        plan, printed = self.plan(400.0, [], fail=True)
        self.assertEqual(plan["width"], 1)
        self.assertEqual(plan["status"], "channel_widths_unreadable")
        self.assertIn("keep their native channels", " ".join(printed))


class SplitTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.outms = self.root / "TARGET.ms"
        self.native = self.root / "TARGET.native.ms"

    def tearDown(self):
        self.temporary.cleanup()

    def calibrator(self, limit_khz, statwt_enabled=True):
        value = calibrator(self.root, statwt_enabled=statwt_enabled)
        value.splitdir = self.root
        value.ms = str(self.root / "calibrated.ms")
        value.band_config = {"target_max_channel_width_khz": limit_khz}
        value.checkpoint_exists = lambda name: False
        value.get_targets = lambda: ["TARGET"]
        value.write_atomic_json = lambda path, data: Path(path).write_text(
            json.dumps(data)
        )
        value.mark_checkpoint = lambda name: None
        value.channel_weight_spread = lambda outms: None
        return value

    def run_split(self, value):
        splits, weights = [], []
        seen_native = []

        def split(**kwargs):
            splits.append(kwargs)
            if kwargs["vis"] == str(self.native):
                # The averaging split reads a weighted, complete native copy.
                seen_native.append(
                    value.statwt_marker(self.native).is_file()
                    if value.statistical_weights_settings()["enabled"]
                    else True
                )
            bare_ms(Path(kwargs["outputvis"]))

        def statwt(**kwargs):
            weights.append(kwargs)
            Path(kwargs["vis"] + ".flagversions").mkdir()

        statwt.casa_task_available = True
        table = SpectralWindowTable([(4096, NATIVE_HZ)])
        with mock.patch.object(CALIBRATION, "split", split), \
                mock.patch.object(CALIBRATION, "statwt", statwt), \
                mock.patch.object(CALIBRATION, "table", lambda: table), \
                mock.patch.object(
                    CALIBRATION.IMAGING_SAFETY,
                    "measurement_set_instance_identity",
                    lambda path: {"id": 1},
                ), mock.patch("builtins.print"):
            value.split_targets()
        return splits, weights, seen_native

    def test_band_4_targets_are_weighted_then_averaged(self):
        value = self.calibrator(400.0)
        splits, weights, seen_native = self.run_split(value)

        self.assertEqual(len(splits), 2)
        native, averaged = splits
        self.assertEqual(native["vis"], value.ms)
        self.assertEqual(native["field"], "TARGET")
        self.assertEqual(native["datacolumn"], "corrected")
        self.assertEqual(native["outputvis"], str(self.native))

        # statwt measures the native channels, with its configured bin.
        self.assertEqual(len(weights), 1)
        self.assertEqual(weights[0]["vis"], str(self.native))
        self.assertEqual(weights[0]["chanbin"], 16)

        self.assertEqual(averaged["vis"], str(self.native))
        self.assertEqual(averaged["outputvis"], str(self.outms))
        self.assertEqual(averaged["datacolumn"], "data")
        self.assertEqual(averaged["width"], 4)
        self.assertTrue(averaged["keepflags"])
        self.assertEqual(seen_native, [True])

        marker = json.loads(value.statwt_marker(self.outms).read_text())
        self.assertEqual(marker["status"], "applied")
        self.assertEqual(marker["channel_average"]["width"], 4)
        self.assertAlmostEqual(
            marker["channel_average"]["output_channel_width_khz"], 390.625
        )

        # The full-resolution copy and what statwt left with it are gone.
        self.assertFalse(self.native.exists())
        self.assertFalse(Path(str(self.native) + ".flagversions").exists())
        self.assertFalse(value.statwt_marker(self.native).exists())

        plan = json.loads(value.target_channel_average_path().read_text())
        self.assertEqual(plan["width"], 4)
        self.assertEqual(
            json.loads((self.root / "split_targets.json").read_text()),
            {"TARGET": str(self.outms)},
        )

    def test_without_the_setting_the_split_is_unchanged(self):
        value = self.calibrator(None)
        value.band_config = {}
        splits, weights, _ = self.run_split(value)

        self.assertEqual(len(splits), 1)
        self.assertEqual(splits[0]["outputvis"], str(self.outms))
        self.assertNotIn("width", splits[0])
        self.assertEqual(weights[0]["vis"], str(self.outms))
        self.assertFalse(value.target_channel_average_path().exists())

    def test_an_interrupted_split_is_redone_from_the_start(self):
        # A run stopped after the native split, before the averaged MS was
        # weighted and marked, left both behind.
        bare_ms(self.native)
        bare_ms(self.outms)
        value = self.calibrator(400.0)
        splits, _, _ = self.run_split(value)

        self.assertEqual(len(splits), 2)
        self.assertTrue(
            (self.root / "TARGET.ms.incomplete.1").is_dir(),
            "an unmarked target MS is preserved, as before",
        )
        self.assertTrue(value.statwt_marker(self.outms).is_file())
        self.assertFalse(self.native.exists())

    def test_a_weighted_averaged_target_is_reused(self):
        value = self.calibrator(400.0)
        bare_ms(self.outms)
        value.statwt_marker(self.outms).write_text(json.dumps(
            {"status": "applied"}
        ))
        splits, weights, _ = self.run_split(value)

        self.assertEqual(splits, [])
        self.assertEqual(weights, [])

    def test_without_statwt_the_target_is_still_averaged(self):
        value = self.calibrator(400.0, statwt_enabled=False)
        splits, weights, _ = self.run_split(value)

        self.assertEqual([item.get("width") for item in splits], [None, 4])
        self.assertEqual(weights, [])
        self.assertFalse(value.statwt_marker(self.outms).exists())
        self.assertFalse(self.native.exists())

    def test_only_a_native_copy_can_be_discarded(self):
        bare_ms(self.outms)
        with self.assertRaisesRegex(ValueError, "native-resolution"):
            CALIBRATION.CalibrationManager.discard_native_target(self.outms)
        self.assertTrue(self.outms.is_dir())


class ConfigurationTests(unittest.TestCase):

    def test_only_band_4_is_averaged_by_default(self):
        config = load_config("test_radio_continuum_config_channel_default")
        self.assertEqual(
            config.CAL_CONFIG["B4"]["target_max_channel_width_khz"], 400.0
        )
        for band, settings in config.CAL_CONFIG.items():
            if band != "B4":
                with self.subTest(band=band):
                    self.assertNotIn("target_max_channel_width_khz", settings)

    def test_the_fallback_table_agrees(self):
        self.assertEqual(
            CALIBRATION.CalibrationManager.CAL_CONFIG["B4"][
                "target_max_channel_width_khz"
            ],
            400.0,
        )

    def test_config_yaml_sets_every_band(self):
        config = load_config("test_radio_continuum_config_channel_yaml")
        config.CONFIG_YAML.clear()
        config.CONFIG_YAML["calibration_target_max_channel_width_khz"] = 0
        config._apply_config_yaml_overrides()
        for band, settings in config.CAL_CONFIG.items():
            with self.subTest(band=band):
                self.assertEqual(settings["target_max_channel_width_khz"], 0.0)

    def test_config_yaml_rejects_a_negative_width(self):
        config = load_config("test_radio_continuum_config_channel_bad")
        config.CONFIG_YAML.clear()
        config.CONFIG_YAML["calibration_target_max_channel_width_khz"] = -1
        with self.assertRaisesRegex(ValueError, "0 or more"):
            config._apply_config_yaml_overrides()


if __name__ == "__main__":
    unittest.main()
