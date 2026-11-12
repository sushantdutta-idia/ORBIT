"""Amplitude self-calibration is kept only when the image says it helped.

Before this, an enabled amp+phase solve was applied and committed without a
look at the result. Phase loops have always been made to earn their place
with an image; the amplitude solve, which can also move the flux scale, was
the one step that was not.
"""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_selfcal_resume_safety import load_selfcal_module


SOURCE_DIR = Path(__file__).resolve().parents[1]
SELFCAL = load_selfcal_module()

# The phase-only image of the uGMRT target that prompted this.
REFERENCE = {"max": 23.06e-3, "rms": 86.4e-6, "dynamic_range": 267.0}


def image(dynamic_range, peak=23.3e-3):
    return {
        "max": peak,
        "rms": peak / dynamic_range,
        "dynamic_range": dynamic_range,
    }


class ColumnTable:

    def __init__(self, columns=("DATA", "MODEL_DATA", "CORRECTED_DATA")):
        self.columns = list(columns)

    def open(self, path):
        return None

    def close(self):
        return None

    def colnames(self):
        return list(self.columns)


class Stop(Exception):
    pass


class AmpSelfcalTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        self.msdir = root / "ms"
        self.imagedir = root / "images"
        self.caltabledir = root / "caltables"
        self.checkdir = root / "checkpoints"
        self.reportdir = root / "reports"
        for directory in (self.msdir, self.imagedir, self.caltabledir,
                          self.checkdir, self.reportdir):
            directory.mkdir()
        self.phase_ms = self.msdir / "phase3.ms"
        self.phase_ms.mkdir()
        self.candidate = self.msdir / "amp_phase.candidate.ms"
        self.final = self.msdir / "final_selfcal.ms"

        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.msdir = self.msdir
        manager.imagedir = self.imagedir
        manager.caltabledir = self.caltabledir
        manager.checkdir = self.checkdir
        manager.reportdir = self.reportdir
        manager.task_config = {
            "enable_amp_selfcal": True,
            "amp_start_dynamic_range": 100.0,
            "phase_improvement_stop": 0.02,
            "applymode": "calonly",
            "calwt": False,
            "parang": False,
            "minblperant": 4,
            "gain_minsnr": 3.0,
            "amp_gaintype": "G",
        }
        manager.band_config = {"amp_solint": "inf", "min_snr_amp": 40.0}
        manager.selfcal_strategy = {"amp_solint": "inf", "amp_selfcal": True}
        manager.phase_tables = ["phase_1.cal", "phase_2.cal", "phase_3.cal"]
        manager.current_ms = str(self.phase_ms)
        manager.reference_antenna = "C00"
        manager.phase_reference_antennas = lambda ms: "C00,W05,W03"
        manager.ensure_model_column = mock.Mock()
        manager.prepare_caltable = lambda caltable: False
        manager.caltable_status = lambda caltable: (True, "")
        self.manager = manager

    def tearDown(self):
        self.temporary.cleanup()

    def run_amp(self, candidate_image, reference=None, imaging_error=None):
        calls = {"gaincal": [], "split": [], "applycal": [], "imaged": []}

        def split(**kwargs):
            calls["split"].append(kwargs)
            Path(kwargs["outputvis"]).mkdir()

        def image_candidate(label):
            calls["imaged"].append((label, self.manager.current_ms))
            if imaging_error:
                raise imaging_error
            return candidate_image

        self.manager.image_candidate = image_candidate
        with mock.patch.object(
            SELFCAL, "gaincal", lambda **kw: calls["gaincal"].append(kw)
        ), mock.patch.object(SELFCAL, "split", split), mock.patch.object(
            SELFCAL, "applycal", lambda **kw: calls["applycal"].append(kw)
        ), mock.patch.object(
            SELFCAL, "tb", ColumnTable(), create=True
        ), mock.patch.object(
            SELFCAL.CASA_COMPAT, "peak_and_noise",
            lambda imstat, path: dict(reference or REFERENCE),
        ), mock.patch("builtins.print"):
            self.manager.run_amp_phase_selfcal()
        return calls

    def rejection(self):
        return json.loads(
            (self.checkdir / "amp_selfcal.rejected.json").read_text()
        )

    def test_a_solution_that_improves_the_image_is_kept(self):
        calls = self.run_amp(image(300.0))

        # The solve is against the last accepted phase loop's model.
        self.manager.ensure_model_column.assert_called_once_with(
            str(self.phase_ms), self.imagedir / "phase_3"
        )
        solve = calls["gaincal"][0]
        self.assertEqual(solve["calmode"], "ap")
        self.assertEqual(solve["solint"], "inf")
        # Not the band's min_snr_amp of 40, which flagged 74 of 81 solutions
        # on the field this was first tried on.
        self.assertEqual(solve["minsnr"], 10.0)
        # Both polarizations together, for the signal-to-noise.
        self.assertEqual(solve["gaintype"], "T")
        # A list: with one reference antenna and no alternate, CASA flags
        # every antenna wherever that antenna's own solution fails.
        self.assertEqual(solve["refant"], "C00,W05,W03")
        self.assertEqual(solve["refantmode"], "flex")
        # Normalised, so the flux scale stays the calibrators'.
        self.assertTrue(solve["solnorm"])
        self.assertEqual(solve["normtype"], "median")

        # Applied to a copy, which is what gets imaged.
        copy, commit = calls["split"]
        self.assertEqual(copy["vis"], str(self.phase_ms))
        self.assertEqual(copy["outputvis"], str(self.candidate))
        self.assertEqual(copy["datacolumn"], "all")
        self.assertEqual(calls["applycal"][0]["vis"], str(self.candidate))
        self.assertEqual(calls["applycal"][0]["applymode"], "calonly")
        # Each visibility takes its own scan's solution. Interpolating
        # linearly towards a failed scan's solution leaves the half of the
        # scan before it that faces it uncorrected.
        self.assertEqual(calls["applycal"][0]["interp"], ["nearest"])
        self.assertEqual(calls["imaged"], [("amp_phase", str(self.candidate))])

        self.assertEqual(commit["vis"], str(self.candidate))
        self.assertEqual(commit["outputvis"], str(self.final))
        self.assertEqual(commit["datacolumn"], "corrected")
        self.assertFalse(self.candidate.exists())
        self.assertEqual(self.manager.current_ms, str(self.final))
        self.assertTrue(self.manager.amp_table.endswith("amp_phase.cal"))
        self.assertTrue((self.checkdir / "amp_selfcal.done").exists())

        record = json.loads((self.reportdir / "amp_selfcal.json").read_text())
        self.assertTrue(record["accepted"])
        self.assertAlmostEqual(record["improvement"], 33.0 / 267.0)
        self.assertEqual(record["measurement_set"], str(self.final))

    def test_the_threshold_and_gain_type_can_be_configured(self):
        self.manager.task_config["amp_minsnr"] = 25.0
        self.manager.task_config["amp_solve_gaintype"] = "G"
        calls = self.run_amp(image(300.0))
        self.assertEqual(calls["gaincal"][0]["minsnr"], 25.0)
        self.assertEqual(calls["gaincal"][0]["gaintype"], "G")

    def test_a_solve_with_no_usable_table_is_recorded_not_applied(self):
        self.manager.caltable_status = lambda caltable: (
            False, "unusable calibration-solution fraction 0.913580"
        )
        calls = self.run_amp(image(300.0))

        self.assertEqual(calls["split"], [])
        self.assertEqual(calls["applycal"], [])
        self.assertEqual(self.manager.current_ms, str(self.phase_ms))
        history = self.manager.amp_history
        self.assertFalse(history["accepted"])
        self.assertTrue(history["attempted"])
        self.assertIn("0.913580", history["rejection_reason"])
        self.assertEqual(history["solve_parameters"]["minsnr"], 10.0)
        # It is not a verdict on the image, so a later attempt may solve
        # again.
        self.assertFalse((self.checkdir / "amp_selfcal.rejected.json").exists())

    def test_an_accepted_solve_makes_an_earlier_final_image_stale(self):
        # The stage was run before, phase-only, and made its final image.
        (self.checkdir / "final_image.done").touch()
        self.run_amp(image(300.0))
        self.assertFalse((self.checkdir / "final_image.done").exists())

    def test_a_rejected_solve_leaves_an_earlier_final_image_in_place(self):
        (self.checkdir / "final_image.done").touch()
        self.run_amp(image(268.0))
        self.assertTrue((self.checkdir / "final_image.done").exists())

    def test_a_solution_that_does_not_improve_the_image_is_discarded(self):
        calls = self.run_amp(image(268.0))

        self.assertEqual(len(calls["split"]), 1)
        self.assertFalse(self.candidate.exists())
        self.assertFalse(self.final.exists())
        self.assertEqual(self.manager.current_ms, str(self.phase_ms))
        self.assertFalse(hasattr(self.manager, "amp_table"))
        self.assertFalse((self.checkdir / "amp_selfcal.done").exists())
        record = self.rejection()
        self.assertFalse(record["accepted"])
        self.assertIn("Dynamic range", record["rejection_reason"])
        self.assertEqual(self.manager.amp_history["accepted"], False)

    def test_a_solution_that_moves_the_flux_scale_is_discarded(self):
        # A fifth more dynamic range, bought with a peak 15% higher.
        self.run_amp(image(320.0, peak=23.06e-3 * 1.15))

        self.assertEqual(self.manager.current_ms, str(self.phase_ms))
        self.assertFalse(self.final.exists())
        record = self.rejection()
        self.assertIn("flux scale", record["rejection_reason"])
        self.assertAlmostEqual(record["peak_change"], 0.15)

    def test_the_dynamic_range_test_is_made_after_the_phase_loops(self):
        # The initial image was below the threshold, so no amplitude interval
        # was chosen then; the phase-calibrated image is well above it.
        self.manager.selfcal_strategy = {"amp_solint": None, "amp_selfcal": True}
        calls = self.run_amp(image(300.0))

        self.assertEqual(calls["gaincal"][0]["solint"], "inf")
        self.assertEqual(self.manager.current_ms, str(self.final))

    def test_a_poorly_modelled_field_is_not_attempted(self):
        calls = self.run_amp(
            image(300.0),
            reference={"max": 5.0e-3, "rms": 8.0e-5, "dynamic_range": 62.5},
        )

        self.assertEqual(calls["gaincal"], [])
        self.assertEqual(calls["split"], [])
        record = self.rejection()
        self.assertFalse(record["attempted"])
        self.assertIn("below the 100", record["rejection_reason"])

    def test_a_failure_while_imaging_leaves_the_phase_only_data_in_use(self):
        with self.assertRaises(Stop):
            self.run_amp(image(300.0), imaging_error=Stop())

        self.assertFalse(self.candidate.exists())
        self.assertEqual(self.manager.current_ms, str(self.phase_ms))
        self.assertFalse((self.checkdir / "amp_selfcal.done").exists())
        self.assertFalse((self.checkdir / "amp_selfcal.rejected.json").exists())

    def test_a_recorded_rejection_is_not_solved_again(self):
        self.run_amp(image(268.0))
        calls = self.run_amp(image(400.0))

        self.assertEqual(calls["gaincal"], [])
        self.assertEqual(self.manager.current_ms, str(self.phase_ms))
        self.assertFalse(self.manager.amp_history["accepted"])

    def test_an_accepted_solution_is_restored_on_resume(self):
        self.run_amp(image(300.0))
        self.manager.current_ms = str(self.phase_ms)
        del self.manager.amp_table
        calls = self.run_amp(image(100.0))

        self.assertEqual(calls["gaincal"], [])
        self.assertEqual(self.manager.current_ms, str(self.final))
        self.assertTrue(self.manager.amp_table.endswith("amp_phase.cal"))
        self.assertTrue(self.manager.amp_history["accepted"])

    def test_disabled_it_does_nothing(self):
        self.manager.task_config["enable_amp_selfcal"] = False
        calls = self.run_amp(image(300.0))
        self.assertEqual(calls["gaincal"], [])
        self.manager.ensure_model_column.assert_not_called()

    def test_without_an_accepted_phase_loop_it_does_nothing(self):
        self.manager.phase_tables = []
        calls = self.run_amp(image(300.0))
        self.assertEqual(calls["gaincal"], [])

    def test_residual_flagging_then_subtracts_the_amp_phase_model(self):
        self.run_amp(image(300.0))
        self.manager.phase_history = [{"loop": 1, "accepted": True}]
        for name in ("amp_phase.model.tt0", "phase_1.model.tt0"):
            (self.imagedir / name).mkdir()
        self.assertEqual(
            self.manager.latest_image_base(), self.imagedir / "amp_phase"
        )


class ConfigurationTests(unittest.TestCase):

    def load(self, name, **yaml):
        spec = importlib.util.spec_from_file_location(
            name, SOURCE_DIR / "config.py"
        )
        config = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(config)
        config.CONFIG_YAML.clear()
        config.CONFIG_YAML.update(yaml)
        if yaml:
            config._apply_config_yaml_overrides()
        return config

    def test_it_stays_off_unless_asked_for(self):
        config = self.load("test_radio_continuum_config_amp_default")
        self.assertFalse(config.SELFCAL_TASK_CONFIG["enable_amp_selfcal"])

    def test_config_yaml_switches_it_on(self):
        config = self.load(
            "test_radio_continuum_config_amp_on",
            selfcal_enable_amp=True,
            selfcal_amp_max_peak_change=0.05,
        )
        self.assertTrue(config.SELFCAL_TASK_CONFIG["enable_amp_selfcal"])
        self.assertEqual(
            config.SELFCAL_TASK_CONFIG["amp_max_peak_change"], 0.05
        )

    def test_the_solve_settings_default_in_the_code_not_the_recipe(self):
        config = self.load("test_radio_continuum_config_amp_recipe")
        for key in ("amp_minsnr", "amp_solve_gaintype"):
            self.assertNotIn(key, config.SELFCAL_TASK_CONFIG)
        self.assertEqual(SELFCAL.DEFAULT_AMP_SELFCAL_MINSNR, 10.0)
        self.assertEqual(SELFCAL.DEFAULT_AMP_SELFCAL_GAINTYPE, "T")

    def test_config_yaml_sets_the_solve(self):
        config = self.load(
            "test_radio_continuum_config_amp_solve",
            selfcal_amp_minsnr=20,
            selfcal_amp_gaintype="g",
        )
        self.assertEqual(config.SELFCAL_TASK_CONFIG["amp_minsnr"], 20.0)
        self.assertEqual(
            config.SELFCAL_TASK_CONFIG["amp_solve_gaintype"], "G"
        )


    def test_selfcal_images_are_point_source_cleans_by_default(self):
        config = self.load("test_radio_continuum_config_scales_default")
        for band, settings in config.SELFCAL_IMAGING_CONFIG.items():
            with self.subTest(band=band):
                self.assertNotIn("scales", settings)

    def test_config_yaml_gives_selfcal_a_multi_scale_model(self):
        config = self.load(
            "test_radio_continuum_config_scales", selfcal_scales=[0, 5, 15]
        )
        for band, settings in config.SELFCAL_IMAGING_CONFIG.items():
            with self.subTest(band=band):
                self.assertEqual(settings["scales"], [0, 5, 15])


class SelfcalScalesTests(unittest.TestCase):

    def controls(self, image_config):
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.task_config = {}
        manager.image_config = image_config
        return manager.clean_control_kwargs()

    def test_configured_scales_reach_every_selfcal_clean(self):
        controls = self.controls({"scales": [0, 5, 15]})
        self.assertEqual(controls["scales"], [0, 5, 15])
        # The minor-cycle controls are unchanged.
        self.assertEqual(controls["gain"], 0.1)
        self.assertEqual(controls["cyclefactor"], 1.5)

    def test_without_them_tclean_is_called_as_before(self):
        self.assertNotIn("scales", self.controls({}))
        self.assertNotIn("scales", self.controls({"scales": []}))


if __name__ == "__main__":
    unittest.main()
