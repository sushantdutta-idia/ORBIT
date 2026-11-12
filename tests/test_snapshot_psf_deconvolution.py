"""A short observation's point-spread function must not stop the clean.

The run that prompted these tests was a uGMRT Band 4 target observed for 33
minutes in two usable scans. Its PSF has a largest sidelobe of 0.44, so
auto-multithresh's sidelobe criterion (2.5 x 0.44 x peak) sat above the
residual peak: CASA reported "zero mask", and the published image was only
ever cleaned to the 5 sigma of the first, unmasked pass.
"""

import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_selfcal_resume_safety import load_selfcal_module
from test_imaging_atomic_safety import load_imaging_module_without_casa


SOURCE_DIR = Path(__file__).resolve().parents[1]
SELFCAL = load_selfcal_module()
IMAGING = load_imaging_module_without_casa()

_spec = importlib.util.spec_from_file_location(
    "snapshot_psf_casa_compat", SOURCE_DIR / "casa_compat.py"
)
COMPAT = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(COMPAT)

SNAPSHOT_SIDELOBE = 0.44

AFTER_BOOTSTRAP = {
    "model_nonzero": True,
    "residual_rms_Jy_beam": 9.0e-5,
    "residual_peak_abs_Jy_beam": 6.3e-4,
    "residual_peak_snr": 7.0,
}
AFTER_MASKED = {
    "model_nonzero": True,
    "residual_rms_Jy_beam": 4.0e-5,
    "residual_peak_abs_Jy_beam": 1.3e-4,
    "residual_peak_snr": 3.2,
}
AFTER_FLOOR = {
    "model_nonzero": True,
    "residual_rms_Jy_beam": 7.0e-5,
    "residual_peak_abs_Jy_beam": 3.4e-4,
    "residual_peak_snr": 4.9,
}


class SidelobeThresholdTests(unittest.TestCase):

    def test_a_snapshot_psf_lowers_the_threshold_below_the_peak(self):
        value, lowered = COMPAT.adaptive_sidelobethreshold(
            2.5, SNAPSHOT_SIDELOBE
        )
        self.assertTrue(lowered)
        # The criterion must end up below the residual peak, or no mask can
        # ever open.
        self.assertLess(value * SNAPSHOT_SIDELOBE, 1.0)
        self.assertAlmostEqual(value * SNAPSHOT_SIDELOBE, 0.6)

    def test_a_good_psf_keeps_the_configured_threshold(self):
        self.assertEqual(
            COMPAT.adaptive_sidelobethreshold(2.5, 0.10), (2.5, False)
        )

    def test_the_threshold_is_never_raised(self):
        value, lowered = COMPAT.adaptive_sidelobethreshold(1.25, 0.10)
        self.assertEqual((value, lowered), (1.25, False))

    def test_it_does_not_go_below_one(self):
        value, lowered = COMPAT.adaptive_sidelobethreshold(2.5, 0.9)
        self.assertEqual((value, lowered), (1.0, True))

    def test_an_unmeasured_psf_changes_nothing(self):
        for sidelobe in (None, 0.0, float("nan"), "unknown"):
            self.assertEqual(
                COMPAT.adaptive_sidelobethreshold(2.5, sidelobe),
                (2.5, False),
            )

    def test_zero_mask_is_recognised_from_code_or_description(self):
        self.assertTrue(COMPAT.is_zero_mask_stop({"stopcode": 7}))
        self.assertTrue(COMPAT.is_zero_mask_stop(
            {"stopcode": None, "casa_stop_description": "zero mask"}
        ))
        self.assertFalse(COMPAT.is_zero_mask_stop({"stopcode": 8}))
        self.assertFalse(COMPAT.is_zero_mask_stop(None))


class SnapshotCleanTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.imagename = str(Path(self.temporary.name) / "final")

    def tearDown(self):
        self.temporary.cleanup()

    def kwargs(self, pblimit=0.05, savemodel="modelcolumn"):
        return {
            "imagename": self.imagename,
            "vis": "target.ms",
            "niter": 100000,
            "nsigma": 3.0,
            "usemask": "auto-multithresh",
            "sidelobethreshold": 2.5,
            "noisethreshold": 4.25,
            "pblimit": pblimit,
            "savemodel": savemodel,
            "restart": False,
            "calcpsf": True,
            "calcres": True,
        }

    def run_selfcal(self, returns, sequence, kwargs=None, **task_config):
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.task_config = dict(task_config)
        manager.clean_diagnostics = {}
        diagnostics = list(sequence)
        manager.deconvolution_diagnostics = lambda base: dict(
            diagnostics.pop(0)
        )
        manager.automask_kwargs = lambda relaxed=False: {
            "usemask": "auto-multithresh",
            "sidelobethreshold": 1.5 if relaxed else 2.5,
        }
        results = list(returns)
        calls = []

        def fake_tclean(**arguments):
            calls.append(arguments)
            return results.pop(0)

        with mock.patch.object(SELFCAL, "tclean", fake_tclean), \
                mock.patch("builtins.print"):
            result = manager.run_tclean_with_recovery(
                kwargs or self.kwargs(), "phase_1"
            )
        return calls, result, manager

    def test_the_masked_pass_uses_a_threshold_the_psf_allows(self):
        calls, _, manager = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
            ],
            sequence=[AFTER_BOOTSTRAP, AFTER_MASKED],
        )

        self.assertEqual(len(calls), 2)
        masked = calls[1]
        self.assertEqual(masked["usemask"], "auto-multithresh")
        self.assertLess(masked["sidelobethreshold"] * SNAPSHOT_SIDELOBE, 1.0)
        # The depth asked of the masked clean is unchanged.
        self.assertEqual(masked["nsigma"], 3.0)
        self.assertEqual(masked["noisethreshold"], 4.25)
        self.assertEqual(
            manager.automask_adaptations["phase_1"][
                "configured_sidelobethreshold"
            ],
            2.5,
        )

    def test_a_full_synthesis_psf_is_cleaned_as_configured(self):
        calls, _, manager = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": 0.08},
                {"stopcode": 8, "maxpsfsidelobe": 0.08},
            ],
            sequence=[AFTER_BOOTSTRAP, AFTER_MASKED],
        )

        self.assertEqual(calls[1]["sidelobethreshold"], 2.5)
        self.assertFalse(hasattr(manager, "automask_adaptations"))

    def test_an_empty_mask_is_followed_by_a_threshold_clean(self):
        calls, result, manager = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 7, "iterdone": 0},
                {"stopcode": 8},
            ],
            sequence=[AFTER_BOOTSTRAP, AFTER_BOOTSTRAP, AFTER_FLOOR],
        )

        self.assertEqual(
            [item["name"] for item in manager.clean_passes["phase_1"]],
            ["bootstrap", "auto-multithresh", "unmasked floor"],
        )
        floor = calls[2]
        self.assertEqual(floor["usemask"], "pb")
        self.assertEqual(floor["nsigma"], 5.0)
        self.assertNotIn("sidelobethreshold", floor)
        # It continues from what was already cleaned ...
        self.assertTrue(floor["restart"])
        self.assertFalse(floor["calcpsf"])
        self.assertFalse(floor["calcres"])
        # ... and the gain solve still gets its model column.
        self.assertEqual(floor["savemodel"], "modelcolumn")
        self.assertEqual(result["floor_clean"]["nsigma"], 5.0)
        self.assertEqual(result["residual_peak_snr"], 4.9)
        self.assertTrue(result["final_stop"]["converged"])

    def test_a_residual_already_below_the_floor_is_left_alone(self):
        quiet = dict(AFTER_BOOTSTRAP, residual_peak_snr=4.6)
        calls, result, _ = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 7, "iterdone": 0},
            ],
            sequence=[quiet, quiet],
        )

        self.assertEqual(len(calls), 2)
        self.assertNotIn("floor_clean", result)

    def test_a_masked_clean_that_worked_needs_no_floor(self):
        calls, result, _ = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 8},
            ],
            sequence=[AFTER_BOOTSTRAP, AFTER_MASKED],
        )

        self.assertEqual(len(calls), 2)
        self.assertNotIn("floor_clean", result)

    def test_a_negative_pblimit_cleans_the_whole_image(self):
        calls, _, _ = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 7, "iterdone": 0},
                {"stopcode": 8},
            ],
            sequence=[AFTER_BOOTSTRAP, AFTER_BOOTSTRAP, AFTER_FLOOR],
            kwargs=self.kwargs(pblimit=-0.01),
        )

        for unmasked in (calls[0], calls[2]):
            # CASA's own beam model is narrower than the measured uGMRT
            # beam, so a beam mask would leave real outer-beam sources
            # undeconvolved.
            self.assertEqual(unmasked["usemask"], "user")
            self.assertEqual(unmasked["mask"], "")
            self.assertNotIn("pbmask", unmasked)
            self.assertEqual(unmasked["pblimit"], -0.01)

    def test_an_explicit_bootstrap_mask_still_wins(self):
        calls, _, _ = self.run_selfcal(
            returns=[{"stopcode": 8}, {"stopcode": 8}],
            sequence=[AFTER_BOOTSTRAP, AFTER_MASKED],
            kwargs=self.kwargs(pblimit=-0.01),
            clean_bootstrap_pbmask=0.1,
        )

        self.assertEqual(calls[0]["usemask"], "pb")
        self.assertEqual(calls[0]["pbmask"], 0.1)

    def test_final_imaging_follows_the_same_path(self):
        manager = IMAGING.FinalImagingManager.__new__(
            IMAGING.FinalImagingManager
        )
        manager.config = {}
        manager.clean_diagnostics = {}
        diagnostics = [AFTER_BOOTSTRAP, AFTER_BOOTSTRAP, AFTER_FLOOR]
        manager.deconvolution_diagnostics = lambda base: dict(
            diagnostics.pop(0)
        )
        results = [
            {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
            {"stopcode": 7, "iterdone": 0},
            {"stopcode": 8},
        ]
        calls = []

        def fake_tclean(**arguments):
            calls.append(arguments)
            return results.pop(0)

        with mock.patch.object(IMAGING, "tclean", fake_tclean), \
                mock.patch("builtins.print"):
            result = manager.run_tclean_with_recovery(
                self.kwargs(pblimit=-0.01, savemodel="none"),
                "final_continuum",
            )

        self.assertEqual(
            [item["name"] for item in manager.clean_passes["final_continuum"]],
            ["bootstrap", "auto-multithresh", "unmasked floor"],
        )
        self.assertLess(
            calls[1]["sidelobethreshold"] * SNAPSHOT_SIDELOBE, 1.0
        )
        self.assertEqual(calls[2]["usemask"], "user")
        self.assertEqual(calls[2]["savemodel"], "none")
        self.assertIn("floor_clean", result)


STALLED = {
    "model_nonzero": True,
    "residual_rms_Jy_beam": 7.5e-5,
    # The artefact beside the brightest source: ten sigma, outside the mask,
    # and not something more cleaning removes.
    "residual_peak_abs_Jy_beam": 7.9e-4,
    "residual_peak_snr": 10.5,
}


class MaskedMajorCycleLimitTests(unittest.TestCase):
    """A masked pass that oscillates around its threshold is bounded."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.imagename = str(Path(self.temporary.name) / "phase_3")

    def tearDown(self):
        self.temporary.cleanup()

    def run_selfcal(self, returns, sequence, **task_config):
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.task_config = dict(task_config)
        manager.clean_diagnostics = {}
        diagnostics = list(sequence)
        manager.deconvolution_diagnostics = lambda base: dict(
            diagnostics.pop(0)
        )
        results = list(returns)
        calls = []

        def fake_tclean(**arguments):
            calls.append(arguments)
            return results.pop(0)

        kwargs = {
            "imagename": self.imagename,
            "vis": "phase2.ms",
            "niter": 100000,
            "nsigma": 3.0,
            "usemask": "auto-multithresh",
            "sidelobethreshold": 2.5,
            "noisethreshold": 4.25,
            "pblimit": -0.01,
            "savemodel": "modelcolumn",
            "restart": False,
            "calcpsf": True,
            "calcres": True,
        }
        with mock.patch.object(SELFCAL, "tclean", fake_tclean), \
                mock.patch("builtins.print"):
            manager.run_tclean_with_recovery(kwargs, "phase_3")
        return calls, manager

    def test_only_the_masked_pass_is_bounded(self):
        calls, _ = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 8},
            ],
            sequence=[AFTER_BOOTSTRAP, AFTER_MASKED],
            clean_max_major_cycles=20,
        )

        bootstrap, masked = calls
        self.assertNotIn("nmajor", bootstrap)
        self.assertEqual(masked["usemask"], "auto-multithresh")
        self.assertEqual(masked["nmajor"], 20)
        # The depth asked of it is unchanged.
        self.assertEqual(masked["nsigma"], 3.0)
        self.assertEqual(masked["niter"], 100000)

    def test_the_built_in_default_applies_when_nothing_is_configured(self):
        calls, _ = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 8},
            ],
            sequence=[AFTER_BOOTSTRAP, AFTER_MASKED],
        )
        self.assertEqual(
            calls[1]["nmajor"], SELFCAL.DEFAULT_CLEAN_MAX_MAJOR_CYCLES
        )

    def test_zero_leaves_casa_unbounded(self):
        calls, _ = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 8},
            ],
            sequence=[AFTER_BOOTSTRAP, AFTER_MASKED],
            clean_max_major_cycles=0,
        )
        self.assertNotIn("nmajor", calls[1])

    def test_a_pass_stopped_by_the_bound_is_not_continued(self):
        calls, manager = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 9, "nmajordone": 20},
            ],
            sequence=[AFTER_BOOTSTRAP, STALLED],
            clean_max_major_cycles=20,
        )

        self.assertEqual(len(calls), 2)
        self.assertEqual(
            [item["name"] for item in manager.clean_passes["phase_3"]],
            ["bootstrap", "auto-multithresh"],
        )
        self.assertEqual(
            manager.clean_passes["phase_3"][1]["stop_reason"],
            "reached the major-cycle limit (nmajor)",
        )

    def test_a_pass_out_of_iterations_is_still_continued_and_bounded(self):
        calls, _ = self.run_selfcal(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 1},
                {"stopcode": 8},
            ],
            sequence=[AFTER_BOOTSTRAP, STALLED, AFTER_MASKED],
            clean_max_major_cycles=20,
            clean_max_continuations=3,
        )

        self.assertEqual(len(calls), 3)
        continuation = calls[2]
        self.assertEqual(continuation["nmajor"], 20)
        self.assertTrue(continuation["restart"])
        self.assertFalse(continuation["calcres"])


class FinalImageMajorCycleLimitTests(unittest.TestCase):
    """The final image's masked pass stalled the same way, and is bounded."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.imagename = str(Path(self.temporary.name) / "final_continuum")

    def tearDown(self):
        self.temporary.cleanup()

    def run_imaging(self, returns, sequence, **task_config):
        manager = IMAGING.FinalImagingManager.__new__(
            IMAGING.FinalImagingManager
        )
        manager.config = {}
        manager.clean_diagnostics = {}
        diagnostics = list(sequence)
        manager.deconvolution_diagnostics = lambda base: dict(
            diagnostics.pop(0)
        )
        results = list(returns)
        calls = []

        def fake_tclean(**arguments):
            calls.append(arguments)
            return results.pop(0)

        kwargs = {
            "imagename": self.imagename,
            "vis": "phase3.ms",
            "niter": 200000,
            "nsigma": 3.0,
            "usemask": "auto-multithresh",
            "sidelobethreshold": 2.5,
            "noisethreshold": 4.25,
            "pblimit": -0.01,
            "savemodel": "none",
            "restart": False,
            "calcpsf": True,
            "calcres": True,
        }
        with mock.patch.object(IMAGING, "tclean", fake_tclean), \
                mock.patch.dict(IMAGING.IMAGING_TASK_CONFIG, task_config), \
                mock.patch("builtins.print"):
            manager.run_tclean_with_recovery(kwargs, "final_continuum")
        return calls, manager

    def test_the_masked_pass_is_bounded_by_default(self):
        calls, _ = self.run_imaging(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 8},
            ],
            sequence=[AFTER_BOOTSTRAP, AFTER_MASKED],
        )

        bootstrap, masked = calls
        self.assertNotIn("nmajor", bootstrap)
        self.assertEqual(
            masked["nmajor"], IMAGING.DEFAULT_CLEAN_MAX_MAJOR_CYCLES
        )
        self.assertEqual(masked["nsigma"], 3.0)
        self.assertEqual(masked["niter"], 200000)

    def test_zero_leaves_casa_unbounded(self):
        calls, _ = self.run_imaging(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 8},
            ],
            sequence=[AFTER_BOOTSTRAP, AFTER_MASKED],
            clean_max_major_cycles=0,
        )
        self.assertNotIn("nmajor", calls[1])

    def test_a_pass_stopped_by_the_bound_is_not_continued(self):
        calls, manager = self.run_imaging(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 9, "nmajordone": 20},
            ],
            sequence=[AFTER_BOOTSTRAP, STALLED],
        )

        self.assertEqual(len(calls), 2)
        final = manager.clean_passes["final_continuum"][-1]
        self.assertEqual(final["name"], "auto-multithresh")
        self.assertFalse(final["converged"])

    def test_a_pass_out_of_iterations_is_still_continued(self):
        calls, _ = self.run_imaging(
            returns=[
                {"stopcode": 8, "maxpsfsidelobe": SNAPSHOT_SIDELOBE},
                {"stopcode": 1},
                {"stopcode": 8},
            ],
            sequence=[AFTER_BOOTSTRAP, STALLED, AFTER_MASKED],
            clean_max_continuations=3,
        )

        self.assertEqual(len(calls), 3)
        self.assertEqual(
            calls[2]["nmajor"], IMAGING.DEFAULT_CLEAN_MAX_MAJOR_CYCLES
        )


class MajorCycleLimitConfigurationTests(unittest.TestCase):

    def load_config(self, name):
        spec = importlib.util.spec_from_file_location(
            name, SOURCE_DIR / "config.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.CONFIG_YAML.clear()
        return module

    def test_the_default_sits_just_below_the_longest_observed_pass(self):
        # The run that set it needed 13, 19 and 21 major cycles to converge.
        self.assertEqual(SELFCAL.DEFAULT_CLEAN_MAX_MAJOR_CYCLES, 20)

    def test_the_default_leaves_the_selfcal_recipe_alone(self):
        # A self-calibration under way keeps its accepted loops only while
        # its recipe, which includes SELFCAL_TASK_CONFIG, is unchanged.
        config = self.load_config("test_radio_continuum_config_nmajor")
        self.assertNotIn(
            "clean_max_major_cycles", config.SELFCAL_TASK_CONFIG
        )

    def test_the_imaging_default_is_the_same_and_also_off_the_recipe(self):
        self.assertEqual(IMAGING.DEFAULT_CLEAN_MAX_MAJOR_CYCLES, 20)
        config = self.load_config("test_radio_continuum_config_nmajor_img")
        self.assertNotIn(
            "clean_max_major_cycles", config.FINAL_IMAGING_TASK_CONFIG
        )

    def test_config_yaml_sets_it(self):
        config = self.load_config("test_radio_continuum_config_nmajor_yaml")
        config.CONFIG_YAML["selfcal_clean_max_major_cycles"] = 12
        config._apply_config_yaml_overrides()
        self.assertEqual(
            config.SELFCAL_TASK_CONFIG["clean_max_major_cycles"], 12
        )
        self.assertNotIn(
            "clean_max_major_cycles", config.FINAL_IMAGING_TASK_CONFIG
        )

    def test_one_config_yaml_key_sets_both_stages(self):
        config = self.load_config("test_radio_continuum_config_nmajor_both")
        config.CONFIG_YAML["clean_max_major_cycles"] = 15
        config._apply_config_yaml_overrides()
        self.assertEqual(
            config.SELFCAL_TASK_CONFIG["clean_max_major_cycles"], 15
        )
        self.assertEqual(
            config.FINAL_IMAGING_TASK_CONFIG["clean_max_major_cycles"], 15
        )

    def test_the_imaging_key_sets_only_imaging(self):
        config = self.load_config("test_radio_continuum_config_nmajor_i")
        config.CONFIG_YAML["imaging_clean_max_major_cycles"] = 30
        config._apply_config_yaml_overrides()
        self.assertEqual(
            config.FINAL_IMAGING_TASK_CONFIG["clean_max_major_cycles"], 30
        )
        self.assertNotIn("clean_max_major_cycles", config.SELFCAL_TASK_CONFIG)

    def test_config_yaml_rejects_a_negative_limit(self):
        config = self.load_config("test_radio_continuum_config_nmajor_bad")
        config.CONFIG_YAML["selfcal_clean_max_major_cycles"] = -1
        with self.assertRaisesRegex(ValueError, "0 \\(no limit\\)"):
            config._apply_config_yaml_overrides()


if __name__ == "__main__":
    unittest.main()
