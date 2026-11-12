"""The deconvolution has to actually happen, and be seen to have happened."""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_selfcal_resume_safety import load_selfcal_module


SOURCE_DIR = Path(__file__).resolve().parents[1]
SELFCAL = load_selfcal_module()


from test_imaging_atomic_safety import load_imaging_module_without_casa


IMAGING = load_imaging_module_without_casa()

DIRTY = {
    # A dirty map of a field with several comparable sources: the rms is
    # sidelobe confusion, so auto-multithresh's threshold sits at the peak.
    "model_nonzero": False,
    "residual_rms_Jy_beam": 3.0e-4,
    "residual_peak_abs_Jy_beam": 1.5e-3,
    "residual_peak_snr": 5.0,
}
CLEANED = {
    "model_nonzero": True,
    "model_flux_Jy": 0.0045,
    "residual_rms_Jy_beam": 3.1e-5,
    "residual_peak_abs_Jy_beam": 1.2e-4,
    "residual_peak_snr": 3.9,
}


class BootstrapCleanTests(unittest.TestCase):
    """The bright sources come out before auto-multithresh measures noise."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.imagename = str(Path(self.temporary.name) / "final")
        for suffix in (".pb.tt0",):
            path = Path(self.imagename + suffix)
            path.mkdir(parents=True)
            (path / "table.dat").touch()

    def tearDown(self):
        self.temporary.cleanup()

    def kwargs(self):
        return {
            "imagename": self.imagename,
            "vis": "target.ms",
            "niter": 100000,
            "nsigma": 3.0,
            "usemask": "auto-multithresh",
            "sidelobethreshold": 2.5,
            "noisethreshold": 4.25,
            "savemodel": "modelcolumn",
            "restart": False,
            "calcpsf": True,
            "calcres": True,
        }

    def selfcal_manager(self, sequence, **task_config):
        manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        manager.task_config = dict(task_config)
        manager.clean_diagnostics = {}
        self.sequence = list(sequence)
        manager.deconvolution_diagnostics = lambda base: dict(
            self.sequence.pop(0)
        )
        manager.automask_kwargs = lambda relaxed=False: {
            "usemask": "auto-multithresh",
            "sidelobethreshold": 1.5 if relaxed else 2.5,
        }
        return manager

    def run_selfcal(self, sequence, **task_config):
        manager = self.selfcal_manager(sequence, **task_config)
        calls = []
        with mock.patch.object(
            SELFCAL, "tclean", lambda **kwargs: calls.append(kwargs)
        ):
            diagnostics = manager.run_tclean_with_recovery(
                self.kwargs(), "initial"
            )
        return calls, diagnostics

    def test_the_first_pass_cleans_inside_the_beam_not_behind_a_mask(self):
        calls, _ = self.run_selfcal([DIRTY, CLEANED])

        bootstrap = calls[0]
        self.assertEqual(bootstrap["usemask"], "pb")
        self.assertEqual(bootstrap["pbmask"], 0.2)
        # Shallow on purpose: it only has to make the rms mean something.
        # auto-multithresh then cleans every source down to nsigma inside
        # its mask, which a 5-sigma unmasked pass would have pre-empted.
        self.assertEqual(bootstrap["nsigma"], 8.0)
        self.assertEqual(bootstrap["niter"], 20000)
        # auto-multithresh's own parameters must not survive into a pass
        # that is not using it, and the bootstrap writes no model column.
        self.assertNotIn("sidelobethreshold", bootstrap)
        self.assertNotIn("noisethreshold", bootstrap)
        self.assertEqual(bootstrap["savemodel"], "none")
        # It is the pass that makes the PSF.
        self.assertTrue(bootstrap["calcpsf"])

    def test_the_masked_pass_continues_from_the_bootstrap(self):
        calls, diagnostics = self.run_selfcal([DIRTY, CLEANED])

        self.assertEqual(len(calls), 2)
        deep = calls[1]
        self.assertEqual(deep["usemask"], "auto-multithresh")
        self.assertEqual(deep["niter"], 100000)
        self.assertEqual(deep["nsigma"], 3.0)
        # Components add to the ones already found rather than starting over.
        self.assertTrue(deep["restart"])
        self.assertFalse(deep["calcpsf"])
        self.assertFalse(deep["calcres"])
        # The gain solve reads the model column, so the deep pass writes it.
        self.assertEqual(deep["savemodel"], "modelcolumn")
        self.assertEqual(diagnostics["bootstrap"]["residual_rms_Jy_beam"],
                         DIRTY["residual_rms_Jy_beam"])

    def test_no_recovery_is_needed_once_the_bootstrap_has_a_model(self):
        calls, diagnostics = self.run_selfcal([DIRTY, CLEANED])

        self.assertFalse(diagnostics["recovery_attempted"])
        self.assertNotIn("fallback_attempted", diagnostics)
        self.assertEqual(len(calls), 2)

    def test_it_can_be_switched_off_or_sized(self):
        calls, _ = self.run_selfcal([CLEANED], clean_bootstrap_enabled=False)
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["usemask"], "auto-multithresh")

        calls, _ = self.run_selfcal([CLEANED], clean_bootstrap_niter=0)
        self.assertEqual(len(calls), 1)

        calls, _ = self.run_selfcal(
            [DIRTY, CLEANED],
            clean_bootstrap_nsigma=7.0,
            clean_bootstrap_niter=1000,
            clean_bootstrap_pbmask=0.35,
        )
        self.assertEqual(calls[0]["nsigma"], 7.0)
        self.assertEqual(calls[0]["niter"], 1000)
        self.assertEqual(calls[0]["pbmask"], 0.35)

    def test_the_recovery_paths_still_guard_a_bootstrap_that_found_nothing(self):
        # A genuinely empty field: no pass produces a model, and the relaxed
        # retry is still reached.
        calls, diagnostics = self.run_selfcal(
            [DIRTY, DIRTY, DIRTY, DIRTY],
            automask_retry_enabled=True,
            automask_retry_min_peak_snr=4.0,
        )
        # Bootstrap, auto-multithresh, relaxed retry, primary-beam clean.
        self.assertEqual(len(calls), 4)
        self.assertEqual(calls[0]["usemask"], "pb")
        self.assertEqual(calls[3]["usemask"], "pb")
        self.assertTrue(diagnostics["fallback_attempted"])
        self.assertFalse(diagnostics["fallback_succeeded"])
        self.assertTrue(diagnostics["before_fallback"]["recovery_attempted"])

    def test_final_imaging_bootstraps_the_same_way(self):
        manager = IMAGING.FinalImagingManager.__new__(IMAGING.FinalImagingManager)
        manager.clean_diagnostics = {}
        sequence = [DIRTY, CLEANED]
        manager.deconvolution_diagnostics = lambda base: dict(sequence.pop(0))
        calls = []
        with mock.patch.dict(
            IMAGING.IMAGING_TASK_CONFIG, {}, clear=False
        ), mock.patch.object(
            IMAGING, "tclean", lambda **kwargs: calls.append(kwargs)
        ), mock.patch.object(
            IMAGING, "adaptive_automask_kwargs",
            lambda config, relaxed=False: {"usemask": "auto-multithresh"},
        ):
            manager.run_tclean_with_recovery(self.kwargs(), "final_continuum")

        self.assertEqual(len(calls), 2)
        self.assertEqual(calls[0]["usemask"], "pb")
        self.assertEqual(calls[1]["usemask"], "auto-multithresh")


class UndeconvolvedPublicationTests(unittest.TestCase):
    """A dirty map is not a final image."""

    def manager(self, diagnostics):
        value = IMAGING.FinalImagingManager.__new__(IMAGING.FinalImagingManager)
        value.clean_diagnostics = {"final_continuum": dict(diagnostics)}
        return value

    def verify(self, diagnostics):
        self.manager(diagnostics).verify_deconvolution("final_continuum")

    def test_real_emission_that_was_never_cleaned_stops_the_stage(self):
        with self.assertRaises(RuntimeError) as raised:
            self.verify(DIRTY)

        message = str(raised.exception)
        self.assertIn("dirty map", message)
        self.assertIn("5.0-sigma", message)
        # It says what to do about it rather than only what went wrong.
        self.assertIn("clean_bootstrap_niter", message)
        self.assertIn("undeconvolved_max_peak_snr", message)

    def test_a_field_with_nothing_in_it_is_published_and_labelled(self):
        empty = dict(DIRTY, residual_peak_snr=2.1)
        printed = []
        with mock.patch("builtins.print", lambda *a, **k: printed.append(a)):
            self.verify(empty)

        text = " ".join(str(part) for call in printed for part in call)
        self.assertIn("nothing was deconvolved", text)
        self.assertIn("point-spread function", text)

    def test_a_deconvolved_image_reports_what_it_took_out(self):
        printed = []
        with mock.patch("builtins.print", lambda *a, **k: printed.append(a)):
            self.verify(CLEANED)

        text = " ".join(str(part) for call in printed for part in call)
        self.assertIn("0.0045 Jy in the CLEAN model", text)

    def test_an_unmeasurable_residual_does_not_block_publication(self):
        printed = []
        with mock.patch("builtins.print", lambda *a, **k: printed.append(a)):
            self.verify(dict(DIRTY, residual_peak_snr=None))

        text = " ".join(str(part) for call in printed for part in call)
        self.assertIn("nothing was deconvolved", text)

    def test_the_threshold_is_configurable(self):
        with mock.patch.dict(
            IMAGING.IMAGING_TASK_CONFIG,
            {"undeconvolved_max_peak_snr": 50.0},
        ):
            printed = []
            with mock.patch("builtins.print", lambda *a, **k: printed.append(a)):
                self.verify(DIRTY)
            self.assertTrue(printed)


if __name__ == "__main__":
    unittest.main()
