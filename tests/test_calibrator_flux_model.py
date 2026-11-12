"""A calibrator given its flux density outright becomes an absolute reference.

The case is uGMRT Band 5 of 2021-02-22: the transfer from 3C48 made the
phase calibrator 0323+055 2.30 Jy where RACS-mid, five months later, had
3.08 Jy at 1250 MHz (finding 154). Named as a secondary flux calibrator
with that flux, it sets the targets' scale itself.
"""

import importlib.util
import json
import shutil
import unittest
from pathlib import Path
from unittest import mock

import test_flux_model_safety as SAFETY

CALIBRATION = SAFETY.CALIBRATION
setjy_result = SAFETY.setjy_result


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_calibrator_flux_model_config", SOURCE_DIR / "config.py")
CONFIG = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CONFIG)

MODEL = {"0323+055": {"flux_jy": 3.08, "reffreq_mhz": 1250.0, "spix": [-0.76]}}


class ParseTests(unittest.TestCase):

    def test_the_flat_form_as_config_yaml_gives_it(self):
        self.assertEqual(
            CONFIG.parse_calibrator_flux_models("0323+055=3.08@1250,-0.76"),
            MODEL)
        both = CONFIG.parse_calibrator_flux_models(
            "0323+055=3.08@1250,-0.76; J0240-2309 = 5.1 @ 1400")
        self.assertEqual(both["J0240-2309"],
                         {"flux_jy": 5.1, "reffreq_mhz": 1400.0, "spix": []})

    def test_a_mapping_as_pyyaml_gives_it(self):
        self.assertEqual(CONFIG.parse_calibrator_flux_models({
            "0323+055": {"flux_jy": 3.08, "reffreq_mhz": 1250, "spix": -0.76},
        }), MODEL)

    def test_what_cannot_be_a_flux_model_is_refused(self):
        for bad in ("0323+055=3.08", "0323+055=-1@1250", "=3@1250",
                    "0323+055=x@1250", "0323+055=3@1250,nan"):
            with self.subTest(value=bad), self.assertRaises(ValueError):
                CONFIG.parse_calibrator_flux_models(bad)
        self.assertEqual(CONFIG.parse_calibrator_flux_models(""), {})


class SetjyTests(unittest.TestCase):

    # The stand-in manager of the flux-model tests, without their tests.
    setUp = SAFETY.FluxModelSafetyTests.setUp
    tearDown = SAFETY.FluxModelSafetyTests.tearDown
    manager = SAFETY.FluxModelSafetyTests.manager

    def two_references(self):
        manager = self.manager(bandpass="3C48")
        manager.task_config["calibrator_flux_models"] = dict(MODEL)
        manager.get_flux_calibrators = lambda: ["3C48", "0323+055"]
        manager.get_flux_calibrator = lambda: "3C48"
        return manager

    def test_the_stated_calibrator_gets_a_manual_model_and_3c48_its_standard(self):
        manager = self.two_references()
        calls = []

        def fake_setjy(**kwargs):
            calls.append(kwargs)
            return setjy_result(3.08 if kwargs["field"] == "0323+055" else 17.8,
                                2.9 if kwargs["field"] == "0323+055" else 16.0)

        CALIBRATION.setjy = fake_setjy
        manager.run_setjy()

        by_field = {call["field"]: call for call in calls}
        self.assertEqual(by_field["3C48"]["standard"], "Perley-Butler 2017")
        stated = by_field["0323+055"]
        self.assertEqual(stated["standard"], "manual")
        self.assertEqual(stated["fluxdensity"], [3.08, 0.0, 0.0, 0.0])
        self.assertEqual(stated["spix"], [-0.76])
        self.assertEqual(stated["reffreq"], "1250MHz")
        self.assertTrue(stated["scalebychan"])
        self.assertEqual(len(calls), 2)
        provenance = json.loads(
            (self.caldir / "flux_model_provenance.json").read_text())
        methods = {app["field"]: app["method"] for app in provenance["applications"]}
        self.assertEqual(methods, {"3C48": "setjy", "0323+055": "setjy-manual"})
        self.assertEqual(provenance["identity"]["calibrator_flux_models"], MODEL)
        self.assertTrue((self.caldir / "setjy.done").is_file())

    def test_a_changed_stated_flux_changes_the_identity(self):
        manager = self.two_references()
        before = manager.flux_model_identity("3C48", "3C48")
        manager.task_config["calibrator_flux_models"] = {
            "0323+055": dict(MODEL["0323+055"], flux_jy=3.0)}
        self.assertNotEqual(before, manager.flux_model_identity("3C48", "3C48"))
        manager.task_config.pop("calibrator_flux_models")
        self.assertNotIn("calibrator_flux_models",
                         manager.flux_model_identity("3C48", "3C48"))

    def test_a_manual_model_setjy_reports_nothing_for_is_still_kept(self):
        manager = self.two_references()

        def fake_setjy(**kwargs):
            if kwargs["field"] == "0323+055":
                return {}
            return setjy_result(17.8, 16.0)

        CALIBRATION.setjy = fake_setjy
        manager.run_setjy()
        provenance = json.loads(
            (self.caldir / "flux_model_provenance.json").read_text())
        stated = [app for app in provenance["applications"]
                  if app["field"] == "0323+055"][0]
        self.assertEqual(stated["flux_qa"]["status"], "not-returned")
        self.assertEqual(provenance["model_data_qa"]["0323+055"]["status"], "passed")


class FluxscaleTests(unittest.TestCase):

    def test_with_the_phase_calibrator_a_reference_nothing_is_bootstrapped(self):
        import tempfile

        root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, root)
        manager = CALIBRATION.CalibrationManager.__new__(
            CALIBRATION.CalibrationManager)
        manager.caldir = root
        (root / "amp.cal").mkdir()
        (root / "amp.cal" / "table.dat").write_text("amplitudes")
        manager.get_flux_calibrators = lambda: ["3C48", "0323+055"]
        manager.get_phase_calibrators = lambda: ["0323+055"]
        manager.caltable_checkpoint_ready = lambda *args: False
        manager.require_caltable = lambda name: None
        marked = []
        manager.mark_checkpoint = marked.append
        manager.archive_uncheckpointed_caltable = lambda name: None
        called = []
        with mock.patch.object(CALIBRATION, "fluxscale",
                               lambda **kwargs: called.append(kwargs)), \
                mock.patch("builtins.print"):
            manager.run_fluxscale()

        self.assertEqual(called, [])
        self.assertEqual((root / "flux.cal" / "table.dat").read_text(), "amplitudes")
        report = json.loads((root / "fluxscale_report.json").read_text())
        self.assertEqual(report["status"], "not_required")
        self.assertEqual(report["reference"], ["3C48", "0323+055"])
        self.assertEqual(marked, ["fluxscale"])


if __name__ == "__main__":
    unittest.main()
