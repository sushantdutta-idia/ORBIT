"""A science change made in config.yaml is not adopted as a code-only upgrade.

The case is Band 5 (2026-10-25): calibration_calibrator_flux_models and
secondary_flux_calibrator, set in config.yaml, changed the calibration
recipe and left its command line as it was, and the completed calibration
was adopted as though only the code had changed.
"""

import importlib.util
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_pipeline():
    spec = importlib.util.spec_from_file_location(
        "test_recipe_change_pipeline", SOURCE_DIR / "pipeline.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PIPELINE = load_pipeline()


class CompatibilityTests(unittest.TestCase):

    def runner(self, current_recipe="new-recipe"):
        runner = PIPELINE.PipelineRunner.__new__(PIPELINE.PipelineRunner)
        runner.recorded_recipe_matches_after_provenance_refreshes = (
            lambda stage, recorded: False)
        runner.recorded_recipe_matches_apart_from_failure_policy = (
            lambda stage, recorded, previous: False)
        runner.normalized_checkpoint_resume_command = lambda command: command
        runner.command_for_stage = lambda stage: ["calibration.py", "a.ms"]
        self.current_recipe = current_recipe
        return runner

    def previous(self, recipe):
        return {
            "status": "completed",
            "product_stage_fingerprint": "old-implementation",
            "product_recipe_fingerprint": recipe,
            "command": ["calibration.py", "a.ms"],
        }

    def test_a_changed_calibration_recipe_is_not_adopted(self):
        compatible, reason = self.runner().completed_stage_implementation_is_compatible(
            "calibration", self.previous("old-recipe"), "new-recipe")
        self.assertFalse(compatible)
        self.assertIn("science-recipe fingerprint changed", reason)

    def test_an_unchanged_recipe_is_adopted(self):
        compatible, reason = self.runner().completed_stage_implementation_is_compatible(
            "calibration", self.previous("new-recipe"), "new-recipe")
        self.assertTrue(compatible)

    def test_a_legacy_record_whose_recipe_was_its_implementation_uses_the_command(self):
        previous = self.previous("old-implementation")
        compatible, reason = self.runner().completed_stage_implementation_is_compatible(
            "flagging", previous, "new-recipe")
        self.assertTrue(compatible)
        self.assertIn("normalized science command matches", reason)

    def test_stages_without_a_recipe_of_their_own_still_use_the_command(self):
        compatible, _ = self.runner().completed_stage_implementation_is_compatible(
            "inspection", self.previous("old-recipe"), "new-recipe")
        self.assertTrue(compatible)


if __name__ == "__main__":
    unittest.main()
