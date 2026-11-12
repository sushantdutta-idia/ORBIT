"""A field too faint for amplitudes is self-calibrated in phase only.

The case is Band 5 (2026-10-24): in G0231-0533 and G0216-0604's target
stage the amplitude rounds lowered the compact sources' peaks by 5 to 14%.
"""

import importlib.util
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load(name, module_name):
    spec = importlib.util.spec_from_file_location(module_name, SOURCE_DIR / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RECIPE = load("selfcal_recipe.py", "test_phase_only_recipe")
TARGETS = load("run_targets.py", "test_phase_only_targets")


def plan(*extra):
    captured = []
    with mock.patch.object(
        RECIPE, "run", lambda args, plan=None: captured.append((args, plan))
    ):
        RECIPE.main([
            "--target", "G1", "--ms", "data.ms", "--imsize", "6144",
            "--cell", "0.5", "--channel-range", "31", "993",
            "--casa", "python3", *extra,
        ])
    return captured[0]


class RecipeTests(unittest.TestCase):

    def test_no_round_solves_an_amplitude(self):
        args, rounds = plan("--phase-only", "--time-phase-solint", "1min")

        self.assertEqual([stage["id"] for stage in rounds],
                         ["phase", "phase_time", "delay", "phase_field"])
        for stage in rounds:
            with self.subTest(stage=stage["id"]):
                self.assertEqual(stage["solve"]["calmode"], "p")
                self.assertFalse(stage["solve"].get("bandpass"))
                command = RECIPE.solve_command(args, "a.ms", stage["solve"])
                self.assertEqual(
                    command[command.index("--calmode") + 1], "p")
                self.assertNotIn("--detrend-radius", command)
        field = rounds[-1]
        self.assertTrue(field["model"]["full"])
        self.assertEqual(field["until"], {"phase_deg": 5.0})
        self.assertEqual(field["most"], 3)

    def test_the_cores_amplitude_rounds_go_too(self):
        _, rounds = plan("--phase-only", "--amplitude-cores")
        self.assertNotIn("amplitude_cores", [stage["id"] for stage in rounds])

    def test_after_peeling_the_whole_field_is_solved_in_phase(self):
        _, rounds = plan("--phase-only", "--whole-field-only")
        self.assertEqual([stage["id"] for stage in rounds], ["phase_field"])

    def test_without_the_option_the_plan_is_as_it_was(self):
        _, rounds = plan()
        self.assertEqual([stage["id"] for stage in rounds],
                         ["phase", "delay", "amplitude", "bandpass"])


class TargetsTests(unittest.TestCase):

    def args(self, *extra):
        parsed = []
        with mock.patch.object(TARGETS, "run", parsed.append):
            TARGETS.main(["--channel-range", "122", "3974",
                          "--average-channels", "4", *extra])
        return parsed[0]

    def test_both_stages_are_told(self):
        args = self.args("--phase-only")
        field = TARGETS.field_command(args, "G1", "kept.ms", Path("w/G1"))
        target = TARGETS.target_command(args, "G1", "peeled.ms", Path("w/G1"))
        self.assertIn("--phase-only", field)
        self.assertNotIn("--amplitude-cores", field)
        self.assertIn("--phase-only", target)

    def test_after_peeling_only(self):
        args = self.args("--phase-only-after-peel")
        field = TARGETS.field_command(args, "G1", "kept.ms", Path("w/G1"))
        target = TARGETS.target_command(args, "G1", "peeled.ms", Path("w/G1"))
        self.assertIn("--amplitude-cores", field)
        self.assertNotIn("--phase-only", field)
        self.assertIn("--phase-only", target)

    def test_without_it_amplitudes_are_solved_as_before(self):
        args = self.args()
        field = TARGETS.field_command(args, "G1", "kept.ms", Path("w/G1"))
        target = TARGETS.target_command(args, "G1", "peeled.ms", Path("w/G1"))
        self.assertIn("--amplitude-cores", field)
        self.assertNotIn("--phase-only", field + target)


if __name__ == "__main__":
    unittest.main()
