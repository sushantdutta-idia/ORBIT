import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_tool():
    spec = importlib.util.spec_from_file_location(
        "test_peel_sources_module", SOURCE_DIR / "peel_sources.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TOOL = load_tool()


def option(command, name, count=1):
    index = command.index(name)
    values = command[index + 1:index + 1 + count]
    return values[0] if count == 1 else values


class StartTests(unittest.TestCase):

    def test_it_starts_beside_the_pipelines_own_inspect_module(self):
        source = (SOURCE_DIR / "peel_sources.py").read_text()
        self.assertIn("sys.path.pop(0)", source.split("import argparse")[0])
        done = subprocess.run(
            [sys.executable, "peel_sources.py", "--help"],
            cwd=SOURCE_DIR, capture_output=True, text=True,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("--source", done.stdout)


class ChooseTests(unittest.TestCase):

    def field(self):
        import numpy as np

        plane = np.zeros((5120, 5120))
        plane[2562, 2567] = 39e-3        # the target, at the centre
        plane[2554, 2464] = 11.6e-3      # 103 arcsec from it: stays
        plane[3162, 4227] = 51e-3
        plane[1641, 2103] = 33e-3
        plane[4148, 2906] = 19e-3
        plane[1820, 3216] = 11.4e-3
        plane[3000, 2895] = 6.4e-3       # too faint
        return plane

    def test_bright_sources_away_from_the_centre_are_taken_brightest_first(self):
        chosen = TOOL.choose_sources(self.field(), 1.0, 10.0, 4.0, 8, 90)

        self.assertEqual(
            [(source["x"], source["y"]) for source in chosen],
            [(3162, 4227), (1641, 2103), (4148, 2906), (1820, 3216)],
        )
        self.assertEqual(chosen[0]["offset_arcmin"], 29.5)

    def test_how_much_is_solved_follows_from_how_bright_a_source_is(self):
        chosen = TOOL.choose_sources(self.field(), 1.0, 10.0, 4.0, 8, 90)

        # An eighth of the peak for the cores; the bandpass in 30 channels
        # for 51 mJy, in 60 for 33, and none below 20.
        self.assertEqual([source["cores"] for source in chosen],
                         [6.38, 4.12, 2.38, 1.43])
        self.assertEqual([source["bins"] for source in chosen],
                         [30, 60, 0, 0])

    def test_the_number_taken_can_be_limited_and_the_threshold_moved(self):
        self.assertEqual(
            len(TOOL.choose_sources(self.field(), 1.0, 10.0, 4.0, 2, 90)), 2)
        self.assertEqual(
            len(TOOL.choose_sources(self.field(), 1.0, 5.0, 4.0, 8, 90)), 5)
        self.assertEqual(
            TOOL.choose_sources(self.field(), 1.0, 60.0, 4.0, 8, 90), [])

    def test_a_source_beside_a_brighter_one_is_never_taken_on_its_own(self):
        # Even with almost nothing kept by radius: the 11.6 mJy source is
        # within two discs of the target, and its disc would cut into it.
        chosen = TOOL.choose_sources(self.field(), 1.0, 10.0, 0.5, 8, 90)

        self.assertNotIn((2554, 2464),
                         [(source["x"], source["y"]) for source in chosen])
        self.assertEqual(len(chosen), 4)


    def test_a_brighter_source_is_taken_out_with_a_larger_disc(self):
        # 90 pixels up to 30 mJy/beam, by the square root of the peak above.
        self.assertEqual(TOOL.disc_radius(30.0, 2000, 90, 240), 90)
        self.assertEqual(TOOL.disc_radius(11.0, 2000, 90, 240), 90)
        self.assertEqual(TOOL.disc_radius(120.0, 2000, 90, 240), 180)
        self.assertEqual(TOOL.disc_radius(257.0, 2000, 90, 240), 263)
        # Never more than three times the disc given.
        self.assertEqual(TOOL.disc_radius(5000.0, 2000, 90, 240), 270)

    def test_a_disc_never_reaches_into_the_kept_centre(self):
        # 10 arcmin out with 4 kept: 360 pixels to spare, enough.
        self.assertEqual(TOOL.disc_radius(257.0, 600, 90, 240), 263)
        # 5 arcmin out: 60 to spare, and the disc given is the least.
        self.assertEqual(TOOL.disc_radius(257.0, 300, 90, 240), 90)
        self.assertEqual(TOOL.disc_radius(257.0, 450, 90, 240), 210)

    def test_a_source_inside_a_bright_ones_disc_goes_out_with_it(self):
        plane = self.field()
        plane[3162, 4227] = 257e-3       # a disc of 263 pixels
        plane[3162, 4027] = 20e-3        # 200 pixels from it: inside
        plane[3462, 4227] = 15e-3        # 300 pixels from it: outside

        chosen = TOOL.choose_sources(plane, 1.0, 10.0, 4.0, 8, 90)

        places = [(source["x"], source["y"]) for source in chosen]
        self.assertEqual(chosen[0]["radius"], 263)
        self.assertNotIn((3162, 4027), places)
        self.assertIn((3462, 4227), places)
        # The others, fainter than 30 mJy/beam or nearly, keep about 90.
        self.assertEqual([source["radius"] for source in chosen[1:]],
                         [94, 90, 90, 90])


class FlowTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.products = self.root / "uGMRT_Products"
        self.products.mkdir()
        self.ms = self.root / "data" / "field.ms"
        self.ms.mkdir(parents=True)
        self.work = self.root / "peel"
        self.calls = []
        self.commands = []
        self.fail_at = None
        test = self

        def run_command(command, log_path):
            script = Path(command[2]).name
            if script == "experiment_extract_source.py":
                kind = "peel" if "--peel" in command else "extract"
                name = option(command, "--name")
            elif script == "selfcal_recipe.py":
                kind = "recipe"
                name = Path(option(command, "--workdir")).name[len("recipe_"):]
            else:
                kind, name = "image", "field"
            test.calls.append((kind, name))
            test.commands.append(command)
            test.assertTrue(Path(log_path).parent.is_dir())
            if test.fail_at == len(test.calls):
                raise TOOL.IMAGER.CommandFailed(f"{kind} exited with 1", 1)
            folder = test.work / name
            if kind == "extract":
                test.assertTrue(Path(option(command, "--ms")).is_dir())
                (folder / f"{name}.ms").mkdir(parents=True, exist_ok=True)
            elif kind == "recipe":
                test.assertTrue(Path(option(command, "--ms")).is_dir())
                report = Path(option(command, "--workdir")) / "G1"
                report.mkdir(parents=True, exist_ok=True)
                (report / "recipe_report.json").write_text(json.dumps({
                    "before": {"peak_Jy_beam": 0.050,
                               "most_negative_Jy_beam": -0.00075},
                    "final": {"peak_Jy_beam": 0.0525,
                              "most_negative_Jy_beam": -0.0007},
                    "reading": ["phase1: largest phase 8 deg (S02)."],
                }))
            elif kind == "peel":
                test.assertTrue(
                    (Path(option(command, "--peel")) / "recipe_report.json")
                    .is_file()
                )
                (folder / f"{name}_peeled.ms").mkdir(parents=True)
            else:
                test.assertTrue(Path(option(command, "--ms")).is_dir())

        patch = mock.patch.object(TOOL.IMAGER, "run_command", run_command)
        patch.start()
        self.addCleanup(patch.stop)

    def run_tool(self, *extra):
        with mock.patch("builtins.print"):
            return TOOL.main([
                "--target", "G1", "--ms", str(self.ms),
                "--models", "images/run_I",
                "--workdir", str(self.work), "--outdir", str(self.products),
                "--channel-range", "51", "768",
                "--source", "3162", "4227", "6", "30",
                "--source", "1641", "2103", "4", "60", *extra,
            ])

    def report(self):
        return json.loads((self.work / TOOL.REPORT_NAME).read_text())

    def test_each_source_is_taken_from_what_the_last_one_left(self):
        self.run_tool()

        self.assertEqual(self.calls, [
            ("extract", "source1"), ("recipe", "source1"), ("peel", "source1"),
            ("extract", "source2"), ("recipe", "source2"), ("peel", "source2"),
            ("image", "field"),
        ])
        first, second = self.commands[0], self.commands[3]
        self.assertEqual(option(first, "--ms"), str(self.ms))
        self.assertEqual(option(first, "--pixel", 2), ["3162", "4227"])
        self.assertNotIn("--also-without", first)
        self.assertEqual(
            option(second, "--ms"),
            str(self.work / "source1" / "source1_peeled.ms"),
        )
        # The second source is taken with the first already gone, and
        # peeled so too, each with the disc the first went out with.
        self.assertEqual(option(second, "--also-without", 3),
                         ["3162", "4227", "90"])
        self.assertEqual(
            option(self.commands[5], "--also-without", 3),
            ["3162", "4227", "90"],
        )
        self.assertEqual(
            option(self.commands[6], "--ms"),
            str(self.work / "source2" / "source2_peeled.ms"),
        )

    def test_each_source_gets_its_own_threshold_and_bandpass_bins(self):
        self.run_tool()

        first, second = self.commands[1], self.commands[4]
        self.assertEqual(option(first, "--core-threshold"), "6")
        self.assertEqual(option(first, "--bandpass-channels"), "30")
        self.assertEqual(option(second, "--core-threshold"), "4")
        self.assertEqual(option(second, "--bandpass-channels"), "60")
        self.assertEqual(option(first, "--imsize"), "1024")
        self.assertIn("--image-first", first)
        self.assertEqual(
            option(first, "--ms"),
            str(self.work / "source1" / "source1.ms"),
        )
        self.assertEqual(option(first, "--channel-range", 2), ["51", "768"])

    def test_no_source_is_solved_against_a_cores_only_model(self):
        self.run_tool()

        # The field is calibrated: only rounds against the source's whole
        # image, which is also what is subtracted.
        self.assertIn("--whole-field-only", self.commands[1])
        self.assertIn("--whole-field-only", self.commands[4])

    def plan_of(self, command):
        """The rounds the recipe makes of a source's command."""

        spec = importlib.util.spec_from_file_location(
            "test_peel_recipe", SOURCE_DIR / "selfcal_recipe.py")
        recipe = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(recipe)
        self.assertEqual(Path(command[2]).name, "selfcal_recipe.py")
        plans = []
        with mock.patch.object(
            recipe, "run", lambda args, plan=None: plans.append(plan)
        ):
            recipe.main(command[3:])
        return [stage["id"] for stage in plans[0]]

    def test_the_recipe_takes_what_a_bright_source_is_given(self):
        self.run_tool()

        self.assertEqual(self.plan_of(self.commands[1]),
                         ["amplitude", "bandpass"])

    def test_the_reference_antennas_reach_each_sources_solves(self):
        # MeerKAT: the recipe's own default names uGMRT antennas.
        self.run_tool("--refant", "m059,m037")

        recipes = [command for (kind, _), command
                   in zip(self.calls, self.commands) if kind == "recipe"]
        self.assertEqual(len(recipes), 2)
        for command in recipes:
            self.assertEqual(option(command, "--refant"), "m059,m037")

    def test_without_reference_antennas_the_recipe_keeps_its_own(self):
        self.run_tool()

        self.assertNotIn("--refant", self.commands[1])

    def test_the_recipe_takes_what_a_faint_source_is_given(self):
        with mock.patch("builtins.print"):
            TOOL.main([
                "--target", "G1", "--ms", str(self.ms),
                "--models", "images/run_I", "--workdir", str(self.work),
                "--outdir", str(self.products),
                "--source", "1820", "3216", "1.5", "0", "--no-image",
            ])

        self.assertEqual(self.plan_of(self.commands[1]), ["amplitude"])

    def test_no_solve_towards_one_source_may_flag_the_fields_data(self):
        self.run_tool()

        self.assertIn("--keep-unsolved", self.commands[1])
        self.assertIn("--keep-unsolved", self.commands[4])

    def test_a_source_too_faint_gets_no_delays_and_no_bandpass(self):
        with mock.patch("builtins.print"):
            TOOL.main([
                "--target", "G1", "--ms", str(self.ms),
                "--models", "images/run_I", "--workdir", str(self.work),
                "--outdir", str(self.products),
                "--source", "1820", "3216", "1.5", "0", "--no-image",
            ])

        recipe = self.commands[1]
        self.assertIn("--phase-amplitude-only", recipe)
        self.assertNotIn("--bandpass-channels", recipe)
        self.assertIn("--keep-unsolved", recipe)

    def test_the_field_is_imaged_as_it_was_before(self):
        self.run_tool()

        image = self.commands[-1]
        self.assertEqual(option(image, "--imsize"), "5120")
        self.assertEqual(option(image, "--channels-out"), "12")
        self.assertIn("--no-local-rms", image)
        self.assertNotIn("--no-negative", image)
        self.assertNotIn("--write-model", image)

    def test_the_report_says_what_became_of_each_source(self):
        self.run_tool()

        report = self.report()
        self.assertEqual(sorted(report["sources"]), ["source1", "source2"])
        entry = report["sources"]["source2"]
        self.assertEqual((entry["x"], entry["y"]), (1641, 2103))
        self.assertEqual(entry["before"]["peak_Jy_beam"], 0.050)
        self.assertEqual(entry["after"]["most_negative_Jy_beam"], -0.0007)
        self.assertEqual(
            report["measurement_set"],
            str(self.work / "source2" / "source2_peeled.ms"),
        )
        self.assertIsNone(report["running"])
        self.assertFalse((self.work / TOOL.IMAGER.LOCK_NAME).exists())

    def test_a_failed_step_stops_and_the_same_command_goes_on(self):
        self.fail_at = 5        # the second source's recipe

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("source2", str(stopped.exception))
        self.assertEqual(self.report()["failed"]["step"], "source2: recipe")
        before = len(self.calls)
        self.fail_at = None

        self.run_tool()

        # The first source is not taken out again.
        self.assertEqual(self.calls[before:], [
            ("extract", "source2"), ("recipe", "source2"), ("peel", "source2"),
            ("image", "field"),
        ])
        self.assertNotIn("failed", self.report())

    def test_a_finished_run_given_again_does_nothing(self):
        self.run_tool()
        count = len(self.calls)

        self.run_tool()

        self.assertEqual(len(self.calls), count)

    def test_sources_taken_out_earlier_are_named_and_numbered_on(self):
        self.run_tool("--already-out", "100", "200")

        self.assertEqual(self.calls[0], ("extract", "source2"))
        self.assertEqual(
            option(self.commands[0], "--also-without", 2), ["100", "200"]
        )
        self.assertEqual(self.commands[3].count("--also-without"), 2)

    def test_sources_can_be_taken_from_an_image_of_the_field(self):
        image = self.root / "field-MFS-image.fits"
        image.write_bytes(b"")
        chosen = [{"x": 3162, "y": 4227, "cores": 6.38, "bins": 30,
                   "peak_mJy": 51.0, "offset_arcmin": 29.5}]
        with mock.patch("builtins.print"), \
                mock.patch.object(TOOL.IMAGER, "fits_plane", lambda path: None), \
                mock.patch.object(TOOL, "choose_sources",
                                  lambda *a: chosen):
            TOOL.main([
                "--target", "G1", "--ms", str(self.ms),
                "--models", "images/run_I", "--workdir", str(self.work),
                "--outdir", str(self.products), "--no-image",
                "--auto-sources", str(image),
            ])

        self.assertEqual(self.calls[0], ("extract", "source1"))
        self.assertEqual(option(self.commands[1], "--core-threshold"), "6.38")
        self.assertEqual(self.report()["chosen"], chosen)

    def test_a_field_with_nothing_to_peel_is_left_as_it_is(self):
        image = self.root / "field-MFS-image.fits"
        image.write_bytes(b"")
        with mock.patch("builtins.print"), \
                mock.patch.object(TOOL.IMAGER, "fits_plane", lambda path: None), \
                mock.patch.object(TOOL, "choose_sources", lambda *a: []):
            TOOL.main([
                "--target", "G1", "--ms", str(self.ms),
                "--models", "images/run_I", "--workdir", str(self.work),
                "--outdir", str(self.products), "--no-image",
                "--auto-sources", str(image),
            ])

        self.assertEqual(self.calls, [])
        report = self.report()
        self.assertEqual(report["chosen"], [])
        self.assertEqual(report["measurement_set"], str(self.ms))

    def test_sources_have_to_come_from_somewhere(self):
        with self.assertRaises(SystemExit) as stopped, \
                mock.patch("builtins.print"):
            TOOL.main([
                "--target", "G1", "--ms", str(self.ms),
                "--models", "images/run_I", "--workdir", str(self.work),
            ])

        self.assertIn("--source or --auto-sources", str(stopped.exception))

    def test_tidying_keeps_only_the_field_as_the_last_source_left_it(self):
        self.run_tool("--tidy")

        self.assertFalse((self.work / "source1" / "source1.ms").exists())
        self.assertFalse(
            (self.work / "source1" / "source1_peeled.ms").exists())
        self.assertTrue(
            (self.work / "source2" / "source2_peeled.ms").is_dir())
        self.assertTrue(self.ms.is_dir())
        self.assertIn("--tidy", self.commands[1])

    def test_a_tidied_run_given_again_does_not_take_the_sources_out_again(self):
        self.run_tool("--tidy")
        count = len(self.calls)

        self.run_tool("--tidy")

        self.assertEqual(len(self.calls), count)
        self.assertEqual(
            self.report()["measurement_set"],
            str(self.work / "source2" / "source2_peeled.ms"),
        )

    def test_the_image_can_be_left_out(self):
        self.run_tool("--no-image")

        self.assertNotIn(("image", "field"), self.calls)

    def test_a_folder_inside_the_products_is_refused(self):
        self.work = self.products / "peel"

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("inside the pipeline's products", str(stopped.exception))
        self.assertEqual(self.calls, [])

    def test_a_folder_of_another_field_is_not_continued(self):
        self.run_tool()
        other = self.root / "data" / "other.ms"
        other.mkdir()
        self.ms = other

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("holds a run on", str(stopped.exception))

    def test_a_dry_run_prints_the_steps_and_starts_nothing(self):
        printed = []
        with mock.patch("builtins.print", lambda *a, **k: printed.append(a[0])):
            TOOL.main([
                "--target", "G1", "--ms", str(self.ms),
                "--models", "images/run_I", "--workdir", str(self.work),
                "--source", "3162", "4227", "6", "30", "--dry-run",
            ])

        self.assertEqual(self.calls, [])
        self.assertFalse(self.work.exists())
        text = "\n".join(printed)
        self.assertIn("experiment_extract_source.py", text)
        self.assertIn("--peel", text)
        self.assertIn("selfcal_recipe.py", text)
        self.assertIn("experiment_wsclean_image.py", text)


if __name__ == "__main__":
    unittest.main()
