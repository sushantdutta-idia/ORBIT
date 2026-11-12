import importlib.util
import json
import os
import socket
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_tool():
    spec = importlib.util.spec_from_file_location(
        "test_run_targets_module", SOURCE_DIR / "run_targets.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TOOL = load_tool()


def option(command, name, count=1):
    index = command.index(name)
    values = command[index + 1:index + 1 + count]
    return values[0] if count == 1 else values


def fits_header(**cards):
    """The first block of a FITS file with these cards."""

    lines = ["SIMPLE  =                    T", "BITPIX  =                  -32"]
    lines += [f"{key:<8}= {value:>20}" for key, value in cards.items()]
    lines.append("END")
    text = "".join(f"{line:<80}" for line in lines)
    return text.ljust(2880).encode("ascii")


class StartTests(unittest.TestCase):

    def test_it_starts_beside_the_pipelines_own_inspect_module(self):
        source = (SOURCE_DIR / "run_targets.py").read_text()
        self.assertIn("sys.path.pop(0)", source.split("import argparse")[0])
        done = subprocess.run(
            [sys.executable, "run_targets.py", "--help"],
            cwd=SOURCE_DIR, capture_output=True, text=True,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("--peel-above", done.stdout)


class FlowTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.products = self.root / "uGMRT_Products"
        calibrated = self.products / "calibration" / "calibrated"
        calibrated.mkdir(parents=True)
        record = {}
        for target in ("G1", "G2", "G3"):
            (calibrated / f"{target}.ms").mkdir()
            record[target] = str(calibrated / f"{target}.ms")
        (self.products / "calibration" / "split_targets.json").write_text(
            json.dumps(record)
        )
        self.work = self.root / "run"
        self.calls = []
        self.commands = {}
        self.fail = set()
        # How many sources the image of each field gives to peel.
        self.sources = {"G1": 2, "G2": 0, "G3": 1}
        # The scans each target has, and those its data say to leave out.
        self.scans = {"G1": [3, 15, 29], "G2": [17], "G3": [5, 18]}
        self.left_out = {}
        test = self

        def scans(command):
            target = option(command, "--target")
            folder = Path(option(command, "--workdir"))
            folder.mkdir(parents=True, exist_ok=True)
            source = option(command, "--ms")
            forced = None
            if "--scans" in command:
                rest = command[command.index("--scans") + 1:]
                forced = []
                while rest and not rest[0].startswith("--"):
                    forced.append(int(rest.pop(0)))
            kept = forced or [
                scan for scan in test.scans[target]
                if scan not in test.left_out.get(target, [])
            ]
            average = (
                int(option(command, "--average-channels"))
                if "--average-channels" in command else 1
            )
            rule = {
                "least_share": float(option(command, "--least-share")),
                "most_excess": float(option(command, "--most-excess")),
                "scans": forced,
            }
            if average > 1:
                rule["average_channels"] = average
            data = source
            if kept != test.scans[target] or average > 1:
                data = str(folder / f"{target}.ms")
                Path(data).mkdir(exist_ok=True)
            (folder / "scans_report.json").write_text(json.dumps({
                "input": source,
                "rule": rule,
                "kept": kept,
                "left_out": [
                    scan for scan in test.scans[target] if scan not in kept],
                "measurement_set": data,
                "reading": [
                    f"scan {scan}: 22.0 min; kept."
                    for scan in test.scans[target]
                ] + [f"{target}: every scan is kept."],
            }))

        def recipe(command):
            target = option(command, "--target")
            folder = Path(option(command, "--workdir")) / target
            final = folder / "images" / "final" / target
            final.mkdir(parents=True, exist_ok=True)
            (final / "wsclean_robust_0_I-MFS-image.fits").write_bytes(b"")
            (folder / "round5.ms").mkdir(exist_ok=True)
            (folder / "recipe_report.json").write_text(json.dumps({
                "input": option(command, "--ms"),
                "done": ["copy", "final:image"],
                "calibrated_measurement_set": str(folder / "round5.ms"),
                "final": {"image_rms_Jy_beam": 25.2e-6,
                          "quietest_tenth_Jy_beam": 22.6e-6,
                          "data_noise_Jy_beam": 21.5e-6},
                "reading": ["final: image 25.2 uJy/beam"],
            }))

        def image(command):
            target = option(command, "--target")
            label = TOOL.robust_label(
                option(command, "--robust"),
                option(command, "--taper") if "--taper" in command else 0)
            folder = Path(option(command, "--workdir")) / target
            folder.mkdir(parents=True, exist_ok=True)
            name = folder / f"wsclean_{label}_I"
            Path(str(name) + "-MFS-image.fits").write_bytes(fits_header(
                BMAJ=3.1 / 3600, BMIN=2.4 / 3600))
            # What WSClean leaves for each sub-band: the image and model
            # are kept for the primary beam, the rest is removed.
            for kind in ("image", "model", "residual", "psf", "dirty"):
                Path(f"{name}-0003-{kind}.fits").write_bytes(b"")
            (folder / f"wsclean_{label}.json").write_text(json.dumps({
                "stokes_i": {"robust_rms_Jy_beam": 17.0e-6,
                             "grid": {"quietest_tenth_rms": 16.2e-6}},
                "stokes_v": {"robust_rms_Jy_beam": 15.9e-6},
            }))

        def run_command(command, log_path):
            script = Path(command[2]).name
            test.assertTrue(Path(log_path).parent.is_dir())
            if script == "selfcal_recipe.py":
                target = option(command, "--target")
                stage = Path(option(command, "--workdir")).name
                if stage == TOOL.TARGET_FOLDER:
                    stage = "target"
            elif script == "peel_sources.py":
                target, stage = option(command, "--target"), "peel"
            elif script == "select_scans.py":
                target, stage = option(command, "--target"), "scans"
            elif script == "experiment_wsclean_image.py":
                target = option(command, "--target")
                stage = TOOL.robust_label(
                    option(command, "--robust"),
                    option(command, "--taper") if "--taper" in command else 0)
            else:
                target = Path(option(command, "--name")).parts[
                    len(test.work.parts)]
                stage = "primary_beam"
                found = Path(option(command, "--name")).name
                if found != "wsclean_robust_0_I":
                    stage = found[len("wsclean_"):-len("_I")] + "_primary_beam"
            test.calls.append((target, stage))
            test.commands[(target, stage)] = command
            test.assertTrue(Path(option(command, "--ms")).is_dir()
                            if "--ms" in command else True)
            if (target, stage) in test.fail:
                raise TOOL.IMAGER.CommandFailed(f"{script} exited with 1", 1)
            if script == "selfcal_recipe.py":
                recipe(command)
            elif script == "select_scans.py":
                scans(command)
            elif script == "experiment_wsclean_image.py":
                image(command)
            elif script == "peel_sources.py":
                folder = Path(option(command, "--workdir"))
                folder.mkdir(parents=True, exist_ok=True)
                count = test.sources[target]
                peeled = folder / "source9" / "source9_peeled.ms"
                if count:
                    peeled.mkdir(parents=True)
                (folder / "peel_report.json").write_text(json.dumps({
                    "chosen": [{"x": 1, "y": 2}] * count,
                    "measurement_set": (
                        str(peeled) if count else option(command, "--ms")),
                }))
            else:
                meerkat = (
                    "--telescope" in command
                    and option(command, "--telescope") == "meerkat"
                )
                Path(option(command, "--name")
                     + ("-MFS-image-pb.fits" if meerkat
                        else "-MFS-image.pbcor.fits")).write_bytes(b"")

        patch = mock.patch.object(TOOL.IMAGER, "run_command", run_command)
        patch.start()
        self.addCleanup(patch.stop)

    def run_tool(self, *extra):
        with mock.patch("builtins.print"):
            return TOOL.main([
                "--products", str(self.products), "--workdir", str(self.work),
                "--channel-range", "51", "768", *extra,
            ])

    def report(self):
        return json.loads((self.work / TOOL.REPORT_NAME).read_text())

    def test_a_faint_field_can_go_without_the_rounds_within_the_scans(self):
        self.run_tool("--targets", "G1", "--time-phase-solint", "none")

        field = self.commands[("G1", "field")]
        self.assertNotIn("--time-phase-solint", field)
        self.assertIn("--amplitude-cores", field)

    def test_averaged_channels_are_counted_in_the_averaged_data(self):
        # Band 5: data in channels 122-3973 of 4096, averaged by 4.
        self.run_tool("--channel-range", "122", "3974",
                      "--average-channels", "4")

        scans = self.commands[("G1", "scans")]
        self.assertEqual(option(scans, "--channel-range", 2), ["122", "3974"])
        self.assertEqual(option(scans, "--average-channels"), "4")
        for key, command in self.commands.items():
            if key[1] in ("scans", "primary_beam") or "--channel-range" not in command:
                continue
            with self.subTest(step=key):
                self.assertEqual(
                    option(command, "--channel-range", 2), ["31", "993"])
        # Every target goes on with its averaged copy, also G2 and G3,
        # which keep all their scans.
        field = self.commands[("G2", "field")]
        self.assertEqual(
            option(field, "--ms"), str(self.work / "G2" / "scans" / "G2.ms"))

    def test_the_averaged_range_keeps_whole_channels_only(self):
        args = TOOL.argparse.Namespace(channel_range=[122, 3974],
                                       average_channels=4)
        self.assertEqual(TOOL.averaged_channel_range(args), (31, 993))
        args.average_channels = 1
        self.assertEqual(TOOL.averaged_channel_range(args), (122, 3974))
        args.channel_range, args.average_channels = [1, 4], 4
        with self.assertRaises(SystemExit):
            TOOL.averaged_channel_range(args)

    def test_no_averaging_leaves_the_scans_rule_as_it_was(self):
        args = TOOL.argparse.Namespace(least_share=0.1, most_excess=1.4,
                                       forced={}, average_channels=1)
        self.assertEqual(TOOL.scans_rule(args, "G1"), {
            "least_share": 0.1, "most_excess": 1.4, "scans": None})
        args.average_channels = 4
        self.assertEqual(TOOL.scans_rule(args, "G1")["average_channels"], 4)

    def test_averaging_must_be_one_or_more(self):
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.run_tool("--average-channels", "0")

    def test_every_target_goes_through_the_four_steps_in_turn(self):
        self.run_tool()

        self.assertEqual(self.calls, [
            # Every target's scans, before the first image is started.
            ("G1", "scans"), ("G2", "scans"), ("G3", "scans"),
            ("G1", "field"), ("G1", "peel"), ("G1", "target"),
            ("G1", "primary_beam"),
            # A field with nothing to peel is imaged once.
            ("G2", "field"), ("G2", "peel"), ("G2", "primary_beam"),
            ("G3", "field"), ("G3", "peel"), ("G3", "target"),
            ("G3", "primary_beam"),
        ])
        report = self.report()
        self.assertEqual(report["targets"]["G1"]["sources_out"], 2)
        self.assertEqual(report["targets"]["G2"]["sources_out"], 0)
        self.assertTrue(
            report["targets"]["G1"]["corrected"].endswith(
                "run/G1/target_whole_field/G1/images/final/G1/"
                "wsclean_robust_0_I-MFS-image.pbcor.fits")
        )
        self.assertTrue(
            report["targets"]["G2"]["corrected"].endswith(
                "run/G2/field/G2/images/final/G2/"
                "wsclean_robust_0_I-MFS-image.pbcor.fits")
        )
        self.assertIn("G1: 25.2 uJy/beam (22.6 in the quietest tenth) "
                      "against 21.5 for the data; 2 sources taken out",
                      report["reading"][0])
        self.assertEqual(report["running"], [])
        # Each target's own record, from which the summary is made.
        own = json.loads((self.work / "G1" / TOOL.TARGET_REPORT).read_text())
        self.assertEqual(own["sources_out"], 2)
        self.assertEqual(own, report["targets"]["G1"])

    def test_images_of_other_weightings_follow_the_primary_beam(self):
        self.run_tool("--targets", "G1", "--extra-robust", "0.5")

        self.assertEqual(self.calls[-3:], [
            ("G1", "primary_beam"), ("G1", "robust_0p5"),
            ("G1", "robust_0p5_primary_beam"),
        ])
        command = self.commands[("G1", "robust_0p5")]
        calibrated = self.work / "G1" / TOOL.TARGET_FOLDER / "G1" / "round5.ms"
        images = (self.work / "G1" / TOOL.TARGET_FOLDER / "G1" / "images"
                  / "final_robust_0p5")
        self.assertEqual(option(command, "--ms"), str(calibrated))
        self.assertEqual(option(command, "--workdir"), str(images))
        self.assertEqual(option(command, "--robust"), "0.5")
        self.assertEqual(option(command, "--imsize"), "5120")
        self.assertEqual(option(command, "--channel-range", 2), ["51", "768"])
        self.assertEqual(option(command, "--scales", 4), ["0", "5", "15", "30"])
        self.assertIn("--no-local-rms", command)
        beam = self.commands[("G1", "robust_0p5_primary_beam")]
        self.assertEqual(option(beam, "--name"),
                         str(images / "G1" / "wsclean_robust_0p5_I"))

        entry = self.report()["targets"]["G1"]
        extra = entry["extra_images"]["robust_0p5"]
        self.assertEqual(extra["beam_arcsec"], [3.1, 2.4])
        self.assertAlmostEqual(extra["image_rms_Jy_beam"], 17.0e-6)
        self.assertTrue(extra["corrected"].endswith(
            "final_robust_0p5/G1/wsclean_robust_0p5_I-MFS-image.pbcor.fits"))
        # The robust-0 image stays the target's product.
        self.assertTrue(entry["corrected"].endswith(
            "images/final/G1/wsclean_robust_0_I-MFS-image.pbcor.fits"))
        self.assertIn('robust 0.5: 17.0 uJy/beam (16.2 in the quietest '
                      'tenth), beam 3.10" x 2.40"', self.report()["reading"][0])
        # Of the sub-bands only what the primary beam needs is kept.
        left = sorted(path.name for path in (images / "G1").glob("*-0003-*"))
        self.assertEqual(left, ["wsclean_robust_0p5_I-0003-image.fits",
                                "wsclean_robust_0p5_I-0003-model.fits"])

    def test_a_finished_target_given_the_option_later_makes_only_that(self):
        self.run_tool("--targets", "G2")
        self.calls.clear()
        self.run_tool("--targets", "G2", "--extra-robust", "0.5", "-0.5")

        self.assertEqual(self.calls, [
            ("G2", "robust_0p5"), ("G2", "robust_0p5_primary_beam"),
            ("G2", "robust_m0p5"), ("G2", "robust_m0p5_primary_beam"),
        ])
        self.calls.clear()
        self.run_tool("--targets", "G2", "--extra-robust", "0.5", "-0.5")
        self.assertEqual(self.calls, [])

    def test_tapered_images_are_made_on_a_coarser_grid_of_the_same_field(self):
        self.run_tool("--targets", "G2", "--extra-robust", "0.5",
                      "--extra-taper", "6", "15")

        self.assertEqual(self.calls[-6:], [
            ("G2", "robust_0p5"), ("G2", "robust_0p5_primary_beam"),
            ("G2", "robust_0_taper6"), ("G2", "robust_0_taper6_primary_beam"),
            ("G2", "robust_0_taper15"), ("G2", "robust_0_taper15_primary_beam"),
        ])
        wide = self.commands[("G2", "robust_0_taper15")]
        self.assertEqual(option(wide, "--taper"), "15")
        self.assertEqual(option(wide, "--robust"), "0")
        # Five pixels to the beam, and the field of the 5120 x 1" image.
        self.assertEqual(option(wide, "--cell"), "3")
        self.assertEqual(option(wide, "--imsize"), "1706")
        narrow = self.commands[("G2", "robust_0_taper6")]
        self.assertEqual(option(narrow, "--cell"), "1.2")
        self.assertEqual(option(narrow, "--imsize"), "4266")
        self.assertTrue(option(wide, "--workdir").endswith(
            "G2/field/G2/images/final_robust_0_taper15"))
        entry = self.report()["targets"]["G2"]
        self.assertEqual(sorted(entry["extra_images"]),
                         ["robust_0_taper15", "robust_0_taper6", "robust_0p5"])
        self.assertEqual(entry["extra_images"]["robust_0_taper15"]["taper_arcsec"],
                         15.0)
        self.assertIn('taper 15": 17.0 uJy/beam', self.report()["reading"][0])

    def test_wider_tapered_fields_go_to_folders_of_their_own(self):
        self.run_tool("--targets", "G2", "--extra-taper", "15",
                      "--taper-field-scale", "3")

        self.assertEqual(self.calls[-2:], [
            ("G2", "robust_0_taper15"), ("G2", "robust_0_taper15_primary_beam"),
        ])
        wide = self.commands[("G2", "robust_0_taper15")]
        # Three times the 5120 x 1" field, in 3" cells.
        self.assertEqual(option(wide, "--cell"), "3")
        self.assertEqual(option(wide, "--imsize"), "5120")
        self.assertTrue(option(wide, "--workdir").endswith(
            "G2/field/G2/images/final_robust_0_taper15_field3"))
        entry = self.report()["targets"]["G2"]["extra_images"]
        self.assertEqual(list(entry), ["robust_0_taper15_field3"])
        self.assertEqual(entry["robust_0_taper15_field3"]["field_arcmin"], 256.0)
        self.assertTrue(entry["robust_0_taper15_field3"]["corrected"].endswith(
            "final_robust_0_taper15_field3/G2/"
            "wsclean_robust_0_taper15_I-MFS-image.pbcor.fits"))

    def test_a_field_narrower_than_the_main_one_is_refused(self):
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.run_tool("--extra-taper", "15", "--taper-field-scale", "0.5")

    def test_the_report_lists_the_extra_images_this_run_was_asked_for(self):
        self.run_tool("--targets", "G2", "--extra-robust", "0.5")
        self.run_tool("--targets", "G2", "--extra-taper", "15")
        entry = self.report()["targets"]["G2"]
        self.assertEqual(list(entry["extra_images"]), ["robust_0_taper15"])

    def test_untapered_robust_above_one_is_refused(self):
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.run_tool("--extra-robust", "2")
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.run_tool("--extra-taper", "0")
        self.assertEqual(self.calls, [])

    def test_the_other_weightings_are_cleaned_as_the_recipe_cleans(self):
        spec = importlib.util.spec_from_file_location(
            "test_run_targets_recipe", SOURCE_DIR / "selfcal_recipe.py")
        recipe = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(recipe)
        captured = []
        with mock.patch.object(recipe, "run",
                               lambda args, plan=None: captured.append(args)):
            recipe.main(["--target", "G1", "--ms", "data.ms", "--imsize",
                         "5120", "--cell", "0.5", "--channel-range", "31",
                         "993", "--casa", "python3"])
        args = captured[0]
        self.assertEqual(tuple(args.scales), TOOL.FINAL_SCALES)
        self.assertEqual(args.auto_threshold, TOOL.FINAL_AUTO_THRESHOLD)
        self.assertEqual(args.nmiter, TOOL.FINAL_NMITER)

    def test_the_field_starts_from_the_pipelines_calibrated_data(self):
        self.run_tool("--targets", "G1")

        field = self.commands[("G1", "field")]
        self.assertEqual(
            option(field, "--ms"),
            str(self.products / "calibration" / "calibrated" / "G1.ms"),
        )
        self.assertEqual(option(field, "--imsize"), "5120")
        # Data with no self-calibration yet: phases within the scans too,
        # and cores down to a fraction of this field's own peak.
        self.assertEqual(option(field, "--time-phase-solint"), "1min")
        # ... and amplitudes against cores, which a usual clean of data
        # with no amplitude self-calibration cannot find.
        self.assertIn("--amplitude-cores", field)
        self.assertEqual(option(field, "--core-fraction"), "0.13")
        self.assertEqual(option(field, "--channels-out"), "12")
        self.assertEqual(option(field, "--outdir"), str(self.products))
        self.assertIn("--tidy", field)
        self.assertNotIn("--keep-unsolved", field)

    def test_peeling_uses_the_fields_own_image_and_models(self):
        self.run_tool("--targets", "G1")

        peel = self.commands[("G1", "peel")]
        final = self.work / "G1" / "field" / "G1" / "images" / "final" / "G1"
        self.assertEqual(option(peel, "--models"),
                         str(final / "wsclean_robust_0_I"))
        self.assertEqual(option(peel, "--auto-sources"),
                         str(final / "wsclean_robust_0_I-MFS-image.fits"))
        self.assertEqual(option(peel, "--ms"),
                         str(self.work / "G1" / "field" / "G1" / "round5.ms"))
        self.assertEqual(option(peel, "--peel-above"), "10")
        self.assertEqual(option(peel, "--keep-radius"), "4")
        self.assertIn("--no-image", peel)
        self.assertIn("--tidy", peel)

    def test_what_is_left_is_calibrated_on_the_inner_field_and_imaged_whole(self):
        self.run_tool("--targets", "G1")

        target = self.commands[("G1", "target")]
        self.assertEqual(
            option(target, "--ms"),
            str(self.work / "G1" / "peel" / "source9" / "source9_peeled.ms"),
        )
        self.assertEqual(option(target, "--imsize"), "3072")
        self.assertEqual(option(target, "--final-imsize"), "5120")
        self.assertIn("--keep-unsolved", target)
        self.assertNotIn("--time-phase-solint", target)
        self.assertNotIn("--amplitude-cores", target)
        # Calibrated data with their bright sources out: no round against
        # a cores-only model, which would hold a tenth of the field.
        self.assertIn("--whole-field-only", target)
        self.assertEqual(
            option(target, "--workdir"),
            str(self.work / "G1" / "target_whole_field"),
        )

    def test_what_an_earlier_run_solved_against_cores_is_done_again(self):
        self.run_tool("--targets", "G1")
        # As an earlier version left it: the stage finished in "target".
        new = self.work / "G1" / "target_whole_field"
        old = self.work / "G1" / "target"
        new.rename(old)
        before = len(self.calls)

        self.run_tool("--targets", "G1")

        self.assertEqual(self.calls[before:],
                         [("G1", "target"), ("G1", "primary_beam")])
        # The earlier result is neither used nor removed.
        self.assertTrue((old / "G1" / "recipe_report.json").is_file())
        self.assertTrue(self.report()["targets"]["G1"]["corrected"].endswith(
            "run/G1/target_whole_field/G1/images/final/G1/"
            "wsclean_robust_0_I-MFS-image.pbcor.fits"))
        beam = self.commands[("G1", "primary_beam")]
        self.assertEqual(option(beam, "--band"), "B4")
        self.assertEqual(option(beam, "--pbmin"), "0.3")

    def test_once_one_has_finished_a_failure_does_not_stop_the_others(self):
        self.fail = {("G2", "peel")}

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("1 of 3 targets stopped: G2", str(stopped.exception))
        self.assertEqual(self.calls[-1], ("G3", "primary_beam"))
        report = self.report()
        self.assertEqual(report["targets"]["G2"]["failed"]["stage"], "peel")
        self.assertIn("G2: stopped at 'peel'.", report["reading"])
        self.assertIn("corrected", report["targets"]["G3"])

    def test_a_failure_before_any_target_has_finished_stops_the_run(self):
        self.fail = {("G1", "peel")}

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        # What stopped the first would as likely stop the next, hours in.
        self.assertEqual(self.calls[3:], [("G1", "field"), ("G1", "peel")])
        self.assertIn("1 of 3 targets stopped: G1", str(stopped.exception))
        self.assertEqual(self.report()["targets"]["G2"]["stages"].keys(),
                         {"scans"})

    def test_the_same_command_goes_on_with_what_is_not_done(self):
        self.fail = {("G1", "peel")}
        with self.assertRaises(SystemExit):
            self.run_tool()
        before = len(self.calls)
        self.fail = set()

        self.run_tool()

        # Neither the scans nor the field of G1 are made again; then the
        # others follow.
        self.assertEqual(self.calls[before:before + 4], [
            ("G1", "peel"), ("G1", "target"), ("G1", "primary_beam"),
            ("G2", "field"),
        ])
        self.assertNotIn("failed", self.report()["targets"]["G1"])
        self.assertEqual(self.calls[-1], ("G3", "primary_beam"))

    def test_a_finished_run_given_again_does_nothing(self):
        self.run_tool()
        count = len(self.calls)

        self.run_tool()

        self.assertEqual(len(self.calls), count)

    def test_missing_data_stop_only_that_target(self):
        (self.products / "calibration" / "calibrated" / "G2.ms").rmdir()

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("G2", str(stopped.exception))
        self.assertIn("are missing",
                      self.report()["targets"]["G2"]["failed"]["error"])
        self.assertIn(("G3", "primary_beam"), self.calls)

    def test_targets_can_be_named_and_ordered(self):
        self.run_tool("--targets", "G3", "G1")

        self.assertEqual(self.calls[:3],
                         [("G3", "scans"), ("G1", "scans"), ("G3", "field")])
        self.assertNotIn(("G2", "field"), self.calls)
        self.assertNotIn(("G2", "scans"), self.calls)

    def test_scans_left_out_give_the_field_a_measurement_set_of_the_rest(self):
        self.left_out = {"G1": [3, 29]}

        self.run_tool("--targets", "G1", "G3")

        kept = self.work / "G1" / "scans" / "G1.ms"
        self.assertEqual(
            option(self.commands[("G1", "field")], "--ms"), str(kept))
        # A target whose scans are all kept is imaged from the products.
        self.assertEqual(
            option(self.commands[("G3", "field")], "--ms"),
            str(self.products / "calibration" / "calibrated" / "G3.ms"),
        )
        report = self.report()
        self.assertEqual(report["targets"]["G1"]["scans"]["kept"], [15])
        self.assertEqual(report["targets"]["G1"]["scans"]["left_out"], [3, 29])
        self.assertIn("without scans 3, 29; ", report["reading"][0])
        self.assertNotIn("without scan", report["reading"][1])
        scans = self.commands[("G1", "scans")]
        self.assertEqual(option(scans, "--least-share"), "0.1")
        self.assertEqual(option(scans, "--most-excess"), "1.4")
        self.assertEqual(option(scans, "--channel-range", 2), ["51", "768"])
        self.assertEqual(option(scans, "--outdir"), str(self.products))
        self.assertNotIn("--scans", scans)

    def test_scans_can_be_named_for_a_target(self):
        self.run_tool("--targets", "G1", "G3", "--scans", "G1=15")

        scans = self.commands[("G1", "scans")]
        self.assertEqual(option(scans, "--scans"), "15")
        self.assertNotIn("--scans", self.commands[("G3", "scans")])
        self.assertEqual(self.report()["targets"]["G1"]["scans"]["kept"], [15])

    def test_scans_named_wrongly_are_refused_before_anything_starts(self):
        for value, said in (("G1", "give a target and its scans"),
                            ("G1=a", "give a target and its scans"),
                            ("G9=15", "G9 is not among the calibrated")):
            with self.assertRaises(SystemExit) as stopped:
                self.run_tool("--scans", value)
            self.assertIn(said, str(stopped.exception))
        self.assertEqual(self.calls, [])

    def test_the_scans_are_printed_before_the_first_image(self):
        printed = []
        order = []
        original = TOOL.IMAGER.run_command

        def run_command(command, log_path):
            order.append(len(printed))
            original(command, log_path)

        with mock.patch.object(TOOL.IMAGER, "run_command", run_command), \
                mock.patch("builtins.print",
                           lambda *a, **k: printed.append(a[0])):
            TOOL.main([
                "--products", str(self.products), "--workdir", str(self.work),
                "--channel-range", "51", "768",
            ])

        table = printed.index("\n##### The scans #####")
        # Three scan steps before the table, the first field after it.
        self.assertLessEqual(order[2], table)
        self.assertGreater(order[3], table)
        self.assertIn("G1: scan 3: 22.0 min; kept.", printed[table + 1])
        self.assertIn("G1: every scan is kept.", printed[table + 1])

    def test_a_change_of_the_scans_does_not_go_on_with_the_old_field(self):
        self.fail = {("G1", "peel")}
        with self.assertRaises(SystemExit):
            self.run_tool("--targets", "G1")
        self.fail = set()
        before = len(self.calls)

        with self.assertRaises(SystemExit):
            self.run_tool("--targets", "G1", "--scans", "G1=15")

        # The scans are decided again; the field made from all of them is
        # neither used nor removed.
        self.assertEqual(self.calls[before:], [("G1", "scans")])
        failed = self.report()["targets"]["G1"]["failed"]
        self.assertEqual(failed["stage"], "field")
        self.assertIn("the scans to use have changed", failed["error"])
        self.assertTrue((self.work / "G1" / "field" / "G1").is_dir())

    def test_a_failure_of_the_first_scans_stops_the_run(self):
        self.fail = {("G1", "scans")}

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertEqual(self.calls, [("G1", "scans")])
        self.assertIn("1 of 3 targets stopped: G1", str(stopped.exception))
        self.assertEqual(
            self.report()["targets"]["G1"]["failed"]["stage"], "scans")

    def test_an_unknown_target_is_said_with_those_there_are(self):
        with self.assertRaises(SystemExit) as stopped:
            self.run_tool("--targets", "G9")

        self.assertIn("Available: G1, G2, G3", str(stopped.exception))
        self.assertEqual(self.calls, [])

    def test_a_folder_inside_the_products_is_refused(self):
        self.work = self.products / "run"

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("inside the pipeline's products", str(stopped.exception))

    def test_products_without_calibrated_targets_are_said_so(self):
        (self.products / "calibration" / "split_targets.json").unlink()

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("No calibrated targets", str(stopped.exception))

    def test_nothing_in_the_products_is_written_or_removed(self):
        before = sorted(
            str(path) for path in self.products.rglob("*"))

        self.run_tool()

        self.assertEqual(
            sorted(str(path) for path in self.products.rglob("*")), before)

    def held(self, folder):
        """A claim on a folder by a run that is still going: this one."""

        folder.mkdir(parents=True, exist_ok=True)
        (folder / "running.lock").write_text(json.dumps({
            "pid": os.getpid(), "host": socket.gethostname(),
            "since": "2026-10-16 08:00:00 UTC",
        }))

    def test_a_target_another_run_holds_is_left_to_it(self):
        self.held(self.work / "G2")

        self.run_tool()

        self.assertNotIn("G2", {target for target, _ in self.calls})
        self.assertIn(("G3", "primary_beam"), self.calls)
        # Its claim is not taken or removed.
        self.assertTrue((self.work / "G2" / "running.lock").is_file())

    def test_runs_of_one_target_each_share_the_folder(self):
        self.run_tool("--targets", "G1")
        self.run_tool("--targets", "G3")

        report = self.report()
        self.assertIn("corrected", report["targets"]["G1"])
        self.assertIn("corrected", report["targets"]["G3"])
        self.assertTrue((self.work / "G1" / TOOL.TARGET_REPORT).is_file())
        self.assertTrue((self.work / "G3" / TOOL.TARGET_REPORT).is_file())
        # Each run held only its own target, and let it go.
        self.assertFalse((self.work / "G1" / "running.lock").exists())
        self.assertFalse((self.work / "running.lock").exists())

    def test_the_summary_takes_what_another_run_has_recorded(self):
        self.run_tool("--targets", "G1")
        record = json.loads((self.work / "G1" / TOOL.TARGET_REPORT).read_text())
        record["sources_out"] = 7      # as another run would leave it
        (self.work / "G1" / TOOL.TARGET_REPORT).write_text(json.dumps(record))

        self.run_tool("--targets", "G3")

        self.assertEqual(self.report()["targets"]["G1"]["sources_out"], 7)

    def test_a_whole_folder_held_by_an_earlier_version_is_refused(self):
        self.held(self.work)

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("held as a whole by another run", str(stopped.exception))
        self.assertEqual(self.calls, [])

    def test_a_whole_folder_claim_that_is_gone_is_passed_over(self):
        self.work.mkdir(parents=True)
        (self.work / "running.lock").write_text(json.dumps({
            "pid": 999999999, "host": socket.gethostname(),
            "since": "2026-10-11 20:00:20 UTC",
        }))

        self.run_tool("--targets", "G1")

        self.assertIn(("G1", "primary_beam"), self.calls)

    def test_threads_memory_and_parallel_cleaning_reach_every_script(self):
        with mock.patch.dict(os.environ, {"SLURM_CPUS_PER_TASK": "16"}):
            self.run_tool("--targets", "G1", "--memory-gb", "100",
                          "--parallel-deconvolution", "2560")

        for stage in ("field", "peel", "target"):
            command = self.commands[("G1", stage)]
            self.assertEqual(option(command, "--threads"), "16", stage)
            self.assertEqual(option(command, "--memory-gb"), "100", stage)
            self.assertEqual(
                option(command, "--parallel-deconvolution"), "2560", stage)

    def test_without_asking_nothing_is_cleaned_in_parallel(self):
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("SLURM_CPUS_PER_TASK", None)
            self.run_tool("--targets", "G1")

        field = self.commands[("G1", "field")]
        self.assertNotIn("--parallel-deconvolution", field)
        self.assertEqual(option(field, "--threads"), "6")
        self.assertEqual(option(field, "--memory-gb"), "20")

    def test_a_dry_run_prints_the_steps_and_starts_nothing(self):
        printed = []
        with mock.patch("builtins.print", lambda *a, **k: printed.append(a[0])):
            TOOL.main([
                "--products", str(self.products), "--workdir", str(self.work),
                "--channel-range", "51", "768", "--dry-run",
            ])

        self.assertEqual(self.calls, [])
        self.assertFalse(self.work.exists())
        text = "\n".join(printed)
        for name in ("# G1", "# G3", "select_scans.py", "selfcal_recipe.py",
                     "peel_sources.py", "correct_primary_beam.py"):
            self.assertIn(name, text)

    def test_ugmrt_still_needs_its_channel_range(self):
        with mock.patch("sys.stderr"), mock.patch("builtins.print"), \
                self.assertRaises(SystemExit):
            TOOL.main(["--products", str(self.products),
                       "--workdir", str(self.work)])
        self.assertEqual(self.calls, [])


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


class MeerKATTests(unittest.TestCase):
    """A MeerKAT field takes its settings from what the pipeline recorded.

    As 1708779077's products have them: MeerKAT L, XX and YY, the
    calibration's reference antenna m059, and 1024 channels of which 0-52
    and 974-1023 lie outside 900-1670 MHz.
    """

    PREFERRED = "m059,m037,m031,m027,m020,m001"

    def setUp(self):
        self.flow = FlowTests("test_every_target_goes_through_the_four_steps_in_turn")
        self.flow.setUp()
        self.addCleanup(self.flow.doCleanups)
        products = self.products = self.flow.products
        write_json(products / "inspection" / "inspection_report.json", {
            "telescope": {"name": "meerkat"},
            "band": {"name": "MEERKAT_L"},
            "correlations": {"setups": {"0": ["XX", "YY"]},
                             "basis": "linear"},
            "reference_antenna": {"recommended": "m059"},
            "spws": {"0": {"id": 0, "nchan": 1024}},
        })
        write_json(products / "calibration" / "calibration_report.json", {
            "telescope": "meerkat", "band": "MEERKAT_L",
            "reference_antenna": "m059",
        })
        # As 1708779077's flagging recorded them: outside 900-1670 MHz, and
        # 1% of each end of the band.
        self.guaranteed = (
            products / "flagging" / "guaranteed_range_mask_report.json")
        write_json(self.guaranteed, {"spws": [{
            "spw": "0", "total_channels": 1024,
            "channel_ranges": [[0, 52], [974, 1023]],
        }]})
        self.edges = products / "flagging" / "edge_flagging_report.json"
        write_json(self.edges, {
            "telescope": "meerkat", "band": "MEERKAT_L",
            "spw_string": "0:0~52;974~1023,0:0~9;1014~1023",
        })

    def run_tool(self, *extra):
        with mock.patch("builtins.print"):
            return TOOL.main(["--products", str(self.products),
                              "--workdir", str(self.flow.work),
                              "--targets", "G1", *extra])

    def test_the_field_is_made_on_the_bands_own_grid_and_channels(self):
        self.run_tool()

        commands = self.flow.commands
        for stage in ("field", "peel"):
            with self.subTest(stage=stage):
                command = commands[("G1", stage)]
                self.assertEqual(option(command, "--cell"), "1.5")
                self.assertEqual(
                    option(command, "--channel-range", 2), ["53", "974"])
        self.assertEqual(option(commands[("G1", "field")], "--imsize"), "6144")
        self.assertEqual(
            option(commands[("G1", "peel")], "--field-imsize"), "6144")
        self.assertEqual(
            option(commands[("G1", "scans")], "--channel-range", 2),
            ["53", "974"])

    def test_every_solve_is_told_meerkat_antennas(self):
        self.run_tool()

        for stage in ("field", "peel", "target"):
            with self.subTest(stage=stage):
                self.assertEqual(
                    option(self.flow.commands[("G1", stage)], "--refant"),
                    self.PREFERRED)

    def test_the_primary_beam_is_katbeams_with_the_products_masks(self):
        self.run_tool()

        beam = self.flow.commands[("G1", "primary_beam")]
        self.assertEqual(option(beam, "--telescope"), "meerkat")
        self.assertEqual(option(beam, "--band"), "MEERKAT_L")
        self.assertEqual(option(beam, "--products"), str(self.products))
        report = json.loads(
            (self.flow.work / TOOL.REPORT_NAME).read_text())
        self.assertTrue(report["targets"]["G1"]["corrected"].endswith(
            "wsclean_robust_0_I-MFS-image-pb.fits"))

    def test_options_given_win_over_the_products(self):
        self.run_tool("--imsize", "4096", "--cell", "2",
                      "--channel-range", "60", "900", "--refant", "m000")

        field = self.flow.commands[("G1", "field")]
        self.assertEqual(option(field, "--imsize"), "4096")
        self.assertEqual(option(field, "--cell"), "2")
        self.assertEqual(option(field, "--channel-range", 2), ["60", "900"])
        self.assertEqual(option(field, "--refant"), "m000")

    def test_the_calibrations_antenna_goes_first(self):
        write_json(self.products / "calibration" / "calibration_report.json",
                   {"telescope": "meerkat", "band": "MEERKAT_L",
                    "reference_antenna": "m010"})
        self.run_tool()

        self.assertEqual(
            option(self.flow.commands[("G1", "field")], "--refant"),
            "m010," + self.PREFERRED)

    def test_a_telescope_the_products_contradict_is_refused(self):
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.run_tool("--telescope", "ugmrt")
        self.assertEqual(self.flow.calls, [])

    def test_a_band_without_a_measured_beam_is_refused(self):
        write_json(self.products / "calibration" / "calibration_report.json",
                   {"telescope": "meerkat", "band": "MEERKAT_S0"})
        write_json(self.products / "inspection" / "inspection_report.json",
                   {"telescope": {"name": "meerkat"},
                    "band": {"name": "MEERKAT_S0"}})
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.run_tool()
        self.assertEqual(self.flow.calls, [])

    def test_without_the_flagging_records_a_channel_range_is_asked_for(self):
        self.guaranteed.unlink()
        self.edges.unlink()
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.run_tool()
        self.run_tool("--channel-range", "53", "974")
        self.assertIn(("G1", "field"), self.flow.calls)

    def test_the_channels_are_read_from_the_flagging_records(self):
        self.assertEqual(TOOL.recorded_channel_range(self.products),
                         (53, 974))
        # Either record alone gives the same edges here.
        self.edges.unlink()
        self.assertEqual(TOOL.recorded_channel_range(self.products),
                         (53, 974))
        write_json(self.guaranteed, {"spws": [
            {"spw": "0", "total_channels": 8, "channel_ranges": []},
            {"spw": "1", "total_channels": 8, "channel_ranges": []},
        ]})
        self.assertIsNone(TOOL.recorded_channel_range(self.products))
        self.assertIsNone(TOOL.recorded_channel_range(self.flow.root))

    def test_s_band_is_imaged_only_without_the_primary_beam_asked_for(self):
        for path, band in (
            (self.products / "inspection" / "inspection_report.json",
             {"name": "MEERKAT_S2"}),
            (self.products / "calibration" / "calibration_report.json",
             "MEERKAT_S2"),
        ):
            record = json.loads(path.read_text())
            record["band"] = band
            write_json(path, record)
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.run_tool()

        self.run_tool("--no-primary-beam")

        self.assertNotIn(("G1", "primary_beam"), self.flow.calls)
        field = self.flow.commands[("G1", "field")]
        self.assertEqual(option(field, "--imsize"), "8192")
        self.assertEqual(option(field, "--cell"), "0.6")
        entry = json.loads(
            (self.flow.work / TOOL.REPORT_NAME).read_text())["targets"]["G1"]
        self.assertIsNone(entry["corrected"])
        self.assertIn("not corrected for the primary beam",
                      TOOL.summary_line("G1", entry))

    def test_a_band_the_products_contradict_is_refused(self):
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            self.run_tool("--band", "MEERKAT_UHF")
        self.assertEqual(self.flow.calls, [])


class UgmrtBandTests(unittest.TestCase):
    """A uGMRT field takes its band, grid and channels from the products.

    Built as Band5_Products_1103 recorded them: 4096 channels of 97.7 kHz,
    3% of each end masked; and Band 4's, averaged by four into the targets
    after its guaranteed range was masked.
    """

    def setUp(self):
        self.flow = FlowTests("test_every_target_goes_through_the_four_steps_in_turn")
        self.flow.setUp()
        self.addCleanup(self.flow.doCleanups)
        self.products = self.flow.products

    def record(self, band, edges, guaranteed=None, width=None):
        products = self.products
        write_json(products / "inspection" / "inspection_report.json", {
            "telescope": {"name": "ugmrt"}, "band": {"name": band},
            "correlations": {"setups": {"0": ["RR", "LL"]},
                             "basis": "circular"},
            "spws": {"0": {"id": 0, "nchan": 4096}},
        })
        write_json(products / "calibration" / "calibration_report.json",
                   {"telescope": "ugmrt", "band": band})
        write_json(products / "flagging" / "edge_flagging_report.json",
                   {"telescope": "ugmrt", "band": band, "spw_string": edges})
        if guaranteed is not None:
            write_json(
                products / "flagging" / "guaranteed_range_mask_report.json",
                {"spws": [{"spw": "0", "total_channels": 4096,
                           "channel_ranges": guaranteed}]})
        if width is not None:
            write_json(
                products / "calibration" / "target_channel_average.json",
                {"schema_version": 1, "width": width})

    def run_tool(self, *extra):
        with mock.patch("builtins.print"):
            return TOOL.main(["--products", str(self.products),
                              "--workdir", str(self.flow.work),
                              "--targets", "G1", *extra])

    def test_band_5_is_made_as_its_run_was(self):
        self.record("B5", "0:0~121;3974~4095")
        self.run_tool("--average-channels", "4")

        field = self.flow.commands[("G1", "field")]
        self.assertEqual(option(field, "--imsize"), "6144")
        self.assertEqual(option(field, "--cell"), "0.5")
        # 122-3974 by hand then, averaged by four: 31-993.
        self.assertEqual(option(field, "--channel-range", 2), ["31", "993"])
        self.assertEqual(
            option(self.flow.commands[("G1", "scans")], "--channel-range", 2),
            ["122", "3974"])
        beam = self.flow.commands[("G1", "primary_beam")]
        self.assertEqual(option(beam, "--band"), "B5")
        self.assertNotIn("--telescope", beam)
        # The recipe keeps its own uGMRT reference antennas.
        self.assertNotIn("--refant", field)

    def test_band_4_counts_its_targets_averaged_channels(self):
        self.record("B4", "0:0~204;3891~4095", guaranteed=[[3072, 4095]],
                    width=4)
        self.run_tool()

        field = self.flow.commands[("G1", "field")]
        self.assertEqual(option(field, "--imsize"), "5120")
        self.assertEqual(option(field, "--cell"), "1")
        self.assertEqual(option(field, "--channel-range", 2), ["52", "768"])

    def test_the_low_bands_get_fields_for_their_beams(self):
        for band, imsize, cell in (("B2", "8192", "3"), ("B3", "8192", "1.5")):
            with self.subTest(band=band):
                self.flow.commands.clear()
                self.record(band, "0:0~286;3809~4095")
                shutil_rmtree(self.flow.work)
                self.run_tool()
                field = self.flow.commands[("G1", "field")]
                self.assertEqual(option(field, "--imsize"), imsize)
                self.assertEqual(option(field, "--cell"), cell)

    def test_a_given_channel_range_still_wins(self):
        self.record("B5", "0:0~121;3974~4095")
        self.run_tool("--channel-range", "200", "3800")

        self.assertEqual(
            option(self.flow.commands[("G1", "field")], "--channel-range", 2),
            ["200", "3800"])


def shutil_rmtree(path):
    import shutil

    shutil.rmtree(path, ignore_errors=True)

    def test_a_field_without_a_noise_image_says_nothing_of_one(self):
        line = TOOL.summary_line("J2344", {"result": {
            "image_rms_Jy_beam": 9.1e-6, "quietest_tenth_Jy_beam": 8.2e-6,
            "data_noise_Jy_beam": None,
        }, "sources_out": 2, "corrected": "x-MFS-image-pb.fits"})

        self.assertIn("J2344: 9.1 uJy/beam (8.2 in the quietest tenth); "
                      "2 sources taken out", line)
        self.assertNotIn("for the data", line)


if __name__ == "__main__":
    unittest.main()
