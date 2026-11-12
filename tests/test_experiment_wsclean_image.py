import importlib.util
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_tool():
    spec = importlib.util.spec_from_file_location(
        "test_experiment_wsclean_image_module",
        SOURCE_DIR / "experiment_wsclean_image.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TOOL = load_tool()


def arguments(*extra):
    parser_args = []

    def capture(args):
        parser_args.append(args)

    with mock.patch.object(TOOL, "run", capture):
        TOOL.main(["--target", "G1", *extra])
    return parser_args[0]


class StartTests(unittest.TestCase):
    """Both tools start from the pipeline folder with a plain Python."""

    def test_they_start_beside_the_pipelines_own_inspect_module(self):
        import subprocess

        self.assertTrue((SOURCE_DIR / "inspect.py").is_file())
        for script in ("experiment_wsclean_image.py", "diagnose_noise_floor.py"):
            source = (SOURCE_DIR / script).read_text()
            # The standard library's inspect must not be shadowed by it.
            self.assertIn("sys.path.pop(0)", source.split("import argparse")[0])
            done = subprocess.run(
                [sys.executable, script, "--help"],
                cwd=SOURCE_DIR, capture_output=True, text=True,
            )
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertIn("--target", done.stdout)


class ChannelRangeTests(unittest.TestCase):

    @staticmethod
    def summary(flagged_from):
        return {
            f"0:{channel}": {
                "flagged": 1000.0 if channel >= flagged_from else 400.0,
                "total": 1000.0,
            }
            for channel in range(1024)
        }

    def test_a_flagged_band_edge_is_left_out(self):
        # uGMRT Band 4 averaged to 1024 channels: the top quarter of the
        # recorded band is beyond the receiver's usable range.
        self.assertEqual(
            TOOL.unflagged_channel_range(self.summary(768)), (0, 768)
        )

    def test_flagged_channels_inside_the_band_do_not_split_it(self):
        summary = self.summary(768)
        for channel in range(300, 340):
            summary[f"0:{channel}"]["flagged"] = 1000.0

        self.assertEqual(TOOL.unflagged_channel_range(summary), (0, 768))

    def test_a_band_flagged_at_its_start_begins_where_the_data_do(self):
        summary = {
            f"0:{channel}": {
                "flagged": 1000.0 if channel < 256 else 300.0,
                "total": 1000.0,
            }
            for channel in range(1024)
        }

        self.assertEqual(TOOL.unflagged_channel_range(summary), (256, 1024))

    def test_other_spectral_windows_and_bad_entries_are_ignored(self):
        summary = self.summary(768)
        summary["1:900"] = {"flagged": 0.0, "total": 1000.0}
        summary["nonsense"] = {"flagged": 0.0, "total": 1000.0}
        summary["0:5"] = {"flagged": "x"}

        self.assertEqual(TOOL.unflagged_channel_range(summary), (0, 768))

    def test_nothing_usable_gives_no_range(self):
        self.assertIsNone(TOOL.unflagged_channel_range(self.summary(0)))
        self.assertIsNone(TOOL.unflagged_channel_range(None))


class CommandTests(unittest.TestCase):

    def command(self, stokes, *extra, channel_range=(0, 768)):
        return TOOL.wsclean_command(
            arguments(*extra), "final_selfcal.ms", "work/img", 5120, 1.0,
            channel_range, stokes, "work/tmp",
        )

    @staticmethod
    def value(command, option, count=1):
        index = command.index(option)
        values = command[index + 1:index + 1 + count]
        return values[0] if count == 1 else values

    def test_the_stokes_i_command_is_the_field_the_pipeline_images(self):
        command = self.command("I")

        self.assertEqual(command[0], "wsclean")
        self.assertEqual(command[-1], "final_selfcal.ms")
        self.assertEqual(self.value(command, "-size", 2), ["5120", "5120"])
        self.assertEqual(self.value(command, "-scale"), "1.0000asec")
        self.assertEqual(self.value(command, "-weight", 2), ["briggs", "0"])
        self.assertEqual(self.value(command, "-pol"), "I")
        self.assertEqual(self.value(command, "-gridder"), "wgridder")
        # Sub-bands deconvolved together, over the channels that hold data.
        self.assertIn("-join-channels", command)
        self.assertEqual(self.value(command, "-channels-out"), "6")
        self.assertEqual(self.value(command, "-fit-spectral-pol"), "3")
        self.assertEqual(
            self.value(command, "-channel-range", 2), ["0", "768"]
        )
        self.assertEqual(self.value(command, "-auto-mask"), "4")
        self.assertEqual(self.value(command, "-auto-threshold"), "1.5")
        self.assertIn("-local-rms", command)
        self.assertEqual(self.value(command, "-multiscale-scales"), "0,5,15")

    def test_the_measurement_set_is_never_written_to(self):
        for stokes in ("I", "V"):
            command = self.command(stokes)
            self.assertIn("-no-update-model-required", command)
            self.assertEqual(self.value(command, "-data-column"), "DATA")
            self.assertEqual(self.value(command, "-temp-dir"), "work/tmp")

    def test_wsclean_is_kept_inside_the_session_memory(self):
        command = self.command("I", "--memory-gb", "12", "--threads", "4")

        self.assertEqual(self.value(command, "-abs-mem"), "12")
        self.assertEqual(self.value(command, "-j"), "4")

    def test_the_stokes_v_image_is_dirty_with_the_same_weighting(self):
        noise, sky = self.command("V"), self.command("I")

        self.assertEqual(self.value(noise, "-niter"), "0")
        self.assertNotIn("-auto-mask", noise)
        self.assertNotIn("-multiscale", noise)
        for option, count in (
            ("-weight", 2), ("-size", 2), ("-scale", 1),
            ("-channels-out", 1), ("-channel-range", 2), ("-gridder", 1),
        ):
            self.assertEqual(
                self.value(noise, option, count),
                self.value(sky, option, count),
            )
        self.assertIn("-join-channels", noise)

    def test_a_container_command_is_split_into_its_words(self):
        # On ilifu WSClean is reached through a container, not on PATH.
        wrapper = "singularity exec /software/astro/containers/oxkat-latest.sif wsclean"

        command = self.command("I", "--wsclean", wrapper)

        self.assertEqual(
            command[:5],
            ["singularity", "exec",
             "/software/astro/containers/oxkat-latest.sif", "wsclean",
             "-name"],
        )
        self.assertEqual(command[-1], "final_selfcal.ms")

    def test_an_older_wsclean_is_given_the_options_it_has(self):
        old = {
            "version": (3, 1),
            "options": {"-use-wgridder", "-local-rms", "-auto-mask"}
            | {f"-o{n}" for n in range(30)},
        }

        command = TOOL.wsclean_command(
            arguments(), "x.ms", "work/img", 5120, 1.0, (0, 768), "I",
            "work/tmp", old,
        )

        self.assertIn("-use-wgridder", command)
        self.assertNotIn("-gridder", command)
        self.assertIn("-local-rms", command)
        # Not in this release's help, so not passed.
        self.assertNotIn("-parallel-deconvolution", command)

    def test_parallel_deconvolution_is_left_off_where_it_crashes(self):
        # WSClean 3.1.0 on ilifu: "double free or corruption" as
        # deconvolution began, with the image split for parallel cleaning.
        def command(version, *extra):
            return TOOL.wsclean_command(
                arguments(*extra), "x.ms", "work/img", 5120, 1.0, (51, 768),
                "I", "work/tmp", {"version": version, "options": set()},
            )

        self.assertNotIn("-parallel-deconvolution", command((3, 1)))
        # What it is meant to keep: the mask on local noise, and scales.
        self.assertIn("-local-rms", command((3, 1)))
        self.assertIn("-multiscale", command((3, 1)))
        self.assertEqual(
            self.value(command((3, 6)), "-parallel-deconvolution"), "2048"
        )
        self.assertEqual(
            self.value(
                command((3, 1), "--parallel-deconvolution", "1024"),
                "-parallel-deconvolution",
            ),
            "1024",
        )
        self.assertNotIn(
            "-parallel-deconvolution",
            command((3, 6), "--parallel-deconvolution", "0"),
        )

    def test_a_positive_model_can_be_written_for_calibration(self):
        sky = self.command("I", "--no-negative", "--write-model")
        noise = self.command("V", "--no-negative", "--write-model")

        self.assertIn("-no-negative", sky)
        self.assertNotIn("-no-update-model-required", sky)
        # The mask, the depth and the scales stay those of the image.
        self.assertIn("-auto-mask", sky)
        self.assertIn("-multiscale", sky)
        # The noise image writes nothing whatever is asked.
        self.assertIn("-no-update-model-required", noise)
        self.assertNotIn("-no-negative", noise)

    def test_a_negative_robust_and_no_scales_are_passed_through(self):
        command = self.command(
            "I", "--robust", "-0.3", "--scales", channel_range=None
        )

        self.assertEqual(self.value(command, "-weight", 2), ["briggs", "-0.3"])
        self.assertNotIn("-multiscale", command)
        self.assertNotIn("-channel-range", command)


class CapabilityTests(unittest.TestCase):

    HELP_3_6 = "\n".join(
        ["Syntax: wsclean [options] <input-ms>", "", "  ** GENERAL **"]
        + [f"-{name} <x>\n   Description mentioning -use-wgridder."
           for name in (
               "j", "abs-mem", "name", "size", "scale", "weight", "pol",
               "gridder", "niter", "nmiter", "mgain", "auto-mask",
               "auto-threshold", "local-rms", "multiscale",
               "multiscale-scales", "join-channels", "channels-out",
               "channel-range", "fit-spectral-pol", "data-column",
               "temp-dir", "no-update-model-required",
               "parallel-deconvolution",
           )]
    )

    def capabilities(self, version_text, help_text):
        def run(command, **kwargs):
            text = version_text if command[-1] == "--version" else help_text
            return mock.Mock(stdout=text)

        with mock.patch.object(TOOL.subprocess, "run", run):
            return TOOL.wsclean_capabilities(["wsclean"])

    def test_version_and_options_are_read_from_wsclean_itself(self):
        found = self.capabilities(
            "INFO:    Mounting image with FUSE.\n\n"
            "WSClean version 3.6 (2025-02-07)\n", self.HELP_3_6,
        )

        self.assertEqual(found["version"], (3, 6))
        self.assertIn("-gridder", found["options"])
        self.assertIn("-local-rms", found["options"])
        # Named in a description, not offered as an option.
        self.assertNotIn("-use-wgridder", found["options"])
        self.assertEqual(
            TOOL.gridder_arguments(found), ["-gridder", "wgridder"]
        )

    def test_a_release_without_gridder_uses_the_older_switch(self):
        old_help = self.HELP_3_6.replace("-gridder <x>", "-use-wgridder")
        found = self.capabilities(
            "WSClean version 3.1.0 (2022-04-01)\n", old_help
        )

        self.assertEqual(found["version"], (3, 1))
        self.assertEqual(TOOL.gridder_arguments(found), ["-use-wgridder"])

    def test_without_a_help_text_the_version_decides(self):
        self.assertEqual(
            TOOL.gridder_arguments({"version": (3, 1), "options": set()}),
            ["-use-wgridder"],
        )
        self.assertEqual(
            TOOL.gridder_arguments({"version": (3, 6), "options": set()}),
            ["-gridder", "wgridder"],
        )
        self.assertEqual(
            TOOL.gridder_arguments(None), ["-gridder", "wgridder"]
        )

    def test_a_command_that_is_not_wsclean_reports_nothing(self):
        found = self.capabilities("bash: wsclean: command not found\n", "")

        self.assertEqual(found, {"version": None, "options": set()})

    def test_a_command_that_cannot_start_reports_nothing(self):
        with mock.patch.object(
            TOOL.subprocess, "run", side_effect=OSError("no such file")
        ):
            found = TOOL.wsclean_capabilities(["nowhere"])

        self.assertEqual(found, {"version": None, "options": set()})


class LadderTests(unittest.TestCase):

    def test_each_rung_drops_one_feature(self):
        ladder = TOOL.feature_ladder(arguments())

        self.assertEqual(
            [TOOL.describe(rung) for rung in ladder],
            [
                "multi-scale, local-rms mask, spectral fit",
                "local-rms mask, spectral fit",
                "spectral fit",
                "plain joined-channel clean",
            ],
        )

    def test_features_already_off_are_not_rungs(self):
        ladder = TOOL.feature_ladder(
            arguments("--scales", "--spectral-terms", "0")
        )

        self.assertEqual(
            [TOOL.describe(rung) for rung in ladder],
            ["local-rms mask", "plain joined-channel clean"],
        )


class CalmerTests(unittest.TestCase):

    def test_the_largest_scale_goes_first(self):
        rung = TOOL.feature_ladder(arguments("--scales", "0", "5", "15", "30"))[0]

        gentler = TOOL.calmer(rung, [0, 5, 15, 30])

        self.assertEqual(gentler["scales"], [0, 5, 15])
        self.assertTrue(gentler["multiscale"])
        self.assertEqual(
            TOOL.describe(gentler),
            "multi-scale to 15 px, local-rms mask, spectral fit",
        )
        # The asked-for rung is left as it was.
        self.assertNotIn("scales", rung)

    def test_with_fewer_than_three_scales_left_multiscale_goes(self):
        rung = dict(TOOL.feature_ladder(arguments())[0], scales=[0, 5, 15])

        gentler = TOOL.calmer(rung, [0, 5, 15, 30])

        self.assertFalse(gentler["multiscale"])
        self.assertNotIn("scales", gentler)
        self.assertEqual(gentler, TOOL.feature_ladder(arguments())[1])

    def test_without_multiscale_there_is_nothing_calmer(self):
        rung = TOOL.feature_ladder(arguments("--scales"))[0]

        self.assertIsNone(TOOL.calmer(rung, []))

    def test_the_command_takes_the_rung_s_scales(self):
        args = arguments("--scales", "0", "5", "15", "30")
        rung = TOOL.calmer(TOOL.feature_ladder(args)[0], args.scales)

        command = TOOL.wsclean_command(
            args, "x.ms", Path("w/n"), 200, 1.0, (51, 768), "I", Path("t"),
            None, rung,
        )

        self.assertEqual(
            command[command.index("-multiscale-scales") + 1], "0,5,15"
        )


class DivergenceWatchTests(unittest.TestCase):

    @staticmethod
    def line(value):
        return f"Estimated standard deviation of background noise: {value}\n"

    def test_a_converging_clean_is_left_alone(self):
        # G0218-0245's robust-0.5 image, which converged at 2.1 times the
        # noise of the data: the residual's noise falls at every iteration.
        watch = TOOL.DivergenceWatch()

        for value in ("88.54 \u00b5Jy", "66.12 \u00b5Jy", "44.82 \u00b5Jy",
                      "36.12 \u00b5Jy", "34.85 \u00b5Jy"):
            self.assertIsNone(watch(self.line(value)))
        self.assertIsNone(watch("Iteration 202, scale 0 px : 2.38 mJy\n"))

    def test_a_diverging_clean_is_stopped_at_twice_its_lowest_noise(self):
        watch = TOOL.DivergenceWatch()
        said = [
            watch(self.line(value))
            for value in ("118.38 \u00b5Jy", "72.13 \u00b5Jy", "58.62 \u00b5Jy",
                          "62.87 \u00b5Jy", "109.38 \u00b5Jy", "257.56 \u00b5Jy")
        ]

        self.assertEqual(said[:5], [None] * 5)
        self.assertIn("from 58.6 to 257.6 uJy/beam", said[5])

    def test_every_unit_wsclean_gives_is_read(self):
        watch = TOOL.DivergenceWatch()

        self.assertIsNone(watch(self.line("40 uJy")))
        self.assertIsNone(watch(self.line("41 \u03bcJy")))
        self.assertIsNotNone(watch(self.line("1.92 mJy")))
        watch = TOOL.DivergenceWatch()
        self.assertIsNone(watch(self.line("2 Jy")))
        self.assertIsNone(watch(self.line("900000000 nJy")))
        self.assertIsNotNone(watch(self.line("4.1e+00 Jy")))


class ProductTests(unittest.TestCase):

    def test_the_joined_image_is_found_under_either_naming(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            for name in (
                "a-0000-image.fits", "a-MFS-image.fits", "a-MFS-residual.fits",
                "b-MFS-I-image.fits", "c-image.fits",
            ):
                (root / name).touch()

            self.assertEqual(
                TOOL.product(root / "a", "image").name, "a-MFS-image.fits"
            )
            self.assertEqual(
                TOOL.product(root / "a", "residual").name,
                "a-MFS-residual.fits",
            )
            self.assertEqual(
                TOOL.product(root / "b", "image").name, "b-MFS-I-image.fits"
            )
            self.assertEqual(
                TOOL.product(root / "c", "image").name, "c-image.fits"
            )
            self.assertIsNone(TOOL.product(root / "a", "model"))


class RunCommandTests(unittest.TestCase):

    def test_output_is_shown_and_kept_and_failure_is_raised(self):
        with tempfile.TemporaryDirectory() as root:
            log = Path(root) / "run.log"
            with mock.patch.object(TOOL.sys, "stdout") as shown:
                TOOL.run_command(
                    [sys.executable, "-c", "print('major iteration 1')"], log
                )
            self.assertIn("major iteration 1", log.read_text())
            self.assertIn(sys.executable, log.read_text().splitlines()[0])
            shown.write.assert_any_call("major iteration 1\n")

            with mock.patch.object(TOOL.sys, "stdout"):
                with self.assertRaises(RuntimeError) as raised:
                    TOOL.run_command(
                        [sys.executable, "-c", "raise SystemExit(3)"], log
                    )
            self.assertIn("status 3", str(raised.exception))

    def test_a_watched_command_is_stopped_when_the_watch_says_so(self):
        script = (
            "import sys, time\n"
            "for value in (60, 50, 120):\n"
            "    print('Estimated standard deviation of background noise: '"
            " + str(value) + ' uJy', flush=True)\n"
            "time.sleep(60)\n"
        )
        with tempfile.TemporaryDirectory() as root:
            log = Path(root) / "run.log"
            with mock.patch.object(TOOL.sys, "stdout"):
                with self.assertRaises(TOOL.Diverged) as raised:
                    TOOL.run_command(
                        [sys.executable, "-c", script], log,
                        watch=TOOL.DivergenceWatch(),
                    )
            text = log.read_text()
        self.assertIn("from 50.0 to 120.0 uJy/beam", str(raised.exception))
        self.assertIn("Stopped: the residual's noise rose", text)
        self.assertFalse(raised.exception.crashed)


def write_fits(path, plane):
    """A minimal single-HDU FITS image from an array indexed [x, y]."""

    nx, ny = plane.shape
    cards = [
        "SIMPLE  =                    T",
        "BITPIX  =                  -32",
        "NAXIS   =                    4",
        f"NAXIS1  = {nx:20d}",
        f"NAXIS2  = {ny:20d}",
        "NAXIS3  =                    1",
        "NAXIS4  =                    1",
        "BSCALE  =                   1. / scale",
        "BZERO   =                   0.",
        "BUNIT   = 'JY/BEAM '",
        "CDELT1  =  -2.777777777778E-04",
        "END",
    ]
    header = "".join(card.ljust(80) for card in cards)
    header = header.ljust(2880 * (len(header) // 2880 + 1))
    data = np.asarray(plane.T, dtype=">f4").tobytes()
    data += b"\0" * (-len(data) % 2880)
    Path(path).write_bytes(header.encode("ascii") + data)


class FitsTests(unittest.TestCase):

    def test_a_wsclean_image_is_read_with_numpy_alone(self):
        plane = np.arange(12, dtype=float).reshape(4, 3) * 1e-5
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "x-MFS-image.fits"
            write_fits(path, plane)
            with mock.patch.dict(
                sys.modules, {"casatools": None, "astropy": None}
            ):
                read = TOOL.fits_plane(path)

        self.assertEqual(read.shape, (4, 3))
        # Indexed [x, y], as the CASA tools return an image.
        np.testing.assert_allclose(read, plane, rtol=1e-6)

    def test_a_truncated_file_is_refused(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "bad.fits"
            path.write_bytes(b"SIMPLE  =                    T".ljust(100))

            with self.assertRaises(ValueError):
                TOOL.fits_plane(path)


def field(random, scales, size=400, sigma=150e-6):
    """Interval images of the same sources with a gain factor each."""

    sources = [((60, 80), 0.038), ((200, 200), 0.023), ((320, 110), 0.020),
               ((110, 330), 0.015), ((300, 320), 0.002)]
    planes = []
    for scale in scales:
        plane = random.normal(0.0, sigma, (size, size))
        for (x, y), flux in sources:
            plane[x * size // 400, y * size // 400] += flux * scale
        planes.append(plane)
    return planes


class TimeSplitTests(unittest.TestCase):

    def setUp(self):
        self.random = np.random.default_rng(21)

    def test_the_brightest_sources_are_found_apart_from_each_other(self):
        (plane,) = field(self.random, [1.0])
        plane[61, 81] = 0.030       # the same source, not another one

        found = TOOL.brightest_sources(plane, 4)

        self.assertEqual(found, [(60, 80), (200, 200), (320, 110), (110, 330)])

    def test_a_scale_that_moves_shows_in_every_source(self):
        table = TOOL.time_split_table(
            field(self.random, [0.85, 0.95, 1.05, 1.15])
        )

        # The faint source is left out: one interval's noise is too large.
        self.assertEqual(len(table["sources"]), 4)
        medians = [entry["median_ratio"] for entry in table["intervals"]]
        for measured, expected in zip(medians, [0.85, 0.95, 1.05, 1.15]):
            self.assertAlmostEqual(measured, expected, delta=0.02)
        text = " ".join(TOOL.time_split_reading(table))
        self.assertIn("change together by 30%", text)
        self.assertIn("one gain per scan cannot follow", text)

    def test_a_steady_scale_is_reported_as_steady(self):
        table = TOOL.time_split_table(field(self.random, [1.0] * 4))

        text = " ".join(TOOL.time_split_reading(table))
        self.assertIn("amplitude scale is steady", text)
        self.assertNotIn("change together", text)

    def test_without_bright_sources_nothing_is_concluded(self):
        planes = [self.random.normal(0.0, 150e-6, (200, 200)) for _ in range(4)]

        table = TOOL.time_split_table(planes)

        self.assertEqual(table["sources"], [])
        self.assertIn(
            "No source is bright enough", TOOL.time_split_reading(table)[0]
        )

    def test_the_command_makes_dirty_images_of_consecutive_stretches(self):
        command = TOOL.time_split_command(
            arguments("--time-split", "4"), "scan15.ms", "work/ts", 5120, 1.0,
            (51, 768), "work/tmp", {"version": (3, 1), "options": set()},
        )

        index = command.index("-intervals-out")
        self.assertEqual(command[index + 1], "4")
        self.assertEqual(command[command.index("-niter") + 1], "0")
        self.assertIn("-no-update-model-required", command)
        self.assertIn("-use-wgridder", command)
        self.assertEqual(command[command.index("-pol") + 1], "I")
        self.assertNotIn("-join-channels", command)
        self.assertEqual(command[-1], "scan15.ms")


class AsymmetryTests(unittest.TestCase):

    def setUp(self):
        self.random = np.random.default_rng(55)

    def image(self, lopsided):
        """Sources with a symmetric beam, and optionally a +/- pair of lobes."""

        plane = self.random.normal(0.0, 30e-6, (900, 900))
        for (x, y), flux in (((200, 250), 0.037), ((450, 450), 0.023),
                             ((700, 300), 0.020), ((300, 700), 0.015)):
            plane[x, y] += flux
            for dx, dy in ((14, 0), (-14, 0)):          # the beam's own lobes
                plane[x + dx, y + dy] += 0.05 * flux
            if lopsided:
                plane[x + 8, y - 9] += 0.13 * flux
                plane[x - 8, y + 9] -= 0.13 * flux
        return plane

    def test_a_lopsided_pattern_common_to_the_sources_is_measured(self):
        found = TOOL.source_asymmetry(self.image(True))

        self.assertEqual(found["sources"], 4)
        self.assertAlmostEqual(found["antisymmetric_peak"], 0.13, delta=0.01)
        self.assertIn(found["antisymmetric_peak_offset"], ([8, -9], [-8, 9]))
        self.assertAlmostEqual(
            found["deepest_negative_median"], -0.13, delta=0.02
        )

    def test_a_symmetric_beam_is_not_called_lopsided(self):
        found = TOOL.source_asymmetry(self.image(False))

        self.assertLess(found["antisymmetric_peak"], 0.02)
        self.assertGreater(found["deepest_negative_median"], -0.02)

    def test_an_image_without_bright_sources_has_no_measure(self):
        plane = self.random.normal(0.0, 30e-6, (400, 400))

        self.assertIsNone(TOOL.source_asymmetry(plane))

    def test_sources_at_the_edge_are_left_out(self):
        plane = self.image(True)
        plane[10, 10] += 0.5

        found = TOOL.source_asymmetry(plane)

        self.assertNotIn([10, 10], found["pixels"])
        self.assertEqual(found["sources"], 4)


class CoreModelTests(unittest.TestCase):

    def command(self, *extra):
        return TOOL.core_model_command(
            arguments("--core-model", "6", *extra), "copy/scan15.ms",
            "work/core", 5120, 1.0, (51, 768), "work/tmp",
            {"version": (3, 1), "options": set()},
        )

    def test_only_bright_positive_emission_goes_into_the_model(self):
        command = self.command()

        self.assertEqual(command[command.index("-threshold") + 1], "0.006")
        self.assertIn("-no-negative", command)
        # No mask grown down to the noise, and no scales: cores only.
        self.assertNotIn("-auto-mask", command)
        self.assertNotIn("-auto-threshold", command)
        self.assertNotIn("-multiscale", command)
        self.assertIn("-join-channels", command)
        self.assertEqual(command[-1], "copy/scan15.ms")

    def test_the_cores_can_be_cut_at_a_fraction_of_the_brightest_pixel(self):
        args = arguments("--core-model-fraction", "0.13")

        command = TOOL.core_model_command(
            args, "copy.ms", "work/coremodel", 5120, 1.0, (51, 768),
            "work/tmp", {"version": (3, 1), "options": set()},
        )

        # One major iteration that ends at (1 - mgain) of the first peak.
        self.assertEqual(command[command.index("-nmiter") + 1], "1")
        self.assertEqual(command[command.index("-mgain") + 1], "0.87")
        self.assertNotIn("-threshold", command)
        self.assertIn("-no-negative", command)
        self.assertNotIn("-no-update-model-required", command)
        # Never into the noise, however faint the field's brightest pixel.
        self.assertEqual(command[command.index("-auto-threshold") + 1], "5")

    def test_the_model_column_is_written(self):
        # Every other mode forbids it; here it is the product.
        self.assertNotIn("-no-update-model-required", self.command())

    def test_a_path_inside_the_products_is_recognised(self):
        with tempfile.TemporaryDirectory() as root:
            products = Path(root) / "uGMRT_Products"
            inside = products / "targets" / "G1" / "final_selfcal.ms"
            outside = Path(root) / "wsclean_test" / "G1_scan15.ms"
            inside.mkdir(parents=True)
            outside.mkdir(parents=True)

            self.assertTrue(TOOL.inside(inside, products))
            self.assertFalse(TOOL.inside(outside, products))


class UvSplitTests(unittest.TestCase):

    EDGES = [0.6, 3.0, 8.0, 20.0, 35.0, 70.0]

    def setUp(self):
        self.random = np.random.default_rng(33)

    def table(self, scales):
        planes = field(self.random, scales)
        return TOOL.time_split_table(
            planes, minimum_ratio=8.0, locate=planes[-1], typical_noise=True
        )

    def test_long_baselines_brighter_than_short_ones_is_reported(self):
        # Central-square antennas low and arm antennas high: flux rises
        # with baseline length for every source.
        table = self.table([0.80, 0.93, 1.0, 1.07, 1.20])

        text = " ".join(TOOL.uv_split_reading(table, self.EDGES))

        self.assertIn("3-8 klambda 0.9", text)
        self.assertIn("35-70 klambda 1.1", text)
        self.assertIn("brighter on the longest baselines", text)
        self.assertIn("depends on baseline length", text)

    def test_sources_fainter_on_long_baselines_are_resolved_not_miscalibrated(self):
        # G0212-0335, scan 15: the phase-centre source fell from 80 to 13
        # mJy/beam between 3-8 and 35-70 klambda and its neighbours by a
        # factor of two. The first wording blamed the antennas for it.
        planes = field(self.random, [1.0] * 5)
        fall = {(60, 80): [1.7, 1.3, 0.9, 0.75, 0.42],
                (200, 200): [2.0, 1.5, 0.84, 0.46, 0.25],
                (320, 110): [1.3, 1.2, 1.0, 0.88, 0.64],
                (110, 330): [1.1, 1.5, 1.15, 0.73, 0.53]}
        flux = {(60, 80): 0.038, (200, 200): 0.023, (320, 110): 0.020,
                (110, 330): 0.015}
        for index, plane in enumerate(planes):
            for pixel, factors in fall.items():
                plane[pixel] += flux[pixel] * (factors[index] - 1.0)
        table = TOOL.time_split_table(
            planes, minimum_ratio=8.0, locate=planes[-1], typical_noise=True
        )

        text = " ".join(TOOL.uv_split_reading(table, self.EDGES))

        self.assertIn("they are resolved", text)
        self.assertIn("multi-scale clean", text)
        self.assertNotIn("antennas on different scales", text)

    def test_a_flat_scale_is_reported_as_flat(self):
        # The shortest range blends neighbours and is not judged.
        table = self.table([1.3, 1.0, 1.0, 1.0, 1.0])

        text = " ".join(TOOL.uv_split_reading(table, self.EDGES))

        self.assertIn("does not depend on baseline length", text)
        self.assertIn("compact", text)

    def test_fainter_on_long_baselines_is_said_so(self):
        table = self.table([1.0, 1.1, 1.03, 0.97, 0.9])

        text = " ".join(TOOL.uv_split_reading(table, self.EDGES))

        self.assertIn("fainter on the longest baselines", text)
        self.assertIn("they are resolved", text)

    def test_each_range_is_a_natural_weighted_dirty_image(self):
        command = TOOL.uv_split_command(
            arguments(), "scan15.ms", "work/b2", 5120, 1.0, (51, 768),
            "work/tmp", {"version": (3, 1), "options": set()}, 8.0, 20.0,
        )

        self.assertEqual(command[command.index("-minuv-l") + 1], "8000")
        self.assertEqual(command[command.index("-maxuv-l") + 1], "20000")
        self.assertEqual(command[command.index("-weight") + 1], "natural")
        self.assertEqual(command[command.index("-niter") + 1], "0")
        self.assertIn("-no-update-model-required", command)
        self.assertEqual(command[-1], "scan15.ms")


class FlowTests(unittest.TestCase):
    """The whole run against stand-ins for CASA and WSClean."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        selfcal = self.root / "out" / "targets" / "G1" / "selfcal"
        self.ms = selfcal / "selfcal_ms" / "final_selfcal.ms"
        self.reference = selfcal / "images" / "final.image.tt0"
        for path in (self.ms, self.reference):
            path.mkdir(parents=True)

        random = np.random.default_rng(11)
        self.reference_plane = random.normal(0.0, 80e-6, (200, 200))
        self.commands = []
        self.fail = None
        self.aborts_with = ()
        self.abort_status = -6
        self.diverges_with = ()
        self.flag_call = None
        self.capabilities = {"version": (3, 6), "options": set()}
        test = self

        class Image:
            def open(self, path):
                assert str(path) == str(test.reference)

            def getchunk(self):
                return test.reference_plane[:, :, None, None]

            def summary(self):
                return {"incr": [-1.0 / 206264.80624709636,
                                 1.0 / 206264.80624709636, 1.0, 1.0]}

            def close(self):
                pass

        def flagdata(**kwargs):
            test.flag_call = kwargs
            return {"spw:channel": ChannelRangeTests.summary(768)}

        def run_command(command, log_path, watch=None):
            test.commands.append(command)
            stokes = command[command.index("-pol") + 1]
            Path(log_path).write_text("log of " + " ".join(command[-3:]))
            if "-no-negative" in command and "-auto-mask" not in command:
                # The cores-only model: a threshold, no mask.
                name = command[command.index("-name") + 1]
                model = np.zeros((200, 200))
                model[60, 80], model[100, 100] = 0.030, 0.018
                write_fits(f"{name}-MFS-model.fits", model)
                dirty = random.normal(0.0, 30e-6, (200, 200))
                for x, y, flux in ((70, 80, 0.030), (120, 110, 0.018)):
                    dirty[x, y] += flux
                    dirty[x + 8, y - 9] += 0.13 * flux
                    dirty[x - 8, y + 9] -= 0.13 * flux
                write_fits(f"{name}-MFS-dirty.fits", dirty)
                return
            if "-minuv-l" in command:
                name = command[command.index("-name") + 1]
                low = float(command[command.index("-minuv-l") + 1]) / 1000.0
                (plane,) = field(random, [0.9 + 0.006 * low], size=200)
                write_fits(f"{name}-image.fits", plane)
                return
            if "-intervals-out" in command:
                name = command[command.index("-name") + 1]
                count = int(command[command.index("-intervals-out") + 1])
                planes = field(
                    random, np.linspace(0.8, 1.2, count), size=200
                )
                for index, plane in enumerate(planes):
                    write_fits(f"{name}-t{index:04d}-image.fits", plane)
                return
            if test.fail == stokes:
                raise RuntimeError("wsclean exited with status 1")
            name = Path(command[command.index("-name") + 1])
            if stokes == "I" and any(
                option in command for option in test.aborts_with
            ):
                # What an abort leaves: the first products, and no image.
                (name.parent / f"{name.name}-0000-dirty.fits").write_bytes(b"x")
                raise TOOL.CommandFailed(
                    "singularity exited with status -6", test.abort_status
                )
            if stokes == "I" and watch is not None and any(
                option in command for option in test.diverges_with
            ):
                # The residual's noise as G0213-0251's robust-0.5 clean
                # gave it, until the watch stops the run.
                (name.parent / f"{name.name}-0000-model.fits").write_bytes(b"x")
                for value in ("118.38 \u00b5Jy", "72.13 \u00b5Jy",
                              "58.62 \u00b5Jy", "62.87 \u00b5Jy",
                              "109.38 \u00b5Jy", "257.56 \u00b5Jy",
                              "769.48 \u00b5Jy", "1.92 mJy"):
                    reason = watch(
                        "Estimated standard deviation of background noise: "
                        + value)
                    if reason:
                        raise TOOL.Diverged(f"wsclean was stopped: {reason}")
                raise AssertionError("the watch did not stop the clean")
            size = int(command[command.index("-size") + 1])
            sigma = {"V": 22e-6, "Q": 22e-6, "I": 35e-6}[stokes]
            for kind in ("image", "residual"):
                write_fits(
                    name.parent / f"{name.name}-MFS-{kind}.fits",
                    random.normal(0.0, sigma, (size, size)),
                )
            Path(log_path).write_text("log")

        self.casatools = types.ModuleType("casatools")
        self.casatools.image = Image
        self.casatasks = types.ModuleType("casatasks")
        self.casatasks.flagdata = flagdata
        self.use_casa(True)
        for patch in (
            mock.patch.object(TOOL, "run_command", run_command),
            mock.patch.object(TOOL.shutil, "which", lambda name: "/bin/" + name),
            mock.patch.object(
                TOOL, "wsclean_capabilities",
                lambda prefix: test.capabilities,
            ),
        ):
            patch.start()
            self.addCleanup(patch.stop)

    def use_casa(self, available):
        """Make CASA importable, or not, as in the two ilifu terminals."""

        if getattr(self, "_casa_patch", None) is not None:
            self._casa_patch.stop()
        modules = (
            {"casatools": self.casatools, "casatasks": self.casatasks}
            if available else {"casatools": None, "casatasks": None}
        )
        self._casa_patch = mock.patch.dict(sys.modules, modules)
        self._casa_patch.start()
        self.addCleanup(self._casa_patch.stop)

    def report_path(self):
        return self.root / "work" / "G1" / "wsclean_robust_0.json"

    def run_tool(self, *extra):
        with mock.patch("builtins.print"):
            TOOL.main([
                "--outdir", str(self.root / "out"), "--target", "G1",
                "--workdir", str(self.root / "work"), *extra,
            ])
        return json.loads(self.report_path().read_text())

    def test_a_folder_another_run_holds_is_refused(self):
        import os
        import socket

        folder = self.root / "work" / "G1"
        folder.mkdir(parents=True)
        lock = folder / TOOL.LOCK_NAME
        # This test's own process stands for a run that is still going.
        lock.write_text(json.dumps({
            "pid": os.getpid(), "host": socket.gethostname(),
            "since": "2026-10-08 20:03:00 UTC",
        }))

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("in use by another run", str(stopped.exception))
        self.assertIn(str(lock), str(stopped.exception))
        self.assertEqual(self.commands, [])
        self.assertTrue(lock.is_file())

    def test_a_claim_left_by_a_run_that_is_gone_is_taken_over(self):
        import socket

        folder = self.root / "work" / "G1"
        folder.mkdir(parents=True)
        lock = folder / TOOL.LOCK_NAME
        lock.write_text(json.dumps(
            {"pid": 2 ** 22 + 12345, "host": socket.gethostname()}
        ))

        with mock.patch.object(TOOL, "process_alive", lambda pid: False):
            self.run_tool()

        self.assertTrue(self.report_path().is_file())
        self.assertFalse(lock.exists())

    def test_a_claim_from_another_machine_is_not_taken_over(self):
        folder = self.root / "work" / "G1"
        folder.mkdir(parents=True)
        lock = folder / TOOL.LOCK_NAME
        lock.write_text(json.dumps({"pid": 1, "host": "some-other-node"}))

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("some-other-node", str(stopped.exception))
        self.assertIn("remove", str(stopped.exception))

    def test_a_claim_of_a_slurm_job_that_has_ended_is_taken_over(self):
        folder = self.root / "work" / "G1"
        folder.mkdir(parents=True)
        lock = folder / TOOL.LOCK_NAME
        lock.write_text(json.dumps(
            {"pid": 1, "host": "compute-201", "slurm_job": "4242"}
        ))

        with mock.patch.object(TOOL, "slurm_job_gone",
                               lambda job: job == "4242"):
            self.run_tool()

        self.assertTrue(self.report_path().is_file())
        self.assertFalse(lock.exists())

    def test_a_claim_of_a_slurm_job_still_going_is_refused(self):
        folder = self.root / "work" / "G1"
        folder.mkdir(parents=True)
        lock = folder / TOOL.LOCK_NAME
        lock.write_text(json.dumps(
            {"pid": 1, "host": "compute-201", "slurm_job": "4242"}
        ))

        with mock.patch.object(TOOL, "slurm_job_gone", lambda job: False), \
                self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("SLURM job 4242", str(stopped.exception))
        self.assertTrue(lock.is_file())

    def test_a_claim_made_in_a_slurm_job_says_which(self):
        with mock.patch.dict("os.environ", {"SLURM_JOB_ID": "777"}):
            mine = TOOL.claim_of_this_run()

        self.assertEqual(mine["slurm_job"], "777")

    def test_squeue_says_whether_a_job_has_ended(self):
        def answer(code, out="", err=""):
            return mock.Mock(returncode=code, stdout=out, stderr=err)

        with mock.patch.object(TOOL.shutil, "which", lambda name: "/bin/x"):
            for reply, gone in (
                (answer(0, "RUNNING\n"), False),
                (answer(0, "PENDING\n"), False),
                (answer(0, ""), True),
                (answer(1, err="slurm_load_jobs error: Invalid job id "
                               "specified"), True),
                (answer(1, err="Unable to contact slurm controller"), False),
            ):
                with mock.patch.object(TOOL.subprocess, "run",
                                       lambda *a, **k: reply):
                    self.assertEqual(TOOL.slurm_job_gone("12"), gone, reply)
        with mock.patch.object(TOOL.shutil, "which", lambda name: None):
            self.assertIsNone(TOOL.slurm_job_gone("12"))
        self.assertIsNone(TOOL.slurm_job_gone(None))

    def test_the_folder_is_held_while_wsclean_runs_and_released_after(self):
        folder = self.root / "work" / "G1"
        seen = []
        real = TOOL.run_command

        def watching(command, log_path, watch=None):
            seen.append((folder / TOOL.LOCK_NAME).is_file())
            return real(command, log_path, watch)

        with mock.patch.object(TOOL, "run_command", watching):
            self.run_tool()

        self.assertTrue(seen and all(seen))
        self.assertFalse((folder / TOOL.LOCK_NAME).exists())

    def test_a_failed_run_does_not_keep_the_folder(self):
        self.fail = "I"

        with self.assertRaises(Exception):
            self.run_tool()

        self.assertFalse(
            (self.root / "work" / "G1" / TOOL.LOCK_NAME).exists()
        )

    def test_noise_first_then_the_image_and_their_ratio(self):
        report = self.run_tool()

        self.assertEqual(
            [c[c.index("-pol") + 1] for c in self.commands], ["V", "I"]
        )
        # The grid of the image it is compared with, the channels with data.
        sky = self.commands[1]
        self.assertEqual(sky[sky.index("-size") + 1], "200")
        self.assertEqual(sky[sky.index("-scale") + 1], "1.0000asec")
        self.assertEqual(sky[sky.index("-channel-range") + 2], "768")
        self.assertEqual(self.flag_call["mode"], "summary")
        self.assertEqual(self.flag_call["action"], "calculate")
        self.assertEqual(sky[-1], str(self.ms))

        self.assertEqual(report["channel_range"], [0, 768])
        self.assertAlmostEqual(
            report["stokes_v"]["robust_rms_Jy_beam"], 22e-6, delta=1e-6
        )
        self.assertAlmostEqual(
            report["stokes_i"]["robust_rms_Jy_beam"], 35e-6, delta=2e-6
        )
        self.assertAlmostEqual(
            report["reference"]["robust_rms_Jy_beam"], 80e-6, delta=4e-6
        )
        self.assertIn("stokes_i_residual", report)
        self.assertEqual(report["notes"], [])
        text = " ".join(report["reading"])
        self.assertIn("1.6 times the noise", text)
        self.assertIn("tclean image it is compared with", text)
        self.assertFalse((self.root / "work" / "G1" / "tmp").exists())

    def write_correlations(self, correlations, basis):
        path = self.root / "out" / "inspection" / "inspection_report.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"correlations": {
            "setups": {"0": correlations}, "basis": basis}}))

    def test_data_without_stokes_v_are_imaged_without_a_noise_image(self):
        # MeerKAT: XX and YY only. WSClean forms no V from them, and Q
        # carries the bright sources' off-axis leakage.
        self.write_correlations(["XX", "YY"], "linear")

        report = self.run_tool()

        self.assertEqual(
            [c[c.index("-pol") + 1] for c in self.commands], ["I"])
        self.assertIsNone(report["noise_stokes"])
        self.assertNotIn("stokes_v", report)
        self.assertNotIn("data_noise", report)
        self.assertIn("No noise image", " ".join(report["reading"]))
        self.assertAlmostEqual(
            report["stokes_i"]["robust_rms_Jy_beam"], 35e-6, delta=2e-6)

    def test_circular_hands_keep_their_stokes_v_noise_image(self):
        self.write_correlations(["RR", "LL"], "circular")

        report = self.run_tool()

        self.assertEqual(
            [c[c.index("-pol") + 1] for c in self.commands], ["V", "I"])
        self.assertEqual(report["data_noise"]["stokes"], "V")
        self.assertEqual(report["stokes_v"], report["data_noise"])
        self.assertAlmostEqual(
            report["stokes_v"]["robust_rms_Jy_beam"], 22e-6, delta=1e-6)

    def test_a_q_noise_image_can_be_asked_for_by_name(self):
        self.write_correlations(["XX", "YY"], "linear")

        report = self.run_tool("--noise-stokes", "Q")

        self.assertEqual(
            [c[c.index("-pol") + 1] for c in self.commands], ["Q", "I"])
        self.assertEqual(report["data_noise"]["stokes"], "Q")
        self.assertNotIn("stokes_v", report)
        self.assertIn("(dirty Stokes Q)", " ".join(report["reading"]))

    def test_it_runs_where_casa_is_not_installed(self):
        # The terminal that reaches the WSClean container has a Python
        # without casatools.
        self.use_casa(False)

        report = self.run_tool(
            "--imsize", "200", "--cell", "1.0", "--channel-range", "51", "768",
        )

        self.assertEqual(len(self.commands), 2)
        sky = self.commands[1]
        self.assertEqual(sky[sky.index("-size") + 1], "200")
        self.assertEqual(
            sky[sky.index("-channel-range") + 1:][:2], ["51", "768"]
        )
        self.assertAlmostEqual(
            report["stokes_i"]["robust_rms_Jy_beam"], 35e-6, delta=2e-6
        )
        self.assertNotIn("reference", report)
        self.assertIn("CASA is not importable", report["notes"][0])
        self.assertIn("1.6 times the noise", " ".join(report["reading"]))

    def test_without_casa_the_grid_has_to_be_given(self):
        self.use_casa(False)

        with self.assertRaises(SystemExit) as raised:
            self.run_tool("--channel-range", "51", "768")

        self.assertIn("--imsize 5120 --cell 1.0", str(raised.exception))
        self.assertEqual(self.commands, [])

    def test_without_casa_the_channels_have_to_be_given(self):
        self.use_casa(False)

        with self.assertRaises(SystemExit) as raised:
            self.run_tool("--imsize", "200", "--cell", "1.0")

        self.assertIn("--channel-range", str(raised.exception))
        self.assertEqual(self.commands, [])

    def test_the_comparison_is_added_later_from_a_python_with_casa(self):
        self.use_casa(False)
        self.run_tool(
            "--imsize", "200", "--cell", "1.0", "--channel-range", "51", "768",
        )
        self.commands.clear()
        self.use_casa(True)

        report = self.run_tool("--measure-only")

        self.assertEqual(self.commands, [])
        self.assertAlmostEqual(
            report["reference"]["robust_rms_Jy_beam"], 80e-6, delta=4e-6
        )
        self.assertEqual(report["notes"], [])
        self.assertEqual(report["channel_range"], [51, 768])
        self.assertIn(
            "tclean image it is compared with", " ".join(report["reading"])
        )

    def test_measuring_needs_an_earlier_run(self):
        with self.assertRaises(SystemExit) as raised:
            self.run_tool("--measure-only")

        self.assertIn("Run without --measure-only first", str(raised.exception))

    def test_without_numpy_the_images_are_still_made(self):
        with mock.patch.dict(sys.modules, {"numpy": None}):
            report = self.run_tool(
                "--imsize", "200", "--cell", "1.0",
                "--channel-range", "51", "768",
            )

        self.assertEqual(len(self.commands), 2)
        self.assertNotIn("stokes_i", report)
        self.assertIn("--measure-only", report["notes"][0])
        self.assertTrue(
            (self.root / "work" / "G1" / "wsclean_robust_0_I-MFS-image.fits")
            .is_file()
        )

    def sky_commands(self):
        return [c for c in self.commands if c[c.index("-pol") + 1] == "I"]

    def test_an_abort_is_retried_with_a_simpler_deconvolution(self):
        # WSClean 3.1.0 on ilifu: "munmap_chunk(): invalid pointer" as
        # deconvolution started, twice.
        self.aborts_with = ("-multiscale",)

        report = self.run_tool()

        first, second = self.sky_commands()
        self.assertIn("-multiscale", first)
        self.assertNotIn("-multiscale", second)
        # Only what was dropped differs.
        self.assertIn("-local-rms", second)
        self.assertIn("-fit-spectral-pol", second)
        self.assertEqual(
            [a["outcome"] for a in report["stokes_i_attempts"]],
            ["failed: singularity exited with status -6", "completed"],
        )
        self.assertEqual(
            report["stokes_i_attempts"][1]["deconvolution"],
            "local-rms mask, spectral fit",
        )
        self.assertEqual(report["commands"]["I"], " ".join(second))
        self.assertAlmostEqual(
            report["stokes_i"]["robust_rms_Jy_beam"], 35e-6, delta=2e-6
        )
        self.assertIn("aborted 1 time(s)", " ".join(report["reading"]))
        work = self.root / "work" / "G1"
        # The crash's log is kept; what it left of the image is not.
        self.assertTrue((work / "wsclean_robust_0_I.crash1.log").is_file())
        self.assertFalse((work / "wsclean_robust_0_I-0000-dirty.fits").exists())

    def test_the_ladder_drops_one_thing_at_a_time(self):
        self.aborts_with = ("-multiscale", "-local-rms")

        report = self.run_tool()

        commands = self.sky_commands()
        self.assertEqual(len(commands), 3)
        self.assertIn("-fit-spectral-pol", commands[2])
        self.assertNotIn("-local-rms", commands[2])
        self.assertEqual(
            [a["deconvolution"] for a in report["stokes_i_attempts"]],
            ["multi-scale, local-rms mask, spectral fit",
             "local-rms mask, spectral fit", "spectral fit"],
        )

    def test_when_every_rung_aborts_the_failure_is_raised(self):
        self.aborts_with = ("-join-channels",)

        with self.assertRaises(TOOL.CommandFailed):
            self.run_tool()

        commands = self.sky_commands()
        self.assertEqual(len(commands), 4)
        self.assertNotIn("-fit-spectral-pol", commands[3])
        report = json.loads(self.report_path().read_text())
        self.assertEqual(len(report["stokes_i_attempts"]), 4)
        self.assertIn("robust_rms_Jy_beam", report["stokes_v"])

    def test_a_refused_command_is_not_retried(self):
        # An option WSClean rejects exits with a status, not a signal.
        self.aborts_with = ("-multiscale",)
        self.abort_status = 255

        with self.assertRaises(TOOL.CommandFailed):
            self.run_tool()

        self.assertEqual(len(self.sky_commands()), 1)

    def test_a_signal_reported_through_a_wrapper_counts_as_an_abort(self):
        self.aborts_with = ("-multiscale",)
        self.abort_status = 134

        self.run_tool()

        self.assertEqual(len(self.sky_commands()), 2)

    def test_a_diverging_clean_is_redone_without_its_largest_scale(self):
        # G0213-0251, robust 0.5, Band 4 (2026-10-29): 30-pixel components
        # of tens of Jy, an image 2200 times the noise of the data.
        self.diverges_with = ("0,5,15,30",)

        report = self.run_tool("--scales", "0", "5", "15", "30")

        first, second = self.sky_commands()
        scales = second[second.index("-multiscale-scales") + 1]
        self.assertEqual(scales, "0,5,15")
        # Only the scales differ.
        self.assertEqual(
            [part for part in first if part != "0,5,15,30"],
            [part for part in second if part != "0,5,15"],
        )
        outcomes = [a["outcome"] for a in report["stokes_i_attempts"]]
        self.assertTrue(outcomes[0].startswith("diverged: "), outcomes)
        self.assertIn("58.6 to 257.6 uJy/beam", outcomes[0])
        self.assertEqual(outcomes[1], "completed")
        self.assertEqual(
            report["stokes_i_attempts"][1]["deconvolution"],
            "multi-scale to 15 px, local-rms mask, spectral fit",
        )
        self.assertEqual(report["commands"]["I"], " ".join(second))
        reading = " ".join(report["reading"])
        self.assertIn("diverged 1 time(s)", reading)
        self.assertNotIn("aborted", reading)
        work = self.root / "work" / "G1"
        # The log of what diverged is kept; what it left of the image is not.
        self.assertTrue((work / "wsclean_robust_0_I.diverged1.log").is_file())
        self.assertFalse((work / "wsclean_robust_0_I-0000-model.fits").exists())

    def test_with_too_few_scales_left_multiscale_goes(self):
        self.diverges_with = ("-multiscale",)

        report = self.run_tool()

        commands = self.sky_commands()
        self.assertEqual(len(commands), 2)
        self.assertNotIn("-multiscale", commands[1])
        self.assertIn("-local-rms", commands[1])
        self.assertEqual(
            report["stokes_i_attempts"][1]["deconvolution"],
            "local-rms mask, spectral fit",
        )

    def test_a_point_clean_after_divergence_still_has_the_crash_ladder(self):
        # The point clean that follows a divergence is the ladder's next
        # rung already: it is not tried twice when it aborts.
        self.diverges_with = ("-multiscale",)
        self.aborts_with = ("-local-rms",)

        report = self.run_tool()

        commands = self.sky_commands()
        self.assertEqual(len(commands), 3)
        self.assertEqual(
            [a["deconvolution"] for a in report["stokes_i_attempts"]],
            ["multi-scale, local-rms mask, spectral fit",
             "local-rms mask, spectral fit", "spectral fit"],
        )

    def test_a_clean_that_diverges_with_nothing_left_to_drop_fails(self):
        self.diverges_with = ("-join-channels",)

        with self.assertRaises(TOOL.Diverged):
            self.run_tool("--scales", "--no-local-rms", "--spectral-terms", "0")

        report = json.loads(self.report_path().read_text())
        self.assertEqual(len(report["stokes_i_attempts"]), 1)
        self.assertNotIn("stokes_i", report)

    def test_the_first_rung_is_what_was_asked_for(self):
        report = self.run_tool("--scales", "--no-local-rms")

        (command,) = self.sky_commands()
        self.assertNotIn("-multiscale", command)
        self.assertNotIn("-local-rms", command)
        self.assertEqual(
            report["stokes_i_attempts"][0]["deconvolution"], "spectral fit"
        )

    def test_the_time_split_compares_stretches_of_the_scan(self):
        self.use_casa(False)

        with mock.patch("builtins.print"):
            TOOL.main([
                "--outdir", str(self.root / "out"), "--target", "G1",
                "--workdir", str(self.root / "work"), "--imsize", "200",
                "--cell", "1.0", "--channel-range", "51", "768",
                "--time-split", "4",
            ])

        (command,) = self.commands
        self.assertIn("-intervals-out", command)
        work = self.root / "work" / "G1"
        report = json.loads((work / "timesplit_robust_0.json").read_text())
        self.assertEqual(report["interval_count"], 4)
        self.assertEqual(len(report["images"]), 4)
        medians = [entry["median_ratio"] for entry in report["intervals"]]
        self.assertLess(medians[0], 0.9)
        self.assertGreater(medians[3], 1.1)
        self.assertIn("change together", " ".join(report["reading"]))
        # The image run's own report and products are not touched.
        self.assertFalse((work / "wsclean_robust_0.json").exists())

    def test_the_core_model_is_written_into_a_copy_of_the_data(self):
        self.use_casa(False)
        copy = self.root / "copies" / "G1_scan15.ms"
        copy.mkdir(parents=True)

        with mock.patch("builtins.print"):
            TOOL.main([
                "--outdir", str(self.root / "out"), "--target", "G1",
                "--ms", str(copy), "--workdir", str(self.root / "work"),
                "--imsize", "200", "--cell", "1.0",
                "--channel-range", "51", "768", "--core-model", "6",
            ])

        (command,) = self.commands
        self.assertEqual(command[-1], str(copy))
        self.assertNotIn("-no-update-model-required", command)
        report = json.loads(
            (self.root / "work" / "G1" / "coremodel_robust_0.json").read_text()
        )
        self.assertAlmostEqual(report["model"]["total_flux_Jy"], 0.048)
        self.assertEqual(report["model"]["pixels_with_flux"], 2)
        # The state of the calibration, read off the uncleaned image.
        self.assertAlmostEqual(
            report["dirty_asymmetry"]["antisymmetric_peak"], 0.13, delta=0.01
        )
        self.assertIn("Before any cleaning", report["reading"][0])
        self.assertIn("13.0% of the peak", report["reading"][0])
        self.assertIn("none negative", report["reading"][1])
        self.assertIn("experiment_phase_offsets.py", report["reading"][2])

    def test_writing_a_model_refuses_a_pipeline_product(self):
        with self.assertRaises(SystemExit) as raised:
            self.run_tool("--no-negative", "--write-model")

        self.assertIn("give it a copy", str(raised.exception))
        self.assertEqual(self.commands, [])

    def test_writing_a_model_into_a_copy_runs_the_ordinary_image(self):
        copy = self.root / "copies" / "G1_scan15.ms"
        copy.mkdir(parents=True)

        report = self.run_tool(
            "--ms", str(copy), "--no-negative", "--write-model"
        )

        noise, sky = self.commands
        self.assertIn("-no-update-model-required", noise)
        self.assertNotIn("-no-update-model-required", sky)
        self.assertIn("-no-negative", sky)
        self.assertIn("robust_rms_Jy_beam", report["stokes_i"])

    def test_the_core_model_refuses_a_pipeline_product(self):
        # The default Measurement Set is the pipeline's own.
        with self.assertRaises(SystemExit) as raised:
            self.run_tool("--core-model", "6")

        self.assertIn("give it a copy", str(raised.exception))
        self.assertEqual(self.commands, [])

    def test_the_uv_split_compares_ranges_of_baseline_length(self):
        self.use_casa(False)

        with mock.patch("builtins.print"):
            TOOL.main([
                "--outdir", str(self.root / "out"), "--target", "G1",
                "--workdir", str(self.root / "work"), "--imsize", "200",
                "--cell", "1.0", "--channel-range", "51", "768",
                "--uv-split", "0.6", "3", "8", "20", "35", "70",
            ])

        self.assertEqual(len(self.commands), 5)
        self.assertEqual(
            [c[c.index("-minuv-l") + 1] for c in self.commands],
            ["600", "3000", "8000", "20000", "35000"],
        )
        work = self.root / "work" / "G1"
        report = json.loads((work / "uvsplit_robust_0.json").read_text())
        self.assertEqual(report["edges_klambda"], [0.6, 3, 8, 20, 35, 70])
        medians = [entry["median_ratio"] for entry in report["intervals"]]
        self.assertEqual(medians, sorted(medians))
        self.assertIn(
            "brighter on the longest baselines", " ".join(report["reading"])
        )

    def test_the_uv_split_needs_at_least_two_ranges(self):
        with self.assertRaises(SystemExit) as raised:
            self.run_tool("--uv-split", "3", "8")

        self.assertIn("at least three increasing edges", str(raised.exception))
        self.assertEqual(self.commands, [])

    def test_a_given_channel_range_skips_the_measurement(self):
        report = self.run_tool("--channel-range", "10", "700")

        self.assertIsNone(self.flag_call)
        self.assertEqual(report["channel_range"], [10, 700])

    def test_the_noise_is_kept_when_the_image_fails(self):
        self.fail = "I"

        with self.assertRaises(RuntimeError):
            self.run_tool()

        report = json.loads(self.report_path().read_text())
        self.assertIn("robust_rms_Jy_beam", report["stokes_v"])
        self.assertNotIn("stokes_i", report)

    def test_existing_images_are_not_overwritten_unasked(self):
        self.run_tool()
        self.commands.clear()

        with self.assertRaises(SystemExit) as raised:
            self.run_tool()
        self.assertIn("--overwrite", str(raised.exception))
        self.assertEqual(self.commands, [])

        self.run_tool("--overwrite")
        self.assertEqual(len(self.commands), 2)

    def test_a_dry_run_prints_the_commands_and_runs_nothing(self):
        with mock.patch("builtins.print") as printed:
            TOOL.main([
                "--outdir", str(self.root / "out"), "--target", "G1",
                "--workdir", str(self.root / "work"), "--dry-run",
                "--channel-range", "0", "768",
            ])

        self.assertEqual(self.commands, [])
        self.assertFalse((self.root / "work").exists())
        lines = [call.args[0] for call in printed.call_args_list]
        self.assertTrue(lines[-1].startswith("wsclean "))
        self.assertIn("-pol I", lines[-1])

    def test_nothing_is_written_beside_the_pipeline_products(self):
        before = sorted(str(p) for p in (self.root / "out").rglob("*"))

        self.run_tool()

        self.assertEqual(
            before, sorted(str(p) for p in (self.root / "out").rglob("*"))
        )

    def test_a_missing_wsclean_is_reported_before_any_work(self):
        with mock.patch.object(TOOL.shutil, "which", lambda name: None):
            with self.assertRaises(SystemExit) as raised:
                self.run_tool()

        self.assertIn("not on PATH", str(raised.exception))
        self.assertIn("singularity exec", str(raised.exception))
        self.assertEqual(self.commands, [])

    def test_wsclean_in_a_container_is_run_through_it(self):
        self.capabilities = {"version": (3, 1), "options": set()}
        wrapper = "singularity exec /software/astro/containers/oxkat-latest.sif wsclean"

        report = self.run_tool("--wsclean", wrapper)

        for command in self.commands:
            self.assertEqual(command[:4], wrapper.split())
            self.assertIn("-use-wgridder", command)
        self.assertEqual(report["wsclean_version"], [3, 1])
        self.assertIn("singularity exec", report["commands"]["I"])

    def test_a_command_that_is_not_wsclean_stops_before_any_work(self):
        self.capabilities = {"version": None, "options": set()}

        with self.assertRaises(SystemExit) as raised:
            self.run_tool("--wsclean", "singularity exec missing.sif wsclean")

        self.assertIn("did not answer as WSClean", str(raised.exception))
        self.assertEqual(self.commands, [])


class NoiseStokesTests(unittest.TestCase):
    """Which dirty image gives the data's noise, from the correlations."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def choice(self, correlations=None, asked="auto"):
        if correlations is not None:
            path = self.root / "inspection" / "inspection_report.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(correlations))
        return TOOL.noise_stokes(TOOL.argparse.Namespace(
            noise_stokes=asked, outdir=str(self.root)))

    def test_without_a_record_it_is_stokes_v_as_before(self):
        self.assertEqual(self.choice(), "V")

    def test_the_parallel_linear_hands_alone_give_none(self):
        self.assertIsNone(self.choice(
            {"correlations": {"setups": {"0": ["XX", "YY"]}}}))

    def test_circular_or_all_four_linear_correlations_give_v(self):
        for correlations in (["RR", "LL"], ["RR", "RL", "LR", "LL"],
                             ["XX", "XY", "YX", "YY"]):
            with self.subTest(correlations=correlations):
                self.assertEqual(self.choice(
                    {"correlations": {"setups": {"0": correlations}}}), "V")

    def test_an_older_list_of_correlations_is_read_too(self):
        self.assertIsNone(self.choice({"correlations": ["XX", "YY"]}))

    def test_a_choice_given_is_kept(self):
        record = {"correlations": {"setups": {"0": ["XX", "YY"]}}}
        self.assertEqual(self.choice(record, "Q"), "Q")
        self.assertEqual(self.choice(record, "V"), "V")
        self.assertIsNone(self.choice({"correlations": ["RR", "LL"]}, "none"))


if __name__ == "__main__":
    unittest.main()
