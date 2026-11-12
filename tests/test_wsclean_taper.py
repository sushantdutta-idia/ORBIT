"""A Gaussian taper for extended emission of low surface brightness."""

import importlib.util
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_wsclean_taper_imager", SOURCE_DIR / "experiment_wsclean_image.py")
IMAGER = importlib.util.module_from_spec(spec)
spec.loader.exec_module(IMAGER)


def arguments(*extra):
    parser_args = []
    captured = []
    original = IMAGER.run
    IMAGER.run = captured.append
    try:
        IMAGER.main(["--target", "G1", "--ms", "x.ms", *extra])
    finally:
        IMAGER.run = original
    return captured[0]


class TaperTests(unittest.TestCase):

    def test_the_taper_goes_to_both_polarizations(self):
        args = arguments("--taper", "15", "--robust", "0")
        for stokes in ("I", "V"):
            command = IMAGER.wsclean_command(
                args, "x.ms", "out/name", 1024, 3.0, (31, 993), stokes, "tmp")
            index = command.index("-taper-gaussian")
            self.assertEqual(command[index + 1], "15asec")
            self.assertEqual(command[command.index("-weight") + 2], "0")

    def test_without_a_taper_the_command_is_as_it_was(self):
        command = IMAGER.wsclean_command(
            arguments(), "x.ms", "out/name", 1024, 3.0, (31, 993), "I", "tmp")
        self.assertNotIn("-taper-gaussian", command)

    def test_the_taper_is_named_in_the_images(self):
        self.assertEqual(IMAGER.image_label(0), "robust_0")
        self.assertEqual(IMAGER.image_label(0.5), "robust_0p5")
        self.assertEqual(IMAGER.image_label(0, 15), "robust_0_taper15")
        self.assertEqual(IMAGER.image_label(-0.5, 7.5), "robust_m0p5_taper7p5")


if __name__ == "__main__":
    unittest.main()
