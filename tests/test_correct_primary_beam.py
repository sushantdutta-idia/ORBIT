import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_tool():
    spec = importlib.util.spec_from_file_location(
        "test_correct_primary_beam_module",
        SOURCE_DIR / "correct_primary_beam.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TOOL = load_tool()


def write_image(path, plane, frequency_hz, weight=None, cell_arcsec=30.0):
    """A FITS file as WSClean writes its images, with its image weight."""

    ny, nx = plane.shape
    cards = [
        "SIMPLE  =                    T",
        "BITPIX  =                  -32",
        "NAXIS   =                    4",
        f"NAXIS1  = {nx:20d}",
        f"NAXIS2  = {ny:20d}",
        "NAXIS3  =                    1",
        "NAXIS4  =                    1",
        "CTYPE1  = 'RA---SIN'",
        f"CRPIX1  = {nx // 2 + 1:20d}",
        f"CDELT1  = {-cell_arcsec / 3600:20.12f}",
        "CTYPE2  = 'DEC--SIN'",
        f"CRPIX2  = {ny // 2 + 1:20d}",
        f"CDELT2  = {cell_arcsec / 3600:20.12f}",
        "CTYPE3  = 'FREQ'",
        f"CRVAL3  = {frequency_hz:20.4f}",
    ]
    if weight is not None:
        cards.append(f"WSCIMGWG= {weight:20.6f}")
    cards.append("END")
    header = "".join(card.ljust(80) for card in cards)
    header = header.ljust(-(-len(header) // 2880) * 2880).encode("ascii")
    data = plane.astype(">f4").tobytes()
    path.write_bytes(header + data + b"\0" * (-len(data) % 2880))


class StartTests(unittest.TestCase):

    def test_it_starts_beside_the_pipelines_own_inspect_module(self):
        source = (SOURCE_DIR / "correct_primary_beam.py").read_text()
        self.assertIn("sys.path.pop(0)", source.split("import argparse")[0])
        done = subprocess.run(
            [sys.executable, "correct_primary_beam.py", "--help"],
            cwd=SOURCE_DIR, capture_output=True, text=True,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("--pbmin", done.stdout)

    def test_the_beam_of_band_4_is_the_pipelines_own(self):
        coefficients, valid = TOOL.band_coefficients("B4")

        self.assertEqual(coefficients, [-3.263, 42.618, -25.580, 5.823])
        self.assertEqual(valid, [550.0, 850.0])
        with self.assertRaises(SystemExit):
            TOOL.band_coefficients("B9")


class BeamTests(unittest.TestCase):

    def setUp(self):
        self.coefficients, _ = TOOL.band_coefficients("B4")
        self.bands = [
            {"frequency_hz": 600e6, "weight": 3.0},
            {"frequency_hz": 800e6, "weight": 1.0},
        ]

    def beam(self, pbmin=0.3):
        return TOOL.mean_beam(
            (201, 201), 0.5, self.bands, self.coefficients, pbmin
        )

    def test_it_is_one_at_the_pointing_centre_and_falls_outwards(self):
        beam, limit, frequency = self.beam()

        self.assertAlmostEqual(beam[100, 100], 1.0)
        self.assertLess(beam[100, 130], beam[100, 110])
        self.assertAlmostEqual(beam[100, 130], beam[130, 100])
        self.assertAlmostEqual(frequency, 650e6)

    def test_it_is_the_weighted_mean_of_the_sub_bands_beams(self):
        beam, _, _ = self.beam()
        radius = 15.0                      # arcmin, 30 pixels of half one
        low = TOOL.SAFETY.ugmrt_primary_beam_response(
            radius * 0.6, self.coefficients)
        high = TOOL.SAFETY.ugmrt_primary_beam_response(
            radius * 0.8, self.coefficients)

        self.assertAlmostEqual(beam[100, 130], 0.75 * low + 0.25 * high)
        # A wider beam at the lower frequency, which has the weight.
        self.assertGreater(beam[100, 130], high)

    def test_outside_the_main_lobe_of_the_narrowest_beam_it_is_blank(self):
        beam, limit, _ = self.beam()
        cut = TOOL.SAFETY.first_outward_beam_cutoff(self.coefficients, 0.3)

        self.assertAlmostEqual(limit, cut / 0.8)
        edge = int(limit / 0.5)
        self.assertGreater(beam[100, 100 + edge], 0.0)
        self.assertEqual(beam[100, 100 + edge + 2], 0.0)
        self.assertEqual(beam[0, 0], 0.0)
        self.assertGreaterEqual(beam[beam > 0].min(), 0.3)

    def test_a_lower_cut_off_reaches_farther(self):
        _, near, _ = self.beam(0.5)
        _, far, _ = self.beam(0.2)

        self.assertLess(near, far)


class FlowTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.name = self.root / "wsclean_robust_0_I"
        self.plane = np.full((121, 121), 2.0e-3)
        write_image(Path(str(self.name) + "-MFS-image.fits"), self.plane, 700e6)
        for index, (frequency, weight) in enumerate(
            [(600e6, 3.0), (800e6, 1.0), (840e6, 0.0)]
        ):
            write_image(
                Path(f"{self.name}-{index:04d}-image.fits"),
                self.plane, frequency, weight,
            )

    def run_tool(self, *extra):
        with mock.patch("builtins.print"):
            TOOL.main(["--name", str(self.name), *extra])
        return json.loads(Path(str(self.name) + "-MFS-pbcor.json").read_text())

    def test_the_image_is_divided_by_its_beam_and_blanked_beyond_it(self):
        report = self.run_tool()

        _, _, corrected, _ = TOOL.read_fits(
            Path(str(self.name) + "-MFS-image.pbcor.fits"))
        _, _, beam, _ = TOOL.read_fits(Path(str(self.name) + "-MFS-pb.fits"))
        self.assertAlmostEqual(corrected[60, 60], 2.0e-3, places=9)
        self.assertAlmostEqual(
            corrected[60, 80], 2.0e-3 / beam[60, 80], places=8)
        self.assertGreater(corrected[60, 80], 2.0e-3)
        self.assertTrue(np.isnan(corrected[0, 0]))
        self.assertTrue(np.isnan(beam[0, 0]))
        # A sub-band without data has no say in the beam.
        self.assertEqual(report["sub_bands"], 2)
        self.assertAlmostEqual(report["mean_frequency_mhz"], 650.0)
        self.assertGreaterEqual(report["beam_at_edge"], 0.3)
        self.assertEqual(report["coefficients"],
                         [-3.263, 42.618, -25.580, 5.823])

    def test_the_image_itself_is_left_as_it_was(self):
        before = Path(str(self.name) + "-MFS-image.fits").read_bytes()

        self.run_tool()

        self.assertEqual(
            Path(str(self.name) + "-MFS-image.fits").read_bytes(), before)
        self.assertEqual(
            Path(str(self.name) + "-MFS-image.pbcor.fits").stat().st_size,
            len(before),
        )

    def test_an_earlier_correction_is_not_overwritten_unasked(self):
        self.run_tool()

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("already exists", str(stopped.exception))
        self.run_tool("--overwrite")

    def test_without_the_sub_band_images_there_is_no_beam(self):
        for path in self.root.glob("*-000?-image.fits"):
            path.unlink()

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("No sub-band images", str(stopped.exception))

    def test_a_sub_band_outside_the_measured_range_is_refused(self):
        write_image(Path(f"{self.name}-0003-image.fits"), self.plane, 900e6, 1.0)

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("outside 550 to 850 MHz", str(stopped.exception))

    def test_an_image_without_weights_counts_each_sub_band_alike(self):
        for path in self.root.glob("*-000?-image.fits"):
            path.unlink()
        write_image(Path(f"{self.name}-0000-image.fits"), self.plane, 600e6)
        write_image(Path(f"{self.name}-0001-image.fits"), self.plane, 800e6)

        report = self.run_tool()

        self.assertAlmostEqual(report["mean_frequency_mhz"], 700.0)

    def test_a_missing_image_is_said_at_once(self):
        with self.assertRaises(SystemExit) as stopped:
            with mock.patch("builtins.print"):
                TOOL.main(["--name", str(self.root / "absent")])

        self.assertIn("No image", str(stopped.exception))


class MeerKATTests(unittest.TestCase):
    """The katbeam correction the MeerKAT target chain asks for.

    Built as 1708779077 is: 856-1712 MHz in eight sub-bands, imaging weights
    that rise with frequency, the 1163-1299 MHz mask emptying the fourth
    sub-band (no weight, no finite pixel) and the other two masks cutting
    the first and the last two.
    """

    EDGES = np.linspace(856e6, 1712e6, 9)
    MASKS_MHZ = [[933.0, 960.0], [1163.0, 1299.0], [1524.0, 1630.0]]
    WEIGHTS = [0.8, 3.6, 4.1, 0.0, 6.4, 9.9, 2.8, 10.1]
    FLUX = 1.5

    def setUp(self):
        from test_meerkat_beam import (
            FakeJimBeam, NX, NY, geometry, wsclean_header, write_image as write,
        )

        self.FakeJimBeam = FakeJimBeam
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.products = self.root / "MeerKAT_Products"
        report = self.products / "flagging" / "known_frequency_mask_report.json"
        report.parent.mkdir(parents=True)
        report.write_text(json.dumps({"applied_intersections": [
            {"spw": "0", "intersection_MHz": item} for item in self.MASKS_MHZ
        ]}))
        images = self.root / "images"
        images.mkdir()
        self.name = images / "wsclean_robust_0_I"
        beam = TOOL.katbeam_module()
        sky = FakeJimBeam("m")
        value = geometry()
        masks = [(low * 1e6, high * 1e6) for low, high in self.MASKS_MHZ]
        weighted = np.zeros((NY, NX))
        for index, (low, high) in enumerate(zip(self.EDGES[:-1],
                                                self.EDGES[1:])):
            centre, width = 0.5 * (low + high), high - low
            effective, _ = beam.effective_band(
                centre, width, (900e6, 1670e6), masks)
            weight = self.WEIGHTS[index]
            if weight > 0:
                plane = self.FLUX * beam.beam_plane(sky, value, effective / 1e6)
                weighted += weight * plane
            else:
                plane = np.full((NY, NX), np.nan)
            write(Path(f"{self.name}-{index:04d}-image.fits"), plane,
                  wsclean_header(centre, width, WSCIMGWG=weight))
        write(Path(f"{self.name}-MFS-image.fits"),
              weighted / sum(self.WEIGHTS), wsclean_header(1284e6, 856e6))

    def correct(self, *extra, band="MEERKAT_L"):
        with mock.patch.object(
            TOOL.katbeam_module(), "load_katbeam",
            return_value=(self.FakeJimBeam, "0.1.2"),
        ), mock.patch("builtins.print"):
            TOOL.main(["--name", str(self.name), "--telescope", "meerkat",
                       "--band", band, "--products", str(self.products),
                       *extra])
        return json.loads(
            Path(f"{self.name}-MFS-pbcor.json").read_text())

    @staticmethod
    def plane(path):
        from astropy.io import fits

        with fits.open(str(path)) as handle:
            return handle[0].data[0, 0].copy()

    def test_the_field_comes_back_at_its_true_flux(self):
        report = self.correct()

        self.assertTrue(report["corrected"].endswith(
            "wsclean_robust_0_I-MFS-image-pb.fits"))
        corrected = self.plane(report["corrected"])
        beam = self.plane(report["beam"])
        inside = beam >= 0.3
        self.assertGreater(inside.sum(), 300)
        np.testing.assert_allclose(corrected[inside], self.FLUX, rtol=2e-4)
        self.assertEqual(report["mfs_weighting"], "wsclean_image_weight")
        self.assertEqual(report["katbeam_model"], "MKAT-AA-L-JIM-2020")
        self.assertTrue(Path(report["provenance"]).is_file())

    def test_the_flagging_masks_are_the_products_own(self):
        report = self.correct()

        self.assertEqual(report["masked_ranges_hz"], [
            [low * 1e6, high * 1e6] for low, high in self.MASKS_MHZ])
        self.assertEqual(report["sub_bands"], 7)
        self.assertEqual(len(report["skipped_sub_bands"]), 1)
        self.assertTrue(report["skipped_sub_bands"][0]["image"].endswith(
            "-0003-image.fits"))

    def test_without_the_mask_record_no_mask_is_assumed(self):
        self.assertEqual(TOOL.flagging_masks_hz(None), [])
        self.assertEqual(TOOL.flagging_masks_hz(self.root / "absent"), [])
        self.assertEqual(len(TOOL.flagging_masks_hz(self.products)), 3)

    def test_a_band_without_a_measured_beam_is_refused(self):
        with self.assertRaises(SystemExit) as stopped:
            self.correct(band="MEERKAT_S0")

        self.assertIn("no measured katbeam model", str(stopped.exception))

    def test_meerkat_needs_its_band(self):
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            TOOL.main(["--name", str(self.name), "--telescope", "meerkat"])

    def test_an_earlier_correction_is_not_overwritten_unasked(self):
        self.correct()

        with self.assertRaises(SystemExit) as stopped:
            self.correct()

        self.assertIn("already exists", str(stopped.exception))
        self.correct("--overwrite")


if __name__ == "__main__":
    unittest.main()
