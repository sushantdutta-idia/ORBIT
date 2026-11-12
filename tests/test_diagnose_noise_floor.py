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
        "test_diagnose_noise_floor_module",
        SOURCE_DIR / "diagnose_noise_floor.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TOOL = load_tool()


class StatisticsTests(unittest.TestCase):

    def setUp(self):
        self.random = np.random.default_rng(7)

    def test_robust_rms_is_the_noise_not_the_sources(self):
        noise = self.random.normal(0.0, 20e-6, 200000)
        noise[:200] = 5e-3

        self.assertAlmostEqual(TOOL.robust_rms(noise), 20e-6, delta=0.5e-6)
        self.assertGreater(TOOL.plain_rms(noise), 150e-6)

    def test_blanked_pixels_are_left_out(self):
        noise = self.random.normal(0.0, 20e-6, 20000)
        blanked = np.concatenate([noise, np.zeros(60000), [np.nan] * 100])

        self.assertAlmostEqual(
            TOOL.robust_rms(blanked), TOOL.robust_rms(noise)
        )
        self.assertIsNone(TOOL.robust_rms(np.zeros(100)))

    def test_the_grid_finds_the_noisy_corner(self):
        plane = self.random.normal(0.0, 20e-6, (500, 500))
        plane[400:, 400:] *= 6.0

        grid = TOOL.rms_grid(plane, 5)
        summary = TOOL.summarize_grid(grid, plane.shape)

        self.assertEqual(len(grid), 5)
        self.assertAlmostEqual(grid[0][0], 20e-6, delta=1.5e-6)
        self.assertAlmostEqual(grid[4][4], 120e-6, delta=8e-6)
        self.assertEqual(
            summary["noisiest_boxes"][0]["centre_pixel"], [450, 450]
        )
        self.assertAlmostEqual(
            summary["quietest_tenth_rms"], 20e-6, delta=1.5e-6
        )
        self.assertGreater(summary["noisiest_over_quietest"], 5.0)

    def test_the_brightest_pixel_is_located(self):
        plane = self.random.normal(0.0, 20e-6, (64, 64))
        plane[40, 12] = 0.025

        self.assertEqual(
            TOOL.brightest_pixel(plane), {"pixel": [40, 12], "value": 0.025}
        )

    def test_the_reading_reports_each_weighting_against_the_goal(self):
        concentrated = {
            "quietest_tenth_rms": 22e-6, "noisiest_over_quietest": 6.0,
        }
        weightings = [
            {
                "robust": 0.5,
                "stokes_v_dirty": {"robust_rms_Jy_beam": 23.6e-6},
                "stokes_v_cleaned": {"robust_rms_Jy_beam": 17.0e-6},
                "psf": {
                    "width_arcsec": 4.9, "largest_sidelobe": 0.44,
                    "far_sidelobe_rms": 0.0061,
                },
            },
            {
                "robust": -0.3,
                "stokes_v_dirty": {"robust_rms_Jy_beam": 31.0e-6},
                "stokes_v_cleaned": {"error": "RuntimeError: no"},
            },
        ]

        lines = TOOL.reading(weightings, concentrated, 20e-6)

        self.assertIn("robust 0.5", lines[0])
        self.assertIn("17.0 uJy/beam in Stokes V (cleaned)", lines[0])
        self.assertIn("at or below the 20", lines[0])
        self.assertIn("largest sidelobe 0.44", lines[0])
        # A clean that failed leaves the dirty image, called what it is.
        self.assertIn("robust -0.3", lines[1])
        self.assertIn("31.0 uJy/beam in Stokes V (dirty, so an upper limit)",
                      lines[1])
        self.assertIn("above the 20", lines[1])
        self.assertIn("concentrated", lines[2])

    def test_a_uniform_image_is_read_as_such(self):
        uniform = {
            "quietest_tenth_rms": 80e-6, "noisiest_over_quietest": 1.3,
        }

        self.assertIn("nearly uniform", TOOL.reading([], uniform, 20e-6)[0])

    def test_psf_statistics_measure_width_and_sidelobes(self):
        size = 801
        x = np.arange(size)[:, None] - 400
        y = np.arange(size)[None, :] - 400
        radius = np.hypot(x, y)
        fwhm = 5.0
        psf = np.exp(-4.0 * np.log(2.0) * (radius / fwhm) ** 2)
        psf[400, 430] = 0.30          # a near sidelobe
        far = radius > 150
        psf[far] += self.random.normal(0.0, 0.01, int(far.sum()))

        stats = TOOL.psf_statistics(psf * 3.0, cell_arcsec=2.0)

        self.assertAlmostEqual(stats["width_pixels"], fwhm, delta=0.6)
        self.assertAlmostEqual(stats["width_arcsec"], 2 * fwhm, delta=1.2)
        self.assertAlmostEqual(stats["largest_sidelobe"], 0.30, delta=0.01)
        self.assertAlmostEqual(stats["far_sidelobe_rms"], 0.01, delta=0.001)

    def test_an_empty_psf_has_no_statistics(self):
        self.assertIsNone(TOOL.psf_statistics(np.zeros((32, 32))))


class ProductTests(unittest.TestCase):

    def test_the_amplitude_corrected_data_are_preferred(self):
        with tempfile.TemporaryDirectory() as root:
            selfcal = Path(root) / "targets" / "G1" / "selfcal"
            for name in ("phase2.ms", "phase10.ms", "final_selfcal.ms"):
                (selfcal / "selfcal_ms" / name).mkdir(parents=True)
            (selfcal / "images" / "final.image.tt0").mkdir(parents=True)

            ms, image = TOOL.resolve_products(root, "G1")
            self.assertEqual(ms.name, "final_selfcal.ms")
            self.assertEqual(image.name, "final.image.tt0")

            (selfcal / "selfcal_ms" / "final_selfcal.ms").rmdir()
            ms, _ = TOOL.resolve_products(root, "G1")
            self.assertEqual(ms.name, "phase10.ms")


class FlowTests(unittest.TestCase):
    """The whole run against stand-ins for CASA."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        selfcal = self.root / "out" / "targets" / "G1" / "selfcal"
        self.ms = selfcal / "selfcal_ms" / "final_selfcal.ms"
        self.image = selfcal / "images" / "final.image.tt0"
        self.residual = selfcal / "images" / "final.residual.tt0"
        for path in (self.ms, self.image, self.residual):
            path.mkdir(parents=True)

        random = np.random.default_rng(3)
        stokes_i = random.normal(0.0, 25e-6, (200, 200))
        stokes_i[150:, 150:] *= 8.0
        stokes_i[170, 170] = 0.025
        self.planes = {
            str(self.image): stokes_i,
            str(self.residual): random.normal(0.0, 25e-6, (200, 200)),
        }
        self.calls = []
        test = self

        class Image:
            def open(self, path):
                self.path = str(path)

            def getchunk(self):
                return test.planes[self.path][:, :, None, None]

            def summary(self):
                return {"incr": [-1.0 / TOOL.ARCSEC_PER_RADIAN,
                                 1.0 / TOOL.ARCSEC_PER_RADIAN, 1.0, 1.0]}

            def close(self):
                pass

        class Metadata:
            def open(self, path):
                pass

            def scannumbers(self):
                return np.asarray([3, 15, 29])

            def close(self):
                pass

        def tclean(**kwargs):
            test.calls.append(kwargs)
            base = kwargs["imagename"]
            if kwargs["stokes"] == "I":
                if test.fail_residual:
                    raise RuntimeError("startmodel refused")
                # What the model leaves: noise, and errors in the first scan.
                scan = kwargs.get("scan", "")
                sigma = {"": 45e-6, "3": 150e-6, "15": 30e-6, "29": 95e-6}[
                    scan
                ]
                Path(base + ".residual.tt0").mkdir(parents=True)
                test.planes[base + ".residual.tt0"] = random.normal(
                    0.0, sigma, (200, 200)
                )
                return
            if kwargs["niter"] > 0:
                if test.fail_clean:
                    raise RuntimeError("clean refused")
                # Cleaning removes the leaked peaks and their sidelobes.
                test.planes[base + ".residual"] = random.normal(
                    0.0, test.floor[kwargs["robust"]], (200, 200)
                )
                return
            scan = kwargs.get("scan", "")
            if scan == "29" and test.fail_last_scan:
                raise RuntimeError("no valid data")
            sigma = {"": 1.4 * test.floor[kwargs["robust"]], "3": 30e-6,
                     "15": 21e-6, "29": 90e-6}[scan]
            for suffix in (".residual", ".psf"):
                Path(base + suffix).mkdir(parents=True, exist_ok=True)
            dirty = random.normal(0.0, sigma, (200, 200))
            dirty[20, 20] = 2.4e-3
            test.planes[base + ".residual"] = dirty
            psf = np.zeros((200, 200))
            psf[100, 100] = 1.0
            psf[100, 120] = 0.44 if kwargs["robust"] > 0 else 0.2
            test.planes[base + ".psf"] = psf

        self.floor = {0.5: 17e-6, 0.0: 21e-6, -0.3: 27e-6}
        self.fail_clean = False
        self.fail_residual = False
        for term in ("tt0", "tt1"):
            (selfcal / "images" / f"final.model.{term}").mkdir()
        self.fail_last_scan = False
        casatools = types.ModuleType("casatools")
        casatools.image = Image
        casatools.msmetadata = Metadata
        casatasks = types.ModuleType("casatasks")
        casatasks.tclean = tclean
        patch = mock.patch.dict(
            sys.modules, {"casatools": casatools, "casatasks": casatasks}
        )
        patch.start()
        self.addCleanup(patch.stop)
        self.addCleanup(self.temporary.cleanup)

    def run_tool(self, *extra):
        arguments = [
            "--outdir", str(self.root / "out"), "--target", "G1",
            "--workdir", str(self.root / "work"), *extra,
        ]
        with mock.patch("builtins.print"):
            TOOL.main(arguments)
        return json.loads(
            (self.root / "work" / "G1" / "noise_floor.json").read_text()
        )

    def test_it_measures_the_floor_and_where_the_image_is_noisy(self):
        report = self.run_tool()

        dirty, clean = self.calls[0], self.calls[1]
        self.assertEqual(dirty["vis"], str(self.ms))
        self.assertEqual(dirty["stokes"], "V")
        self.assertEqual(dirty["niter"], 0)
        # The grid and weighting of the image being judged.
        self.assertEqual(dirty["imsize"], 200)
        self.assertEqual(dirty["cell"], "1.0000arcsec")
        self.assertEqual(dirty["robust"], 0.5)
        self.assertNotIn("savemodel", dirty)
        # The clean continues from the dirty image instead of remaking it.
        self.assertEqual(clean["imagename"], dirty["imagename"])
        self.assertTrue(clean["restart"])
        self.assertFalse(clean["calcpsf"])
        self.assertFalse(clean["calcres"])
        self.assertEqual(
            (clean["niter"], clean["nsigma"], clean["nmajor"]), (5000, 5.0, 4)
        )
        self.assertEqual(
            [call.get("scan") for call in self.calls],
            [None, None, "3", "15", "29"],
        )

        (weighting,) = report["weightings"]
        self.assertAlmostEqual(
            weighting["stokes_v_dirty"]["robust_rms_Jy_beam"],
            1.4 * 17e-6, delta=1.5e-6,
        )
        self.assertAlmostEqual(
            weighting["stokes_v_dirty"]["max_abs_Jy_beam"], 2.4e-3
        )
        self.assertAlmostEqual(
            weighting["stokes_v_cleaned"]["robust_rms_Jy_beam"],
            17e-6, delta=1e-6,
        )
        self.assertAlmostEqual(weighting["psf"]["largest_sidelobe"], 0.44)
        self.assertAlmostEqual(
            report["stokes_v_scans"]["29"]["robust_rms_Jy_beam"],
            90e-6, delta=5e-6,
        )
        self.assertEqual(report["stokes_i"]["brightest"]["pixel"], [170, 170])
        self.assertGreater(
            report["stokes_i"]["grid"]["noisiest_over_quietest"], 5.0
        )
        self.assertIn("stokes_i_residual", report)
        self.assertIn("(cleaned), at or below the 20", report["reading"][0])
        self.assertIn("concentrated", report["reading"][1])

    def test_several_weightings_are_measured_side_by_side(self):
        report = self.run_tool("--robust", "0.5", "0", "-0.3", "--skip-scans")

        self.assertEqual(
            [entry["robust"] for entry in report["weightings"]],
            [0.5, 0.0, -0.3],
        )
        self.assertEqual(
            sorted({call["robust"] for call in self.calls}), [-0.3, 0.0, 0.5]
        )
        names = {Path(call["imagename"]).name for call in self.calls}
        self.assertEqual(len(names), 3)
        floors = [
            entry["stokes_v_cleaned"]["robust_rms_Jy_beam"]
            for entry in report["weightings"]
        ]
        self.assertEqual(floors, sorted(floors))
        self.assertAlmostEqual(
            report["weightings"][2]["psf"]["largest_sidelobe"], 0.2
        )
        self.assertNotIn("stokes_v_scans", report)
        self.assertIn("robust -0.3", report["reading"][2])
        self.assertIn("above the 20", report["reading"][2])

    def test_the_model_residual_is_set_against_the_noise_of_each_scan(self):
        report = self.run_tool("--model-residual")

        residual_calls = [c for c in self.calls if c["stokes"] == "I"]
        first = residual_calls[0]
        images = self.image.parent
        self.assertEqual(
            first["startmodel"],
            [str(images / "final.model.tt0"), str(images / "final.model.tt1")],
        )
        self.assertEqual(first["niter"], 0)
        self.assertEqual(first["nterms"], 2)
        # The gridder the model was made with, and no write to the data.
        self.assertEqual(first["gridder"], "wproject")
        self.assertEqual(first["wprojplanes"], 128)
        self.assertEqual(first["savemodel"], "none")
        self.assertEqual(first["robust"], 0.5)
        self.assertEqual(
            [call.get("scan") for call in residual_calls],
            [None, "3", "15", "29"],
        )

        weighting = report["weightings"][0]
        self.assertAlmostEqual(
            weighting["stokes_i_model_residual"]["robust_rms_Jy_beam"],
            45e-6, delta=3e-6,
        )
        self.assertIn("grid", weighting["stokes_i_model_residual"])
        self.assertEqual(list(report["scans"]), ["3", "15", "29"])
        self.assertAlmostEqual(
            report["scans"]["3"]["stokes_i"], 150e-6, delta=8e-6
        )
        self.assertAlmostEqual(
            report["scans"]["3"]["stokes_v"], 30e-6, delta=2e-6
        )
        text = " ".join(report["reading"])
        self.assertIn("Stokes I residual of the existing model is", text)
        self.assertIn("scan 3: Stokes I residual", text)
        self.assertIn("(5.0 times)", text)

    def test_without_the_option_no_model_is_imaged(self):
        report = self.run_tool("--skip-scans")

        self.assertFalse([c for c in self.calls if c["stokes"] == "I"])
        self.assertNotIn("stokes_i_model_residual", report["weightings"][0])
        self.assertEqual(report["scans"], {})

    def test_a_residual_that_fails_is_recorded_and_the_rest_kept(self):
        self.fail_residual = True

        report = self.run_tool("--model-residual")

        weighting = report["weightings"][0]
        self.assertIn(
            "startmodel refused", weighting["stokes_i_model_residual"]["error"]
        )
        self.assertIn("robust_rms_Jy_beam", weighting["stokes_v_cleaned"])
        self.assertIn(
            "startmodel refused",
            report["stokes_i_model_residual_scans"]["15"]["error"],
        )
        self.assertIn("robust_rms_Jy_beam", report["stokes_v_scans"]["15"])

    def test_the_model_residual_needs_model_images(self):
        for term in ("tt0", "tt1"):
            (self.image.parent / f"final.model.{term}").rmdir()

        with self.assertRaises(SystemExit) as raised:
            self.run_tool("--model-residual")

        self.assertIn("none were found", str(raised.exception))
        self.assertEqual(self.calls, [])

    def test_a_clean_that_fails_leaves_the_dirty_measurement(self):
        self.fail_clean = True

        report = self.run_tool("--skip-scans")

        weighting = report["weightings"][0]
        self.assertIn("clean refused", weighting["stokes_v_cleaned"]["error"])
        self.assertIn("robust_rms_Jy_beam", weighting["stokes_v_dirty"])
        self.assertIn("dirty, so an upper limit", report["reading"][0])

    def test_the_clean_can_be_left_out(self):
        report = self.run_tool("--skip-scans", "--v-niter", "0")

        self.assertEqual(len(self.calls), 1)
        self.assertNotIn("stokes_v_cleaned", report["weightings"][0])

    def test_a_scan_that_cannot_be_imaged_is_recorded_not_fatal(self):
        self.fail_last_scan = True

        report = self.run_tool()

        self.assertIn(
            "no valid data", report["stokes_v_scans"]["29"]["error"]
        )
        self.assertIn("robust_rms_Jy_beam", report["stokes_v_scans"]["15"])

    def test_a_rerun_replaces_its_own_images(self):
        self.run_tool("--skip-scans")
        stale = self.root / "work" / "G1" / "stokes_v_robust_0p5.mask"
        stale.mkdir()

        self.run_tool("--skip-scans")

        self.assertFalse(stale.exists())
        self.assertEqual(len(self.calls), 4)

    def test_nothing_is_written_beside_the_pipeline_products(self):
        before = sorted(
            str(path) for path in (self.root / "out").rglob("*")
        )

        self.run_tool()

        self.assertEqual(
            before,
            sorted(str(path) for path in (self.root / "out").rglob("*")),
        )


if __name__ == "__main__":
    unittest.main()
