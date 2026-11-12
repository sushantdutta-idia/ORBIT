import importlib.util
import json
import math
import subprocess
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
        "test_experiment_extract_source_module",
        SOURCE_DIR / "experiment_extract_source.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TOOL = load_tool()
IMAGER = TOOL.load_imager()
CASA = "singularity exec /containers/casa.sif python3"
WSCLEAN = "singularity exec /containers/oxkat.sif wsclean"


def write_model(path, plane, crval=(33.0, -3.5), cell_arcsec=1.0):
    """A FITS file as WSClean writes its images: four axes, sine projection."""

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
        f"CRVAL1  = {crval[0]:20.10f}",
        f"CDELT1  = {-cell_arcsec / 3600:20.12f}",
        "CTYPE2  = 'DEC--SIN'",
        f"CRPIX2  = {ny // 2 + 1:20d}",
        f"CRVAL2  = {crval[1]:20.10f}",
        f"CDELT2  = {cell_arcsec / 3600:20.12f} / degrees",
        "END",
    ]
    header = "".join(card.ljust(80) for card in cards)
    header = header.ljust(-(-len(header) // 2880) * 2880).encode("ascii")
    data = plane.astype(">f4").tobytes()
    path.write_bytes(header + data + b"\0" * (-len(data) % 2880))


def arguments(*extra):
    captured = []
    with mock.patch.object(TOOL, "run", captured.append):
        TOOL.main([
            "--ms", "data.ms", "--models", "images/run_I",
            "--pixel", "30", "20", "--casa", CASA, "--wsclean", WSCLEAN,
            *extra,
        ])
    return captured[0]


class StartTests(unittest.TestCase):

    def test_it_starts_beside_the_pipelines_own_inspect_module(self):
        source = (SOURCE_DIR / "experiment_extract_source.py").read_text()
        self.assertIn("sys.path.pop(0)", source.split("import argparse")[0])
        done = subprocess.run(
            [sys.executable, "experiment_extract_source.py", "--help"],
            cwd=SOURCE_DIR, capture_output=True, text=True,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("--pixel", done.stdout)
        # What it passes to itself under CASA is not offered to the user.
        self.assertNotIn("--casa-part", done.stdout)

    def test_the_source_and_the_models_have_to_be_named(self):
        with self.assertRaises(SystemExit), mock.patch("sys.stderr"):
            TOOL.main(["--ms", "data.ms", "--models", "run_I"])


class ModelTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def test_an_image_is_read_and_written_back_byte_for_byte(self):
        plane = np.arange(48 * 64, dtype=float).reshape(48, 64)
        path = self.root / "run_I-0000-model.fits"
        write_model(path, plane)

        header, cards, read, tail = TOOL.read_fits(path)
        TOOL.write_fits(self.root / "copy.fits", header, read, tail)

        self.assertEqual(read.shape, (48, 64))
        self.assertEqual(read[20, 30], plane[20, 30])
        self.assertEqual(cards["CTYPE1"], "RA---SIN")
        self.assertEqual(float(cards["CDELT2"]), round(1 / 3600, 12))
        self.assertEqual(
            (self.root / "copy.fits").read_bytes(), path.read_bytes()
        )

    def test_a_file_cut_short_is_refused(self):
        path = self.root / "short.fits"
        path.write_bytes(b"SIMPLE  =                    T".ljust(100))

        with self.assertRaises(ValueError):
            TOOL.read_fits(path)

    def test_only_the_disc_about_the_source_is_taken_out(self):
        plane = np.ones((60, 80))
        plane[20, 30] = 5.0

        rest, taken = TOOL.without_disc(plane, 30, 20, 4)

        self.assertEqual(rest[20, 30], 0.0)
        self.assertEqual(rest[20, 34], 0.0)
        self.assertEqual(rest[20, 35], 1.0)
        self.assertEqual(rest[23, 33], 1.0)     # outside the circle
        self.assertEqual(plane[20, 30], 5.0)    # the original is kept
        self.assertAlmostEqual(taken + rest.sum(), plane.sum())
        self.assertAlmostEqual(taken, 49 - 1 + 5.0)   # 49 pixels in the disc

    def test_a_source_at_the_edge_loses_only_what_is_in_the_image(self):
        plane = np.ones((20, 20))

        rest, taken = TOOL.without_disc(plane, 0, 0, 3)

        self.assertEqual(rest[0, 0], 0.0)
        self.assertAlmostEqual(taken + rest.sum(), 400.0)

    def test_the_centre_pixel_is_the_phase_centre(self):
        write_model(self.root / "m.fits", np.zeros((100, 100)))
        _, cards, _, _ = TOOL.read_fits(self.root / "m.fits")

        ra, dec = TOOL.pixel_direction(cards, 50, 50)

        self.assertAlmostEqual(ra, 33.0)
        self.assertAlmostEqual(dec, -3.5)

    def test_north_is_up_and_east_is_to_the_left(self):
        write_model(self.root / "m.fits", np.zeros((100, 100)),
                    cell_arcsec=36.0)
        _, cards, _, _ = TOOL.read_fits(self.root / "m.fits")

        ra_n, dec_n = TOOL.pixel_direction(cards, 50, 60)   # ten pixels up
        ra_e, dec_e = TOOL.pixel_direction(cards, 40, 50)   # ten to the left

        self.assertAlmostEqual(ra_n, 33.0, places=6)
        self.assertAlmostEqual(dec_n, -3.4, places=4)
        self.assertAlmostEqual(
            (ra_e - 33.0) * math.cos(math.radians(-3.5)), 0.1, places=4
        )
        self.assertAlmostEqual(dec_e, -3.5, places=3)

    def test_the_sub_band_models_are_found_in_order(self):
        for index in (2, 0, 1):
            (self.root / f"run_I-{index:04d}-model.fits").write_bytes(b"")
        (self.root / "run_I-MFS-model.fits").write_bytes(b"")
        (self.root / "run_I-0000-image.fits").write_bytes(b"")

        found = TOOL.model_images(self.root / "run_I")

        self.assertEqual(
            [path.name for path in found],
            [f"run_I-{index:04d}-model.fits" for index in range(3)],
        )

    def test_a_name_with_no_models_is_said_so(self):
        with self.assertRaises(SystemExit) as stopped:
            TOOL.model_images(self.root / "absent")

        self.assertIn("-NNNN-model.fits", str(stopped.exception))

    def test_each_sub_band_gets_a_model_of_the_rest_of_the_field(self):
        images = []
        for index in range(3):
            plane = np.zeros((60, 80))
            plane[20, 30] = 0.050 / (index + 1)   # the source, fading
            plane[45, 60] = 0.010                 # another source
            path = self.root / f"run_I-{index:04d}-model.fits"
            write_model(path, plane)
            images.append(path)
        folder = self.root / "out"
        folder.mkdir()

        cards, table = TOOL.split_models(images, folder, 30, 20, 5)

        self.assertEqual(cards["NAXIS1"], "80")
        self.assertEqual([row["sub_band"] for row in table], [0, 1, 2])
        self.assertAlmostEqual(table[0]["source_model_Jy"], 0.050)
        self.assertAlmostEqual(table[2]["source_model_Jy"], 0.050 / 3)
        self.assertAlmostEqual(table[1]["rest_model_Jy"], 0.010)
        _, _, rest, _ = TOOL.read_fits(folder / "rest-0001-model.fits")
        self.assertEqual(rest[20, 30], 0.0)
        self.assertAlmostEqual(float(rest[45, 60]), 0.010)
        self.assertEqual(
            (folder / "rest-0001-model.fits").stat().st_size,
            images[1].stat().st_size,
        )

    def test_a_pixel_outside_the_image_is_refused(self):
        path = self.root / "run_I-0000-model.fits"
        write_model(path, np.zeros((60, 80)))

        with self.assertRaises(SystemExit) as stopped:
            TOOL.split_models([path], self.root, 90, 20, 5)

        self.assertIn("outside", str(stopped.exception))


class CommandTests(unittest.TestCase):

    def test_wsclean_predicts_the_rest_of_the_field_into_the_copy(self):
        command = TOOL.predict_command(
            arguments("--channel-range", "51", "768"), "work/field.ms",
            "work/rest", 12, {"version": (3, 1), "options": set()}, IMAGER,
        )

        self.assertEqual(command[:4], WSCLEAN.split())
        self.assertIn("-predict", command)
        self.assertEqual(command[command.index("-name") + 1], "work/rest")
        self.assertEqual(command[command.index("-channels-out") + 1], "12")
        index = command.index("-channel-range")
        self.assertEqual(command[index + 1:index + 3], ["51", "768"])
        self.assertIn("-use-wgridder", command)
        self.assertEqual(command[-1], "work/field.ms")

    def test_casa_is_given_this_script_and_where_the_source_is(self):
        command = TOOL.casa_command(
            arguments(), "w/field.ms", "w/source1.ms", (33.1234, -3.0425),
            "w/subtracted.ms",
        )

        self.assertEqual(command[:4], CASA.split())
        self.assertTrue(command[4].endswith("experiment_extract_source.py"))
        self.assertEqual(
            command[command.index("--casa-part") + 1], "subtract-shift"
        )
        self.assertEqual(
            command[command.index("--subtracted-ms") + 1], "w/subtracted.ms"
        )
        self.assertEqual(
            command[command.index("--phasecenter") + 1],
            "J2000 33.12340000deg -3.04250000deg",
        )
        self.assertEqual(command[command.index("--out-ms") + 1],
                         "w/source1.ms")
        self.assertNotIn("--original-ms", command)
        command = TOOL.casa_command(
            arguments(), "w/field.ms", "w/source1.ms", (33.1234, -3.0425),
            "w/subtracted.ms", "w/source1/source1.ms",
        )
        self.assertEqual(command[command.index("--original-ms") + 1],
                         "w/source1/source1.ms")

    def test_under_casa_the_model_is_subtracted_and_the_centre_moved(self):
        calls = []
        root = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(root))
        field = root / "field.ms"
        field.mkdir()
        (root / "source1.ms").mkdir()      # left by an earlier attempt

        def task(name):
            def run(**kwargs):
                calls.append((name, kwargs))
                if "outputvis" in kwargs:
                    Path(kwargs["outputvis"]).mkdir()
            return run

        casatasks = types.ModuleType("casatasks")
        casatasks.uvsub = task("uvsub")
        casatasks.split = task("split")
        casatasks.phaseshift = task("phaseshift")
        with mock.patch.dict(sys.modules, {"casatasks": casatasks}), \
                mock.patch("builtins.print"):
            status = TOOL.main([
                "--casa-part", "subtract-shift", "--field-ms", str(field),
                "--subtracted-ms", str(root / "subtracted.ms"),
                "--out-ms", str(root / "source1.ms"),
                "--phasecenter", "J2000 33.1deg -3.0deg",
            ])

        self.assertEqual(status, 0)
        self.assertEqual([name for name, _ in calls],
                         ["uvsub", "split", "phaseshift"])
        self.assertEqual(calls[0][1], {"vis": str(field)})
        self.assertEqual(calls[1][1]["datacolumn"], "corrected")
        self.assertEqual(calls[2][1]["vis"], str(root / "subtracted.ms"))
        self.assertEqual(calls[2][1]["phasecenter"], "J2000 33.1deg -3.0deg")
        self.assertTrue((root / "source1.ms").is_dir())


class AsCalibratedTests(unittest.TestCase):
    """The source taken out of the field's own data, through its gains."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        self.original, self.calibrated = root / "source1.ms", root / "round2.ms"
        self.out = root / "subtracted.ms"
        for path in (self.original, self.calibrated):
            path.mkdir()
            (path / "table.dat").write_text("rows")
        rng = np.random.default_rng(7)
        rows, shape = 50, (2, 6, 50)
        self.data = (rng.normal(size=shape) + 1j * rng.normal(size=shape))
        # Gains as a noisy faint-source solve gives them: 0.4 to 1.4.
        self.gain = rng.uniform(0.4, 1.4, shape) * np.exp(
            1j * rng.uniform(-1, 1, shape))
        self.model = 0.01 * np.exp(1j * rng.uniform(-3, 3, shape))
        seen = self.data / self.gain
        seen[0, 0, 3] = 0.0           # a visibility the solve left at zero
        layout = {"ANTENNA1": np.arange(rows) % 7,
                  "ANTENNA2": np.arange(rows) % 5,
                  "TIME": np.linspace(0, 700, rows)}
        self.tables = {
            str(self.original): {**layout, "DATA": self.data.copy()},
            str(self.calibrated): {**layout, "DATA": seen,
                                   "MODEL_DATA": self.model.copy()},
        }
        test = self

        class Table:
            def open(self, path, nomodify=True):
                path = str(path)
                if path not in test.tables:
                    test.tables[path] = {
                        key: value.copy() for key, value in
                        test.tables[str(test.original)].items()}
                self.columns = test.tables[path]

            def nrows(self):
                return len(self.columns["TIME"])

            def getcol(self, name, start=0, count=-1):
                values = self.columns[name]
                end = None if count < 0 else start + count
                return values[..., start:end].copy()

            def putcol(self, name, values, start=0, count=-1):
                self.columns[name][..., start:start + values.shape[-1]] = values

            def close(self):
                pass

        casatools = types.ModuleType("casatools")
        casatools.table = Table
        patch = mock.patch.dict(sys.modules, {"casatools": casatools})
        patch.start()
        self.addCleanup(patch.stop)

    def test_the_model_goes_out_as_the_field_sees_it(self):
        undefined, median = TOOL.subtract_as_calibrated(
            self.calibrated, self.original, self.out, chunk=16)

        result = self.tables[str(self.out)]["DATA"]
        expected = self.data - self.model * self.gain
        expected[0, 0, 3] = self.data[0, 0, 3] - self.model[0, 0, 3]
        np.testing.assert_allclose(result, expected, atol=1e-12)
        self.assertAlmostEqual(undefined, 1 / self.data.size)
        self.assertGreater(median, 0.4)
        self.assertLess(median, 1.4)
        # The data given are not written to.
        np.testing.assert_array_equal(
            self.tables[str(self.original)]["DATA"], self.data)
        self.assertTrue((self.out / "table.dat").is_file())

    def test_data_not_the_same_row_for_row_are_refused(self):
        self.tables[str(self.calibrated)]["ANTENNA2"] = (
            self.tables[str(self.calibrated)]["ANTENNA2"][::-1].copy())
        with self.assertRaises(SystemExit) as stopped:
            TOOL.subtract_as_calibrated(
                self.calibrated, self.original, self.out)
        self.assertIn("ANTENNA2", str(stopped.exception))

    def test_under_casa_the_original_data_are_used_and_the_centre_moved(self):
        calls = []

        def task(name):
            def run(**kwargs):
                calls.append((name, kwargs))
                if "outputvis" in kwargs:
                    Path(kwargs["outputvis"]).mkdir()
            return run

        casatasks = types.ModuleType("casatasks")
        casatasks.uvsub = task("uvsub")
        casatasks.split = task("split")
        casatasks.phaseshift = task("phaseshift")
        out = self.original.parent / "residual.ms"
        with mock.patch.dict(sys.modules, {"casatasks": casatasks}), \
                mock.patch("builtins.print"):
            status = TOOL.main([
                "--casa-part", "subtract-shift",
                "--field-ms", str(self.calibrated),
                "--original-ms", str(self.original),
                "--subtracted-ms", str(self.out),
                "--out-ms", str(out),
                "--phasecenter", "J2000 33.1deg -3.0deg",
            ])

        self.assertEqual(status, 0)
        self.assertEqual([name for name, _ in calls], ["phaseshift"])
        self.assertEqual(calls[0][1]["vis"], str(self.out))
        np.testing.assert_allclose(
            self.tables[str(self.out)]["DATA"][:, :, 10],
            (self.data - self.model * self.gain)[:, :, 10], atol=1e-12)


class AddBackTests(unittest.TestCase):

    def test_without_a_direction_casa_is_asked_to_add_the_model(self):
        command = TOOL.casa_command(arguments(), "w/residual.ms", "w/out.ms")

        self.assertEqual(command[command.index("--casa-part") + 1], "add")
        self.assertNotIn("--phasecenter", command)
        self.assertNotIn("--subtracted-ms", command)

    def test_under_casa_the_model_column_is_added_and_the_result_kept(self):
        calls = []
        root = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(root))
        field = root / "residual.ms"
        field.mkdir()

        def task(name):
            def run(**kwargs):
                calls.append((name, kwargs))
                if "outputvis" in kwargs:
                    Path(kwargs["outputvis"]).mkdir()
            return run

        casatasks = types.ModuleType("casatasks")
        casatasks.uvsub = task("uvsub")
        casatasks.split = task("split")
        with mock.patch.dict(sys.modules, {"casatasks": casatasks}), \
                mock.patch("builtins.print"):
            TOOL.main([
                "--casa-part", "add", "--field-ms", str(field),
                "--out-ms", str(root / "peeled.ms"),
            ])

        self.assertEqual([name for name, _ in calls], ["uvsub", "split"])
        self.assertEqual(calls[0][1], {"vis": str(field), "reverse": True})
        self.assertEqual(calls[1][1]["datacolumn"], "corrected")
        self.assertTrue((root / "peeled.ms").is_dir())

    def test_sources_already_out_are_left_out_of_the_rest_as_well(self):
        root = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(root))
        plane = np.zeros((60, 80))
        plane[20, 30] = 0.050      # the source
        plane[45, 60] = 0.010      # one taken out earlier
        plane[10, 70] = 0.002      # the rest
        write_model(root / "run_I-0000-model.fits", plane)

        _, table = TOOL.split_models(
            [root / "run_I-0000-model.fits"], root, 30, 20, 5,
            others=[(60, 45)],
        )

        _, _, rest, _ = TOOL.read_fits(root / "rest-0000-model.fits")
        self.assertEqual(rest[45, 60], 0.0)
        self.assertEqual(rest[20, 30], 0.0)
        self.assertAlmostEqual(float(rest[10, 70]), 0.002)
        self.assertAlmostEqual(table[0]["source_model_Jy"], 0.050)
        self.assertAlmostEqual(table[0]["rest_model_Jy"], 0.002)


    def test_each_source_already_out_keeps_its_own_disc(self):
        root = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(root))
        plane = np.zeros((60, 80))
        plane[20, 30] = 0.050      # the source, taken with a disc of 5
        plane[45, 60] = 0.010      # one taken out earlier, with a disc of 12
        plane[45, 69] = 0.004      # 9 pixels from it: inside that disc
        plane[10, 70] = 0.002      # the rest
        write_model(root / "run_I-0000-model.fits", plane)

        _, table = TOOL.split_models(
            [root / "run_I-0000-model.fits"], root, 30, 20, 5,
            others=[(60, 45, 12)],
        )

        _, _, rest, _ = TOOL.read_fits(root / "rest-0000-model.fits")
        self.assertEqual(rest[45, 69], 0.0)
        self.assertAlmostEqual(float(rest[10, 70]), 0.002)
        self.assertAlmostEqual(table[0]["rest_model_Jy"], 0.002)


class FlowTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.products = self.root / "uGMRT_Products"
        self.products.mkdir()
        self.ms = self.root / "data" / "G1_bpa.ms"
        self.ms.mkdir(parents=True)
        (self.ms / "table.dat").write_text("visibilities")
        (self.root / "images").mkdir()
        for index in range(3):
            plane = np.zeros((60, 80))
            plane[20, 30] = 0.050
            plane[45, 60] = 0.010
            write_model(
                self.root / "images" / f"run_I-{index:04d}-model.fits", plane
            )
        self.calls = []
        self.directions = []
        self.casa_writes = True
        test = self

        def run_command(command, log_path):
            folder = Path(log_path).parent
            if "--casa-part" not in command:
                name = Path(command[command.index("-name") + 1])
                test.calls.append(f"predict {name.name}")
                test.assertTrue(Path(command[-1]).is_dir())
                # The model images named are there to predict.
                test.assertEqual(len(list(
                    name.parent.glob(name.name + "-*-model.fits"))), 3)
                return
            operation = command[command.index("--casa-part") + 1]
            test.calls.append(f"casa {operation}")
            test.assertTrue(
                Path(command[command.index("--field-ms") + 1]).is_dir()
            )
            if operation == "subtract-shift":
                Path(command[command.index("--subtracted-ms") + 1]).mkdir()
                test.directions.append(
                    command[command.index("--phasecenter") + 1]
                )
            if test.casa_writes:
                out = Path(command[command.index("--out-ms") + 1])
                out.mkdir()
                (out / "table.dat").write_text(operation)

        patches = [
            mock.patch.object(IMAGER, "run_command", run_command),
            mock.patch.object(
                IMAGER, "wsclean_capabilities",
                lambda prefix: {"version": (3, 1), "options": set()},
            ),
            mock.patch.object(TOOL, "load_imager", lambda: IMAGER),
            mock.patch.object(
                TOOL.shutil, "which", lambda name: f"/usr/bin/{name}"
            ),
        ]
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)

    def run_tool(self, *extra):
        with mock.patch("builtins.print"):
            return TOOL.main([
                "--ms", str(self.ms),
                "--models", str(self.root / "images" / "run_I"),
                "--pixel", "30", "20", "--radius", "5", "--name", "source1",
                "--workdir", str(self.root / "work"),
                "--outdir", str(self.products),
                "--channel-range", "51", "768", *extra,
            ])

    def test_the_source_is_left_alone_in_a_new_measurement_set(self):
        self.run_tool()

        folder = self.root / "work" / "source1"
        self.assertEqual(self.calls, ["predict rest", "casa subtract-shift"])
        self.assertTrue((folder / "source1.ms").is_dir())
        report = json.loads((folder / "extract.json").read_text())
        self.assertEqual(report["pixel"], [30, 20])
        self.assertEqual(report["measurement_set"],
                         str(folder / "source1.ms"))
        self.assertAlmostEqual(report["sub_bands"][0]["source_model_Jy"], 0.05)
        self.assertAlmostEqual(report["sub_bands"][2]["rest_model_Jy"], 0.01)
        # Ten pixels to the left, which is east, and ten south of the
        # centre, an arcsecond each.
        self.assertAlmostEqual(report["dec_deg"], -3.5 - 10 / 3600, places=6)
        self.assertAlmostEqual(
            (report["ra_deg"] - 33.0) * math.cos(math.radians(-3.5)),
            10 / 3600, places=6,
        )

    def test_the_data_given_are_only_copied(self):
        self.run_tool()

        self.assertEqual(
            sorted(path.name for path in self.ms.iterdir()), ["table.dat"]
        )
        self.assertEqual(
            sorted(path.name for path in self.ms.parent.iterdir()),
            ["G1_bpa.ms"],
        )
        self.assertEqual(list(self.products.iterdir()), [])

    def test_what_has_served_is_removed_and_the_claim_released(self):
        self.run_tool()

        folder = self.root / "work" / "source1"
        self.assertEqual(
            sorted(path.name for path in folder.iterdir()),
            ["extract.json", "source1.ms"],
        )

    def test_the_intermediate_products_can_be_kept(self):
        self.run_tool("--keep-intermediate")

        folder = self.root / "work" / "source1"
        self.assertTrue((folder / "field.ms").is_dir())
        self.assertTrue((folder / "rest-0000-model.fits").is_file())

    def test_an_earlier_result_is_not_overwritten_unasked(self):
        self.run_tool()
        self.calls.clear()

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("already exists", str(stopped.exception))
        self.assertEqual(self.calls, [])
        self.run_tool("--overwrite")
        self.assertEqual(self.calls, ["predict rest", "casa subtract-shift"])

    def test_a_folder_inside_the_products_is_refused(self):
        with self.assertRaises(SystemExit) as stopped, \
                mock.patch("builtins.print"):
            TOOL.main([
                "--ms", str(self.ms),
                "--models", str(self.root / "images" / "run_I"),
                "--pixel", "30", "20",
                "--workdir", str(self.products / "work"),
                "--outdir", str(self.products),
            ])

        self.assertIn("inside the pipeline's products", str(stopped.exception))
        self.assertEqual(self.calls, [])

    def test_a_casa_step_that_writes_nothing_is_reported(self):
        self.casa_writes = False

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("No Measurement Set was written", str(stopped.exception))
        self.assertFalse(
            (self.root / "work" / "source1" / IMAGER.LOCK_NAME).exists()
        )

    def test_a_terminal_without_singularity_is_said_so(self):
        with mock.patch.object(TOOL.shutil, "which", lambda name: None), \
                self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("not on PATH", str(stopped.exception))
        self.assertEqual(self.calls, [])

    def test_a_missing_measurement_set_is_said_at_once(self):
        self.ms = self.root / "data" / "absent.ms"

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("No Measurement Set", str(stopped.exception))

    def recipe(self, finished=True, original=None):
        """What selfcal_recipe.py leaves for the source it calibrated."""

        folder = self.root / "recipe_source1" / "G1"
        final = folder / "images" / "final" / "G1"
        final.mkdir(parents=True)
        for index in range(3):
            plane = np.zeros((40, 40))
            plane[20, 20] = 0.052
            write_model(final / f"wsclean_robust_0_I-{index:04d}-model.fits",
                        plane)
        calibrated = folder / "round5.ms"
        calibrated.mkdir()
        (calibrated / "table.dat").write_text("calibrated towards the source")
        record = {
            "done": ["copy", "final:image"] if finished else ["copy"],
            "calibrated_measurement_set": str(calibrated),
        }
        if original is not None:
            record["input"] = str(original)
        (folder / "recipe_report.json").write_text(json.dumps(record))
        return folder

    def peel(self, recipe, *extra):
        with mock.patch("builtins.print"):
            return TOOL.main([
                "--peel", str(recipe),
                "--models", str(self.root / "images" / "run_I"),
                "--pixel", "30", "20", "--radius", "5", "--name", "source1",
                "--workdir", str(self.root / "work"),
                "--outdir", str(self.products),
                "--channel-range", "51", "768", *extra,
            ])

    def test_a_calibrated_source_is_taken_out_and_the_field_given_back(self):
        recipe = self.recipe()

        self.peel(recipe)

        folder = self.root / "work" / "source1"
        # Its own model out of its calibrated data, back to the field's
        # centre, and the rest of the field's model added again.
        self.assertEqual(self.calls, [
            "predict wsclean_robust_0_I", "casa subtract-shift",
            "predict rest", "casa add",
        ])
        self.assertEqual(self.directions,
                         ["J2000 33.00000000deg -3.50000000deg"])
        self.assertEqual(
            (folder / "source1_peeled.ms" / "table.dat").read_text(), "add"
        )
        report = json.loads((folder / "peel.json").read_text())
        self.assertEqual(report["measurement_set"],
                         str(folder / "source1_peeled.ms"))
        self.assertEqual(report["field_centre_deg"], [33.0, -3.5])
        self.assertAlmostEqual(report["sub_bands"][0]["rest_model_Jy"], 0.01)

    def test_the_field_keeps_its_own_calibration(self):
        original = self.root / "work" / "source1" / "source1.ms"
        original.mkdir(parents=True)
        seen = []
        real = TOOL.casa_command

        def watching(*args, **kwargs):
            command = real(*args, **kwargs)
            seen.append(command)
            return command

        with mock.patch.object(TOOL, "casa_command", watching):
            self.peel(self.recipe(original=original))

        subtract = seen[0]
        self.assertEqual(subtract[subtract.index("--original-ms") + 1],
                         str(original))
        report = json.loads(
            (self.root / "work" / "source1" / "peel.json").read_text())
        self.assertEqual(report["calibration_kept"], "the field's")
        self.assertEqual(report["original"], str(original))

    def test_without_the_recipes_input_the_old_way_is_said(self):
        self.peel(self.recipe(original=self.root / "gone.ms"))

        report = json.loads(
            (self.root / "work" / "source1" / "peel.json").read_text())
        self.assertIn("the source's", report["calibration_kept"])
        self.assertIsNone(report["original"])

    def test_the_recipes_data_are_only_copied(self):
        recipe = self.recipe()

        self.peel(recipe)

        self.assertEqual(
            sorted(path.name for path in (recipe / "round5.ms").iterdir()),
            ["table.dat"],
        )
        self.assertEqual(
            sorted(path.name for path in (self.root / "work" / "source1")
                   .iterdir()),
            ["peel.json", "source1_peeled.ms"],
        )

    def test_sources_already_out_stay_out_of_what_is_added_back(self):
        recipe = self.recipe()
        seen = []
        real = TOOL.split_models

        def watching(images, folder, x, y, radius, others=()):
            seen.append(list(others))
            return real(images, folder, x, y, radius, others)

        with mock.patch.object(TOOL, "split_models", watching):
            self.peel(recipe, "--also-without", "60", "45")

        self.assertEqual(seen, [[(60, 45, 5)]])
        report = json.loads(
            (self.root / "work" / "source1" / "peel.json").read_text()
        )
        # With the disc it went out with: this one's, none being given.
        self.assertEqual(report["already_out"], [[60, 45, 5]])
        self.assertAlmostEqual(report["sub_bands"][0]["rest_model_Jy"], 0.0)

    def test_a_recipe_that_has_not_finished_is_not_peeled(self):
        recipe = self.recipe(finished=False)

        with self.assertRaises(SystemExit) as stopped:
            self.peel(recipe)

        self.assertIn("has not finished", str(stopped.exception))
        self.assertEqual(self.calls, [])

    def test_a_folder_that_is_not_a_recipes_is_said_so(self):
        with self.assertRaises(SystemExit) as stopped:
            self.peel(self.root / "nowhere")

        self.assertIn("No recipe report", str(stopped.exception))

    def test_models_with_different_sub_bands_are_refused(self):
        recipe = self.recipe()
        (recipe / "images" / "final" / "G1"
         / "wsclean_robust_0_I-0002-model.fits").unlink()

        with self.assertRaises(SystemExit) as stopped:
            self.peel(recipe)

        self.assertIn("have to be the same", str(stopped.exception))
        self.assertEqual(self.calls, [])

    def test_extracting_the_next_source_leaves_out_the_one_already_gone(self):
        self.run_tool("--also-without", "60", "45")

        report = json.loads(
            (self.root / "work" / "source1" / "extract.json").read_text()
        )
        # With the disc it went out with: this one's, none being given.
        self.assertEqual(report["already_out"], [[60, 45, 5]])
        # Nothing is left to subtract where the earlier source was.
        self.assertAlmostEqual(report["sub_bands"][0]["rest_model_Jy"], 0.0)


if __name__ == "__main__":
    unittest.main()
