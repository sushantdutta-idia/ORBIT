import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_imaging_module_without_casa():
    casatasks = types.ModuleType("casatasks")

    def unused_task(**kwargs):
        return None

    for name in [
        "exportfits",
        "imhead",
        "impbcor",
        "imsmooth",
        "imstat",
        "rmtables",
        "tclean",
        "widebandpbcor",
    ]:
        setattr(casatasks, name, unused_task)

    casatools = types.ModuleType("casatools")

    class UnusedTool:
        pass

    casatools.image = UnusedTool
    casatools.msmetadata = UnusedTool
    casatools.quanta = UnusedTool
    casatools.table = UnusedTool

    path = SOURCE_DIR / "imaging.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_imaging_atomic_safety",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    with mock.patch.dict(
        sys.modules,
        {
            "casatasks": casatasks,
            "casatools": casatools,
        },
    ):
        spec.loader.exec_module(module)
    return module


IMAGING = load_imaging_module_without_casa()


def make_casa_image(path, marker=None):
    path.mkdir(parents=True, exist_ok=True)
    (path / "table.dat").touch()
    if marker:
        (path / marker).touch()
    return path


class ImagingAtomicSafetyTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def manager(self):
        manager = object.__new__(IMAGING.FinalImagingManager)
        manager.pbmin = 0.2
        manager.stale_pbcor_products = []
        return manager

    def test_atomic_impbcor_keeps_old_output_until_new_one_validates(self):
        manager = self.manager()
        restored = make_casa_image(self.root / "restored.image")
        beam = make_casa_image(self.root / "beam.image")
        corrected = make_casa_image(
            self.root / "restored.image.pbcor",
            marker="old-generation",
        )

        def successful_impbcor(**kwargs):
            self.assertTrue((corrected / "old-generation").exists())
            make_casa_image(Path(kwargs["outfile"]), marker="new-generation")

        with mock.patch.object(IMAGING, "impbcor", successful_impbcor):
            published = manager.call_impbcor_atomic(
                restored,
                beam,
                corrected,
            )

        self.assertEqual(Path(published), corrected.resolve())
        self.assertTrue((corrected / "new-generation").exists())
        archived = list(self.root.glob("restored.image.pbcor.stale.*"))
        self.assertEqual(len(archived), 1)
        self.assertTrue((archived[0] / "old-generation").exists())

    def test_failed_impbcor_leaves_old_output_at_canonical_path(self):
        manager = self.manager()
        restored = make_casa_image(self.root / "restored.image")
        beam = make_casa_image(self.root / "beam.image")
        corrected = make_casa_image(
            self.root / "restored.image.pbcor",
            marker="old-generation",
        )

        def failed_impbcor(**kwargs):
            self.assertTrue((corrected / "old-generation").exists())
            raise RuntimeError("simulated CASA failure")

        with mock.patch.object(IMAGING, "impbcor", failed_impbcor):
            with self.assertRaisesRegex(RuntimeError, "simulated CASA failure"):
                manager.call_impbcor_atomic(restored, beam, corrected)

        self.assertTrue((corrected / "old-generation").exists())
        self.assertEqual(
            list(self.root.glob("restored.image.pbcor.stale.*")),
            [],
        )

    def test_uncommitted_beam_is_preserved_before_retry(self):
        manager = self.manager()
        uncommitted = make_casa_image(
            self.root / "target.ugmrt.pb.tt0",
            marker="old-uncommitted-generation",
        )

        archived = manager.archive_incomplete_casa_product(
            uncommitted,
            "uGMRT primary-beam image",
        )

        self.assertFalse(uncommitted.exists())
        archived_path = Path(archived)
        self.assertTrue(archived_path.is_dir())
        self.assertTrue(
            (archived_path / "old-uncommitted-generation").exists()
        )
        self.assertIn(str(archived_path), manager.stale_pbcor_products)

    def test_live_geometry_normalizes_frequency_and_retains_casa_record(self):
        manager = self.manager()
        image_path = make_casa_image(self.root / "restored.image")
        coordinate_record = {
            "direction0": {
                "crpix": IMAGING.np.asarray([2.0, 2.0]),
            },
            "spectral1": {
                "wcs": {"crval": 1.284, "unit": "GHz"},
            },
        }

        class FakeCoords:
            def findcoordinate(self, kind):
                self_kind = kind
                if self_kind != "direction":
                    raise AssertionError(self_kind)
                return {"pixel": [0, 1]}

            def increment(self, **kwargs):
                return {"numeric": IMAGING.np.asarray([1e-5, 1e-5])}

            def referencepixel(self, **kwargs):
                return {"numeric": IMAGING.np.asarray([2.0, 2.0])}

            def referencevalue(self, type=None, **kwargs):
                if type == "direction":
                    return {"numeric": IMAGING.np.asarray([0.1, -0.2])}
                if type == "spectral":
                    return {"numeric": IMAGING.np.asarray([1.284])}
                raise AssertionError(type)

            def units(self, **kwargs):
                return "GHz"

            def torecord(self):
                return coordinate_record

            def done(self):
                return None

        class FakeImageTool:
            def open(self, path):
                self.path = path

            def shape(self):
                return [4, 4, 1, 1]

            def coordsys(self):
                return FakeCoords()

            def close(self):
                return None

        class FakeQuanta:
            def convert(self, quantity, output_unit):
                self_output_unit = output_unit
                if quantity["unit"] != "GHz" or self_output_unit != "Hz":
                    raise AssertionError((quantity, self_output_unit))
                return {"value": quantity["value"] * 1e9, "unit": "Hz"}

        with (
            mock.patch.object(IMAGING, "image_tool", FakeImageTool),
            mock.patch.object(IMAGING, "quanta_tool", FakeQuanta),
        ):
            descriptor = manager.image_geometry_descriptor(
                image_path,
                include_record=True,
            )

        self.assertEqual(descriptor["geometry"]["frequency_hz"], 1.284e9)
        self.assertIs(descriptor["coordinate_record"], coordinate_record)

    def test_ugmrt_pb_accepts_casa_normalized_equivalent_geometry(self):
        manager = self.manager()
        source = make_casa_image(self.root / "target.image.tt0")
        manager.image_tt0 = str(source)
        manager.imagename = str(self.root / "target")
        manager.band = "B4"
        manager.config = {"imsize": 4, "cell": "1arcsec"}
        manager.ugmrt_pb_coefficients = [-3.263, 42.618, -25.580, 5.823]
        manager.ugmrt_pb_valid_frequency_mhz = [550.0, 850.0]
        source_geometry = {
            "shape": [4, 4, 1, 1],
            "pixel_axes": [0, 1],
            "direction_increment_rad": [-1.0e-5, 1.0e-5],
            "direction_reference_pixel": [2.0, 2.0],
            "direction_reference_value_rad": [0.2, -0.3],
            "frequency_hz": 650.0e6,
            "coordsys_sha256": "source-record",
        }
        normalized_geometry = dict(source_geometry)
        normalized_geometry["coordsys_sha256"] = "casa-normalized-record"

        def descriptor(path, include_record=False):
            resolved = Path(path).resolve()
            if resolved == source.resolve():
                value = {
                    "geometry": source_geometry,
                    "stat_sha256": "source-stat",
                }
                if include_record:
                    value["coordinate_record"] = {"coords": "source"}
                return value
            return {
                "geometry": normalized_geometry,
                "stat_sha256": "beam-stat",
            }

        class FakeOutput:
            def fromarray(self, outfile, **kwargs):
                make_casa_image(Path(outfile))

            def done(self):
                return None

        manager.image_geometry_descriptor = descriptor
        with mock.patch.object(IMAGING, "image_tool", FakeOutput):
            beam = manager.build_ugmrt_primary_beam_image()

        beam_path = Path(beam)
        self.assertTrue((beam_path / "table.dat").is_file())
        provenance = IMAGING.load_json(
            Path(str(beam_path) + ".provenance.json"),
            {},
        )
        self.assertEqual(provenance["pb_geometry"], normalized_geometry)
        self.assertEqual(provenance["source_geometry"], source_geometry)

    def build_ugmrt_beam_fixture(self):
        manager = self.manager()
        source = make_casa_image(self.root / "target.image.tt0")
        (source / "table.f0").write_bytes(b"restored pixels")
        manager.image_tt0 = str(source)
        manager.imagename = str(self.root / "target")
        manager.band = "B4"
        manager.config = {"imsize": 4, "cell": "1arcsec"}
        manager.ugmrt_pb_coefficients = [-3.263, 42.618, -25.580, 5.823]
        manager.ugmrt_pb_valid_frequency_mhz = [550.0, 850.0]
        geometry = {
            "shape": [4, 4, 1, 1],
            "pixel_axes": [0, 1],
            "direction_increment_rad": [-1.0e-5, 1.0e-5],
            "direction_reference_pixel": [2.0, 2.0],
            "direction_reference_value_rad": [0.2, -0.3],
            "frequency_hz": 650.0e6,
            "coordsys_sha256": "record",
        }

        def descriptor(path, include_record=False):
            value = {
                "geometry": geometry,
                "stat_sha256": IMAGING.IMAGING_SAFETY.path_stat_fingerprint(
                    path
                ),
            }
            if include_record:
                value["coordinate_record"] = {"coords": "source"}
            return value

        class FakeOutput:
            generation = 0

            def fromarray(self, outfile, **kwargs):
                FakeOutput.generation += 1
                image = make_casa_image(Path(outfile))
                (image / "table.f0").write_bytes(b"beam pixels")

            def done(self):
                return None

        manager.image_geometry_descriptor = descriptor
        return manager, FakeOutput

    def test_beam_with_new_timestamps_is_reused_by_content(self):
        manager, output = self.build_ugmrt_beam_fixture()
        with mock.patch.object(IMAGING, "image_tool", output):
            beam = Path(manager.build_ugmrt_primary_beam_image())
            data = beam / "table.f0"
            stat = data.stat()
            IMAGING.os.utime(
                data,
                ns=(stat.st_atime_ns, stat.st_mtime_ns + 9 * 10**9),
            )
            again = manager.build_ugmrt_primary_beam_image()

        self.assertEqual(Path(again), beam)
        self.assertEqual(output.generation, 1)
        self.assertEqual(manager.stale_pbcor_products, [])

    def test_beam_with_a_stale_stat_only_record_is_rebuilt(self):
        manager, output = self.build_ugmrt_beam_fixture()
        with mock.patch.object(IMAGING, "image_tool", output):
            beam = Path(manager.build_ugmrt_primary_beam_image())
            provenance_path = Path(str(beam) + ".provenance.json")
            legacy = IMAGING.load_json(provenance_path, {})
            legacy.pop("pb_content_identity")
            legacy.pop("source_content_identity")
            legacy["pb_stat_sha256"] = "recorded-before-timestamps-moved"
            manager.save_json(legacy, provenance_path)

            rebuilt = manager.build_ugmrt_primary_beam_image()

        self.assertEqual(Path(rebuilt), beam)
        self.assertEqual(output.generation, 2)
        self.assertEqual(len(manager.stale_pbcor_products), 1)
        self.assertTrue(Path(manager.stale_pbcor_products[0]).is_dir())
        refreshed = IMAGING.load_json(provenance_path, {})
        self.assertIn("pb_content_identity", refreshed)
        self.assertIn("source_content_identity", refreshed)

    def test_katbeam_branch_never_reuses_a_corrected_image(self):
        manager = self.manager()
        manager.primary_beam_provider = "katbeam"
        manager.checkdir = self.root / "checkpoints"
        manager.checkdir.mkdir()
        calls = []
        manager.run_katbeam_correction = lambda: calls.append("corrected")

        manager.run_primary_beam_correction()

        self.assertEqual(calls, ["corrected"])
        self.assertTrue((manager.checkdir / "pbcor.done").exists())

    def test_ugmrt_pbcor_provenance_uses_live_pbmin_and_is_force_scoped(self):
        manager = self.manager()
        manager.imagename = str(self.root / "target_final_continuum")
        manager.image_tt0 = str(
            make_casa_image(self.root / "target_final_continuum.image.tt0")
        )
        beam = make_casa_image(
            self.root / "target_final_continuum.ugmrt.pb.tt0"
        )
        corrected = make_casa_image(
            self.root / "target_final_continuum.image.tt0.pbcor"
        )
        manager.band = "B3"
        manager.config = {}  # pbmin is deliberately not stored here.
        manager.ugmrt_pb_coefficients = [-3.129, 38.816, -21.608, 4.483]
        manager.ugmrt_pb_valid_frequency_mhz = [250.0, 500.0]
        manager.pbcor_status = {}
        manager.build_ugmrt_primary_beam_image = lambda: str(beam)
        manager.clear_pbcor_dependent_checkpoints = lambda: None
        manager.call_impbcor_atomic = lambda *args: str(corrected)
        manager.find_pbcor_image = lambda: str(corrected)

        manager.run_ugmrt_primary_beam_correction()

        provenance = IMAGING.load_json(
            Path(str(corrected) + ".provenance.json"),
            {},
        )
        self.assertEqual(provenance["pbmin"], 0.2)
        self.assertEqual(
            manager.pbcor_method,
            "ugmrt_polynomial_impbcor_tt0",
        )
        candidates = set(manager.product_candidates())
        self.assertIn(beam, candidates)
        self.assertIn(Path(str(beam) + ".provenance.json"), candidates)

    def test_ugmrt_subband_uses_measured_polynomial_not_tclean_pb(self):
        manager = self.manager()
        manager.make_subbands = True
        manager.primary_beam_provider = "ugmrt_polynomial"
        manager.band = "B3"
        manager.target = "target"
        manager.reportdir = self.root / "reports"
        manager.subbanddir = self.root / "subbands"
        manager.checkdir = self.root / "checkpoints"
        for directory in [
            manager.reportdir,
            manager.subbanddir,
            manager.checkdir,
        ]:
            directory.mkdir()
        manager.require_pbcor = True
        manager.clean_diagnostics = {}
        manager.stale_pbcor_workdirectories = []
        manager.make_subband_plan = lambda: [{
            "index": 1,
            "name": "subband_01",
            "spw": "0:0~3",
            "nchan": 4,
            "freq_min_MHz": 300.0,
            "freq_max_MHz": 320.0,
            "centre_MHz": 310.0,
        }]
        root = manager.subbanddir / "target_subband_01"
        image = make_casa_image(Path(str(root) + ".image.tt0"))
        beam = make_casa_image(Path(str(root) + ".ugmrt.pb.tt0"))
        corrected = make_casa_image(Path(str(image) + ".pbcor"))
        manager.find_subband_image = lambda imagename: str(image)
        manager.find_subband_pbcor = lambda imagename: str(corrected)
        calls = []

        def measured_pb(source_image=None, image_root=None):
            calls.append((source_image, image_root))
            return str(beam)

        manager.build_ugmrt_primary_beam_image = measured_pb
        manager.call_impbcor_atomic = lambda *args: str(corrected)
        manager.deconvolution_diagnostics = lambda imagename: {}

        manager.run_subband_imaging()

        self.assertEqual(calls, [(str(image), str(root))])
        product = manager.subband_products[0]
        self.assertEqual(product["pb_image"], str(beam))
        self.assertEqual(
            product["pbcor_method"],
            "ugmrt_polynomial_impbcor",
        )
        provenance = IMAGING.load_json(
            Path(str(corrected) + ".provenance.json"),
            {},
        )
        self.assertEqual(provenance["provider"], "ugmrt_polynomial")
        self.assertEqual(provenance["pbmin"], 0.2)

    def test_replacing_pbcor_clears_only_dependent_checkpoints(self):
        manager = self.manager()
        manager.checkdir = self.root / "checkpoints"
        manager.checkdir.mkdir()
        for stage in [
            *IMAGING.IMAGING_SAFETY.PBCOR_DEPENDENT_CHECKPOINTS,
            "tclean",
            "subband_imaging",
        ]:
            (manager.checkdir / f"{stage}.done").touch()

        manager.clear_pbcor_dependent_checkpoints()

        for stage in IMAGING.IMAGING_SAFETY.PBCOR_DEPENDENT_CHECKPOINTS:
            self.assertFalse((manager.checkdir / f"{stage}.done").exists())
        self.assertTrue((manager.checkdir / "tclean.done").exists())
        self.assertTrue((manager.checkdir / "subband_imaging.done").exists())

    def test_force_cleanup_removes_only_image_root_hidden_partials(self):
        manager = self.manager()
        manager.imagename = str(self.root / "target_final_continuum")
        manager.product_candidates = lambda: []
        manager.subband_product_candidates = lambda: []
        scoped_dir = self.root / (
            ".target_final_continuum.katbeam.pb.tt0.partial.17"
        )
        scoped_dir.mkdir()
        scoped_file = self.root / (
            ".target_final_continuum.katbeam.pb.tt0.provenance.json.tmp.17"
        )
        scoped_file.touch()
        unrelated = self.root / ".another_target.image.tt0.pbcor.partial.17"
        unrelated.mkdir()

        manager.remove_existing_products()

        self.assertFalse(scoped_dir.exists())
        self.assertFalse(scoped_file.exists())
        self.assertTrue(unrelated.exists())


if __name__ == "__main__":
    unittest.main()
