import importlib.util
import os
import tempfile
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_imaging_safety_module():
    path = SOURCE_DIR / "imaging_safety.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_imaging_safety",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SAFETY = load_imaging_safety_module()


class ImagingSafetyTests(unittest.TestCase):

    @staticmethod
    def converter(quantity, output_unit):
        scales = {
            "Hz": 1.0,
            "kHz": 1e3,
            "MHz": 1e6,
            "GHz": 1e9,
        }
        if output_unit != "Hz":
            raise AssertionError(output_unit)
        return {
            "value": quantity["value"] * scales[quantity["unit"]],
            "unit": output_unit,
        }

    def test_frequency_is_normalized_from_the_native_casa_axis_unit(self):
        self.assertEqual(
            SAFETY.normalize_frequency_hz(
                1.284,
                "GHz",
                self.converter,
            ),
            1.284e9,
        )

    def test_geometry_comparison_uses_live_grid_not_record_hash(self):
        source = {
            "shape": [4096, 4096, 1, 1],
            "pixel_axes": [0, 1],
            "direction_increment_rad": [-1.0e-5, 1.0e-5],
            "direction_reference_pixel": [2047.0, 2047.0],
            "direction_reference_value_rad": [6.2, -0.3],
            "frequency_hz": 650.0e6,
            "coordsys_sha256": "source-serialized-record",
        }
        casa_normalized = dict(source)
        casa_normalized.update({
            "direction_reference_value_rad": [
                6.2 - 2.0 * 3.141592653589793,
                -0.3 + 1.0e-13,
            ],
            "frequency_hz": 650.0e6 + 1.0e-4,
            "coordsys_sha256": "casa-normalized-record",
        })

        self.assertTrue(
            SAFETY.image_geometries_match(source, casa_normalized)
        )
        self.assertIsNone(
            SAFETY.image_geometry_mismatch(source, casa_normalized)
        )

    def test_geometry_comparison_rejects_real_pixel_grid_changes(self):
        source = {
            "shape": [4, 4, 1, 1],
            "pixel_axes": [0, 1],
            "direction_increment_rad": [-1.0e-5, 1.0e-5],
            "direction_reference_pixel": [2.0, 2.0],
            "direction_reference_value_rad": [0.2, -0.3],
            "frequency_hz": 650.0e6,
            "coordsys_sha256": "one",
        }
        shifted = dict(source)
        shifted["direction_reference_pixel"] = [2.01, 2.0]

        reason = SAFETY.image_geometry_mismatch(source, shifted)
        self.assertIn("reference pixels differ", reason)
        self.assertFalse(SAFETY.image_geometries_match(source, shifted))

    def test_invalid_frequency_or_converter_result_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "greater than zero"):
            SAFETY.normalize_frequency_hz(0.0, "Hz", self.converter)

        with self.assertRaisesRegex(ValueError, "quantity value"):
            SAFETY.normalize_frequency_hz(
                1.0,
                "GHz",
                lambda quantity, output_unit: None,
            )

    def test_force_cleanup_patterns_are_scoped_to_the_image_root(self):
        patterns = SAFETY.image_root_partial_patterns(
            "/products/images/target_final_continuum"
        )

        self.assertEqual(
            patterns,
            (
                ".target_final_continuum*.partial.*",
                ".target_final_continuum*.tmp.*",
            ),
        )

    def test_pbcor_dependency_order_includes_every_derived_checkpoint(self):
        self.assertEqual(
            SAFETY.PBCOR_DEPENDENT_CHECKPOINTS,
            (
                "pbcor",
                "common_beam",
                "metrics",
                "qa_plots",
                "fits",
                "imaging",
            ),
        )

    def test_ugmrt_primary_beam_uses_the_ncRA_2023_polynomial(self):
        coefficients = [-3.129, 38.816, -21.608, 4.483]

        self.assertEqual(
            SAFETY.ugmrt_primary_beam_response(0.0, coefficients),
            1.0,
        )
        self.assertAlmostEqual(
            SAFETY.ugmrt_primary_beam_response(30.0, coefficients),
            0.04690243,
            places=8,
        )
        with self.assertRaisesRegex(ValueError, "four coefficients"):
            SAFETY.ugmrt_primary_beam_response(1.0, [1.0, 2.0])

    def test_ugmrt_beam_cutoff_uses_first_main_lobe_crossing(self):
        coefficients = [-3.263, 42.618, -25.580, 5.823]
        cutoff = SAFETY.first_outward_beam_cutoff(coefficients, 0.3)

        self.assertGreater(cutoff, 18.0)
        self.assertLess(cutoff, 18.5)
        self.assertAlmostEqual(
            SAFETY.ugmrt_primary_beam_response(cutoff, coefficients),
            0.3,
            places=8,
        )
        # The unconstrained fit rises again far outside the main lobe; that
        # disconnected response must remain beyond the chosen cutoff.
        self.assertGreater(
            SAFETY.ugmrt_primary_beam_response(45.0, coefficients),
            0.3,
        )

    def test_casa_table_fingerprint_detects_durable_changes_but_ignores_locks(self):
        with tempfile.TemporaryDirectory() as temporary:
            table = Path(temporary) / "image.table"
            table.mkdir()
            (table / "table.dat").write_text("original")
            first = SAFETY.path_stat_fingerprint(table)

            (table / "table.lock").write_text("transient")
            self.assertEqual(SAFETY.path_stat_fingerprint(table), first)

            (table / "table.f0").write_text("new pixels")
            self.assertNotEqual(SAFETY.path_stat_fingerprint(table), first)


class ContentIdentityTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.image = Path(self.temporary.name) / "target.ugmrt.pb.tt0"
        self.image.mkdir()
        (self.image / "table.dat").write_bytes(b"header")
        (self.image / "table.f0").write_bytes(b"pixels-v1")
        (self.image / "table.lock").write_bytes(b"lock")
        (self.image / "logtable").mkdir()
        (self.image / "logtable" / "table.dat").write_bytes(b"history")

    def tearDown(self):
        self.temporary.cleanup()

    def shift_timestamps(self, seconds):
        for item in self.image.rglob("*"):
            if item.is_file():
                stat = item.stat()
                os.utime(
                    item,
                    ns=(stat.st_atime_ns, stat.st_mtime_ns + seconds * 10**9),
                )

    def provenance(self):
        return {
            "pb_stat_sha256": SAFETY.path_stat_fingerprint(self.image),
            "pb_content_identity": SAFETY.casa_table_content_identity(
                self.image
            ),
        }

    def test_new_timestamps_do_not_change_a_content_identity(self):
        record = self.provenance()
        self.shift_timestamps(7)

        self.assertNotEqual(
            SAFETY.path_stat_fingerprint(self.image),
            record["pb_stat_sha256"],
        )
        self.assertEqual(
            SAFETY.artifact_matches_provenance(record, "pb", self.image),
            (True, ""),
        )

    def test_metadata_rewrites_are_not_pixel_changes(self):
        record = self.provenance()
        (self.image / "table.dat").write_bytes(b"header rewritten by casa")
        (self.image / "logtable" / "table.dat").write_bytes(b"one more open")

        self.assertEqual(
            SAFETY.artifact_matches_provenance(record, "pb", self.image),
            (True, ""),
        )

    def test_a_table_without_data_files_is_compared_whole(self):
        bare = self.image.parent / "bare.tt0"
        bare.mkdir()
        (bare / "table.dat").write_bytes(b"only metadata")
        record = {
            "pb_content_identity": SAFETY.casa_table_content_identity(bare)
        }
        self.assertTrue(
            SAFETY.artifact_matches_provenance(record, "pb", bare)[0]
        )

        (bare / "table.dat").write_bytes(b"metadata rewritten")
        matches, detail = SAFETY.artifact_matches_provenance(
            record,
            "pb",
            bare,
        )
        self.assertFalse(matches)
        self.assertIn("rewritten table.dat", detail)

    def test_rewritten_or_added_files_are_named(self):
        record = self.provenance()
        (self.image / "table.f0").write_bytes(b"pixels-v2")
        (self.image / "table.f1").write_bytes(b"extra")

        matches, detail = SAFETY.artifact_matches_provenance(
            record,
            "pb",
            self.image,
        )

        self.assertFalse(matches)
        self.assertIn("rewritten table.f0", detail)
        self.assertIn("added table.f1", detail)

    def test_lock_files_are_not_content(self):
        record = self.provenance()
        (self.image / "table.lock").write_bytes(b"another reader")
        self.assertTrue(
            SAFETY.artifact_matches_provenance(record, "pb", self.image)[0]
        )

    def test_older_stat_only_provenance_is_still_compared_by_stat(self):
        record = {"pb_stat_sha256": SAFETY.path_stat_fingerprint(self.image)}
        self.assertTrue(
            SAFETY.artifact_matches_provenance(record, "pb", self.image)[0]
        )

        self.shift_timestamps(7)
        matches, detail = SAFETY.artifact_matches_provenance(
            record,
            "pb",
            self.image,
        )
        self.assertFalse(matches)
        self.assertIn("before content identities", detail)

    def test_fits_identity_ignores_modification_time_only(self):
        fits = Path(self.temporary.name) / "target.fits"
        fits.write_bytes(b"SIMPLE")
        recorded = {str(fits.resolve()): SAFETY.file_content_identity(fits)}

        stat = fits.stat()
        os.utime(fits, ns=(stat.st_atime_ns, stat.st_mtime_ns + 5 * 10**9))
        self.assertTrue(
            SAFETY.validate_file_content_identities([fits], recorded)[0]
        )

        fits.write_bytes(b"SIMPLF")
        self.assertFalse(
            SAFETY.validate_file_content_identities([fits], recorded)[0]
        )


if __name__ == "__main__":
    unittest.main()
