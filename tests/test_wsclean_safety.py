import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_wsclean_module():
    path = SOURCE_DIR / "wsclean_imaging.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_wsclean_imaging",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


WSCLEAN = load_wsclean_module()


class WSCleanFieldSelectionTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.outdir = self.root / "Products"
        self.ms = self.root / "calibrated-target.ms"
        self.ms.mkdir()
        (self.ms / "table.dat").touch()
        inspection = self.outdir / "inspection"
        inspection.mkdir(parents=True)
        (inspection / "inspection_report.json").write_text(
            json.dumps(
                {
                    "measurement_set": str(self.ms),
                    "telescope": {"name": "lofar"},
                    "band": {"name": "LOFAR_HBA"},
                    "data_columns": ["DATA", "CORRECTED_DATA"],
                }
            )
        )

    def tearDown(self):
        self.temporary.cleanup()

    def write_fields(self, values):
        path = self.outdir / "inspection" / "fields.json"
        path.write_text(json.dumps(values))

    def manager(self, target, beam_mode="none"):
        return WSCLEAN.WSCleanImagingManager(
            self.ms,
            outdir=self.outdir,
            target=target,
            telescope="lofar",
            receiver_band="HBA",
            executable="wsclean",
            threads=2,
            profile="science",
            beam_mode=beam_mode,
            datacolumn="CORRECTED_DATA",
            dry_run=True,
        )

    def test_named_target_becomes_the_exact_numeric_field_selector(self):
        self.write_fields(
            [
                {"id": 0, "name": "calibrator"},
                {"id": 7, "name": "science_target"},
                {"id": 11, "name": "other_target"},
            ]
        )

        manager = self.manager("science_target")
        command = manager.build_command()

        self.assertEqual(manager.field_id, 7)
        self.assertEqual(command[command.index("-field") + 1], "7")
        self.assertNotIn("science_target", command)

    def test_unresolved_target_in_a_multifield_ms_fails_closed(self):
        self.write_fields(
            [
                {"id": 2, "name": "target_a"},
                {"id": 5, "name": "target_b"},
            ]
        )

        with self.assertRaisesRegex(
            RuntimeError,
            "could not be mapped.*multiple MS fields",
        ):
            self.manager("misspelled_target")

    def test_beam_mode_is_required_by_the_backend_too(self):
        self.write_fields([{"id": 0, "name": "science_target"}])

        with self.assertRaisesRegex(ValueError, "beam state must be explicit"):
            self.manager("science_target", beam_mode=None)

    def test_reusable_report_rejects_changed_fits_content(self):
        self.write_fields([{"id": 0, "name": "science_target"}])
        manager = self.manager("science_target")
        restored = manager.fitsdir / "science_target-MFS-image.fits"
        restored.write_bytes(b"original image")
        identity = WSCLEAN.IMAGING_SAFETY.file_content_identity(restored)
        manager.report_path.write_text(json.dumps({
            "imaging_backend": "wsclean",
            "provenance_fingerprint": manager.provenance_fingerprint,
            "measurement_set": manager.ms,
            "measurement_set_stat_sha256": (
                WSCLEAN.IMAGING_SAFETY.path_stat_fingerprint(manager.ms)
            ),
            "telescope": manager.telescope,
            "band": manager.band,
            "datacolumn": manager.datacolumn,
            "target": manager.target,
            "execution_profile": manager.profile,
            "imaging_parameters": {
                "beam_mode": manager.beam_mode,
                "command": manager.build_command(),
            },
            "products": {
                "restored_image": str(restored),
                "pbcor_image": None,
                "fits": [str(restored)],
                "fits_content_identities": {
                    str(restored.resolve()): identity,
                },
            },
        }))

        self.assertIsNotNone(manager.reusable_report())
        restored.write_bytes(b"replacement image")
        self.assertIsNone(manager.reusable_report())


if __name__ == "__main__":
    unittest.main()
