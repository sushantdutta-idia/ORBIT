import importlib.util
import tempfile
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_source_detection_module():
    path = SOURCE_DIR / "source_detection.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_source_detection",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SOURCE_DETECTION = load_source_detection_module()


class SourceDetectionImageSelectionTests(unittest.TestCase):

    def test_pb_preference_rejects_report_without_pb_product(self):
        with tempfile.TemporaryDirectory() as temporary:
            outdir = Path(temporary) / "Products"
            fitsdir = outdir / "imaging" / "fits"
            fitsdir.mkdir(parents=True)
            restored = fitsdir / "target.image.tt0.fits"
            restored.touch()
            report = {
                "products": {
                    "restored_image": str(restored),
                    "pbcor_image": None,
                    # A real restored FITS candidate must not be used as a
                    # silent fallback for PB-preferred science cataloguing.
                    "fits": [str(restored)],
                }
            }

            with self.assertRaisesRegex(
                RuntimeError,
                "records no PB-corrected product",
            ):
                SOURCE_DETECTION.select_input_image(
                    None,
                    report,
                    outdir,
                    "pbcor",
                )

    def test_restored_preference_can_select_the_recorded_restored_fits(self):
        with tempfile.TemporaryDirectory() as temporary:
            outdir = Path(temporary) / "Products"
            restored = outdir / "target-restored.fits"
            restored.parent.mkdir(parents=True)
            restored.touch()
            report = {
                "products": {
                    "restored_image": str(restored),
                    "pbcor_image": None,
                    "fits": [str(restored)],
                }
            }

            selected = SOURCE_DETECTION.select_input_image(
                None,
                report,
                outdir,
                "restored",
            )
            self.assertEqual(selected, restored.resolve())


if __name__ == "__main__":
    unittest.main()
