import csv
import json
import tempfile
import unittest
from pathlib import Path

import importlib.util

SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_ddcal():
    spec = importlib.util.spec_from_file_location(
        "test_ddcal_workflow_module", SOURCE_DIR / "ddcal.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ddcal = load_ddcal()

try:
    import numpy as np
except ImportError:
    np = None


class DDCalWorkflowTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.outdir = self.root / "Products"
        self.ms = self.root / "target.ms"
        self.ms.mkdir()
        # A valid CASA table has durable metadata even in this CASA-free
        # workflow fixture.  DD checkpoint provenance hashes the materialized
        # table inventory, so keep the fixture structurally representative.
        (self.ms / "table.dat").write_bytes(b"mock-main-table")
        self.image = self.root / "target.pbcor.fits"
        self.image.touch()

        inspection = {
            "band": "B4",
            "target_fields": ["target"],
            "spws": {
                "0": {
                    "fmin_MHz": 550.0,
                    "fmax_MHz": 750.0,
                }
            },
        }
        imaging = {
            "measurement_set": str(self.ms),
            "datacolumn": "CORRECTED_DATA",
            "target": "target",
            "band": "B4",
            "imaging_parameters": {
                "imsize": 4096,
                "cell": "0.8arcsec",
                "robust": 0.5,
                "pbmin": 0.3,
            },
            "products": {
                "fits": [str(self.image)],
                "casa_images": [],
            },
        }

        self.write_json(
            self.outdir / "inspection" / "inspection_report.json",
            inspection,
        )
        self.write_json(
            self.outdir / "imaging" / "reports" / "imaging_report.json",
            imaging,
        )

        self.catalog = (
            self.outdir /
            "source_detection" /
            "catalogs" /
            "target_sources.csv"
        )
        self.catalog.parent.mkdir(parents=True)
        with open(self.catalog, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "Source_id",
                    "RA",
                    "DEC",
                    "Total_flux",
                    "Peak_flux",
                ],
            )
            writer.writeheader()
            for index in range(5):
                writer.writerow(
                    {
                        "Source_id": f"S{index + 1}",
                        "RA": 150.0 + 0.2 * index,
                        "DEC": -30.0 + 0.2 * index,
                        "Total_flux": 0.5 - 0.05 * index,
                        "Peak_flux": 0.4 - 0.04 * index,
                    }
                )

        self.write_json(
            self.outdir /
            "source_detection" /
            "reports" /
            "source_detection_report.json",
            {
                "status": "success",
                "products": {
                    "source_catalog_csv": str(self.catalog),
                },
            },
        )

    def tearDown(self):
        self.temporary.cleanup()

    @staticmethod
    def write_json(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))

    def manager(self, mode="auto"):
        return ddcal.DDCalManager(
            ms=str(self.ms),
            outdir=str(self.outdir),
            mode=mode,
            engine="auto",
            execute=False,
            strict=True,
            ncpu=8,
            skip_imaging=False,
        )

    def test_quantity_parsing(self):
        self.assertEqual(ddcal.parse_quantity("5min", "min"), 5.0)
        self.assertEqual(ddcal.parse_quantity("16MHz", "mhz"), 16.0)
        self.assertEqual(ddcal.parse_quantity("0.8arcsec", "arcsec"), 0.8)

    def test_catalogue_constructor_and_units(self):
        sources = ddcal.read_source_catalogue(self.catalog)
        self.assertEqual(len(sources), 5)
        self.assertAlmostEqual(sources[0].flux_jy, 0.5)
        self.assertAlmostEqual(sources[0].dec_deg, -30.0)

    def test_auto_survey_prepares_documented_solver_sequence(self):
        manager = self.manager(mode="auto")
        manager.run()

        report = json.loads(
            (
                self.outdir /
                "ddcal" /
                "reports" /
                "ddcal_report.json"
            ).read_text()
        )
        commands = report["commands"]

        self.assertEqual(report["status"], "prepared")
        self.assertEqual(report["mode"], "survey")
        self.assertEqual(report["engine"], "killms-ddfacet")
        self.assertEqual(report["solver_direction_count"], 5)
        self.assertEqual(
            [item["name"] for item in commands],
            [
                "ddfacet_initial",
                "make_model_clusters",
                "killms_solve",
                "ddfacet_dd_image",
            ],
        )
        self.assertIn(
            "--DDESolutions-DDSols=",
            commands[-1]["command"],
        )
        self.assertIn(
            "--Predict-InitDicoModel=",
            commands[-1]["command"],
        )
        self.assertIn(
            "--NChanSols 13",
            commands[2]["command"],
        )
        pb_script, _pb_fits, _metrics = (
            manager.write_dd_pbcor_script(
                self.image,
                self.image,
            )
        )
        compile(
            pb_script.read_text(),
            str(pb_script),
            "exec",
        )

    def test_single_mode_limits_solver_clusters(self):
        manager = self.manager(mode="single")
        manager.run()

        report = json.loads(
            (
                self.outdir /
                "ddcal" /
                "reports" /
                "ddcal_report.json"
            ).read_text()
        )

        self.assertEqual(report["mode"], "single")
        self.assertEqual(report["solver_direction_count"], 3)
        make_model = next(
            item
            for item in report["commands"]
            if item["name"] == "make_model_clusters"
        )
        self.assertIn("--NCluster 0", make_model["command"])
        self.assertIn("--ds9PreClusterFile", make_model["command"])
        solver_regions = (
            self.outdir /
            "ddcal" /
            "directions" /
            "dd_solver_directions.reg"
        ).read_text()
        self.assertEqual(solver_regions.count("circle("), 3)

    @unittest.skipIf(np is None, "NumPy is not installed in this interpreter")
    def test_killms_solution_validation_rejects_unusable_gains(self):
        manager = self.manager()
        valid = self.root / "valid_solutions.npz"
        invalid = self.root / "invalid_solutions.npz"
        np.savez(valid, G=np.ones((2, 3, 4), dtype=np.complex128))
        np.savez(
            invalid,
            G=np.full((2, 3, 4), np.nan + 1j * np.nan),
        )

        qa = manager.validate_killms_solutions(valid)
        self.assertEqual(qa["usable_fraction"], 1.0)
        self.assertTrue(
            (
                self.outdir /
                "ddcal" /
                "solutions" /
                "killms_solution_qa.json"
            ).is_file()
        )

        with self.assertRaisesRegex(
            RuntimeError,
            "usable finite/non-zero solution fraction",
        ):
            manager.validate_killms_solutions(invalid)


if __name__ == "__main__":
    unittest.main()
