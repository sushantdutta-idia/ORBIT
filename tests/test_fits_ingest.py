import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def _load_fits_ingest():
    path = SOURCE_DIR / "fits_ingest.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_fits_ingest",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fits_ingest = _load_fits_ingest()


def _card(keyword, value=None):
    if value is None:
        text = keyword
    elif isinstance(value, bool):
        text = f"{keyword:<8}= {'T' if value else 'F':>20}"
    elif isinstance(value, str):
        escaped = value.replace("'", "''")
        text = f"{keyword:<8}= '{escaped}'"
    else:
        text = f"{keyword:<8}= {value:>20}"
    return text[:80].ljust(80).encode("ascii")


def _hdu(cards, data=b""):
    header = b"".join(_card(key, value) for key, value in cards)
    header += _card("END")
    header += b" " * ((-len(header)) % fits_ingest.FITS_BLOCK_SIZE)
    data += b"\0" * ((-len(data)) % fits_ingest.FITS_BLOCK_SIZE)
    return header + data


def _uvfits_bytes(telescope=None):
    cards = [
        ("SIMPLE", True),
        ("BITPIX", -32),
        ("NAXIS", 3),
        ("NAXIS1", 0),
        ("NAXIS2", 3),
        ("NAXIS3", 1),
        ("GROUPS", True),
        ("PCOUNT", 5),
        ("GCOUNT", 2),
    ]
    if telescope:
        cards.append(("TELESCOP", telescope))
    # 4 bytes * 2 groups * (5 parameters + 3 data values).
    return _hdu(cards, b"\0" * 64)


def _fitsidi_extension():
    return _hdu(
        [
            ("XTENSION", "BINTABLE"),
            ("BITPIX", 8),
            ("NAXIS", 2),
            ("NAXIS1", 16),
            ("NAXIS2", 1),
            ("PCOUNT", 0),
            ("GCOUNT", 1),
            ("TFIELDS", 1),
            ("EXTNAME", "UV_DATA"),
        ],
        b"\0" * 16,
    )


def _fitsidi_bytes():
    return _hdu(
        [
            ("SIMPLE", True),
            ("BITPIX", 8),
            ("NAXIS", 0),
            ("EXTEND", True),
        ]
    ) + _fitsidi_extension()


def _image_fits_bytes():
    return _hdu(
        [
            ("SIMPLE", True),
            ("BITPIX", -32),
            ("NAXIS", 2),
            ("NAXIS1", 2),
            ("NAXIS2", 2),
            ("PCOUNT", 0),
            ("GCOUNT", 1),
        ],
        b"\0" * 16,
    )


class FITSDetectionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, name, value):
        path = self.root / name
        path.write_bytes(value)
        return path

    def test_detects_random_groups_by_content_not_extension(self):
        path = self.write("visibilities.notfits", _uvfits_bytes())
        self.assertEqual(fits_ingest.detect_visibility_fits(path), "uvfits")
        summary = fits_ingest.inspect_fits_structure(path)
        self.assertTrue(summary["random_groups_uvfits"])
        self.assertFalse(summary["fitsidi_uv_data"])

    def test_detects_fitsidi_uv_data_table(self):
        path = self.write("idi.fits", _fitsidi_bytes())
        self.assertEqual(fits_ingest.detect_visibility_fits(path), "fitsidi")
        self.assertEqual(
            fits_ingest.detect_visibility_fits(path, "fits-idi"),
            "fitsidi",
        )

    def test_rejects_ordinary_image_even_with_uvfits_override(self):
        path = self.write("continuum.fits", _image_fits_bytes())
        with self.assertRaisesRegex(
            fits_ingest.FITSDetectionError,
            "ordinary image FITS",
        ):
            fits_ingest.detect_visibility_fits(path)
        with self.assertRaisesRegex(
            fits_ingest.FITSDetectionError,
            "does not match",
        ):
            fits_ingest.detect_visibility_fits(path, "uvfits")

    def test_auto_rejects_ambiguous_but_override_selects_signature(self):
        path = self.write("ambiguous.fits", _uvfits_bytes() + _fitsidi_extension())
        with self.assertRaisesRegex(fits_ingest.FITSDetectionError, "both"):
            fits_ingest.detect_visibility_fits(path)
        self.assertEqual(
            fits_ingest.detect_visibility_fits(path, "uvfits"),
            "uvfits",
        )
        self.assertEqual(
            fits_ingest.detect_visibility_fits(path, "fitsidi"),
            "fitsidi",
        )

    def test_truncated_payload_fails_before_casa(self):
        path = self.write("truncated.fits", _uvfits_bytes()[:-fits_ingest.FITS_BLOCK_SIZE])
        with self.assertRaises(fits_ingest.FITSDetectionError):
            fits_ingest.detect_visibility_fits(path)


class FITSProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.source = self.root / "raw.fits"
        self.source.write_bytes(_uvfits_bytes())
        self.output = self.root / "Products" / "ingest" / "raw.ms"
        self.outdir = self.root / "Products"
        self.calls = []

    def tearDown(self):
        self.temporary.cleanup()

    @staticmethod
    def validator(path):
        path = Path(path)
        if not (path / "table.dat").is_file():
            raise fits_ingest.MeasurementSetValidationError("missing fake MS")
        return {
            "main_rows": 2,
            "main_columns": sorted(
                fits_ingest.REQUIRED_MAIN_COLUMNS | {"DATA"}
            ),
            "data_columns": ["DATA"],
            "subtables": {},
        }

    def loader(self, name, required=True):
        self.calls.append((name, required))

        def task(**kwargs):
            destination = Path(kwargs["vis"])
            destination.mkdir()
            (destination / "table.dat").touch()

        return task

    def manager(self, **values):
        settings = {
            "input_fits": self.source,
            "output_ms": self.output,
            "outdir": self.outdir,
            "task_loader": self.loader,
            "ms_validator": self.validator,
        }
        settings.update(values)
        return fits_ingest.FITSIngestor(**settings)

    def test_import_report_checkpoint_and_matching_reuse(self):
        result = self.manager().run()
        self.assertEqual(result, self.output.resolve())
        self.assertEqual(self.calls, [("importuvfits", True)])

        report_path = self.outdir / "ingest" / "ingest_report.json"
        marker = self.outdir / "ingest" / "ingest.done"
        report = json.loads(report_path.read_text())
        self.assertTrue(marker.is_file())
        self.assertEqual(report["status"], "completed")
        self.assertEqual(report["fits_format"], "uvfits")
        self.assertEqual(report["casa_task"], "importuvfits")
        self.assertEqual(report["measurement_set"], str(self.output.resolve()))
        self.assertEqual(report["source_fits"]["path"], str(self.source.resolve()))
        self.assertEqual(len(report["source_fits"]["sha256"]), 64)

        def forbidden_loader(*args, **kwargs):
            raise AssertionError("matching output should have been reused")

        reused = fits_ingest.FITSIngestor(
            self.source,
            self.output,
            outdir=self.outdir,
            task_loader=forbidden_loader,
            ms_validator=self.validator,
        ).run()
        self.assertEqual(reused, self.output.resolve())

    def test_existing_output_is_not_reused_after_source_identity_changes(self):
        self.manager().run()
        stat = self.source.stat()
        os.utime(
            self.source,
            ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000),
        )
        with self.assertRaisesRegex(RuntimeError, "without a matching"):
            self.manager().run()
        self.assertEqual(len(self.calls), 1)

    def test_ugmrt_uvfits_uses_specific_importer(self):
        self.manager(telescope="ugmrt").run()
        self.assertEqual(self.calls, [("importgmrt", True)])
        report = json.loads(
            (self.outdir / "ingest" / "ingest_report.json").read_text()
        )
        self.assertEqual(report["telescope"], "ugmrt")
        self.assertEqual(report["casa_task"], "importgmrt")

    def test_fits_telescope_conflict_fails_before_task_loading(self):
        self.source.write_bytes(_uvfits_bytes("MeerKAT"))
        with self.assertRaisesRegex(fits_ingest.FITSDetectionError, "conflicts"):
            self.manager(telescope="ugmrt").run()
        self.assertEqual(self.calls, [])

    def test_task_failure_removes_partial_and_writes_no_checkpoint(self):
        def failing_loader(name, required=True):
            def task(**kwargs):
                destination = Path(kwargs["vis"])
                destination.mkdir()
                (destination / "partial").touch()
                raise RuntimeError("synthetic CASA failure")

            return task

        manager = fits_ingest.FITSIngestor(
            self.source,
            self.output,
            outdir=self.outdir,
            task_loader=failing_loader,
            ms_validator=self.validator,
        )
        with self.assertRaisesRegex(RuntimeError, "synthetic CASA failure"):
            manager.run()
        self.assertFalse(self.output.exists())
        self.assertFalse((self.outdir / "ingest" / "ingest.done").exists())
        self.assertEqual(list(self.output.parent.glob(".*.partial.*")), [])

    def test_output_cannot_contain_ingest_provenance_directory(self):
        with self.assertRaisesRegex(ValueError, "cannot contain"):
            fits_ingest.FITSIngestor(
                self.source,
                self.outdir,
                outdir=self.outdir,
                task_loader=self.loader,
                ms_validator=self.validator,
            )


if __name__ == "__main__":
    unittest.main()
