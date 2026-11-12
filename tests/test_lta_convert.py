import importlib.util
import json
import os
import stat
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]


def _load(name, module_name):
    spec = importlib.util.spec_from_file_location(module_name, SOURCE_DIR / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


lta_convert = _load("lta_convert.py", "test_radio_continuum_lta_convert")
fits_ingest = _load("fits_ingest.py", "test_radio_continuum_lta_fits_ingest")


def _card(keyword, value=None):
    if value is None:
        text = keyword
    elif isinstance(value, bool):
        text = f"{keyword:<8}= {'T' if value else 'F':>20}"
    elif isinstance(value, str):
        text = f"{keyword:<8}= '{value}'"
    else:
        text = f"{keyword:<8}= {value:>20}"
    return text[:80].ljust(80).encode("ascii")


def uvfits_bytes(telescope="GMRT"):
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
        ("TELESCOP", telescope),
    ]
    header = b"".join(_card(key, value) for key, value in cards) + _card("END")
    header += b" " * ((-len(header)) % 2880)
    return header + b"\0" * 2880


def image_fits_bytes():
    cards = [
        ("SIMPLE", True),
        ("BITPIX", -32),
        ("NAXIS", 2),
        ("NAXIS1", 2),
        ("NAXIS2", 2),
    ]
    header = b"".join(_card(key, value) for key, value in cards) + _card("END")
    header += b" " * ((-len(header)) % 2880)
    return header + b"\0" * 2880


class FakeTools:
    """Stand-ins for NCRA's listscan and gvfits, written as Python scripts.

    listscan writes ``<name>.log`` naming TEST.FITS into its working folder;
    gvfits copies a prepared FITS file to the name in the log. ``mode``
    switches in the failures under test.
    """

    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.template = self.root / "template.fits"
        self.template.write_bytes(uvfits_bytes())
        self.calls = self.root / "calls.txt"
        self.listscan = self.write("listscan", self.listscan_source())
        self.gvfits = self.write("gvfits", self.gvfits_source())

    def write(self, name, body):
        path = self.root / name
        path.write_text(f"#!{sys.executable}\n" + textwrap.dedent(body))
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        return path

    def listscan_source(self, mode="normal"):
        return f"""
        import pathlib, sys
        with open({str(self.calls)!r}, "a") as calls:
            calls.write("listscan " + " ".join(sys.argv[1:]) + "\\n")
        lta = pathlib.Path(sys.argv[1])
        if not lta.is_file():
            sys.exit("cannot open " + str(lta))
        mode = {mode!r}
        if mode == "fail":
            print("listscan: not an LTA file")
            sys.exit(3)
        if mode == "nolib":
            print(sys.argv[0] + ": error while loading shared libraries: "
                  "libgfortran.so.3: cannot open shared object file: "
                  "No such file or directory")
            sys.exit(127)
        import os
        print("LD_LIBRARY_PATH=" + os.environ.get("LD_LIBRARY_PATH", ""))
        stem = lta.name.split(".")[0]
        if mode == "beside":
            target = lta.resolve().with_name(stem + ".log")
        else:
            target = pathlib.Path.cwd() / (stem + ".log")
        name = "/tmp/elsewhere/OUT.FITS" if mode == "absolute" else "TEST.FITS"
        if mode != "nolog":
            target.write_text(
                "LTAFILE  " + lta.name + "\\n"
                "FITS     " + name + "\\n"
                "SCAN 1   3C48\\n"
            )
        print("listed 1 scan")
        """

    def gvfits_source(self, mode="normal"):
        return f"""
        import os, pathlib, re, shutil, sys
        with open({str(self.calls)!r}, "a") as calls:
            calls.write("gvfits " + " ".join(sys.argv[1:]) + "\\n")
        log = pathlib.Path(sys.argv[1]).read_text()
        lta = re.search(r"LTAFILE\\s+(\\S+)", log).group(1)
        if not pathlib.Path(lta).is_file():
            sys.exit("cannot open " + lta)
        name = re.search(r"FITS\\s+(\\S+)", log).group(1)
        mode = {mode!r}
        if mode == "fail":
            print("gvfits: bad record")
            sys.exit(4)
        if mode == "nothing":
            sys.exit(0)
        if mode == "touch":
            target = pathlib.Path(lta).resolve()
            os.utime(target, ns=(0, target.stat().st_mtime_ns + 10**9))
        shutil.copyfile({str(self.template)!r}, name)
        pathlib.Path("gvfits.log").write_text("converted " + lta + "\\n")
        if mode == "two":
            shutil.copyfile({str(self.template)!r}, "OTHER.FITS")
        print("wrote " + name)
        """

    def set_listscan(self, mode):
        self.write("listscan", self.listscan_source(mode))

    def set_gvfits(self, mode):
        self.write("gvfits", self.gvfits_source(mode))

    def call_count(self, name):
        if not self.calls.exists():
            return 0
        return sum(
            1 for line in self.calls.read_text().splitlines()
            if line.startswith(name + " ")
        )


class LTANamesTests(unittest.TestCase):
    def test_lta_names(self):
        for name in ["obs.lta", "obs.LTA", "obs.ltb", "obs.lta_1", "28_066_gwb.lta_12"]:
            with self.subTest(name=name):
                self.assertTrue(lta_convert.is_lta_path(name))
                self.assertTrue(fits_ingest.is_lta_path(name))
        for name in ["obs.fits", "obs.ms", "obs.lta.fits", "obs.ltc", "lta"]:
            with self.subTest(name=name):
                self.assertFalse(lta_convert.is_lta_path(name))
                self.assertFalse(fits_ingest.is_lta_path(name))

    def test_stems_keep_parts_apart(self):
        self.assertEqual(lta_convert.lta_stem("/a/28_066_22NOV2020_gwb.lta"), "28_066_22NOV2020_gwb")
        self.assertEqual(lta_convert.lta_stem("obs.ltb"), "obs_ltb")
        self.assertEqual(lta_convert.lta_stem("obs.LTA_2"), "obs_lta_2")
        self.assertEqual(lta_convert.lta_stem("my obs.v2.lta"), "my_obs_v2")

    def test_fits_names_in_log(self):
        text = "FITS  TEST.FITS\nOTHER 'b.fits'\nFITS  TEST.FITS\n"
        self.assertEqual(lta_convert.fits_names_in_log(text), ["TEST.FITS", "b.fits"])


class LTAConversionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.archive = self.root / "archive"
        self.archive.mkdir()
        self.lta = self.archive / "28_066_22NOV2020_gwb.lta"
        self.lta.write_bytes(b"GMRT LTA stand-in " * 100)
        self.tools = FakeTools(self.root / "tools")
        self.output = self.root / "Products" / "ingest" / "lta" / "28_066_22NOV2020_gwb.FITS"
        self.messages = []

    def tearDown(self):
        self.temporary.cleanup()

    def convert(self, **values):
        settings = {
            "listscan": str(self.tools.listscan),
            "gvfits": str(self.tools.gvfits),
            "echo": self.messages.append,
        }
        settings.update(values)
        return lta_convert.convert(self.lta, self.output, **settings)

    def leftovers(self):
        return sorted(path.name for path in self.output.parent.glob(".*partial*"))

    def test_converts_and_records_without_writing_beside_the_lta(self):
        before = sorted(os.listdir(self.archive))
        record = self.convert()

        self.assertFalse(record["reused"])
        self.assertTrue(self.output.is_file())
        self.assertEqual(self.output.read_bytes(), self.tools.template.read_bytes())
        self.assertEqual(sorted(os.listdir(self.archive)), before)
        self.assertEqual(self.leftovers(), [])
        self.assertEqual(record["lta"]["path"], str(self.lta.resolve()))
        self.assertEqual(len(record["lta"]["sha256"]), 64)
        self.assertEqual(record["fits"]["path"], str(self.output))
        self.assertEqual(record["fits_name_in_log"], "TEST.FITS")
        self.assertEqual(record["listscan"]["command"], [str(self.tools.listscan.resolve())])
        self.assertEqual(len(record["gvfits"]["sha256"]), 64)
        log = Path(record["kept_files"]["listscan.log"])
        self.assertIn("TEST.FITS", log.read_text())
        self.assertIn("wrote TEST.FITS", Path(record["kept_files"]["gvfits.out"]).read_text())
        self.assertIn("converted", Path(record["kept_files"]["gvfits.log"]).read_text())
        self.assertTrue(self.output.with_name(self.output.name + lta_convert.RECORD_SUFFIX).is_file())
        # The output is GMRT random-groups UVFITS, as importgmrt expects.
        summary = fits_ingest.inspect_fits_structure(self.output)
        self.assertEqual(fits_ingest.resolve_visibility_format(summary), "uvfits")
        self.assertEqual(summary["detected_telescopes"], ["ugmrt"])

    def test_second_conversion_of_the_same_file_is_reused(self):
        self.convert()
        self.tools.set_listscan("fail")
        record = self.convert()
        self.assertTrue(record["reused"])
        self.assertEqual(self.tools.call_count("listscan"), 1)

    def test_changed_lta_is_converted_again(self):
        self.convert()
        stat_now = self.lta.stat()
        os.utime(self.lta, ns=(stat_now.st_atime_ns, stat_now.st_mtime_ns + 10**9))
        record = self.convert()
        self.assertFalse(record["reused"])
        self.assertEqual(self.tools.call_count("listscan"), 2)

    def test_listscan_failure_is_reported_and_leaves_nothing(self):
        self.tools.set_listscan("fail")
        with self.assertRaisesRegex(lta_convert.LTAConversionError, r"status 3(.|\n)*not an LTA file"):
            self.convert()
        self.assertFalse(self.output.exists())
        self.assertEqual(self.leftovers(), [])
        self.assertEqual(self.tools.call_count("gvfits"), 0)

    def test_missing_log_is_an_error(self):
        self.tools.set_listscan("nolog")
        with self.assertRaisesRegex(lta_convert.LTAConversionError, "no .log file"):
            self.convert()
        self.assertEqual(self.tools.call_count("gvfits"), 0)

    def test_log_written_beside_the_lta_is_named_in_the_error(self):
        self.tools.set_listscan("beside")
        with self.assertRaisesRegex(lta_convert.LTAConversionError, "beside the LTA file"):
            self.convert()

    def test_output_outside_the_working_folder_is_refused_before_gvfits(self):
        self.tools.set_listscan("absolute")
        with self.assertRaisesRegex(lta_convert.LTAConversionError, "outside the working"):
            self.convert()
        self.assertEqual(self.tools.call_count("gvfits"), 0)

    def test_gvfits_failure_and_missing_output(self):
        for mode, pattern in [("fail", r"status 4(.|\n)*bad record"), ("nothing", "wrote no FITS")]:
            with self.subTest(mode=mode):
                self.tools.set_gvfits(mode)
                with self.assertRaisesRegex(lta_convert.LTAConversionError, pattern):
                    self.convert()
                self.assertFalse(self.output.exists())
                self.assertEqual(self.leftovers(), [])

    def test_the_file_named_in_the_log_is_taken(self):
        self.tools.set_gvfits("two")
        record = self.convert()
        self.assertEqual(record["fits_name_in_log"], "TEST.FITS")

    def test_lta_changed_during_conversion_is_refused(self):
        self.tools.set_gvfits("touch")
        with self.assertRaisesRegex(lta_convert.LTAConversionError, "changed while"):
            self.convert()
        self.assertFalse(self.output.exists())

    def test_missing_and_unexecutable_tools(self):
        with self.assertRaisesRegex(FileNotFoundError, "CAPTURE-CASA6"):
            self.convert(listscan=str(self.root / "nowhere" / "listscan"))
        plain = self.root / "tools" / "plain"
        plain.write_text("not a program")
        plain.chmod(0o644)
        with self.assertRaisesRegex(FileNotFoundError, "chmod \\+x"):
            self.convert(gvfits=str(plain))

    def test_relative_tool_and_wrapper_command(self):
        # The programs run inside the staging folder, so a relative path
        # must be resolved first, and a command with arguments must work.
        previous = os.getcwd()
        os.chdir(self.tools.root)
        try:
            record = self.convert(
                listscan="./listscan",
                gvfits=f"{sys.executable} {self.tools.gvfits}",
            )
        finally:
            os.chdir(previous)
        self.assertEqual(record["listscan"]["command"], [str(self.tools.listscan.resolve())])
        self.assertEqual(record["gvfits"]["command"][1], str(self.tools.gvfits))
        self.assertTrue(self.output.is_file())

    def test_libraries_beside_the_programs_are_found_first(self):
        (self.tools.root / "lib").mkdir()
        with mock.patch.dict(os.environ, {"LD_LIBRARY_PATH": "/opt/other"}):
            record = self.convert()
        said = Path(record["kept_files"]["listscan.out"]).read_text()
        self.assertIn(
            f"LD_LIBRARY_PATH={self.tools.root}:{self.tools.root / 'lib'}:/opt/other",
            said,
        )

    def test_a_missing_library_is_explained(self):
        self.tools.set_listscan("nolib")
        with self.assertRaises(lta_convert.LTAConversionError) as raised:
            self.convert()
        message = str(raised.exception)
        self.assertIn("status 127", message)
        self.assertIn("listscan needs libgfortran.so.3", message)
        self.assertIn(f"{self.tools.root / 'lib'}/", message)
        self.assertIn("libgfortran-3.0.0", message)

    def test_too_little_space_stops_before_conversion(self):
        usage = mock.Mock(free=10)
        with mock.patch.object(lta_convert.shutil, "disk_usage", return_value=usage):
            with self.assertRaisesRegex(lta_convert.LTAConversionError, "will not fit"):
                self.convert()
        self.assertEqual(self.tools.call_count("listscan"), 0)

    def test_command_line(self):
        code = lta_convert.main([
            str(self.lta),
            str(self.output),
            "--listscan", str(self.tools.listscan),
            "--gvfits", str(self.tools.gvfits),
        ])
        self.assertEqual(code, 0)
        self.assertTrue(self.output.is_file())



class LTAIngestTests(unittest.TestCase):
    """fits_ingest.py converts an LTA file, then imports it with importgmrt."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.lta = self.root / "archive" / "obs_gwb.lta"
        self.lta.parent.mkdir()
        self.lta.write_bytes(b"GMRT LTA stand-in " * 100)
        self.tools = FakeTools(self.root / "tools")
        self.outdir = self.root / "Products"
        self.output = self.outdir / "ingest" / "obs_gwb.ms"
        self.converted = self.outdir / "ingest" / "lta" / "obs_gwb.FITS"
        self.calls = []
        self.fail_import = False

    def tearDown(self):
        self.temporary.cleanup()

    @staticmethod
    def validator(path):
        path = Path(path)
        if not (path / "table.dat").is_file():
            raise fits_ingest.MeasurementSetValidationError("missing fake MS")
        return {
            "main_rows": 2,
            "main_columns": sorted(fits_ingest.REQUIRED_MAIN_COLUMNS | {"DATA"}),
            "data_columns": ["DATA"],
            "subtables": {},
        }

    def loader(self, name, required=True):
        def task(**kwargs):
            self.calls.append((name, dict(kwargs)))
            # importgmrt must be given the converted FITS, which exists.
            if not Path(kwargs["fitsfile"]).is_file():
                raise RuntimeError("the FITS to import is missing")
            if self.fail_import:
                raise RuntimeError("synthetic importgmrt failure")
            destination = Path(kwargs["vis"])
            destination.mkdir()
            (destination / "table.dat").touch()

        return task

    def ingestor(self, **values):
        settings = {
            "input_fits": self.lta,
            "output_ms": self.output,
            "outdir": self.outdir,
            "task_loader": self.loader,
            "ms_validator": self.validator,
            "listscan": str(self.tools.listscan),
            "gvfits": str(self.tools.gvfits),
        }
        settings.update(values)
        return fits_ingest.FITSIngestor(**settings)

    def report(self):
        return json.loads((self.outdir / "ingest" / "ingest_report.json").read_text())

    def test_lta_is_converted_then_imported_with_importgmrt(self):
        result = self.ingestor().run()

        self.assertEqual(result, self.output)
        self.assertEqual([name for name, _ in self.calls], ["importgmrt"])
        self.assertEqual(self.calls[0][1]["fitsfile"], str(self.converted))
        report = self.report()
        self.assertEqual(report["source_format"], "lta")
        self.assertEqual(report["source_fits"]["path"], str(self.lta))
        self.assertEqual(len(report["source_fits"]["sha256"]), 64)
        self.assertEqual(report["fits_format"], "uvfits")
        self.assertEqual(report["casa_task"], "importgmrt")
        self.assertEqual(report["telescope"], "ugmrt")
        self.assertEqual(report["lta_conversion"]["lta"], report["source_fits"])
        self.assertFalse(report["converted_fits_kept"])
        self.assertEqual(report["fits_summary"]["detected_telescopes"], ["ugmrt"])
        self.assertTrue((self.outdir / "ingest" / "ingest.done").is_file())
        # The large FITS is gone once the MS is validated; the logs stay.
        self.assertFalse(self.converted.exists())
        self.assertFalse(
            self.converted.with_name(self.converted.name + lta_convert.RECORD_SUFFIX).exists()
        )
        self.assertTrue(self.converted.with_name(self.converted.name + ".listscan.log").is_file())

    def test_finished_import_is_reused_without_converting(self):
        self.ingestor().run()
        self.tools.set_listscan("fail")
        self.calls.clear()
        self.assertEqual(self.ingestor().run(), self.output)
        self.assertEqual(self.calls, [])
        self.assertEqual(self.tools.call_count("listscan"), 1)

    def test_converted_fits_can_be_kept(self):
        self.ingestor(keep_converted_fits=True).run()
        self.assertTrue(self.converted.is_file())
        self.assertTrue(self.report()["converted_fits_kept"])

    def test_failed_import_keeps_the_fits_for_the_retry(self):
        self.fail_import = True
        with self.assertRaisesRegex(RuntimeError, "synthetic importgmrt failure"):
            self.ingestor().run()
        self.assertTrue(self.converted.is_file())
        self.assertFalse(self.output.exists())
        self.assertFalse((self.outdir / "ingest" / "ingest.done").exists())

        self.fail_import = False
        self.ingestor().run()
        self.assertEqual(self.tools.call_count("listscan"), 1)
        self.assertTrue(self.output.is_dir())
        self.assertTrue(self.report()["lta_conversion"]["reused"])

    def test_other_telescope_is_refused_before_converting(self):
        with self.assertRaisesRegex(fits_ingest.FITSDetectionError, "GMRT's own format"):
            self.ingestor(telescope="meerkat").run()
        self.assertEqual(self.tools.call_count("listscan"), 0)
        self.assertEqual(self.calls, [])

    def test_explicit_lta_format_for_any_name(self):
        renamed = self.lta.with_name("raw.data")
        self.lta.rename(renamed)
        self.lta = renamed
        self.converted = self.outdir / "ingest" / "lta" / "raw_data.FITS"
        self.ingestor(fits_format="lta").run()
        self.assertEqual(self.calls[0][1]["fitsfile"], str(self.converted))
        self.assertEqual(self.report()["requested_fits_format"], "lta")

    def test_lta_name_with_uvfits_format_is_read_as_fits(self):
        with self.assertRaises(fits_ingest.FITSDetectionError):
            self.ingestor(fits_format="uvfits").run()
        self.assertEqual(self.tools.call_count("listscan"), 0)

    def test_gvfits_output_that_is_not_uvfits_is_refused(self):
        self.tools.template.write_bytes(image_fits_bytes())
        with self.assertRaisesRegex(fits_ingest.FITSDetectionError, "random-groups UVFITS"):
            self.ingestor().run()
        self.assertEqual(self.calls, [])
        self.assertFalse(self.output.exists())

    def test_command_line_options(self):
        args = fits_ingest.build_parser().parse_args([
            str(self.lta), str(self.output),
            "--fits-format", "lta",
            "--listscan", "/opt/gmrt/listscan",
            "--gvfits", "singularity exec /x.sif /opt/gmrt/gvfits",
            "--keep-lta-fits",
        ])
        self.assertEqual(args.fits_format, "lta")
        self.assertEqual(args.listscan, "/opt/gmrt/listscan")
        self.assertEqual(args.gvfits, "singularity exec /x.sif /opt/gmrt/gvfits")
        self.assertTrue(args.keep_lta_fits)


if __name__ == "__main__":
    unittest.main()
