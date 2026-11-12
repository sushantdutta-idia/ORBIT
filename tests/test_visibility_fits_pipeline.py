import contextlib
import hashlib
import importlib.util
import io
import json
import os
import tempfile
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_pipeline_module():
    path = SOURCE_DIR / "pipeline.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_pipeline_visibility_fits",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PIPELINE = load_pipeline_module()


def fits_card(keyword, value=None):
    """Return one simple, standards-sized FITS header card."""

    if keyword == "END":
        text = "END"
    elif isinstance(value, bool):
        text = f"{keyword:<8}= {'T' if value else 'F':>20}"
    elif isinstance(value, str):
        text = f"{keyword:<8}= '{value}'"
    else:
        text = f"{keyword:<8}= {value:>20}"
    return text.ljust(80).encode("ascii")


def write_random_groups_uvfits(path):
    """Write enough of a random-groups header for safe type detection."""

    cards = [
        fits_card("SIMPLE", True),
        fits_card("BITPIX", -32),
        fits_card("NAXIS", 3),
        fits_card("NAXIS1", 0),
        fits_card("NAXIS2", 1),
        fits_card("NAXIS3", 1),
        fits_card("GROUPS", True),
        fits_card("PCOUNT", 5),
        fits_card("GCOUNT", 1),
        fits_card("END"),
    ]
    header = b"".join(cards)
    # BITPIX=-32, GCOUNT=1, PCOUNT=5 and one array element declare 24
    # payload bytes; every FITS HDU payload is padded to 2880 bytes.
    payload = bytes(24).ljust(2880, b"\0")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header.ljust(2880, b" ") + payload)
    return path


def write_image_fits(path):
    """Write a valid-looking image FITS header, never visibility data."""

    cards = [
        fits_card("SIMPLE", True),
        fits_card("BITPIX", -32),
        fits_card("NAXIS", 2),
        fits_card("NAXIS1", 2),
        fits_card("NAXIS2", 2),
        fits_card("EXTEND", True),
        fits_card("END"),
    ]
    header = b"".join(cards).ljust(2880, b" ")
    data = bytes(16).ljust(2880, b"\0")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + data)
    return path


def make_measurement_set(path):
    path.mkdir(parents=True, exist_ok=True)
    (path / "table.dat").touch()
    return path


class VisibilityFITSPipelineTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.outdir = self.root / "Products"

    def tearDown(self):
        self.temporary.cleanup()

    def parse(self, input_path, *extra):
        argv = [
            str(input_path),
            "--outdir",
            str(self.outdir),
            "--dry-run",
            "--mpi",
            "off",
            "--no-copy-input-ms",
            "--telescope",
            "ugmrt",
            "--stages",
            "inspection",
            *extra,
        ]
        return PIPELINE.build_parser().parse_args(argv)

    def runner(self, input_path, *extra):
        return PIPELINE.PipelineRunner(self.parse(input_path, *extra))

    def test_measurement_set_keeps_the_existing_stage_path(self):
        measurement_set = make_measurement_set(self.root / "raw.ms")
        runner = self.runner(measurement_set)

        self.assertEqual(runner.selected_stages(), ["inspection"])
        command = runner.command_for_stage("inspection")
        self.assertIn(str(measurement_set.resolve()), command)
        self.assertNotIn("fits_ingest.py", command)

    def test_uvfits_and_generic_fits_plan_ingest_before_inspection(self):
        for filename in ["raw.uvfits", "raw.fits"]:
            with self.subTest(filename=filename):
                source = write_random_groups_uvfits(self.root / filename)
                runner = self.runner(source)

                self.assertEqual(
                    runner.selected_stages(),
                    ["ingest", "inspection"],
                )
                expected_ms = (
                    self.outdir / "ingest" / f"{source.stem}.ms"
                ).resolve()
                ingest = runner.command_for_stage("ingest")
                inspection = runner.command_for_stage("inspection")

                self.assertIn(str(SOURCE_DIR / "fits_ingest.py"), ingest)
                self.assertIn(str(source.resolve()), ingest)
                self.assertIn(str(expected_ms), ingest)
                self.assertEqual(
                    ingest[ingest.index("--telescope") + 1],
                    "ugmrt",
                )
                self.assertIn(str(expected_ms), inspection)
                self.assertNotIn(str(source.resolve()), inspection)

    def test_explicit_fits_format_is_accepted_and_forwarded(self):
        source = write_random_groups_uvfits(self.root / "raw.data")
        runner = self.runner(source, "--fits-format", "uvfits")

        self.assertEqual(runner.args.fits_format, "uvfits")
        command = runner.command_for_stage("ingest")
        self.assertEqual(
            command[command.index("--fits-format") + 1],
            "uvfits",
        )

        fitsidi_args = self.parse(source, "--fits-format", "fitsidi")
        self.assertEqual(fitsidi_args.fits_format, "fitsidi")

    def test_imported_ms_override_becomes_the_downstream_input(self):
        source = write_random_groups_uvfits(self.root / "raw.uvfits")
        imported_ms = self.root / "scratch" / "converted.ms"
        runner = self.runner(
            source,
            "--imported-ms",
            str(imported_ms),
        )

        ingest = runner.command_for_stage("ingest")
        inspection = runner.command_for_stage("inspection")
        expected = str(imported_ms.resolve())
        self.assertIn(expected, ingest)
        self.assertIn(expected, inspection)

    def test_explicit_downstream_ms_does_not_trigger_unrelated_fits_ingest(self):
        source = write_random_groups_uvfits(self.root / "archival.uvfits")
        imaging_ms = self.root / "already_calibrated.ms"
        runner = self.runner(
            source,
            "--data-state",
            "calibrated",
            "--stages",
            "imaging",
            "--imaging-ms",
            str(imaging_ms),
        )

        self.assertEqual(runner.selected_stages(), ["imaging"])
        command = runner.command_for_stage("imaging")
        self.assertIn(str(imaging_ms), command)
        self.assertNotIn("fits_ingest.py", command)

    def test_imported_ms_cannot_contain_ingest_provenance(self):
        source = write_random_groups_uvfits(self.root / "raw.uvfits")
        runner = self.runner(
            source,
            "--imported-ms",
            str(self.outdir),
        )

        with self.assertRaisesRegex(RuntimeError, "cannot contain"):
            runner.preflight(runner.selected_stages())

    def test_generated_sbatch_keeps_fits_source_and_import_options(self):
        source = write_random_groups_uvfits(self.root / "raw.uvfits")
        imported_ms = self.root / "scratch" / "raw.ms"
        sbatch = self.root / "fits_pipeline.sbatch"
        args = self.parse(
            source,
            "--fits-format",
            "uvfits",
            "--imported-ms",
            str(imported_ms),
            "--write-sbatch",
            str(sbatch),
        )

        with contextlib.redirect_stdout(io.StringIO()):
            PIPELINE.write_sbatch_from_config(args)

        content = sbatch.read_text()
        self.assertIn(str(source), content)
        self.assertIn("--fits-format uvfits", content)
        self.assertIn(f"--imported-ms {imported_ms}", content)

    def test_image_fits_is_rejected_before_any_casa_stage_is_planned(self):
        image = write_image_fits(self.root / "restored_image.fits")

        with self.assertRaisesRegex(
            (RuntimeError, ValueError),
            r"(?i)(image|visibility).*(FITS|visibility)|FITS.*(image|visibility)",
        ):
            runner = self.runner(image)
            stages = runner.selected_stages()
            runner.preflight(stages)

    def test_lofar_visibility_fits_is_not_treated_as_a_linc_product(self):
        source = write_random_groups_uvfits(self.root / "lofar.uvfits")
        args = self.parse(
            source,
            "--telescope", "lofar",
            "--receiver-band", "HBA",
            "--data-state", "calibrated",
            "--imager", "wsclean",
            "--lofar-beam-mode", "full",
            "--stages", "inspection,imaging",
        )
        runner = PIPELINE.PipelineRunner(args)

        with self.assertRaisesRegex(RuntimeError, "finalized.*Measurement Set"):
            runner.preflight(runner.selected_stages())

    def test_ingest_fingerprint_binds_resume_to_the_source_input(self):
        first = write_random_groups_uvfits(
            self.root / "first" / "same_name.uvfits"
        )
        second = write_random_groups_uvfits(
            self.root / "second" / "same_name.uvfits"
        )

        first_runner = self.runner(first)
        second_runner = self.runner(second)

        # Both inputs intentionally resolve to the same default output MS.
        self.assertEqual(
            first_runner.imported_ms,
            second_runner.imported_ms,
        )
        self.assertNotEqual(
            first_runner.stage_fingerprint("ingest"),
            second_runner.stage_fingerprint("ingest"),
        )

    def test_ingest_fingerprint_does_not_change_after_auto_inspection(self):
        source = write_random_groups_uvfits(self.root / "raw.uvfits")
        runner = self.runner(source, "--telescope", "auto")
        before = runner.stage_fingerprint("ingest")
        inspection = runner.report_path("inspection")
        inspection.parent.mkdir(parents=True, exist_ok=True)
        inspection.write_text(
            json.dumps(
                {
                    "measurement_set": str(runner.imported_ms),
                    "telescope": {"name": "ugmrt"},
                    "band": {"name": "UGMRT_B4"},
                }
            )
        )
        runner._stage_fingerprints.clear()

        after = runner.stage_fingerprint("ingest")

        self.assertEqual(before, after)

    def test_custom_stage_order_cannot_put_inspection_before_ingest(self):
        source = write_random_groups_uvfits(self.root / "raw.uvfits")
        runner = self.runner(source, "--stages", "inspection,ingest")

        self.assertEqual(runner.selected_stages(), ["ingest", "inspection"])

    def test_ingest_resume_rejects_a_report_for_another_source(self):
        first = write_random_groups_uvfits(
            self.root / "first" / "same_name.uvfits"
        )
        second = write_random_groups_uvfits(
            self.root / "second" / "same_name.uvfits"
        )
        runner = self.runner(second)
        imported_ms = make_measurement_set(runner.imported_ms)
        ingest_dir = self.outdir / "ingest"
        (ingest_dir / "ingest.done").touch()
        first_stat = first.stat()
        (ingest_dir / "ingest_report.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "stage": "ingest",
                    "status": "completed",
                    "source_fits": {
                        "path": str(first.resolve()),
                        "size_bytes": first_stat.st_size,
                        "mtime_ns": first_stat.st_mtime_ns,
                        "sha256": "not-used-for-the-path-mismatch",
                    },
                    "requested_fits_format": "auto",
                    "fits_format": "uvfits",
                    "casa_task": "importgmrt",
                    "requested_telescope": "ugmrt",
                    "telescope": "ugmrt",
                    "gmrt_flagfiles": [],
                    "measurement_set": str(imported_ms.resolve()),
                }
            )
        )

        valid, reason = runner.validate_stage_completion("ingest")

        self.assertFalse(valid)
        self.assertRegex(reason, r"(?i)different.*FITS")

    def test_ingest_resume_checks_checksum_not_only_size_and_mtime(self):
        source = write_random_groups_uvfits(self.root / "raw.uvfits")
        original_stat = source.stat()
        original_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
        runner = self.runner(source)
        imported_ms = make_measurement_set(runner.imported_ms)
        ingest_dir = self.outdir / "ingest"
        (ingest_dir / "ingest.done").touch()
        (ingest_dir / "ingest_report.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "stage": "ingest",
                    "status": "completed",
                    "source_fits": {
                        "path": str(source.resolve()),
                        "size_bytes": original_stat.st_size,
                        "mtime_ns": original_stat.st_mtime_ns,
                        "sha256": original_sha256,
                    },
                    "requested_fits_format": "auto",
                    "fits_format": "uvfits",
                    "casa_task": "importgmrt",
                    "requested_telescope": "ugmrt",
                    "telescope": "ugmrt",
                    "gmrt_flagfiles": [],
                    "measurement_set": str(imported_ms.resolve()),
                }
            )
        )

        changed = bytearray(source.read_bytes())
        changed[-1] = ord("X")
        source.write_bytes(changed)
        os.utime(
            source,
            ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns),
        )
        self.assertEqual(source.stat().st_size, original_stat.st_size)
        self.assertEqual(source.stat().st_mtime_ns, original_stat.st_mtime_ns)

        valid, reason = runner.validate_stage_completion("ingest")

        self.assertFalse(valid)
        self.assertRegex(reason, r"(?i)(checksum|sha-?256|content)")

    def test_dry_run_prints_ingest_before_inspection(self):
        source = write_random_groups_uvfits(self.root / "raw.uvfits")
        runner = self.runner(source)
        output = io.StringIO()

        with contextlib.redirect_stdout(output):
            runner.run()

        text = output.getvalue()
        self.assertLess(text.index("Stage: ingest"), text.index("Stage: inspection"))
        self.assertIn(str(source.resolve()), text)
        self.assertIn(str((self.outdir / "ingest" / "raw.ms").resolve()), text)



def write_tool(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\nexit 0\n")
    path.chmod(0o755)
    return path


class LTAPipelineTests(unittest.TestCase):
    """A GMRT LTA file is planned through the ingest stage like FITS."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.outdir = self.root / "Products"
        self.lta = self.root / "archive" / "28_066_22NOV2020_gwb.lta"
        self.lta.parent.mkdir()
        self.lta.write_bytes(b"GMRT LTA stand-in")
        self.listscan = write_tool(self.root / "gmrt" / "listscan")
        self.gvfits = write_tool(self.root / "gmrt" / "gvfits")

    def tearDown(self):
        self.temporary.cleanup()

    def parse(self, input_path, *extra):
        argv = [
            str(input_path),
            "--outdir", str(self.outdir),
            "--dry-run",
            "--mpi", "off",
            "--no-copy-input-ms",
            "--telescope", "ugmrt",
            "--stages", "inspection",
            "--listscan", str(self.listscan),
            "--gvfits", str(self.gvfits),
            *extra,
        ]
        return PIPELINE.build_parser().parse_args(argv)

    def runner(self, input_path, *extra):
        return PIPELINE.PipelineRunner(self.parse(input_path, *extra))

    def test_lta_plans_ingest_with_the_converters(self):
        runner = self.runner(self.lta)

        self.assertTrue(runner.input_is_lta())
        self.assertEqual(runner.input_kind, "fits")
        self.assertEqual(runner.selected_stages(), ["ingest", "inspection"])
        expected_ms = (self.outdir / "ingest" / "28_066_22NOV2020_gwb.ms").resolve()
        self.assertEqual(runner.imported_ms, expected_ms)
        ingest = runner.command_for_stage("ingest")
        self.assertIn(str(SOURCE_DIR / "fits_ingest.py"), ingest)
        self.assertIn(str(self.lta), ingest)
        self.assertIn(str(expected_ms), ingest)
        self.assertEqual(ingest[ingest.index("--listscan") + 1], str(self.listscan))
        self.assertEqual(ingest[ingest.index("--gvfits") + 1], str(self.gvfits))
        self.assertNotIn("--keep-lta-fits", ingest)
        self.assertIn(str(expected_ms), runner.command_for_stage("inspection"))

        kept = self.runner(self.lta, "--keep-lta-fits").command_for_stage("ingest")
        self.assertIn("--keep-lta-fits", kept)

    def test_imported_ms_names_for_other_lta_files(self):
        for name, expected in [
            ("obs.ltb", "obs_ltb.ms"),
            ("obs.lta_1", "obs_lta_1.ms"),
            ("OBS.LTA", "OBS.ms"),
        ]:
            with self.subTest(name=name):
                source = self.root / "archive" / name
                source.write_bytes(b"x")
                runner = self.runner(source)
                self.assertTrue(runner.input_is_lta())
                self.assertEqual(runner.imported_ms.name, expected)

    def test_lta_not_yet_present_is_still_planned_as_lta(self):
        runner = self.runner(self.root / "later" / "obs.lta")
        self.assertTrue(runner.input_is_lta())
        self.assertEqual(runner.selected_stages(), ["ingest", "inspection"])

    def test_fits_command_is_unchanged_by_lta_support(self):
        source = write_random_groups_uvfits(self.root / "raw.uvfits")
        runner = self.runner(source)
        self.assertFalse(runner.input_is_lta())
        ingest = runner.command_for_stage("ingest")
        self.assertNotIn("--listscan", ingest)
        self.assertNotIn("--gvfits", ingest)

    def test_preflight_accepts_lta_without_reading_it_as_fits(self):
        runner = self.runner(self.lta)
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            runner.preflight(runner.selected_stages())
        self.assertNotIn("WARNING", output.getvalue())

    def test_missing_converters_warn_in_a_dry_run_and_stop_a_real_run(self):
        runner = self.runner(
            self.lta, "--listscan", str(self.root / "nowhere" / "listscan")
        )
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            runner.check_lta_tools()
        self.assertIn("CAPTURE-CASA6", output.getvalue())

        runner.args.dry_run = False
        with self.assertRaisesRegex(FileNotFoundError, "listscan.*not found"):
            runner.check_lta_tools()

    def test_lta_with_another_telescope_is_refused(self):
        runner = self.runner(self.lta, "--telescope", "meerkat")
        with self.assertRaisesRegex(RuntimeError, "GMRT data"):
            runner.check_lta_tools()

    def test_converter_options_do_not_count_as_a_changed_command(self):
        normalize = PIPELINE.PipelineRunner.normalized_checkpoint_resume_command
        first = ["casa", "-c", "fits_ingest.py", "obs.lta", "obs.ms",
                 "--listscan", "/a/listscan", "--gvfits", "/a/gvfits"]
        second = ["casa", "-c", "fits_ingest.py", "obs.lta", "obs.ms",
                  "--listscan", "/b/listscan", "--gvfits", "/b/gvfits",
                  "--keep-lta-fits"]
        self.assertEqual(normalize(first), normalize(second))
        self.assertEqual(normalize(first), ["casa", "-c", "fits_ingest.py", "obs.lta", "obs.ms"])

    def test_lta_run_needs_the_new_ingest_files(self):
        runner = self.runner(self.lta)
        required = runner.required_source_versions(runner.selected_stages())
        self.assertEqual(required["fits_ingest.py"], 1)
        self.assertEqual(required["lta_convert.py"], 2)
        self.assertEqual(PIPELINE.outdated_source_files(SOURCE_DIR, required), [])

        fits_runner = self.runner(write_random_groups_uvfits(self.root / "raw.uvfits"))
        self.assertNotIn(
            "lta_convert.py",
            fits_runner.required_source_versions(fits_runner.selected_stages()),
        )

    def test_ingest_fingerprint_covers_the_converter(self):
        runner = self.runner(self.lta)
        before = runner.stage_fingerprint("ingest")
        runner._stage_fingerprints.clear()
        runner.root = self.root / "copy"
        runner.root.mkdir()
        for name in ["fits_ingest.py", "config.py", "config.yaml",
                     "telescope_profiles.py", "casa_compat.py", "qa_products.py"]:
            (runner.root / name).write_bytes((SOURCE_DIR / name).read_bytes())
        (runner.root / "lta_convert.py").write_text("# another converter\n")
        self.assertNotEqual(before, runner.stage_fingerprint("ingest"))

    def test_generated_sbatch_keeps_the_converters(self):
        sbatch = self.root / "lta_pipeline.sbatch"
        args = self.parse(self.lta, "--keep-lta-fits", "--write-sbatch", str(sbatch))
        with contextlib.redirect_stdout(io.StringIO()):
            PIPELINE.write_sbatch_from_config(args)
        content = sbatch.read_text()
        self.assertIn(str(self.lta), content)
        self.assertIn(f"--listscan {self.listscan}", content)
        self.assertIn(f"--gvfits {self.gvfits}", content)
        self.assertIn("--keep-lta-fits", content)

    def test_completed_lta_ingest_is_accepted_on_resume(self):
        runner = self.runner(self.lta)
        imported_ms = make_measurement_set(runner.imported_ms)
        ingest_dir = self.outdir / "ingest"
        (ingest_dir / "ingest.done").touch()
        lta_stat = self.lta.stat()
        (ingest_dir / "ingest_report.json").write_text(json.dumps({
            "schema_version": 1,
            "stage": "ingest",
            "status": "completed",
            "source_fits": {
                "path": str(self.lta),
                "size_bytes": lta_stat.st_size,
                "mtime_ns": lta_stat.st_mtime_ns,
                "sha256": hashlib.sha256(self.lta.read_bytes()).hexdigest(),
            },
            "source_format": "lta",
            "requested_fits_format": "auto",
            "fits_format": "uvfits",
            "casa_task": "importgmrt",
            "requested_telescope": "ugmrt",
            "telescope": "ugmrt",
            "gmrt_flagfiles": [],
            "measurement_set": str(imported_ms.resolve()),
        }))

        valid, reason = runner.validate_stage_completion("ingest")
        self.assertTrue(valid, reason)

        self.lta.write_bytes(b"GMRT LTA stand-iX")
        os.utime(self.lta, ns=(lta_stat.st_atime_ns, lta_stat.st_mtime_ns))
        valid, reason = runner.validate_stage_completion("ingest")
        self.assertFalse(valid)
        self.assertRegex(reason, r"(?i)sha-?256")


if __name__ == "__main__":
    unittest.main()
