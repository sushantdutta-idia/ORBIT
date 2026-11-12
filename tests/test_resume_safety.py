import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_pipeline_module():
    path = SOURCE_DIR / "pipeline.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_pipeline_resume_safety",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PIPELINE = load_pipeline_module()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2))


def make_casa_table(path):
    path.mkdir(parents=True, exist_ok=True)
    (path / "table.dat").touch()
    return path


def make_measurement_set(path):
    make_casa_table(path)
    for name in [
        "ANTENNA",
        "DATA_DESCRIPTION",
        "OBSERVATION",
        "POLARIZATION",
        "SPECTRAL_WINDOW",
    ]:
        subtable = make_casa_table(path / name)
        (subtable / "table.dat").write_bytes(name.encode("ascii"))
    return path


class PipelineResumeSafetyTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.ms = self.root / "input.ms"
        make_casa_table(self.ms)

    def tearDown(self):
        self.temporary.cleanup()

    def runner(self, outdir=None, *extra):
        outdir = outdir or self.root / "Products"
        argv = [
            str(self.ms),
            "--outdir",
            str(outdir),
            "--dry-run",
            "--mpi",
            "off",
            "--no-copy-input-ms",
            *extra,
        ]
        args = PIPELINE.build_parser().parse_args(argv)
        return PIPELINE.PipelineRunner(args)

    def prepare_meerkat_calibration_completion(
        self,
        report_updates=None,
        telescope="meerkat",
        receiver_band="L",
    ):
        checkpoint_prefix = telescope
        runner = self.runner(
            None,
            "--telescope",
            telescope,
            "--receiver-band",
            receiver_band,
            "--data-state",
            "raw",
            "--stages",
            "calibration",
        )
        caldir = runner.outdir / "calibration"
        caldir.mkdir(parents=True, exist_ok=True)
        for checkpoint in [
            "calibration.done",
            f"{checkpoint_prefix}_postcal_flagging.done",
            f"{checkpoint_prefix}_round2_complete.done",
        ]:
            (caldir / checkpoint).touch()

        write_json(
            caldir / f"{checkpoint_prefix}_postcal_flagging_report.json",
            {
                "status": "completed",
                "correlation": "ABS_XX,YY",
                "parameters": {
                    "datacolumn": "CORRECTED",
                    "ntime": "scan",
                    "combinescans": False,
                    "extendflags": False,
                },
            },
        )
        write_json(
            caldir / f"{checkpoint_prefix}_round2_report.json",
            {"status": "completed"},
        )
        for directory in [caldir, caldir / "round1"]:
            for table in [
                "delay.cal",
                "prebp.cal",
                "bandpass.cal",
                "phase.cal",
                "phase_scan.cal",
                "amp.cal",
                "flux.cal",
            ]:
                make_casa_table(directory / table)

        target_ms = make_measurement_set(caldir / "targets" / "target.ms")
        write_json(
            runner.report_path("split_targets"),
            {"target": str(target_ms)},
        )
        target_instance_identity = (
            PIPELINE.IMAGING_SAFETY.measurement_set_instance_identity(
                target_ms
            )
        )
        report = {
            "measurement_set": str(self.ms),
            "measurement_set_used": str(self.ms),
            "measurement_set_stat_sha256": (
                PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(self.ms)
            ),
            "original_measurement_set": str(self.ms),
            "telescope": telescope,
            "band": (
                f"MEERKAT_{receiver_band}"
                if telescope == "meerkat"
                else receiver_band
            ),
            "calibration_parameters": {
                "solution_rounds": 2,
                "postcal_flagging": True,
                "prebp_gaintype": "G",
            },
            "checkpoints": {
                f"{checkpoint_prefix}_postcal_flagging": True,
                f"{checkpoint_prefix}_round2_complete": True,
            },
            "split_target_identities": {
                "schema_version": 1,
                "targets": {
                    "target": {
                        "path": str(target_ms.resolve()),
                        "instance_identity": target_instance_identity,
                    },
                },
            },
        }
        report["provenance_fingerprint"] = runner.stage_fingerprint(
            "calibration"
        )
        if report_updates:
            for key, value in report_updates.items():
                if isinstance(value, dict) and isinstance(report.get(key), dict):
                    report[key].update(value)
                else:
                    report[key] = value
        write_json(runner.report_path("calibration"), report)
        return runner, report

    def test_ugmrt_calibration_resume_requires_two_round_evidence(self):
        runner, _ = self.prepare_meerkat_calibration_completion(
            telescope="ugmrt",
            receiver_band="B3",
        )

        valid, reason = runner.validate_stage_completion("calibration")

        self.assertTrue(valid, reason)

    def rewrite_postcal_parameters(self, runner, telescope, **updates):
        path = (
            runner.outdir / "calibration" /
            f"{telescope}_postcal_flagging_report.json"
        )
        report = json.loads(path.read_text())
        report["parameters"].update(updates)
        path.write_text(json.dumps(report))

    def test_calibration_accepts_recorded_memory_bounded_postcal_chunk(self):
        runner, _ = self.prepare_meerkat_calibration_completion(
            telescope="ugmrt",
            receiver_band="B4",
        )
        self.rewrite_postcal_parameters(
            runner,
            "ugmrt",
            ntime="919s",
            flagdata_chunk_plan={"status": "bounded", "ntime": "919s"},
        )

        valid, reason = runner.validate_stage_completion("calibration")

        self.assertTrue(valid, reason)

    def test_calibration_requires_scan_local_flag_extension(self):
        runner, _ = self.prepare_meerkat_calibration_completion(
            telescope="ugmrt",
            receiver_band="B4",
        )
        self.rewrite_postcal_parameters(
            runner,
            "ugmrt",
            extend={
                "applied": True,
                "combinescans": False,
                "growtime": 80.0,
                "growfreq": 80.0,
            },
        )
        valid, reason = runner.validate_stage_completion("calibration")
        self.assertTrue(valid, reason)

        self.rewrite_postcal_parameters(
            runner,
            "ugmrt",
            extend={"applied": True, "combinescans": True},
        )
        valid, reason = runner.validate_stage_completion("calibration")
        self.assertFalse(valid)
        self.assertIn("extension is not scan-local", reason)

    def test_calibration_rejects_unexplained_or_cross_scan_postcal_chunk(self):
        cases = [
            {"ntime": "919s"},
            {
                "ntime": "919s",
                "flagdata_chunk_plan": {"status": "scan_fits", "ntime": "scan"},
            },
            {
                "ntime": "919s",
                "flagdata_chunk_plan": {"status": "bounded", "ntime": "600s"},
            },
            {"ntime": "0s", "flagdata_chunk_plan": {"status": "bounded",
                                                    "ntime": "0s"}},
            {"ntime": "scan", "combinescans": True},
        ]
        for updates in cases:
            with self.subTest(updates=updates):
                runner, _ = self.prepare_meerkat_calibration_completion(
                    telescope="ugmrt",
                    receiver_band="B4",
                )
                self.rewrite_postcal_parameters(runner, "ugmrt", **updates)

                valid, reason = runner.validate_stage_completion(
                    "calibration"
                )

                self.assertFalse(valid)
                self.assertIn("scan-local", reason)

    def test_calibration_resume_rejects_replaced_split_target_ms(self):
        runner, _ = self.prepare_meerkat_calibration_completion()
        split_report = json.loads(
            runner.report_path("split_targets").read_text()
        )
        target_ms = Path(split_report["target"])
        target_ms.rename(target_ms.with_name("target.original.ms"))
        make_measurement_set(target_ms)

        valid, reason = runner.validate_stage_completion("calibration")

        self.assertFalse(valid)
        self.assertRegex(reason.lower(), r"target|replaced|identity|bound")

    def test_calibration_resume_rebuilds_only_replaced_split_target_ms(self):
        runner, _ = self.prepare_meerkat_calibration_completion()
        split_report = json.loads(
            runner.report_path("split_targets").read_text()
        )
        target_ms = Path(split_report["target"])
        target_ms.rename(target_ms.with_name("target.original.ms"))
        make_measurement_set(target_ms)
        runner.args.resume = True
        current_recipe = runner.stage_recipe_fingerprint("calibration")
        runner.resume_state = {
            "schema_version": 1,
            "created_utc": PIPELINE.utc_now(),
            "stages": {
                "calibration": {
                    "status": "completed",
                    "product_stage_fingerprint": "0" * 64,
                    "product_recipe_fingerprint": current_recipe,
                    "command": runner.command_for_stage("calibration"),
                },
            },
        }

        self.assertFalse(runner.should_skip_stage("calibration"))
        entry = runner.resume_state["stages"]["calibration"]
        self.assertEqual(entry["status"], "stale")
        self.assertIn("target", entry["validation_error"].lower())

    def test_meerkat_calibration_resume_rejects_a_one_round_report(self):
        runner, report = self.prepare_meerkat_calibration_completion(
            {"calibration_parameters": {"solution_rounds": 1}}
        )

        valid, reason = runner.validate_stage_completion("calibration")

        self.assertFalse(valid)
        self.assertRegex(reason.lower(), r"round|solution")

    def test_meerkat_calibration_resume_requires_reported_postcal_checkpoints(self):
        for checkpoint in [
            "meerkat_postcal_flagging",
            "meerkat_round2_complete",
        ]:
            with self.subTest(checkpoint=checkpoint):
                runner, report = self.prepare_meerkat_calibration_completion()
                report["checkpoints"][checkpoint] = False
                write_json(runner.report_path("calibration"), report)

                valid, reason = runner.validate_stage_completion("calibration")

                self.assertFalse(valid)
                self.assertRegex(reason.lower(), r"postcal|round|flag|checkpoint")

    def test_meerkat_calibration_resume_requires_durable_postcal_markers(self):
        for checkpoint in [
            "meerkat_postcal_flagging",
            "meerkat_round2_complete",
        ]:
            with self.subTest(checkpoint=checkpoint):
                runner, _ = self.prepare_meerkat_calibration_completion()
                (runner.outdir / "calibration" / f"{checkpoint}.done").unlink()

                valid, reason = runner.validate_stage_completion("calibration")

                self.assertFalse(valid)
                self.assertIn(checkpoint, reason)

    def prepare_meerkat_imaging_completion(self, band):
        canonical_band = f"MEERKAT_{band}"
        model = {
            "UHF": "MKAT-AA-UHF-JIM-2020",
            "L": "MKAT-AA-L-JIM-2020",
        }[band]
        outdir = self.root / f"Products-{band}"
        runner = self.runner(
            outdir,
            "--telescope",
            "meerkat",
            "--receiver-band",
            band,
            "--data-state",
            "raw",
            "--target",
            "target_b",
            "--stages",
            "imaging",
            # These tests validate the CASA backend's katbeam contract.
            "--imager",
            "casa",
        )
        image_root = outdir / "imaging" / "images" / "target_b"
        restored = make_casa_table(Path(str(image_root) + ".image.tt0"))
        pbcor = make_casa_table(Path(str(image_root) + ".image.tt0.pbcor"))
        katbeam_pb = make_casa_table(
            Path(str(image_root) + ".katbeam.pb.tt0")
        )
        fits = outdir / "imaging" / "fits" / "target_b.pbcor.fits"
        fits.parent.mkdir(parents=True, exist_ok=True)
        fits.touch()
        checkpoint = outdir / "imaging" / "checkpoints" / "imaging.done"
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        checkpoint.touch()

        report = {
            "measurement_set": str(self.ms),
            "measurement_set_stat_sha256": (
                PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(self.ms)
            ),
            "target": "target_b",
            "telescope": "meerkat",
            "band": canonical_band,
            "imaging_backend": "casa",
            "imaging_parameters": {
                "primary_beam_provider": "katbeam",
                "katbeam_model": model,
                "imsize": 4096,
                "cell": "1.5arcsec",
                "pbmin": 0.3,
            },
            "products": {
                "restored_image": str(restored),
                "pbcor_image": str(pbcor),
                "pbcor_method": "katbeam_impbcor_tt0",
                "pbcor_status": {
                    "katbeam_model": model,
                    "katbeam_pb_image": str(katbeam_pb),
                    "katbeam_frequency_mhz": 1284.0,
                    "intensity_taylor_term": 0,
                    "spectral_terms_corrected": False,
                },
                "fits": [str(fits)],
                "fits_content_identities": {
                    str(fits.resolve()): (
                        PIPELINE.IMAGING_SAFETY.file_content_identity(fits)
                    ),
                },
            },
        }
        report["provenance_fingerprint"] = runner.stage_fingerprint(
            "imaging"
        )
        write_json(runner.report_path("imaging"), report)
        provenance = Path(str(katbeam_pb) + ".provenance.json")
        return runner, report, provenance, katbeam_pb, restored

    def write_matching_katbeam_provenance(
        self,
        runner,
        report,
        provenance,
        katbeam_pb,
        restored,
    ):
        parameters = report["imaging_parameters"]
        products = report["products"]
        pbcor = Path(products["pbcor_image"])
        geometry = {
            "shape": [4096, 4096, 1, 1],
            "pixel_axes": [0, 1],
            "direction_increment_rad": [-1.0e-5, 1.0e-5],
            "direction_reference_pixel": [2047.0, 2047.0],
            "direction_reference_value_rad": [0.1, -0.2],
            "coordsys_sha256": "test-coordinate-record",
            "frequency_hz": 1.284e9,
        }
        write_json(
            provenance,
            {
                "schema_version": 2,
                "provider": "katbeam",
                "katbeam_version": "test-version",
                "model": parameters["katbeam_model"],
                "frequency_mhz": 1284.0,
                "source_image": str(restored.resolve()),
                "source_stat_sha256": (
                    PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(restored)
                ),
                "source_geometry": geometry,
                "imsize": parameters["imsize"],
                "cell": parameters["cell"],
                "pb_image": str(katbeam_pb.resolve()),
                "pb_stat_sha256": (
                    PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(katbeam_pb)
                ),
                "pb_geometry": geometry,
            },
        )
        corrected_provenance = Path(str(pbcor) + ".provenance.json")
        write_json(
            corrected_provenance,
            {
                "schema_version": 1,
                "provider": "katbeam",
                "method": "katbeam_impbcor_tt0",
                "model": parameters["katbeam_model"],
                "source_image": str(restored.resolve()),
                "source_stat_sha256": (
                    PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(restored)
                ),
                "pb_image": str(katbeam_pb.resolve()),
                "pb_stat_sha256": (
                    PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(katbeam_pb)
                ),
                "pbcor_image": str(pbcor.resolve()),
                "pbcor_stat_sha256": (
                    PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(pbcor)
                ),
                "pbmin": parameters["pbmin"],
            },
        )
        products["pbcor_status"]["pbcor_provenance"] = str(
            corrected_provenance
        )
        write_json(runner.report_path("imaging"), report)
        return corrected_provenance

    def test_meerkat_imaging_resume_requires_pbcor_product(self):
        runner, report, provenance, katbeam_pb, restored = (
            self.prepare_meerkat_imaging_completion("L")
        )
        report["products"]["pbcor_image"] = None
        write_json(runner.report_path("imaging"), report)

        valid, reason = runner.validate_stage_completion("imaging")

        self.assertFalse(valid)
        self.assertIn("pbcor", reason.lower())

    def test_meerkat_uhf_and_l_resume_require_katbeam_pb_provenance(self):
        for band in ["UHF", "L"]:
            with self.subTest(band=band):
                runner, report, provenance, katbeam_pb, restored = (
                    self.prepare_meerkat_imaging_completion(band)
                )
                self.assertFalse(provenance.exists())

                valid, reason = runner.validate_stage_completion("imaging")

                self.assertFalse(valid)
                self.assertRegex(reason.lower(), r"katbeam|provenance")

    def test_meerkat_imaging_resume_accepts_matching_katbeam_provenance(self):
        runner, report, provenance, katbeam_pb, restored = (
            self.prepare_meerkat_imaging_completion("L")
        )
        self.write_matching_katbeam_provenance(
            runner,
            report,
            provenance,
            katbeam_pb,
            restored,
        )

        valid, reason = runner.validate_stage_completion("imaging")

        self.assertTrue(valid, reason)

    def test_meerkat_imaging_validation_binds_images_by_content(self):
        runner, report, provenance, katbeam_pb, restored = (
            self.prepare_meerkat_imaging_completion("L")
        )
        corrected_provenance = self.write_matching_katbeam_provenance(
            runner,
            report,
            provenance,
            katbeam_pb,
            restored,
        )
        safety = PIPELINE.IMAGING_SAFETY
        pbcor = Path(json.loads(corrected_provenance.read_text())["pbcor_image"])
        for path, prefixes in (
            (provenance, [("source", restored), ("pb", katbeam_pb)]),
            (
                corrected_provenance,
                [("source", restored), ("pb", katbeam_pb), ("pbcor", pbcor)],
            ),
        ):
            contents = json.loads(path.read_text())
            for prefix, artifact in prefixes:
                contents[f"{prefix}_content_identity"] = (
                    safety.casa_table_content_identity(artifact)
                )
            write_json(path, contents)

        # Timestamps read back differently, as on the network filesystem.
        for item in katbeam_pb.rglob("*"):
            if item.is_file():
                stat = item.stat()
                os.utime(item, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10**9))
        valid, reason = runner.validate_stage_completion("imaging")
        self.assertTrue(valid, reason)

        (katbeam_pb / "table.dat").write_bytes(b"rewritten beam")
        valid, reason = runner.validate_stage_completion("imaging")
        self.assertFalse(valid)
        self.assertIn("rewritten table.dat", reason)

    def test_meerkat_imaging_resume_rejects_mismatched_katbeam_provenance(self):
        runner, report, provenance, katbeam_pb, restored = (
            self.prepare_meerkat_imaging_completion("L")
        )
        self.write_matching_katbeam_provenance(
            runner,
            report,
            provenance,
            katbeam_pb,
            restored,
        )
        contents = json.loads(provenance.read_text())
        contents["model"] = "a-different-katbeam-model"
        write_json(provenance, contents)

        valid, reason = runner.validate_stage_completion("imaging")

        self.assertFalse(valid)
        self.assertRegex(reason.lower(), r"katbeam|model|provenance")

    def test_meerkat_imaging_resume_requires_katbeam_pb_artifact(self):
        runner, report, provenance, katbeam_pb, restored = (
            self.prepare_meerkat_imaging_completion("L")
        )
        self.write_matching_katbeam_provenance(
            runner,
            report,
            provenance,
            katbeam_pb,
            restored,
        )
        (katbeam_pb / "table.dat").unlink()

        valid, reason = runner.validate_stage_completion("imaging")

        self.assertFalse(valid)
        self.assertRegex(reason.lower(), r"katbeam|beam|artifact")

    def test_meerkat_imaging_resume_rejects_mutated_image_tables(self):
        for artifact_name in ["restored", "beam", "pbcor"]:
            with self.subTest(artifact=artifact_name):
                runner, report, provenance, katbeam_pb, restored = (
                    self.prepare_meerkat_imaging_completion("L")
                )
                self.write_matching_katbeam_provenance(
                    runner,
                    report,
                    provenance,
                    katbeam_pb,
                    restored,
                )
                artifacts = {
                    "restored": restored,
                    "beam": katbeam_pb,
                    "pbcor": Path(report["products"]["pbcor_image"]),
                }
                (artifacts[artifact_name] / "mutated-pixels").write_text(
                    "changed"
                )

                valid, reason = runner.validate_stage_completion("imaging")

                self.assertFalse(valid)
                self.assertRegex(reason.lower(), r"match|provenance|fingerprint")

    def prepare_source_detection_completion(self):
        runner = self.runner(
            None,
            "--telescope",
            "meerkat",
            "--receiver-band",
            "L",
            "--data-state",
            "raw",
            "--target",
            "target_b",
            "--stages",
            "source_detection",
        )
        root = runner.outdir / "source_detection"
        checkpoint = root / "checkpoints" / "source_detection.done"
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        checkpoint.touch()
        products = {}
        for key, suffix in [
            ("source_catalog_csv", "catalog.csv"),
            ("source_catalog_fits", "catalog.fits"),
            ("source_regions", "sources.reg"),
        ]:
            path = root / suffix
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
            products[key] = str(path)
        report = {
            "status": "success",
            "target": "target_b",
            "telescope": "meerkat",
            "band": "MEERKAT_L",
            "imaging_backend": "casa",
            "input_image": None,
            "parameters": {"image_preference": "pbcor"},
            "primary_beam": {"method": "katbeam_impbcor_tt0"},
            "products": products,
            "product_content_identities": {
                str(Path(value).resolve()): (
                    PIPELINE.IMAGING_SAFETY.file_content_identity(value)
                )
                for value in products.values()
            },
        }
        source_image = (
            runner.outdir / "imaging" / "fits" / "target_b.pbcor.fits"
        )
        source_image.parent.mkdir(parents=True, exist_ok=True)
        source_image.touch()
        report["input_image"] = str(source_image)
        write_json(
            runner.report_path("imaging"),
            {
                "target": "target_b",
                "telescope": "meerkat",
                "band": "MEERKAT_L",
                "imaging_backend": "casa",
                "products": {
                    "restored_image": str(
                        runner.outdir /
                        "imaging" /
                        "fits" /
                        "target_b.restored.fits"
                    ),
                    "pbcor_image": str(source_image),
                    "pbcor_method": "katbeam_impbcor_tt0",
                    "fits": [str(source_image)],
                },
            },
        )
        write_json(runner.report_path("source_detection"), report)
        return runner, report

    def test_source_detection_resume_requires_matching_target(self):
        runner, report = self.prepare_source_detection_completion()
        report["target"] = "target_a"
        write_json(runner.report_path("source_detection"), report)

        valid, reason = runner.validate_stage_completion("source_detection")

        self.assertFalse(valid)
        self.assertIn("target", reason.lower())

    def test_source_detection_resume_requires_matching_telescope(self):
        runner, report = self.prepare_source_detection_completion()
        report["telescope"] = "lofar"
        write_json(runner.report_path("source_detection"), report)

        valid, reason = runner.validate_stage_completion("source_detection")

        self.assertFalse(valid)
        self.assertIn("telescope", reason.lower())

    def test_source_detection_resume_requires_matching_receiver_band(self):
        runner, report = self.prepare_source_detection_completion()
        report["band"] = "MEERKAT_UHF"
        write_json(runner.report_path("source_detection"), report)

        valid, reason = runner.validate_stage_completion("source_detection")

        self.assertFalse(valid)
        self.assertIn("band", reason.lower())

    def test_source_detection_resume_accepts_matching_context(self):
        runner, _ = self.prepare_source_detection_completion()

        valid, reason = runner.validate_stage_completion("source_detection")

        self.assertTrue(valid, reason)

    def test_raw_lofar_imaging_is_rejected_before_execution(self):
        runner = self.runner(
            None,
            "--telescope",
            "lofar",
            "--receiver-band",
            "HBA",
            "--data-state",
            "raw",
            "--imager",
            "auto",
            "--lofar-beam-mode",
            "full",
            "--stages",
            "imaging",
        )

        with self.assertRaisesRegex(ValueError, "Raw LOFAR"):
            runner.preflight(["imaging"])

    def test_target_b_imaging_does_not_reuse_target_a_selfcal_ms(self):
        runner = self.runner(
            None,
            "--telescope",
            "meerkat",
            "--receiver-band",
            "L",
            "--data-state",
            "raw",
            "--target",
            "target_b",
            "--stages",
            "imaging",
        )
        target_a_selfcal = make_casa_table(
            runner.outdir / "selfcal" / "ms" / "target_a.final.ms"
        )
        target_a_split = make_casa_table(
            runner.outdir / "calibration" / "targets" / "target_a.ms"
        )
        target_b_split = make_casa_table(
            runner.outdir / "calibration" / "targets" / "target_b.ms"
        )
        write_json(
            runner.report_path("selfcal"),
            {
                "target": "target_a",
                "telescope": "meerkat",
                "band": "MEERKAT_L",
                "final_ms": str(target_a_selfcal),
            },
        )
        write_json(
            runner.report_path("split_targets"),
            {
                "target_a": str(target_a_split),
                "target_b": str(target_b_split),
            },
        )

        selected = Path(runner.resolve_imaging_ms()).resolve()
        command = runner.command_for_stage("imaging")

        self.assertEqual(selected, target_b_split.resolve())
        self.assertIn(str(target_b_split), command)
        self.assertNotIn(str(target_a_selfcal), command)

    def test_selfcal_resume_binds_the_calibrated_target_input_and_recipe(self):
        runner = self.runner(
            None,
            "--telescope", "meerkat",
            "--receiver-band", "L",
            "--data-state", "raw",
            "--target", "target_a",
            "--stages", "selfcal",
        )
        target = make_casa_table(
            runner.outdir / "calibration" / "targets" / "target_a.ms"
        )
        final_ms = make_casa_table(
            runner.outdir / "selfcal" / "selfcal_ms" / "final.ms"
        )
        final_image = make_casa_table(
            runner.outdir / "selfcal" / "images" / "final.image.tt0"
        )
        write_json(runner.report_path("split_targets"), {
            "target_a": str(target),
        })
        (runner.outdir / "selfcal" / "checkpoints").mkdir(
            parents=True, exist_ok=True
        )
        (runner.outdir / "selfcal" / "checkpoints" / "selfcal.done").touch()
        report = {
            "telescope": "meerkat",
            "band": "MEERKAT_L",
            "target": "target_a",
            "measurement_set": str(target),
            "input_measurement_set": str(target),
            "input_measurement_set_stat_sha256": (
                PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(target)
            ),
            "final_ms": str(final_ms),
            "final_ms_stat_sha256": (
                PIPELINE.IMAGING_SAFETY.path_stat_fingerprint(final_ms)
            ),
            "final_image": str(final_image),
            "provenance_fingerprint": runner.stage_fingerprint("selfcal"),
        }
        write_json(runner.report_path("selfcal"), report)

        valid, reason = runner.validate_stage_completion("selfcal")
        self.assertTrue(valid, reason)

        report["input_measurement_set"] = str(
            runner.outdir / "calibration" / "targets" / "other.ms"
        )
        write_json(runner.report_path("selfcal"), report)
        valid, reason = runner.validate_stage_completion("selfcal")
        self.assertFalse(valid)
        self.assertIn("input MS", reason)


if __name__ == "__main__":
    unittest.main()
