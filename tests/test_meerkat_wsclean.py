"""MeerKAT imaging with WSClean, end to end, and the telescope defaults."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from astropy.io import fits

from test_meerkat_beam import (
    BEAM, FakeJimBeam, NX, NY, geometry, wsclean_header, write_image,
)
from test_resume_safety import PIPELINE, make_casa_table, write_json
from test_wsclean_safety import WSCLEAN


MODEL = "MKAT-AA-L-JIM-2020"
# MeerKAT L-band: 856-1712 MHz in eight equal output channels.
EDGES = np.linspace(856.0, 1712.0, 9)
CENTRES = 0.5 * (EDGES[:-1] + EDGES[1:])
WIDTH = float(EDGES[1] - EDGES[0])
FLUX = 1.7


class FakeSubprocessResult:
    returncode = 0


def write_wsclean_products(prefix):
    """Write what a joined-channel WSClean run leaves behind."""

    value = geometry()
    beam = FakeJimBeam("m")
    for index, freq in enumerate(CENTRES):
        plane = FLUX * BEAM.beam_plane(beam, value, freq)
        write_image(Path(f"{prefix}-{index:04d}-image.fits"), plane,
                    wsclean_header(freq * 1e6, WIDTH * 1e6))
        write_image(Path(f"{prefix}-{index:04d}-dirty.fits"), plane,
                    wsclean_header(freq * 1e6, WIDTH * 1e6))
    mean = np.mean([BEAM.beam_plane(beam, value, f) for f in CENTRES], axis=0)
    for kind in ("image", "residual", "model", "dirty"):
        write_image(Path(f"{prefix}-MFS-{kind}.fits"), FLUX * mean,
                    wsclean_header(1284e6, 856e6))
    write_image(Path(f"{prefix}-psf.fits"), mean,
                wsclean_header(1284e6, 856e6))


class MeerkatFixture(unittest.TestCase):
    """A run root with a self-calibrated MeerKAT target ready to image."""

    TARGET = "target_b"

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.outdir = self.root / "Products"
        self.raw = self.root / "raw.ms"
        make_casa_table(self.raw)
        (self.raw / "table.dat").write_bytes(b"DATA FLAG")
        self.final = self.outdir / "selfcal" / "selfcal_ms" / "final_selfcal.ms"
        make_casa_table(self.final)
        (self.final / "table.dat").write_bytes(b"DATA FLAG WEIGHT")
        self.split = self.outdir / "calibration" / "targets" / f"{self.TARGET}.ms"
        make_casa_table(self.split)
        (self.split / "table.dat").write_bytes(b"DATA CORRECTED_DATA FLAG")

        write_json(self.outdir / "inspection" / "inspection_report.json", {
            "measurement_set": str(self.raw),
            "telescope": {"name": "meerkat"},
            "band": {"name": "MEERKAT_L"},
            "data_columns": ["DATA", "FLAG"],
        })
        write_json(self.outdir / "inspection" / "fields.json", [
            {"id": 0, "name": "calibrator"}, {"id": 5, "name": self.TARGET},
        ])
        write_json(self.outdir / "calibration" / "calibration_report.json", {
            "telescope": "meerkat", "band": "MEERKAT_L",
            "targets": [self.TARGET],
        })
        write_json(self.outdir / "calibration" / "split_targets.json",
                   {self.TARGET: str(self.split)})
        write_json(self.outdir / "selfcal" / "reports" / "selfcal_report.json", {
            "telescope": "meerkat", "band": "MEERKAT_L",
            "target": self.TARGET, "final_ms": str(self.final),
        })

        executable = self.root / "wsclean"
        executable.write_text("#!/bin/sh\n")
        self.executable = str(executable)

    def tearDown(self):
        self.temporary.cleanup()

    def manager(self, ms=None, **kwargs):
        defaults = dict(
            outdir=self.outdir, target=None, telescope="meerkat",
            receiver_band="L", executable=self.executable, threads=4,
            profile="science", beam_mode=None, datacolumn="auto",
            dry_run=False, provenance_fingerprint="f" * 64,
        )
        defaults.update(kwargs)
        return WSCLEAN.WSCleanImagingManager(ms or self.final, **defaults)

    def run_manager(self, manager):
        def fake_wsclean(command, check=False):
            write_wsclean_products(manager.prefix)
            return FakeSubprocessResult()

        with mock.patch.object(WSCLEAN.subprocess, "run", fake_wsclean), \
                mock.patch.object(WSCLEAN, "_MEERKAT_BEAM", BEAM), \
                mock.patch.object(
                    BEAM, "load_katbeam",
                    return_value=(FakeJimBeam, "0.1.2")), \
                mock.patch("builtins.print"):
            return manager.run()


class ManagerInputTests(MeerkatFixture):

    def test_the_self_calibrated_ms_is_accepted_without_being_inspected(self):
        manager = self.manager()
        self.assertEqual(manager.telescope, "meerkat")
        self.assertEqual(manager.beam_mode, "katbeam")
        self.assertEqual(manager.target, self.TARGET)
        self.assertFalse(manager.inspected_input)

    def test_a_pipeline_made_ms_is_one_field_so_no_field_is_selected(self):
        # The inspection report numbers the fields of the raw observation,
        # where this target is field 5. Selecting it in the split MS, which
        # holds one field, would select nothing.
        manager = self.manager()
        self.assertIsNone(manager.field_id)
        self.assertNotIn("-field", manager.build_command())

    def test_the_data_column_is_chosen_as_the_casa_backend_chooses_it(self):
        self.assertEqual(self.manager(self.final).datacolumn, "DATA")
        # Any other MS uses CORRECTED_DATA when it declares one.
        self.assertEqual(self.manager(self.split).datacolumn, "CORRECTED_DATA")
        (self.split / "table.dat").write_bytes(b"DATA FLAG")
        self.assertEqual(self.manager(self.split).datacolumn, "DATA")

    def test_an_explicit_column_must_exist(self):
        self.assertEqual(
            self.manager(self.final, datacolumn="DATA").datacolumn, "DATA")
        with self.assertRaisesRegex(RuntimeError, "not present"):
            self.manager(self.final, datacolumn="CORRECTED_DATA")

    def test_an_ms_the_pipeline_did_not_make_is_refused(self):
        stranger = self.root / "elsewhere.ms"
        make_casa_table(stranger)
        with self.assertRaisesRegex(RuntimeError, "neither the Measurement Set"):
            self.manager(stranger)

    def test_a_target_that_contradicts_the_ms_is_refused(self):
        with self.assertRaisesRegex(RuntimeError, "does not match the target"):
            self.manager(target="someone_else")

    def test_reports_come_from_the_run_root_for_a_per_target_scope(self):
        scope = self.outdir / "targets" / self.TARGET
        write_json(scope / "selfcal" / "reports" / "selfcal_report.json", {
            "target": self.TARGET, "final_ms": str(self.final),
        })
        manager = self.manager(outdir=scope, upstream_outdir=self.outdir)
        self.assertEqual(manager.target, self.TARGET)
        self.assertEqual(manager.report_path.parent.parent.parent.resolve(),
                         scope.resolve())


class BeamModeTests(MeerkatFixture):

    def test_meerkat_defaults_to_katbeam_and_has_no_lofar_modes(self):
        self.assertEqual(self.manager().beam_mode, "katbeam")
        self.assertEqual(self.manager(beam_mode="none").beam_mode, "none")
        for mode in ("full", "differential"):
            with self.assertRaisesRegex(ValueError, "LOFAR EveryBeam"):
                self.manager(beam_mode=mode)

    def test_ugmrt_is_never_imaged_here(self):
        write_json(self.outdir / "inspection" / "inspection_report.json", {
            "measurement_set": str(self.final),
            "telescope": {"name": "ugmrt"}, "band": {"name": "B4"},
            "data_columns": ["DATA"],
        })
        for report in ("calibration", "selfcal"):
            path = (self.outdir / "calibration" / "calibration_report.json"
                    if report == "calibration" else
                    self.outdir / "selfcal" / "reports" / "selfcal_report.json")
            data = json.loads(path.read_text())
            data["telescope"], data["band"] = "ugmrt", "B4"
            write_json(path, data)
        with self.assertRaisesRegex(RuntimeError, "tclean"):
            self.manager(telescope="ugmrt", receiver_band="B4")

    def test_a_band_without_a_measured_beam_is_refused(self):
        for report in ("calibration", "selfcal"):
            path = (self.outdir / "calibration" / "calibration_report.json"
                    if report == "calibration" else
                    self.outdir / "selfcal" / "reports" / "selfcal_report.json")
            data = json.loads(path.read_text())
            data["band"] = "MEERKAT_S0"
            write_json(path, data)
        inspection = self.outdir / "inspection" / "inspection_report.json"
        data = json.loads(inspection.read_text())
        data["band"] = {"name": "MEERKAT_S0"}
        write_json(inspection, data)
        with self.assertRaisesRegex(RuntimeError, "No validated MeerKAT"):
            self.manager(receiver_band="S0")


class CommandTests(MeerkatFixture):

    def test_the_command_images_the_band_with_joined_channels(self):
        command = self.manager().build_command()
        joined = command.index("-join-channels")
        self.assertEqual(command[command.index("-channels-out") + 1], "8")
        self.assertEqual(command[command.index("-fit-spectral-pol") + 1], "2")
        self.assertGreater(command.index("-channels-out"), joined)
        self.assertEqual(command[command.index("-data-column") + 1], "DATA")
        self.assertEqual(command[command.index("-weight") + 2], "-0.5")
        self.assertIn("-multiscale", command)
        self.assertEqual(command[command.index("-size") + 1], "6144")

    def test_the_beam_is_not_left_to_everybeam(self):
        command = self.manager().build_command()
        self.assertNotIn("-apply-primary-beam", command)
        self.assertNotIn("-primary-beam-limit", command)

    def test_a_job_with_a_memory_limit_caps_wsclean(self):
        snapshot = {"job_limit": {"limit_bytes": 28 * 1024 ** 3}}
        fake = mock.Mock(snapshot=lambda: snapshot)
        with mock.patch.object(WSCLEAN, "load_local_module",
                               return_value=fake):
            command = self.manager().build_command()
        self.assertEqual(command[command.index("-abs-mem") + 1], "22")

    def test_without_a_job_limit_nothing_is_added(self):
        fake = mock.Mock(snapshot=lambda: {"job_limit": None})
        with mock.patch.object(WSCLEAN, "load_local_module",
                               return_value=fake):
            self.assertNotIn("-abs-mem", self.manager().build_command())

    def test_the_cap_can_be_disabled_or_overridden(self):
        snapshot = {"job_limit": {"limit_bytes": 28 * 1024 ** 3}}
        fake = mock.Mock(snapshot=lambda: snapshot)
        with mock.patch.object(WSCLEAN, "load_local_module",
                               return_value=fake):
            with mock.patch.dict(WSCLEAN.IMAGING_TASK_CONFIG,
                                 {"wsclean_memory_fraction": 0}):
                self.assertNotIn("-abs-mem", self.manager().build_command())
            command = self.manager(extra_args=["-mem", "50"]).build_command()
        self.assertNotIn("-abs-mem", command)

    def test_the_cap_does_not_decide_whether_an_image_is_reused(self):
        first = ["wsclean", "-j", "4", "-abs-mem", "22", "x.ms"]
        second = ["wsclean", "-j", "4", "-abs-mem", "60", "x.ms"]
        self.assertEqual(WSCLEAN.comparable_command(first),
                         WSCLEAN.comparable_command(second))
        self.assertNotEqual(WSCLEAN.comparable_command(first),
                            WSCLEAN.comparable_command(
                                ["wsclean", "-j", "8", "x.ms"]))


class RunTests(MeerkatFixture):

    def test_the_report_records_the_katbeam_correction(self):
        manager = self.manager()
        report = self.run_manager(manager)

        parameters, products = report["imaging_parameters"], report["products"]
        self.assertEqual(report["imaging_backend"], "wsclean")
        self.assertEqual(report["datacolumn"], "DATA")
        self.assertEqual(parameters["primary_beam_provider"], "katbeam")
        self.assertEqual(parameters["beam_mode"], "katbeam")
        self.assertEqual(parameters["mfs_beam_mode"], "band_average")
        self.assertEqual(parameters["pbmin"], 0.3)
        self.assertEqual(products["pbcor_method"], "katbeam_wsclean_mfs")

        status = products["pbcor_status"]
        self.assertEqual(status["katbeam_model"], MODEL)
        self.assertEqual(status["katbeam_version"], "0.1.2")
        self.assertIs(status["spectral_terms_corrected"], False)
        self.assertAlmostEqual(status["katbeam_frequency_mhz"], 1284.0)
        self.assertTrue(Path(status["pbcor_provenance"]).is_file())

    def test_the_published_images_are_the_mfs_pair(self):
        report = self.run_manager(self.manager())
        products = report["products"]
        self.assertTrue(products["restored_image"].endswith("-MFS-image.fits"))
        self.assertTrue(products["pbcor_image"].endswith("-MFS-image-pb.fits"))
        self.assertEqual(len(products["subbands"]), 8)
        corrected = [item for item in products["subbands"]
                     if item.get("pbcor_image")]
        self.assertGreater(len(corrected), 0)

    def test_the_published_flux_scale_is_right(self):
        report = self.run_manager(self.manager())
        with fits.open(report["products"]["pbcor_image"]) as handle:
            corrected = handle[0].data[0, 0]
        with fits.open(
            report["products"]["pbcor_status"]["katbeam_pb_image"]
        ) as handle:
            beam = handle[0].data[0, 0]
        inside = beam >= 0.3
        self.assertGreater(inside.sum(), 300)
        # A flat-spectrum field seen through the real (frequency-dependent)
        # beam comes back at its true flux, to the noise-free precision of
        # the fake beam, over the whole usable field.
        np.testing.assert_allclose(corrected[inside], FLUX, rtol=0.02)

    def test_every_published_fits_file_is_inventoried(self):
        report = self.run_manager(self.manager())
        products = report["products"]
        names = {Path(path).name for path in products["fits"]}
        prefix = "target_b_final_continuum"
        for expected in ("-MFS-image.fits", "-MFS-image-pb.fits",
                         "-MFS-katbeam-pb.fits", "-0000-image-pb.fits",
                         "-0004-katbeam-pb.fits"):
            self.assertIn(prefix + expected, names)
        valid, reason = PIPELINE.IMAGING_SAFETY.validate_file_content_identities(
            products["fits"], products["fits_content_identities"])
        self.assertTrue(valid, reason)

    def test_a_finished_run_is_reused(self):
        manager = self.manager()
        self.run_manager(manager)
        with mock.patch.object(
            WSCLEAN.subprocess, "run",
            side_effect=AssertionError("WSClean must not run again"),
        ), mock.patch("builtins.print"):
            again = self.manager().run()
        self.assertEqual(again["imaging_backend"], "wsclean")

    def test_reuse_survives_a_different_memory_allocation(self):
        snapshot = {"job_limit": {"limit_bytes": 28 * 1024 ** 3}}
        with mock.patch.object(WSCLEAN, "load_local_module",
                               return_value=mock.Mock(snapshot=lambda: snapshot)):
            self.run_manager(self.manager())
        bigger = {"job_limit": {"limit_bytes": 64 * 1024 ** 3}}
        with mock.patch.object(WSCLEAN, "load_local_module",
                               return_value=mock.Mock(snapshot=lambda: bigger)):
            self.assertIsNotNone(self.manager().reusable_report())

    def test_changed_beam_settings_do_not_reuse_the_old_correction(self):
        self.run_manager(self.manager())
        for key, value in (("meerkat_mfs_beam_mode", "reference_frequency"),
                           ("pbmin", 0.5)):
            with self.subTest(key=key):
                with mock.patch.dict(WSCLEAN.IMAGING_TASK_CONFIG, {key: value}):
                    self.assertIsNone(self.manager().reusable_report())

    def test_a_tampered_corrected_image_is_not_reused(self):
        report = self.run_manager(self.manager())
        pbcor = Path(report["products"]["pbcor_image"])
        with fits.open(pbcor, mode="update") as handle:
            self.assertTrue(np.isfinite(handle[0].data[0, 0, NY // 2, NX // 2]))
            handle[0].data[0, 0, NY // 2, NX // 2] += 5.0
        self.assertIsNone(self.manager().reusable_report())

    def test_a_missing_provenance_file_is_not_reused(self):
        report = self.run_manager(self.manager())
        Path(report["products"]["pbcor_status"]["pbcor_provenance"]).unlink()
        self.assertIsNone(self.manager().reusable_report())

    def test_a_failed_wsclean_run_says_so(self):
        manager = self.manager()
        with mock.patch.object(
            WSCLEAN.subprocess, "run",
            return_value=mock.Mock(returncode=1),
        ), mock.patch.object(WSCLEAN, "_MEERKAT_BEAM", BEAM), \
                mock.patch.object(
                    BEAM, "load_katbeam",
                    return_value=(FakeJimBeam, "0.1.2")), \
                mock.patch("builtins.print"):
            with self.assertRaisesRegex(RuntimeError, "status 1"):
                manager.run()

    def test_preflight_names_a_missing_wsclean_and_a_missing_katbeam(self):
        manager = self.manager(executable=str(self.root / "no_such_wsclean"))
        with self.assertRaisesRegex(FileNotFoundError, "pass --wsclean"):
            manager.preflight()
        manager = self.manager()
        with mock.patch.object(WSCLEAN, "_MEERKAT_BEAM", BEAM), \
                mock.patch.dict("sys.modules", {"katbeam": None}):
            with self.assertRaisesRegex(RuntimeError, "--imager casa"):
                manager.preflight()


class RetryTests(MeerkatFixture):
    """A correction that fails must not cost the WSClean run before it."""

    def failing_correction(self, manager):
        def fake_wsclean(command, check=False):
            write_wsclean_products(manager.prefix)
            return FakeSubprocessResult()

        def broken(*args, **kwargs):
            raise BEAM.BeamCorrectionError("the beam went wrong")

        with mock.patch.object(WSCLEAN.subprocess, "run", fake_wsclean), \
                mock.patch.object(WSCLEAN, "_MEERKAT_BEAM", BEAM), \
                mock.patch.object(BEAM, "load_katbeam",
                                  return_value=(FakeJimBeam, "0.1.2")), \
                mock.patch.object(BEAM, "correct_products", broken), \
                mock.patch("builtins.print"):
            with self.assertRaisesRegex(RuntimeError, "went wrong"):
                manager.run()

    def retry(self, manager):
        with mock.patch.object(
            WSCLEAN.subprocess, "run",
            side_effect=AssertionError("WSClean must not run again"),
        ), mock.patch.object(WSCLEAN, "_MEERKAT_BEAM", BEAM), \
                mock.patch.object(BEAM, "load_katbeam",
                                  return_value=(FakeJimBeam, "0.1.2")), \
                mock.patch("builtins.print"):
            return manager.run()

    def test_a_failed_correction_is_retried_without_wsclean(self):
        manager = self.manager()
        self.failing_correction(manager)
        self.assertTrue(manager.finished_marker_path().is_file())
        self.assertFalse(manager.report_path.exists())

        report = self.retry(self.manager())
        self.assertTrue(report["products"]["pbcor_image"].endswith(
            "-MFS-image-pb.fits"))

    def test_a_changed_ms_or_command_or_output_runs_wsclean_again(self):
        def reran(prepare):
            manager = self.manager()
            self.failing_correction(manager)
            prepare(manager)
            return self.manager().finished_wsclean_products() is None

        self.assertFalse(reran(lambda manager: None))
        self.assertTrue(reran(
            lambda manager: (self.final / "table.dat").write_bytes(b"DATA X")))
        (self.final / "table.dat").write_bytes(b"DATA FLAG WEIGHT")

        def edit_raw_output(manager):
            path = Path(f"{manager.prefix}-MFS-image.fits")
            with fits.open(path, mode="update") as handle:
                handle[0].data[0, 0, NY // 2, NX // 2] += 1.0
        self.assertTrue(reran(edit_raw_output))

    def test_a_different_command_is_a_different_run(self):
        manager = self.manager()
        self.failing_correction(manager)
        self.assertIsNotNone(self.manager().finished_wsclean_products())
        changed = self.manager(threads=9)
        self.assertIsNone(changed.finished_wsclean_products())

    def test_the_memory_cap_does_not_discard_a_finished_run(self):
        manager = self.manager()
        self.failing_correction(manager)
        snapshot = {"job_limit": {"limit_bytes": 28 * 1024 ** 3}}
        with mock.patch.object(WSCLEAN, "load_local_module",
                               return_value=mock.Mock(snapshot=lambda: snapshot)):
            self.assertIsNotNone(self.manager().finished_wsclean_products())

    def test_lofar_keeps_its_behaviour_and_writes_no_marker(self):
        write_json(self.outdir / "inspection" / "inspection_report.json", {
            "measurement_set": str(self.final),
            "telescope": {"name": "lofar"}, "band": {"name": "LOFAR_HBA"},
            "data_columns": ["DATA", "CORRECTED_DATA"],
        })
        for path in (self.outdir / "calibration" / "calibration_report.json",
                     self.outdir / "selfcal" / "reports" / "selfcal_report.json"):
            data = json.loads(path.read_text())
            data["telescope"], data["band"] = "lofar", "LOFAR_HBA"
            write_json(path, data)
        manager = WSCLEAN.WSCleanImagingManager(
            self.final, outdir=self.outdir, telescope="lofar",
            receiver_band="HBA", executable=self.executable, threads=2,
            beam_mode="full", datacolumn="CORRECTED_DATA",
            provenance_fingerprint="f" * 64, target=self.TARGET)

        def fake_wsclean(command, check=False):
            for suffix in ("-MFS-image.fits", "-MFS-image-pb.fits"):
                write_image(Path(f"{manager.prefix}{suffix}"),
                            np.ones((NY, NX)), wsclean_header(150e6, 30e6))
            return FakeSubprocessResult()

        with mock.patch.object(WSCLEAN.subprocess, "run", fake_wsclean), \
                mock.patch("builtins.print"):
            report = manager.run()
        self.assertEqual(report["imaging_parameters"]["primary_beam_provider"],
                         "wsclean_everybeam")
        self.assertFalse(manager.finished_marker_path().exists())


class BiasEstimateIsInformationalTests(MeerkatFixture):

    def test_a_failing_estimate_does_not_stop_the_correction(self):
        class AwkwardBeam(FakeJimBeam):
            def I(self, x, y, freq_mhz):
                if np.asarray(x).shape[0] == 1:
                    raise ValueError("only two-dimensional grids")
                return super().I(x, y, freq_mhz)

        manager = self.manager()

        def fake_wsclean(command, check=False):
            write_wsclean_products(manager.prefix)
            return FakeSubprocessResult()

        with mock.patch.object(WSCLEAN.subprocess, "run", fake_wsclean), \
                mock.patch.object(WSCLEAN, "_MEERKAT_BEAM", BEAM), \
                mock.patch.object(BEAM, "load_katbeam",
                                  return_value=(AwkwardBeam, "0.1.2")), \
                mock.patch("builtins.print") as printed:
            report = manager.run()

        bias = report["products"]["pbcor_status"]["single_frequency_flux_bias"]
        self.assertIn("error", bias)
        self.assertTrue(Path(report["products"]["pbcor_image"]).is_file())
        shown = " ".join(str(call) for call in printed.call_args_list)
        self.assertIn("correction itself is unaffected", shown)


class FlaggingMaskTests(MeerkatFixture):
    """Finding 165: the masks the flagging stage applied reach the beam."""

    MASKS_MHZ = [[933.0, 960.0], [1163.0, 1299.0], [1524.0, 1630.0]]

    def write_mask_report(self):
        write_json(
            self.outdir / "flagging" / "known_frequency_mask_report.json", {
                "scheme": "spw_intersected_static_frequency_masks",
                "applied_intersections": [
                    {"spw": "0", "intersection_MHz": list(item)}
                    for item in self.MASKS_MHZ
                ],
            })

    def forbid_wsclean(self, manager):
        with mock.patch.object(
            WSCLEAN.subprocess, "run",
            side_effect=AssertionError("WSClean must not run again"),
        ), mock.patch.object(WSCLEAN, "_MEERKAT_BEAM", BEAM), \
                mock.patch.object(BEAM, "load_katbeam",
                                  return_value=(FakeJimBeam, "0.1.2")), \
                mock.patch("builtins.print"):
            return manager.run()

    def test_without_a_report_no_mask_is_assumed(self):
        self.assertEqual(self.manager().known_frequency_masks_hz(), [])

    def test_the_masks_reach_the_correction(self):
        self.write_mask_report()
        report = self.run_manager(self.manager())
        status = report["products"]["pbcor_status"]
        self.assertEqual(status["masked_ranges_hz"],
                         [[low * 1e6, high * 1e6]
                          for low, high in self.MASKS_MHZ])
        provenance = json.loads(
            Path(status["pbcor_provenance"]).read_text())
        # 1177-1284 MHz lies inside the 1163-1299 MHz mask; 1498-1605 MHz
        # keeps only 1498-1524 MHz.
        self.assertIn("masked", provenance["subbands"][3]["skipped"])
        self.assertAlmostEqual(
            provenance["subbands"][6]["frequency_mhz"], 1511.0)
        self.assertEqual(provenance["mfs_weighting"], "covered_fraction")

    def test_an_old_correction_is_redone_without_imaging_again(self):
        self.run_manager(self.manager())
        manager = self.manager()
        old = json.loads(manager.report_path.read_text())
        old["products"]["pbcor_status"].pop("mfs_weighting")
        write_json(manager.report_path, old)
        self.assertIsNone(manager.reusable_report())
        again = self.forbid_wsclean(self.manager())
        self.assertEqual(
            again["products"]["pbcor_status"]["mfs_weighting"],
            "covered_fraction")

    def test_new_masks_redo_the_correction_without_imaging_again(self):
        self.run_manager(self.manager())
        self.write_mask_report()
        self.assertIsNone(self.manager().reusable_report())
        again = self.forbid_wsclean(self.manager())
        self.assertEqual(len(again["products"]["pbcor_status"][
            "masked_ranges_hz"]), 3)
        self.assertIsNotNone(self.manager().reusable_report())


class OrchestratorValidationTests(MeerkatFixture):
    """The report the stage writes is what the orchestrator checks."""

    def runner(self):
        argv = [
            str(self.raw), "--outdir", str(self.outdir), "--dry-run",
            "--mpi", "off", "--no-copy-input-ms", "--telescope", "meerkat",
            "--receiver-band", "L", "--data-state", "raw",
            "--target", self.TARGET, "--stages", "imaging",
        ]
        return PIPELINE.PipelineRunner(PIPELINE.build_parser().parse_args(argv))

    def finished(self):
        runner = self.runner()
        manager = self.manager(
            provenance_fingerprint=runner.stage_fingerprint("imaging"))
        self.run_manager(manager)
        return runner, manager

    def rewrite_report(self, manager, change):
        report = json.loads(manager.report_path.read_text())
        change(report)
        write_json(manager.report_path, report)
        return report

    def test_a_report_the_stage_wrote_is_valid(self):
        runner, manager = self.finished()
        self.assertEqual(runner.imaging_backend(), "wsclean")
        valid, reason = runner.validate_stage_completion("imaging")
        self.assertTrue(valid, reason)

    def test_the_orchestrator_expects_the_backend_of_this_run(self):
        runner, manager = self.finished()
        self.rewrite_report(
            manager, lambda r: r.update(imaging_backend="casa"))
        valid, reason = runner.validate_stage_completion("imaging")
        self.assertFalse(valid)
        self.assertIn("backend", reason)

    def test_provenance_mismatches_are_each_caught(self):
        cases = {
            "model": lambda r: r["products"]["pbcor_status"].update(
                katbeam_model="MKAT-AA-UHF-JIM-2020"),
            "method": lambda r: r["products"].update(pbcor_method="other"),
            "version": lambda r: r["products"]["pbcor_status"].update(
                katbeam_version="unknown"),
            "frequency": lambda r: r["products"]["pbcor_status"].update(
                katbeam_frequency_mhz=0),
            "spectral": lambda r: r["products"]["pbcor_status"].update(
                spectral_terms_corrected=True),
            "mode": lambda r: r["products"]["pbcor_status"].update(
                mfs_beam_mode="reference_frequency"),
            "no pbcor": lambda r: r["products"].update(pbcor_image=None),
            # A report written before finding 165 records no weighting.
            "weighting": lambda r: r["products"]["pbcor_status"].pop(
                "mfs_weighting"),
            "other weighting": lambda r: r["products"]["pbcor_status"].update(
                mfs_weighting="uniform_band_samples"),
        }
        for label, change in cases.items():
            with self.subTest(label):
                runner, manager = self.finished()
                self.rewrite_report(manager, change)
                valid, reason = runner.validate_stage_completion("imaging")
                self.assertFalse(valid, label)
                self.assertTrue(reason)

    def test_a_changed_corrected_image_no_longer_validates(self):
        runner, manager = self.finished()
        report = json.loads(manager.report_path.read_text())
        pbcor = Path(report["products"]["pbcor_image"])
        with fits.open(pbcor, mode="update") as handle:
            # The centre is inside the field; a blanked corner is NaN, and
            # NaN plus one leaves the file unchanged.
            self.assertTrue(np.isfinite(handle[0].data[0, 0, NY // 2, NX // 2]))
            handle[0].data[0, 0, NY // 2, NX // 2] += 1.0
        valid, reason = runner.validate_stage_completion("imaging")
        self.assertFalse(valid)
        self.assertRegex(reason, "content changed|no longer matches")

    def test_a_changed_beam_image_no_longer_validates(self):
        runner, manager = self.finished()
        report = json.loads(manager.report_path.read_text())
        beam = Path(report["products"]["pbcor_status"]["katbeam_pb_image"])
        with fits.open(beam, mode="update") as handle:
            handle[0].data[0, 0, 3, 3] *= 0.5
        valid, reason = runner.validate_stage_completion("imaging")
        self.assertFalse(valid)

    def test_a_provenance_file_that_names_other_images_is_rejected(self):
        runner, manager = self.finished()
        report = json.loads(manager.report_path.read_text())
        path = Path(report["products"]["pbcor_status"]["pbcor_provenance"])
        provenance = json.loads(path.read_text())
        provenance["source_image"] = str(self.root / "other-image.fits")
        write_json(path, provenance)
        valid, reason = runner.validate_stage_completion("imaging")
        self.assertFalse(valid)
        self.assertIn("source_image", reason)


class OrchestratorCommandContractTests(MeerkatFixture):
    """The command the orchestrator builds is one the stage can parse.

    The WSClean stage once failed to start for every telescope because the
    orchestrator passed --upstream-outdir, which the stage's parser did not
    know. Building the real command and parsing it with the real parser is
    what catches that.
    """

    def parse_stage_command(self, argv):
        args = PIPELINE.build_parser().parse_args(argv)
        runner = PIPELINE.PipelineRunner(args)
        command = runner.command_for_stage("imaging")
        self.assertEqual(command[1], str(Path(WSCLEAN.__file__)))
        return WSCLEAN.build_parser().parse_args(command[2:]), runner

    def test_meerkat_command_is_accepted_by_the_stage_parser(self):
        parsed, _ = self.parse_stage_command([
            str(self.raw), "--outdir", str(self.outdir), "--dry-run",
            "--no-copy-input-ms", "--telescope", "meerkat",
            "--receiver-band", "L", "--stages", "imaging",
            "--target", self.TARGET,
        ])
        self.assertEqual(Path(parsed.ms).resolve(), self.final.resolve())
        self.assertEqual(Path(parsed.upstream_outdir).resolve(),
                         self.outdir.resolve())
        self.assertEqual(parsed.datacolumn, "auto")
        self.assertIsNone(parsed.beam_mode)

    def test_lofar_command_is_accepted_by_the_stage_parser(self):
        parsed, _ = self.parse_stage_command([
            str(self.final), "--outdir", str(self.outdir), "--dry-run",
            "--no-copy-input-ms", "--telescope", "lofar",
            "--receiver-band", "HBA", "--data-state", "calibrated",
            "--stages", "imaging", "--imaging-ms", str(self.final),
            "--lofar-beam-mode", "full",
        ])
        self.assertEqual(parsed.beam_mode, "full")
        self.assertEqual(parsed.datacolumn, "CORRECTED_DATA")

    def test_a_per_target_scope_is_accepted_by_the_stage_parser(self):
        args = PIPELINE.build_parser().parse_args([
            str(self.raw), "--outdir", str(self.outdir), "--dry-run",
            "--no-copy-input-ms", "--telescope", "meerkat",
            "--receiver-band", "L", "--stages", "imaging",
        ])
        runner = PIPELINE.PipelineRunner(args)
        scopes = runner.science_target_scopes(["imaging"])
        self.assertEqual([target for target, _ in scopes], [self.TARGET])
        target, root = scopes[0]
        with runner.target_scope(target, root):
            command = runner.command_for_stage("imaging")
        parsed = WSCLEAN.build_parser().parse_args(command[2:])
        self.assertEqual(Path(parsed.outdir).resolve(), Path(root).resolve())
        self.assertEqual(Path(parsed.upstream_outdir).resolve(),
                         self.outdir.resolve())
        self.assertEqual(parsed.target, self.TARGET)


class ProfileDefaultTests(unittest.TestCase):

    def test_the_defaults_the_user_asked_for(self):
        profiles = PIPELINE.TELESCOPE_PROFILE_MODULE
        self.assertEqual(profiles.default_imager("meerkat"), "wsclean")
        self.assertEqual(profiles.default_imager("lofar"), "wsclean")
        self.assertEqual(profiles.default_imager("ugmrt"), "casa")


class RecipeStabilityTests(unittest.TestCase):
    """Settings that cannot affect a run must not change its recipe."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.ms = self.root / "in.ms"
        make_casa_table(self.ms)

    def tearDown(self):
        self.temporary.cleanup()

    def recipe(self, *extra):
        args = PIPELINE.build_parser().parse_args([
            str(self.ms), "--outdir", str(self.root / "P"), "--dry-run",
            "--no-copy-input-ms", "--stages", "imaging", *extra,
        ])
        runner = PIPELINE.PipelineRunner(args)
        return runner.computed_stage_recipe_fingerprint("imaging")

    def changed_by(self, settings, band_settings, *extra):
        before = self.recipe(*extra)
        config = PIPELINE.CONFIG
        task = dict(config.FINAL_IMAGING_TASK_CONFIG)
        task.update(settings)
        band = dict(config.FINAL_IMAGING_CONFIG["MEERKAT_L"])
        band.update(band_settings)
        bands = dict(config.FINAL_IMAGING_CONFIG, MEERKAT_L=band)
        with mock.patch.dict(PIPELINE.IMAGING_TASK_CONFIG, task), \
                mock.patch.object(config, "FINAL_IMAGING_CONFIG", bands):
            return self.recipe(*extra) != before

    def test_the_memory_cap_is_never_science(self):
        for extra in (
            ("--telescope", "meerkat", "--receiver-band", "L"),
            ("--telescope", "meerkat", "--receiver-band", "L",
             "--imager", "casa"),
            ("--telescope", "ugmrt", "--receiver-band", "B4"),
        ):
            with self.subTest(extra=extra):
                self.assertFalse(self.changed_by(
                    {"wsclean_memory_fraction": 0.3}, {}, *extra))

    def test_wsclean_only_settings_leave_casa_recipes_alone(self):
        for extra in (
            ("--telescope", "meerkat", "--receiver-band", "L",
             "--imager", "casa"),
            ("--telescope", "ugmrt", "--receiver-band", "B4"),
        ):
            with self.subTest(extra=extra):
                self.assertFalse(self.changed_by(
                    {"meerkat_mfs_beam_mode": "reference_frequency"},
                    {"wsclean_channels_out": 3}, *extra))

    def test_they_do_change_a_meerkat_wsclean_recipe(self):
        extra = ("--telescope", "meerkat", "--receiver-band", "L")
        self.assertTrue(self.changed_by(
            {"meerkat_mfs_beam_mode": "reference_frequency"}, {}, *extra))
        self.assertTrue(self.changed_by({}, {"wsclean_channels_out": 3},
                                        *extra))

    def test_the_flagging_masks_enter_only_a_meerkat_wsclean_recipe(self):
        # Finding 165: the katbeam correction reads them, so a MeerKAT
        # WSClean image is redone when they change, and nothing else is.
        report = self.root / "P" / "flagging" / "known_frequency_mask_report.json"
        for extra, expected in (
            (("--telescope", "meerkat", "--receiver-band", "L"), True),
            (("--telescope", "meerkat", "--receiver-band", "L",
              "--imager", "casa"), False),
            (("--telescope", "ugmrt", "--receiver-band", "B4"), False),
        ):
            with self.subTest(extra=extra):
                if report.exists():
                    report.unlink()
                before = self.recipe(*extra)
                write_json(report, {"applied_intersections": [
                    {"spw": "0", "intersection_MHz": [1163.0, 1299.0]}]})
                self.assertEqual(self.recipe(*extra) != before, expected)


if __name__ == "__main__":
    unittest.main()
