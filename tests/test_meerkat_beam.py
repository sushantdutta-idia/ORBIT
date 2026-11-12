"""MeerKAT primary-beam correction of WSClean FITS images."""

import hashlib
import importlib.util
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from astropy.io import fits


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_radio_continuum_meerkat_beam", SOURCE_DIR / "meerkat_beam.py"
)
BEAM = importlib.util.module_from_spec(spec)
spec.loader.exec_module(BEAM)

NX, NY = 64, 48
CDELT = 0.02  # degrees per pixel
CRPIX_X, CRPIX_Y = 33.0, 25.0


class FakeJimBeam:
    """A circular beam whose width scales as 1/frequency, like a real one."""

    def __init__(self, model):
        self.model = model

    def I(self, x, y, freq_mhz):
        sigma = (1000.0 / float(freq_mhz)) / 2.3548
        return np.exp(-0.5 * ((np.asarray(x) / sigma) ** 2 +
                              (np.asarray(y) / sigma) ** 2))


class AsymmetricBeam:
    """Distinguishes x from y, so a transposed grid cannot pass."""

    def __init__(self, model):
        pass

    def I(self, x, y, freq_mhz):
        return np.exp(-0.5 * (((np.asarray(x) - 0.1) / 0.5) ** 2 +
                              (np.asarray(y) / 0.9) ** 2)) / np.exp(
            -0.5 * (0.1 / 0.5) ** 2)


def wsclean_header(frequency_hz, bandwidth_hz=None, **overrides):
    """A header laid out the way WSClean writes it (RA, Dec, freq, Stokes)."""

    header = fits.Header()
    header["SIMPLE"] = True
    header["BITPIX"] = -32
    header["NAXIS"] = 4
    header["NAXIS1"], header["NAXIS2"], header["NAXIS3"], header["NAXIS4"] = (
        NX, NY, 1, 1
    )
    header["CTYPE1"], header["CRPIX1"], header["CRVAL1"] = (
        "RA---SIN", CRPIX_X, 40.0)
    header["CDELT1"], header["CUNIT1"] = -CDELT, "deg"
    header["CTYPE2"], header["CRPIX2"], header["CRVAL2"] = (
        "DEC--SIN", CRPIX_Y, -30.0)
    header["CDELT2"], header["CUNIT2"] = CDELT, "deg"
    header["CTYPE3"], header["CRPIX3"] = "FREQ", 1.0
    header["CRVAL3"], header["CDELT3"] = frequency_hz, bandwidth_hz or 0.0
    header["CTYPE4"], header["CRPIX4"] = "STOKES", 1.0
    header["CRVAL4"], header["CDELT4"] = 1.0, 1.0
    header["BUNIT"] = "JY/BEAM"
    header["BMAJ"], header["BMIN"], header["BPA"] = 0.001, 0.001, 0.0
    for key, value in overrides.items():
        header[key] = value
    return header


def write_image(path, plane, header):
    fits.PrimaryHDU(
        data=np.asarray(plane, dtype=np.float32).reshape(1, 1, NY, NX),
        header=header,
    ).writeto(str(path), overwrite=True)


def geometry(**overrides):
    return BEAM.image_geometry(wsclean_header(1.28e9, 8e8, **overrides))


def identity(path):
    data = Path(path).read_bytes()
    return {"path": str(Path(path).resolve()), "size_bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


class GeometryTests(unittest.TestCase):

    def test_a_wsclean_header_is_read(self):
        value = geometry()
        self.assertEqual((value["nx"], value["ny"]), (NX, NY))
        self.assertEqual(value["crpix_x"], CRPIX_X)
        self.assertAlmostEqual(value["cdelt_x_deg"], -CDELT)
        self.assertAlmostEqual(value["frequency_hz"], 1.28e9)
        self.assertAlmostEqual(value["bandwidth_hz"], 8e8)

    def test_units_are_converted(self):
        value = BEAM.image_geometry(wsclean_header(
            1280.0, 800.0, CUNIT3="MHz", CDELT1=-CDELT * 3600.0,
            CUNIT1="arcsec"))
        self.assertAlmostEqual(value["frequency_hz"], 1.28e9)
        self.assertAlmostEqual(value["cdelt_x_deg"], -CDELT)

    def test_unsupported_layouts_are_refused(self):
        cases = {
            "not RA/Dec": dict(CTYPE1="GLON-SIN"),
            "projection": dict(CTYPE1="RA---CAR", CTYPE2="DEC--CAR"),
            "several planes": dict(NAXIS3=4),
        }
        for label, overrides in cases.items():
            with self.subTest(label):
                with self.assertRaises(BEAM.BeamCorrectionError):
                    BEAM.image_geometry(wsclean_header(1.28e9, **overrides))
        header = wsclean_header(1.28e9)
        header["CTYPE3"] = "VELO"
        with self.assertRaisesRegex(BEAM.BeamCorrectionError, "frequency"):
            BEAM.image_geometry(header)


class BeamPlaneTests(unittest.TestCase):

    def test_pixels_get_the_beam_at_their_own_offsets(self):
        # The orientation the CASA backend uses: offsets from the reference
        # pixel with the image's own increments, signs included. An
        # asymmetric beam makes a swapped or mirrored grid visible.
        value = geometry()
        plane = BEAM.beam_plane(AsymmetricBeam("m"), value, 1280.0)
        for iy, ix in [(0, 0), (10, 50), (24, 32), (47, 63), (5, 3)]:
            x = (ix - (CRPIX_X - 1.0)) * -CDELT
            y = (iy - (CRPIX_Y - 1.0)) * CDELT
            self.assertAlmostEqual(
                float(plane[iy, ix]),
                float(AsymmetricBeam("m").I(x, y, 1280.0)),
                places=5,
            )

    def test_chunking_does_not_change_the_result(self):
        value = geometry()
        whole = BEAM.beam_plane(FakeJimBeam("m"), value, 1100.0, 256)
        pieces = BEAM.beam_plane(FakeJimBeam("m"), value, 1100.0, 7)
        np.testing.assert_array_equal(whole, pieces)

    def test_a_beam_with_nan_is_refused(self):
        class Broken(FakeJimBeam):
            def I(self, x, y, freq_mhz):
                out = super().I(x, y, freq_mhz)
                out[0, 0] = np.nan
                return out
        with self.assertRaisesRegex(BEAM.BeamCorrectionError, "non-finite"):
            BEAM.beam_plane(Broken("m"), geometry(), 1280.0)

    def test_a_wrong_shape_is_refused(self):
        class Wrong(FakeJimBeam):
            def I(self, x, y, freq_mhz):
                return np.ones(3)
        with self.assertRaisesRegex(BEAM.BeamCorrectionError, "invalid"):
            BEAM.beam_plane(Wrong("m"), geometry(), 1280.0)

    def test_a_beam_that_is_not_unity_at_the_centre_is_refused(self):
        with self.assertRaisesRegex(BEAM.BeamCorrectionError, "unity"):
            BEAM.check_peak(np.full((4, 4), 0.5, dtype=np.float32))
        self.assertAlmostEqual(
            BEAM.check_peak(np.array([[1.0, 0.2]], dtype=np.float32)), 1.0)


class FrequencyTests(unittest.TestCase):

    def test_a_subband_inside_the_range_is_whole(self):
        centre, fraction = BEAM.effective_band(
            1.2e9, 1e8, (9e8, 1.67e9))
        self.assertEqual((centre, fraction), (1.2e9, 1.0))

    def test_an_edge_subband_uses_only_its_unflagged_part(self):
        # 856-963 MHz with everything below 900 MHz flagged.
        centre, fraction = BEAM.effective_band(
            909.5e6, 107e6, (900e6, 1.67e9))
        self.assertAlmostEqual(fraction, (963e6 - 900e6) / 107e6)
        self.assertAlmostEqual(centre, 0.5 * (900e6 + 963e6))

    def test_a_subband_outside_the_range_has_no_weight(self):
        self.assertEqual(
            BEAM.effective_band(800e6, 50e6, (900e6, 1.67e9))[1], 0.0)

    def test_without_a_bandwidth_the_nominal_frequency_is_used(self):
        self.assertEqual(BEAM.effective_band(1e9, None, None), (1e9, 1.0))
        self.assertEqual(BEAM.band_samples_mhz(1e9, None, None), [1000.0])

    def test_band_samples_stay_inside_the_science_range(self):
        samples = BEAM.band_samples_mhz(1.284e9, 0.856e9, (900e6, 1.67e9))
        self.assertAlmostEqual(min(samples), 900.0)
        self.assertAlmostEqual(max(samples), 1670.0)
        self.assertEqual(len(samples), BEAM.SAMPLES_PER_BAND)


class CorrectionTests(unittest.TestCase):
    """A uniform flat-spectrum sky seen through a frequency-dependent beam."""

    SUBBANDS = [1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0]
    WIDTH = 100.0
    FLUX = 2.5

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.dir = Path(self.temporary.name)
        self.prefix = "TARGET_final_continuum"
        self.beam = FakeJimBeam("m")
        self.value = geometry()
        self.sub_paths = []
        weights = []
        beams = []
        for index, freq in enumerate(self.SUBBANDS):
            path = self.dir / f"{self.prefix}-{index:04d}-image.fits"
            plane = self.FLUX * BEAM.beam_plane(self.beam, self.value, freq)
            write_image(path, plane, wsclean_header(
                freq * 1e6, self.WIDTH * 1e6))
            self.sub_paths.append(path)
            beams.append(BEAM.beam_plane(self.beam, self.value, freq))
            weights.append(1.0)
        self.mean_beam = np.mean(beams, axis=0)
        self.mfs = self.dir / f"{self.prefix}-MFS-image.fits"
        write_image(self.mfs, self.FLUX * self.mean_beam, wsclean_header(
            1250e6, 600e6))

    def tearDown(self):
        self.temporary.cleanup()

    def correct(self, **kwargs):
        defaults = dict(model="MKAT-AA-L-JIM-2020", pbmin=0.3,
                        science_limits_hz=(900e6, 1.67e9), identity=identity)
        defaults.update(kwargs)
        with mock.patch.object(
            BEAM, "load_katbeam", return_value=(FakeJimBeam, "0.1.2")
        ):
            return BEAM.correct_products(self.mfs, self.sub_paths, **defaults)

    @staticmethod
    def read(path):
        with fits.open(str(path)) as handle:
            return handle[0].data[0, 0].copy(), handle[0].header.copy()

    def test_the_band_averaged_correction_recovers_a_flat_spectrum_flux(self):
        record = self.correct()
        corrected, _ = self.read(record["pbcor_image"])
        beam, _ = self.read(record["pb_image"])
        inside = beam >= 0.3
        self.assertGreater(inside.sum(), 500)
        np.testing.assert_allclose(corrected[inside], self.FLUX, rtol=2e-4)
        self.assertEqual(record["mfs_beam_mode"], "band_average")
        self.assertIs(record["spectral_terms_corrected"], False)

    def test_each_subband_is_corrected_at_its_own_frequency(self):
        record = self.correct()
        self.assertEqual(len(record["subbands"]), len(self.SUBBANDS))
        for item, freq in zip(record["subbands"], self.SUBBANDS):
            self.assertAlmostEqual(item["frequency_mhz"], freq)
            corrected, _ = self.read(item["pbcor_image"])
            beam, _ = self.read(item["pb_image"])
            inside = beam >= 0.3
            np.testing.assert_allclose(
                corrected[inside], self.FLUX, rtol=2e-4)

    def test_a_single_frequency_correction_carries_the_predicted_bias(self):
        record = self.correct(mode="reference_frequency", pbmin=0.25)
        corrected, _ = self.read(record["pbcor_image"])
        beam, _ = self.read(record["pb_image"])
        inside = beam >= 0.25

        # Dividing by the 1250 MHz beam leaves, at every pixel, exactly
        # mean(B_k) / B_ref - 1 of the true flux: the physics the estimator
        # summarises. It is a few percent by the edge of the field.
        reference = BEAM.beam_plane(self.beam, self.value, 1250.0)
        expected = self.mean_beam / reference - 1.0
        measured = corrected / self.FLUX - 1.0
        np.testing.assert_allclose(
            measured[inside], expected[inside], atol=5e-4)
        self.assertGreater(float(np.max(measured[inside])), 0.02)
        self.assertEqual(record["mfs_beam_mode"], "reference_frequency")

        # And the number the report carries is that same quantity.
        bias = record["single_frequency_flux_bias"]["0.5"]
        x = bias["radius_deg"]
        analytic = np.mean([
            float(self.beam.I(x, 0.0, f)) for f in self.SUBBANDS
        ]) / float(self.beam.I(x, 0.0, 1250.0)) - 1.0
        self.assertAlmostEqual(
            bias["flux_bias_of_single_frequency_correction"], analytic,
            places=6)

    def test_pixels_below_pbmin_are_blank_and_unit_fields_are_right(self):
        record = self.correct(pbmin=0.4)
        corrected, header = self.read(record["pbcor_image"])
        beam, pb_header = self.read(record["pb_image"])
        self.assertTrue(np.all(np.isnan(corrected[beam < 0.4])))
        self.assertTrue(np.all(np.isfinite(corrected[beam >= 0.4])))
        self.assertEqual(header["BUNIT"], "JY/BEAM")
        # A beam is dimensionless and has no restoring beam.
        for key in ("BUNIT", "BMAJ", "BMIN", "BPA"):
            self.assertNotIn(key, pb_header)
        self.assertEqual(corrected.shape, (NY, NX))

    def test_products_follow_wsclean_naming(self):
        record = self.correct()
        self.assertEqual(
            Path(record["pbcor_image"]).name,
            f"{self.prefix}-MFS-image-pb.fits")
        self.assertEqual(
            Path(record["pb_image"]).name,
            f"{self.prefix}-MFS-katbeam-pb.fits")
        self.assertEqual(
            Path(record["subbands"][2]["pbcor_image"]).name,
            f"{self.prefix}-0002-image-pb.fits")

    def test_an_edge_subband_is_weighted_by_what_survived_flagging(self):
        record = self.correct(science_limits_hz=(1.05e9, 1.45e9))
        by_freq = {round(item["nominal_frequency_mhz"]): item
                   for item in record["subbands"]}
        # 1000 MHz spans 950-1050: only the top 0 MHz... lies outside, i.e.
        # none of it is inside 1050+, so it is skipped; 1100 spans 1050-1150.
        self.assertIn("skipped", by_freq[1000])
        self.assertAlmostEqual(by_freq[1100]["covered_fraction"], 1.0)
        self.assertIn("skipped", by_freq[1500])
        self.assertEqual(len(record["effective_frequencies_mhz"]), 4)

    def test_partial_edge_changes_the_effective_frequency(self):
        record = self.correct(science_limits_hz=(1.02e9, 1.67e9))
        first = record["subbands"][0]
        self.assertAlmostEqual(first["covered_fraction"], 0.30, places=2)
        self.assertAlmostEqual(first["frequency_mhz"], 0.5 * (1020 + 1050))

    def test_without_subbands_the_mfs_beam_is_averaged_over_the_band(self):
        with mock.patch.object(
            BEAM, "load_katbeam", return_value=(FakeJimBeam, "0.1.2")
        ):
            record = BEAM.correct_products(
                self.mfs, [], model="m", pbmin=0.3,
                science_limits_hz=(1.0e9, 1.5e9), identity=identity)
        self.assertEqual(record["subbands"], [])
        self.assertEqual(len(record["effective_frequencies_mhz"]),
                         BEAM.SAMPLES_PER_BAND)
        self.assertAlmostEqual(min(record["effective_frequencies_mhz"]),
                               1000.0)

    def test_identities_describe_the_files_on_disk(self):
        record = self.correct()
        for key in ("source", "pb", "pbcor"):
            path = record[f"{key}_image"]
            self.assertEqual(
                record[key]["sha256"],
                hashlib.sha256(Path(path).read_bytes()).hexdigest())
        item = record["subbands"][0]
        self.assertEqual(
            item["pbcor"]["sha256"],
            hashlib.sha256(Path(item["pbcor_image"]).read_bytes()).hexdigest())

    def test_provenance_round_trips_and_leaves_no_partial_files(self):
        record = self.correct()
        path = BEAM.write_provenance(self.dir / "p.provenance.json", record)
        self.assertEqual(json.loads(path.read_text())["model"], record["model"])
        self.assertEqual(
            [p.name for p in self.dir.iterdir() if ".partial." in p.name], [])

    def test_the_bias_estimate_matches_a_direct_calculation(self):
        freqs = [1000.0, 1250.0, 1500.0]
        weights = [1.0, 1.0, 1.0]
        result = BEAM.flat_spectrum_bias(
            FakeJimBeam("m"), 1250.0, freqs, weights, levels=[0.5])
        item = result["0.5"]
        r = item["radius_deg"]
        beams = [float(FakeJimBeam("m").I(r, 0.0, f)) for f in freqs]
        self.assertAlmostEqual(
            item["flux_bias_of_single_frequency_correction"],
            np.mean(beams) / beams[1] - 1.0, places=6)

    def test_an_unknown_mode_is_refused(self):
        with self.assertRaisesRegex(BEAM.BeamCorrectionError, "mode"):
            self.correct(mode="average")


class KatbeamAvailabilityTests(unittest.TestCase):

    def test_a_missing_package_says_how_to_proceed(self):
        with mock.patch.dict(sys.modules, {"katbeam": None}):
            with self.assertRaisesRegex(
                BEAM.BeamCorrectionError, "--imager casa"
            ):
                BEAM.load_katbeam()

    def test_an_unversioned_package_is_refused(self):
        fake = types.SimpleNamespace(JimBeam=FakeJimBeam)
        with mock.patch.dict(sys.modules, {"katbeam": fake}), \
                mock.patch("importlib.metadata.version",
                           side_effect=Exception("none")):
            with self.assertRaisesRegex(
                BEAM.BeamCorrectionError, "no discoverable version"
            ):
                BEAM.load_katbeam()

    def test_the_version_is_recorded(self):
        fake = types.SimpleNamespace(JimBeam=FakeJimBeam, __version__="0.1.2")
        with mock.patch.dict(sys.modules, {"katbeam": fake}):
            cls, version = BEAM.load_katbeam()
        self.assertIs(cls, FakeJimBeam)
        self.assertEqual(version, "0.1.2")


class ClassifyOutputsTests(unittest.TestCase):

    def paths(self, prefix, *suffixes):
        return [Path("/o") / f"{prefix}{suffix}" for suffix in suffixes]

    def test_wsclean_outputs_are_sorted_into_mfs_and_subbands(self):
        prefix = "G0231-0533_final_continuum"
        paths = self.paths(
            prefix, "-MFS-image.fits", "-0000-image.fits", "-0001-image.fits",
            "-MFS-residual.fits", "-0000-model.fits", "-MFS-image-pb.fits",
            "-MFS-katbeam-pb.fits", "-MFS-dirty.fits", "-psf.fits")
        mfs, subbands = BEAM.classify_outputs(paths, prefix)
        self.assertEqual(mfs.name, f"{prefix}-MFS-image.fits")
        self.assertEqual([p.name for p in subbands],
                         [f"{prefix}-0000-image.fits",
                          f"{prefix}-0001-image.fits"])

    def test_a_target_name_cannot_pose_as_a_product_type(self):
        prefix = "Model_beam_0001_final_continuum"
        paths = self.paths(prefix, "-MFS-image.fits", "-0000-image.fits",
                           "-MFS-residual.fits")
        mfs, subbands = BEAM.classify_outputs(paths, prefix)
        self.assertEqual(mfs.name, f"{prefix}-MFS-image.fits")
        self.assertEqual(len(subbands), 1)

    def test_without_joined_channels_the_plain_image_is_used(self):
        prefix = "T_final_continuum"
        mfs, subbands = BEAM.classify_outputs(
            self.paths(prefix, "-image.fits", "-residual.fits"), prefix)
        self.assertEqual(mfs.name, f"{prefix}-image.fits")
        self.assertEqual(subbands, [])


class MaskedFrequencyTests(unittest.TestCase):
    """Finding 165: a sub-band sits where its unflagged data lie."""

    def test_a_mask_moves_a_subband_to_its_unflagged_part(self):
        # 1350-1450 MHz less 1420-1470 MHz leaves 1350-1420 MHz.
        centre, fraction = BEAM.effective_band(
            1.4e9, 1e8, (9e8, 1.67e9), [(1.42e9, 1.47e9)])
        self.assertAlmostEqual(centre, 1.385e9)
        self.assertAlmostEqual(fraction, 0.7)

    def test_masks_inside_a_subband_weigh_each_surviving_piece(self):
        # 900-1100 MHz less 950-1000 and 1050-1200 leaves 900-950 and
        # 1000-1050: centre 975 MHz, half the width.
        centre, fraction = BEAM.effective_band(
            1.0e9, 2e8, None, [(1.05e9, 1.2e9), (0.95e9, 1.0e9)])
        self.assertAlmostEqual(centre, 0.975e9)
        self.assertAlmostEqual(fraction, 0.5)

    def test_the_meerkat_l_masks_empty_one_of_eight_subbands(self):
        # 1708779077: 856-1712 MHz in eight sub-bands, the 900-1670 MHz
        # science range and the three L-band masks.
        edges = np.linspace(856e6, 1712e6, 9)
        masks = [(933e6, 960e6), (1163e6, 1299e6), (1524e6, 1630e6)]
        bands = [BEAM.effective_band(0.5 * (low + high), high - low,
                                     (900e6, 1670e6), masks)
                 for low, high in zip(edges[:-1], edges[1:])]
        self.assertEqual(bands[3][1], 0.0)
        # 1498-1605 MHz keeps only 1498-1524 MHz.
        self.assertAlmostEqual(bands[6][0], 1511e6)
        self.assertAlmostEqual(bands[6][1], 26e6 / 107e6)
        # 1605-1712 MHz keeps 1630-1670 MHz, not 1605-1670 MHz.
        self.assertAlmostEqual(bands[7][0], 1650e6)

    def test_without_masks_nothing_changes(self):
        self.assertEqual(
            BEAM.effective_band(1.4e9, 1e8, (9e8, 1.67e9), []),
            BEAM.effective_band(1.4e9, 1e8, (9e8, 1.67e9)))

    def test_band_samples_avoid_the_masks(self):
        samples = BEAM.band_samples_mhz(
            1.25e9, 5e8, None, count=10, masked_hz=[(1.1e9, 1.3e9)])
        self.assertEqual(len(samples), 10)
        self.assertTrue(all(not 1100.0 < value < 1300.0 for value in samples))
        self.assertTrue(all(1000.0 <= value <= 1500.0 for value in samples))

    def test_the_wsclean_image_weight_is_read_and_checked(self):
        self.assertEqual(BEAM.wsclean_image_weight({"WSCIMGWG": 3.5e8}), 3.5e8)
        self.assertEqual(BEAM.wsclean_image_weight({"WSCIMGWG": 0.0}), 0.0)
        for bad in ({}, {"WSCIMGWG": -1.0}, {"WSCIMGWG": float("nan")},
                    {"WSCIMGWG": "x"}, {"WSCIMGWG": True}):
            with self.subTest(header=bad):
                self.assertIsNone(BEAM.wsclean_image_weight(bad))


class ImageWeightedCorrectionTests(unittest.TestCase):
    """The MFS beam WSClean's own weights imply (finding 165).

    Built like 1708779077: the imaging weights rise with frequency, one
    sub-band is wholly masked (no weight, no finite pixels) and one is cut
    by a mask. Each sub-band image is a flat-spectrum sky seen through the
    beam where its data lie, and the MFS image is WSClean's weighted mean.
    """

    SUBBANDS = [1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0]
    WEIGHTS = [1.0, 2.0, 0.0, 4.0, 6.0, 9.0]
    WIDTH = 100.0
    FLUX = 2.0
    MASKS = [(1.15e9, 1.25e9), (1.42e9, 1.47e9)]

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.dir = Path(self.temporary.name)
        self.prefix = "TARGET_final_continuum"
        self.beam = FakeJimBeam("m")
        self.value = geometry()
        self.sub_paths = []
        weighted = np.zeros((NY, NX))
        self.data_frequencies = []
        for index, (freq, weight) in enumerate(zip(self.SUBBANDS,
                                                    self.WEIGHTS)):
            path = self.dir / f"{self.prefix}-{index:04d}-image.fits"
            centre, fraction = BEAM.effective_band(
                freq * 1e6, self.WIDTH * 1e6, None, self.MASKS)
            if weight > 0:
                plane = self.FLUX * BEAM.beam_plane(
                    self.beam, self.value, centre / 1e6)
                weighted += weight * plane
                self.data_frequencies.append(centre / 1e6)
            else:
                plane = np.full((NY, NX), np.nan)
            write_image(path, plane, wsclean_header(
                freq * 1e6, self.WIDTH * 1e6, WSCIMGWG=weight))
            self.sub_paths.append(path)
        self.mfs = self.dir / f"{self.prefix}-MFS-image.fits"
        write_image(self.mfs, weighted / sum(self.WEIGHTS), wsclean_header(
            1250e6, 600e6))

    def tearDown(self):
        self.temporary.cleanup()

    def correct(self, paths=None, **kwargs):
        defaults = dict(model="MKAT-AA-L-JIM-2020", pbmin=0.3,
                        science_limits_hz=None, identity=identity,
                        masked_ranges_hz=self.MASKS)
        defaults.update(kwargs)
        with mock.patch.object(
            BEAM, "load_katbeam", return_value=(FakeJimBeam, "0.1.2")
        ):
            return BEAM.correct_products(
                self.mfs, paths or self.sub_paths, **defaults)

    @staticmethod
    def read(path):
        with fits.open(str(path)) as handle:
            return handle[0].data[0, 0].copy()

    def test_the_mfs_image_comes_back_at_its_true_flux(self):
        record = self.correct()
        corrected = self.read(record["pbcor_image"])
        beam = self.read(record["pb_image"])
        inside = beam >= 0.3
        self.assertGreater(inside.sum(), 500)
        np.testing.assert_allclose(corrected[inside], self.FLUX, rtol=2e-4)
        self.assertEqual(record["mfs_weighting"], "wsclean_image_weight")
        self.assertEqual(record["schema_version"], 2)

        # The equal-coverage beam used before misses by more than 3% at the
        # edge of the field, so the test can tell the two apart.
        coverage = []
        stack = []
        for freq in self.SUBBANDS:
            centre, fraction = BEAM.effective_band(
                freq * 1e6, self.WIDTH * 1e6, None, self.MASKS)
            if fraction > 0:
                coverage.append(fraction)
                stack.append(BEAM.beam_plane(self.beam, self.value,
                                             centre / 1e6))
        old = np.tensordot(coverage, stack, axes=1) / sum(coverage)
        mfs = self.read(self.mfs)
        self.assertGreater(
            float(np.max(np.abs(mfs[inside] / old[inside] / self.FLUX - 1))),
            0.03)

    def test_the_mfs_frequency_is_the_weighted_one(self):
        record = self.correct()
        weights = [weight for weight in self.WEIGHTS if weight > 0]
        expected = sum(
            w * f for w, f in zip(weights, self.data_frequencies)
        ) / sum(weights)
        self.assertAlmostEqual(record["mfs_effective_frequency_mhz"], expected)
        self.assertEqual(record["frequency_weights"], weights)
        self.assertEqual(record["masked_ranges_hz"],
                         [list(mask) for mask in sorted(self.MASKS)])

    def test_a_cut_subband_is_corrected_where_its_data_are(self):
        record = self.correct()
        item = record["subbands"][4]
        self.assertAlmostEqual(item["frequency_mhz"], 1385.0)
        self.assertAlmostEqual(item["covered_fraction"], 0.7)
        self.assertEqual(item["image_weight"], 6.0)
        corrected = self.read(item["pbcor_image"])
        beam = self.read(item["pb_image"])
        inside = beam >= 0.3
        np.testing.assert_allclose(corrected[inside], self.FLUX, rtol=2e-4)

    def test_a_masked_subband_is_skipped_and_its_old_products_go(self):
        stale = [self.dir / f"{self.prefix}-0002-katbeam-pb.fits",
                 self.dir / f"{self.prefix}-0002-image-pb.fits"]
        for path in stale:
            write_image(path, np.ones((NY, NX)), wsclean_header(1.2e9, 1e8))
        record = self.correct()
        item = record["subbands"][2]
        self.assertIn("masked", item["skipped"])
        self.assertNotIn("pbcor_image", item)
        self.assertEqual(sorted(item["removed_stale_products"]),
                         sorted(str(path.resolve()) for path in stale))
        for path in stale:
            self.assertFalse(path.exists())
        self.assertEqual(len(record["effective_frequencies_mhz"]), 5)

    def test_a_subband_without_weight_is_skipped_even_unmasked(self):
        record = self.correct(masked_ranges_hz=[(1.42e9, 1.47e9)])
        self.assertIn("no imaging weight", record["subbands"][2]["skipped"])
        corrected = self.read(record["pbcor_image"])
        beam = self.read(record["pb_image"])
        inside = beam >= 0.3
        np.testing.assert_allclose(corrected[inside], self.FLUX, rtol=2e-4)

    def test_one_missing_weight_falls_back_to_coverage_for_all(self):
        path = self.sub_paths[0]
        with fits.open(str(path), mode="update") as handle:
            del handle[0].header["WSCIMGWG"]
        record = self.correct()
        self.assertEqual(record["mfs_weighting"], "covered_fraction")
        self.assertIsNone(record["subbands"][0]["image_weight"])
        self.assertEqual(record["subbands"][1]["mfs_weight"], 1.0)

    def test_an_all_blank_image_without_a_weight_is_skipped(self):
        for path in self.sub_paths:
            with fits.open(str(path), mode="update") as handle:
                del handle[0].header["WSCIMGWG"]
        record = self.correct(masked_ranges_hz=[])
        self.assertIn("no finite pixels", record["subbands"][2]["skipped"])


if __name__ == "__main__":
    unittest.main()
