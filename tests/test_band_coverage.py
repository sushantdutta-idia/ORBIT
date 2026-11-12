"""Every uGMRT and MeerKAT band has what a run needs, sized for its beam.

One rule sizes every band's images (finding 167): four or five pixels across
the beam of the longest baselines at the band's centre, and a field that
reaches out to where the primary beam, at the band's lowest frequency, has
fallen to about a tenth. It is the rule that set Band 4's 5120 x 1" (beam
0.12 at the edge at 550 MHz), and the runs that were made meet it: uGMRT Band
5 at 6144 x 0.5" and MeerKAT L at 6144 x 1.5". It is checked here for the
per-target chain's grid and for the pipeline's own self-calibration and final
images, for every band.
"""

import importlib.util
import math
import unittest
from pathlib import Path

from test_resume_safety import PIPELINE


SOURCE_DIR = Path(__file__).resolve().parents[1]
CONFIG = PIPELINE.CONFIG
SAFETY = PIPELINE.IMAGING_SAFETY
PROFILES = PIPELINE.TELESCOPE_PROFILE_MODULE.TELESCOPE_PROFILES


def load_run_targets():
    spec = importlib.util.spec_from_file_location(
        "test_band_coverage_run_targets", SOURCE_DIR / "run_targets.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RUN_TARGETS = load_run_targets()

BANDS = (
    [("ugmrt", band) for band in ("B2", "B3", "B4", "B5")]
    + [("meerkat", band) for band in (
        "MEERKAT_UHF", "MEERKAT_L", "MEERKAT_S0", "MEERKAT_S1",
        "MEERKAT_S2", "MEERKAT_S3", "MEERKAT_S4")]
)
# The longest baselines: the GMRT's arms, MeerKAT's outer antennas.
LONGEST_BASELINE_M = {"ugmrt": 25000.0, "meerkat": 7700.0}
LIGHT = 299792458.0
ARCSEC_PER_RADIAN = 206264.80624709636
EDGE_BEAM = 0.15
PIXELS_PER_BEAM = (3.0, 8.0)


def science_range_hz(profile):
    """Where a band holds data worth imaging: science, beam or recorded."""

    valid = profile.get("primary_beam_valid_frequency_mhz")
    low = profile.get("science_freq_min") or (
        valid[0] * 1e6 if valid else profile["freq_min"])
    high = profile.get("science_freq_max") or (
        valid[1] * 1e6 if valid else profile["freq_max"])
    return float(low), float(high)


def arcsec(value):
    text = str(value).strip().lower()
    return float(text[:-len("arcsec")]) if text.endswith("arcsec") else float(text)


def meerkat_beam(radius_arcmin, frequency_ghz):
    """The cosine-tapered MeerKAT beam of Mauch et al. (2020).

    FWHM 57.5' at 1.5 GHz, scaling as 1/frequency: an approximation to
    katbeam, good enough to size a field, used here because katbeam itself
    is not installed where the tests run.
    """

    width = 57.5 * 1.5 / float(frequency_ghz)
    x = 1.189 * float(radius_arcmin) / width
    if abs(1.0 - 4.0 * x * x) < 1e-9:
        return (math.pi / 4.0) ** 2
    return (math.cos(math.pi * x) / (1.0 - 4.0 * x * x)) ** 2


def radius_to_edge_beam(telescope, profile, frequency_hz):
    """Arcmin out to where the beam at this frequency falls to EDGE_BEAM."""

    ghz = frequency_hz / 1e9
    if telescope == "ugmrt":
        return SAFETY.first_outward_beam_cutoff(
            profile["primary_beam_coefficients"], EDGE_BEAM) / ghz
    radius = 0.0
    while meerkat_beam(radius, ghz) > EDGE_BEAM:
        radius += 0.05
    return radius


def problems(telescope, band, imsize, cell):
    """What is wrong with a grid for a band, as sentences; none if right."""

    profile = PROFILES[telescope]["bands"][band]
    low, high = science_range_hz(profile)
    centre = 0.5 * (low + high)
    beam = LIGHT / centre / LONGEST_BASELINE_M[telescope] * ARCSEC_PER_RADIAN
    found = []
    pixels = beam / float(cell)
    if not PIXELS_PER_BEAM[0] <= pixels <= PIXELS_PER_BEAM[1]:
        found.append(f"{pixels:.1f} pixels to the {beam:.1f}\" beam")
    half = float(imsize) * float(cell) / 120.0
    needed = radius_to_edge_beam(telescope, profile, low)
    if half < needed:
        found.append(
            f"the field ends {half:.1f}' out, where the beam at "
            f"{low / 1e6:.0f} MHz has not yet fallen to {EDGE_BEAM:g} "
            f"({needed:.1f}')")
    return found


class SettingsTests(unittest.TestCase):

    def test_every_band_has_what_each_stage_reads(self):
        for telescope, band in BANDS:
            with self.subTest(band=band):
                for table in ("CAL_CONFIG", "SELFCAL_CONFIG",
                              "SELFCAL_IMAGING_CONFIG",
                              "FINAL_IMAGING_CONFIG", "DD_CONFIG"):
                    self.assertIn(band, getattr(CONFIG, table), table)
                self.assertIn(band, CONFIG.FLAGGING_CONFIG["edge_fraction"])
                self.assertIn((telescope, band), RUN_TARGETS.TARGET_GRIDS)

    def test_a_primary_beam_is_measured_wherever_one_is_applied(self):
        for telescope, band in BANDS:
            profile = PROFILES[telescope]["bands"][band]
            with self.subTest(band=band):
                if telescope == "ugmrt":
                    low, high = science_range_hz(profile)
                    valid = profile["primary_beam_valid_frequency_mhz"]
                    self.assertTrue(RUN_TARGETS.has_primary_beam(
                        telescope, profile))
                    self.assertLessEqual(valid[0] * 1e6, low)
                    self.assertGreaterEqual(valid[1] * 1e6, high)
                elif band in ("MEERKAT_UHF", "MEERKAT_L"):
                    self.assertTrue(profile.get("katbeam_model"))
                else:
                    # katbeam publishes no S-band model: the chain asks for
                    # --no-primary-beam rather than guess one.
                    self.assertFalse(RUN_TARGETS.has_primary_beam(
                        telescope, profile))


class GeometryTests(unittest.TestCase):

    def test_the_rule_is_met_by_the_runs_that_were_made(self):
        for telescope, band, imsize, cell in (
            ("ugmrt", "B4", 5120, 1.0),
            ("ugmrt", "B5", 6144, 0.5),
            ("meerkat", "MEERKAT_L", 6144, 1.5),
        ):
            with self.subTest(band=band):
                self.assertEqual(problems(telescope, band, imsize, cell), [])

    def test_the_grids_given_before_fail_it(self):
        # What the pipeline held for the bands not yet run: so the rule can
        # tell a field too small or a pixel too fine.
        for telescope, band, imsize, cell in (
            ("ugmrt", "B2", 8192, 1.0),
            ("ugmrt", "B3", 8192, 0.8),
            ("ugmrt", "B5", 8192, 0.3),
            ("meerkat", "MEERKAT_UHF", 6144, 2.0),
            ("meerkat", "MEERKAT_S0", 6144, 0.6),
        ):
            with self.subTest(band=band):
                self.assertNotEqual(problems(telescope, band, imsize, cell), [])

    def test_the_target_chain_sizes_every_band_by_it(self):
        for telescope, band in BANDS:
            imsize, cell = RUN_TARGETS.TARGET_GRIDS[(telescope, band)]
            with self.subTest(band=band):
                self.assertEqual(problems(telescope, band, imsize, cell), [])

    def test_the_pipelines_own_images_are_sized_by_it(self):
        for telescope, band in BANDS:
            for table in ("SELFCAL_IMAGING_CONFIG", "FINAL_IMAGING_CONFIG"):
                settings = getattr(CONFIG, table)[band]
                with self.subTest(band=band, table=table):
                    self.assertEqual(problems(
                        telescope, band, settings["imsize"],
                        arcsec(settings["cell"])), [])

    def test_the_chain_and_the_pipeline_agree_on_each_band(self):
        for telescope, band in BANDS:
            imsize, cell = RUN_TARGETS.TARGET_GRIDS[(telescope, band)]
            final = CONFIG.FINAL_IMAGING_CONFIG[band]
            with self.subTest(band=band):
                self.assertEqual(arcsec(final["cell"]), cell)
                self.assertEqual(int(final["imsize"]), imsize)


if __name__ == "__main__":
    unittest.main()
