import cmath
import importlib.util
import json
import math
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_tool():
    spec = importlib.util.spec_from_file_location(
        "test_experiment_phase_offsets_module",
        SOURCE_DIR / "experiment_phase_offsets.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TOOL = load_tool()
NAMES = ["C00", "C01", "C02", "E02", "W03"]


def columns(phases_deg, flagged=()):
    """Gain-table columns for one G solution per antenna."""

    rows = len(phases_deg)
    values = np.ones((2, 1, rows), dtype=complex)
    flags = np.zeros((2, 1, rows), dtype=bool)
    for row, (right, left) in enumerate(phases_deg):
        values[0, 0, row] = cmath.exp(1j * math.radians(right))
        values[1, 0, row] = cmath.exp(1j * math.radians(left))
    for row, pol in flagged:
        flags[pol, 0, row] = True
    return np.arange(rows), values, flags


class PhaseTableTests(unittest.TestCase):

    def test_phases_are_read_per_antenna_and_polarization(self):
        antennas, values, flags = columns(
            [(0, 0), (12, -8), (-20, 15), (3, 2), (-1, 1)], flagged=[(3, 1)]
        )

        table = TOOL.phase_table(antennas, values, flags, NAMES)

        self.assertEqual([entry["antenna"] for entry in table], NAMES)
        self.assertAlmostEqual(table[1]["phase_deg"][0], 12.0)
        self.assertAlmostEqual(table[2]["phase_deg"][1], 15.0)
        self.assertIsNone(table[3]["phase_deg"][1])

    def test_only_the_first_solution_of_an_antenna_is_taken(self):
        antennas, values, flags = columns([(5, 5), (9, 9), (40, 40)])
        antennas = np.asarray([0, 1, 0])

        table = TOOL.phase_table(antennas, values, flags, NAMES)

        self.assertEqual(len(table), 2)
        self.assertAlmostEqual(table[0]["phase_deg"][0], 5.0)

    def test_offsets_of_several_degrees_are_called_offsets(self):
        table = TOOL.phase_table(
            *columns([(0, 0), (12, -8), (-20, 15), (9, 7), (-11, 6)]), NAMES
        )

        summary = TOOL.phase_summary(table)
        text = " ".join(TOOL.reading(summary))

        self.assertEqual(summary["solved_antennas"], 5)
        self.assertEqual(summary["beyond_5_deg"], 4)
        self.assertAlmostEqual(summary["largest_deg"], 20.0)
        self.assertGreater(summary["rms_deg"], 8.0)
        self.assertIn("fixed offsets", text)

    def test_phases_near_zero_are_not(self):
        table = TOOL.phase_table(
            *columns([(0, 0), (1, -1), (-2, 1), (0.5, 0.2), (-1, 1)]), NAMES
        )

        text = " ".join(TOOL.reading(TOOL.phase_summary(table)))

        self.assertIn("do not explain", text)

    def test_a_table_with_nothing_solved_applies_nothing(self):
        antennas, values, flags = columns(
            [(0, 0), (5, 5)], flagged=[(0, 0), (0, 1), (1, 0), (1, 1)]
        )
        summary = TOOL.phase_summary(
            TOOL.phase_table(antennas, values, flags, NAMES)
        )

        self.assertEqual(summary["solved_antennas"], 0)
        self.assertIn("nothing was applied", TOOL.reading(summary)[0])

    def test_the_corrected_data_get_their_own_name(self):
        self.assertEqual(
            str(TOOL.output_path("wsclean_test/G1_scan15.ms/", "_ph")),
            "wsclean_test/G1_scan15_ph.ms",
        )


def series_columns(tracks, amplitudes=None, failed=()):
    """Gain-table columns for several solutions per antenna.

    ``tracks`` gives, per antenna, the phase in degrees at each time.
    """

    steps = len(tracks[0])
    rows = len(tracks) * steps
    antennas = np.repeat(np.arange(len(tracks)), steps)
    times = np.tile(np.arange(steps) * 60.0, len(tracks))
    values = np.ones((2, 1, rows), dtype=complex)
    flags = np.zeros((2, 1, rows), dtype=bool)
    for antenna, track in enumerate(tracks):
        gain = 1.0 if amplitudes is None else amplitudes[antenna]
        for step, phase in enumerate(track):
            row = antenna * steps + step
            values[:, 0, row] = gain * cmath.exp(1j * math.radians(phase))
            if (antenna, step) in failed:
                flags[:, 0, row] = True
    return antennas, times, values, flags


class SeriesTests(unittest.TestCase):

    def test_an_offset_and_a_wander_are_told_apart(self):
        steady = [0, 1, -1, 0, 1, -1]
        offset = [30, 31, 29, 30, 31, 29]
        wandering = [-40, -10, 25, 50, 20, -45]
        columns_ = series_columns([steady, offset, wandering])

        series = TOOL.gain_series(*columns_, NAMES)

        first = [entry["polarizations"][0] for entry in series]
        self.assertAlmostEqual(first[0]["mean_phase_deg"], 0.0, delta=0.5)
        self.assertAlmostEqual(first[1]["mean_phase_deg"], 30.0, delta=0.5)
        self.assertLess(first[1]["phase_rms_deg"], 2.0)
        self.assertGreater(first[2]["phase_rms_deg"], 30.0)
        text = " ".join(TOOL.series_reading(series, amplitude=False))
        self.assertIn("offset by more than 5 deg: C01", text)
        self.assertIn("more than 10 deg: C02", text)
        self.assertIn("cannot follow", text)

    def test_the_largest_single_solution_is_given_beside_the_mean(self):
        # One antenna at +80 degrees in one scan and -80 in the next.
        antennas, times, values, flags = series_columns([[0, 0], [80, -80]])
        values[:, 0, 2] *= 1.3         # its first scan
        values[:, 0, 3] *= 0.9         # its second

        series = TOOL.gain_series(antennas, times, values, flags, NAMES)

        pol = series[1]["polarizations"][0]
        self.assertAlmostEqual(pol["mean_phase_deg"], 0.0, places=6)
        self.assertAlmostEqual(pol["largest_phase_deg"], 80.0)
        self.assertAlmostEqual(pol["largest_amplitude_step"], 0.3)

    def test_phases_that_wrap_are_averaged_on_the_circle(self):
        series = TOOL.gain_series(
            *series_columns([[175, -175, 178, -178]]), NAMES
        )

        pol = series[0]["polarizations"][0]
        self.assertAlmostEqual(abs(pol["mean_phase_deg"]), 180.0, delta=1.0)
        self.assertLess(pol["phase_rms_deg"], 5.0)

    def test_failed_solutions_are_counted_not_averaged(self):
        series = TOOL.gain_series(
            *series_columns([[10, 10, 90, 10]], failed={(0, 2)}), NAMES
        )

        pol = series[0]["polarizations"][0]
        self.assertEqual(pol["failed"], 1)
        self.assertAlmostEqual(pol["mean_phase_deg"], 10.0, delta=0.1)

    def test_amplitudes_far_from_the_median_are_named(self):
        series = TOOL.gain_series(
            *series_columns([[0, 0], [0, 0], [0, 0]],
                            amplitudes=[1.0, 0.8, 1.03]),
            NAMES,
        )

        text = " ".join(TOOL.series_reading(series, amplitude=True))

        self.assertIn("Amplitudes more than 10% from the median: C01.", text)

    def test_an_antenna_with_no_solution_is_reported_as_failed(self):
        series = TOOL.gain_series(
            *series_columns([[5, 5], [0, 0]], failed={(0, 0), (0, 1)}), NAMES
        )

        self.assertEqual(series[0]["polarizations"][0], {"failed": 2})
        self.assertIn(
            "1 of 2 antennas solved",
            TOOL.series_reading(series, amplitude=False)[0],
        )


class DelayTests(unittest.TestCase):

    def table(self, delays):
        values = np.zeros((2, 1, len(delays)))
        values[:, 0, :] = delays
        flags = np.zeros((2, 1, len(delays)), dtype=bool)
        return TOOL.delay_table(np.arange(len(delays)), values, flags, NAMES)

    def test_a_delay_is_stated_as_what_it_does_across_the_band(self):
        # 0.3 ns over 280 MHz: 0.36 * 0.3 * 280 = 30 degrees.
        text = " ".join(TOOL.delay_reading(self.table([0, 0.3, -0.05]), 280.0))

        self.assertIn("+0.30 ns on C01", text)
        self.assertIn("30 deg across 280 MHz", text)
        self.assertIn("C01 (+30)", text)
        self.assertNotIn("C02 (", text)

    def test_small_delays_are_said_to_be_small(self):
        text = " ".join(TOOL.delay_reading(self.table([0, 0.02, -0.03]), 280.0))

        self.assertIn("delays are not what is left", text)

    def test_flagged_delays_are_left_out(self):
        values = np.zeros((2, 1, 2))
        flags = np.ones((2, 1, 2), dtype=bool)
        table = TOOL.delay_table(np.arange(2), values, flags, NAMES)

        self.assertEqual(table[0]["delay_ns"], [None, None])
        self.assertEqual(
            TOOL.delay_reading(table, 280.0),
            ["No antenna has a delay solution."],
        )


def band_columns(bins=12, bow_deg=20.0, slope_deg=0.0, flagged=(0, 11)):
    """Bandpass-table columns: C02 bows across the band, W03 has a ripple.

    Indexed [polarization, bin, antenna]. The first and last bins are
    flagged, as the band's edges are.
    """

    values = np.ones((2, bins, len(NAMES)), dtype=complex)
    flags = np.zeros((2, bins, len(NAMES)), dtype=bool)
    centre = (bins - 1) / 2.0
    for chan in range(bins):
        offset = (chan - centre) / centre
        bow = bow_deg * (1.0 - offset * offset) + slope_deg * offset
        values[:, chan, 2] = cmath.exp(1j * math.radians(bow))
        values[:, chan, 4] = 1.0 + 0.2 * (-1) ** chan
    for chan in flagged:
        flags[:, chan, :] = True
    return values, flags


class BandpassTests(unittest.TestCase):

    def statistics(self, **kwargs):
        values, flags = band_columns(**kwargs)
        table = TOOL.bandpass_table(
            np.arange(len(NAMES)), values, flags, NAMES
        )
        return table, TOOL.bandpass_statistics(table)

    def test_each_bin_of_each_antenna_is_read_and_flagged_ones_left_out(self):
        table, statistics = self.statistics()

        self.assertEqual([entry["antenna"] for entry in table], NAMES)
        bins = table[2]["polarizations"][0]["phase_deg"]
        self.assertIsNone(bins[0])
        self.assertIsNone(bins[-1])
        self.assertAlmostEqual(max(b for b in bins if b is not None), 19.8,
                               places=0)
        self.assertEqual(statistics[2]["polarizations"][0]["bins"], 10)

    def test_a_phase_that_bows_is_not_taken_for_a_delay(self):
        _, statistics = self.statistics(bow_deg=20.0)

        bowed = statistics[2]["polarizations"][0]
        self.assertGreater(bowed["phase_rms_deg"], 3.0)
        # A bow is symmetric about the band centre: no line fits it.
        self.assertAlmostEqual(bowed["slope_deg_across_band"], 0.0, places=6)
        self.assertAlmostEqual(
            bowed["phase_rms_after_line_deg"], bowed["phase_rms_deg"]
        )
        flat = statistics[1]["polarizations"][0]
        self.assertAlmostEqual(flat["phase_rms_deg"], 0.0, places=6)

    def test_a_straight_line_is_what_a_delay_would_have_taken(self):
        _, statistics = self.statistics(bow_deg=0.0, slope_deg=30.0)

        sloped = statistics[2]["polarizations"][1]
        self.assertGreater(sloped["phase_rms_deg"], 10.0)
        self.assertAlmostEqual(sloped["phase_rms_after_line_deg"], 0.0,
                               places=6)
        self.assertAlmostEqual(
            sloped["slope_deg_across_band"], 2 * 30.0 * 9 / 11, places=4
        )

    def test_amplitudes_that_ripple_are_measured(self):
        _, statistics = self.statistics()

        self.assertAlmostEqual(
            statistics[4]["polarizations"][0]["amplitude_rms"], 0.2
        )
        self.assertAlmostEqual(
            statistics[0]["polarizations"][0]["amplitude_rms"], 0.0
        )

    def test_the_reading_names_what_bends_and_what_ripples(self):
        _, statistics = self.statistics()

        lines = TOOL.bandpass_reading(statistics, 30)

        self.assertIn("5 antennas solved, a gain for every 30 channels",
                      lines[0])
        self.assertIn("(C02)", lines[0])
        self.assertEqual(lines[1], "More than 3 deg rms: C02.")
        self.assertIn("more than 3 deg: C02.", lines[2])
        self.assertIn("more than 5%: W03.", lines[3])

    def test_an_antenna_with_too_few_bins_has_no_statistics(self):
        values, flags = band_columns()
        flags[:, :, 3] = True
        flags[:, 5, 3] = False
        table = TOOL.bandpass_table(
            np.arange(len(NAMES)), values, flags, NAMES
        )

        statistics = TOOL.bandpass_statistics(table)

        self.assertEqual(statistics[3]["polarizations"][0], {"bins": 1})
        self.assertNotIn("E02", " ".join(TOOL.bandpass_reading(statistics, 30)))

    def test_nothing_solved_is_said_so(self):
        self.assertEqual(
            TOOL.bandpass_reading([], 30),
            ["No antenna has a solution: nothing was applied."],
        )

    def test_the_shared_amplitude_is_the_median_over_antennas(self):
        values, flags = band_columns()
        values[:, 5, :] *= 1.03
        table = TOOL.bandpass_table(
            np.arange(len(NAMES)), values, flags, NAMES
        )

        common = TOOL.common_amplitude(table)

        # W03's ripple of 20% does not move the median of five antennas.
        self.assertAlmostEqual(common[0][5], 1.03)
        self.assertAlmostEqual(common[1][4], 1.0)
        self.assertIsNone(common[0][0])
        line = TOOL.common_reading(common)[0]
        self.assertIn("1.000 to 1.030", line)
        self.assertEqual(TOOL.common_amplitude([]), [])
        self.assertEqual(TOOL.common_reading([[None, None]]), [])


class FlowTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.ms = self.root / "G1_scan15.ms"
        (self.ms / "ANTENNA").mkdir(parents=True)
        self.calls = []
        self.ms_columns = ["DATA", "MODEL_DATA"]
        self.phases = [(0, 0), (12, -8), (-20, 15), (9, 7), (-11, 6)]
        self.flagged = []
        self.series = None
        self.scan = 15
        self.delays = [0.0, 0.02, -0.31, 0.05, 0.12]
        self.band = band_columns()
        self.written = []
        test = self

        class Table:
            def open(self, path, nomodify=True):
                self.path = Path(path)
                self.nomodify = nomodify

            def putcol(self, name, values):
                assert not self.nomodify, "table was opened read-only"
                test.written.append((self.path.name, name, values))

            def close(self):
                pass

            def colnames(self):
                return list(test.ms_columns)

            def getcol(self, name):
                if self.path.name == "ANTENNA":
                    return np.asarray(NAMES)
                if self.path.name.endswith(".bandpass.cal"):
                    values, flags = test.band
                    return {"CPARAM": values, "FLAG": flags,
                            "ANTENNA1": np.arange(len(NAMES))}[name]
                if self.path.name.endswith(".delay.cal"):
                    delays = np.zeros((2, 1, len(NAMES)))
                    delays[:, 0, :] = test.delays
                    flags = np.zeros((2, 1, len(NAMES)), dtype=bool)
                    return {"FPARAM": delays, "FLAG": flags,
                            "ANTENNA1": np.arange(len(NAMES))}[name]
                if test.series is not None:
                    antennas, times, values, flags = test.series
                else:
                    antennas, values, flags = columns(
                        test.phases, test.flagged
                    )
                    times = np.zeros(len(antennas))
                return {"CPARAM": values, "FLAG": flags,
                        "ANTENNA1": antennas, "TIME": times,
                        "SCAN_NUMBER": np.full(len(antennas), test.scan),
                        }[name]

        def gaincal(**kwargs):
            test.calls.append(("gaincal", kwargs))
            if kwargs.get("gaintype") == "K":
                if test.delays is not None:
                    Path(kwargs["caltable"]).mkdir()
                return
            if test.phases is not None:
                Path(kwargs["caltable"]).mkdir()

        def bandpass(**kwargs):
            test.calls.append(("bandpass", kwargs))
            if test.band is not None:
                Path(kwargs["caltable"]).mkdir()

        def applycal(**kwargs):
            test.calls.append(("applycal", kwargs))

        def split(**kwargs):
            test.calls.append(("split", kwargs))
            Path(kwargs["outputvis"]).mkdir()

        def flagdata(**kwargs):
            test.calls.append(("flagdata", kwargs))

        casatasks = types.ModuleType("casatasks")
        casatasks.flagdata = flagdata
        casatasks.gaincal, casatasks.applycal, casatasks.split = (
            gaincal, applycal, split
        )
        casatasks.bandpass = bandpass
        casatools = types.ModuleType("casatools")
        casatools.table = Table
        patch = mock.patch.dict(
            sys.modules, {"casatasks": casatasks, "casatools": casatools}
        )
        patch.start()
        self.addCleanup(patch.stop)

    def run_tool(self, *extra):
        with mock.patch("builtins.print"):
            TOOL.main(["--ms", str(self.ms), *extra])
        return json.loads(
            (self.root / "G1_scan15_ph.json").read_text()
        )

    def test_one_phase_per_antenna_is_solved_applied_and_split(self):
        report = self.run_tool()

        names = [name for name, _ in self.calls]
        self.assertEqual(names, ["gaincal", "applycal", "split"])
        solve = self.calls[0][1]
        self.assertEqual(solve["calmode"], "p")
        self.assertEqual(solve["solint"], "inf")
        self.assertEqual(solve["gaintype"], "G")
        self.assertEqual(solve["refant"], "C00,W03,E05,E02,E06")
        self.assertEqual(solve["refantmode"], "flex")
        apply = self.calls[1][1]
        self.assertEqual(apply["vis"], str(self.ms))
        self.assertEqual(apply["applymode"], "calonly")
        self.assertEqual(apply["interp"], ["nearest"])
        cut = self.calls[2][1]
        self.assertEqual(cut["datacolumn"], "corrected")
        self.assertEqual(cut["outputvis"], str(self.root / "G1_scan15_ph.ms"))

        self.assertEqual(report["summary"]["beyond_5_deg"], 4)
        self.assertEqual(
            report["corrected_measurement_set"],
            str(self.root / "G1_scan15_ph.ms"),
        )
        self.assertEqual(report["phases"][2]["antenna"], "C02")

    def test_a_shorter_interval_reports_and_interpolates(self):
        self.series = series_columns(
            [[0, 1, -1, 0], [30, 31, 29, 30], [-40, 10, 50, -30],
             [2, 2, 2, 2], [0, 0, 0, 0]]
        )

        report = self.run_tool("--solint", "1min")

        self.assertEqual(self.calls[0][1]["solint"], "1min")
        # Solutions within a scan are joined, not stepped.
        self.assertEqual(self.calls[1][1]["interp"], ["linear"])
        self.assertEqual(report["series"][2]["solutions"], 4)
        self.assertIn("C02", " ".join(report["reading"]))
        self.assertIn("cannot follow", " ".join(report["reading"]))

    def test_one_solution_per_scan_is_not_joined_to_the_next_scans(self):
        # Three scans, one solution each: C01 is offset in the second only.
        self.series = series_columns(
            [[0, 0, 0], [5, -80, 10], [0, 0, 0], [2, 2, 2], [0, 0, 0]]
        )

        self.run_tool()

        self.assertEqual(self.calls[0][1]["solint"], "inf")
        self.assertEqual(self.calls[1][1]["interp"], ["nearest"])

    def test_amplitudes_are_solved_without_moving_the_flux_scale(self):
        report = self.run_tool("--calmode", "ap", "--gaintype", "T")

        solve = self.calls[0][1]
        self.assertEqual(solve["calmode"], "ap")
        self.assertEqual(solve["gaintype"], "T")
        self.assertTrue(solve["solnorm"])
        self.assertEqual(solve["normtype"], "median")
        self.assertIn("Amplitudes", " ".join(report["reading"]))

    def test_short_baselines_can_be_kept_out_of_the_gain_solve(self):
        report = self.run_tool(
            "--calmode", "ap", "--uvrange", ">3klambda", "--delays"
        )

        delay, gain = self.calls[0][1], self.calls[1][1]
        self.assertEqual(gain["uvrange"], ">3klambda")
        # Not the delay solve: the central antennas' baselines to the
        # reference antenna are all short.
        self.assertNotIn("uvrange", delay)
        self.assertEqual(report["solve"]["uvrange"], ">3klambda")

    def test_by_default_every_baseline_is_used(self):
        self.run_tool()

        self.assertNotIn("uvrange", self.calls[0][1])

    def test_weights_follow_an_amplitude_correction(self):
        self.run_tool("--calmode", "ap", "--delays")

        apply = [kw for name, kw in self.calls if name == "applycal"][0]
        # Not for the delays, which change no amplitude.
        self.assertEqual(apply["calwt"], [False, True])

    def test_weights_are_left_alone_by_a_phase_correction(self):
        self.run_tool()

        apply = [kw for name, kw in self.calls if name == "applycal"][0]
        self.assertEqual(apply["calwt"], [False])

    def test_a_phase_solve_does_not_normalise(self):
        self.run_tool()

        self.assertNotIn("solnorm", self.calls[0][1])

    def test_delays_are_solved_first_and_the_phases_on_top(self):
        report = self.run_tool("--delays")

        delay, phase, apply, _ = self.calls
        self.assertEqual(delay[1]["gaintype"], "K")
        self.assertEqual(delay[1]["solint"], "inf")
        table = str(self.root / "G1_scan15_ph.delay.cal")
        self.assertEqual(delay[1]["caltable"], table)
        # The phase solve sees the data with the delays taken out.
        self.assertEqual(phase[1]["gaintable"], [table])
        self.assertEqual(phase[1]["calmode"], "p")
        # Both are applied, delays first.
        self.assertEqual(
            apply[1]["gaintable"],
            [table, str(self.root / "G1_scan15_ph.phase.cal")],
        )
        self.assertEqual(apply[1]["interp"], ["nearest", "nearest"])
        self.assertAlmostEqual(report["delays"][2]["delay_ns"][0], -0.31)
        text = " ".join(report["delay_reading"])
        self.assertIn("-0.31 ns on C02", text)
        self.assertIn("C02 (-31)", text)

    def test_a_delay_too_large_to_be_a_residual_is_not_applied(self):
        # C02 comes out at 17.85 ns, which turns its phase five times
        # round across the band; the others are a few tenths.
        self.delays = [0.0, 0.02, 17.85, 0.05, 0.12]

        report = self.run_tool("--delays")

        table, column, values = self.written[0]
        self.assertEqual((table, column),
                         ("G1_scan15_ph.delay.cal", "FPARAM"))
        np.testing.assert_allclose(values[:, 0, 2], 0.0)
        np.testing.assert_allclose(values[:, 0, 4], 0.12)
        self.assertEqual(
            [entry["antenna"] for entry in report["rejected_delays"]], ["C02"]
        )
        self.assertEqual(report["rejected_delays"][0]["delay_ns"],
                         [17.85, 17.85])
        self.assertEqual(report["delays"][2]["delay_ns"], [0.0, 0.0])
        # The antenna still gets its phase: the table is applied.
        self.assertEqual(self.calls[-2][0], "applycal")

    def test_the_limit_on_a_delay_can_be_raised(self):
        self.delays = [0.0, 0.02, 17.85, 0.05, 0.12]

        report = self.run_tool("--delays", "--max-delay-ns", "30")

        self.assertEqual(report["rejected_delays"], [])
        self.assertEqual(self.written, [])

    def test_an_amplitude_far_from_one_is_not_applied(self):
        antennas, values, flags = columns(self.phases)
        values[:, 0, 1] *= 0.2       # C01: a solve that went wrong
        values[:, 0, 3] *= 1.31      # E02: large, and believable
        self.series = (antennas, np.zeros(len(antennas)), values, flags)

        report = self.run_tool("--calmode", "ap")

        table, column, written = self.written[0]
        self.assertEqual((table, column),
                         ("G1_scan15_ph.phase.cal", "FLAG"))
        self.assertTrue(written[:, 0, 1].all())
        self.assertFalse(written[:, 0, 3].any())
        self.assertEqual(report["rejected_amplitudes"], ["C01"])
        self.assertEqual(
            report["series"][1]["polarizations"][0], {"failed": 1}
        )
        self.assertAlmostEqual(
            report["series"][3]["polarizations"][0]["mean_amplitude"], 1.31
        )

    def test_an_antenna_far_too_high_in_a_scan_is_mended_not_refused(self):
        # S04 on G0216-0604: beyond the old limit of 2.5 in one of two
        # scans, and real. Refused, it left a hole beside every source.
        antennas, values, flags = columns(self.phases)
        values[:, 0, 2] *= 3.2
        self.series = (antennas, np.zeros(len(antennas)), values, flags)

        report = self.run_tool("--calmode", "ap")

        self.assertEqual(report["rejected_amplitudes"], [])
        self.assertEqual(report["flagged_rejected"], [])
        self.assertEqual(self.written, [])
        self.assertNotIn("flagdata", [name for name, _ in self.calls])
        self.assertAlmostEqual(
            report["series"][2]["polarizations"][0]["mean_amplitude"], 3.2
        )

    def test_data_whose_amplitude_is_refused_are_flagged_for_that_scan(self):
        antennas, values, flags = columns(self.phases)
        values[:, 0, 2] *= 5.5
        self.series = (antennas, np.zeros(len(antennas)), values, flags)
        self.scan = 5

        report = self.run_tool("--calmode", "ap")

        self.assertEqual(report["rejected_amplitudes"], [NAMES[2]])
        self.assertEqual(report["flagged_rejected"],
                         [{"antenna": NAMES[2], "scan": 5}])
        # After the corrected data are written, and in them, not in the
        # data given.
        names = [name for name, _ in self.calls]
        self.assertEqual(names[-2:], ["split", "flagdata"])
        flag = self.calls[-1][1]
        self.assertEqual(flag["vis"], self.calls[-2][1]["outputvis"])
        self.assertEqual((flag["mode"], flag["antenna"], flag["scan"]),
                         ("manual", NAMES[2], "5"))
        self.assertFalse(flag["flagbackup"])

    def test_a_solve_towards_one_source_flags_none_of_the_fields_data(self):
        antennas, values, flags = columns(self.phases)
        values[:, 0, 2] *= 5.5
        self.series = (antennas, np.zeros(len(antennas)), values, flags)

        report = self.run_tool("--calmode", "ap", "--keep-unsolved")

        self.assertEqual(report["rejected_amplitudes"], [NAMES[2]])
        self.assertEqual(report["flagged_rejected"], [])
        self.assertNotIn("flagdata", [name for name, _ in self.calls])

    def test_within_a_scan_a_refused_amplitude_flags_nothing(self):
        # One solution of several in a scan is not the scan.
        antennas, values, flags = columns(self.phases)
        values[:, 0, 2] *= 5.5
        self.series = (antennas, np.zeros(len(antennas)), values, flags)

        report = self.run_tool("--calmode", "ap", "--solint", "1min")

        self.assertEqual(report["rejected_amplitudes"], [NAMES[2]])
        self.assertEqual(report["flagged_rejected"], [])

    def test_a_phase_solve_rejects_no_amplitude(self):
        report = self.run_tool()

        self.assertEqual(report["rejected_amplitudes"], [])
        self.assertEqual(self.written, [])

    def test_without_the_option_no_delay_is_solved(self):
        report = self.run_tool()

        self.assertEqual(
            [kw.get("gaintype") for name, kw in self.calls
             if name == "gaincal"], ["G"],
        )
        self.assertNotIn("gaintable", self.calls[0][1])
        self.assertIsNone(report["delays"])

    def test_a_delay_solve_that_writes_no_table_is_reported(self):
        self.delays = None

        with self.assertRaises(SystemExit) as raised:
            self.run_tool("--delays")

        self.assertIn("no delay table", str(raised.exception))

    def test_it_needs_the_model_column(self):
        self.ms_columns = ["DATA"]

        with self.assertRaises(SystemExit) as raised:
            self.run_tool()

        self.assertIn("--core-model", str(raised.exception))
        self.assertEqual(self.calls, [])

    def test_an_earlier_result_is_not_overwritten_unasked(self):
        self.run_tool()
        self.calls.clear()

        with self.assertRaises(SystemExit) as raised:
            self.run_tool()
        self.assertIn("--overwrite", str(raised.exception))
        self.assertEqual(self.calls, [])

        self.run_tool("--overwrite")
        self.assertEqual(len(self.calls), 3)

    def test_nothing_is_applied_when_nothing_was_solved(self):
        self.flagged = [(row, pol) for row in range(5) for pol in range(2)]

        report = self.run_tool()

        self.assertEqual([name for name, _ in self.calls], ["gaincal"])
        self.assertNotIn("corrected_measurement_set", report)

    def test_a_solve_that_writes_no_table_is_reported(self):
        self.phases = None

        with self.assertRaises(SystemExit) as raised:
            self.run_tool()

        self.assertIn("wrote no table", str(raised.exception))

    def run_band(self, *extra):
        with mock.patch("builtins.print"):
            TOOL.main(["--ms", str(self.ms), "--bandpass", "30",
                       "--suffix", "_bp", *extra])
        return json.loads((self.root / "G1_scan15_bp.json").read_text())

    def test_a_gain_per_channel_bin_is_solved_and_its_phases_applied(self):
        report = self.run_band()

        names = [name for name, _ in self.calls]
        self.assertEqual(names, ["bandpass", "applycal", "split"])
        solve = self.calls[0][1]
        self.assertEqual(solve["solint"], "inf,30ch")
        self.assertEqual(solve["bandtype"], "B")
        self.assertEqual(solve["combine"], "")
        self.assertFalse(solve["solnorm"])
        self.assertNotIn("gaintable", solve)
        # Only the phases go into the data: every gain is given unit size.
        self.assertEqual(len(self.written), 1)
        table, column, values = self.written[0]
        self.assertEqual((table, column), ("G1_scan15_bp.bandpass.cal",
                                           "CPARAM"))
        np.testing.assert_allclose(abs(values), 1.0)
        np.testing.assert_allclose(
            np.angle(values), np.angle(self.band[0])
        )
        apply = self.calls[1][1]
        self.assertEqual(
            apply["gaintable"], [str(self.root / "G1_scan15_bp.bandpass.cal")]
        )
        self.assertEqual(apply["interp"], ["nearest,linear"])
        self.assertEqual(apply["calwt"], [False])
        self.assertEqual(apply["applymode"], "calonly")
        self.assertEqual(
            self.calls[2][1]["outputvis"], str(self.root / "G1_scan15_bp.ms")
        )
        self.assertFalse(report["solve"]["amplitudes_applied"])
        self.assertIn("phases only", report["reading"][-1])
        # What was solved is reported as solved, ripple included.
        self.assertAlmostEqual(
            report["statistics"][4]["polarizations"][0]["amplitude_rms"], 0.2
        )

    def test_amplitudes_across_the_band_are_applied_only_when_asked(self):
        report = self.run_band("--bandpass-amplitude")

        apply = self.calls[1][1]
        self.assertEqual(apply["calwt"], [True])
        self.assertTrue(report["solve"]["amplitudes_applied"])
        self.assertIn("phases and amplitudes", report["reading"][-1])
        # What a whole scan could not solve is flagged, not guessed.
        self.assertEqual(apply["interp"], ["nearest,linearflag"])
        self.assertEqual(apply["applymode"], "calflag")

    def test_what_all_antennas_share_is_the_model_and_is_not_applied(self):
        values, flags = self.band
        # Every antenna 4% high in bin 3 and 4% low in bin 4: a step of
        # the model. C01 alone is at twice the others from bin 8 on.
        values[:, 3, :] *= 1.04
        values[:, 4, :] *= 0.96
        values[:, 8:, 1] *= 2.0

        report = self.run_band("--bandpass-amplitude")

        self.assertAlmostEqual(report["common_amplitude"][0][3], 1.04)
        self.assertAlmostEqual(report["common_amplitude"][0][4], 0.96)
        self.assertIsNone(report["common_amplitude"][0][0])
        self.assertIn("common to all antennas", " ".join(report["reading"]))
        self.assertIn("never", " ".join(report["reading"]))
        table, column, written = self.written[0]
        self.assertEqual(column, "CPARAM")
        # The shared step is gone from what is applied ...
        np.testing.assert_allclose(abs(written[0, 3, 0]), 1.0)
        np.testing.assert_allclose(abs(written[0, 4, 0]), 1.0)
        # ... the antenna's own factor of two is kept, with its phase.
        np.testing.assert_allclose(abs(written[0, 9, 1]), 2.0)
        np.testing.assert_allclose(abs(written[0, 5, 1]), 1.0)
        np.testing.assert_allclose(
            np.angle(written[:, 5, 2]), np.angle(values[:, 5, 2])
        )

    def test_towards_one_source_unsolved_channels_are_not_flagged(self):
        report = self.run_band("--bandpass-amplitude", "--keep-unsolved")

        apply = self.calls[1][1]
        self.assertEqual(apply["interp"], ["nearest,linear"])
        self.assertEqual(apply["applymode"], "calonly")
        # The amplitudes are still corrected, and the weights follow.
        self.assertEqual(apply["calwt"], [True])
        self.assertTrue(report["solve"]["amplitudes_applied"])

    def test_a_bin_far_from_what_the_antennas_share_is_flagged(self):
        values, flags = self.band
        values[:, 5, 1] *= 3.0       # one bin of C01, three times the rest

        report = self.run_band("--bandpass-amplitude")

        flagged = [w for w in self.written if w[1] == "FLAG"]
        self.assertEqual(len(flagged), 1)
        self.assertTrue(flagged[0][2][:, 5, 1].all())
        self.assertFalse(flagged[0][2][:, 6, 1].any())
        self.assertEqual(report["rejected_bins"], {"C01": 2})
        # It is not counted in that antenna's run across the band either.
        self.assertIsNone(
            report["bandpass"][1]["polarizations"][0]["amplitude"][5]
        )

    def test_the_phase_only_correction_leaves_unsolved_data_alone(self):
        self.run_band()

        apply = self.calls[1][1]
        self.assertEqual(apply["interp"], ["nearest,linear"])
        self.assertEqual(apply["applymode"], "calonly")

    def test_an_earlier_solve_is_taken_out_first_and_applied_with_it(self):
        earlier = self.root / "G1_scan15_t1.phase.cal"
        earlier.mkdir()

        report = self.run_band("--prior", str(earlier))

        solve = self.calls[0][1]
        self.assertEqual(solve["gaintable"], [str(earlier)])
        self.assertEqual(solve["interp"], ["linear"])
        apply = self.calls[1][1]
        self.assertEqual(
            apply["gaintable"],
            [str(earlier), str(self.root / "G1_scan15_bp.bandpass.cal")],
        )
        self.assertEqual(apply["interp"], ["linear", "nearest,linear"])
        # The earlier table carries amplitudes, and weights follow those.
        self.assertEqual(apply["calwt"], [True, False])
        self.assertEqual(report["prior"], [str(earlier)])

    def test_an_earlier_table_that_is_not_there_is_said_at_once(self):
        with self.assertRaises(SystemExit) as stopped:
            self.run_band("--prior", str(self.root / "absent.cal"))

        self.assertIn("No calibration table", str(stopped.exception))
        self.assertEqual(self.calls, [])

    def test_a_bandpass_solve_that_writes_no_table_is_reported(self):
        self.band = None

        with self.assertRaises(SystemExit) as stopped:
            self.run_band()

        self.assertIn("no solution was possible", str(stopped.exception))

    def test_the_bandpass_solve_needs_the_model_column(self):
        self.ms_columns = ["DATA"]

        with self.assertRaises(SystemExit) as stopped:
            self.run_band()

        self.assertIn("MODEL_DATA", str(stopped.exception))

    def test_a_second_copy_of_the_same_solve_is_refused(self):
        import os
        import socket

        lock = self.root / "G1_scan15_bp.lock"
        # This test's own process stands for a solve that is still going.
        lock.write_text(json.dumps({
            "pid": os.getpid(), "host": socket.gethostname(),
            "since": "2026-10-09 16:12:41 UTC",
        }))

        with self.assertRaises(SystemExit) as stopped:
            self.run_band("--overwrite")

        self.assertIn("Another run of this solve", str(stopped.exception))
        self.assertEqual(self.calls, [])
        self.assertTrue(lock.is_file())

    def test_a_solve_with_another_suffix_is_not_held_up(self):
        import os
        import socket

        (self.root / "G1_scan15_bp.lock").write_text(json.dumps(
            {"pid": os.getpid(), "host": socket.gethostname()}
        ))

        self.run_tool()

        self.assertEqual(self.calls[0][0], "gaincal")

    def test_a_claim_left_by_a_solve_that_is_gone_is_taken_over(self):
        import socket

        lock = self.root / "G1_scan15_bp.lock"
        lock.write_text(json.dumps(
            {"pid": 2 ** 22 + 4321, "host": socket.gethostname()}
        ))

        with mock.patch.object(TOOL, "process_alive", lambda pid: False):
            self.run_band()

        self.assertEqual(self.calls[0][0], "bandpass")
        self.assertFalse(lock.exists())

    def test_the_claim_is_released_when_the_solve_fails(self):
        self.band = None

        with self.assertRaises(SystemExit):
            self.run_band()

        self.assertFalse((self.root / "G1_scan15_bp.lock").exists())

    def test_bandpass_options_are_refused_without_the_bandpass(self):
        with self.assertRaises(SystemExit), mock.patch("sys.stderr"):
            TOOL.main(["--ms", str(self.ms), "--bandpass-amplitude"])

        self.assertEqual(self.calls, [])


if __name__ == "__main__":
    unittest.main()
