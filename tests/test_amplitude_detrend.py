"""The cores-only amplitude round keeps antenna faults but not the model's bias.

The case is G0212-0335 in Band 5 (2026-10-22): the first cores-only
amplitude solve put the central square at 1.03-1.38 and the far arm
antennas at 0.33-0.6. The target's core rose by 84%, and 40% of its flux was
lost to a bowl of negatives around it.
"""

import importlib.util
import math
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load(name, module_name):
    spec = importlib.util.spec_from_file_location(module_name, SOURCE_DIR / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SOLVE = load("experiment_phase_offsets.py", "test_detrend_phase_offsets")
RECIPE = load("selfcal_recipe.py", "test_detrend_selfcal_recipe")


def gmrt_like_array():
    """14 antennas within 1 km of the middle, 12 along arms out to 14 km."""

    positions, names = [], []
    rng = np.random.default_rng(4)
    for index in range(14):
        angle = rng.uniform(0, 2 * math.pi)
        radius = rng.uniform(150.0, 900.0)
        positions.append([radius * math.cos(angle), radius * math.sin(angle), 0.0])
        names.append(f"C{index:02d}")
    for arm, angle in (("E", 0.3), ("S", 2.4), ("W", 4.4)):
        for index, radius in enumerate((2000.0, 5000.0, 9000.0, 14000.0)):
            positions.append(
                [radius * math.cos(angle), radius * math.sin(angle), 0.0]
            )
            names.append(f"{arm}{index + 2:02d}")
    return np.array(positions), names


def gains(amplitude_of, antennas, phases=None):
    """CPARAM of two hands, one channel, one row per antenna."""

    values = np.zeros((2, 1, len(antennas)), dtype=complex)
    for row, antenna in enumerate(antennas):
        phase = 0.0 if phases is None else phases[row]
        values[:, 0, row] = amplitude_of(antenna) * np.exp(1j * phase)
    return values


class RadialTrendTests(unittest.TestCase):

    def setUp(self):
        self.positions, self.names = gmrt_like_array()
        self.antennas = np.arange(len(self.names))
        centre = np.median(self.positions, axis=0)
        self.distance = np.linalg.norm(self.positions - centre, axis=1)
        self.flags = np.zeros((2, 1, len(self.names)), dtype=bool)

    def biased(self, antenna):
        # The model's doing: 1.3 at the middle, falling with log distance.
        value = 1.3 - 0.45 * math.log10(max(self.distance[antenna], 100.0) / 300.0)
        # One antenna's own fault, as S02 had in Band 4.
        return value * (1.6 if self.names[antenna] == "C05" else 1.0)

    def test_a_run_with_distance_is_found_and_divided_out(self):
        phases = np.linspace(-1.0, 1.0, len(self.names))
        values = gains(self.biased, self.antennas, phases)

        trend = SOLVE.radial_trend(
            self.antennas, values, self.flags, self.positions, self.names
        )
        self.assertGreater(trend["ratio"], 1.3)
        self.assertLess(trend["slope_per_decade"], -0.3)

        fixed = SOLVE.detrend_amplitudes(self.antennas, values, trend)
        size = np.abs(fixed[0, 0])
        inner = np.median(size[self.distance <= np.median(self.distance)])
        outer = np.median(size[self.distance > np.median(self.distance)])
        self.assertAlmostEqual(inner / outer, 1.0, delta=0.05)
        # The faulty antenna keeps its own excess over its neighbours.
        fault = self.names.index("C05")
        self.assertAlmostEqual(size[fault] / inner, 1.6, delta=0.15)
        # The phases are not touched.
        np.testing.assert_allclose(np.angle(fixed), np.angle(values), atol=1e-12)

    def test_antenna_faults_alone_make_no_run(self):
        values = gains(
            lambda antenna: 1.6 if self.names[antenna] in ("C05", "W11") else 1.0,
            self.antennas,
        )
        trend = SOLVE.radial_trend(
            self.antennas, values, self.flags, self.positions, self.names
        )
        self.assertAlmostEqual(trend["ratio"], 1.0, delta=0.02)

    def test_flagged_solutions_and_too_few_antennas(self):
        values = gains(self.biased, self.antennas)
        flags = self.flags.copy()
        flags[:, :, 20:] = True
        trend = SOLVE.radial_trend(
            self.antennas, values, flags, self.positions, self.names
        )
        self.assertEqual(len(trend["antennas"]), 20)
        self.assertIsNone(SOLVE.radial_trend(
            self.antennas[:5], values[:, :, :5], self.flags[:, :, :5],
            self.positions, self.names,
        ))


# G0212-0335, Band 5, amplitude_cores1 of 2026-10-23: each antenna's median
# amplitude as CASA normalised it, and the line radial_trend drew through
# them. The median of the line is 0.85, of the antennas 1.07.
G0212_CORES1 = {
    "C00": (1.294, 1.117), "C01": (1.135, 1.25), "C02": (1.078, 1.477),
    "C04": (1.241, 1.154), "C06": (1.379, 1.349), "C08": (1.214, 1.184),
    "C09": (1.106, 1.415), "C10": (1.055, 1.145), "C11": (1.201, 1.128),
    "C13": (1.1, 1.003), "C14": (1.104, 1.098), "E02": (1.089, 0.848),
    "E03": (0.927, 0.753), "E04": (1.079, 0.659), "E05": (0.578, 0.612),
    "E06": (0.462, 0.578), "S01": (0.843, 0.849), "S02": (0.899, 0.774),
    "S03": (0.636, 0.7), "S04": (0.717, 0.636), "S06": (0.422, 0.563),
    "W01": (1.104, 0.953), "W02": (0.785, 0.824), "W03": (0.572, 0.72),
    "W04": (0.392, 0.65), "W05": (0.352, 0.6),
}


class FluxScaleTests(unittest.TestCase):
    """Dividing the line out must not move the flux scale (2026-10-23).

    Scaled to the median of the line, G0212-0335's antennas came out at a
    median of 0.83: every visibility was turned up by 1.6, and its field's
    Stokes V noise rose from 22 to 28 uJy/beam.
    """

    def setUp(self):
        names = list(G0212_CORES1)
        self.antennas = np.arange(len(names))
        amplitude = np.array([G0212_CORES1[name][0] for name in names])
        self.values = np.zeros((2, 1, len(names)), dtype=complex)
        self.values[:, 0, :] = amplitude
        self.flags = np.zeros(self.values.shape, dtype=bool)
        self.trend = {"line": {
            index: G0212_CORES1[name][1] for index, name in enumerate(names)
        }}

    def test_the_median_antenna_stays_at_one(self):
        fixed = SOLVE.detrend_amplitudes(self.antennas, self.values, self.trend)
        fixed, medians = SOLVE.normalise_amplitudes(fixed, self.flags)

        size = np.abs(fixed)
        self.assertAlmostEqual(float(np.median(size[0])), 1.0, places=6)
        self.assertAlmostEqual(float(np.median(size[1])), 1.0, places=6)
        self.assertEqual(len(medians), 2)
        # Each antenna keeps its own departure from the line.
        own = np.abs(self.values[0, 0]) / np.array(
            [self.trend["line"][i] for i in self.antennas])
        np.testing.assert_allclose(size[0, 0] * medians[0], own, rtol=1e-3)
        # Scaled to the median of the line, as before, the median antenna
        # came out at 0.83: the visibilities were turned up by 1/0.83^2.
        line = np.median(list(self.trend["line"].values()))
        before = SOLVE.detrend_amplitudes(
            self.antennas, self.values, self.trend) * line
        self.assertLess(float(np.median(np.abs(before[0]))), 0.87)

    def test_each_hand_and_only_unflagged_solutions_count(self):
        values = self.values.copy()
        values[1] *= 1.2
        flags = self.flags.copy()
        values[0, 0, :5] = 40.0
        flags[0, 0, :5] = True
        fixed, medians = SOLVE.normalise_amplitudes(values, flags)

        size = np.abs(fixed)
        self.assertAlmostEqual(float(np.median(size[0][~flags[0]])), 1.0, places=6)
        self.assertAlmostEqual(float(np.median(size[1])), 1.0, places=6)
        self.assertAlmostEqual(medians[1] / medians[0],
                               1.2 * np.median(np.abs(self.values[0, 0])) /
                               np.median(np.abs(self.values[0, 0, 5:])), places=3)
        np.testing.assert_allclose(np.angle(fixed), np.angle(values), atol=1e-12)

    def test_the_round_says_when_it_was_set_back(self):
        found = RECIPE.solve_findings({
            "solve": {"calmode": "ap"}, "summary": {"solved_antennas": 26},
            "corrected_measurement_set": "x.ms", "series": [],
            "amplitude_medians": [0.83, 0.84],
        })
        self.assertEqual(found["rescaled"], [0.83, 0.84])
        line = RECIPE.reading({"rounds": [{"id": "amplitude_cores1",
                                           "solve": found}]})[0]
        self.assertIn("median amplitude 0.83/0.84 set back to 1", line)
        quiet = RECIPE.solve_findings({
            "solve": {"calmode": "ap"}, "summary": {"solved_antennas": 26},
            "series": [], "amplitude_medians": [1.004, 0.995],
        })
        self.assertNotIn("rescaled", quiet)
        band = RECIPE.bandpass_findings({
            "statistics": [], "amplitude_medians": [1.08, 1.07],
        })
        self.assertEqual(band["rescaled"], [1.08, 1.07])


class RecipeTests(unittest.TestCase):

    def plan(self, *extra):
        captured = []
        with mock.patch.object(
            RECIPE, "run", lambda args, plan=None: captured.append((args, plan))
        ):
            RECIPE.main([
                "--target", "G1", "--ms", "data.ms", "--imsize", "5120",
                "--cell", "0.5", "--channel-range", "31", "993",
                "--casa", "python3", *extra,
            ])
        return captured[0]

    def test_every_amplitude_solve_and_the_bandpass_divide_the_run_out(self):
        # The whole-field rounds brought G0212-0335's bias back (central
        # square over arms 1.32, 1.19, 1.09 in turn): finding 148.
        args, plan = self.plan("--amplitude-cores")
        rounds = {stage["id"]: stage for stage in plan}
        for name in ("amplitude_cores", "amplitude", "bandpass"):
            with self.subTest(round=name):
                command = RECIPE.solve_command(args, "a.ms", rounds[name]["solve"])
                self.assertEqual(
                    command[command.index("--detrend-radius") + 1], "0.05")
        phase = RECIPE.solve_command(args, "a.ms", rounds["phase"]["solve"])
        self.assertNotIn("--detrend-radius", phase)

    def test_the_correction_can_be_switched_off(self):
        args, plan = self.plan("--amplitude-detrend", "0")
        rounds = {stage["id"]: stage for stage in plan}
        for name in ("amplitude", "bandpass"):
            command = RECIPE.solve_command(args, "a.ms", rounds[name]["solve"])
            self.assertNotIn("--detrend-radius", command)

    def test_the_bandpass_round_says_what_was_divided_out(self):
        found = RECIPE.bandpass_findings({
            "statistics": [], "corrected_measurement_set": "x.ms",
            "radial_trend": {"divided_out": True, "inner_median": 0.97,
                             "outer_median": 1.04},
        })
        self.assertEqual(found["radial_trend"], [0.97, 1.04])

    def test_the_round_says_what_was_divided_out(self):
        found = RECIPE.solve_findings({
            "solve": {"calmode": "ap"},
            "summary": {"solved_antennas": 26},
            "corrected_measurement_set": "x.ms",
            "series": [],
            "radial_trend": {"divided_out": True, "inner_median": 1.15,
                             "outer_median": 0.7},
        })
        self.assertEqual(found["radial_trend"], [1.15, 0.7])
        quiet = RECIPE.solve_findings({
            "solve": {"calmode": "ap"}, "summary": {"solved_antennas": 26},
            "series": [], "radial_trend": {"divided_out": False},
        })
        self.assertNotIn("radial_trend", quiet)



class SolveFlowTests(unittest.TestCase):
    """The solve script, end to end, with stand-ins for CASA."""

    def setUp(self):
        import json
        import sys
        import tempfile
        import types

        self.json = json
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.ms = self.root / "round8.ms"
        (self.ms / "ANTENNA").mkdir(parents=True)
        self.positions, self.names = gmrt_like_array()
        count = len(self.names)
        centre = np.median(self.positions, axis=0)
        distance = np.linalg.norm(self.positions - centre, axis=1)
        amplitude = np.array([
            1.3 - 0.45 * math.log10(max(d, 100.0) / 300.0) for d in distance
        ])
        self.values = np.zeros((2, 1, count), dtype=complex)
        self.values[:, 0, :] = amplitude
        self.calls, self.written = [], []
        test = self

        class Table:
            def open(self, path, nomodify=True):
                self.path, self.nomodify = Path(path), nomodify

            def close(self):
                pass

            def colnames(self):
                return ["DATA", "MODEL_DATA"]

            def putcol(self, name, values):
                test.written.append((name, np.array(values)))

            def getcol(self, name):
                if self.path.name == "ANTENNA":
                    return {"NAME": np.array(test.names),
                            "POSITION": test.positions.T}[name]
                return {
                    "CPARAM": test.values.copy(),
                    "FLAG": np.zeros((2, 1, count), dtype=bool),
                    "ANTENNA1": np.arange(count),
                    "TIME": np.zeros(count),
                    "SCAN_NUMBER": np.full(count, 3),
                }[name]

        def gaincal(**kwargs):
            test.calls.append(("gaincal", kwargs))
            Path(kwargs["caltable"]).mkdir()

        def applycal(**kwargs):
            test.calls.append(("applycal", kwargs))

        def split(**kwargs):
            test.calls.append(("split", kwargs))
            Path(kwargs["outputvis"]).mkdir()

        def flagdata(**kwargs):
            test.calls.append(("flagdata", kwargs))

        casatasks = types.ModuleType("casatasks")
        casatasks.gaincal, casatasks.applycal = gaincal, applycal
        casatasks.split, casatasks.flagdata = split, flagdata
        casatools = types.ModuleType("casatools")
        casatools.table = Table
        patch = mock.patch.dict(
            sys.modules, {"casatasks": casatasks, "casatools": casatools}
        )
        patch.start()
        self.addCleanup(patch.stop)

    def run_solve(self, *extra):
        with mock.patch("builtins.print"):
            SOLVE.main(["--ms", str(self.ms), "--calmode", "ap",
                        "--amplitude-limits", "0.6", "2.5", *extra])
        return self.json.loads((self.root / "round8_ph.json").read_text())

    def test_the_run_is_divided_out_before_the_limits_are_judged(self):
        report = self.run_solve("--detrend-radius", "0.1")

        self.assertTrue(report["radial_trend"]["divided_out"])
        self.assertGreater(report["radial_trend"]["ratio"], 1.3)
        written = [values for name, values in self.written if name == "CPARAM"]
        self.assertEqual(len(written), 1)
        size = np.abs(written[0][0, 0])
        self.assertLess(size.max() / size.min(), 1.05)
        self.assertAlmostEqual(float(np.median(size)), 1.0, places=6)
        self.assertEqual(len(report["amplitude_medians"]), 2)
        # The far antennas, near 0.55 before, are no longer out of limits.
        self.assertEqual(report["rejected_amplitudes"], [])
        self.assertIn("applycal", [name for name, _ in self.calls])

    def test_without_the_option_the_far_antennas_are_refused_as_before(self):
        report = self.run_solve()

        self.assertIsNone(report["radial_trend"])
        self.assertTrue(report["rejected_amplitudes"])
        # The run is kept; only the level is set, as CASA's solnorm would.
        written = [values for name, values in self.written if name == "CPARAM"]
        size = np.abs(written[0][0, 0])
        self.assertGreater(size.max() / size.min(), 2.0)
        self.assertAlmostEqual(float(np.median(size)), 1.0, places=6)

    def test_a_solve_casa_normalised_is_left_as_it_is(self):
        self.values /= np.median(np.abs(self.values[0, 0]))
        report = self.run_solve("--amplitude-limits", "0.25", "4")

        self.assertNotIn("CPARAM", [name for name, _ in self.written])
        self.assertEqual(report["amplitude_medians"], [1.0, 1.0])



class BandpassFlowTests(unittest.TestCase):
    """The bandpass round, end to end: its antennas' run is divided out."""

    def setUp(self):
        import json
        import sys
        import tempfile
        import types

        self.json = json
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.ms = self.root / "round14.ms"
        (self.ms / "ANTENNA").mkdir(parents=True)
        self.positions, self.names = gmrt_like_array()
        count, bins = len(self.names), 8
        centre = np.median(self.positions, axis=0)
        distance = np.linalg.norm(self.positions - centre, axis=1)
        # The far antennas high, as G0212-0335's bandpass round had them.
        level = np.array([
            0.95 + 0.12 * math.log10(max(d, 100.0) / 300.0) for d in distance
        ])
        self.values = np.ones((2, bins, count), dtype=complex) * level
        self.written = []
        test = self

        class Table:
            def open(self, path, nomodify=True):
                self.path = Path(path)

            def close(self):
                pass

            def colnames(self):
                return ["DATA", "MODEL_DATA"]

            def putcol(self, name, values):
                test.written.append((name, np.array(values)))

            def getcol(self, name):
                if self.path.name == "ANTENNA":
                    return {"NAME": np.array(test.names),
                            "POSITION": test.positions.T}[name]
                return {
                    "CPARAM": test.values.copy(),
                    "FLAG": np.zeros(test.values.shape, dtype=bool),
                    "ANTENNA1": np.arange(count),
                }[name]

        def bandpass(**kwargs):
            Path(kwargs["caltable"]).mkdir()

        def applycal(**kwargs):
            pass

        def split(**kwargs):
            Path(kwargs["outputvis"]).mkdir()

        casatasks = types.ModuleType("casatasks")
        casatasks.bandpass, casatasks.applycal, casatasks.split = (
            bandpass, applycal, split
        )
        casatasks.flagdata = lambda **kwargs: None
        casatasks.gaincal = lambda **kwargs: None
        casatools = types.ModuleType("casatools")
        casatools.table = Table
        patch = mock.patch.dict(
            sys.modules, {"casatasks": casatasks, "casatools": casatools}
        )
        patch.start()
        self.addCleanup(patch.stop)

    def run_round(self, *extra):
        with mock.patch("builtins.print"):
            SOLVE.main(["--ms", str(self.ms), "--bandpass", "60",
                        "--bandpass-amplitude", *extra])
        reports = sorted(self.root.glob("*.json"))
        return self.json.loads(reports[0].read_text())

    def test_the_run_is_divided_out_of_the_antennas_bandpasses(self):
        report = self.run_round("--detrend-radius", "0.05")

        self.assertTrue(report["radial_trend"]["divided_out"])
        self.assertLess(report["radial_trend"]["ratio"], 0.95)
        written = [values for name, values in self.written if name == "CPARAM"]
        size = np.abs(written[-1])
        self.assertLess(size.max() / size.min(), 1.03)
        self.assertAlmostEqual(float(np.median(size)), 1.0, places=6)

    def test_without_the_option_the_bandpasses_keep_the_run(self):
        report = self.run_round()

        self.assertIsNone(report["radial_trend"])
        written = [values for name, values in self.written if name == "CPARAM"]
        size = np.abs(written[-1])
        self.assertGreater(size.max() / size.min(), 1.1)


if __name__ == "__main__":
    unittest.main()
