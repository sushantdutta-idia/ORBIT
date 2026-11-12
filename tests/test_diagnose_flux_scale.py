"""The calibrators' flux density scan by scan, against Perley-Butler 2017."""

import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_diagnose_flux_scale_tool", SOURCE_DIR / "diagnose_flux_scale.py")
TOOL = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TOOL)


class ModelTests(unittest.TestCase):

    def test_3c48_as_setjy_gave_it(self):
        # setjy reported 15.80 Jy at 1450 MHz, the band's top edge.
        self.assertAlmostEqual(TOOL.model_flux("3C48", 1.45e9), 15.80, delta=0.02)
        self.assertAlmostEqual(TOOL.model_flux("3C48", 1.25e9), 17.80, delta=0.03)
        self.assertIsNone(TOOL.model_flux("0323+055", 1.25e9))


class FlowTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        self.ms = root / "Band5.ms"
        self.ms.mkdir()
        self.out = root / "flux.json"
        nchan = 8
        freqs = np.linspace(1.05e9, 1.45e9, nchan)
        rows = []
        # 3C48 in scan 1 at 0.8 of its model; 0323+055 in scans 2 and 6 at
        # 2.16 and 2.30 Jy; one long baseline at 10 Jy that the cut leaves
        # out; autocorrelations at 50 and 60.
        for scan, field, flux in ((1, 0, None), (2, 1, 2.16), (6, 1, 2.30)):
            for k in range(6):
                rows.append((scan, field, 0, k + 1, 300.0, flux))
            rows.append((scan, field, 0, 9, 20000.0, 10.0))
            rows.append((scan, field, 0, 0, 0.0, "auto"))
        self.rows = rows
        test = self

        def column(name, picked):
            out = []
            for scan, field, a1, a2, length, flux in picked:
                if name == "SCAN_NUMBER":
                    out.append(scan)
                elif name == "TIME":
                    out.append(1000.0 * scan)
                elif name == "UVW":
                    out.append([length, 0.0, 0.0])
                elif name in ("CORRECTED_DATA", "DATA"):
                    if flux == "auto":
                        level = 50.0 if field == 0 else 60.0
                        out.append(np.full((2, nchan), level, dtype=complex))
                    elif flux is None:
                        model = np.array([TOOL.model_flux("3C48", f) for f in freqs])
                        out.append(np.vstack([0.8 * model, 0.8 * model]).astype(complex))
                    else:
                        out.append(np.full((2, nchan), flux, dtype=complex))
                elif name == "FLAG":
                    out.append(np.zeros((2, nchan), dtype=bool))
            values = np.array(out)
            if name == "UVW":
                return values.T
            if values.ndim == 3:
                return np.moveaxis(values, 0, 2)
            return values

        class Selection:
            def __init__(self, picked):
                self.picked = picked

            def nrows(self):
                return len(self.picked)

            def getcol(self, name, start=0, count=-1):
                picked = self.picked[start:None if count < 0 else start + count]
                return column(name, picked)

            def close(self):
                pass

        class Table:
            def open(self, path):
                self.path = Path(path)

            def close(self):
                pass

            def colnames(self):
                return ["DATA", "CORRECTED_DATA", "FLAG", "UVW"]

            def getcol(self, name):
                if self.path.name == "FIELD":
                    return np.array(["3C48", "0323+055"])
                if self.path.name == "SPECTRAL_WINDOW":
                    return freqs.reshape(nchan, 1)
                raise KeyError(name)

            def query(self, condition):
                picked = []
                for row in test.rows:
                    scan, field, a1, a2, _, _ = row
                    env = {"FIELD_ID": field, "SCAN_NUMBER": scan,
                           "ANTENNA1": a1, "ANTENNA2": a2}
                    expression = condition.replace("&&", " and ").replace(
                        "==", " == ").replace("!=", " != ")
                    if eval(expression, {}, env):
                        picked.append(row)
                return Selection(picked)

        self.casatools = types.SimpleNamespace(table=Table)

    def run_tool(self, *extra):
        with mock.patch.dict(sys.modules, {"casatools": self.casatools}), \
                contextlib.redirect_stdout(io.StringIO()) as printed:
            TOOL.main(["--ms", str(self.ms), "--out", str(self.out),
                       "--chunk", "3", *extra])
        return json.loads(self.out.read_text()), printed.getvalue()

    def test_each_scan_against_the_model_and_the_autocorrelations(self):
        report, printed = self.run_tool()

        flux_cal = report["fields"]["3C48"][0]
        self.assertEqual(flux_cal["scan"], 1)
        for band in flux_cal["sub_bands"]:
            self.assertAlmostEqual(band["ratio"], 0.8, places=3)
        self.assertEqual(flux_cal["autocorrelation_median"], 50.0)
        phase = report["fields"]["0323+055"]
        self.assertEqual([entry["scan"] for entry in phase], [2, 6])
        # The 10 Jy long baseline is outside the cut.
        self.assertAlmostEqual(phase[0]["flux_jy"], 2.16, places=3)
        self.assertAlmostEqual(phase[1]["flux_jy"], 2.30, places=3)
        self.assertIsNone(phase[0]["sub_bands"][0]["ratio"])
        self.assertEqual(phase[1]["autocorrelation_median"], 60.0)
        self.assertIn("0.800 of the model", printed)

    def test_a_longer_cut_takes_the_long_baseline_in(self):
        report, _ = self.run_tool("--max-klambda", "200")
        self.assertGreater(report["fields"]["0323+055"][0]["flux_jy"], 3.0)

    def test_flux_against_baseline_length(self):
        report, printed = self.run_tool(
            "--uv-bins", "0", "5", "200", "--model", "0323+055=2.2@1250")
        bins = report["uv_profiles"]["0323+055"]["bins"]
        self.assertEqual([entry["klambda"] for entry in bins],
                         [[0.0, 5.0], [5.0, 200.0]])
        # Both scans' short baselines, 2.16 and 2.30 Jy, against 2.2 Jy.
        self.assertAlmostEqual(bins[0]["flux_jy"], 2.23, places=3)
        for band in bins[0]["sub_bands"]:
            self.assertAlmostEqual(band["ratio"], 2.23 / 2.2, places=3)
        self.assertAlmostEqual(bins[1]["flux_jy"], 10.0, places=3)
        # 3C48 keeps its standard model.
        short = report["uv_profiles"]["3C48"]["bins"][0]["sub_bands"]
        self.assertAlmostEqual(short[0]["ratio"], 0.8, places=3)
        self.assertIn("against baseline length", printed)
        # The per-scan report now uses the stated model too.
        self.assertAlmostEqual(
            report["fields"]["0323+055"][0]["sub_bands"][0]["ratio"],
            2.16 / 2.2, places=3)

    def test_a_scan_can_be_left_out_of_the_profile(self):
        report, _ = self.run_tool("--uv-bins", "0", "5", "200",
                                  "--profile-skip-scans", "6")
        profile = report["uv_profiles"]["0323+055"]
        self.assertEqual(profile["skipped_scans"], [6])
        self.assertAlmostEqual(profile["bins"][0]["flux_jy"], 2.16, places=3)
        # The scan stays in the per-scan report.
        self.assertEqual(
            [entry["scan"] for entry in report["fields"]["0323+055"]], [2, 6])

    def test_an_uncalibrated_ms_serves_only_the_time_profile(self):
        with mock.patch.object(TOOL, "time_profiles",
                               return_value={"1": {
                                   "field": "3C48", "step_start_s": [0.0],
                                   "array_median": [1.0], "antennas": {}}}) \
                as profiles:
            self.casatools.table.colnames = lambda self: ["DATA", "FLAG",
                                                          "UVW"]
            report, printed = self.run_tool("--time-profile", "1")
        self.assertEqual(report["fields"], {})
        self.assertEqual(report["column"], None)
        self.assertIn("1", report["time_profiles"])
        self.assertIn("only the time profile", printed)
        profiles.assert_called_once()

    def test_an_uncalibrated_ms_without_a_time_profile_is_refused(self):
        self.casatools.table.colnames = lambda self: ["DATA", "FLAG", "UVW"]
        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()
        self.assertIn("not calibrated", str(stopped.exception))

    def test_a_model_that_cannot_be_read_is_refused(self):
        with self.assertRaises(SystemExit):
            self.run_tool("--model", "0323+055=3.08")

    def test_an_unknown_field_is_said_so(self):
        with self.assertRaises(SystemExit) as stopped:
            self.run_tool("--fields", "3C286")
        self.assertIn("3C286", str(stopped.exception))


class StatedModelTests(unittest.TestCase):

    def test_a_power_law_about_its_reference(self):
        stated = TOOL.parse_models(["0323+055=3.08@1250,-0.76"])
        self.assertAlmostEqual(TOOL.model_flux("0323+055", 1.25e9, stated),
                               3.08)
        self.assertAlmostEqual(
            TOOL.model_flux("0323+055", 1.3675e9, stated),
            3.08 * (1367.5 / 1250) ** -0.76, places=6)


class TimeProfileTests(unittest.TestCase):

    def test_an_antenna_that_settles_late(self):
        # Three antennas, 20 one-minute steps of 10 s rows on each baseline.
        # Antenna 2 is at 20% of its level for the first five minutes.
        times, first, second, levels = [], [], [], []
        for t in np.arange(0.0, 1200.0, 10.0):
            for a1, a2 in ((0, 1), (0, 2), (1, 2)):
                low = t < 300 and 2 in (a1, a2)
                times.append(t)
                first.append(a1)
                second.append(a2)
                levels.append(0.2 if low else 1.0)
        steps, antennas = TOOL.time_profile(
            times, first, second, levels, ["C00", "C01", "S04"], 60.0)
        self.assertEqual(len(steps), 20)
        self.assertEqual(steps[5], 300.0)
        self.assertEqual(antennas["S04"]["settled_step"], 5)
        self.assertAlmostEqual(antennas["S04"]["ratios"][0], 0.2)
        # The others share a baseline with it, but the median over their
        # two baselines is not dragged below 0.6.
        self.assertEqual(antennas["C00"]["settled_step"], 5)
        self.assertAlmostEqual(antennas["C00"]["ratios"][0], 0.6)

    def test_an_antenna_that_never_settles(self):
        times = np.repeat(np.arange(0.0, 600.0, 10.0), 1)
        levels = np.where(np.arange(times.size) % 12 < 6, 0.5, 1.5)
        steps, antennas = TOOL.time_profile(
            times, np.zeros(times.size, dtype=int),
            np.ones(times.size, dtype=int), levels, ["A", "B"], 60.0)
        self.assertIsNone(antennas["A"]["settled_step"])


if __name__ == "__main__":
    unittest.main()
