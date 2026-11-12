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
    "test_diagnose_delay_drift_tool", SOURCE_DIR / "diagnose_delay_drift.py"
)
TOOL = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TOOL)


class DriftTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        self.ms = root / "Band5.ms"
        self.ms.mkdir()
        self.caldir = root / "calibration"
        for name in ("delay.cal", "bandpass.cal"):
            (self.caldir / name).mkdir(parents=True)
        (self.caldir / "delay_quarantined_antennas.json").write_text(json.dumps({
            "antennas": [{"id": 3, "name": "C03"}, {"id": 29, "name": "W06"}],
        }))
        self.out = root / "drift" / "Band5.K.cal"
        self.calls = []
        # Antenna 1 drifts in one hand; antenna 2 stays put; one solution fails.
        antennas = [0, 1, 2, 0, 1, 2, 0, 1, 2]
        scans = [2, 2, 2, 10, 10, 10, 23, 23, 23]
        fields = [1, 1, 1, 1, 1, 1, 1, 1, 1]
        r = [0.0, 0.2, 0.1, 0.0, 2.6, 0.1, 0.0, 5.1, 0.2]
        l = [0.0, 0.1, 0.0, 0.0, 0.2, 0.1, 0.0, 0.3, 0.1]
        fparam = np.array([[r], [l]], dtype=float).reshape(2, 1, 9)
        flag = np.zeros((2, 1, 9), dtype=bool)
        flag[1, 0, 8] = True
        self.tables = {
            "ANTENNA": {"NAME": np.array(["C00", "C01", "C02"])},
            "FIELD": {"NAME": np.array(["3C48", "0323+055"])},
            "K": {"ANTENNA1": np.array(antennas), "SCAN_NUMBER": np.array(scans),
                  "FIELD_ID": np.array(fields), "FPARAM": fparam, "FLAG": flag},
        }
        test = self

        class Table:
            def open(self, path):
                path = Path(path)
                self.columns = test.tables["K" if path.suffix == ".cal" else path.name]

            def getcol(self, name):
                return self.columns[name]

            def close(self):
                pass

        def gaincal(**kwargs):
            test.calls.append(kwargs)
            Path(kwargs["caltable"]).mkdir(parents=True)

        self.modules = {
            "casatasks": types.SimpleNamespace(gaincal=gaincal),
            "casatools": types.SimpleNamespace(table=Table),
        }

    def run_tool(self, *extra):
        with mock.patch.dict(sys.modules, self.modules), \
                contextlib.redirect_stdout(io.StringIO()) as output:
            TOOL.main([
                "--ms", str(self.ms), "--caldir", str(self.caldir),
                "--fields", "0323+055,3C48", "--out", str(self.out), *extra,
            ])
        return output.getvalue()

    def test_a_delay_per_scan_with_the_pipelines_tables_applied(self):
        printed = self.run_tool()

        call = self.calls[0]
        self.assertEqual(call["gaintype"], "K")
        self.assertEqual(call["solint"], "inf")
        self.assertEqual(call["combine"], "")
        self.assertEqual(call["antenna"], "!C03,W06")
        self.assertEqual(call["gaintable"], [
            str(self.caldir.resolve() / "delay.cal"),
            str(self.caldir.resolve() / "bandpass.cal"),
        ])
        self.assertEqual(call["gainfield"], ["3C48", "3C48"])
        self.assertEqual(call["vis"], str(self.ms.resolve()))

        report = json.loads(self.out.with_suffix(".json").read_text())
        self.assertEqual(report["scans"], [2, 10, 23])
        by_name = {entry["antenna"]: entry for entry in report["antennas"]}
        self.assertEqual(by_name["C01"]["spread_ns"], [4.9, 0.2])
        self.assertEqual(by_name["C02"]["by_scan"]["23"]["delay_ns"], [0.2, None])
        self.assertEqual(by_name["C02"]["spread_ns"], [0.1, 0.1])
        self.assertIn("C01", printed)
        self.assertIn("+5.1/ +0.3", printed)
        self.assertIn("4.9/0.2", printed)

    def test_antennas_to_leave_out_can_be_named(self):
        self.run_tool("--exclude")
        self.assertEqual(self.calls[0]["antenna"], "")

    def test_a_missing_table_stops_it_before_solving(self):
        (self.caldir / "bandpass.cal").rmdir()
        with self.assertRaises(SystemExit):
            self.run_tool()
        self.assertEqual(self.calls, [])


if __name__ == "__main__":
    unittest.main()
