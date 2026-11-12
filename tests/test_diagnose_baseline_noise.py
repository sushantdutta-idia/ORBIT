"""Noise by baseline length against what the weights say."""

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
    "test_diagnose_baseline_noise_tool", SOURCE_DIR / "diagnose_baseline_noise.py")
TOOL = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TOOL)


class FlowTests(unittest.TestCase):

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        self.ms = root / "round3.ms"
        self.ms.mkdir()
        self.out = root / "noise.json"
        rng = np.random.default_rng(3)
        rows, nchan = 400, 64
        # Half the baselines 0.3 klambda long, half 20; each hand with
        # weight 4 (sigma 0.5 per part), the short ones three times as
        # noisy as that, as interference would make them; a sky that is
        # the same in every channel; one autocorrelation in each ten rows.
        length_m = np.where(np.arange(rows) % 2 == 0, 300 * 0.24, 20000 * 0.24)
        self.uvw = np.vstack([length_m, np.zeros(rows), np.zeros(rows)])
        noise = np.where(np.arange(rows) % 2 == 0, 1.5, 0.5)
        sky = 0.2 + 0.1j
        data = np.empty((2, nchan, rows), dtype=complex)
        for corr in range(2):
            data[corr] = sky + noise * (rng.normal(size=(nchan, rows)) +
                                        1j * rng.normal(size=(nchan, rows)))
        self.columns = {
            "ANTENNA1": np.zeros(rows, dtype=int),
            "ANTENNA2": np.where(np.arange(rows) % 10 == 0, 0, 1),
            "UVW": self.uvw,
            "DATA": data,
            "FLAG": np.zeros((2, nchan, rows), dtype=bool),
            "WEIGHT": np.full((2, rows), 4.0),
        }
        freqs = np.full(nchan, 1.25e9)
        test = self

        class Table:
            def open(self, path):
                self.path = Path(path)

            def close(self):
                pass

            def colnames(self):
                return ["DATA", "FLAG", "WEIGHT", "UVW"]

            def iscelldefined(self, name, row):
                return False

            def nrows(self):
                return rows

            def getcol(self, name, start=0, count=-1):
                if self.path.name == "SPECTRAL_WINDOW":
                    return freqs.reshape(nchan, 1)
                values = test.columns[name]
                return values[..., start:start + count]

        self.casatools = types.SimpleNamespace(table=Table)

    def test_the_noisier_short_baselines_show_against_their_weights(self):
        with mock.patch.dict(sys.modules, {"casatools": self.casatools}), \
                contextlib.redirect_stdout(io.StringIO()):
            TOOL.main(["--ms", str(self.ms), "--out", str(self.out),
                       "--chunk", "97", "--edges", "0", "1", "100"])
        report = json.loads(self.out.read_text())
        short, long = report["bins"]
        for name in ("I", "V"):
            self.assertAlmostEqual(long[f"{name}_over_weights"], 1.0, delta=0.08)
            self.assertAlmostEqual(short[f"{name}_over_weights"], 3.0, delta=0.25)
            # The sky is gone from the channel differences.
            self.assertAlmostEqual(long[f"{name}_jy"], 0.5, delta=0.05)
        self.assertAlmostEqual(short["share_of_data"] + long["share_of_data"], 1.0)
        self.assertAlmostEqual(short["share_of_data"], 0.5, delta=0.06)


if __name__ == "__main__":
    unittest.main()
