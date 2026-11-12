import importlib.util
import json
import math
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load(filename, name):
    spec = importlib.util.spec_from_file_location(name, SOURCE_DIR / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TOOL = load("select_scans.py", "test_select_scans_module")
RECIPE = load("selfcal_recipe.py", "test_select_scans_recipe")


def option(command, name, count=1):
    index = command.index(name)
    values = command[index + 1:index + 1 + count]
    return values[0] if count == 1 else values


def scan(number, amount, sigma, weight=1.0, minutes=None, robust=None):
    """A scan of ``amount`` visibilities of scatter ``sigma`` per component."""

    return {
        "scan": number, "rows": int(amount), "integrations": 10,
        "minutes": amount if minutes is None else minutes,
        "antennas": 27, "flagged": 0.4, "samples": int(amount),
        "sum_w": weight * amount, "sum_w2": weight ** 2 * amount,
        "sum_w2v2": weight ** 2 * amount * 2 * sigma ** 2,
        "sum_v2": amount * 2 * sigma ** 2,
        "robust_jy": sigma if robust is None else robust,
    }


def g0212():
    """G0212-0335 as its images gave it: 61.0, 23.9 and 110.9 uJy/beam."""

    return [scan(3, 9300, 1.65e-2), scan(15, 22200, 1.0e-2),
            scan(29, 1300, 1.12e-2)]


class StartTests(unittest.TestCase):

    def test_it_starts_beside_the_pipelines_own_inspect_module(self):
        source = (SOURCE_DIR / "select_scans.py").read_text()
        self.assertIn("sys.path.pop(0)", source.split("import argparse")[0])
        done = subprocess.run(
            [sys.executable, "select_scans.py", "--help"],
            cwd=SOURCE_DIR, capture_output=True, text=True,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("--least-share", done.stdout)
        self.assertNotIn("--casa-part", done.stdout)

    def test_a_target_and_a_work_folder_are_asked_for(self):
        done = subprocess.run(
            [sys.executable, "select_scans.py", "--ms", "x.ms",
             "--channel-range", "51", "768"],
            cwd=SOURCE_DIR, capture_output=True, text=True,
        )
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("--target and --workdir are required", done.stderr)


class HandsTests(unittest.TestCase):

    def test_the_parallel_hands_are_found_among_the_correlations(self):
        self.assertEqual(TOOL.parallel_hands([5, 8]), (0, 1))
        self.assertEqual(TOOL.parallel_hands([5, 6, 7, 8]), (0, 3))
        self.assertEqual(TOOL.parallel_hands([9, 10, 11, 12]), (0, 3))

    def test_data_without_both_hands_are_refused(self):
        with self.assertRaises(SystemExit) as stopped:
            TOOL.parallel_hands([5])
        self.assertIn("do not hold both parallel hands",
                      str(stopped.exception))

    def test_corrected_data_are_read_if_there_are_any(self):
        self.assertEqual(TOOL.data_column(["DATA", "FLAG"]), "DATA")
        self.assertEqual(
            TOOL.data_column(["DATA", "CORRECTED_DATA"]), "CORRECTED_DATA")


class DerivedTests(unittest.TestCase):

    def test_scatter_noise_and_share_follow_from_the_sums(self):
        scans = g0212()

        main = TOOL.derived(scans)

        self.assertEqual(main, 15)
        first, whole, last = scans
        self.assertAlmostEqual(whole["rms_jy"], 1.0e-2)
        self.assertAlmostEqual(whole["noise_jy"], 1.0e-2 / math.sqrt(22200))
        self.assertAlmostEqual(first["excess"], 1.65)
        self.assertAlmostEqual(whole["excess"], 1.0)
        # The shares the images of the three scans gave: 13, 86 and 4%.
        self.assertAlmostEqual(first["share"], 0.128, places=3)
        self.assertAlmostEqual(whole["share"], 0.833, places=3)
        self.assertAlmostEqual(last["share"], 0.039, places=3)
        self.assertAlmostEqual(sum(s["share"] for s in scans), 1.0)
        self.assertAlmostEqual(whole["outliers"], 1.0)

    def test_the_weights_the_data_carry_do_not_change_a_scans_own_noise(self):
        plain, weighted = [scan(1, 1000, 2e-2)], [scan(1, 1000, 2e-2, 7.0)]
        TOOL.derived(plain)
        TOOL.derived(weighted)

        self.assertAlmostEqual(plain[0]["noise_jy"], weighted[0]["noise_jy"])
        self.assertAlmostEqual(weighted[0]["rms_jy"], 2e-2)

    def test_interference_left_in_shows_as_scatter_beyond_the_bulk(self):
        scans = [scan(1, 1000, 3e-2, robust=1e-2)]
        TOOL.derived(scans)

        self.assertAlmostEqual(scans[0]["outliers"], 3.0)

    def test_a_scan_without_data_has_no_noise_and_no_share(self):
        empty = scan(7, 0, 1.0)
        empty["robust_jy"] = None
        scans = [scan(5, 1000, 1e-2), empty]

        self.assertEqual(TOOL.derived(scans), 5)

        self.assertIsNone(scans[1]["noise_jy"])
        self.assertIsNone(scans[1]["excess"])
        self.assertEqual(scans[1]["share"], 0.0)
        self.assertAlmostEqual(scans[0]["share"], 1.0)

    def test_the_floor_is_that_of_the_scans_as_an_imager_adds_them(self):
        scans = g0212()
        # Equal weights: each scan counts by its amount of data, not by
        # its noise, so the noisy scan costs more than it adds.
        everything = TOOL.floor(scans)
        whole = TOOL.floor(scans[1:2])

        expected = math.sqrt(
            9300 * 1.65e-2 ** 2 + 22200 * 1.0e-2 ** 2 + 1300 * 1.12e-2 ** 2
        ) / (9300 + 22200 + 1300)
        self.assertAlmostEqual(everything, expected)
        self.assertGreater(everything, whole)
        self.assertIsNone(TOOL.floor([]))


class ChooseTests(unittest.TestCase):

    def choose(self, scans, forced=None, least=0.10, most=1.4):
        TOOL.derived(scans)
        return TOOL.choose(scans, least, most, forced)

    def test_g0212_keeps_its_one_whole_scan(self):
        scans = g0212()

        self.assertEqual(self.choose(scans), [15])

        self.assertIn("1.65 times as noisy as scan 15", scans[0]["why"])
        self.assertIn("4% of the weight", scans[2]["why"])
        self.assertIsNone(scans[1]["why"])

    def test_two_whole_scans_alike_are_both_kept(self):
        scans = [scan(5, 22000, 1.0e-2), scan(18, 22000, 1.15e-2)]

        self.assertEqual(self.choose(scans), [5, 18])

    def test_a_noisy_scan_is_left_out_whatever_its_length(self):
        scans = [scan(5, 22000, 1.5e-2), scan(18, 22000, 1.0e-2)]

        self.assertEqual(self.choose(scans), [18])
        self.assertIn("1.50 times as noisy as scan 18", scans[0]["why"])

    def test_the_scan_with_the_most_weight_is_always_kept(self):
        # A short quiet scan does not make the long one the noisy one.
        scans = [scan(5, 22000, 1.3e-2), scan(29, 1300, 1.0e-2)]

        self.assertEqual(self.choose(scans), [5])
        self.assertIn("too little for solutions of its own", scans[1]["why"])
        self.assertAlmostEqual(scans[1]["excess"], 1.0 / 1.3)
        # Nor does a limit no scan can meet leave a target without data.
        alike = [scan(5, 22000, 1.0e-2), scan(18, 21000, 1.0e-2)]
        self.assertEqual(self.choose(alike, least=0.9), [5])

    def test_a_single_scan_is_kept(self):
        self.assertEqual(self.choose([scan(17, 22000, 1e-2)]), [17])

    def test_the_limits_can_be_moved(self):
        self.assertEqual(self.choose(g0212(), most=2.0), [3, 15])
        self.assertEqual(self.choose(g0212(), least=0.01, most=2.0),
                         [3, 15, 29])

    def test_named_scans_overrule_the_rules_which_are_still_recorded(self):
        scans = g0212()

        self.assertEqual(self.choose(scans, forced=[3, 15]), [3, 15])

        self.assertTrue(scans[0]["kept"])
        self.assertIn("times as noisy", scans[0]["rule"])
        self.assertFalse(scans[2]["kept"])
        self.assertEqual(scans[2]["why"], "it was not named with --scans")

    def test_a_named_scan_that_is_not_there_is_said_with_those_there_are(self):
        with self.assertRaises(SystemExit) as stopped:
            self.choose(g0212(), forced=[16])

        self.assertIn("No scan 16", str(stopped.exception))
        self.assertIn("there are: 3, 15, 29", str(stopped.exception))

    def test_a_scan_without_data_is_left_out_and_cannot_be_named(self):
        empty = scan(7, 0, 1.0)
        empty["robust_jy"] = None
        scans = [scan(5, 1000, 1e-2), empty]

        self.assertEqual(self.choose(scans), [5])
        self.assertEqual(scans[1]["why"], "it holds no data")
        with self.assertRaises(SystemExit):
            self.choose([scan(5, 1000, 1e-2), dict(empty)], forced=[7])


class ReadingTests(unittest.TestCase):

    def test_each_scan_has_a_line_and_the_cost_is_said(self):
        scans = g0212()
        TOOL.derived(scans)
        TOOL.choose(scans, 0.10, 1.4)

        lines = TOOL.reading("G0212-0335", scans)

        self.assertEqual(len(lines), 4)
        self.assertIn("scan 3: 9300.0 min, 27 antennas, 40% flagged", lines[0])
        self.assertIn("13% of the weight, 1.65 times the noise of the main "
                      "scan", lines[0])
        self.assertIn("left out: it is 1.65 times as noisy", lines[0])
        self.assertTrue(lines[1].endswith("; kept."))
        self.assertIn("G0212-0335: kept scan 15; the data then allow",
                      lines[3])
        self.assertIn("with every scan", lines[3])

    def test_every_scan_kept_is_said_so(self):
        scans = [scan(5, 22000, 1.0e-2), scan(18, 22000, 1.1e-2)]
        TOOL.derived(scans)
        TOOL.choose(scans, 0.10, 1.4)

        lines = TOOL.reading("G0216-0604", scans)

        self.assertIn("G0216-0604: every scan is kept", lines[-1])

    def test_interference_and_named_scans_are_said(self):
        scans = [scan(3, 9300, 1.65e-2, robust=0.825e-2),
                 scan(15, 22200, 1e-2)]
        TOOL.derived(scans)
        TOOL.choose(scans, 0.10, 1.4, [3])

        lines = TOOL.reading("G", scans, forced=True)

        self.assertIn("(scatter 2.0 times what its bulk has, against 1.0 in "
                      "the main scan: interference left in)", lines[0])
        self.assertIn("kept, named with --scans though it is 1.65 times",
                      lines[0])
        self.assertIn("not named with --scans, though the rules would keep "
                      "it", lines[1])

    def test_the_scatter_usual_in_whole_scans_is_not_called_interference(self):
        # As measured on the nine fields: 1.2 to 1.6 in every whole scan,
        # 3.35 in the piece of G0212-0335's first.
        scans = [scan(9, 22000, 1.57e-2, robust=1.0e-2),
                 scan(23, 22200, 1.34e-2, robust=1.0e-2),
                 scan(3, 9300, 3.35e-2, robust=1.0e-2)]
        TOOL.derived(scans)
        TOOL.choose(scans, 0.01, 9.0)

        lines = TOOL.reading("G", scans)

        self.assertNotIn("interference", lines[0])
        self.assertNotIn("interference", lines[1])
        self.assertIn("scatter 3.4 times what its bulk has, against 1.3 in "
                      "the main scan", lines[2])

    def test_what_the_rules_say_of_a_scan_named_out_is_still_said(self):
        scans = g0212()
        TOOL.derived(scans)
        TOOL.choose(scans, 0.10, 1.4, [15])

        lines = TOOL.reading("G0212-0335", scans, forced=True)

        self.assertIn("left out: it was not named with --scans, and it is "
                      "1.65 times as noisy as scan 15", lines[0])
        self.assertIn("left out: it was not named with --scans, and it "
                      "holds 4% of the weight", lines[2])
        self.assertTrue(lines[1].endswith("; kept."))


# ---------------------------------------------------------------------------
# The part that runs under CASA, on tables made of arrays
# ---------------------------------------------------------------------------

class FakeTable:
    """The few calls of casatools.table that the measurement makes."""

    tables = {}
    opened = []

    def __init__(self, columns=None):
        self.columns = columns

    def open(self, path, *args, **kwargs):
        FakeTable.opened.append((path, args, kwargs))
        self.columns = FakeTable.tables[path]

    def close(self):
        pass

    def colnames(self):
        return list(self.columns)

    def nrows(self):
        return self.columns["FLAG"].shape[-1] if "FLAG" in self.columns else 1

    def iscelldefined(self, name, row):
        return True

    def getcell(self, name, row):
        return self.columns[name][..., row]

    def getcol(self, name, startrow=0, nrow=-1):
        end = None if nrow < 0 else startrow + nrow
        return self.columns[name][..., startrow:end]

    def getcolslice(self, name, blc, trc, incr, startrow=0, nrow=-1):
        end = None if nrow < 0 else startrow + nrow
        return self.columns[name][
            blc[0]:trc[0] + 1, blc[1]:trc[1] + 1, startrow:end]

    def query(self, text):
        number = int(text.split("==")[1].split()[0])
        assert "ANTENNA1!=ANTENNA2" in text, text
        rows = (self.columns["SCAN_NUMBER"] == number) & (
            self.columns["ANTENNA1"] != self.columns["ANTENNA2"])
        return FakeTable(
            {name: value[..., rows] for name, value in self.columns.items()})


def fake_ms(spectrum=True, seed=3):
    """Two scans, four antennas with their autocorrelations, twelve channels.

    Both hands hold the same bright sky and noise of their own: 0.02 per
    component in scan 3 and 0.01 in scan 15.
    """

    rng = np.random.default_rng(seed)
    pairs = [(a, b) for a in range(4) for b in range(a, 4)]
    scans, times, first, second = [], [], [], []
    for number, steps, start in ((3, 6, 100.0), (15, 14, 5000.0)):
        for step in range(steps):
            for a, b in pairs:
                scans.append(number)
                times.append(start + 10.0 * step)
                first.append(a)
                second.append(b)
    rows = len(scans)
    scans = np.array(scans)
    sigma = np.where(scans == 3, 0.02, 0.01)
    sky = 5.0 * np.exp(1j * rng.uniform(0, 6.28, (12, rows)))
    noise = (rng.normal(size=(2, 12, rows))
             + 1j * rng.normal(size=(2, 12, rows))) * sigma
    data = (sky[None] + noise).astype(np.complex64)
    flag = rng.uniform(size=(2, 12, rows)) < 0.2
    row_flag = rng.uniform(size=rows) < 0.1
    weight = rng.uniform(0.5, 2.0, size=(2, 12, rows)).astype(np.float32)
    columns = {
        "DATA": data, "FLAG": flag, "FLAG_ROW": row_flag,
        "SCAN_NUMBER": scans, "TIME": np.array(times),
        "ANTENNA1": np.array(first), "ANTENNA2": np.array(second),
        "INTERVAL": np.full(rows, 10.0),
        "DATA_DESC_ID": np.zeros(rows, int),
        "WEIGHT": weight[:, 0, :].copy(),
    }
    if spectrum:
        columns["WEIGHT_SPECTRUM"] = weight
    return columns


def expected(columns, number, first, end, spectrum=True):
    """The sums of a scan, worked out the long way."""

    rows = (columns["SCAN_NUMBER"] == number) & (
        columns["ANTENNA1"] != columns["ANTENNA2"])
    data = columns["DATA"][:, first:end, rows]
    flag = columns["FLAG"][:, first:end, rows] | columns["FLAG_ROW"][rows]
    if spectrum:
        weight = columns["WEIGHT_SPECTRUM"][:, first:end, rows]
    else:
        weight = np.repeat(
            columns["WEIGHT"][:, None, rows], end - first, axis=1)
    use = ~flag[0] & ~flag[1]
    mean = (0.5 * (weight[0] + weight[1]).astype(float))[use]
    difference = (0.5 * (data[0] - data[1]))[use]
    power = np.abs(difference).astype(float) ** 2
    return {
        "samples": int(use.sum()),
        "sum_w": float(mean.sum()),
        "sum_w2": float((mean ** 2).sum()),
        "sum_w2v2": float((mean ** 2 * power).sum()),
        "sum_v2": float(power.sum()),
        "flagged": float(flag.mean()),
        "integrations": len(set(
            columns["TIME"][rows][use.any(axis=0)].tolist())),
    }


class MeasureTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.out = Path(self.temporary.name) / "measured.json"
        FakeTable.opened = []

    def measure(self, columns, channel_range=(2, 10), corr_types=(5, 8),
                chunk=None):
        FakeTable.tables = {
            "t.ms": columns,
            "t.ms/POLARIZATION": {
                "CORR_TYPE": np.array(corr_types)[:, None]},
        }
        module = types.SimpleNamespace(table=FakeTable)
        args = types.SimpleNamespace(
            ms="t.ms", channel_range=list(channel_range), out=str(self.out))
        patches = [mock.patch.dict(sys.modules, {"casatools": module}),
                   mock.patch("builtins.print")]
        if chunk:
            patches.append(mock.patch.object(TOOL, "CHUNK_SAMPLES", chunk))
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)
        TOOL.casa_measure(args)
        return json.loads(self.out.read_text())

    def test_each_scan_is_summed_over_its_unflagged_cross_correlations(self):
        columns = fake_ms()

        found = self.measure(columns)

        self.assertEqual(found["column"], "DATA")
        self.assertEqual(found["weights"], "WEIGHT_SPECTRUM")
        self.assertEqual([s["scan"] for s in found["scans"]], [3, 15])
        for measured in found["scans"]:
            want = expected(columns, measured["scan"], 2, 10)
            for key, value in want.items():
                self.assertAlmostEqual(
                    measured[key], value, delta=1e-6 * max(1.0, abs(value)),
                    msg=key)
            self.assertEqual(measured["antennas"], 4)
        first, whole = found["scans"]
        # Six baselines of four antennas: the autocorrelations are not read.
        self.assertEqual(first["rows"], 6 * 6)
        self.assertEqual(whole["rows"], 14 * 6)
        self.assertAlmostEqual(whole["minutes"], 14 * 10.0 / 60.0)

    def test_the_sky_cancels_and_the_noise_of_each_scan_is_what_was_put_in(self):
        found = self.measure(fake_ms())

        scans = found["scans"]
        TOOL.derived(scans)
        # Each hand has the scatter put in; half their difference has that
        # over the square root of two.
        self.assertAlmostEqual(
            scans[0]["rms_jy"], 0.02 / math.sqrt(2), delta=0.0015)
        self.assertAlmostEqual(
            scans[1]["rms_jy"], 0.01 / math.sqrt(2), delta=0.0007)
        self.assertAlmostEqual(scans[0]["excess"], 2.0, delta=0.2)
        self.assertAlmostEqual(
            scans[1]["robust_jy"], 0.01 / math.sqrt(2), delta=0.001)

    def test_reading_in_pieces_gives_the_same_sums(self):
        columns = fake_ms()
        whole = self.measure(columns)

        pieces = self.measure(columns, chunk=8 * 5)

        for one, other in zip(whole["scans"], pieces["scans"]):
            for key in ("samples", "sum_w", "sum_w2", "sum_w2v2", "sum_v2",
                        "flagged", "integrations", "antennas"):
                self.assertAlmostEqual(one[key], other[key], places=6)

    def test_without_weights_per_channel_the_rows_weights_are_used(self):
        columns = fake_ms(spectrum=False)

        found = self.measure(columns)

        self.assertEqual(found["weights"], "WEIGHT")
        want = expected(columns, 15, 2, 10, spectrum=False)
        self.assertAlmostEqual(found["scans"][1]["sum_w"], want["sum_w"],
                               places=3)
        self.assertAlmostEqual(found["scans"][1]["sum_w2v2"],
                               want["sum_w2v2"], places=6)

    def test_the_hands_are_taken_from_among_four_correlations(self):
        two = fake_ms()
        four = dict(two)
        zeros = np.zeros_like(two["DATA"][:1])
        four["DATA"] = np.concatenate(
            [two["DATA"][:1], zeros, zeros, two["DATA"][1:]])
        for name in ("FLAG", "WEIGHT_SPECTRUM"):
            four[name] = np.concatenate(
                [two[name][:1], two[name][:1], two[name][1:], two[name][1:]])
        four["WEIGHT"] = four["WEIGHT_SPECTRUM"][:, 0, :].copy()

        found = self.measure(four, corr_types=(5, 6, 7, 8))

        want = expected(two, 15, 2, 10)
        self.assertEqual(found["scans"][1]["samples"], want["samples"])
        self.assertAlmostEqual(found["scans"][1]["sum_w2v2"],
                               want["sum_w2v2"], places=6)

    def test_a_scan_flagged_throughout_is_reported_empty(self):
        columns = fake_ms()
        columns["FLAG"] = columns["FLAG"] | (columns["SCAN_NUMBER"] == 3)

        found = self.measure(columns)

        first = found["scans"][0]
        self.assertEqual(first["samples"], 0)
        self.assertEqual(first["integrations"], 0)
        self.assertIsNone(first["robust_jy"])
        self.assertEqual(first["flagged"], 1.0)

    def test_the_data_are_opened_only_to_be_read(self):
        self.measure(fake_ms())

        # casatools opens a table for reading unless told otherwise.
        for path, args, kwargs in FakeTable.opened:
            self.assertEqual(args, ())
            self.assertEqual(kwargs, {})

    def test_a_channel_range_beyond_the_data_is_refused(self):
        with self.assertRaises(SystemExit) as stopped:
            self.measure(fake_ms(), channel_range=(2, 13))

        self.assertIn("does not fit 12 channels", str(stopped.exception))

    def test_more_than_one_spectral_window_is_refused(self):
        columns = fake_ms()
        columns["DATA_DESC_ID"][::2] = 1

        with self.assertRaises(SystemExit) as stopped:
            self.measure(columns)

        self.assertIn("more than one spectral window", str(stopped.exception))


class SplitTests(unittest.TestCase):

    def test_the_kept_scans_are_split_from_the_column_that_was_read(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder) / "T.ms"
            out.mkdir()
            (out / "old").write_text("")
            calls = []

            def split(**kwargs):
                calls.append(kwargs)
                Path(kwargs["outputvis"]).mkdir()

            FakeTable.tables = {"t.ms": {"DATA": None, "FLAG": None}}
            modules = {
                "casatasks": types.SimpleNamespace(split=split),
                "casatools": types.SimpleNamespace(table=FakeTable),
            }
            args = types.SimpleNamespace(ms="t.ms", out=str(out),
                                         scans=[15, 29])
            with mock.patch.dict(sys.modules, modules), \
                    mock.patch("builtins.print"):
                TOOL.casa_split(args)

            self.assertEqual(calls, [{
                "vis": "t.ms", "outputvis": str(out), "scan": "15,29",
                "datacolumn": "data",
            }])
            # What an earlier, interrupted split left is not built upon.
            self.assertFalse((out / "old").exists())

    def test_channels_are_averaged_in_the_split_when_asked(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder) / "T.ms"
            calls = []

            def split(**kwargs):
                calls.append(kwargs)
                Path(kwargs["outputvis"]).mkdir()

            FakeTable.tables = {"t.ms": {"DATA": None}}
            modules = {
                "casatasks": types.SimpleNamespace(split=split),
                "casatools": types.SimpleNamespace(table=FakeTable),
            }
            args = types.SimpleNamespace(ms="t.ms", out=str(out), scans=[5, 18],
                                         average_channels=4)
            with mock.patch.dict(sys.modules, modules), \
                    mock.patch("builtins.print"):
                TOOL.casa_split(args)

            self.assertEqual(calls[0]["width"], 4)
            self.assertEqual(calls[0]["scan"], "5,18")

    def test_corrected_data_are_split_if_those_were_read(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder) / "T.ms"
            calls = []

            def split(**kwargs):
                calls.append(kwargs)
                Path(kwargs["outputvis"]).mkdir()

            FakeTable.tables = {"t.ms": {"DATA": None,
                                         "CORRECTED_DATA": None}}
            modules = {
                "casatasks": types.SimpleNamespace(split=split),
                "casatools": types.SimpleNamespace(table=FakeTable),
            }
            args = types.SimpleNamespace(ms="t.ms", out=str(out), scans=[15])
            with mock.patch.dict(sys.modules, modules), \
                    mock.patch("builtins.print"):
                TOOL.casa_split(args)

            self.assertEqual(calls[0]["datacolumn"], "corrected")


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

class RunTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.products = self.root / "uGMRT_Products"
        self.ms = self.products / "calibration" / "calibrated" / "T.ms"
        self.ms.mkdir(parents=True)
        self.work = self.root / "run" / "T" / "scans"
        self.measured = g0212()
        self.commands = []
        self.asked = []
        test = self

        def run_command(command, log_path):
            test.assertTrue(Path(log_path).parent.is_dir())
            test.commands.append(command)
            part = option(command, "--casa-part")
            out = Path(option(command, "--out"))
            if part == "measure":
                out.write_text(json.dumps({
                    "column": "DATA", "weights": "WEIGHT_SPECTRUM",
                    "scans": test.measured,
                }))
            else:
                # As the split itself does: what is there is not built upon.
                if out.exists():
                    shutil.rmtree(out)
                out.mkdir()

        def casa_answers(casa, ms):
            test.asked.append(casa)
            return casa.startswith("good"), "6.7.6"

        imager = types.SimpleNamespace(
            run_command=run_command, inside=RECIPE.IMAGER.inside,
            claimed=RECIPE.IMAGER.claimed,
        )
        recipe = types.SimpleNamespace(
            IMAGER=imager, CASA_CANDIDATES=("bad python3", "good python3"),
            casa_answers=casa_answers,
        )
        for patch in (
            mock.patch.object(TOOL, "load_sibling", lambda *a: recipe),
            mock.patch.object(TOOL.shutil, "which", lambda name: name),
            mock.patch("builtins.print"),
        ):
            patch.start()
            self.addCleanup(patch.stop)

    def run_tool(self, *extra):
        return TOOL.main([
            "--target", "T", "--ms", str(self.ms), "--workdir", str(self.work),
            "--channel-range", "51", "768", "--outdir", str(self.products),
            *extra,
        ])

    def report(self):
        return json.loads((self.work / TOOL.REPORT_NAME).read_text())

    def parts(self):
        return [option(command, "--casa-part") for command in self.commands]

    def test_scans_left_out_give_a_measurement_set_of_those_kept(self):
        self.run_tool()

        self.assertEqual(self.parts(), ["measure", "split"])
        measure, split = self.commands
        # The first command that answers as CASA runs this script again.
        self.assertEqual(measure[:2], ["good", "python3"])
        self.assertEqual(Path(measure[2]).name, "select_scans.py")
        self.assertEqual(option(measure, "--ms"), str(self.ms))
        self.assertEqual(option(measure, "--channel-range", 2), ["51", "768"])
        self.assertEqual(option(split, "--scans"), "15")
        self.assertEqual(option(split, "--out"), str(self.work / "T.ms"))
        report = self.report()
        self.assertEqual(report["kept"], [15])
        self.assertEqual(report["left_out"], [3, 29])
        self.assertEqual(report["main_scan"], 15)
        self.assertEqual(report["measurement_set"], str(self.work / "T.ms"))
        self.assertEqual(report["input"], str(self.ms))
        self.assertEqual(report["rule"], {
            "least_share": 0.1, "most_excess": 1.4, "scans": None})
        self.assertLess(report["noise_kept_jy"], report["noise_every_scan_jy"])
        self.assertEqual(len(report["reading"]), 4)

    def test_every_scan_kept_goes_on_with_the_data_given(self):
        self.measured = [scan(5, 22000, 1.0e-2), scan(18, 22000, 1.1e-2)]

        self.run_tool()

        self.assertEqual(self.parts(), ["measure"])
        report = self.report()
        self.assertEqual(report["measurement_set"], str(self.ms))
        self.assertEqual(report["left_out"], [])
        self.assertFalse((self.work / "T.ms").exists())

    def test_averaging_makes_a_copy_even_with_every_scan_kept(self):
        self.measured = [scan(5, 22000, 1.0e-2), scan(18, 22000, 1.1e-2)]

        self.run_tool("--average-channels", "4")

        self.assertEqual(self.parts(), ["measure", "split"])
        measure, split = self.commands
        self.assertNotIn("--average-channels", measure)
        self.assertEqual(option(split, "--average-channels"), "4")
        self.assertEqual(option(split, "--scans", 2), ["5", "18"])
        report = self.report()
        self.assertEqual(report["measurement_set"], str(self.work / "T.ms"))
        self.assertEqual(report["rule"]["average_channels"], 4)
        self.assertEqual(report["split"], {"scans": [5, 18], "average_channels": 4})

        # Given again, nothing is redone; another averaging splits again.
        self.run_tool("--average-channels", "4")
        self.assertEqual(self.parts(), ["measure", "split"])
        self.run_tool("--average-channels", "2")
        self.assertEqual(self.parts(), ["measure", "split", "split"])

    def test_given_again_nothing_is_measured_or_split_again(self):
        self.run_tool()
        first = self.report()
        asked = len(self.asked)

        self.run_tool()

        self.assertEqual(self.parts(), ["measure", "split"])
        self.assertEqual(len(self.asked), asked)
        self.assertEqual(self.report()["kept"], first["kept"])

    def test_other_scans_named_are_split_without_measuring_again(self):
        self.run_tool()

        self.run_tool("--scans", "3", "15")

        self.assertEqual(self.parts(), ["measure", "split", "split"])
        self.assertEqual(option(self.commands[-1], "--scans", 2), ["3", "15"])
        report = self.report()
        self.assertEqual(report["kept"], [3, 15])
        self.assertEqual(report["rule"]["scans"], [3, 15])
        self.assertIn("named with --scans though", report["reading"][0])

    def test_a_split_that_did_not_finish_is_made_again(self):
        self.run_tool()
        (self.work / "T.ms").rmdir()

        self.run_tool()

        self.assertEqual(self.parts(), ["measure", "split", "split"])

    def test_the_limits_given_decide(self):
        self.run_tool("--least-share", "0.01", "--most-excess", "2")

        self.assertEqual(self.parts(), ["measure"])
        self.assertEqual(self.report()["kept"], [3, 15, 29])

    def test_a_given_casa_is_the_only_one_asked(self):
        self.run_tool("--casa", "good own")

        self.assertEqual(self.asked, ["good own"])
        self.assertEqual(self.commands[0][:2], ["good", "own"])

    def test_no_casa_is_said_with_what_was_tried(self):
        with self.assertRaises(SystemExit) as stopped:
            self.run_tool("--casa", "bad one")

        self.assertIn("no command ran a Python that imports casatasks",
                      str(stopped.exception))
        self.assertEqual(self.commands, [])

    def test_other_data_in_the_same_folder_are_measured_afresh(self):
        self.run_tool()
        other = self.ms.with_name("U.ms")
        other.mkdir()

        TOOL.main([
            "--target", "T", "--ms", str(other), "--workdir", str(self.work),
            "--channel-range", "51", "768", "--outdir", str(self.products),
        ])

        self.assertEqual(self.parts(), ["measure", "split", "measure", "split"])
        self.assertEqual(self.report()["input"], str(other))

    def test_a_folder_inside_the_products_is_refused(self):
        self.work = self.products / "scans"

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("inside the pipeline's products", str(stopped.exception))
        self.assertEqual(self.commands, [])

    def test_missing_data_are_said_so(self):
        self.ms.rmdir()

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("No Measurement Set at", str(stopped.exception))

    def test_nothing_in_the_products_is_written(self):
        before = sorted(str(path) for path in self.products.rglob("*"))

        self.run_tool()

        self.assertEqual(
            sorted(str(path) for path in self.products.rglob("*")), before)


if __name__ == "__main__":
    unittest.main()
