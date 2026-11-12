import contextlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_radio_continuum_diagnose_flag_losses",
    SOURCE_DIR / "diagnose_flag_losses.py",
)
DIAGNOSE = importlib.util.module_from_spec(spec)
spec.loader.exec_module(DIAGNOSE)


def counts(fraction, total=1000):
    return {"flagged": int(round(fraction * total)), "total": total}


def summary(field_fractions, scan_fractions):
    return {
        "flagged": 0,
        "total": 1,
        "field": {name: counts(value) for name, value in field_fractions.items()},
        "scan": {scan: counts(value) for scan, value in scan_fractions.items()},
    }


class HelperTests(unittest.TestCase):

    def test_casa_log_time_round_trips(self):
        seconds = DIAGNOSE.casa_time_to_mjd_seconds("2020/11/22/14:05:38.4")
        self.assertEqual(
            DIAGNOSE.mjd_seconds_to_utc(seconds),
            "2020/11/22/14:05:38.4",
        )

    def test_selector_expansion(self):
        channels = DIAGNOSE.parse_spw_selection("0:0~2;7,1:4~5")
        self.assertEqual(channels["0"], {0, 1, 2, 7})
        self.assertEqual(channels["1"], {4, 5})

    def test_channel_runs(self):
        fractions = {0: 1.0, 1: 1.0, 2: 0.2, 3: 0.995, 4: 1.0, 6: 1.0}
        self.assertEqual(
            DIAGNOSE.channel_runs(fractions),
            [(0, 1), (3, 4), (6, 6)],
        )

    def test_run_classification_allows_hanning_neighbour(self):
        deliberate = set(range(0, 5))
        self.assertEqual(
            DIAGNOSE.classify_run((0, 5), deliberate),
            "deliberate mask",
        )
        self.assertTrue(
            DIAGNOSE.classify_run((20, 24), deliberate).startswith(
                "data-driven"
            )
        )
        self.assertTrue(
            DIAGNOSE.classify_run((3, 9), deliberate).startswith("mixed")
        )

    def test_loss_attribution(self):
        self.assertEqual(
            DIAGNOSE.attribute_loss([
                ("on arrival", 0.95),
                ("after extension", 0.99),
            ]),
            "already flagged on arrival",
        )
        self.assertEqual(
            DIAGNOSE.attribute_loss([
                ("on arrival", 0.05),
                ("after fixed masks", 0.30),
                ("after auto-flag 1", 0.40),
                ("after extension", 0.97),
            ]),
            "after extension",
        )


class OutdirResolutionTests(unittest.TestCase):

    def test_relative_outdir_is_found_next_to_the_script(self):
        with tempfile.TemporaryDirectory() as script_dir, \
                tempfile.TemporaryDirectory() as elsewhere:
            report = Path(script_dir) / "Products" / "inspection"
            report.mkdir(parents=True)
            (report / "inspection_report.json").write_text("{}")

            cwd = Path.cwd()
            try:
                import os
                os.chdir(elsewhere)
                found, searched = DIAGNOSE.resolve_outdir(
                    "Products",
                    script_dir=script_dir,
                )
            finally:
                os.chdir(cwd)

            self.assertEqual(found, (Path(script_dir) / "Products").resolve())
            self.assertEqual(len(searched), 2)

    def test_missing_output_explains_where_it_looked(self):
        with tempfile.TemporaryDirectory() as empty:
            with self.assertRaises(SystemExit) as raised:
                DIAGNOSE.main(["--outdir", str(Path(empty) / "Products")])
        self.assertIn("Looked in:", str(raised.exception))


class ReportTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.outdir = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, relative, value):
        path = self.outdir / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))

    def test_band4_style_losses_are_attributed(self):
        start = DIAGNOSE.casa_time_to_mjd_seconds
        spans = {
            "1": {"fields": [0], "start": start("2020/11/22/13:54:16.6"),
                  "end": start("2020/11/22/14:01:15.3")},
            "2": {"fields": [1], "start": start("2020/11/22/14:05:30.0"),
                  "end": start("2020/11/22/14:10:12.2")},
            "3": {"fields": [2], "start": start("2020/11/22/14:12:00.0"),
                  "end": start("2020/11/22/14:50:00.0")},
        }
        # 100 channels over 550-950 MHz, ascending.
        frequencies = [550.0 + (channel + 0.5) * 4.0 for channel in range(100)]
        self.write("inspection/inspection_report.json", {
            "flux_calibrators": ["3C147"],
            "phase_calibrators": ["J0000+0000"],
            "target_fields": ["TARGET"],
            "spws": {"0": {"nchan": 100,
                           "channel_frequencies_MHz": frequencies}},
            "time_coverage": {"scan_spans": spans},
        })
        self.write("inspection/fields.json", [
            {"id": 0, "name": "3C147"},
            {"id": 1, "name": "J0000+0000"},
            {"id": 2, "name": "TARGET"},
        ])

        fields = {"3C147": 0.1, "J0000+0000": 0.1, "TARGET": 0.1}
        self.write("flagging/preflag_summary.json",
                   summary(fields, {"1": 0.96, "2": 0.05, "3": 0.05}))
        self.write("flagging/iterative_baseline_summary.json",
                   summary(fields, {"1": 0.97, "2": 0.30, "3": 0.30}))
        self.write("flagging/iter_1_summary.json",
                   summary(fields, {"1": 0.97, "2": 0.40, "3": 0.35}))
        self.write("flagging/postflag_summary.json",
                   summary(fields, {"1": 0.98, "2": 0.97, "3": 0.36}))

        # Deliberate: 5% edges and everything above 850 MHz (channel 75+).
        self.write("flagging/edge_flagging_report.json", {
            "spw_string": "0:75~99,0:0~4;95~99",
            "guaranteed_science_range_mask": {"spw_selection": "0:75~99"},
            "percentage_edge_mask": {"spw_selection": "0:0~4;95~99"},
        })
        channel = {}
        for index in range(100):
            flagged = index <= 5 or 60 <= index <= 64 or index >= 75
            channel[f"0:{index}"] = counts(1.0 if flagged else 0.1)
        self.write("flagging/spw_statistics.json",
                   {"by_spw": {}, "by_spw_channel": {"spw:channel": channel}})

        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            DIAGNOSE.main([
                "--outdir", str(self.outdir),
                "--time", "2020/11/22/14:05:38.4",
            ])
        text = output.getvalue()

        scan_one = text.split("scan 1:")[1].split("scan 2:")[0]
        self.assertIn("already flagged on arrival", scan_one)

        scan_two = text.split("scan 2:")[1].split("Fully flagged")[0]
        self.assertIn("contains --time", scan_two)
        self.assertIn("phase cal", scan_two)
        self.assertIn("most data removed: after extension", scan_two)

        self.assertNotIn("scan 3:", text)

        runs = text.split("Fully flagged channel runs")[1]
        self.assertRegex(runs, r"chan\s+0-5\s.*deliberate mask")
        self.assertRegex(runs, r"chan\s+60-64\s.*data-driven")
        self.assertRegex(runs, r"chan\s+75-99\s.*deliberate mask")
        self.assertIn("792.000", runs)


if __name__ == "__main__":
    unittest.main()
