"""Observation-specific manual flags, applied before the automatic flaggers.

The case is uGMRT Band 5 of 2021-02-22: on the run of 2026-10-25 the
amplitude solves on 3C48 gave SNR < 3 for 44-46 of 52 solutions in each of
the first five 2-minute intervals of its 27-minute scan (phases were fine),
and the bandpass was solved over the whole scan. Scan 10 of 0323+055 lost
half its 1300 MHz flux and most of its amplitude solutions in its second
half.
"""

import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_flagging_selectors import FLAGGING


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_manual_flags_config", SOURCE_DIR / "config.py")
CONFIG = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CONFIG)

START = ("mode='manual' field='3C48' "
         "timerange='2021/02/22/07:36:00~2021/02/22/07:47:00' "
         "reason='array settling at the start of the observation'")
SCAN_10 = "mode='manual' scan='10' reason='0323+055 scan 10 corrupted'"


class ParseTests(unittest.TestCase):

    def test_a_command_gives_its_selection_and_reason(self):
        selection, reason = FLAGGING.parse_manual_flag_command(START)
        self.assertEqual(selection, {
            "field": "3C48",
            "timerange": "2021/02/22/07:36:00~2021/02/22/07:47:00",
        })
        self.assertEqual(reason, "array settling at the start of the "
                         "observation")
        selection, reason = FLAGGING.parse_manual_flag_command(SCAN_10)
        self.assertEqual(selection, {"scan": "10"})

    def test_mode_may_be_left_out(self):
        self.assertEqual(
            FLAGGING.parse_manual_flag_command("antenna='W06'"),
            ({"antenna": "W06"}, None),
        )

    def test_what_is_not_a_manual_selection_is_refused(self):
        for bad in (
            "mode='manual'",                      # flags everything
            "mode='manual' reason='bad'",         # still everything
            "mode='tfcrop' field='3C48'",         # not a manual flag
            "mode='manual' fields='3C48'",        # misspelt key
            "mode='manual' scan='10' scan='11'",  # given twice
            "mode='manual' scan=''",              # empty selection
            "mode='manual' scan",                 # not key='value'
            "mode='manual' field='3C48",          # unclosed quote
        ):
            with self.subTest(command=bad), self.assertRaises(ValueError):
                FLAGGING.parse_manual_flag_command(bad)


class ConfigTests(unittest.TestCase):

    def test_config_yaml_gives_a_list(self):
        self.assertEqual(
            CONFIG.parse_manual_flag_commands([START, SCAN_10]),
            [START, SCAN_10])
        self.assertEqual(CONFIG.parse_manual_flag_commands(SCAN_10),
                         [SCAN_10])
        self.assertEqual(CONFIG.parse_manual_flag_commands(None), [])
        self.assertEqual(CONFIG.parse_manual_flag_commands([]), [])
        with self.assertRaises(ValueError):
            CONFIG.parse_manual_flag_commands({"scan": "10"})

    def test_the_one_line_form_survives_the_fallback_parser(self):
        # Without PyYAML a list is split at commas, so commands carry none.
        line = f'["{START}", "{SCAN_10}"]'
        self.assertEqual(
            CONFIG.parse_manual_flag_commands(
                CONFIG._parse_config_yaml_scalar(line)),
            [START, SCAN_10])

    def test_without_commands_the_recipe_has_no_key(self):
        # Earlier runs' flagging recipes hash FLAGGING_CONFIG; an empty
        # default would have changed every one of them.
        self.assertNotIn("manual_commands", CONFIG.FLAGGING_CONFIG)


class Summary:
    """flagdata summary and manual flagging on a pretend observation."""

    def __init__(self, sizes):
        # sizes: selection value -> (flagged, total)
        self.sizes = dict(sizes)
        self.applied = []

    def __call__(self, **kwargs):
        key = tuple(sorted((k, v) for k, v in kwargs.items()
                           if k in FLAGGING.MANUAL_FLAG_SELECTION_KEYS))
        if kwargs["mode"] == "summary":
            flagged, total = self.sizes.get(key, (0, 0))
            return {"flagged": flagged, "total": total}
        if kwargs["mode"] == "manual":
            self.applied.append(dict(key))
            flagged, total = self.sizes.get(key, (0, 0))
            self.sizes[key] = (total, total)
            return None
        raise AssertionError(kwargs["mode"])


class FlagStepTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def flagger(self, commands):
        value = FLAGGING.CASAFlagger.__new__(FLAGGING.CASAFlagger)
        value.ms = str(self.root / "obs.ms")
        value.root = self.root
        value.manual_flag_commands = list(commands)
        value.flagdata_ntime = lambda: "scan"
        return value

    def run_step(self, flagger, fake, total=1000):
        with mock.patch.object(FLAGGING, "flagdata", fake), \
                mock.patch("builtins.print"):
            return flagger.flag_manual_commands(total)

    def test_each_command_is_counted_applied_and_reported(self):
        fake = Summary({
            (("field", "3C48"),
             ("timerange", "2021/02/22/07:36:00~2021/02/22/07:47:00")):
                (10, 100),
            (("scan", "10"),): (30, 60),
        })
        report = self.run_step(self.flagger([START, SCAN_10]), fake)
        self.assertEqual(report["status"], "applied")
        self.assertEqual(fake.applied, [
            {"field": "3C48",
             "timerange": "2021/02/22/07:36:00~2021/02/22/07:47:00"},
            {"scan": "10"},
        ])
        first, second = report["commands"]
        self.assertEqual(first["selected_visibilities"], 100)
        self.assertAlmostEqual(first["share_of_observation"], 0.1)
        self.assertAlmostEqual(first["flagged_before"], 0.1)
        self.assertAlmostEqual(first["flagged_after"], 1.0)
        self.assertAlmostEqual(second["flagged_before"], 0.5)
        self.assertEqual(second["reason"], "0323+055 scan 10 corrupted")
        self.assertTrue((self.root / "manual_flag_report.json").is_file())

    def test_a_command_that_selects_nothing_stops_the_stage(self):
        fake = Summary({})
        with self.assertRaisesRegex(RuntimeError, "selects no data"):
            self.run_step(self.flagger([SCAN_10]), fake)
        self.assertEqual(fake.applied, [])

    def test_a_command_selecting_most_of_the_data_is_not_applied(self):
        fake = Summary({(("field", "0323+055"),): (0, 400)})
        with self.assertRaisesRegex(RuntimeError, "40.0% of the observation"):
            self.run_step(self.flagger(["field='0323+055'"]), fake)
        self.assertEqual(fake.applied, [])

    def test_without_commands_nothing_is_flagged(self):
        fake = Summary({})
        report = self.run_step(self.flagger([]), fake)
        self.assertEqual(report, {"status": "none", "commands": []})
        self.assertEqual(fake.applied, [])


if __name__ == "__main__":
    unittest.main()
