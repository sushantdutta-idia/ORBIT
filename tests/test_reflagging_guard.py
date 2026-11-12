"""An MS already flagged by an earlier run is not flagged again unasked."""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_flagging_selectors import FLAGGING


def flagger(ms):
    value = FLAGGING.CASAFlagger.__new__(FLAGGING.CASAFlagger)
    value.ms = str(ms)
    return value


def version_list(ms, *lines):
    directory = Path(str(ms) + ".flagversions")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "FLAG_VERSION_LIST").write_text("\n".join(lines) + "\n")


class PreviouslyFlaggedInputTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.ms = Path(self.temporary.name) / "working.ms"
        self.ms.mkdir()

    def tearDown(self):
        self.temporary.cleanup()

    def test_a_fresh_copy_is_flagged_normally(self):
        self.assertIsNone(flagger(self.ms).refuse_previously_flagged_input())

    def test_versions_from_other_software_are_not_ours(self):
        version_list(self.ms, "Original : Original flags at import into CASA",
                     "flagdata_1 : Flags autosave on 2026-01-01")
        self.assertEqual(FLAGGING.earlier_pipeline_flag_versions(self.ms), [])
        self.assertIsNone(flagger(self.ms).refuse_previously_flagged_input())

    def test_an_earlier_runs_flags_stop_a_fresh_run(self):
        version_list(
            self.ms,
            "before_flagging : radio continuum pipeline flag state: "
            "before_flagging",
            "iter_1 : radio continuum pipeline flag state: iter_1",
            "before_ugmrt_postcal_4242_1 : Before corrected-data tfcrop/rflag "
            "and round-2 solves",
        )
        self.assertEqual(
            FLAGGING.earlier_pipeline_flag_versions(self.ms),
            ["before_flagging", "iter_1", "before_ugmrt_postcal_4242_1"],
        )
        with self.assertRaisesRegex(RuntimeError, "earlier run") as raised:
            flagger(self.ms).refuse_previously_flagged_input()
        message = str(raised.exception)
        self.assertIn("working_ms", message)
        # Never advise deleting what may be the observer's only copy.
        self.assertIn("do not delete it", message)

    def test_it_can_be_overridden_on_purpose(self):
        version_list(
            self.ms,
            "before_flagging : radio continuum pipeline flag state: "
            "before_flagging",
        )
        with mock.patch.dict(
            FLAGGING.FLAGGING_DEFAULTS,
            {"allow_previously_flagged_input": True},
        ), mock.patch("builtins.print"):
            self.assertEqual(
                flagger(self.ms).refuse_previously_flagged_input(),
                ["before_flagging"],
            )


class AttributionTests(unittest.TestCase):

    def test_the_observers_run(self):
        # 42.7% arrived flagged, 50.5% at the end of the flagging stage.
        attribution = FLAGGING.CASAFlagger.flagging_attribution(
            {"flagged": 427, "total": 1000},
            {"flagged": 505, "total": 1000},
        )
        self.assertAlmostEqual(attribution["input_flagged_fraction"], 0.427)
        self.assertAlmostEqual(
            attribution["this_run_fraction_of_received"],
            (0.505 - 0.427) / (1 - 0.427),
        )

    def test_missing_summaries_give_nothing(self):
        self.assertIsNone(
            FLAGGING.CASAFlagger.flagging_attribution({}, {"flagged": 1,
                                                           "total": 2})
        )


if __name__ == "__main__":
    unittest.main()
