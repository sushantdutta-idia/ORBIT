"""The self-calibration report says what its final image reached.

In the uGMRT Band 4 run that prompted this, self-calibration finished with
three accepted loops and its report held ``"final_metrics": {}``: the
measurement sat after an early return in ``generate_qa_plots``.
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from test_selfcal_resume_safety import load_selfcal_module


SELFCAL = load_selfcal_module()

FINAL = {
    "max": 0.0231,
    "rms": 8.6e-5,
    "dynamic_range": 268.6,
}


class FinalMetricsTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.manager = SELFCAL.SelfCalManager.__new__(SELFCAL.SelfCalManager)
        self.manager.reportdir = self.root
        self.manager.generate_native_qa_plots = lambda: {"plots": ["a.png"]}

    def tearDown(self):
        self.temporary.cleanup()

    def generate(self, measure):
        with mock.patch.object(
            SELFCAL.CASA_COMPAT, "peak_and_noise", measure
        ), mock.patch("builtins.print") as printed:
            status = self.manager.generate_qa_plots()
        return status, " ".join(
            str(argument) for call in printed.call_args_list
            for argument in call.args
        )

    def test_the_final_image_is_measured_before_the_plots(self):
        self.manager.final_image = str(self.root / "final.image.tt0")
        seen = []

        def measure(imstat, image):
            seen.append(image)
            return dict(FINAL)

        status, _ = self.generate(measure)

        self.assertEqual(seen, [self.manager.final_image])
        self.assertEqual(
            self.manager.final_metrics,
            {
                "peak_Jy_beam": 0.0231,
                "rms_Jy_beam": 8.6e-5,
                "dynamic_range": 268.6,
            },
        )
        written = json.loads((self.root / "selfcal_metrics.json").read_text())
        self.assertEqual(written, self.manager.final_metrics)
        # The figures are still made.
        self.assertEqual(status, {"plots": ["a.png"]})

    def test_without_a_final_image_nothing_is_claimed(self):
        def measure(imstat, image):
            raise AssertionError("no image to measure")

        status, _ = self.generate(measure)

        self.assertFalse(hasattr(self.manager, "final_metrics"))
        self.assertFalse((self.root / "selfcal_metrics.json").exists())
        self.assertEqual(status, {"plots": ["a.png"]})

    def test_an_unreadable_image_does_not_stop_the_stage(self):
        self.manager.final_image = str(self.root / "final.image.tt0")

        def measure(imstat, image):
            raise RuntimeError("image is locked")

        status, printed = self.generate(measure)

        self.assertFalse(hasattr(self.manager, "final_metrics"))
        self.assertIn("could not measure", printed)
        self.assertEqual(status, {"plots": ["a.png"]})


if __name__ == "__main__":
    unittest.main()
