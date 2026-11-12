import importlib.util
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_module():
    path = SOURCE_DIR / "casa_compat.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_casa_compat",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


COMPAT = load_module()


class CASACompatibilitySafetyTests(unittest.TestCase):

    def test_signature_mismatch_does_not_drop_science_keyword(self):
        def task(vis=""):
            return vis

        wrapped = COMPAT.compat_task(task, name="flagdata")

        with self.assertRaisesRegex(RuntimeError, "required science keyword"):
            wrapped(vis="input.ms", extendflags=False)

    def test_runtime_unknown_keyword_does_not_retry_weakened_call(self):
        calls = []

        def task(**kwargs):
            calls.append(dict(kwargs))
            raise TypeError("unexpected keyword argument 'calwt'")

        wrapped = COMPAT.compat_task(task, name="applycal")

        with self.assertRaisesRegex(RuntimeError, "calwt"):
            wrapped(vis="input.ms", calwt=True)

        self.assertEqual(len(calls), 1)

    def test_explicit_cosmetic_keyword_may_be_dropped(self):
        def task(vis=""):
            return vis

        wrapped = COMPAT.compat_task(
            task,
            name="plotms",
            optional_keywords={"highres"},
        )

        self.assertEqual(wrapped(vis="input.ms", highres=True), "input.ms")


if __name__ == "__main__":
    unittest.main()
