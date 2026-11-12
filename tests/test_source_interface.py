import importlib.util
import tempfile
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "test_radio_continuum_pipeline_source_interface",
    SOURCE_DIR / "pipeline.py",
)
PIPELINE = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PIPELINE)


class SourceInterfaceTests(unittest.TestCase):

    def test_shipped_files_satisfy_the_orchestrator(self):
        self.assertEqual(PIPELINE.outdated_source_files(SOURCE_DIR), [])

    def test_both_deconvolving_stages_are_covered(self):
        # Files are copied to the cluster one at a time. Both stage scripts
        # that deconvolve must be checked, or half a release can run.
        for name in ("imaging.py", "selfcal.py"):
            with self.subTest(name=name):
                self.assertIn(name, PIPELINE.SOURCE_INTERFACE_VERSIONS)
                declared = [
                    line for line in
                    (SOURCE_DIR / name).read_text().splitlines()
                    if line.startswith("PIPELINE_INTERFACE_VERSION")
                ]
                self.assertEqual(len(declared), 1, name)

    def test_stale_missing_and_unmarked_copies_are_named(self):
        required = {"current.py": 2, "older.py": 2, "unmarked.py": 1,
                    "absent.py": 1}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "current.py").write_text("PIPELINE_INTERFACE_VERSION = 3\n")
            (root / "older.py").write_text("PIPELINE_INTERFACE_VERSION = 1\n")
            (root / "unmarked.py").write_text(
                "def f():\n    PIPELINE_INTERFACE_VERSION = 9\n"
            )
            problems = PIPELINE.outdated_source_files(root, required)

        self.assertEqual(len(problems), 3)
        self.assertIn("older.py is from an earlier release (interface 1, needs 2)",
                      problems)
        self.assertIn(
            "unmarked.py is from an earlier release (interface unmarked, needs 1)",
            problems,
        )
        self.assertIn("absent.py is missing", problems)

    def test_preflight_stops_before_any_stage(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = PIPELINE.build_parser().parse_args([
                "--outdir", str(Path(temporary) / "Products"),
                "--telescope", "ugmrt",
                "--receiver-band", "B4",
                "--stages", "imaging",
                "--dry-run",
            ])
            runner = PIPELINE.PipelineRunner(args)
            runner.root = Path(temporary)
            with self.assertRaisesRegex(RuntimeError, "imaging.py is missing"):
                runner.preflight(["imaging"])


if __name__ == "__main__":
    unittest.main()
