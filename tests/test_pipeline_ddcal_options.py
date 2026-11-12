import tempfile
import unittest
from pathlib import Path

from test_resume_safety import PIPELINE as pipeline


class PipelineDDCalIntegrationTests(unittest.TestCase):

    def test_generated_sbatch_embeds_dd_options(self):
        with tempfile.TemporaryDirectory() as temporary:
            outfile = Path(temporary) / "ddcal.sbatch"
            args = pipeline.build_parser().parse_args(
                [
                    "/cluster/data/target.ms",
                    "--outdir",
                    "/cluster/scratch/target_products",
                    "--include-ddcal",
                    "--dd-mode",
                    "survey",
                    "--dd-engine",
                    "killms-ddfacet",
                    "--execute-ddcal",
                    "--write-sbatch",
                    str(outfile),
                ]
            )

            pipeline.write_sbatch_from_config(args)
            content = outfile.read_text()

            self.assertIn("/cluster/data/target.ms", content)
            self.assertIn(
                "--outdir /cluster/scratch/target_products",
                content,
            )
            self.assertIn("--include-ddcal", content)
            self.assertIn("--dd-mode survey", content)
            self.assertIn("--dd-engine killms-ddfacet", content)
            self.assertIn("--execute-ddcal", content)

    def test_dd_stage_receives_mode_engine_and_execution_policy(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = pipeline.build_parser().parse_args(
                [
                    "--outdir",
                    temporary,
                    "--stages",
                    "ddcal",
                    "--include-ddcal",
                    "--dd-mode",
                    "single",
                    "--dd-engine",
                    "auto",
                    "--execute-ddcal",
                    "--skip-dd-post-imaging",
                    "--dry-run",
                ]
            )
            runner = pipeline.PipelineRunner(args)
            command = runner.command_for_stage("ddcal")

            self.assertIn("--mode", command)
            self.assertIn("single", command)
            self.assertIn("--engine", command)
            self.assertIn("auto", command)
            self.assertIn("--execute", command)
            self.assertIn("--skip-imaging", command)


if __name__ == "__main__":
    unittest.main()
