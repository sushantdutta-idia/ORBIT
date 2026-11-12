import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_tool():
    spec = importlib.util.spec_from_file_location(
        "test_selfcal_recipe_module", SOURCE_DIR / "selfcal_recipe.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


TOOL = load_tool()
CASA = "singularity exec /containers/casa.sif python3"
FOUND = {"wsclean_version": [3, 1], "casa": CASA, "casa_version": "6.6.4"}


def arguments(*extra):
    captured = []
    with mock.patch.object(
        TOOL, "run", lambda args, plan=None: captured.append(args)
    ):
        TOOL.main([
            "--target", "G1", "--ms", "data.ms", "--imsize", "5120",
            "--cell", "1.0", "--channel-range", "51", "768",
            "--casa", CASA, *extra,
        ])
    return captured[0]


def option(command, name, count=1):
    index = command.index(name)
    values = command[index + 1:index + 1 + count]
    return values[0] if count == 1 else values


def solve_report(phase, amplitude=None, delay=None, applied=True):
    """What experiment_phase_offsets.py reports for one offset antenna."""

    quiet = {"mean_phase_deg": 0.5, "mean_amplitude": 1.0}
    offset = {
        "mean_phase_deg": -float(phase),
        "mean_amplitude": 1.0 if amplitude is None else float(amplitude),
    }
    report = {
        "solve": {"calmode": "p" if amplitude is None else "ap"},
        "summary": {"solved_antennas": 2 if applied else 0},
        "series": [
            {"antenna": "C00", "polarizations": [quiet, quiet]},
            {"antenna": "S02", "polarizations": [offset, {"failed": 1}]},
        ],
        "delays": None if delay is None else [
            {"antenna": "C00", "delay_ns": [0.0, None]},
            {"antenna": "S04", "delay_ns": [float(delay), -0.1]},
        ],
    }
    if applied:
        report["corrected_measurement_set"] = "somewhere"
    return report


def band_report(phase=10.5, amplitude=0.326, applied=True):
    """What the solve per channel bin reports: W04 runs across the band."""

    steady = {"bins": 25, "phase_rms_deg": 1.0,
              "phase_rms_after_line_deg": 0.9, "amplitude_rms": 0.02}
    running = {"bins": 22, "phase_rms_deg": float(phase),
               "phase_rms_after_line_deg": float(phase) - 0.2,
               "amplitude_rms": float(amplitude)}
    report = {
        "solve": {"bandpass_channels": 30, "amplitudes_applied": True},
        "statistics": [
            {"antenna": "C00", "polarizations": [steady, steady]},
            {"antenna": "W04", "polarizations": [running, {"bins": 0}]},
        ],
    }
    if applied:
        report["corrected_measurement_set"] = "somewhere"
    return report


IMAGE_REPORT = {
    "stokes_v": {"robust_rms_Jy_beam": 23.5e-6},
    "stokes_i": {
        "image": "final-MFS-image.fits",
        "robust_rms_Jy_beam": 30.1e-6,
        "max_Jy_beam": 39.2e-3,
        "min_Jy_beam": -0.48e-3,
        "grid": {"quietest_tenth_rms": 24.9e-6},
    },
    "source_asymmetry": {
        "deepest_negative_median": -0.012, "antisymmetric_peak": 0.032,
    },
}
CORE_REPORT = {
    "dirty_asymmetry": {"antisymmetric_peak": 0.15, "antisymmetric_rms": 0.039},
    "model": {"total_flux_Jy": 0.21},
}


class StartTests(unittest.TestCase):

    def test_it_starts_beside_the_pipelines_own_inspect_module(self):
        source = (SOURCE_DIR / "selfcal_recipe.py").read_text()
        self.assertIn("sys.path.pop(0)", source.split("import argparse")[0])
        done = subprocess.run(
            [sys.executable, "selfcal_recipe.py", "--help"],
            cwd=SOURCE_DIR, capture_output=True, text=True,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("--channel-range", done.stdout)


class CommandTests(unittest.TestCase):

    def test_a_cores_only_model_is_asked_of_the_imaging_script(self):
        command = TOOL.model_command(
            arguments(), "work/round0.ms", "work/images/phase1",
            {"core": "phase"},
        )

        self.assertEqual(command[0], sys.executable)
        self.assertTrue(command[1].endswith("experiment_wsclean_image.py"))
        self.assertEqual(option(command, "--ms"), "work/round0.ms")
        self.assertEqual(option(command, "--core-model"), "6")
        self.assertEqual(option(command, "--channel-range", 2), ["51", "768"])
        self.assertEqual(option(command, "--imsize"), "5120")
        self.assertEqual(option(command, "--threads"), "6")
        # The container command stays one argument for the script to split.
        self.assertEqual(option(command, "--wsclean"), TOOL.DEFAULT_WSCLEAN)
        self.assertNotIn("--write-model", command)

    def test_the_threshold_of_the_cores_can_be_set_for_another_field(self):
        command = TOOL.model_command(
            arguments("--core-threshold", "2.5"), "a.ms", "images",
            {"core": "phase"},
        )

        self.assertEqual(option(command, "--core-model"), "2.5")

    def test_the_model_of_the_whole_field_is_a_usual_clean_and_written(self):
        command = TOOL.model_command(
            arguments(), "a.ms", "images", {"full": True}
        )

        self.assertIn("--write-model", command)
        # Positive components only would end 15 to 40% above the data
        # around resolved sources, and bias the amplitudes solved on it.
        self.assertNotIn("--no-negative", command)
        # Multi-scale with a local noise map aborts in WSClean 3.1.0.
        self.assertIn("--no-local-rms", command)
        self.assertEqual(
            option(command, "--scales", 4), ["0", "5", "15", "30"]
        )
        self.assertEqual(option(command, "--auto-threshold"), "1")
        self.assertEqual(option(command, "--nmiter"), "20")
        self.assertNotIn("--core-model", command)

    def test_the_final_image_writes_nothing_to_the_data(self):
        command = TOOL.image_command(arguments(), "a.ms", "images")

        self.assertNotIn("--write-model", command)
        self.assertNotIn("--core-model", command)
        self.assertNotIn("--no-negative", command)
        self.assertEqual(option(command, "--channels-out"), "12")

    def test_the_final_image_can_be_larger_than_the_models(self):
        args = arguments("--final-imsize", "5120")
        args.imsize = 3072

        final = TOOL.image_command(args, "a.ms", "images", args.final_imsize)
        model = TOOL.model_command(args, "a.ms", "images", {"full": True})

        self.assertEqual(option(final, "--imsize"), "5120")
        self.assertEqual(option(model, "--imsize"), "3072")
        # Without the option the final image is as large as the models.
        self.assertEqual(
            option(TOOL.image_command(args, "a.ms", "images"), "--imsize"),
            "3072",
        )

    def test_without_a_channel_range_every_channel_is_imaged(self):
        args = arguments()
        args.channel_range = None

        command = TOOL.image_command(args, "a.ms", "images")

        self.assertIn("--all-channels", command)
        self.assertNotIn("--channel-range", command)

    def test_a_solve_is_run_by_the_python_that_has_casa(self):
        command = TOOL.solve_command(
            arguments(), "work/round0.ms", {"calmode": "p"}
        )

        self.assertEqual(
            command[:4],
            ["singularity", "exec", "/containers/casa.sif", "python3"],
        )
        self.assertTrue(command[4].endswith("experiment_phase_offsets.py"))
        self.assertEqual(option(command, "--ms"), "work/round0.ms")
        self.assertEqual(option(command, "--calmode"), "p")
        self.assertEqual(option(command, "--solint"), "inf")
        self.assertEqual(option(command, "--suffix"), "_solved")
        self.assertIn("--overwrite", command)
        self.assertNotIn("--delays", command)
        self.assertNotIn("--uvrange", command)

    def test_delays_and_a_baseline_range_are_passed_on(self):
        command = TOOL.solve_command(
            arguments(), "a.ms",
            {"calmode": "ap", "delays": True, "uvrange": ">3klambda"},
        )

        self.assertIn("--delays", command)
        self.assertEqual(option(command, "--uvrange"), ">3klambda")
        self.assertEqual(option(command, "--calmode"), "ap")

    def test_the_bandpass_round_solves_a_gain_per_channel_bin(self):
        command = TOOL.solve_command(
            arguments("--bandpass-channels", "24"), "a.ms", {"bandpass": True}
        )

        self.assertEqual(option(command, "--bandpass"), "24")
        self.assertIn("--bandpass-amplitude", command)
        self.assertNotIn("--calmode", command)
        self.assertEqual(option(command, "--suffix"), "_solved")
        self.assertIn("--overwrite", command)

    def test_the_bandpass_round_can_leave_unsolved_data_unflagged(self):
        kept = TOOL.solve_command(
            arguments("--keep-unsolved"), "a.ms", {"bandpass": True}
        )
        usual = TOOL.solve_command(arguments(), "a.ms", {"bandpass": True})

        self.assertIn("--keep-unsolved", kept)
        self.assertNotIn("--keep-unsolved", usual)

    def test_the_number_of_sub_bands_reaches_every_model(self):
        args = arguments("--channels-out", "6")

        for model in ({"core": "phase"}, {"full": True}):
            command = TOOL.model_command(args, "a.ms", "images", model)
            self.assertEqual(option(command, "--channels-out"), "6")

    def test_the_cores_can_follow_the_fields_own_brightest_pixel(self):
        command = TOOL.model_command(
            arguments("--core-fraction", "0.13"), "a.ms", "images",
            {"core": "phase"},
        )

        self.assertEqual(option(command, "--core-model-fraction"), "0.13")
        self.assertNotIn("--core-model", command)

    def test_a_solve_that_may_flag_nothing_says_so_in_every_mode(self):
        gains = TOOL.solve_command(
            arguments("--keep-unsolved"), "a.ms", {"calmode": "ap"})
        band = TOOL.solve_command(
            arguments("--keep-unsolved"), "a.ms", {"bandpass": True})
        plain = TOOL.solve_command(arguments(), "a.ms", {"calmode": "ap"})

        self.assertIn("--keep-unsolved", gains)
        self.assertIn("--keep-unsolved", band)
        self.assertNotIn("--keep-unsolved", plain)

    def test_the_amplitudes_cores_go_twice_as_deep_as_the_phases(self):
        # By hand: 6 mJy/beam for the phases, 3 for the amplitudes.
        command = TOOL.model_command(
            arguments(), "a.ms", "images", {"core": "amplitude"})
        self.assertEqual(option(command, "--core-model"), "3")
        self.assertNotIn("--write-model", command)

        command = TOOL.model_command(
            arguments("--core-fraction", "0.13"), "a.ms", "images",
            {"core": "amplitude"},
        )
        self.assertEqual(option(command, "--core-model-fraction"), "0.065")

    def test_finished_rounds_keep_only_the_bands_images(self):
        folder = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(folder))
        names = [
            f"wsclean_robust_0_I-{part}-{kind}.fits"
            for part in ("0000", "0011", "MFS")
            for kind in ("dirty", "psf", "image", "model", "residual")
        ] + ["wsclean_robust_0_V-0000-image.fits",
             "wsclean_robust_0_V-MFS-image.fits", "wsclean_robust_0.json"]
        for name in names:
            (folder / name).write_bytes(b"")

        removed = TOOL.tidy_images(folder)

        self.assertEqual(
            sorted(path.name for path in folder.iterdir()),
            ["wsclean_robust_0.json",
             "wsclean_robust_0_I-MFS-image.fits",
             "wsclean_robust_0_I-MFS-model.fits",
             "wsclean_robust_0_I-MFS-residual.fits",
             "wsclean_robust_0_V-MFS-image.fits"],
        )
        self.assertEqual(removed, 13)

    def test_the_final_image_keeps_what_peeling_and_the_beam_read(self):
        folder = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(folder))
        for part in ("0000", "0011", "MFS"):
            for kind in ("dirty", "psf", "image", "model", "residual"):
                (folder / f"wsclean_robust_0_I-{part}-{kind}.fits").write_bytes(b"")
        (folder / "wsclean_robust_0_V-0000-image.fits").write_bytes(b"")

        TOOL.tidy_images(folder, final=True)

        kept = sorted(path.name for path in folder.iterdir())
        self.assertIn("wsclean_robust_0_I-0000-model.fits", kept)
        self.assertIn("wsclean_robust_0_I-0011-image.fits", kept)
        self.assertNotIn("wsclean_robust_0_I-0000-psf.fits", kept)
        self.assertNotIn("wsclean_robust_0_I-0011-residual.fits", kept)
        self.assertNotIn("wsclean_robust_0_V-0000-image.fits", kept)

    def test_the_solve_leaves_its_products_beside_the_data(self):
        ms, report = TOOL.solved_products("work/G1/round2.ms")

        self.assertEqual(ms, Path("work/G1/round2_solved.ms"))
        self.assertEqual(report, Path("work/G1/round2_solved.json"))


class PlanTests(unittest.TestCase):

    def test_the_plan_is_the_sequence_that_worked_by_hand(self):
        stages = [stage["id"] for stage in TOOL.DEFAULT_PLAN]

        self.assertEqual(
            stages, ["phase", "delay", "amplitude", "bandpass"],
        )
        # Phases and delays against cores; amplitudes only against the
        # whole field, on every baseline.
        self.assertEqual(TOOL.DEFAULT_PLAN[0]["model"], {"core": "phase"})
        self.assertTrue(TOOL.DEFAULT_PLAN[1]["solve"]["delays"])
        for stage in TOOL.DEFAULT_PLAN:
            if stage["solve"].get("calmode") == "ap":
                self.assertEqual(stage["model"], {"full": True})
                self.assertNotIn("uvrange", stage["solve"])

    def test_a_stage_that_may_repeat_numbers_its_rounds(self):
        self.assertEqual(TOOL.rounds_of({"id": "delay"}), ["delay"])
        self.assertEqual(
            TOOL.rounds_of({"id": "phase", "most": 3}),
            ["phase1", "phase2", "phase3"],
        )

    def test_a_round_is_settled_when_every_correction_is_small(self):
        until = {"phase_deg": 5.0, "amplitude": 0.10}

        self.assertTrue(TOOL.settled({"phase_deg": 2, "amplitude": 0.07}, until))
        self.assertFalse(TOOL.settled({"phase_deg": 2, "amplitude": 0.3}, until))
        self.assertFalse(TOOL.settled({"phase_deg": 9, "amplitude": 0.0}, until))

    def test_a_missing_number_does_not_settle_a_round(self):
        self.assertFalse(TOOL.settled({"phase_deg": 1}, {"amplitude": 0.1}))
        self.assertFalse(TOOL.settled({}, {"phase_deg": 10.0}))

    def test_a_stage_without_a_limit_is_done_once(self):
        self.assertFalse(TOOL.settled({"phase_deg": 0.0}, None))


class FindingsTests(unittest.TestCase):

    def test_the_largest_correction_is_named_with_its_antenna(self):
        found = TOOL.solve_findings(solve_report(82, delay=2.26))

        self.assertEqual(found["phase_deg"], 82.0)
        self.assertEqual(found["phase_antenna"], "S02")
        self.assertEqual(found["delay_ns"], 2.26)
        self.assertEqual(found["delay_antenna"], "S04")
        self.assertEqual(found["antennas_solved"], 2)
        self.assertTrue(found["applied"])
        self.assertNotIn("amplitude", found)

    def test_with_a_solution_per_scan_the_largest_one_counts_not_the_mean(self):
        report = solve_report(2)       # a mean of 2 degrees over the scans
        report["series"][1]["polarizations"][0].update(
            largest_phase_deg=80.0, largest_amplitude_step=0.3)

        found = TOOL.solve_findings(report)

        self.assertEqual(found["phase_deg"], 80.0)
        self.assertEqual(found["phase_antenna"], "S02")
        self.assertFalse(TOOL.settled(found, {"phase_deg": 10.0}))

    def test_amplitudes_count_only_in_an_amplitude_solve(self):
        found = TOOL.solve_findings(solve_report(1, amplitude=1.35))

        self.assertAlmostEqual(found["amplitude"], 0.35)
        self.assertEqual(found["amplitude_antenna"], "S02")
        self.assertEqual(found["amplitude_range"], [1.0, 1.35])

    def test_a_report_that_is_missing_gives_nothing_to_settle_on(self):
        found = TOOL.solve_findings(None)

        self.assertFalse(found["applied"])
        self.assertNotIn("phase_deg", found)

    def test_a_solve_per_channel_bin_is_read_by_what_runs_across_the_band(self):
        found = TOOL.solve_findings(band_report())

        self.assertEqual(found["calmode"], "bandpass")
        self.assertEqual(found["band_phase_deg"], 10.5)
        self.assertEqual(found["band_phase_antenna"], "W04")
        self.assertEqual(found["band_amplitude"], 0.326)
        self.assertEqual(found["antennas_solved"], 2)
        self.assertTrue(found["applied"])
        self.assertFalse(TOOL.solve_findings(
            band_report(applied=False))["applied"])

    def test_what_a_solve_did_not_believe_is_passed_on(self):
        report = solve_report(35, delay=0.3)
        report["rejected_delays"] = [{"antenna": "C10", "delay_ns": [17.85]}]
        report["rejected_amplitudes"] = ["C10"]

        found = TOOL.solve_findings(report)
        lines = TOOL.reading({"rounds": [{"id": "delay", "solve": found}]})

        self.assertEqual(
            found["rejected"], ["delay of C10", "amplitude of C10"]
        )
        self.assertIn(
            "not believed and not applied: delay of C10, amplitude of C10",
            lines[0],
        )
        band = band_report()
        band["rejected_bins"] = {"W04": 3}
        self.assertEqual(
            TOOL.solve_findings(band)["rejected"], ["3 bins of W04"]
        )

    def test_an_image_is_read_against_the_noise_of_its_data(self):
        found = TOOL.image_findings(IMAGE_REPORT)

        self.assertEqual(found["image_rms_Jy_beam"], 30.1e-6)
        self.assertEqual(found["data_noise_Jy_beam"], 23.5e-6)
        self.assertEqual(found["quietest_tenth_Jy_beam"], 24.9e-6)
        self.assertEqual(found["negative_beside_sources"], -0.012)
        line = TOOL.image_line("final", found)
        self.assertIn("30.1 uJy/beam", line)
        self.assertIn("24.9 in its quietest tenth", line)
        self.assertIn("23.5 for the data", line)
        self.assertIn("-1.2%", line)

    def test_an_image_that_was_not_measured_is_still_reported(self):
        line = TOOL.image_line("final", TOOL.image_findings({}))

        self.assertIn("? uJy/beam", line)

    def test_a_noise_image_under_its_new_name_is_read(self):
        found = TOOL.image_findings({
            "stokes_i": {"robust_rms_Jy_beam": 30.1e-6},
            "data_noise": {"stokes": "Q", "robust_rms_Jy_beam": 19.0e-6},
        })

        self.assertEqual(found["data_noise_Jy_beam"], 19.0e-6)

    def test_without_a_noise_image_the_line_says_nothing_of_one(self):
        # MeerKAT: XX and YY only, so no Stokes V image is made.
        found = TOOL.image_findings({
            "stokes_i": {"robust_rms_Jy_beam": 9.1e-6,
                         "grid": {"quietest_tenth_rms": 8.2e-6}},
            "noise_stokes": None,
        })
        line = TOOL.image_line("final", found)

        self.assertIsNone(found["data_noise_Jy_beam"])
        self.assertIn("final: image 9.1 uJy/beam (8.2 in its quietest tenth)",
                      line)
        self.assertNotIn("for the data", line)


class PreflightTests(unittest.TestCase):

    def setUp(self):
        self.which = mock.patch.object(
            TOOL.shutil, "which", lambda name: f"/usr/bin/{name}"
        )
        self.which.start()
        self.addCleanup(self.which.stop)
        self.capabilities = mock.patch.object(
            TOOL.IMAGER, "wsclean_capabilities",
            lambda prefix: {"version": (3, 1), "options": set()},
        )
        self.capabilities.start()
        self.addCleanup(self.capabilities.stop)

    def test_the_first_casa_that_answers_is_used(self):
        args = arguments()
        args.casa = None
        tried = []

        def answers(casa, ms):
            tried.append(casa)
            if "py3.12" in casa:
                return False, "No module named 'casatasks'"
            return True, "6.6.4"

        with mock.patch.object(TOOL, "casa_answers", answers):
            problems, found = TOOL.preflight(args, "data.ms")

        self.assertEqual(problems, [])
        self.assertEqual(tried, list(TOOL.CASA_CANDIDATES[:2]))
        self.assertEqual(args.casa, TOOL.CASA_CANDIDATES[1])
        self.assertEqual(found["casa_version"], "6.6.4")
        self.assertEqual(found["wsclean_version"], (3, 1))

    def test_a_casa_that_was_named_is_the_only_one_tried(self):
        args = arguments()
        tried = []

        def answers(casa, ms):
            tried.append((casa, ms))
            return False, "cannot open data.ms"

        with mock.patch.object(TOOL, "casa_answers", answers):
            problems, _ = TOOL.preflight(args, "data.ms")

        self.assertEqual(tried, [(CASA, "data.ms")])
        self.assertEqual(len(problems), 1)
        self.assertIn("cannot open data.ms", problems[0])
        self.assertIn("--casa", problems[0])

    def test_a_command_that_is_not_wsclean_is_a_problem(self):
        with mock.patch.object(
            TOOL.IMAGER, "wsclean_capabilities",
            lambda prefix: {"version": None, "options": set()},
        ), mock.patch.object(TOOL, "casa_answers", lambda *a: (True, "6.6")):
            problems, _ = TOOL.preflight(arguments(), "data.ms")

        self.assertEqual(len(problems), 1)
        self.assertIn("did not answer as WSClean", problems[0])

    def test_a_terminal_without_singularity_is_said_so(self):
        with mock.patch.object(TOOL.shutil, "which", lambda name: None):
            problems, _ = TOOL.preflight(arguments(), "data.ms")

        self.assertEqual(len(problems), 2)
        self.assertIn("not on PATH", problems[0])
        self.assertIn("not on PATH", problems[1])

    def test_casa_is_asked_to_open_the_measurement_set(self):
        calls = []

        def fake_run(command, **kwargs):
            calls.append(command)
            return types.SimpleNamespace(
                stdout="some banner\nCASA_OK 6.6.4.34 65286\n"
            )

        with mock.patch.object(TOOL.subprocess, "run", fake_run):
            ok, detail = TOOL.casa_answers(CASA, "data.ms")

        self.assertTrue(ok)
        self.assertEqual(detail, "6.6.4.34")
        self.assertEqual(calls[0][:4], CASA.split())
        self.assertEqual(calls[0][4], "-c")
        self.assertIn("casatasks", calls[0][5])
        self.assertIn("t.open(sys.argv[1])", calls[0][5])
        self.assertEqual(calls[0][6], "data.ms")

    def test_what_casa_said_instead_is_passed_on(self):
        def fake_run(command, **kwargs):
            return types.SimpleNamespace(
                stdout="ModuleNotFoundError: No module named 'casatasks'\n"
            )

        with mock.patch.object(TOOL.subprocess, "run", fake_run):
            ok, detail = TOOL.casa_answers(CASA, "data.ms")

        self.assertFalse(ok)
        self.assertIn("No module named 'casatasks'", detail)

    def test_a_command_that_cannot_start_is_not_an_answer(self):
        def fake_run(command, **kwargs):
            raise FileNotFoundError("singularity")

        with mock.patch.object(TOOL.subprocess, "run", fake_run):
            ok, detail = TOOL.casa_answers(CASA, "data.ms")

        self.assertFalse(ok)
        self.assertIn("singularity", detail)


class FlowTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.products = self.root / "uGMRT_Products"
        self.products.mkdir()
        self.source = self.root / "input" / "G1_scan15.ms"
        self.source.mkdir(parents=True)
        (self.source / "table.dat").write_text("visibilities")
        self.work = self.root / "work"
        self.calls = []
        # What each solve finds, in order: phases, a delay, amplitudes.
        self.solves = [
            solve_report(82), solve_report(16), solve_report(7),
            solve_report(3, delay=2.26),
            solve_report(4, amplitude=1.35),
            solve_report(2, amplitude=1.04),
            band_report(),
        ]
        self.fail_at = None
        patch = mock.patch.object(TOOL.IMAGER, "run_command", self.runner)
        patch.start()
        self.addCleanup(patch.stop)
        patch = mock.patch.object(
            TOOL, "preflight", lambda args, ms: ([], dict(FOUND))
        )
        patch.start()
        self.addCleanup(patch.stop)

    def runner(self, command, log_path):
        """Stand in for both scripts: leave what each would leave."""

        script = next(
            Path(part).name for part in command if part.endswith(".py")
        )
        ms = Path(option(command, "--ms"))
        self.assertTrue(ms.is_dir(), f"{ms} does not exist for {script}")
        self.assertTrue(Path(log_path).parent.is_dir())
        if script == "experiment_wsclean_image.py":
            kind = (
                "core" if "--core-model" in command
                else "model" if "--write-model" in command else "image"
            )
            self.calls.append((kind, ms.name))
            if self.fail_at == len(self.calls):
                raise TOOL.IMAGER.CommandFailed("wsclean exited with status 1", 1)
            folder = Path(option(command, "--workdir")) / "G1"
            folder.mkdir(parents=True, exist_ok=True)
            if kind == "core":
                (folder / "coremodel_robust_0.json").write_text(
                    json.dumps(CORE_REPORT)
                )
            else:
                (folder / "wsclean_robust_0.json").write_text(
                    json.dumps(IMAGE_REPORT)
                )
            return
        self.calls.append(("solve", ms.name))
        if self.fail_at == len(self.calls):
            raise TOOL.IMAGER.CommandFailed("python3 exited with status 1", 1)
        report = self.solves.pop(0)
        solved, report_path = TOOL.solved_products(ms)
        report_path.write_text(json.dumps(report))
        if report.get("corrected_measurement_set"):
            solved.mkdir()
            (solved / "table.dat").write_text(f"after {ms.name}")

    def run_tool(self, *extra):
        with mock.patch("builtins.print"):
            return TOOL.main([
                "--target", "G1", "--ms", str(self.source),
                "--workdir", str(self.work),
                "--outdir", str(self.products),
                "--imsize", "5120", "--cell", "1.0",
                "--channel-range", "51", "768", "--casa", CASA, *extra,
            ])

    def report(self):
        return json.loads(
            (self.work / "G1" / TOOL.REPORT_NAME).read_text()
        )

    def test_the_rounds_run_in_order_each_on_the_last_ones_result(self):
        self.run_tool()

        self.assertEqual(self.calls, [
            ("core", "round0.ms"), ("solve", "round0.ms"),
            ("core", "round1.ms"), ("solve", "round1.ms"),
            ("core", "round2.ms"), ("solve", "round2.ms"),
            ("core", "round3.ms"), ("solve", "round3.ms"),
            ("model", "round4.ms"), ("solve", "round4.ms"),
            ("model", "round5.ms"), ("solve", "round5.ms"),
            ("model", "round6.ms"), ("solve", "round6.ms"),
            ("image", "round7.ms"),
        ])
        report = self.report()
        self.assertEqual(
            [entry["id"] for entry in report["rounds"]],
            ["phase1", "phase2", "phase3", "delay", "amplitude1",
             "amplitude2", "bandpass"],
        )
        self.assertEqual(
            report["calibrated_measurement_set"],
            str(self.work / "G1" / "round7.ms"),
        )
        self.assertEqual(
            (self.work / "G1" / "round7.ms" / "table.dat").read_text(),
            "after round6.ms",
        )
        self.assertEqual(report["final"]["image_rms_Jy_beam"], 30.1e-6)
        self.assertEqual(report["versions"]["casa_version"], "6.6.4")
        self.assertIsNone(report["running"])

    def test_phases_stop_repeating_once_they_are_small(self):
        # 82 and 16 degrees call for another round; 7 does not.
        self.run_tool()

        report = self.report()
        self.assertEqual(
            [entry["solve"]["phase_deg"] for entry in report["rounds"][:3]],
            [82.0, 16.0, 7.0],
        )
        self.assertNotIn("phase4", [entry["id"] for entry in report["rounds"]])

    def test_phases_that_never_settle_stop_at_the_most_allowed(self):
        self.solves = [solve_report(40)] * 4 + self.solves[3:]

        self.run_tool()

        names = [entry["id"] for entry in self.report()["rounds"]]
        self.assertEqual(names[:5], ["phase1", "phase2", "phase3", "phase4",
                                     "delay"])

    def test_amplitudes_still_far_off_get_another_round(self):
        self.solves[5] = solve_report(2, amplitude=1.2)
        self.solves.insert(6, solve_report(1, amplitude=1.03))

        self.run_tool()

        self.assertEqual(self.calls[-7:], [
            ("model", "round5.ms"), ("solve", "round5.ms"),
            ("model", "round6.ms"), ("solve", "round6.ms"),
            ("model", "round7.ms"), ("solve", "round7.ms"),
            ("image", "round8.ms"),
        ])
        names = [entry["id"] for entry in self.report()["rounds"]]
        self.assertEqual(names[-3:], ["amplitude2", "amplitude3", "bandpass"])

    def test_the_measurement_set_given_is_only_copied(self):
        self.run_tool()

        self.assertEqual(
            sorted(path.name for path in self.source.iterdir()),
            ["table.dat"],
        )
        self.assertEqual(
            sorted(path.name for path in self.source.parent.iterdir()),
            ["G1_scan15.ms"],
        )
        self.assertEqual(list(self.products.iterdir()), [])

    def test_each_rounds_copy_is_removed_once_it_has_served(self):
        self.run_tool()

        kept = sorted(
            path.name for path in (self.work / "G1").glob("round*.ms")
        )
        self.assertEqual(kept, ["round7.ms"])
        # What each solve found stays.
        self.assertTrue((self.work / "G1" / "round0_solved.json").is_file())

    def test_the_copies_can_be_kept(self):
        self.run_tool("--keep-intermediate")

        kept = sorted(
            path.name for path in (self.work / "G1").glob("round*.ms")
        )
        self.assertEqual(kept, [f"round{n}.ms" for n in range(8)])

    def test_the_reading_gives_each_round_and_the_image(self):
        self.run_tool()

        lines = self.report()["reading"]
        self.assertIn("phase1: lopsided part of the dirty image 15.0% of the "
                      "peaks; largest phase 82 deg (S02).", lines)
        self.assertIn("largest delay 2.26 ns (S04)", " ".join(lines))
        self.assertIn("amplitudes 1.00 to 1.35", " ".join(lines))
        self.assertIn(
            "bandpass: across the band, phases up to 10.5 deg rms (W04); "
            "amplitudes up to 32.6% rms (W04).", lines,
        )
        self.assertTrue(lines[-1].startswith("final: image 30.1 uJy/beam"))
        # The model of the whole field is an image too, and is reported.
        self.assertTrue(
            any(line.startswith("amplitude1 model: image") for line in lines)
        )

    def test_a_failed_step_stops_the_run_and_says_how_to_continue(self):
        self.fail_at = 4

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("phase2", str(stopped.exception))
        self.assertIn("continues from this step", str(stopped.exception))
        report = self.report()
        self.assertEqual(
            report["done"],
            ["copy", "phase1:model", "phase1:solve", "phase2:model"],
        )
        self.assertEqual(report["failed"]["step"], "phase2: solve on round1.ms")

    def test_the_same_command_continues_after_the_last_finished_step(self):
        self.fail_at = 4
        with self.assertRaises(SystemExit):
            self.run_tool()
        before = list(self.calls)
        self.fail_at = None
        (self.source / "table.dat").write_text("changed since")

        self.run_tool()

        # The failed solve is the first thing run again; the copy and the
        # three steps before it are not.
        self.assertEqual(before[-1], ("solve", "round1.ms"))
        self.assertEqual(self.calls[len(before)], ("solve", "round1.ms"))
        self.assertEqual(self.calls[len(before) + 1], ("core", "round2.ms"))
        self.assertEqual(self.calls[-1], ("image", "round7.ms"))
        self.assertNotIn("failed", self.report())

    def test_a_finished_run_given_again_runs_nothing(self):
        self.run_tool()
        count = len(self.calls)

        self.run_tool()

        self.assertEqual(len(self.calls), count)

    def test_a_solve_that_applies_nothing_stops_the_run(self):
        self.solves[1] = solve_report(0, applied=False)

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("applied nothing", str(stopped.exception))
        self.assertEqual(self.calls[-1], ("solve", "round1.ms"))
        self.assertIn("NOTHING APPLIED", " ".join(self.report()["reading"]))

    def test_nothing_is_started_when_a_tool_does_not_answer(self):
        with mock.patch.object(
            TOOL, "preflight", lambda args, ms: (["CASA: no command"], {})
        ), self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("Cannot start", str(stopped.exception))
        self.assertEqual(self.calls, [])
        self.assertFalse((self.work / "G1" / "round0.ms").exists())

    def test_a_work_folder_inside_the_products_is_refused(self):
        self.work = self.products / "work"

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("inside the pipeline's products", str(stopped.exception))
        self.assertEqual(self.calls, [])

    def test_a_work_folder_of_another_run_is_not_continued(self):
        self.run_tool()
        other = self.root / "input" / "other.ms"
        other.mkdir()
        self.source = other

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("holds a run on", str(stopped.exception))

    def test_restart_discards_an_earlier_run_of_this_script(self):
        self.run_tool()
        count = len(self.calls)
        self.solves = [solve_report(3), solve_report(1, delay=0.1),
                       solve_report(1, amplitude=1.02), band_report()]

        self.run_tool("--restart")

        self.assertEqual(self.calls[count], ("core", "round0.ms"))
        self.assertEqual(self.calls[-1], ("image", "round4.ms"))

    def test_restart_removes_nothing_it_did_not_make(self):
        folder = self.work / "G1"
        folder.mkdir(parents=True)
        (folder / "notes.txt").write_text("not a run of this script")

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool("--restart")

        self.assertIn("will not remove it", str(stopped.exception))
        self.assertTrue((folder / "notes.txt").is_file())

    def test_a_second_copy_of_the_run_is_refused(self):
        import os
        import socket

        folder = self.work / "G1"
        folder.mkdir(parents=True)
        lock = folder / TOOL.IMAGER.LOCK_NAME
        lock.write_text(json.dumps({
            "pid": os.getpid(), "host": socket.gethostname(),
            "since": "2026-10-08 20:03:00 UTC",
        }))

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("in use by another run", str(stopped.exception))
        self.assertEqual(self.calls, [])
        self.assertFalse((folder / "round0.ms").exists())

    def test_restart_does_not_clear_a_folder_another_run_holds(self):
        import os
        import socket

        self.run_tool()
        folder = self.work / "G1"
        lock = folder / TOOL.IMAGER.LOCK_NAME
        lock.write_text(json.dumps(
            {"pid": os.getpid(), "host": socket.gethostname()}
        ))

        with self.assertRaises(SystemExit):
            self.run_tool("--restart")

        self.assertTrue((folder / "round7.ms").is_dir())
        self.assertTrue((folder / TOOL.REPORT_NAME).is_file())

    def test_the_folder_is_released_when_the_run_ends_or_fails(self):
        self.fail_at = 2
        with self.assertRaises(SystemExit):
            self.run_tool()
        self.assertFalse(
            (self.work / "G1" / TOOL.IMAGER.LOCK_NAME).exists()
        )
        self.fail_at = None

        self.run_tool()

        self.assertFalse(
            (self.work / "G1" / TOOL.IMAGER.LOCK_NAME).exists()
        )
        self.assertEqual(self.calls[-1], ("image", "round7.ms"))

    def test_the_data_can_be_imaged_before_any_round(self):
        self.run_tool("--image-first")

        self.assertEqual(self.calls[0], ("image", "round0.ms"))
        self.assertEqual(self.calls[1], ("core", "round0.ms"))
        report = self.report()
        self.assertEqual(report["before"]["image_rms_Jy_beam"], 30.1e-6)
        self.assertTrue(report["reading"][0].startswith("before: image"))
        self.assertIn("before:image", report["done"])

    def test_a_round_within_the_scan_can_be_added_at_the_end(self):
        within = solve_report(1, amplitude=1.0)
        within["series"][1]["solutions"] = 11
        within["series"][1]["polarizations"][0].update(
            amplitude_rms=0.042, phase_rms_deg=3.5
        )
        self.solves.append(within)

        self.run_tool("--time-solint", "2min")

        self.assertEqual(self.calls[-3:], [
            ("model", "round7.ms"), ("solve", "round7.ms"),
            ("image", "round8.ms"),
        ])
        report = self.report()
        self.assertEqual(report["rounds"][-1]["id"], "time")
        self.assertEqual(report["rounds"][-1]["solve"]["time_amplitude"], 0.042)
        self.assertIn(
            "within the scan, amplitudes move by up to 4.2% rms (S02) and "
            "phases by up to 3.5 deg", " ".join(report["reading"]),
        )

    def test_the_round_within_the_scan_solves_at_the_interval_given(self):
        command = TOOL.solve_command(
            arguments(), "a.ms", {"calmode": "ap", "solint": "2min"}
        )

        self.assertEqual(option(command, "--solint"), "2min")
        self.assertEqual(option(command, "--calmode"), "ap")

    def test_a_faint_source_gets_only_phases_and_amplitudes(self):
        self.solves = [solve_report(3), solve_report(1, amplitude=1.02)]

        self.run_tool("--phase-amplitude-only")

        self.assertEqual(self.calls, [
            ("core", "round0.ms"), ("solve", "round0.ms"),
            ("model", "round1.ms"), ("solve", "round1.ms"),
            ("image", "round2.ms"),
        ])
        self.assertEqual(
            [entry["id"] for entry in self.report()["rounds"]],
            ["phase1", "amplitude1"],
        )

    def test_only_the_final_image_takes_the_larger_size(self):
        seen = []
        real = self.runner

        def watching(command, log_path):
            if "experiment_wsclean_image.py" in " ".join(command):
                seen.append(option(command, "--imsize"))
            return real(command, log_path)

        with mock.patch.object(TOOL.IMAGER, "run_command", watching):
            self.run_tool("--final-imsize", "8192")

        self.assertEqual(seen[-1], "8192")
        self.assertEqual(set(seen[:-1]), {"5120"})

    def test_phases_within_the_scans_follow_the_phases_per_scan(self):
        self.solves = [
            solve_report(7),                       # per scan: small at once
            solve_report(40), solve_report(3),     # within the scans
            solve_report(3, delay=0.3),
            solve_report(2, amplitude=1.04), band_report(),
        ]

        self.run_tool("--time-phase-solint", "1min")

        names = [entry["id"] for entry in self.report()["rounds"]]
        self.assertEqual(names, ["phase1", "phase_time1", "phase_time2",
                                 "delay", "amplitude1", "bandpass"])
        # Cores-only models for them, as for the phases per scan.
        self.assertEqual(self.calls[2:6], [
            ("core", "round1.ms"), ("solve", "round1.ms"),
            ("core", "round2.ms"), ("solve", "round2.ms"),
        ])

    def test_amplitudes_against_cores_come_before_the_whole_fields(self):
        self.solves = [
            solve_report(7), solve_report(3, delay=0.3),
            # Against the cores: far from one, then near it.
            solve_report(4, amplitude=1.35), solve_report(2, amplitude=1.05),
            solve_report(2, amplitude=1.02), band_report(),
        ]
        commands = []
        real = self.runner

        def watching(command, log_path):
            commands.append(command)
            return real(command, log_path)

        with mock.patch.object(TOOL.IMAGER, "run_command", watching):
            self.run_tool("--amplitude-cores")

        names = [entry["id"] for entry in self.report()["rounds"]]
        self.assertEqual(names, [
            "phase1", "delay", "amplitude_cores1", "amplitude_cores2",
            "amplitude1", "bandpass",
        ])
        self.assertEqual(self.calls[4:10], [
            ("core", "round2.ms"), ("solve", "round2.ms"),
            ("core", "round3.ms"), ("solve", "round3.ms"),
            ("model", "round4.ms"), ("solve", "round4.ms"),
        ])
        solves = [command for command in commands
                  if "experiment_phase_offsets.py" in " ".join(command)]
        self.assertEqual(option(solves[2], "--calmode"), "ap")
        self.assertEqual(option(solves[2], "--uvrange"), ">3klambda")
        self.assertEqual(option(solves[2], "--gaintype"), "G")
        # The whole field's rounds stay on every baseline.
        self.assertEqual(option(solves[4], "--calmode"), "ap")
        self.assertNotIn("--uvrange", solves[4])
        cores = [command for command in commands if "--core-model" in command]
        self.assertEqual(
            [option(command, "--core-model") for command in cores],
            ["6", "6", "3", "3"],
        )
        self.assertIn("amplitude_cores1: ", " ".join(self.report()["reading"]))

    def test_the_baselines_of_the_cores_amplitudes_can_be_given(self):
        self.solves = [
            solve_report(7), solve_report(3, delay=0.3),
            solve_report(2, amplitude=1.05),
            solve_report(2, amplitude=1.02), band_report(),
        ]
        commands = []
        real = self.runner

        def watching(command, log_path):
            commands.append(command)
            return real(command, log_path)

        with mock.patch.object(TOOL.IMAGER, "run_command", watching):
            self.run_tool("--amplitude-cores",
                          "--amplitude-cores-uvrange", ">5klambda")

        ranges = [option(command, "--uvrange") for command in commands
                  if "--uvrange" in command]
        self.assertEqual(ranges, [">5klambda"])

    def test_without_asking_no_amplitude_is_solved_against_cores(self):
        self.run_tool()

        names = [entry["id"] for entry in self.report()["rounds"]]
        self.assertFalse([name for name in names if "cores" in name])

    def test_calibrated_data_can_be_solved_against_the_whole_field_only(self):
        self.solves = [
            solve_report(6, amplitude=1.03), solve_report(2, amplitude=1.01),
            band_report(),
        ]

        self.run_tool("--whole-field-only")

        names = [entry["id"] for entry in self.report()["rounds"]]
        self.assertEqual(names, ["amplitude1", "amplitude2", "bandpass"])
        # The first thing made is an image of the data as they came.
        self.assertEqual(self.calls, [
            ("model", "round0.ms"), ("solve", "round0.ms"),
            ("model", "round1.ms"), ("solve", "round1.ms"),
            ("model", "round2.ms"), ("solve", "round2.ms"),
            ("image", "round3.ms"),
        ])
        self.assertFalse([call for call in self.calls if call[0] == "core"])

    def test_rounds_against_cores_cannot_be_asked_for_with_it(self):
        for extra in (("--time-phase-solint", "1min"), ("--amplitude-cores",)):
            with self.assertRaises(SystemExit), \
                    mock.patch("sys.stderr"):
                self.run_tool("--whole-field-only", *extra)
        self.assertEqual(self.calls, [])

    def test_tidying_removes_the_sub_band_images_of_finished_rounds(self):
        seen = []
        with mock.patch.object(
            TOOL, "tidy_images",
            lambda folder, final=False: seen.append(
                (Path(folder).parent.name, final))
        ):
            self.run_tool("--tidy")

        self.assertIn(("phase1", False), seen)
        self.assertIn(("bandpass", False), seen)
        self.assertIn(("final", True), seen)

    def test_without_asking_no_image_is_made_first(self):
        self.run_tool()

        self.assertEqual(self.calls[0], ("core", "round0.ms"))
        self.assertNotIn("before", self.report())

    def test_a_dry_run_prints_the_steps_and_starts_nothing(self):
        printed = []
        with mock.patch("builtins.print", lambda *a, **k: printed.append(a[0])):
            TOOL.main([
                "--target", "G1", "--ms", str(self.source),
                "--workdir", str(self.work),
                "--imsize", "5120", "--cell", "1.0", "--dry-run",
            ])

        text = "\n".join(printed)
        self.assertEqual(self.calls, [])
        self.assertFalse(self.work.exists())
        self.assertIn("--core-model 6", text)
        self.assertNotIn("--uvrange", text)
        self.assertIn("--delays", text)
        self.assertIn("--write-model", text)
        self.assertIn("--bandpass 30 --bandpass-amplitude", text)
        self.assertNotIn("--no-negative", text)
        self.assertIn("up to 4 times, until phase_deg <= 10", text)
        self.assertIn("casa-modular-py3.12-latest.sif", text)

    def test_a_missing_measurement_set_is_said_at_once(self):
        self.source = self.root / "input" / "absent.ms"

        with self.assertRaises(SystemExit) as stopped:
            self.run_tool()

        self.assertIn("No Measurement Set", str(stopped.exception))


if __name__ == "__main__":
    unittest.main()
