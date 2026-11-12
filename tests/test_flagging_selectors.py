import importlib.util
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, SOURCE_DIR / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SELECTORS = load_module(
    "test_radio_continuum_flagging_selectors",
    "flagging_selectors.py",
)


def spw(*frequencies):
    return {
        "nchan": len(frequencies),
        "channel_frequencies_MHz": list(frequencies),
        "fmin_MHz": min(frequencies),
        "fmax_MHz": max(frequencies),
    }


def load_flagging_module_without_casa():
    casatasks = types.ModuleType("casatasks")

    def unused_task(**kwargs):
        return None

    casatasks.flagdata = unused_task
    casatasks.flagmanager = unused_task

    casatools = types.ModuleType("casatools")

    class UnusedTable:
        pass

    casatools.table = UnusedTable
    numpy = types.ModuleType("numpy")
    matplotlib = types.ModuleType("matplotlib")
    pyplot = types.ModuleType("matplotlib.pyplot")

    def unused_use(*args, **kwargs):
        return None

    matplotlib.use = unused_use
    matplotlib.pyplot = pyplot
    path = SOURCE_DIR / "flagging.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_frequency_flagging_integration",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    with mock.patch.dict(
        sys.modules,
        {
            "casatasks": casatasks,
            "casatools": casatools,
            "numpy": numpy,
            "matplotlib": matplotlib,
            "matplotlib.pyplot": pyplot,
        },
    ):
        spec.loader.exec_module(module)
    return module


FLAGGING = load_flagging_module_without_casa()


class FrequencySelectorTests(unittest.TestCase):

    def test_science_mask_selects_only_outside_channels_per_spw(self):
        spectral_windows = {
            "0": spw(544.0, 579.9, 580.0, 700.0, 1015.0, 1015.1, 1088.0),
            "1": spw(600.0, 700.0, 900.0, 1000.0),
            # Descending channel order must be preserved in the indices.
            "2": spw(1088.0, 1015.1, 1015.0, 580.0, 579.9, 544.0),
        }

        report = SELECTORS.build_science_range_channel_masks(
            spectral_windows,
            580.0,
            1015.0,
        )

        self.assertEqual(
            report["spw_selection"],
            "0:0~1;5~6,2:0~1;4~5",
        )
        self.assertNotIn("1:", report["spw_selection"])
        self.assertEqual(report["science_range_MHz"], [580.0, 1015.0])

    def test_static_masks_are_intersected_with_real_spw_channels(self):
        spectral_windows = {
            "0": spw(900.0, 932.0, 933.0, 940.0, 960.0, 961.0),
            "1": spw(1160.0, 1163.0, 1200.0, 1299.0, 1300.0),
            "2": spw(1400.0, 1450.0, 1500.0),
        }

        report = SELECTORS.build_intersected_static_masks(
            spectral_windows,
            ["933~960MHz", "1163~1299MHz", "1524~1630MHz"],
        )

        self.assertEqual(report["spw_selection"], "0:2~4,1:1~3")
        self.assertNotIn("*:", report["spw_selection"])
        self.assertEqual(len(report["applied_intersections"]), 2)
        self.assertTrue(report["skipped_non_overlaps"])
        self.assertTrue(all(
            item["reason"] in {
                "no_spw_overlap",
                "no_channel_centres_in_overlap",
            }
            for item in report["skipped_non_overlaps"]
        ))

    def test_nonoverlapping_static_mask_is_reported_without_selection(self):
        report = SELECTORS.build_intersected_static_masks(
            {"4": spw(1300.0, 1350.0, 1400.0)},
            ["933~960MHz"],
        )

        self.assertEqual(report["spw_selection"], "")
        self.assertEqual(report["status"], "no_overlapping_channels")
        self.assertEqual(
            report["skipped_non_overlaps"][0]["reason"],
            "no_spw_overlap",
        )

    def test_frequency_masking_requires_channel_order_metadata(self):
        with self.assertRaisesRegex(ValueError, "Rerun the inspection"):
            SELECTORS.build_science_range_channel_masks(
                {"0": {"nchan": 4, "fmin_MHz": 544, "fmax_MHz": 1088}},
                580,
                1015,
            )


class FlaggingFrequencyIntegrationTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def manager(self, telescope, band, spectral_windows):
        manager = object.__new__(FLAGGING.CASAFlagger)
        manager.ms = str(self.root / "input.ms")
        manager.root = self.root
        manager.telescope = telescope
        manager.band = band
        manager.inspection = {"spws": spectral_windows}
        manager.edge_fraction = 0.10
        manager.known_frequency_masks = []
        return manager

    def test_meerkat_uhf_edge_step_uses_guaranteed_range_selector(self):
        manager = self.manager(
            "meerkat",
            "MEERKAT_UHF",
            {"0": spw(544.0, 579.0, 580.0, 1015.0, 1016.0, 1088.0)},
        )
        calls = []

        with mock.patch.object(
            FLAGGING,
            "flagdata",
            lambda **kwargs: calls.append(kwargs),
        ):
            report = manager.flag_edge_channels()

        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["spw"], "0:0~1;4~5")
        self.assertEqual(
            report["scheme"],
            "outside_guaranteed_science_range",
        )
        saved = json.loads((self.root / "edge_flagging_report.json").read_text())
        self.assertEqual(saved["telescope"], "meerkat")
        self.assertEqual(saved["science_range_MHz"], [580.0, 1015.0])

    def test_ugmrt_b4_wide_backend_masks_channels_above_receiver(self):
        manager = self.manager(
            "ugmrt",
            "B4",
            {"0": spw(550.0, 700.0, 850.0, 850.1, 900.0, 950.0)},
        )
        calls = []

        with mock.patch.object(
            FLAGGING,
            "flagdata",
            lambda **kwargs: calls.append(kwargs),
        ):
            report = manager.flag_edge_channels()

        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["spw"], "0:3~5")
        self.assertEqual(report["science_range_MHz"], [550.0, 850.0])
        self.assertEqual(
            report["digitized_receiver_range_MHz"],
            [550.0, 950.0],
        )

    def test_static_meerkat_masks_never_leak_to_another_telescope(self):
        manager = self.manager(
            "ugmrt",
            "B5",
            {"0": spw(900.0, 933.0, 940.0, 960.0)},
        )
        manager.known_frequency_masks = ["933~960MHz"]
        calls = []

        with mock.patch.object(
            FLAGGING,
            "flagdata",
            lambda **kwargs: calls.append(kwargs),
        ):
            report = manager.flag_known_frequency_ranges()

        self.assertEqual(calls, [])
        self.assertEqual(report["status"], "skipped_telescope_band_guard")

    def test_meerkat_l_static_step_applies_only_intersected_spws(self):
        manager = self.manager(
            "meerkat",
            "MEERKAT_L",
            {
                "0": spw(900.0, 932.0, 933.0, 940.0, 960.0, 961.0),
                "1": spw(1400.0, 1450.0, 1500.0),
            },
        )
        manager.known_frequency_masks = [
            "933~960MHz",
            "1524~1630MHz",
        ]
        calls = []

        with mock.patch.object(
            FLAGGING,
            "flagdata",
            lambda **kwargs: calls.append(kwargs),
        ):
            report = manager.flag_known_frequency_ranges()

        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["spw"], "0:2~4")
        self.assertTrue(report["skipped_non_overlaps"])
        saved = json.loads(
            (self.root / "known_frequency_mask_report.json").read_text()
        )
        self.assertEqual(saved["status"], "applied")

    def test_completed_flagging_refreshes_only_orchestrator_provenance(self):
        manager = object.__new__(FLAGGING.CASAFlagger)
        manager.root = self.root
        manager.donefile = self.root / "flagging.done"
        manager.donefile.touch()
        manager.ms = str(self.root / "input.ms")
        manager.telescope = "ugmrt"
        manager.band = "B4"
        manager.provenance_fingerprint = "b" * 64
        (self.root / "postflag_summary.json").write_text(json.dumps({
            "measurement_set": str(Path(manager.ms).resolve()),
            "telescope": "ugmrt",
            "band": "B4",
            "provenance_fingerprint": "a" * 64,
            "flagged": 123,
        }))
        qa_calls = []
        manager.regenerate_qa_plots = lambda: qa_calls.append(True)

        manager.run()

        report = json.loads(
            (self.root / "postflag_summary.json").read_text()
        )
        self.assertEqual(report["provenance_fingerprint"], "b" * 64)
        self.assertEqual(report["flagged"], 123)
        self.assertEqual(len(report["provenance_migration_history"]), 1)
        self.assertEqual(qa_calls, [True])

    def test_completed_flagging_never_relabels_another_ms(self):
        manager = object.__new__(FLAGGING.CASAFlagger)
        manager.root = self.root
        manager.donefile = self.root / "flagging.done"
        manager.donefile.touch()
        manager.ms = str(self.root / "current.ms")
        manager.telescope = "ugmrt"
        manager.band = "B4"
        manager.provenance_fingerprint = "b" * 64
        (self.root / "postflag_summary.json").write_text(json.dumps({
            "measurement_set": str(self.root / "different.ms"),
            "telescope": "ugmrt",
            "band": "B4",
            "provenance_fingerprint": "a" * 64,
        }))
        manager.regenerate_qa_plots = lambda: None

        with self.assertRaisesRegex(RuntimeError, "different Measurement Set"):
            manager.run()



GIB = 1024 ** 3


def chunk_group(nant=64, nchan=32768, ncorr=4, integration_s=8.0,
                n_integrations=150, scan=1, spw=0):
    return {
        "scan": scan,
        "spw": spw,
        "rows_per_integration": nant * (nant + 1) // 2,
        "integration_s": integration_s,
        "n_integrations": n_integrations,
        "nchan": nchan,
        "ncorr": ncorr,
    }


class FlagdataChunkPlanningTests(unittest.TestCase):

    def test_whole_scan_kept_when_it_fits(self):
        plan = SELECTORS.plan_flagdata_time_chunk(
            [chunk_group(nchan=4096)],
            232 * GIB,
        )
        self.assertEqual(plan["status"], "scan_fits")
        self.assertEqual(plan["ntime"], "scan")

    def test_bounded_chunk_fits_every_constrained_group(self):
        groups = [
            chunk_group(nchan=32768, integration_s=8.0, scan=1),
            chunk_group(nchan=32768, integration_s=4.0,
                        n_integrations=300, scan=2),
            chunk_group(nchan=4096, scan=3),
        ]
        limit = 232 * GIB
        plan = SELECTORS.plan_flagdata_time_chunk(groups, limit)

        self.assertEqual(plan["status"], "bounded")
        self.assertRegex(plan["ntime"], r"^\d+s$")
        seconds = plan["chunk_s"]
        budget = limit * SELECTORS.FLAGDATA_MEMORY_FRACTION

        for group in groups:
            integrations = min(
                group["n_integrations"],
                seconds / group["integration_s"],
            )
            loaded = (
                group["rows_per_integration"] * group["nchan"] *
                group["ncorr"] * SELECTORS.FLAGDATA_BYTES_PER_ELEMENT *
                integrations
            )
            self.assertLessEqual(loaded, budget)
            self.assertGreaterEqual(
                integrations,
                min(group["n_integrations"],
                    SELECTORS.FLAGDATA_MINIMUM_INTEGRATIONS),
            )

    def test_insufficient_memory_is_reported_not_guessed(self):
        plan = SELECTORS.plan_flagdata_time_chunk(
            [chunk_group()],
            1 * GIB,
        )
        self.assertEqual(plan["status"], "insufficient_memory")
        self.assertEqual(plan["ntime"], "scan")
        self.assertGreater(plan["required_memory_bytes"], 1 * GIB)

    def test_unknown_limit_keeps_whole_scan(self):
        plan = SELECTORS.plan_flagdata_time_chunk([chunk_group()], None)
        self.assertEqual(plan["status"], "no_memory_limit_detected")
        self.assertEqual(plan["ntime"], "scan")

    def test_invalid_fraction_is_rejected(self):
        with self.assertRaises(ValueError):
            SELECTORS.plan_flagdata_time_chunk(
                [chunk_group()], 232 * GIB, memory_fraction=1.5,
            )


class MemoryLimitDetectionTests(unittest.TestCase):

    @staticmethod
    def reader(files):
        return lambda path: files.get(path)

    def test_smallest_source_wins(self):
        files = {
            "/proc/self/cgroup": "0::/job_1/step_0\n",
            "/sys/fs/cgroup/job_1/step_0/memory.max": "max\n",
            "/sys/fs/cgroup/job_1/memory.max": str(300 * GIB),
        }
        result = SELECTORS.detect_memory_limit_bytes(
            environ={"SLURM_MEM_PER_NODE": str(232 * 1024)},
            read_text=self.reader(files),
            physical_bytes=500 * GIB,
        )
        self.assertEqual(result["limit_bytes"], 232 * GIB)
        self.assertEqual(result["limiting_source"], "slurm_mem_per_node")
        self.assertEqual(result["sources_bytes"]["cgroup"], 300 * GIB)

    def test_tighter_cgroup_ancestor_is_used(self):
        files = {
            "/proc/self/cgroup": "0::/job_1/step_0\n",
            "/sys/fs/cgroup/job_1/memory.max": str(100 * GIB),
        }
        result = SELECTORS.detect_memory_limit_bytes(
            environ={"SLURM_MEM_PER_NODE": str(232 * 1024)},
            read_text=self.reader(files),
            physical_bytes=500 * GIB,
        )
        self.assertEqual(result["limiting_source"], "cgroup")
        self.assertEqual(result["limit_bytes"], 100 * GIB)

    def test_cgroup_v1_unlimited_sentinel_is_ignored(self):
        files = {
            "/proc/self/cgroup": "5:memory:/slurm/job_1\n",
            "/sys/fs/cgroup/memory/slurm/job_1/memory.limit_in_bytes":
                "9223372036854771712\n",
        }
        result = SELECTORS.detect_memory_limit_bytes(
            environ={},
            read_text=self.reader(files),
            physical_bytes=64 * GIB,
        )
        self.assertNotIn("cgroup", result["sources_bytes"])
        self.assertEqual(result["limit_bytes"], 64 * GIB)

    def test_configured_limit_can_only_tighten(self):
        result = SELECTORS.detect_memory_limit_bytes(
            environ={"SLURM_MEM_PER_NODE": str(232 * 1024)},
            read_text=self.reader({}),
            physical_bytes=500 * GIB,
            configured_gb=900,
        )
        self.assertEqual(result["limit_bytes"], 232 * GIB)


class FakeMetadata:

    def __init__(self, scans):
        self.scans = scans

    def open(self, path):
        pass

    def close(self):
        pass

    def scannumbers(self):
        return list(self.scans)

    def timesforscan(self, scan):
        return self.scans[scan]["times"]

    def antennasforscan(self, scan):
        return list(range(self.scans[scan]["nant"]))

    def spwsforscan(self, scan):
        return self.scans[scan]["spws"]

    def nchan(self, spw):
        return {0: 4096, 1: 32768}[spw]

    def exposuretime(self, scan, spwid):
        return {"value": 7.5, "unit": "s"}


class FakePolarizationTable:

    def open(self, path):
        pass

    def getcol(self, name):
        return [2, 4]

    def close(self):
        pass


class ScanChunkSummaryTests(unittest.TestCase):

    def test_geometry_from_metadata(self):
        scans = {
            1: {"times": [0.0, 8.0, 16.0, 16.0, 24.0], "nant": 64,
                "spws": [0, 1]},
            2: {"times": [100.0], "nant": 30, "spws": [0]},
        }
        groups = SELECTORS.summarize_scan_chunks(
            "fake.ms",
            lambda: FakeMetadata(scans),
            FakePolarizationTable,
        )
        by_key = {(item["scan"], item["spw"]): item for item in groups}

        wide = by_key[(1, 1)]
        self.assertEqual(wide["rows_per_integration"], 64 * 65 // 2)
        self.assertEqual(wide["n_integrations"], 4)
        self.assertEqual(wide["integration_s"], 8.0)
        self.assertEqual(wide["nchan"], 32768)
        self.assertEqual(wide["ncorr"], 4)

        single = by_key[(2, 0)]
        self.assertEqual(single["n_integrations"], 1)
        self.assertEqual(single["integration_s"], 7.5)


class FlaggerChunkPlanTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def manager(self):
        manager = object.__new__(FLAGGING.CASAFlagger)
        manager.ms = str(self.root / "input.ms")
        manager.root = self.root
        manager.telescope = "meerkat"
        manager.band = "MEERKAT_UHF"
        manager.inspection = {
            "spws": {"0": spw(544.0, 579.0, 580.0, 1015.0, 1016.0, 1088.0)}
        }
        manager.edge_fraction = 0.10
        manager.known_frequency_masks = []
        return manager

    def patched(self, limit_bytes):
        casatools = types.ModuleType("casatools")
        casatools.msmetadata = lambda: FakeMetadata({
            1: {"times": [8.0 * index for index in range(150)],
                "nant": 64, "spws": [1]},
        })
        limit = {
            "limit_bytes": limit_bytes,
            "limiting_source": "slurm_mem_per_node",
            "sources_bytes": {"slurm_mem_per_node": limit_bytes},
        }
        return [
            mock.patch.dict(sys.modules, {"casatools": casatools}),
            mock.patch.object(FLAGGING, "table", FakePolarizationTable),
            mock.patch.object(
                FLAGGING.FLAGGING_SELECTORS,
                "detect_memory_limit_bytes",
                lambda **kwargs: limit,
            ),
        ]

    def test_bounded_chunk_reaches_flagdata_and_is_recorded(self):
        manager = self.manager()
        calls = []
        patches = self.patched(232 * GIB)
        for item in patches:
            item.start()
        try:
            with mock.patch.object(
                FLAGGING, "flagdata", lambda **kwargs: calls.append(kwargs),
            ):
                manager.flag_edge_channels()
        finally:
            for item in reversed(patches):
                item.stop()

        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["ntime"], "304s")
        saved = json.loads(
            (self.root / "flagdata_chunk_plan.json").read_text()
        )
        self.assertEqual(saved["status"], "bounded")

    def test_insufficient_memory_fails_before_any_flag_is_written(self):
        manager = self.manager()
        calls = []
        patches = self.patched(1 * GIB)
        for item in patches:
            item.start()
        try:
            with mock.patch.object(
                FLAGGING, "flagdata", lambda **kwargs: calls.append(kwargs),
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "No flags were written",
                ):
                    manager.flag_edge_channels()
        finally:
            for item in reversed(patches):
                item.stop()

        self.assertEqual(calls, [])


if __name__ == "__main__":
    unittest.main()
