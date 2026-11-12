import importlib.util
import unittest
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_profiles_module():
    path = SOURCE_DIR / "telescope_profiles.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_telescope_profiles",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PROFILES = load_profiles_module()


def load_config_module():
    path = SOURCE_DIR / "config.py"
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_band_keyed_config",
        path,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONFIG = load_config_module()


class TelescopeIdentityTests(unittest.TestCase):

    def test_telescope_aliases_normalize_to_canonical_names(self):
        cases = {
            "GMRT": "ugmrt",
            "uGMRT": "ugmrt",
            "Giant Metrewave Radio Telescope": "ugmrt",
            "MeerKAT": "meerkat",
            "mkat": "meerkat",
            "LOFAR": "lofar",
            "ILT": "lofar",
            "International LOFAR Telescope": "lofar",
        }

        for value, expected in cases.items():
            with self.subTest(value=value):
                self.assertEqual(
                    PROFILES.normalize_telescope(value),
                    expected,
                )

    def test_auto_resolution_prefers_observation_metadata(self):
        result = PROFILES.resolve_telescope(
            "auto",
            observatory_names=["MeerKAT"],
            antenna_names=["m000", "m001"],
        )

        self.assertEqual(result["name"], "meerkat")
        self.assertEqual(result["method"], "observation_metadata")
        self.assertEqual(result["requested"], "auto")

    def test_antenna_names_are_a_supported_fallback(self):
        cases = [
            (["C00", "C05", "W01"], "ugmrt"),
            (["m000", "m031", "m063"], "meerkat"),
            (["CS001HBA0", "CS002HBA1", "DE601HBA"], "lofar"),
        ]

        for antennas, expected in cases:
            with self.subTest(antennas=antennas):
                result = PROFILES.resolve_telescope(
                    "auto",
                    antenna_names=antennas,
                )
                self.assertEqual(result["name"], expected)
                self.assertEqual(result["method"], "antenna_names")

        with self.assertRaisesRegex(ValueError, "Could not identify"):
            PROFILES.resolve_telescope(
                "auto",
                antenna_names=["m001", "UNKNOWN_A", "UNKNOWN_B"],
            )

    def test_explicit_telescope_is_verified_against_metadata(self):
        verified = PROFILES.resolve_telescope(
            "lofar",
            observatory_names=["LOFAR"],
        )
        self.assertEqual(verified["method"], "explicit_verified")

        with self.assertRaisesRegex(ValueError, "conflicts"):
            PROFILES.resolve_telescope(
                "ugmrt",
                observatory_names=["MeerKAT"],
            )

        with self.assertRaisesRegex(ValueError, "antenna names identify"):
            PROFILES.resolve_telescope(
                "auto",
                observatory_names=["MeerKAT"],
                antenna_names=["CS001HBA0", "CS002HBA1"],
            )

    def test_auto_resolution_fails_closed_for_unknown_or_mixed_metadata(self):
        with self.assertRaisesRegex(ValueError, "Could not identify"):
            PROFILES.resolve_telescope(
                "auto",
                observatory_names=["UNKNOWN"],
            )

        with self.assertRaisesRegex(ValueError, "multiple supported telescopes"):
            PROFILES.resolve_telescope(
                "auto",
                observatory_names=["MeerKAT", "LOFAR"],
            )


class ReceiverBandTests(unittest.TestCase):

    def test_profile_schema_and_ugmrt_band_edges(self):
        self.assertEqual(PROFILES.PROFILE_SCHEMA_VERSION, 4)
        self.assertEqual(
            PROFILES.TELESCOPE_PROFILES["ugmrt"]["bands"]["B4"]["freq_max"],
            950e6,
        )
        b4 = PROFILES.TELESCOPE_PROFILES["ugmrt"]["bands"]["B4"]
        self.assertEqual(
            (b4["guaranteed_freq_min"], b4["guaranteed_freq_max"]),
            (550e6, 850e6),
        )
        self.assertEqual(
            PROFILES.TELESCOPE_PROFILES["ugmrt"]["bands"]["B2"]["freq_min"],
            125e6,
        )
        self.assertEqual(
            PROFILES.TELESCOPE_PROFILES["ugmrt"]["bands"]["B5"]["freq_max"],
            1460e6,
        )
        self.assertEqual(PROFILES.detect_band("ugmrt", 850e6)["name"], "B4")
        self.assertEqual(PROFILES.detect_band("ugmrt", 1460e6)["name"], "B5")

        self.assertEqual(PROFILES.detect_band("ugmrt", 900e6)["name"], "B4")
        with self.assertRaisesRegex(ValueError, "outside the configured"):
            PROFILES.detect_band("ugmrt", 950e6 + 1.0)
        with self.assertRaisesRegex(ValueError, "outside the configured"):
            PROFILES.detect_band("ugmrt", 1460e6 + 1.0)
        with self.assertRaisesRegex(ValueError, "outside the configured"):
            PROFILES.detect_band("ugmrt", 125e6 - 1.0)

    def test_same_frequency_is_resolved_inside_telescope_namespace(self):
        ugmrt = PROFILES.detect_band("ugmrt", 150e6)
        lofar = PROFILES.detect_band("lofar", 150e6)

        self.assertEqual(ugmrt["name"], "B2")
        self.assertEqual(lofar["name"], "LOFAR_HBA")
        self.assertEqual(ugmrt["telescope"], "ugmrt")
        self.assertEqual(lofar["telescope"], "lofar")

    def test_band_aliases_are_telescope_scoped(self):
        cases = [
            ("ugmrt", "Band-4", "B4"),
            ("meerkat", "UHF", "MEERKAT_UHF"),
            ("meerkat", "L-band", "MEERKAT_L"),
            ("meerkat", "S0", "MEERKAT_S0"),
            ("meerkat", "S-band-4", "MEERKAT_S4"),
            ("lofar", "HBA", "LOFAR_HBA"),
            ("lofar", "low_band", "LOFAR_LBA"),
        ]

        for telescope, value, expected in cases:
            with self.subTest(telescope=telescope, value=value):
                self.assertEqual(
                    PROFILES.normalize_band(telescope, value),
                    expected,
                )

        with self.assertRaisesRegex(ValueError, "Unknown lofar receiver band"):
            PROFILES.normalize_band("lofar", "B4")

        with self.assertRaisesRegex(ValueError, "ambiguous.*S0.*S4"):
            PROFILES.normalize_band("meerkat", "S")

    def test_meerkat_s_tunings_and_science_ranges_are_explicit(self):
        bands = PROFILES.TELESCOPE_PROFILES["meerkat"]["bands"]
        expected_s_ranges = {
            "MEERKAT_S0": (1750e6, 2625e6),
            "MEERKAT_S1": (1968.75e6, 2843.75e6),
            "MEERKAT_S2": (2187.5e6, 3062.5e6),
            "MEERKAT_S3": (2406.25e6, 3281.25e6),
            "MEERKAT_S4": (2625e6, 3500e6),
        }

        self.assertNotIn("MEERKAT_S", bands)
        for band, frequency_range in expected_s_ranges.items():
            with self.subTest(band=band):
                self.assertEqual(
                    (bands[band]["freq_min"], bands[band]["freq_max"]),
                    frequency_range,
                )

        for band, receiver_range, science_range in [
            ("MEERKAT_UHF", (544e6, 1088e6), (580e6, 1015e6)),
            ("MEERKAT_L", (856e6, 1712e6), (900e6, 1670e6)),
        ]:
            with self.subTest(band=band):
                profile = bands[band]
                self.assertEqual(
                    (profile["freq_min"], profile["freq_max"]),
                    receiver_range,
                )
                self.assertEqual(
                    (
                        profile["science_freq_min"],
                        profile["science_freq_max"],
                    ),
                    science_range,
                )
                self.assertEqual(
                    (
                        profile["guaranteed_freq_min"],
                        profile["guaranteed_freq_max"],
                    ),
                    science_range,
                )

        detected = PROFILES.detect_band(
            "meerkat",
            2200e6,
            fmin_hz=1750e6,
            fmax_hz=2625e6,
        )
        self.assertEqual(detected["name"], "MEERKAT_S0")
        with self.assertRaisesRegex(ValueError, "ambiguous"):
            PROFILES.detect_band(
                "meerkat",
                2500e6,
                fmin_hz=2490e6,
                fmax_hz=2510e6,
            )

    def test_every_band_keyed_config_map_contains_s0_through_s4(self):
        expected = {f"MEERKAT_S{index}" for index in range(5)}
        mappings = {
            "edge_fraction": CONFIG.FLAGGING_CONFIG["edge_fraction"],
            "calibration": CONFIG.CAL_CONFIG,
            "selfcal": CONFIG.SELFCAL_CONFIG,
            "selfcal_imaging": CONFIG.SELFCAL_IMAGING_CONFIG,
            "final_imaging": CONFIG.FINAL_IMAGING_CONFIG,
            "ddcal": CONFIG.DD_CONFIG,
        }

        for name, mapping in mappings.items():
            with self.subTest(mapping=name):
                self.assertTrue(expected <= set(mapping))
                self.assertNotIn("MEERKAT_S", mapping)

        with self.assertRaisesRegex(ValueError, "ambiguous"):
            CONFIG.normalize_band("S", telescope="meerkat")

    def test_ugmrt_profiles_publish_measured_polynomial_beams(self):
        expected = {
            "B2": ([-3.089, 39.314, -23.011, 5.037], [125.0, 250.0]),
            "B3": ([-3.129, 38.816, -21.608, 4.483], [250.0, 500.0]),
            "B4": ([-3.263, 42.618, -25.580, 5.823], [550.0, 850.0]),
            "B5": ([-2.614, 27.594, -13.268, 2.395], [1050.0, 1450.0]),
        }

        bands = PROFILES.TELESCOPE_PROFILES["ugmrt"]["bands"]
        self.assertEqual(bands["B2"]["freq_min"], 125e6)
        for band, (coefficients, valid_range) in expected.items():
            with self.subTest(band=band):
                self.assertEqual(
                    bands[band]["primary_beam_provider"],
                    "ugmrt_polynomial",
                )
                self.assertEqual(
                    bands[band]["primary_beam_coefficients"],
                    coefficients,
                )
                self.assertEqual(
                    bands[band]["primary_beam_valid_frequency_mhz"],
                    valid_range,
                )

    def test_ugmrt_is_imaged_only_with_tclean(self):
        for imager in ("wsclean", "other"):
            with self.subTest(imager=imager):
                with self.assertRaisesRegex(
                    ValueError,
                    "validated only with the CASA backend",
                ):
                    PROFILES.validate_stage_capabilities(
                        "ugmrt",
                        "raw",
                        ["imaging"],
                        imager=imager,
                        receiver_band="B4",
                    )
        # Leaving the choice to the pipeline reaches the same backend.
        PROFILES.validate_stage_capabilities(
            "ugmrt", "raw", ["imaging"], imager="auto", receiver_band="B4",
        )
        PROFILES.validate_stage_capabilities(
            "ugmrt", "raw", ["imaging"], imager="casa", receiver_band="B4",
        )

    def test_meerkat_accepts_either_validated_backend(self):
        for imager in ("auto", "wsclean", "casa"):
            with self.subTest(imager=imager):
                PROFILES.validate_stage_capabilities(
                    "meerkat", "raw", ["imaging"], imager=imager,
                    receiver_band="L",
                )
        with self.assertRaisesRegex(ValueError, "Unknown imager"):
            PROFILES.validate_stage_capabilities(
                "meerkat", "raw", ["imaging"], imager="other",
                receiver_band="L",
            )

    def test_meerkat_overlap_requires_an_override_when_coverage_is_ambiguous(self):
        with self.assertRaisesRegex(ValueError, "ambiguous"):
            PROFILES.detect_band(
                "meerkat",
                900e6,
                fmin_hz=900e6,
                fmax_hz=1000e6,
            )

        explicit = PROFILES.detect_band(
            "meerkat",
            900e6,
            fmin_hz=900e6,
            fmax_hz=1000e6,
            requested_band="L",
        )
        self.assertEqual(explicit["name"], "MEERKAT_L")
        self.assertEqual(explicit["method"], "explicit")

    def test_frequency_coverage_can_disambiguate_meerkat_uhf(self):
        result = PROFILES.detect_band(
            "meerkat",
            900e6,
            fmin_hz=700e6,
            fmax_hz=1000e6,
        )

        self.assertEqual(result["name"], "MEERKAT_UHF")
        self.assertEqual(result["method"], "frequency")

    def test_explicit_band_must_contain_the_centre_frequency(self):
        with self.assertRaisesRegex(ValueError, "does not contain"):
            PROFILES.detect_band(
                "lofar",
                150e6,
                requested_band="LBA",
            )

    def test_auto_band_must_contain_the_complete_observed_range(self):
        with self.assertRaisesRegex(ValueError, "complete observed range"):
            PROFILES.detect_band(
                "lofar",
                150e6,
                fmin_hz=80e6,
                fmax_hz=180e6,
            )

        with self.assertRaisesRegex(ValueError, "does not contain"):
            PROFILES.detect_band(
                "meerkat",
                900e6,
                fmin_hz=500e6,
                fmax_hz=1000e6,
                requested_band="UHF",
            )


class TelescopeProfileContractTests(unittest.TestCase):

    def test_profile_getters_return_isolated_deep_copies(self):
        first = PROFILES.get_telescope_profile("meerkat")
        second = PROFILES.get_telescope_profile("meerkat")

        first["preferred_refants"].append("changed")
        first["bands"]["MEERKAT_L"]["cell"] = "99arcsec"

        self.assertNotIn("changed", second["preferred_refants"])
        self.assertEqual(
            second["bands"]["MEERKAT_L"]["cell"],
            "1.5arcsec",
        )
        self.assertNotIn(
            "changed",
            PROFILES.TELESCOPE_PROFILES["meerkat"]["preferred_refants"],
        )

    def test_registry_has_required_capability_and_band_fields(self):
        required_profile_keys = {
            "label",
            "metadata_patterns",
            "antenna_pattern",
            "preferred_refants",
            "raw_calibration_backend",
            "default_imager",
            "casa_selfcal_supported",
            "raw_pipeline_supported",
            "bands",
        }
        required_band_keys = {
            "aliases",
            "label",
            "freq_min",
            "freq_max",
            "cell",
            "imsize",
            "primary_beam_provider",
        }

        self.assertEqual(
            set(PROFILES.TELESCOPE_PROFILES),
            {"ugmrt", "meerkat", "lofar"},
        )

        for telescope, profile in PROFILES.TELESCOPE_PROFILES.items():
            with self.subTest(telescope=telescope):
                self.assertTrue(required_profile_keys <= set(profile))
                self.assertTrue(profile["bands"])
                for band, values in profile["bands"].items():
                    with self.subTest(telescope=telescope, band=band):
                        self.assertTrue(required_band_keys <= set(values))
                        self.assertLess(values["freq_min"], values["freq_max"])
                        self.assertGreater(values["imsize"], 0)
                        if (
                            values["primary_beam_provider"] ==
                            "ugmrt_polynomial"
                        ):
                            self.assertEqual(
                                len(values["primary_beam_coefficients"]),
                                4,
                            )
                            valid_range = values[
                                "primary_beam_valid_frequency_mhz"
                            ]
                            self.assertEqual(len(valid_range), 2)
                            self.assertLess(valid_range[0], valid_range[1])

    def test_backend_defaults_and_capability_guards(self):
        # uGMRT is always tclean; MeerKAT and LOFAR default to WSClean.
        self.assertEqual(PROFILES.default_imager("ugmrt"), "casa")
        self.assertEqual(PROFILES.default_imager("meerkat"), "wsclean")
        self.assertEqual(PROFILES.default_imager("lofar"), "wsclean")

        PROFILES.validate_stage_capabilities(
            "lofar",
            "calibrated",
            ["imaging"],
            imager="auto",
        )

        with self.assertRaisesRegex(ValueError, "Raw LOFAR calibration"):
            PROFILES.validate_stage_capabilities(
                "lofar",
                "raw",
                ["calibration"],
            )

        with self.assertRaisesRegex(ValueError, "self-calibration stage"):
            PROFILES.validate_stage_capabilities(
                "lofar",
                "calibrated",
                ["selfcal"],
            )

        with self.assertRaisesRegex(ValueError, "validated only for uGMRT"):
            PROFILES.validate_stage_capabilities(
                "meerkat",
                "raw",
                ["ddcal"],
            )

        with self.assertRaisesRegex(ValueError, "requires WSClean"):
            PROFILES.validate_stage_capabilities(
                "lofar",
                "calibrated",
                ["imaging"],
                imager="casa",
            )

    def test_lofar_flagging_is_always_rejected(self):
        for state in ("raw", "calibrated"):
            with self.subTest(state=state):
                with self.assertRaisesRegex(
                    ValueError,
                    "flagging stage is disabled for LOFAR",
                ):
                    PROFILES.validate_stage_capabilities(
                        "lofar",
                        state,
                        ["flagging"],
                    )

    def test_known_meerkat_s_band_blocks_pb_dependent_stages(self):
        for stage in ("imaging", "source_detection"):
            with self.subTest(stage=stage):
                with self.assertRaisesRegex(
                    ValueError,
                    "MeerKAT S0-band has no validated primary-beam model",
                ):
                    PROFILES.validate_stage_capabilities(
                        "meerkat",
                        "raw",
                        [stage],
                        receiver_band="S0",
                    )

        # Unknown bands are resolved after inspection; UHF/L have validated
        # katbeam models and remain supported by these profile-level gates.
        PROFILES.validate_stage_capabilities(
            "meerkat",
            "raw",
            ["imaging", "source_detection"],
        )
        for band in ("UHF", "L"):
            with self.subTest(band=band):
                PROFILES.validate_stage_capabilities(
                    "meerkat",
                    "raw",
                    ["imaging", "source_detection"],
                    receiver_band=band,
                )


if __name__ == "__main__":
    unittest.main()
