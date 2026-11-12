"""Field roles from observing intents, in the conventions the exports use.

MeerKAT's katdal export tags a gain calibrator CALIBRATE_PHASE,CALIBRATE_AMPLI
and the flux calibrator CALIBRATE_BANDPASS,CALIBRATE_FLUX. On 1708779077
(2026-11-11) ORBIT read the gain calibrator J2341-3506 as a second flux
standard, which the MeerKAT flux policy refuses at setjy. CALIBRATE_AMPLI on
its own, or with no CALIBRATE_FLUX field anywhere, still marks the flux
calibrator, as ALMA and the VLA use it.
"""

import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest import mock


SOURCE_DIR = Path(__file__).resolve().parents[1]


def load_inspect_module_without_casa():
    casatools = types.ModuleType("casatools")

    class Unused:
        def __init__(self, *args, **kwargs):
            pass

    casatools.msmetadata = Unused
    casatools.table = Unused
    spec = importlib.util.spec_from_file_location(
        "test_radio_continuum_inspection", SOURCE_DIR / "inspect.py")
    module = importlib.util.module_from_spec(spec)
    with mock.patch.dict(sys.modules, {"casatools": casatools}):
        spec.loader.exec_module(module)
    return module


INSPECT = load_inspect_module_without_casa()


class Metadata:
    """The part of casatools.msmetadata that classify_fields reads."""

    def __init__(self, intents):
        self.names = list(intents)
        self.intents = dict(intents)

    def fieldnames(self):
        return list(self.names)

    def fieldsforintent(self, pattern):
        key = pattern.strip("*").upper()
        return [index for index, name in enumerate(self.names)
                if any(key in item.upper() for item in self.intents[name])]

    def intentsforfield(self, field_id):
        return list(self.intents[self.names[field_id]])


def inspector(intents):
    value = INSPECT.MSInspector.__new__(INSPECT.MSInspector)
    value.msmd = Metadata(intents)
    value.flux_catalogue = set()
    value.phase_catalogue = set()
    return value


class IntentTests(unittest.TestCase):

    def test_a_katdal_gain_calibrator_is_the_phase_calibrator(self):
        fields = inspector({
            "J1939-6342": ["CALIBRATE_BANDPASS", "CALIBRATE_FLUX"],
            "J2341-3506": ["CALIBRATE_AMPLI", "CALIBRATE_PHASE"],
            "J2344": ["TARGET"],
        })
        roles = fields.classify_fields()
        self.assertEqual(roles, {
            "flux_calibrators": ["J1939-6342"],
            "bandpass_calibrators": ["J1939-6342"],
            "phase_calibrators": ["J2341-3506"],
            "targets": ["J2344"],
        })
        strategy = fields.calibration_strategy(roles)
        self.assertEqual(strategy["primary_flux_calibrator"], "J1939-6342")
        self.assertIsNone(strategy["secondary_flux_calibrator"])
        self.assertEqual(strategy["phase_calibrator"], "J2341-3506")
        self.assertEqual(strategy["target_fields"], ["J2344"])

    def test_an_amplitude_calibrator_on_its_own_sets_the_flux_scale(self):
        # The ALMA and VLA reading, with no CALIBRATE_FLUX anywhere.
        roles = inspector({
            "J0137+3309": ["CALIBRATE_AMPLI", "CALIBRATE_BANDPASS"],
            "J0238+1636": ["CALIBRATE_PHASE"],
            "NGC1068": ["OBSERVE_TARGET"],
        }).classify_fields()
        self.assertEqual(roles["flux_calibrators"], ["J0137+3309"])
        self.assertEqual(roles["phase_calibrators"], ["J0238+1636"])
        self.assertEqual(roles["targets"], ["NGC1068"])

    def test_amplitude_and_phase_without_a_flux_field_stays_the_flux_scale(self):
        roles = inspector({
            "3C286": ["CALIBRATE_AMPLI", "CALIBRATE_PHASE"],
            "J1331+3030": ["CALIBRATE_PHASE"],
            "FIELD": ["TARGET"],
        }).classify_fields()
        self.assertEqual(roles["flux_calibrators"], ["3C286"])
        self.assertEqual(roles["phase_calibrators"], ["J1331+3030"])

    def test_a_flux_field_that_also_carries_ampli_keeps_its_role(self):
        roles = inspector({
            "J1939-6342": ["CALIBRATE_AMPLI", "CALIBRATE_BANDPASS",
                           "CALIBRATE_FLUX"],
            "J2341-3506": ["CALIBRATE_AMPLI", "CALIBRATE_PHASE"],
            "J2344": ["TARGET"],
        }).classify_fields()
        self.assertEqual(roles["flux_calibrators"], ["J1939-6342"])
        self.assertEqual(roles["phase_calibrators"], ["J2341-3506"])


if __name__ == "__main__":
    unittest.main()
