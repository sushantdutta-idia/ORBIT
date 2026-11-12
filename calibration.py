#!/usr/bin/env python3
#=============================================================
# calibration.py
#
# CASA-native uGMRT and MeerKAT Stokes-I calibration pipeline.
#
# Author: Sushant Dutta
# Compatible with CASA 6+
#
# Run as: casa --nologger --nogui -c calibration.py MS_File.ms
#==============================================================

import importlib.util
import hashlib


# pipeline.py checks this before any stage runs, so a copy of this file from
# an earlier release is reported by name instead of failing inside a stage.
PIPELINE_INTERFACE_VERSION = 5
import itertools
import math
import os
import re
import shutil
import sys
import json
import time
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = os.path.dirname(
    os.path.abspath(__file__)
)

if (
    sys.path and
    os.path.abspath(sys.path[0]) == SCRIPT_DIR
):

    sys.path.pop(0)

import argparse


def load_casa_compat_module():

    compat_path = Path(SCRIPT_DIR) / "casa_compat.py"
    spec = importlib.util.spec_from_file_location(
        "ugmrt_casa_compat",
        str(compat_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


CASA_COMPAT = load_casa_compat_module()
setjy = CASA_COMPAT.import_casa_task("setjy")
ft = CASA_COMPAT.import_casa_task("ft")
gaincal = CASA_COMPAT.import_casa_task("gaincal")
bandpass = CASA_COMPAT.import_casa_task("bandpass")
fluxscale = CASA_COMPAT.import_casa_task("fluxscale")
applycal = CASA_COMPAT.import_casa_task("applycal")
split = CASA_COMPAT.import_casa_task("split")
flagdata = CASA_COMPAT.import_casa_task("flagdata")
flagmanager = CASA_COMPAT.import_casa_task("flagmanager")
# plotms is deliberately disabled by default.  The ilifu CASA 6.7 build
# may expose a partial/broken plotting installation, so merely importing or
# invoking plotms can terminate an otherwise valid reduction.  QA products
# are generated directly from CASA tables by qa_products.py instead.
plotms = None
hanningsmooth = CASA_COMPAT.import_casa_task("hanningsmooth")
statwt = CASA_COMPAT.import_casa_task("statwt", required=False)
from casatools import table


tb = table()


def load_config_module():

    config_path = Path(SCRIPT_DIR) / "config.py"

    if not config_path.exists():

        return None

    removed_script_dir = False

    if (
        sys.path and
        os.path.abspath(sys.path[0]) == SCRIPT_DIR
    ):

        sys.path.pop(0)
        removed_script_dir = True

    try:

        spec = importlib.util.spec_from_file_location(
            "ugmrt_pipeline_config",
            str(config_path),
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        return module

    finally:

        if removed_script_dir:

            sys.path.insert(0, SCRIPT_DIR)


def load_qa_products_module():

    module_path = Path(SCRIPT_DIR) / "qa_products.py"
    spec = importlib.util.spec_from_file_location(
        "ugmrt_qa_products",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_imaging_safety_module():

    module_path = Path(SCRIPT_DIR) / "imaging_safety.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_imaging_safety_for_calibration",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_telescope_profiles_module():

    module_path = Path(SCRIPT_DIR) / "telescope_profiles.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_telescope_profiles",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_field_list(value):

    if value is None:

        return None

    if isinstance(value, (list, tuple)):

        fields = [str(item).strip() for item in value]

    else:

        fields = [item.strip() for item in str(value).split(",")]

    fields = [field for field in fields if field]

    return fields or None


def parse_target_phase_map(value):

    if value is None:

        return None

    mapping = {}

    def add_pair(target, calibrator):

        target_text = str(target).strip()
        calibrator_text = str(calibrator).strip()

        if not target_text or not calibrator_text:
            raise ValueError(
                "Target/phase mappings must contain both a target and a "
                "phase calibrator."
            )

        existing = mapping.get(target_text)
        if existing is not None and existing != calibrator_text:
            raise ValueError(
                "Conflicting phase-calibrator mappings were provided for "
                f"target {target_text!r}."
            )

        mapping[target_text] = calibrator_text

    def parse_item(item):

        if isinstance(item, dict):
            for target, calibrator in item.items():
                add_pair(target, calibrator)
            return

        text = str(item).strip()
        if not text:
            return

        for separator in (":", "="):
            if separator in text:
                target, calibrator = text.split(separator, 1)
                add_pair(target, calibrator)
                return

        raise ValueError(
            "Target/phase mappings must use TARGET:PHASE pairs."
        )

    if isinstance(value, dict):
        parse_item(value)
    elif isinstance(value, (list, tuple)):
        for item in value:
            parse_item(item)
    else:
        for item in str(value).split(","):
            parse_item(item)

    return mapping or None


CONFIG = load_config_module()
RUN_CONFIG = (
    getattr(CONFIG, "PIPELINE_RUN_CONFIG", {})
    if CONFIG is not None
    else {}
)
CENTRAL_CAL_CONFIG = (
    getattr(CONFIG, "CAL_CONFIG", {})
    if CONFIG is not None
    else {}
)
CENTRAL_CAL_TASK_CONFIG = (
    getattr(CONFIG, "CALIBRATION_TASK_CONFIG", {})
    if CONFIG is not None
    else {}
)
TELESCOPE_CAL_TASK_CONFIG = (
    getattr(CONFIG, "TELESCOPE_CALIBRATION_TASK_CONFIG", {})
    if CONFIG is not None
    else {}
)
QA_PRODUCTS = load_qa_products_module()
IMAGING_SAFETY = load_imaging_safety_module()
TELESCOPE_PROFILE_MODULE = load_telescope_profiles_module()


def load_flagging_selectors_module():

    module_path = Path(SCRIPT_DIR) / "flagging_selectors.py"
    spec = importlib.util.spec_from_file_location(
        "radio_continuum_flagging_selectors_for_calibration",
        str(module_path),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FLAGGING_SELECTORS = load_flagging_selectors_module()


def _iter_casa_array_values(value):
    """Yield scalar values from CASA/NumPy arrays without importing NumPy."""

    flat = getattr(value, "flat", None)
    if flat is not None:
        try:
            yield from flat
            return
        except TypeError:
            pass

    if isinstance(value, (list, tuple)):
        for item in value:
            yield from _iter_casa_array_values(item)
        return

    yield value


def sha256_file_or_directory(path):
    """Hash an external CASA model deterministically, without NumPy."""

    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Flux-calibrator model not found: {path}")

    digest = hashlib.sha256()

    def update_file(filename):
        with open(filename, "rb") as handle:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)

    if path.is_file():
        digest.update(b"external-model-file\0")
        update_file(path)
        return digest.hexdigest()

    if not path.is_dir():
        raise RuntimeError(
            f"Flux-calibrator model is neither a file nor directory: {path}"
        )

    digest.update(b"external-model-directory-v1\0")
    transient_names = {".DS_Store", "table.lock"}

    for entry in sorted(path.rglob("*"), key=lambda item: item.as_posix()):
        relative = entry.relative_to(path).as_posix()
        if entry.name in transient_names or entry.name.endswith(".lock"):
            continue

        relative_bytes = relative.encode("utf-8", errors="surrogateescape")
        if entry.is_symlink():
            digest.update(b"L\0" + relative_bytes + b"\0")
            digest.update(os.readlink(entry).encode(
                "utf-8",
                errors="surrogateescape",
            ))
        elif entry.is_dir():
            digest.update(b"D\0" + relative_bytes + b"\0")
        elif entry.is_file():
            digest.update(b"F\0" + relative_bytes + b"\0")
            update_file(entry)

    return digest.hexdigest()


def to_builtin_json(value):
    """Convert a CASA/NumPy result into a JSON-serializable structure."""

    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        try:
            value = tolist()
        except (TypeError, ValueError):
            pass

    if isinstance(value, dict):
        return {str(key): to_builtin_json(item) for key, item in value.items()}

    if isinstance(value, (list, tuple)):
        return [to_builtin_json(item) for item in value]

    try:
        json.dumps(value)
        return value
    except (TypeError, ValueError):
        return str(value)


def _normalized_spw_id(value):

    try:
        return str(int(value))
    except (TypeError, ValueError, OverflowError):
        text = str(value).strip()
        return text or None


def validate_setjy_flux_result(result, expected_spws=None):
    """Validate CASA setjy's returned per-SPW Stokes flux record."""

    if not isinstance(result, dict):
        raise RuntimeError(
            "CASA setjy returned no per-SPW flux-density record"
        )

    records = []

    def visit(node, path=()):
        if not isinstance(node, dict):
            return

        if "fluxd" in node:
            raw_flux = list(_iter_casa_array_values(node["fluxd"]))
            if len(raw_flux) < 4:
                raise RuntimeError(
                    "CASA setjy returned an incomplete Stokes flux vector"
                )

            flux = []
            for value in raw_flux[:4]:
                try:
                    numeric = complex(value)
                except (TypeError, ValueError, OverflowError) as exc:
                    raise RuntimeError(
                        "CASA setjy returned a non-numeric flux density"
                    ) from exc
                if (
                    not math.isfinite(numeric.real) or
                    not math.isfinite(numeric.imag) or
                    numeric.imag != 0.0
                ):
                    raise RuntimeError(
                        "CASA setjy returned a non-finite or complex-valued "
                        "Stokes flux density"
                    )
                flux.append(float(numeric.real))

            if flux[0] <= 0.0:
                raise RuntimeError(
                    "CASA setjy returned a non-positive Stokes-I flux density"
                )

            num_solutions = node.get("numSol")
            if num_solutions is not None:
                try:
                    if int(num_solutions) < 1:
                        raise RuntimeError(
                            "CASA setjy returned zero fitted flux solutions"
                        )
                except (TypeError, ValueError, OverflowError) as exc:
                    raise RuntimeError(
                        "CASA setjy returned an invalid numSol value"
                    ) from exc

            spw_value = (
                node.get("spwID")
                if "spwID" in node
                else node.get("spw_id")
                if "spw_id" in node
                else path[-1]
                if path
                else None
            )
            spw_id = _normalized_spw_id(spw_value)
            if spw_id is None:
                raise RuntimeError(
                    "CASA setjy returned a flux record without an SPW ID"
                )

            records.append({
                "spw_id": spw_id,
                "flux_density_jy": flux,
                "num_solutions": (
                    int(num_solutions)
                    if num_solutions is not None
                    else None
                ),
            })

        for key, value in node.items():
            if isinstance(value, dict):
                visit(value, path + (str(key),))

    visit(result)

    if not records:
        raise RuntimeError(
            "CASA setjy returned no usable per-SPW flux-density records"
        )

    returned_spws = {record["spw_id"] for record in records}
    expected_spws = {
        _normalized_spw_id(spw)
        for spw in (expected_spws or [])
    }
    expected_spws.discard(None)
    missing_spws = sorted(expected_spws.difference(returned_spws))
    if missing_spws:
        raise RuntimeError(
            "CASA setjy did not return flux densities for SPW(s): "
            + ", ".join(missing_spws)
        )

    if all(
        math.isclose(
            record["flux_density_jy"][0],
            1.0,
            rel_tol=0.0,
            abs_tol=1e-6,
        )
        for record in records
    ):
        raise RuntimeError(
            "CASA setjy returned exactly 1 Jy in every SPW, which is its "
            "unsafe fallback for an unrecognized calibrator/model"
        )

    records.sort(key=lambda record: record["spw_id"])
    return {
        "status": "passed",
        "expected_spws": sorted(expected_spws),
        "returned_spws": sorted(returned_spws),
        "records": records,
    }


# MODEL_DATA validation reads at most this many values, in blocks of at most
# MODEL_DATA_CHECK_BLOCK_VALUES. A calibrator field on high-resolution data
# holds tens of billions of values (150-300 GiB); reading it whole killed
# calibration, and checking it value by value in Python would take days.
MODEL_DATA_CHECK_MAX_VALUES = 1 << 26
MODEL_DATA_CHECK_BLOCK_VALUES = 1 << 23


def _count_model_data_values(values):
    """Count total, non-finite and non-zero values in one block."""

    try:
        import numpy as np

        array = np.asarray(values, dtype=complex)
    except (ImportError, TypeError, ValueError):
        array = None

    if array is not None:
        finite = np.isfinite(array)
        total = int(array.size)
        finite_count = int(np.count_nonzero(finite))
        return {
            "total_values": total,
            "finite_values": finite_count,
            "nonfinite_values": total - finite_count,
            "nonzero_values": int(np.count_nonzero(finite & (array != 0))),
        }

    counts = {
        "total_values": 0,
        "finite_values": 0,
        "nonfinite_values": 0,
        "nonzero_values": 0,
    }
    for value in _iter_casa_array_values(values):
        counts["total_values"] += 1
        try:
            numeric = complex(value)
            finite = (
                math.isfinite(numeric.real) and
                math.isfinite(numeric.imag)
            )
        except (TypeError, ValueError, OverflowError):
            finite = False
            numeric = None

        if not finite:
            counts["nonfinite_values"] += 1
            continue

        counts["finite_values"] += 1
        if numeric != 0:
            counts["nonzero_values"] += 1

    return counts


def summarize_model_data(values, counts=None):
    """Validate finite, non-zero MODEL_DATA.

    Pass ``counts`` to judge totals already accumulated over several blocks.
    """

    summary = {
        "status": "failed",
        "total_values": 0,
        "finite_values": 0,
        "nonfinite_values": 0,
        "nonzero_values": 0,
        "reason": "",
    }
    summary.update(
        counts if counts is not None else _count_model_data_values(values)
    )

    if summary["total_values"] < 1:
        summary["reason"] = "selected MODEL_DATA contains no values"
    elif summary["nonfinite_values"]:
        summary["reason"] = (
            f"selected MODEL_DATA contains "
            f"{summary['nonfinite_values']} non-finite values"
        )
    elif summary["nonzero_values"] < 1:
        summary["reason"] = "selected MODEL_DATA is entirely zero"
    else:
        summary["status"] = "passed"

    return summary


def _json_group_value(value):
    """Convert a CASA/NumPy scalar into a deterministic JSON value."""

    item = getattr(value, "item", None)
    if callable(item):
        try:
            value = item()
        except (TypeError, ValueError):
            pass

    if isinstance(value, float) and not math.isfinite(value):
        return str(value)

    if value is None or isinstance(value, (bool, int, float, str)):
        return value

    return str(value)


def summarize_caltable_groups(row_totals,
                              row_flagged,
                              row_usable,
                              group_columns,
                              max_group_flagged_fraction,
                              ignored_row_indices=None):
    """Aggregate row-aligned calibration QA by CASA identifier columns."""

    metrics = {}
    first_failure = ""
    ignored_row_indices = set(ignored_row_indices or ())

    for column_name, raw_values in group_columns.items():
        aggregates = {}
        order = []

        for row_index, raw_value in enumerate(raw_values):
            if row_index in ignored_row_indices:
                continue
            value = _json_group_value(raw_value)
            token = (type(value).__name__, repr(value))

            if token not in aggregates:
                aggregates[token] = {
                    "value": value,
                    "rows": 0,
                    "total_solutions": 0,
                    "flagged_solutions": 0,
                    "usable_solutions": 0,
                }
                order.append(token)

            group = aggregates[token]
            group["rows"] += 1
            group["total_solutions"] += row_totals[row_index]
            group["flagged_solutions"] += row_flagged[row_index]
            group["usable_solutions"] += row_usable[row_index]

        groups = []
        column_failed = False

        for token in order:
            group = aggregates[token]
            total = group["total_solutions"]
            fraction = (
                group["flagged_solutions"] / total
                if total
                else None
            )
            group["flagged_fraction"] = fraction
            group["usable_fraction"] = (
                group["usable_solutions"] / total
                if total
                else None
            )
            group["status"] = "passed"
            group["reason"] = ""

            if group["usable_solutions"] < 1:
                group["status"] = "failed"
                group["reason"] = "group has no usable solution"
            elif (
                fraction is not None and
                fraction > max_group_flagged_fraction
            ):
                group["status"] = "failed"
                group["reason"] = (
                    f"flagged fraction {fraction:.6f} exceeds the "
                    "configured group maximum "
                    f"{max_group_flagged_fraction:.6f}"
                )

            if group["status"] == "failed":
                column_failed = True
                if not first_failure:
                    first_failure = (
                        f"{column_name}={group['value']!r}: "
                        f"{group['reason']}"
                    )

            groups.append(group)

        metrics[column_name] = {
            "status": "failed" if column_failed else "passed",
            "group_count": len(groups),
            "groups": groups,
        }

    return metrics, first_failure


def summarize_caltable_solutions(parameter_column,
                                 parameters,
                                 flags,
                                 max_flagged_fraction,
                                 row_count=None,
                                 group_columns=None,
                                 max_group_flagged_fraction=0.80,
                                 ignored_group_row_indices=None):
    """Return fail-closed QA metrics for one CASA solution column."""

    summary = {
        "parameter_column": parameter_column,
        "max_flagged_solution_fraction": max_flagged_fraction,
        "max_group_flagged_solution_fraction": (
            max_group_flagged_fraction
        ),
        "total_solutions": 0,
        "flagged_solutions": 0,
        "flagged_fraction": None,
        "usable_solutions": 0,
        "nonfinite_usable_solutions": 0,
        "zero_usable_solutions": 0,
        "group_qa_status": "not_evaluated",
        "group_metrics": {},
        "ignored_group_row_indices": sorted(
            set(ignored_group_row_indices or ())
        ),
        "status": "failed",
        "reason": "",
    }

    group_columns = group_columns or {}
    if group_columns and row_count is None:
        summary["reason"] = (
            "group identifier columns require a calibration-table row count"
        )
        return summary

    if row_count is not None:
        try:
            row_count = int(row_count)
        except (TypeError, ValueError):
            summary["reason"] = "calibration-table row count is invalid"
            return summary

        if row_count < 1:
            summary["reason"] = "calibration-table row count is invalid"
            return summary

        row_totals = [0] * row_count
        row_flagged = [0] * row_count
        row_usable = [0] * row_count

        normalized_groups = {}
        for column_name, values in group_columns.items():
            values = list(_iter_casa_array_values(values))
            if len(values) != row_count:
                summary["reason"] = (
                    f"{column_name} has {len(values)} values for "
                    f"{row_count} calibration-table rows"
                )
                return summary
            normalized_groups[column_name] = values
        group_columns = normalized_groups
    else:
        row_totals = row_flagged = row_usable = None

    sentinel = object()
    pairs = itertools.zip_longest(
        _iter_casa_array_values(parameters),
        _iter_casa_array_values(flags),
        fillvalue=sentinel,
    )

    for solution_index, (parameter, flagged) in enumerate(pairs):
        if parameter is sentinel or flagged is sentinel:
            summary["reason"] = (
                f"{parameter_column} and FLAG have different numbers "
                "of solution elements"
            )
            return summary

        summary["total_solutions"] += 1

        if row_totals is not None:
            row_index = solution_index % row_count
            row_totals[row_index] += 1

        try:
            is_flagged = bool(flagged)
        except Exception:
            summary["reason"] = "FLAG contains an unreadable value"
            return summary

        if is_flagged:
            summary["flagged_solutions"] += 1
            if row_flagged is not None:
                row_flagged[row_index] += 1
            continue

        try:
            numeric = complex(parameter)
            is_finite = (
                math.isfinite(numeric.real) and
                math.isfinite(numeric.imag)
            )
        except (TypeError, ValueError, OverflowError):
            is_finite = False
            numeric = None

        if not is_finite:
            summary["nonfinite_usable_solutions"] += 1
        elif parameter_column == "CPARAM" and numeric == 0:
            summary["zero_usable_solutions"] += 1
        else:
            summary["usable_solutions"] += 1
            if row_usable is not None:
                row_usable[row_index] += 1

    total = summary["total_solutions"]
    if total < 1:
        summary["reason"] = f"{parameter_column} contains no solutions"
        return summary

    flagged_fraction = summary["flagged_solutions"] / total
    summary["flagged_fraction"] = flagged_fraction

    if row_count is not None and total % row_count:
        summary["reason"] = (
            "solution-array size is not divisible by the calibration-table "
            "row count; per-group coverage cannot be aligned safely"
        )
        return summary

    group_failure = ""
    if group_columns:
        group_metrics, group_failure = summarize_caltable_groups(
            row_totals,
            row_flagged,
            row_usable,
            group_columns,
            max_group_flagged_fraction,
            ignored_row_indices=ignored_group_row_indices,
        )
        summary["group_metrics"] = group_metrics
        summary["group_qa_status"] = (
            "failed" if group_failure else "passed"
        )

    if summary["usable_solutions"] < 1:
        summary["reason"] = "all calibration solutions are flagged"
        return summary

    if flagged_fraction > max_flagged_fraction:
        summary["reason"] = (
            "flagged calibration-solution fraction "
            f"{flagged_fraction:.6f} exceeds the configured maximum "
            f"{max_flagged_fraction:.6f}"
        )
        return summary

    if summary["nonfinite_usable_solutions"]:
        summary["reason"] = (
            f"{summary['nonfinite_usable_solutions']} unflagged "
            "calibration solutions are non-finite"
        )
        return summary

    if group_failure:
        summary["reason"] = group_failure
        return summary

    if summary["zero_usable_solutions"]:
        summary["reason"] = (
            f"{summary['zero_usable_solutions']} unflagged CPARAM "
            "solutions are zero and cannot be applied"
        )
        return summary

    summary["status"] = "passed"
    return summary


class CalibrationManager:
    """
    Telescope-aware uGMRT/MeerKAT calibration manager.
    """

    # ---------------------------------------------------------
    # Band-dependent calibration configuration
    # ---------------------------------------------------------

    CAL_CONFIG = {

        "B2": {
        "delay_solint": "inf",
        "phase_solint": "60s",
        "amp_solint": "600s",
        "combine_delay": "scan",
        "hanningsmooth": True
    },

    "B3": {
        "delay_solint": "inf",
        "phase_solint": "30s",
        "amp_solint": "300s",
        "combine_delay": "scan",
        "hanningsmooth": True
    },

    "B4": {
        "delay_solint": "inf",
        "phase_solint": "20s",
        "amp_solint": "180s",
        "combine_delay": "scan",
        "hanningsmooth": True,
        "target_max_channel_width_khz": 400.0
    },

    "B5": {
        "delay_solint": "inf",
        "phase_solint": "10s",
        "amp_solint": "120s",
        "combine_delay": "",
        "hanningsmooth": False
    }
    }

    # ---------------------------------------------------------

    def __init__(self,
                 ms,
                 outdir="Products",
                 flux_calibrator=None,
                 bandpass_calibrator=None,
                 secondary_flux_calibrator=None,
                 phase_calibrator=None,
                 target_phase_map=None,
                 target_fields=None,
                 reference_antenna=None,
                 telescope=None,
                 flux_standard=None,
                 flux_model=None,
                 flux_model_type="casa-image",
                 provenance_fingerprint=None,
                 recipe_fingerprint=None):

        self.ms = ms
        self.input_ms = ms
        self.requested_telescope = telescope
        self.flux_standard_override = flux_standard
        self.flux_model = flux_model
        self.flux_model_type = str(
            flux_model_type or "casa-image"
        ).strip().lower()
        if self.flux_model_type not in {"casa-image", "component-list"}:
            raise ValueError(
                "flux_model_type must be casa-image or component-list"
            )
        self.provenance_fingerprint = (
            str(provenance_fingerprint).strip()
            if provenance_fingerprint
            else None
        )
        self.recipe_fingerprint = (
            str(recipe_fingerprint).strip()
            if recipe_fingerprint
            else self.provenance_fingerprint
        )
        self.field_overrides = {
            "flux_calibrator": parse_field_list(flux_calibrator),
            "bandpass_calibrator": parse_field_list(bandpass_calibrator),
            "secondary_flux_calibrator": parse_field_list(
                secondary_flux_calibrator
            ),
            "phase_calibrator": parse_field_list(phase_calibrator),
            "target_phase_map": parse_target_phase_map(target_phase_map),
            "target_fields": parse_field_list(target_fields),
            "reference_antenna": reference_antenna,
        }

        if CENTRAL_CAL_CONFIG:

            self.CAL_CONFIG = CENTRAL_CAL_CONFIG

        self.task_config = {
            "flux_standard": "Perley-Butler 2017",
            "delay_minsnr": 3.0,
            "gain_minsnr": 3.0,
            "minblperant": 4,
            "gain_solmode": "L1R",
            "prebp_gaintype": "G",
            "bandpass_fillgaps": 8,
            "bandpass_solint": "inf",
            "bandpass_solnorm": True,
            "fluxscale_fitorder": 1,
            "refantmode": "strict",
            "target_phase_transfer_solint": "inf",
            "max_caltable_flagged_fraction": 0.50,
            "max_caltable_group_flagged_fraction": 0.80,
            "auto_quarantine_failed_delay_antennas": True,
            "max_auto_quarantine_delay_antennas": 1,
            "max_auto_quarantine_delay_fraction": 0.10,
            "min_antennas_after_delay_quarantine": 8,
            "applymode": "calflagstrict",
            "calwt": True,
            "parang": False,
        }
        self.task_config.update(
            CENTRAL_CAL_TASK_CONFIG
        )

        self.outdir = Path(outdir).expanduser().resolve()

        #
        # Restore Hanning-smoothed MS
        # if it already exists
        #

        hann_ms = (
            self.outdir /
            "calibration" /
            "hann.ms"
        )

        if (
            hann_ms.is_dir() and
            (hann_ms / "table.dat").is_file()
        ):

            print(
                f"Using existing smoothed MS:\n"
                f"{hann_ms}"
            )
            self.ms = str(hann_ms)

        self.caldir = (
            self.outdir / "calibration"
        )

        self.plotdir = (
            self.caldir / "plots"
        )

        self.splitdir = (
            self.caldir / "calibrated"
        )

        # Create directories

        self.caldir.mkdir(
            parents=True,
            exist_ok=True
        )

        self.plotdir.mkdir(
            exist_ok=True
        )

        self.splitdir.mkdir(
            exist_ok=True
        )

        # Load reports

        self.report = (
            self.load_inspection_report()
        )

        self.telescope = TELESCOPE_PROFILE_MODULE.telescope_name_from_report(
            self.report
        )
        if self.telescope is None:
            raise RuntimeError(
                "Inspection report has no supported telescope identity. "
                "Rerun inspection with --telescope."
            )

        if self.requested_telescope not in (None, "", "auto"):
            requested = TELESCOPE_PROFILE_MODULE.normalize_telescope(
                self.requested_telescope
            )
            if requested != self.telescope:
                raise RuntimeError(
                    f"Requested telescope {requested} conflicts with the "
                    f"inspection report ({self.telescope})."
                )

        if self.telescope == "lofar":
            raise RuntimeError(
                "Raw LOFAR calibration cannot use the CASA uGMRT/MeerKAT "
                "recipe. Run LINC/DP3 first and pass its target MS with "
                "--data-state calibrated."
            )

        self.task_config.update(
            TELESCOPE_CAL_TASK_CONFIG.get(self.telescope, {})
        )
        if self.flux_standard_override:
            self.task_config["flux_standard"] = self.flux_standard_override

        self.flagging = (
            self.load_flagging_outputs()
        )
        self.load_delay_quarantine()

        # Band configuration

        self.band = self.get_band_name()

        self.band_config = (
            self.choose_band_config()
        )

        print(
            f"Loaded calibration "
            f"configuration for {self.telescope}/{self.band}"
        )

    # =========================================================
    # REPORT LOADING
    # =========================================================

    def load_inspection_report(self):

        report_file = (
            self.outdir /
            "inspection" /
            "inspection_report.json"
        )

        if not report_file.exists():

            raise FileNotFoundError(
                f"Inspection report not found:\n"
                f"{report_file}"
            )

        with open(report_file) as f:

            report = json.load(f)

        print(
            "Loaded inspection report."
        )

        return report

    # ---------------------------------------------------------

    def load_flagging_outputs(self):

        flagging = {}

        # Recommended refant

        refant_file = (
            self.outdir /
            "flagging" /
            "recommended_refant.txt"
        )

        if refant_file.exists():

            with open(refant_file) as f:

                flagging["recommended_refant"] = (
                    f.read().strip()
                )

        # Bad antennas

        badant_file = (
            self.outdir /
            "flagging" /
            "bad_antennas.txt"
        )

        if badant_file.exists():

            with open(badant_file) as f:

                bad = []

                for line in f:

                    ant = line.strip()

                    if ant:

                        bad.append(ant)

                flagging["bad_antennas"] = bad

        else:

            flagging["bad_antennas"] = []

        print(
            "Loaded flagging outputs."
        )

        return flagging

    # ---------------------------------------------------------

    def get_antenna_selection(self):

        bad = self.flagging["bad_antennas"]

        if not bad:

            return ""

        return "!" + ",".join(bad)

    # ---------------------------------------------------------

    def delay_quarantine_path(self):

        return self.caldir / "delay_quarantined_antennas.json"

    # ---------------------------------------------------------

    def load_delay_quarantine(self):
        """Restore a previously committed delay-antenna exclusion."""

        path = self.delay_quarantine_path()
        self.delay_quarantine = {
            "schema_version": 1,
            "measurement_set": str(
                Path(self.input_ms).expanduser().resolve()
            ),
            "antennas": [],
        }
        if not path.is_file():
            return

        try:
            with open(path) as handle:
                stored = json.load(handle)
        except (OSError, ValueError, TypeError) as exc:
            raise RuntimeError(
                f"Delay-antenna quarantine is unreadable: {path}"
            ) from exc

        expected_ms = self.delay_quarantine["measurement_set"]
        if (
            not isinstance(stored, dict) or
            stored.get("schema_version") != 1 or
            stored.get("measurement_set") != expected_ms or
            not isinstance(stored.get("antennas"), list)
        ):
            raise RuntimeError(
                "Delay-antenna quarantine does not belong to the current "
                f"Measurement Set: {path}"
            )

        names = []
        for item in stored["antennas"]:
            if not isinstance(item, dict):
                raise RuntimeError(
                    f"Delay-antenna quarantine is malformed: {path}"
                )
            name = str(item.get("name", "")).strip()
            try:
                int(item.get("id"))
            except (TypeError, ValueError):
                name = ""
            if not name:
                raise RuntimeError(
                    f"Delay-antenna quarantine is malformed: {path}"
                )
            names.append(name)

        self.delay_quarantine = stored
        self.flagging["bad_antennas"] = list(dict.fromkeys(
            list(self.flagging.get("bad_antennas", [])) + names
        ))
        if names:
            print(
                "Restored delay-solution antenna quarantine: "
                + ", ".join(names)
            )

    # ---------------------------------------------------------

    def inspection_antenna_id_name_map(self):

        path = self.outdir / "inspection" / "antennas.json"
        try:
            with open(path) as handle:
                entries = json.load(handle)
        except (OSError, ValueError, TypeError) as exc:
            raise RuntimeError(
                "Cannot map a failed calibration antenna ID to a name "
                f"because {path} is missing or unreadable."
            ) from exc

        mapping = {}
        # Inspection products written before schema 2 stored only antenna
        # names. Their list order is the Measurement Set ANTENNA row order,
        # so enumerate is the authoritative backward-compatible ID mapping.
        for list_index, item in enumerate(
            entries if isinstance(entries, list) else []
        ):
            if not isinstance(item, dict):
                continue
            try:
                antenna_id = int(item.get("id", list_index))
            except (TypeError, ValueError):
                continue
            name = str(item.get("name", "")).strip()
            if name:
                mapping[antenna_id] = name

        if not mapping:
            raise RuntimeError(
                f"Inspection antenna inventory is empty: {path}"
            )
        return mapping

    # ---------------------------------------------------------

    def solver_excluded_antenna_ids(self):
        """Return inspection IDs deliberately excluded from gain solves."""

        bad = list(
            getattr(self, "flagging", {}).get("bad_antennas", [])
        )
        if not bad:
            return set()

        mapping = self.inspection_antenna_id_name_map()
        by_name = {name: antenna_id for antenna_id, name in mapping.items()}
        excluded = set()
        unknown = []
        for value in bad:
            text = str(value).strip()
            if text in by_name:
                excluded.add(by_name[text])
                continue
            try:
                antenna_id = int(text)
            except ValueError:
                unknown.append(text)
                continue
            if antenna_id in mapping:
                excluded.add(antenna_id)
            else:
                unknown.append(text)

        if unknown:
            raise RuntimeError(
                "Excluded calibration antenna(s) are absent from the "
                "inspection inventory: " + ", ".join(unknown)
            )
        return excluded

    # ---------------------------------------------------------

    def failed_delay_antenna_ids(self, filename="delay.cal",
                                 max_usable_fraction=0.0):
        """Return a narrowly recoverable set of all-unsolved antennas.

        Recovery is allowed only when ANTENNA1 is the sole failed grouping.
        An SPW, field, observation, non-finite, or globally sparse failure is
        never converted into an antenna exclusion. An antenna counts when it
        has no usable solution, or (``max_usable_fraction`` > 0, used for the
        per-channel bandpass) when no more than that share is usable.
        """

        path = self.caltable_qa_path(filename)
        try:
            with open(path) as handle:
                qa = json.load(handle)
        except (OSError, ValueError, TypeError):
            return []

        if not isinstance(qa, dict):
            return []
        if qa.get("nonfinite_usable_solutions", 0):
            return []

        groups = qa.get("group_metrics", {})
        if not isinstance(groups, dict):
            return []
        antenna_metrics = groups.get("ANTENNA1")
        if not isinstance(antenna_metrics, dict):
            return []

        for group_name, metrics in groups.items():
            if (
                group_name != "ANTENNA1" and
                isinstance(metrics, dict) and
                metrics.get("status") == "failed"
            ):
                return []

        failed = []
        for group in antenna_metrics.get("groups", []):
            if not isinstance(group, dict) or group.get("status") != "failed":
                continue
            none_usable = (
                group.get("reason") == "group has no usable solution" and
                int(group.get("usable_solutions", 0)) == 0
            )
            try:
                usable_fraction = float(group.get("usable_fraction"))
            except (TypeError, ValueError):
                usable_fraction = None
            nearly_none = (
                max_usable_fraction > 0 and
                usable_fraction is not None and
                usable_fraction <= max_usable_fraction
            )
            if not (none_usable or nearly_none):
                return []
            try:
                failed.append(int(group["value"]))
            except (KeyError, TypeError, ValueError):
                return []

        return sorted(set(failed))

    # ---------------------------------------------------------

    def quarantine_failed_delay_antennas(self, antenna_ids):
        """Commit one bounded delay exclusion and return antenna names."""

        return self.quarantine_unsolvable_antennas(
            antenna_ids, "delay.cal", "no usable delay solution"
        )

    def unsolvable_table_message(self, filename, reason, failed_ids):
        """Why a solve failed, and what would let it proceed."""

        message = (
            "CASA did not produce a usable calibration table "
            f"{self.caldir / filename}: {reason}"
        )
        if not failed_ids:
            return message
        try:
            mapping = self.inspection_antenna_id_name_map()
        except Exception:
            mapping = {}
        names = ", ".join(
            str(mapping.get(value, value)) for value in failed_ids
        )
        return (
            message + f". Antennas with no usable solution at all: {names}. "
            "They were not excluded automatically: "
            + self.quarantine_refusal_reasons(failed_ids, mapping)
            + ". If they are not working in this observation, raise that "
            "limit in config.yaml (the run then excludes them and goes on) "
            "or flag them in the flagging stage."
        )

    def quarantine_refusal_reasons(self, failed_ids, mapping):
        """Which quarantine limit stopped these antennas being excluded."""

        config = self.task_config
        maximum = int(config.get("max_auto_quarantine_delay_antennas", 1))
        fraction = float(config.get("max_auto_quarantine_delay_fraction", 0.10))
        minimum = int(config.get("min_antennas_after_delay_quarantine", 8))
        already = {
            int(item["id"]): item.get("name", item["id"])
            for item in self.delay_quarantine.get("antennas", [])
        }
        failed = set(int(value) for value in failed_ids)
        try:
            active = self.active_measurement_set_antenna_ids()
        except Exception:
            active = set()
        reasons = []
        if not config.get("auto_quarantine_failed_delay_antennas", True):
            reasons.append(
                "calibration_auto_quarantine_failed_delay_antennas is off"
            )
        if len(failed) > maximum:
            reasons.append(
                f"{len(failed)} at once is more than "
                f"calibration_max_auto_quarantine_delay_antennas ({maximum})"
            )
        if active:
            total = len(set(already) | failed)
            if total / len(active) > fraction:
                before = ", ".join(str(name) for name in already.values())
                reasons.append(
                    f"with the {len(already)} excluded before"
                    + (f" ({before})" if before else "")
                    + f" that makes {total} of {len(active)} antennas "
                    f"({100 * total / len(active):.0f}%), more than "
                    f"calibration_max_auto_quarantine_delay_fraction "
                    f"({fraction:g}); {math.ceil(100 * total / len(active)) / 100:g} "
                    "would allow it"
                )
            if len(active) - total < minimum:
                reasons.append(
                    f"only {len(active) - total} antennas would remain, fewer "
                    f"than calibration_min_antennas_after_delay_quarantine "
                    f"({minimum})"
                )
        if failed & set(already):
            reasons.append("one of them was excluded already")
        if not reasons:
            reasons.append(
                "one is the reference antenna or not in the measurement set"
            )
        return "; ".join(reasons)

    def quarantine_unsolvable_antennas(self, antenna_ids, filename, reason):
        """Exclude antennas a solve could not use at all, within bounds.

        An antenna with no usable solution on a bright calibrator, over every
        interval, is not a calibration failure but an antenna that is not
        working: a dead or disconnected receiver produces noise that
        interference flagging leaves alone, because nothing in it looks like
        interference. It is excluded from every calibration solve (and its
        target data are flagged by applycal for want of solutions) only
        within the same limits as always: at most
        ``max_auto_quarantine_delay_antennas`` per event, the array's total
        exclusions within ``max_auto_quarantine_delay_fraction``, at least
        ``min_antennas_after_delay_quarantine`` left, and never the
        reference antenna.
        """

        if not self.task_config.get(
            "auto_quarantine_failed_delay_antennas", True
        ):
            return []

        antenna_ids = sorted(set(int(value) for value in antenna_ids))
        if not antenna_ids:
            return []

        maximum = int(self.task_config.get(
            "max_auto_quarantine_delay_antennas", 1
        ))
        active = self.active_measurement_set_antenna_ids()
        maximum_fraction = float(self.task_config.get(
            "max_auto_quarantine_delay_fraction", 0.10
        ))
        minimum_remaining = int(self.task_config.get(
            "min_antennas_after_delay_quarantine", 8
        ))

        if (
            maximum < 0 or
            not 0.0 <= maximum_fraction <= 1.0 or
            minimum_remaining < 2
        ):
            raise RuntimeError(
                "Delay-antenna quarantine limits are invalid in config.yaml."
            )
        if not active or any(value not in active for value in antenna_ids):
            return []
        already = {
            int(item["id"]) for item in self.delay_quarantine.get("antennas", [])
        }
        excluded_total = len(already | set(antenna_ids))
        if (
            len(antenna_ids) > maximum or
            excluded_total / len(active) > maximum_fraction or
            len(active) - excluded_total < minimum_remaining
        ):
            return []

        mapping = self.inspection_antenna_id_name_map()
        try:
            antennas = [
                {
                    "id": value,
                    "name": mapping[value],
                    "source_table": str(self.caldir / filename),
                    "reason": reason,
                }
                for value in antenna_ids
            ]
        except KeyError as exc:
            raise RuntimeError(
                f"Inspection has no name for failed antenna ID {exc.args[0]}."
            ) from exc

        names = [item["name"] for item in antennas]
        reference = self.get_reference_antenna()
        if reference in names:
            raise RuntimeError(
                f"Reference antenna {reference} has {reason} in {filename}. "
                "Select a different --reference-antenna; it will not be "
                "quarantined automatically."
            )

        existing = {
            int(item["id"]): item
            for item in self.delay_quarantine.get("antennas", [])
        }
        if any(value in existing for value in antenna_ids):
            return []
        existing.update({item["id"]: item for item in antennas})
        record = {
            "schema_version": 1,
            "measurement_set": str(
                Path(self.input_ms).expanduser().resolve()
            ),
            "source_table": str(self.caldir / filename),
            "reason": reason,
            "antennas": [existing[key] for key in sorted(existing)],
        }
        self.write_atomic_json(self.delay_quarantine_path(), record)
        self.delay_quarantine = record
        self.flagging["bad_antennas"] = list(dict.fromkeys(
            list(self.flagging.get("bad_antennas", [])) + names
        ))

        print(
            f"WARNING: Quarantining an antenna with {reason} in {filename}: "
            + ", ".join(
                f"{item['name']} (ID {item['id']})" for item in antennas
            )
        )
        print(
            "It will be excluded from every calibration solve, and the "
            "solves it already took part in are recomputed without it."
        )
        return names

    # =========================================================
    # BAND CONFIGURATION
    # =========================================================

    def choose_band_config(self):

        if self.band not in self.CAL_CONFIG:

            raise RuntimeError(
                f"Unsupported band: "
                f"{self.band}"
            )

        return self.CAL_CONFIG[
            self.band
        ]

    # =========================================================
    # REPORT NORMALIZATION
    # =========================================================

    def get_band_name(self):

        band = self.report.get(
            "band",
            {}
        )

        if isinstance(band, dict):

            band = (
                band.get("name") or
                band.get("band")
            )

        if not band:

            raise RuntimeError(
                "Inspection report does not contain a "
                "usable band name."
            )

        band = str(band).upper()
        band = band.replace("BAND-", "B")
        band = band.replace("BAND", "B")

        return band

    # =========================================================
    # CALIBRATOR SELECTION
    # =========================================================

    def get_flux_calibrators(self):

        overrides = self.field_overrides["flux_calibrator"] or []
        secondary_overrides = (
            self.field_overrides["secondary_flux_calibrator"] or []
        )

        if overrides or secondary_overrides:

            fields = []

            if overrides:
                fields.extend(
                    self.validate_field_list(
                        overrides,
                        "flux calibrator",
                    )
                )
            else:
                candidates = list(dict.fromkeys(
                    self.report.get("flux_calibrators") or
                    self.report.get("fields", {}).get(
                        "flux_calibrators",
                        []
                    )
                ))
                # With two standards observed (uGMRT Band 4: 3C48 first,
                # 3C147 last) inspection's strategy names one as the
                # primary; honour it, as get_bandpass_calibrator does, when
                # it is still among the vetted flux calibrators. A named
                # secondary replaces inspection's own.
                strategy_primary = (
                    self.report.get("calibration_strategy", {}).get(
                        "primary_flux_calibrator"
                    )
                )
                if len(candidates) > 1 and strategy_primary in candidates:
                    self.validate_field_name(
                        strategy_primary,
                        "inspection-selected primary flux calibrator",
                    )
                    inferred_primary = strategy_primary
                else:
                    inferred_primary = self.require_unique_field(
                        candidates,
                        "primary flux calibrator",
                        "--flux-calibrator",
                    )
                fields.append(inferred_primary)

            if secondary_overrides:
                fields.extend(
                    self.validate_field_list(
                        secondary_overrides,
                        "secondary flux calibrator",
                    )
                )

            return list(dict.fromkeys(fields))

        flux = (
            self.report.get("flux_calibrators") or
            self.report.get("fields", {}).get(
                "flux_calibrators",
                []
            )
        )

        return self.require_unique_field(
            flux,
            "flux calibrator",
            "--flux-calibrator",
            allow_multiple=True,
        )

    def get_flux_calibrator(self):

        return self.get_flux_calibrators()[0]

    # ---------------------------------------------------------

    def get_secondary_flux_calibrator(self):

        flux = self.get_flux_calibrators()

        if len(flux) > 1:

            return flux[1]

        return None

    # ---------------------------------------------------------

    def get_bandpass_calibrator(self):

        override = self.field_overrides["bandpass_calibrator"] or []

        if override:
            if len(override) != 1:
                raise RuntimeError(
                    "Exactly one bandpass calibrator may be configured. "
                    "Pass a single --bandpass-calibrator field."
                )
            self.validate_field_name(override[0], "bandpass calibrator")
            return override[0]

        candidates = (
            self.report.get("bandpass_calibrators") or
            self.report.get("fields", {}).get(
                "bandpass_calibrators",
                [],
            )
        )

        if not candidates:
            return self.get_flux_calibrator()

        # Inspection records a deterministic calibration strategy even when
        # source-name heuristics classify more than one standard as capable
        # of bandpass calibration. Honour that role assignment only when it
        # is still present in the independently vetted candidate list. This
        # resolves the common inferred 3C147(primary+BP), 3C48(secondary)
        # uGMRT arrangement without choosing an arbitrary list element.
        strategy_bandpass = (
            self.report.get("calibration_strategy", {}).get(
                "bandpass_calibrator"
            )
        )
        if strategy_bandpass in candidates:
            self.validate_field_name(
                strategy_bandpass,
                "inspection-selected bandpass calibrator",
            )
            return strategy_bandpass

        # An explicit primary flux-calibrator choice is also a safe
        # deterministic bandpass choice when inspection independently
        # classified that same field as a bandpass candidate.  This is the
        # common 3C147(primary+BP), 3C48(secondary) uGMRT arrangement.  Keep
        # failing closed when the explicit primary is not in the vetted
        # bandpass candidate set, because that may indicate a genuinely
        # distinct bandpass calibrator or incorrect field metadata.
        if self.field_overrides["flux_calibrator"]:
            primary_flux = self.get_flux_calibrator()
            if primary_flux in candidates:
                return primary_flux

        return self.require_unique_field(
            candidates,
            "bandpass calibrator",
            "--bandpass-calibrator",
        )

    # ---------------------------------------------------------

    def get_phase_calibrators(self):

        override = self.field_overrides["phase_calibrator"] or []

        if override:

            return self.validate_field_list(
                override,
                "phase calibrator",
            )

        phase = (
            self.report.get("phase_calibrators") or
            self.report.get("fields", {}).get(
                "phase_calibrators",
                []
            )
        )

        return self.require_unique_field(
            phase,
            "phase calibrator",
            "--phase-calibrator",
            allow_multiple=True,
        )

    def get_phase_calibrator(self):

        return self.get_phase_calibrators()[0]

    # ---------------------------------------------------------

    def get_target_phase_map(self):

        targets = self.get_targets()
        phases = self.get_phase_calibrators()
        configured = self.field_overrides["target_phase_map"]

        if configured:

            mapping = {}
            for target, calibrator in configured.items():
                self.validate_field_name(
                    target,
                    "target in --target-phase-map",
                )
                self.validate_field_name(
                    calibrator,
                    "phase calibrator in --target-phase-map",
                )
                mapping[str(target).strip()] = str(calibrator).strip()

            missing = [
                target
                for target in targets
                if target not in mapping
            ]
            if missing:
                raise RuntimeError(
                    "The target/phase mapping is incomplete. Missing "
                    f"targets: {', '.join(missing)}."
                )

            unknown_targets = sorted(
                set(mapping) - set(targets)
            )
            if unknown_targets:
                raise RuntimeError(
                    "The target/phase mapping references unknown or "
                    "unselected targets: "
                    + ", ".join(unknown_targets)
                )

            unknown_phases = sorted(
                set(mapping.values()) - set(phases)
            )
            if unknown_phases:
                raise RuntimeError(
                    "The target/phase mapping references phase "
                    "calibrators that are not selected for this run: "
                    + ", ".join(unknown_phases)
                )

            return {
                target: mapping[target]
                for target in targets
            }

        if not targets:

            return {}

        if len(phases) == 1:

            return {
                target: phases[0]
                for target in targets
            }

        raise RuntimeError(
            "Multiple phase calibrators are available: "
            + ", ".join(phases)
            + ". Pass --target-phase-map TARGET:PHASE,... to choose "
              "which phase calibrator should be applied to each target."
        )

    # ---------------------------------------------------------

    def get_targets(self):
        overrides = self.field_overrides["target_fields"]

        if overrides:

            for field in overrides:

                self.validate_field_name(field, "target")

            return overrides

        targets = (
            self.report.get("target_fields") or
            self.report.get("fields", {}).get(
                "targets",
                []
            )
        )

        return self.require_unique_field(targets, "target field", "--targets", allow_multiple=True)

    # ---------------------------------------------------------

    def get_reference_antenna(self):

        override = self.field_overrides["reference_antenna"]

        if override:

            return self.validate_reference_antenna(override)

        # Priority 1:
        # flagging recommendation

        if (
            "recommended_refant"
            in self.flagging
        ):

            recommended = self.flagging["recommended_refant"]

            if recommended:

                return self.validate_reference_antenna(recommended)

        # Priority 2:
        # inspection report

        refant = self.report.get(
            "reference_antenna",
            ""
        )

        if isinstance(refant, dict):

            refant = refant.get("recommended", "")

        if refant:

            return self.validate_reference_antenna(refant)

        raise RuntimeError(
            "No usable reference antenna was selected. Run the flagging "
            "stage to obtain a recommendation or pass --refant explicitly."
        )

    # ---------------------------------------------------------

    def validate_reference_antenna(self, antenna):

        antenna = str(antenna).strip()

        if not antenna or antenna.upper() == "UNKNOWN":

            raise RuntimeError(
                "No usable reference antenna was selected. Run the "
                "flagging stage to obtain a recommendation or pass "
                "--refant explicitly."
            )

        antenna_file = self.outdir / "inspection" / "antennas.json"

        if antenna_file.exists():

            with open(antenna_file) as f:

                known = {
                    item.get("name")
                    for item in json.load(f)
                    if item.get("name")
                }

            if known and antenna not in known:

                raise RuntimeError(
                    f"Reference antenna '{antenna}' does not exist in "
                    "the inspection report."
                )

        return antenna

    # ---------------------------------------------------------

    def known_fields(self):

        fields_file = self.outdir / "inspection" / "fields.json"

        if fields_file.exists():

            with open(fields_file) as f:

                return {
                    item.get("name")
                    for item in json.load(f)
                    if item.get("name")
                }

        groups = self.report.get("fields", {})

        return {
            name
            for names in groups.values()
            if isinstance(names, list)
            for name in names
        }

    # ---------------------------------------------------------

    def validate_field_name(self, field, role):

        known = self.known_fields()

        if known and field not in known:

            raise RuntimeError(
                f"Configured {role} field '{field}' does not exist in "
                "the inspection report."
            )

    # ---------------------------------------------------------

    def validate_field_list(self, fields, role):

        normalized = list(dict.fromkeys(fields or []))

        for field in normalized:

            self.validate_field_name(field, role)

        return normalized

    # ---------------------------------------------------------

    @staticmethod
    def require_unique_field(fields,
                             role,
                             option,
                             allow_multiple=False):

        fields = list(dict.fromkeys(fields or []))

        if len(fields) == 0:

            raise RuntimeError(
                f"No {role} is available in the inspection report. "
                f"Pass {option} explicitly."
            )

        if len(fields) > 1 and not allow_multiple:

            raise RuntimeError(
                f"Multiple {role}s are available: {', '.join(fields)}. "
                f"Pass {option} explicitly."
            )

        return fields if allow_multiple else fields[0]

    # =========================================================
    # CHECKPOINTS
    # =========================================================

    def caltable_qa_path(self, filename):

        return self.caldir / f"{filename}.qa.json"

    # ---------------------------------------------------------

    def write_caltable_qa(self, filename, qa):

        outfile = self.caltable_qa_path(filename)
        outfile.parent.mkdir(parents=True, exist_ok=True)
        temporary = outfile.with_name(
            f".{outfile.name}.tmp.{os.getpid()}"
        )

        with open(temporary, "w") as handle:
            json.dump(qa, handle, indent=2)
            handle.flush()
            os.fsync(handle.fileno())

        os.replace(temporary, outfile)

    # ---------------------------------------------------------

    def checkpoint_exists(self,
                          stage):

        donefile = (
            self.caldir /
            f"{stage}.done"
        )

        return donefile.exists()

    # ---------------------------------------------------------

    def mark_checkpoint(self,
                        stage):

        donefile = (
            self.caldir /
            f"{stage}.done"
        )

        donefile.touch()

    # ---------------------------------------------------------

    def caltable_status(self, filename):

        path = self.caldir / filename

        qa = {
            "schema_version": 2,
            "table": str(path),
            "status": "failed",
            "reason": "",
        }

        def finish(status, reason=""):
            qa["status"] = status
            qa["reason"] = reason
            self.write_caltable_qa(filename, qa)
            return status == "passed", reason

        try:
            max_flagged_fraction = float(
                self.task_config.get(
                    "max_caltable_flagged_fraction",
                    0.50,
                )
            )
        except (TypeError, ValueError):
            return finish(
                "failed",
                "max_caltable_flagged_fraction is not numeric",
            )

        qa["max_flagged_solution_fraction"] = max_flagged_fraction

        if (
            not math.isfinite(max_flagged_fraction) or
            not 0.0 <= max_flagged_fraction <= 1.0
        ):
            return finish(
                "failed",
                "max_caltable_flagged_fraction must be between 0 and 1",
            )

        try:
            max_group_flagged_fraction = float(
                self.task_config.get(
                    "max_caltable_group_flagged_fraction",
                    0.80,
                )
            )
        except (TypeError, ValueError):
            return finish(
                "failed",
                "max_caltable_group_flagged_fraction is not numeric",
            )

        qa["max_group_flagged_solution_fraction"] = (
            max_group_flagged_fraction
        )

        if (
            not math.isfinite(max_group_flagged_fraction) or
            not 0.0 <= max_group_flagged_fraction <= 1.0
        ):
            return finish(
                "failed",
                "max_caltable_group_flagged_fraction must be between 0 and 1",
            )

        if not path.is_dir():

            return finish("failed", "table directory is missing")

        if not (path / "table.dat").is_file():

            return finish("failed", "table.dat is missing")

        opened = False

        try:

            tb.open(str(path))
            opened = True
            rows = int(tb.nrows())
            qa["rows"] = rows
            columns = set(tb.colnames())

            parameter_column = (
                "CPARAM" if "CPARAM" in columns
                else "FPARAM" if "FPARAM" in columns
                else None
            )
            qa["parameter_column"] = parameter_column

            if rows < 1:
                return finish(
                    "failed",
                    "table contains no solution rows",
                )

            if parameter_column is None:
                return finish(
                    "failed",
                    "table has neither CPARAM nor FPARAM",
                )

            if "FLAG" not in columns:
                return finish(
                    "failed",
                    "table has no FLAG column",
                )

            parameters = tb.getcol(parameter_column)
            flags = tb.getcol("FLAG")

            group_columns = {}
            group_column_status = {}
            for group_column in [
                "ANTENNA1",
                "SPECTRAL_WINDOW_ID",
                "FIELD_ID",
                "OBSERVATION_ID",
            ]:
                if group_column not in columns:
                    group_column_status[group_column] = {
                        "status": "missing",
                    }
                    continue

                group_values = list(
                    _iter_casa_array_values(tb.getcol(group_column))
                )
                if len(group_values) != rows:
                    group_column_status[group_column] = {
                        "status": "skipped",
                        "reason": (
                            f"column has {len(group_values)} values for "
                            f"{rows} table rows"
                        ),
                    }
                    continue

                group_columns[group_column] = group_values
                group_column_status[group_column] = {
                    "status": "evaluated",
                }

            qa["group_column_status"] = group_column_status

        except Exception as exc:

            return finish(
                "failed",
                f"CASA cannot read the table ({exc})",
            )

        finally:

            if opened:

                tb.close()

        ignored_group_rows = set()
        if "ANTENNA1" in group_columns and getattr(self, "ms", None):
            try:
                active_antenna_ids = self.active_measurement_set_antenna_ids()
            except Exception as exc:
                return finish(
                    "failed",
                    "could not determine participating Measurement Set "
                    f"antennas ({exc})",
                )

            excluded_nonparticipating_ids = sorted({
                int(value)
                for value in group_columns["ANTENNA1"]
                if int(value) not in active_antenna_ids
            })
            try:
                excluded_solver_ids = self.solver_excluded_antenna_ids()
            except Exception as exc:
                return finish(
                    "failed",
                    f"could not validate excluded solver antennas ({exc})",
                )
            excluded_set = (
                set(excluded_nonparticipating_ids) | excluded_solver_ids
            )
            ignored_group_rows.update(
                row_index
                for row_index, value in enumerate(
                    group_columns["ANTENNA1"]
                )
                if int(value) in excluded_set
            )
            qa["excluded_nonparticipating_antenna_ids"] = (
                excluded_nonparticipating_ids
            )
            qa["excluded_solver_antenna_ids"] = sorted(excluded_solver_ids)
            if excluded_nonparticipating_ids:
                qa["group_column_status"]["ANTENNA1"].update({
                    "excluded_nonparticipating_values": (
                        excluded_nonparticipating_ids
                    ),
                    "reason": (
                        "CASA emitted placeholder solution rows for "
                        "antenna IDs absent from all MAIN-table baselines"
                    ),
                })
            if excluded_solver_ids:
                qa["group_column_status"]["ANTENNA1"].update({
                    "excluded_solver_values": sorted(excluded_solver_ids),
                    "solver_exclusion_reason": (
                        "antennas were explicitly excluded from every "
                        "calibration solve"
                    ),
                })

        qa["rows"] = rows
        solution_qa = summarize_caltable_solutions(
            parameter_column,
            parameters,
            flags,
            max_flagged_fraction,
            row_count=rows,
            group_columns=group_columns,
            max_group_flagged_fraction=max_group_flagged_fraction,
            ignored_group_row_indices=ignored_group_rows,
        )
        qa.update(solution_qa)

        return finish(
            solution_qa["status"],
            solution_qa["reason"],
        )

    # ---------------------------------------------------------

    def active_measurement_set_antenna_ids(self):
        """Return antenna IDs that occur in at least one MAIN-table row.

        CASA gain solvers can publish fully flagged placeholder rows for every
        ANTENNA-subtable entry, including antennas that never occur in the
        Measurement Set. Those placeholders are not failed participating
        antennas and must not trip per-antenna solution coverage QA.
        """

        cached = getattr(self, "_active_ms_antenna_ids", None)
        if cached is not None:
            return set(cached)

        tool = table()
        opened = False
        try:
            tool.open(str(self.ms))
            opened = True
            if int(tool.nrows()) < 1:
                raise RuntimeError("MAIN table has no rows")

            columns = set(tool.colnames())
            missing = {
                column for column in ("ANTENNA1", "ANTENNA2")
                if column not in columns
            }
            if missing:
                raise RuntimeError(
                    "MAIN table is missing " + ", ".join(sorted(missing))
                )

            active = set()
            for column in ("ANTENNA1", "ANTENNA2"):
                active.update(
                    int(value)
                    for value in _iter_casa_array_values(
                        tool.getcol(column)
                    )
                )
        finally:
            if opened:
                tool.close()

        if not active:
            raise RuntimeError("MAIN table contains no antenna IDs")

        self._active_ms_antenna_ids = frozenset(active)
        return set(active)

    # ---------------------------------------------------------

    def archive_uncheckpointed_caltable(self, filename):

        path = self.caldir / filename

        if not path.exists():

            return

        index = 1

        while True:

            archived = path.with_name(
                f"{path.name}.incomplete.{index}"
            )

            if not archived.exists():

                break

            index += 1

        path.rename(archived)
        print(
            "Preserved uncheckpointed calibration table as "
            f"{archived}"
        )

    # ---------------------------------------------------------

    def caltable_checkpoint_ready(self, stage, filename):

        if not self.checkpoint_exists(stage):

            self.archive_uncheckpointed_caltable(filename)
            return False

        valid, reason = self.caltable_status(filename)

        if valid:

            return True

        print(
            f"WARNING: {stage} checkpoint is stale because "
            f"{filename} is unusable ({reason}); recomputing it."
        )
        (self.caldir / f"{stage}.done").unlink()
        self.archive_uncheckpointed_caltable(filename)

        return False

    # ---------------------------------------------------------

    def require_caltable(self, filename):

        valid, reason = self.caltable_status(filename)

        if not valid:

            raise RuntimeError(
                f"CASA did not produce a usable calibration table "
                f"{self.caldir / filename}: {reason}"
            )

    # =========================================================

    def print_summary(self):

        print("\n")

        print("="*60)

        print("Calibration Summary")

        print("="*60)

        print(
            f"MS                : {self.ms}"
        )

        print(
            f"Band              : {self.band}"
        )

        print(
            "Flux Calibrators  : "
            + ", ".join(self.get_flux_calibrators())
        )

        print(
            "Bandpass Cal      : "
            f"{self.get_bandpass_calibrator()}"
        )

        print(
            "Phase Calibrators : "
            + ", ".join(self.get_phase_calibrators())
        )

        print(
            "Target/Phase Map  : "
            + ", ".join(
                f"{target}:{phase}"
                for target, phase in self.get_target_phase_map().items()
            )
        )

        print(
            f"Reference Antenna : "
            f"{self.get_reference_antenna()}"
        )

        print(
            f"Targets           : "
            f"{len(self.get_targets())}"
        )

        print(
            f"Bad Antennas      : "
            f"{self.flagging['bad_antennas']}"
        )

        print("="*60)

        print("\n")

    # =========================================================
    # CALIBRATION TASKS
    # =========================================================

    def hanning_ms_path(self):

        return self.caldir / "hann.ms"

    # ---------------------------------------------------------

    def select_hanning_ms_for_restart(self):

        if not self.band_config["hanningsmooth"]:
            return True

        hann_ms = self.hanning_ms_path()

        if not (
            hann_ms.is_dir() and
            (hann_ms / "table.dat").is_file()
        ):

            print(
                "WARNING: setjy checkpoint is stale because the "
                "Hanning-smoothed Measurement Set is missing or invalid; "
                "recomputing it."
            )
            (self.caldir / "setjy.done").unlink()

            if hann_ms.exists():

                index = 1

                while True:

                    archived = hann_ms.with_name(
                        f"{hann_ms.name}.incomplete.{index}"
                    )

                    if not archived.exists():

                        break

                    index += 1

                hann_ms.rename(archived)

            self.ms = self.input_ms
            return False

        self.ms = str(hann_ms)

        return True

    # ---------------------------------------------------------

    def expected_spw_ids(self):

        spws = self.report.get("spws", {})
        expected = []

        if isinstance(spws, dict):
            for key, value in spws.items():
                spw_id = value.get("id", key) if isinstance(value, dict) else key
                normalized = _normalized_spw_id(spw_id)
                if normalized is not None:
                    expected.append(normalized)

        if not expected:
            count = self.report.get("nspws")
            try:
                expected = [str(index) for index in range(int(count))]
            except (TypeError, ValueError, OverflowError):
                pass

        return sorted(set(expected))

    # ---------------------------------------------------------

    def flux_model_provenance_path(self):

        return self.caldir / "flux_model_provenance.json"

    # ---------------------------------------------------------

    def flux_model_identity(self, fluxcal, bandpasscal):

        fluxcals = self.get_flux_calibrators()

        identity = {
            "schema_version": 2,
            "measurement_set": str(
                Path(self.input_ms).expanduser().resolve()
            ),
            "telescope": self.telescope,
            "band": self.band,
            "primary_field": str(fluxcal),
            "flux_fields": [str(field) for field in fluxcals],
            "bandpass_field": str(bandpasscal),
            "bandpass_is_primary": str(bandpasscal) == str(fluxcal),
            "pipeline_recipe_fingerprint": self.recipe_fingerprint,
        }

        if self.flux_model:
            model_path = Path(self.flux_model).expanduser().resolve()
            identity.update({
                "mode": "external-ft",
                "flux_model_type": self.flux_model_type,
                "flux_model_path": str(model_path),
                "flux_model_sha256": sha256_file_or_directory(model_path),
                "distinct_bandpass_standard": (
                    self.task_config["flux_standard"]
                    if str(bandpasscal) != str(fluxcal)
                    else None
                ),
            })
        else:
            identity.update({
                "mode": "builtin-setjy",
                "flux_standard": self.task_config["flux_standard"],
            })

        manual = self.manual_flux_models()
        if manual:
            identity["calibrator_flux_models"] = {
                field: manual[field] for field in sorted(manual)
            }

        return identity

    def manual_flux_models(self):
        """Calibrators whose flux density the configuration states outright."""

        return dict(self.task_config.get("calibrator_flux_models") or {})

    def setjy_field(self, field, expected_spws):
        """Model one calibrator: its stated flux if it has one, else the standard.

        Returns the provenance record of the application. A stated model is
        CASA setjy's manual one: Stokes I at the reference frequency and a
        spectral index (and curvature) about it, scaled channel by channel.
        """

        model = self.manual_flux_models().get(str(field))
        if model is None:
            result = setjy(
                vis=self.ms,
                field=field,
                standard=self.task_config["flux_standard"],
                scalebychan=True,
                usescratch=True,
            )
            return {
                "field": str(field),
                "method": "setjy",
                "standard": self.task_config["flux_standard"],
                "flux_qa": validate_setjy_flux_result(
                    result,
                    expected_spws=expected_spws,
                ),
            }

        result = setjy(
            vis=self.ms,
            field=field,
            standard="manual",
            fluxdensity=[float(model["flux_jy"]), 0.0, 0.0, 0.0],
            spix=[float(term) for term in model.get("spix") or []] or [0.0],
            reffreq=f"{float(model['reffreq_mhz']):g}MHz",
            scalebychan=True,
            usescratch=True,
        )
        try:
            flux_qa = validate_setjy_flux_result(
                result,
                expected_spws=expected_spws,
            )
        except RuntimeError as exc:
            # The model column is checked either way; what setjy returns
            # for a manual model is only a report of it.
            flux_qa = {"status": "not-returned", "reason": str(exc)}
        return {
            "field": str(field),
            "method": "setjy-manual",
            "model": model,
            "flux_qa": flux_qa,
        }

    # ---------------------------------------------------------

    def load_flux_model_provenance(self):

        path = self.flux_model_provenance_path()
        if not path.is_file():
            return None

        try:
            with open(path) as handle:
                value = json.load(handle)
        except (OSError, ValueError, TypeError):
            return None

        return value if isinstance(value, dict) else None

    # ---------------------------------------------------------

    def flux_model_checkpoint_matches(self, identity):

        stored = self.load_flux_model_provenance()
        if not (
            stored and
            stored.get("status") == "completed" and
            isinstance(stored.get("identity"), dict) and
            stored.get("model_data_qa")
        ):
            return False

        def semantic_identity(value, ignore_recipe=False):
            value = dict(value)
            value.pop("pipeline_provenance_fingerprint", None)
            if ignore_recipe:
                value.pop("pipeline_recipe_fingerprint", None)
            return value

        legacy_identity = not stored["identity"].get(
            "pipeline_recipe_fingerprint"
        )
        if semantic_identity(
            stored["identity"],
            ignore_recipe=legacy_identity,
        ) != semantic_identity(
            identity,
            ignore_recipe=legacy_identity,
        ):
            return False

        # Migrate a checkpoint made before implementation and science-recipe
        # provenance were separated. The semantic setjy/ft identity is exact,
        # including model content hash, fields, band, MS, and flux standard.
        if stored["identity"] != identity:
            history = stored.setdefault("identity_migration_history", [])
            history.append({
                "migrated_utc": datetime.now(timezone.utc).isoformat(),
                "previous_identity": stored["identity"],
                "reason": (
                    "implementation provenance separated from science recipe"
                ),
            })
            stored["identity"] = identity
            self.write_atomic_json(
                self.flux_model_provenance_path(),
                stored,
            )

        return True

    # ---------------------------------------------------------

    def field_ids_for_name(self, field):

        field_text = str(field).strip()
        fields_file = self.outdir / "inspection" / "fields.json"
        entries = []

        if fields_file.is_file():
            try:
                with open(fields_file) as handle:
                    entries = json.load(handle)
            except (OSError, ValueError, TypeError):
                entries = []

        if not isinstance(entries, list) or not entries:
            field_table = Path(self.ms) / "FIELD"
            opened = False
            try:
                tb.open(str(field_table))
                opened = True
                names = list(_iter_casa_array_values(tb.getcol("NAME")))
                entries = [
                    {"id": index, "name": str(name)}
                    for index, name in enumerate(names)
                ]
            except Exception as exc:
                raise RuntimeError(
                    "Could not resolve the selected calibrator field ID"
                ) from exc
            finally:
                if opened:
                    tb.close()

        matches = []
        for entry in entries:
            if not isinstance(entry, dict) or "id" not in entry:
                continue
            if str(entry.get("name", "")) == field_text:
                matches.append(int(entry["id"]))

        if not matches and re.fullmatch(r"\d+", field_text):
            numeric = int(field_text)
            known_ids = {
                int(entry["id"])
                for entry in entries
                if isinstance(entry, dict) and "id" in entry
            }
            if numeric in known_ids:
                matches.append(numeric)

        matches = sorted(set(matches))
        if not matches:
            raise RuntimeError(
                f"Could not resolve calibrator field {field!r} to FIELD_ID"
            )

        return matches

    # ---------------------------------------------------------

    def validate_selected_model_data(self, field):

        field_ids = self.field_ids_for_name(field)
        selection = None
        opened = False

        try:
            tb.open(str(self.ms))
            opened = True
            if "MODEL_DATA" not in set(tb.colnames()):
                raise RuntimeError(
                    "Measurement Set has no physical MODEL_DATA column"
                )

            expression = "FIELD_ID IN [" + ",".join(
                str(field_id) for field_id in field_ids
            ) + "]"
            selection = tb.query(
                expression,
                columns="MODEL_DATA",
            )
            selected_rows = int(selection.nrows())
            if selected_rows < 1:
                raise RuntimeError(
                    f"MODEL_DATA selection for field {field!r} has no rows"
                )

            def read_rows(start, count):
                try:
                    return selection.getcol("MODEL_DATA", start, count)
                except Exception:
                    getvarcol = getattr(selection, "getvarcol", None)
                    if not callable(getvarcol):
                        raise
                    return list(
                        getvarcol("MODEL_DATA", start, count).values()
                    )

            def values_in(item):
                size = getattr(item, "size", None)
                if size is not None:
                    return int(size)
                return sum(1 for _ in _iter_casa_array_values(item))

            # Size blocks from a small probe. getcol returns one array whose
            # last axis is the row axis; getvarcol returns one array per row.
            # Taking the largest row keeps blocks bounded when spectral
            # windows of different widths interleave in time.
            probe_rows = min(64, selected_rows)
            probe = read_rows(0, probe_rows)
            if isinstance(probe, list):
                values_per_row = max(
                    (values_in(row) for row in probe),
                    default=1,
                )
            else:
                values_per_row = values_in(probe) // probe_rows
            values_per_row = max(1, values_per_row)

            rows_per_block = max(
                1, MODEL_DATA_CHECK_BLOCK_VALUES // values_per_row
            )
            row_budget = max(1, MODEL_DATA_CHECK_MAX_VALUES // values_per_row)

            if selected_rows <= row_budget:
                starts = list(range(0, selected_rows, rows_per_block))
                coverage = "complete"
            else:
                # Evenly spaced blocks across the whole selection, always
                # including the first and last rows.
                blocks = max(2, -(-row_budget // rows_per_block))
                last_start = max(0, selected_rows - rows_per_block)
                starts = sorted({
                    round(index * last_start / (blocks - 1))
                    for index in range(blocks)
                })
                coverage = "sampled"

            counts = {
                "total_values": 0,
                "finite_values": 0,
                "nonfinite_values": 0,
                "nonzero_values": 0,
            }
            rows_checked = 0
            for start in starts:
                count = min(rows_per_block, selected_rows - start)
                block = _count_model_data_values(read_rows(start, count))
                for key in counts:
                    counts[key] += block[key]
                rows_checked += count

        finally:
            if selection is not None:
                selection.close()
            if opened:
                tb.close()

        qa = summarize_model_data(None, counts=counts)
        qa.update({
            "field": str(field),
            "field_ids": field_ids,
            "selected_rows": selected_rows,
            "rows_checked": rows_checked,
            "coverage": coverage,
            "values_per_row": values_per_row,
            "value_budget": MODEL_DATA_CHECK_MAX_VALUES,
        })
        if qa["status"] != "passed":
            raise RuntimeError(
                f"MODEL_DATA validation failed for field {field!r}: "
                f"{qa['reason']}"
            )

        return qa

    # ---------------------------------------------------------

    def invalidate_model_dependent_checkpoints(self, reason):

        print(
            "WARNING: Flux-model provenance changed or became invalid; "
            f"invalidating dependent calibration checkpoints ({reason})."
        )

        stages = [
            "setjy",
            "delay",
            "prebp",
            "bandpass",
            "phase",
            "phase_scan",
            "amp",
            "fluxscale",
            "applycal",
            "split",
            "qa",
            "calibration",
            f"{self.telescope}_postcal_flagging",
            f"{self.telescope}_round2_prepared",
            f"{self.telescope}_round2_complete",
        ]
        for stage in stages:
            (self.caldir / f"{stage}.done").unlink(missing_ok=True)

        round1 = self.caldir / "round1"
        if round1.exists():
            index = 1
            while True:
                stale = self.caldir / f"round1.stale.{index}"
                if not stale.exists():
                    break
                index += 1
            round1.rename(stale)

        if self.splitdir.exists() and any(self.splitdir.iterdir()):
            index = 1
            while True:
                stale = self.caldir / f"calibrated.stale.{index}"
                if not stale.exists():
                    break
                index += 1
            self.splitdir.rename(stale)
            self.splitdir.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------

    def validate_flux_model_policy(self, fluxcal):

        normalized_fluxcal = re.sub(
            r"[^A-Z0-9]+",
            "",
            str(fluxcal).upper(),
        )
        if self.telescope == "meerkat" and not self.flux_model:
            is_j0408 = (
                "0408" in normalized_fluxcal and
                (
                    normalized_fluxcal.endswith("65") or
                    normalized_fluxcal.endswith("6545")
                )
            )
            is_sr_primary = any(
                token in normalized_fluxcal
                for token in ("1934638", "19396342")
            )

            if self.band == "MEERKAT_UHF":
                raise RuntimeError(
                    "MeerKAT UHF calibration requires a validated, "
                    "versioned calibrator model. Pass --flux-model and "
                    "--flux-model-type."
                )
            if is_j0408:
                raise RuntimeError(
                    "MeerKAT J0408-6545/PKS0408-65 requires its documented "
                    "broadband model. Pass --flux-model and "
                    "--flux-model-type."
                )
            if not is_sr_primary and not self.flux_standard_override:
                raise RuntimeError(
                    f"MeerKAT flux calibrator {fluxcal!r} is not one of the "
                    "built-in Stevens-Reynolds primary aliases. Pass a "
                    "validated --flux-model or an explicit --flux-standard."
                )

        if (
            not self.flux_model and
            str(self.task_config["flux_standard"]).strip().lower() == "manual"
        ):
            raise RuntimeError(
                "CASA setjy standard='manual' requires explicit spectral "
                "flux parameters; this pipeline instead requires a validated "
                "external --flux-model for manual models."
            )

    # ---------------------------------------------------------

    def run_setjy(self):
        fluxcals = self.get_flux_calibrators()
        fluxcal = fluxcals[0]
        bandpasscal = self.get_bandpass_calibrator()
        for selected_fluxcal in fluxcals:
            self.validate_flux_model_policy(selected_fluxcal)
        identity = self.flux_model_identity(fluxcal, bandpasscal)

        setjy_checkpoint = self.checkpoint_exists("setjy")
        if not setjy_checkpoint:
            downstream_stages = [
                "delay",
                "prebp",
                "bandpass",
                "phase",
                "phase_scan",
                "amp",
                "fluxscale",
                "applycal",
                f"{self.telescope}_round2_complete",
                "split",
                "calibration",
            ]
            if (
                any(self.checkpoint_exists(stage) for stage in downstream_stages)
                and not self.flux_model_checkpoint_matches(identity)
            ):
                self.invalidate_model_dependent_checkpoints(
                    "setjy marker is absent and downstream products do not "
                    "have matching flux-model provenance"
                )

        if setjy_checkpoint:
            hanning_ready = self.select_hanning_ms_for_restart()
            if hanning_ready and self.flux_model_checkpoint_matches(identity):
                try:
                    for selected_fluxcal in fluxcals:
                        self.validate_selected_model_data(selected_fluxcal)
                    if str(bandpasscal) not in {
                        str(field)
                        for field in fluxcals
                    }:
                        self.validate_selected_model_data(bandpasscal)
                except Exception as exc:
                    self.invalidate_model_dependent_checkpoints(str(exc))
                else:
                    print("Flux-model application already completed.")
                    return
            else:
                reason = (
                    "Hanning-smoothed MS is unavailable"
                    if not hanning_ready
                    else "stored flux-model identity does not match"
                )
                self.invalidate_model_dependent_checkpoints(reason)

        print(f"Applying flux model to {fluxcal}")

        hann_ms = str(self.hanning_ms_path())

        if self.band_config["hanningsmooth"]:

            hann_path = Path(hann_ms)

            if (
                hann_path.exists() and
                not (hann_path / "table.dat").is_file()
            ):

                index = 1

                while True:

                    archived = hann_path.with_name(
                        f"{hann_path.name}.incomplete.{index}"
                    )

                    if not archived.exists():

                        break

                    index += 1

                hann_path.rename(archived)
                print(
                    "Preserved incomplete Hanning MS as "
                    f"{archived}"
                )

            if not Path(hann_ms).exists():

                print("Running Hanning smoothing.")

                hanningsmooth(
                    vis=self.ms,
                    outputvis=hann_ms,
                    datacolumn="DATA"
                )

            else:

                print(
                    "Using existing "
                    "Hanning-smoothed MS."
                )

            if not (
                Path(hann_ms).is_dir() and
                (Path(hann_ms) / "table.dat").is_file()
            ):

                raise RuntimeError(
                    "Hanning smoothing did not produce a usable "
                    f"Measurement Set: {hann_ms}"
                )

            self.ms = hann_ms

        provenance = {
            "schema_version": 1,
            "status": "running",
            "identity": identity,
            "applications": [],
            "model_data_qa": {},
        }
        self.write_atomic_json(
            self.flux_model_provenance_path(),
            provenance,
        )

        expected_spws = self.expected_spw_ids()

        try:
            modeled_fields = set()

            if self.flux_model:
                model_path = Path(self.flux_model).expanduser().resolve()
                ft_kwargs = {
                    "vis": self.ms,
                    "field": fluxcal,
                    "usescratch": True,
                    "incremental": False,
                }
                if self.flux_model_type == "component-list":
                    ft_kwargs["complist"] = str(model_path)
                else:
                    ft_kwargs["model"] = str(model_path)

                ft(**ft_kwargs)
                provenance["applications"].append({
                    "field": str(fluxcal),
                    "method": "ft",
                    "flux_model_type": self.flux_model_type,
                    "model_path": str(model_path),
                    "model_sha256": identity["flux_model_sha256"],
                })
                modeled_fields.add(str(fluxcal))
            else:
                provenance["applications"].append(
                    self.setjy_field(fluxcal, expected_spws)
                )
                modeled_fields.add(str(fluxcal))

            provenance["model_data_qa"][str(fluxcal)] = (
                self.validate_selected_model_data(fluxcal)
            )

            for extra_fluxcal in fluxcals[1:]:
                try:
                    extra_application = self.setjy_field(
                        extra_fluxcal, expected_spws
                    )
                    extra_model_qa = self.validate_selected_model_data(
                        extra_fluxcal
                    )
                except Exception as exc:
                    raise RuntimeError(
                        "An additional flux calibrator could not receive a "
                        "validated broadband model. Pass a supported "
                        "--flux-standard/--flux-model combination or reduce "
                        "the explicit --flux-calibrator list "
                        f"({extra_fluxcal!r}: {exc})."
                    ) from exc

                extra_application["role"] = "additional-flux-calibrator"
                provenance["applications"].append(extra_application)
                provenance["model_data_qa"][str(extra_fluxcal)] = (
                    extra_model_qa
                )
                modeled_fields.add(str(extra_fluxcal))

            if str(bandpasscal) not in modeled_fields:
                try:
                    bp_result = setjy(
                        vis=self.ms,
                        field=bandpasscal,
                        standard=self.task_config["flux_standard"],
                        scalebychan=True,
                        usescratch=True,
                    )
                    bp_flux_qa = validate_setjy_flux_result(
                        bp_result,
                        expected_spws=expected_spws,
                    )
                    bp_model_qa = self.validate_selected_model_data(
                        bandpasscal
                    )
                except Exception as exc:
                    raise RuntimeError(
                        "The distinct bandpass calibrator could not receive "
                        "a validated broadband built-in CASA model. Set the "
                        "bandpass calibrator equal to the modeled primary "
                        "flux calibrator, or add a future explicit bandpass "
                        "model facility; refusing to solve K/preBP/B against "
                        f"an unmodeled field ({exc})."
                    ) from exc

                provenance["applications"].append({
                    "field": str(bandpasscal),
                    "method": "setjy",
                    "role": "distinct-bandpass-calibrator",
                    "standard": self.task_config["flux_standard"],
                    "flux_qa": bp_flux_qa,
                })
                provenance["model_data_qa"][str(bandpasscal)] = bp_model_qa

        except Exception as exc:
            provenance["status"] = "failed"
            provenance["error"] = str(exc)
            self.write_atomic_json(
                self.flux_model_provenance_path(),
                provenance,
            )
            raise

        provenance["status"] = "completed"
        self.write_atomic_json(
            self.flux_model_provenance_path(),
            provenance,
        )
        self.mark_checkpoint("setjy")

    # ---------------------------------------------------------

    def run_delaycal(self):

        if self.caltable_checkpoint_ready("delay", "delay.cal"):

            print("Delay calibration already completed.")
            return

        self.solve_delay_table()

        self.require_caltable("delay.cal")

        self.mark_checkpoint("delay")

    def solve_delay_table(self):
        """Solve delay.cal, quarantining an antenna it cannot use at all."""

        fluxcal = self.get_bandpass_calibrator()
        refant = self.get_reference_antenna()
        combine = self.band_config["combine_delay"]
        delay_solint = self.band_config["delay_solint"]

        def solve():
            print(
                f"Running delay calibration "
                f"(solint={delay_solint})"
            )
            gaincal(
                vis=self.ms,
                caltable=str(self.caldir / "delay.cal"),
                field=fluxcal,
                gaintype="K",
                antenna=self.get_antenna_selection(),
                refant=refant,
                refantmode=self.task_config.get("refantmode", "strict"),
                combine=combine,
                solint=delay_solint,
                minsnr=self.task_config["delay_minsnr"],
                parang=self.task_config["parang"],
            )

        solve()
        valid, reason = self.caltable_status("delay.cal")
        if not valid:
            failed_ids = self.failed_delay_antenna_ids("delay.cal")
            quarantined = self.quarantine_failed_delay_antennas(failed_ids)
            if quarantined:
                self.archive_uncheckpointed_caltable("delay.cal")
                solve()
            else:
                raise RuntimeError(
                    self.unsolvable_table_message("delay.cal", reason,
                                                  failed_ids)
                )

    # =========================================================
    # PRE-BANDPASS PHASE CALIBRATION
    # =========================================================

    def run_prebp(self):

        if self.caltable_checkpoint_ready("prebp", "prebp.cal"):

            print(
                "Pre-bandpass phase calibration already completed."
            )

            return

        self.solve_prebp_table()
        valid, reason = self.caltable_status("prebp.cal")
        if not valid:
            failed_ids = self.failed_delay_antenna_ids("prebp.cal")
            quarantined = self.quarantine_unsolvable_antennas(
                failed_ids,
                "prebp.cal",
                "no usable pre-bandpass phase solution",
            )
            if not quarantined:
                raise RuntimeError(
                    self.unsolvable_table_message("prebp.cal", reason,
                                                  failed_ids)
                )
            # The delay solve already included the antenna now excluded;
            # "excluded from every solve" means that one too.
            self.archive_uncheckpointed_caltable("prebp.cal")
            self.archive_uncheckpointed_caltable("delay.cal")
            self.solve_delay_table()
            self.require_caltable("delay.cal")
            self.solve_prebp_table()

        self.finish_prebp()

    def solve_prebp_table(self):

        fluxcal = self.get_bandpass_calibrator()

        refant = self.get_reference_antenna()

        antenna=self.get_antenna_selection()

        combine = self.band_config["combine_delay"]

        print(
            "Running pre-bandpass phase calibration."
        )

        gaincal(
            vis=self.ms,

            caltable=str(
                self.caldir /
                "prebp.cal"
            ),

            field=fluxcal,

            refant=refant,

            refantmode=self.task_config.get("refantmode", "strict"),

            antenna=antenna,

            gaintype=self.task_config.get("prebp_gaintype", "G"),

            calmode="p",

            combine=combine,

            solint="int",

            minsnr=self.task_config[
                "gain_minsnr"
            ],

            minblperant=self.task_config[
                "minblperant"
            ],

            solmode=self.task_config[
                "gain_solmode"
            ],

            parang=self.task_config[
                "parang"
            ],

            gaintable=[

                str(
                    self.caldir /
                    "delay.cal"
                )
            ]
        )

    def finish_prebp(self):

        self.require_caltable("prebp.cal")

        self.mark_checkpoint(
            "prebp"
        )

    # ---------------------------------------------------------

    def bandpass_quarantine_max_usable_fraction(self):
        """The usable share at or below which an antenna is excluded.

        An antenna that fails the per-antenna table check on its own (more
        than ``max_caltable_group_flagged_fraction`` of its bandpass
        solutions flagged) has no signal at the calibrator worth keeping.
        In the Band 5 data of 2021-02-22, C03 and C05 failed SNR 3 in every
        channel on the 27-minute 3C48 scan (C05 kept 1%, by chance), while
        every other antenna kept 94%, the rest being the band edges; on the
        scan's 9.7 settled minutes (finding 161) noise alone kept 10% of
        theirs, and the earlier fixed limit of 5% stopped the run on two
        dead antennas. The limits on how many are excluded still apply.
        """

        try:
            group_maximum = float(self.task_config.get(
                "max_caltable_group_flagged_fraction", 0.80))
        except (TypeError, ValueError):
            group_maximum = 0.80
        return min(1.0, max(0.0, 1.0 - group_maximum))

    def run_bandpass(self):

        if self.caltable_checkpoint_ready("bandpass", "bandpass.cal"):

            print("Bandpass already completed.")
            return

        self.solve_bandpass_table()
        valid, reason = self.caltable_status("bandpass.cal")
        if not valid:
            usable_limit = self.bandpass_quarantine_max_usable_fraction()
            failed_ids = self.failed_delay_antenna_ids(
                "bandpass.cal",
                max_usable_fraction=usable_limit,
            )
            quarantined = self.quarantine_unsolvable_antennas(
                failed_ids,
                "bandpass.cal",
                "no usable bandpass solution (at most "
                f"{usable_limit:.0%} of channels usable)",
            )
            if not quarantined:
                raise RuntimeError(
                    self.unsolvable_table_message("bandpass.cal", reason,
                                                  failed_ids)
                )
            # Excluded from every solve: the delay and pre-bandpass tables
            # were solved with these antennas, so both are solved again.
            self.archive_uncheckpointed_caltable("bandpass.cal")
            self.archive_uncheckpointed_caltable("prebp.cal")
            self.archive_uncheckpointed_caltable("delay.cal")
            self.solve_delay_table()
            self.require_caltable("delay.cal")
            self.solve_prebp_table()
            self.require_caltable("prebp.cal")
            self.solve_bandpass_table()

        self.require_caltable("bandpass.cal")

        self.mark_checkpoint(
            "bandpass"
        )

    def solve_bandpass_table(self):

        fluxcal = self.get_bandpass_calibrator()

        refant = self.get_reference_antenna()

        antenna=self.get_antenna_selection()

        combine = self.band_config["combine_delay"]

        print(
            "Running bandpass calibration."
        )

        bandpass(
            vis=self.ms,

            caltable=str(
                self.caldir /
                "bandpass.cal"
            ),

            field=fluxcal,

            refant=refant,

            antenna=antenna,

            combine=combine,

            bandtype="B",

            fillgaps=self.task_config[
                "bandpass_fillgaps"
            ],

            solint=self.task_config[
                "bandpass_solint"
            ],

            solnorm=self.task_config[
                "bandpass_solnorm"
            ],

            minsnr=self.task_config[
                "gain_minsnr"
            ],

            minblperant=self.task_config[
                "minblperant"
            ],

            parang=self.task_config[
                "parang"
            ],

            gaintable=[
                str(
                    self.caldir /
                    "delay.cal"
                ),

                str(
                    self.caldir /
                    "prebp.cal"
                )
            ]
        )

    # ---------------------------------------------------------

    def run_phasecal(self):

        if self.caltable_checkpoint_ready("phase", "phase.cal"):

            print("Phase calibration already completed.")
            return

        phasecals = self.get_phase_calibrators()
        fluxcals = self.get_flux_calibrators()

        bandpasscal = (
            self.get_bandpass_calibrator()
        )

        fields = sorted(
            {
                field
                for field in (
                    fluxcals +
                    [bandpasscal] +
                    phasecals
                )
                if field
            }
        )

        refant = (
            self.get_reference_antenna()
        )

        antenna = (
            self.get_antenna_selection()
        )

        combine = (
            self.band_config[
                "combine_delay"
            ]
        )

        phase_solint = (
            self.band_config[
                "phase_solint"
            ]
        )

        def solve(gaintype):
            print(
                f"Running phase calibration "
                f"(solint={phase_solint}, gaintype={gaintype})"
            )

            gaincal(
                vis=self.ms,
                caltable=str(self.caldir / "phase.cal"),
                field=",".join(fields),
                antenna=antenna,
                refant=refant,
                refantmode=self.task_config.get("refantmode", "strict"),
                minblperant=self.task_config["minblperant"],
                gaintype=gaintype,
                calmode="p",
                combine="",
                solint=phase_solint,
                solmode=self.task_config["gain_solmode"],
                minsnr=self.task_config["gain_minsnr"],
                parang=self.task_config["parang"],
                gaintable=[
                    str(self.caldir / "delay.cal"),
                    str(self.caldir / "bandpass.cal"),
                ],
                gainfield=[bandpasscal, bandpasscal],
            )

        gaintype = self.calibrator_gaintype("phase")
        solve(gaintype)
        valid, reason = self.caltable_status("phase.cal")
        if not valid and gaintype == "T":
            # T fits one phase to both hands, which holds only while the R-L
            # phase stays where the bandpass left it. Where it drifts, the
            # two hands cancel and solutions fail SNR although each hand
            # alone is strong: in the Band 5 data of 2021-02-22, C11, S02,
            # S03 and the reference C00 kept 9-28% of their T solutions but
            # all of their per-hand pre-bandpass solutions. One phase per
            # hand (G) is the standard solve; it is used for the rest of the
            # run when it passes where T failed.
            self.archive_uncheckpointed_caltable("phase.cal")
            solve("G")
            valid_g, reason_g = self.caltable_status("phase.cal")
            if not valid_g:
                raise RuntimeError(
                    "CASA did not produce a usable calibration table "
                    f"{self.caldir / 'phase.cal'} with one phase for both "
                    f"hands (T: {reason}) or one per hand (G: {reason_g})."
                )
            self.record_gaintype_switch("phase.cal", reason)

        self.require_caltable("phase.cal")

        self.mark_checkpoint(
            "phase"
        )

    # ---------------------------------------------------------

    def gaintype_switch_path(self):

        return self.caldir / "calibrator_gaintype.json"

    def read_gaintype_switch(self):
        """This run's switch of the calibrator gain solves to G, if any."""

        try:
            record = json.loads(self.gaintype_switch_path().read_text())
        except (OSError, ValueError, TypeError):
            return None
        if not isinstance(record, dict) or record.get("gaintype") != "G":
            return None
        expected = str(Path(self.input_ms).expanduser().resolve())
        if record.get("measurement_set") != expected:
            return None
        return record

    def record_gaintype_switch(self, table, reason):
        record = {
            "schema_version": 1,
            "measurement_set": str(Path(self.input_ms).expanduser().resolve()),
            "gaintype": "G",
            "configured_gaintype": {
                "phase": self.task_config.get("phase_gaintype", "T"),
                "amp": self.task_config.get("amp_gaintype", "T"),
            },
            "table": table,
            "reason": (
                "one phase for both hands (T) failed the table check, one "
                "per hand (G) passed: " + str(reason)
            ),
            "created_utc": datetime.now(timezone.utc).isoformat(),
        }
        self.write_atomic_json(self.gaintype_switch_path(), record)
        print(
            "WARNING: the phase solve with one phase for both hands (T) "
            f"failed ({reason}); one phase per hand (G) passed. The "
            "calibrator gain solves of this run use G."
        )

    def calibrator_gaintype(self, kind):
        """Gaintype of the calibrator gain solves, 'phase' or 'amp'.

        The telescope's choice (T for uGMRT, for signal-to-noise) unless this
        run has switched to G because the two hands were not coherent; the
        switch is kept on disk, so a resumed run solves the same way.
        """

        if self.read_gaintype_switch():
            return "G"
        return self.task_config.get(f"{kind}_gaintype", "T")

    # ---------------------------------------------------------

    def run_phase_transfer_cal(self):
        """Solve the scan-averaged phase table transferred to the targets.

        ``phase.cal`` is deliberately short (tens of seconds) so that the
        amplitude solve sees phase-coherent data.  Interpolating that table
        across the minutes-long gap between calibrator scans transfers one
        short-solution noise realisation to the target instead of averaging
        it down.  The scan-averaged table solved here is the one applied to
        science targets, which is the standard calibrator-transfer recipe.
        """

        if self.caltable_checkpoint_ready("phase_scan", "phase_scan.cal"):

            print("Scan-averaged phase calibration already completed.")
            return

        phasecals = self.get_phase_calibrators()
        fluxcals = self.get_flux_calibrators()
        bandpasscal = self.get_bandpass_calibrator()

        fields = sorted({
            field
            for field in (fluxcals + [bandpasscal] + phasecals)
            if field
        })

        solint = str(
            self.task_config.get("target_phase_transfer_solint", "inf")
        )

        print(
            f"Running scan-averaged phase calibration (solint={solint})"
        )

        gaincal(
            vis=self.ms,
            caltable=str(self.caldir / "phase_scan.cal"),
            field=",".join(fields),
            antenna=self.get_antenna_selection(),
            refant=self.get_reference_antenna(),
            refantmode=self.task_config.get("refantmode", "strict"),
            minblperant=self.task_config["minblperant"],
            gaintype=self.calibrator_gaintype("phase"),
            calmode="p",
            combine="",
            solint=solint,
            solmode=self.task_config["gain_solmode"],
            minsnr=self.task_config["gain_minsnr"],
            parang=self.task_config["parang"],
            gaintable=[
                str(self.caldir / "delay.cal"),
                str(self.caldir / "bandpass.cal"),
            ],
            gainfield=[bandpasscal, bandpasscal],
        )

        self.require_caltable("phase_scan.cal")

        self.mark_checkpoint("phase_scan")

    # ---------------------------------------------------------

    def run_ampcal(self):

        if self.caltable_checkpoint_ready("amp", "amp.cal"):

            print("Amplitude calibration already completed.")
            return

        phasecals = self.get_phase_calibrators()
        fluxcals = self.get_flux_calibrators()

        bandpasscal = (
            self.get_bandpass_calibrator()
        )

        fields = sorted(
            {
                field
                for field in (
                    fluxcals +
                    [bandpasscal] +
                    phasecals
                )
                if field
            }
        )

        refant = (
            self.get_reference_antenna()
        )

        antenna = (
            self.get_antenna_selection()
        )

        combine = (
            self.band_config[
                "combine_delay"
            ]
        )

        amp_solint = (
            self.band_config[
                "amp_solint"
            ]
        )

        print(
            f"Running amplitude calibration "
            f"(solint={amp_solint})"
        )

        gaincal(
            vis=self.ms,

            caltable=str(
                self.caldir /
                "amp.cal"
            ),

            field=",".join(fields),

            antenna=antenna,

            refant=refant,

            refantmode=self.task_config.get("refantmode", "strict"),

            gaintype=self.calibrator_gaintype("amp"),

            calmode="ap",

            minblperant=self.task_config[
                "minblperant"
            ],

            combine="",

            solint=amp_solint,

            solmode=self.task_config[
                "gain_solmode"
            ],

            minsnr=self.task_config[
                "gain_minsnr"
            ],

            parang=self.task_config[
                "parang"
            ],

            gaintable=[

                str(
                    self.caldir /
                    "delay.cal"
                ),

                str(
                    self.caldir /
                    "bandpass.cal"
                ),

                str(
                    self.caldir /
                    "phase.cal"
                )
            ],

            gainfield=[bandpasscal, bandpasscal, ""]
        )

        self.require_caltable("amp.cal")

        self.mark_checkpoint(
            "amp"
        )

    # =========================================================
    # FLUX SCALE
    # =========================================================

    def run_fluxscale(self):

        if self.caltable_checkpoint_ready("fluxscale", "flux.cal"):

            print(
                "Fluxscale already completed."
            )

            return

        fluxcals = self.get_flux_calibrators()
        primary_fluxcal = fluxcals[0]
        phasecals = self.get_phase_calibrators()

        # Every field that received a validated absolute model in run_setjy
        # is a fluxscale *reference*.  A reference field must never also be
        # listed as a transfer field: CASA would then be asked to bootstrap a
        # calibrator from itself, which is undefined and silently corrupts
        # the derived flux-density scale of the remaining transfer fields.
        reference_fields = list(dict.fromkeys(fluxcals))
        transfer = sorted(
            {
                field
                for field in phasecals
                if field and field not in set(reference_fields)
            }
        )

        fluxtable = self.caldir / "flux.cal"
        amptable = self.caldir / "amp.cal"

        if not transfer:
            # Degenerate but legal setup: the phase calibrator is itself a
            # modelled flux calibrator, so amp.cal already carries the
            # absolute scale and there is nothing to bootstrap.  Publish it
            # under the name the apply/QA chain expects instead of calling
            # fluxscale with an empty transfer list.
            print(
                "No transfer field is distinct from the modelled flux "
                "calibrator(s); amp.cal already carries the absolute flux "
                "scale and is published as flux.cal."
            )
            if fluxtable.exists():
                self.archive_uncheckpointed_caltable("flux.cal")
            shutil.copytree(amptable, fluxtable)
            self.write_atomic_json(
                self.caldir / "fluxscale_report.json",
                {
                    "status": "not_required",
                    "reason": (
                        "every gain field is a modelled absolute-flux "
                        "reference"
                    ),
                    "reference": reference_fields,
                    "transfer": [],
                },
            )
        else:
            # A spectral fit needs more spectral windows than its order.
            # uGMRT/MeerKAT wideband setups are frequently a single wide SPW,
            # where an unclamped fitorder makes fluxscale fail or return a
            # meaningless spectrum for the transferred calibrator.
            configured_fitorder = int(
                self.task_config.get("fluxscale_fitorder", 1)
            )
            spw_count = len(self.expected_spw_ids())
            fitorder = max(0, min(configured_fitorder, max(0, spw_count - 1)))

            print(
                "Running fluxscale "
                f"(reference={','.join(reference_fields)}, "
                f"transfer={','.join(transfer)}, fitorder={fitorder})"
            )

            fluxscale_result = fluxscale(
                vis=self.ms,

                caltable=str(amptable),

                fluxtable=str(fluxtable),

                reference=",".join(reference_fields),

                transfer=transfer,

                fitorder=fitorder,
            )

            self.write_atomic_json(
                self.caldir / "fluxscale_report.json",
                {
                    "status": "completed",
                    "reference": reference_fields,
                    "transfer": transfer,
                    "configured_fitorder": configured_fitorder,
                    "effective_fitorder": fitorder,
                    "spectral_window_count": spw_count,
                    "result": to_builtin_json(fluxscale_result),
                },
            )

        self.require_caltable("flux.cal")

        self.mark_checkpoint(
            "fluxscale"
        )

    # =========================================================
    # APPLYCAL
    # =========================================================

    # ---------------------------------------------------------
    # TRANSFER TO TARGETS ACROSS AN ISOLATED FAILED SOLUTION
    # ---------------------------------------------------------
    #
    # applycal interpolates a time-dependent table linearly between the two
    # solutions that bracket each visibility, and flags the visibility when
    # either of them is flagged. One failed solution for an antenna at one
    # calibrator scan therefore removes that antenna from the science scans on
    # both sides of it. In the uGMRT Band 4 run that prompted this, five
    # antennas with one to three failed solutions out of twelve kept 6 to 7%
    # of their data on a target, the scan-by-scan flagged fraction stepped
    # down from 86% to 45% across the first ten calibrator scans, and a
    # target observed for three scans was imaged from 26% of its data.
    #
    # A calibrator's own scan stays flagged where its solution failed: that
    # failure is a statement about those data. A target between two good
    # solutions of the same antenna is a different case, and what every
    # interpolating calibration does there is use them. The tables applied
    # to targets are therefore copies in which a failed solution with good
    # ones close by is replaced from them. The originals are untouched and
    # remain what the calibrators receive.

    def target_transfer_bridge_settings(self):
        config = self.task_config
        return {
            "enabled": bool(config.get("target_transfer_bridge", False)),
            "max_gap_s": 60.0 * float(
                config.get("target_transfer_bridge_max_gap_minutes", 120.0)
            ),
            "min_valid_fraction": float(
                config.get("target_transfer_bridge_min_valid_fraction", 0.5)
            ),
        }

    @staticmethod
    def bridge_solution_series(times, gains, flags, max_gap_s,
                               min_valid_fraction=0.5):
        """Fill failed solutions of one antenna from its good neighbours.

        ``times``, ``gains`` and ``flags`` describe one antenna, field,
        spectral window, polarization and channel, in time order. A failed
        solution between two good ones no more than ``max_gap_s`` apart
        takes their linear interpolation, amplitude and phase separately and
        the phase along the shorter arc. One with a good solution on one side
        only, at the start or end of the series, takes that solution if it is
        within half of ``max_gap_s``. Only solutions that were good to begin
        with are interpolated from, and an antenna with fewer than
        ``min_valid_fraction`` of its solutions good is left alone: that is
        a failing antenna, not a failed solution.

        Returns the new gains, the new flags and one record per fill.
        """

        import numpy as np

        times = np.asarray(times, dtype=float)
        gains = np.array(gains, dtype=complex)
        flags = np.array(flags, dtype=bool)
        count = times.size
        valid = (
            ~flags & np.isfinite(gains.real) & np.isfinite(gains.imag) &
            (np.abs(gains) > 0)
        )
        good = int(np.count_nonzero(valid))
        if (
            count < 2 or good == 0 or good == count or
            good < min_valid_fraction * count
        ):
            return gains, flags, []

        valid_indices = np.flatnonzero(valid)
        original = gains.copy()
        fills = []
        for index in np.flatnonzero(~valid):
            before = valid_indices[valid_indices < index]
            after = valid_indices[valid_indices > index]
            if before.size and after.size:
                low, high = int(before[-1]), int(after[0])
                gap = float(times[high] - times[low])
                if not 0.0 <= gap <= max_gap_s:
                    continue
                weight = (
                    float(times[index] - times[low]) / gap if gap > 0 else 0.5
                )
                amplitude = (
                    (1.0 - weight) * abs(original[low]) +
                    weight * abs(original[high])
                )
                arc = np.angle(original[high] * np.conj(original[low]))
                phase = np.angle(original[low]) + weight * arc
                value = amplitude * np.exp(1j * phase)
                kind = "interpolated"
                span = gap
            else:
                neighbour = int(before[-1] if before.size else after[0])
                span = abs(float(times[index] - times[neighbour]))
                if span > 0.5 * max_gap_s:
                    continue
                value = original[neighbour]
                kind = "nearest"
            gains[index] = value
            flags[index] = False
            fills.append({
                "index": int(index),
                "time": float(times[index]),
                "kind": kind,
                "span_s": span,
            })
        return gains, flags, fills

    def bridge_transfer_table(self, source, destination):
        """Write a copy of a gain table with isolated failures bridged.

        Returns a record of what was filled, or None when the copy could not
        be made or read, in which case the caller applies the original.
        """

        import numpy as np

        settings = self.target_transfer_bridge_settings()
        source = Path(source)
        destination = Path(destination)
        if destination.exists():
            shutil.rmtree(destination)

        tool = table()
        opened = False
        try:
            shutil.copytree(source, destination)
            tool.open(str(destination), nomodify=False)
            opened = True
            times = np.asarray(tool.getcol("TIME"), dtype=float)
            antennas = np.asarray(tool.getcol("ANTENNA1"), dtype=int)
            fields = np.asarray(tool.getcol("FIELD_ID"), dtype=int)
            windows = np.asarray(
                tool.getcol("SPECTRAL_WINDOW_ID"), dtype=int
            )
            gains = np.array(tool.getcol("CPARAM"), dtype=complex)
            flags = np.array(tool.getcol("FLAG"), dtype=bool)
            if gains.ndim != 3 or gains.shape != flags.shape:
                raise ValueError(
                    f"unexpected solution shape {gains.shape}"
                )

            fills = []
            keys = sorted(set(zip(fields.tolist(), windows.tolist(),
                                  antennas.tolist())))
            for field, window, antenna in keys:
                rows = np.flatnonzero(
                    (fields == field) & (windows == window) &
                    (antennas == antenna)
                )
                rows = rows[np.argsort(times[rows], kind="stable")]
                for pol in range(gains.shape[0]):
                    for channel in range(gains.shape[1]):
                        new_gains, new_flags, filled = (
                            self.bridge_solution_series(
                                times[rows],
                                gains[pol, channel, rows],
                                flags[pol, channel, rows],
                                settings["max_gap_s"],
                                settings["min_valid_fraction"],
                            )
                        )
                        if not filled:
                            continue
                        gains[pol, channel, rows] = new_gains
                        flags[pol, channel, rows] = new_flags
                        for item in filled:
                            fills.append({
                                "field_id": int(field),
                                "spectral_window": int(window),
                                "antenna_id": int(antenna),
                                "polarization": int(pol),
                                "channel": int(channel),
                                "row": int(rows[item["index"]]),
                                "time": item["time"],
                                "kind": item["kind"],
                                "span_minutes": item["span_s"] / 60.0,
                            })

            if fills:
                tool.putcol("CPARAM", gains)
                tool.putcol("FLAG", flags)
        except Exception as error:
            if opened:
                tool.close()
                opened = False
            if destination.exists():
                shutil.rmtree(destination, ignore_errors=True)
            print(
                f"WARNING: could not bridge failed solutions in {source.name} "
                f"({type(error).__name__}: {error}). Targets receive the "
                "table as solved."
            )
            return None
        finally:
            if opened:
                tool.close()

        return {
            "source": str(source),
            "table": str(destination),
            "solutions": int(flags.size),
            "flagged_before": int(np.count_nonzero(flags)) + len(fills),
            "bridged": len(fills),
            "antennas": sorted({item["antenna_id"] for item in fills}),
            "fills": fills,
        }

    def target_transfer_tables(self, phase_table, flux_table):
        """The phase and amplitude tables science targets are corrected with.

        With bridging on these are the bridged copies; otherwise, or where a
        copy could not be made, the tables as solved.
        """

        tables = {"phase": str(phase_table), "flux": str(flux_table)}
        settings = self.target_transfer_bridge_settings()
        report_path = self.caldir / "target_transfer_bridge.json"
        if not settings["enabled"]:
            return tables

        report = {
            "schema_version": 1,
            "max_gap_minutes": settings["max_gap_s"] / 60.0,
            "min_valid_fraction": settings["min_valid_fraction"],
            "tables": {},
        }
        for role, source in (("phase", phase_table), ("flux", flux_table)):
            source = Path(source)
            destination = source.with_name(f"{source.stem}.target.cal")
            record = self.bridge_transfer_table(source, destination)
            if record is None:
                report["tables"][role] = {
                    "source": str(source),
                    "status": "not_bridged",
                }
                continue
            tables[role] = record["table"]
            report["tables"][role] = dict(record, status="bridged")
            print(
                f"  {source.name}: {record['bridged']} of "
                f"{record['flagged_before']} failed solutions bridged for "
                "the targets"
                + (
                    f" (antenna ids {record['antennas']})"
                    if record["antennas"] else ""
                )
                + "."
            )

        self.write_atomic_json(report_path, report)
        return tables

    def run_applycal(self):

        if self.checkpoint_exists("applycal"):

            opened = False

            try:

                tb.open(self.ms)
                opened = True
                corrected_exists = "CORRECTED_DATA" in tb.colnames()

            finally:

                if opened:

                    tb.close()

            if corrected_exists:

                print(
                    "Applycal already completed."
                )

                return

            print(
                "WARNING: applycal checkpoint is stale because "
                "CORRECTED_DATA is missing; applying solutions again."
            )
            (self.caldir / "applycal.done").unlink()

        phasecals = self.get_phase_calibrators()
        fluxcals = self.get_flux_calibrators()

        bandpasscal = (
            self.get_bandpass_calibrator()
        )

        targets = (
            self.get_targets()
        )
        target_phase_map = self.get_target_phase_map()

        applications = []
        seen_fields = set()

        short_phase_table = str(self.caldir / "phase.cal")
        transfer_phase_table = str(self.caldir / "phase_scan.cal")
        flux_table = str(self.caldir / "flux.cal")

        def add_application(field, gain_field, phase_table,
                            amplitude_table=flux_table):

            if field and field not in seen_fields:

                applications.append(
                    (field, gain_field, phase_table, amplitude_table)
                )
                seen_fields.add(field)

        # A calibrator is corrected with its own short-interval solutions,
        # which is what its own amplitude solve assumed.  A science target
        # instead receives the scan-averaged phase table of its mapped phase
        # calibrator: interpolating the short table across the calibrator
        # cycle would transfer single-interval solution noise to the target.
        # Relying on CASA's implicit field lookup can leave a target
        # uncalibrated, so every mapping is stated explicitly.
        for fluxcal in fluxcals:
            add_application(fluxcal, fluxcal, short_phase_table)
        add_application(bandpasscal, bandpasscal, short_phase_table)
        for phasecal in phasecals:
            add_application(phasecal, phasecal, short_phase_table)

        # Calibrators above keep the tables as solved. Targets receive
        # copies in which an isolated failed solution is bridged from the
        # same antenna's good solutions on either side.
        target_tables = self.target_transfer_tables(
            transfer_phase_table,
            flux_table,
        )

        for target in targets:

            add_application(
                target,
                target_phase_map[target],
                target_tables["phase"],
                target_tables["flux"],
            )

        print(
            "Applying calibration to calibrators and targets using the "
            "configured temporal-gain mapping."
        )

        for field, gain_field, phase_table, amplitude_table in applications:

            print(
                f"  applycal → {field} (temporal gains: {gain_field}, "
                f"phase table: {Path(phase_table).name}, amplitude table: "
                f"{Path(amplitude_table).name})"
            )

            applycal(
                vis=self.ms,

                field=field,

                # prebp.cal is deliberately absent.  It is a per-integration
                # phase table of the bandpass calibrator whose only purpose
                # is to keep that field coherent while bandpass.cal is being
                # solved.  Interpolating it onto another field or onto the
                # target injects an unrelated antenna-phase offset.
                gaintable=[

                    str(
                        self.caldir /
                        "delay.cal"
                    ),

                    str(
                        self.caldir /
                        "bandpass.cal"
                    ),

                    phase_table,

                    amplitude_table,
                ],

                gainfield=[

                    bandpasscal,

                    bandpasscal,

                    gain_field,

                    gain_field,
                ],

                interp=[

                    "nearest",

                    "linear,linearflag",

                    "linear,linearflag",

                    "linear,linearflag"
                ],

                applymode=self.task_config[
                    "applymode"
                ],

                parang=self.task_config[
                    "parang"
                ],

                calwt=self.task_config[
                    "calwt"
                ]
            )

        self.mark_checkpoint(
            "applycal"
        )

    # =========================================================
    # GUARDED CORRECTED-DATA FLAGGING AND FRESH SOLUTIONS
    # =========================================================

    @staticmethod
    def summary_flag_fraction(summary):

        try:
            total = float(summary.get("total", 0.0))
            flagged = float(summary.get("flagged", 0.0))
        except (AttributeError, TypeError, ValueError):
            return None

        if (
            total <= 0 or
            flagged < 0 or
            flagged > total or
            not math.isfinite(total) or
            not math.isfinite(flagged)
        ):
            return None

        fraction = flagged / total
        return fraction if math.isfinite(fraction) else None

    @classmethod
    def summary_group_flag_fractions(cls, summary, category):

        if not isinstance(summary, dict):
            return {}

        groups = summary.get(category, {})
        if not isinstance(groups, dict):
            return {}

        fractions = {}
        for name, values in groups.items():
            fraction = cls.summary_flag_fraction(values)
            if fraction is not None:
                fractions[str(name)] = fraction

        return fractions

    def absolute_parallel_hand_correlation(self):
        """Return a strict absolute-amplitude selector for parallel hands.

        CASA correlation-type codes follow the MeasurementSet definition:
        RR/LL are 5/8 and XX/YY are 9/12.  Corrected-data auto-flagging must
        never fall back to cross-hands or ``ABS_ALL`` when the metadata are
        unreadable, because cross-hand statistics can be very different from
        those of the Stokes-I parallel hands.
        """

        cached = getattr(self, "_postcal_parallel_correlation", None)
        if cached:
            return cached

        tool = table()
        opened = False
        try:
            tool.open(str(Path(self.ms) / "POLARIZATION"))
            opened = True
            values = tool.getcol("CORR_TYPE")
            flattened = values.ravel() if hasattr(values, "ravel") else values
            present = {int(value) for value in flattened}
        except Exception as exc:
            raise RuntimeError(
                "Cannot identify parallel-hand correlations for guarded "
                f"post-calibration flagging: {exc}"
            ) from exc
        finally:
            if opened:
                tool.close()

        circular = [
            name for code, name in ((5, "RR"), (8, "LL"))
            if code in present
        ]
        linear = [
            name for code, name in ((9, "XX"), (12, "YY"))
            if code in present
        ]
        if circular and linear:
            raise RuntimeError(
                "The POLARIZATION table contains both circular- and "
                "linear-feed parallel hands; split the polarization setups "
                "before post-calibration auto-flagging."
            )

        hands = circular or linear
        if not hands:
            raise RuntimeError(
                "The POLARIZATION table has no recognized parallel-hand "
                "correlations (RR/LL or XX/YY)."
            )
        if len(hands) != 2:
            raise RuntimeError(
                "Guarded Stokes-I post-calibration flagging requires both "
                f"parallel hands, but only {','.join(hands)} is present."
            )

        # ABS applies to every comma-separated hand following the operator.
        correlation = f"ABS_{','.join(hands)}"
        self._postcal_parallel_correlation = correlation
        return correlation

    def postcal_flagdata_chunk_plan(self):
        """Bound corrected-data flagdata chunks exactly as flagging.py does.

        The corrected-data pass loads the same whole-scan cubes as the
        pre-calibration pass. Without the same bound, a dataset that only just
        survived flagging would be killed here instead, many hours later.
        """

        cached = getattr(self, "_postcal_flagdata_chunk_plan", None)
        if cached is not None:
            return cached

        limit = FLAGGING_SELECTORS.detect_memory_limit_bytes(
            configured_gb=RUN_CONFIG.get("flagdata_memory_limit_gb"),
        )

        try:
            from casatools import msmetadata

            groups = FLAGGING_SELECTORS.summarize_scan_chunks(
                self.ms,
                msmetadata,
                table,
            )
        except Exception as exc:
            plan = {
                "ntime": "scan",
                "status": "estimation_failed",
                "reason": f"scan geometry could not be read: {exc}",
            }
        else:
            plan = FLAGGING_SELECTORS.plan_flagdata_time_chunk(
                groups,
                limit["limit_bytes"],
                memory_fraction=float(RUN_CONFIG.get(
                    "flagdata_memory_fraction",
                    FLAGGING_SELECTORS.FLAGDATA_MEMORY_FRACTION,
                )),
                bytes_per_element=float(RUN_CONFIG.get(
                    "flagdata_bytes_per_element",
                    FLAGGING_SELECTORS.FLAGDATA_BYTES_PER_ELEMENT,
                )),
            )

        plan["memory_limit"] = limit
        self._postcal_flagdata_chunk_plan = plan

        caldir = getattr(self, "caldir", None)
        if caldir is not None:
            self.write_atomic_json(
                Path(caldir) / "postcal_flagdata_chunk_plan.json",
                plan,
            )

        gib = 1024 ** 3
        if plan["status"] == "insufficient_memory":
            raise RuntimeError(
                "Corrected-data flagdata cannot fit even "
                f"{plan['minimum_integrations']} integrations of the largest "
                "scan in memory: about "
                f"{plan['required_memory_bytes'] / gib:.0f} GiB is needed, "
                f"against {limit['limit_bytes'] / gib:.0f} GiB "
                f"({limit['limiting_source']}). Request more memory. No "
                "flags were written."
            )
        if plan["status"] == "bounded":
            print(
                "Corrected-data flagdata uses "
                f"ntime={plan['ntime']} because whole-scan chunks would need "
                "about "
                f"{plan['largest_scan_chunk']['scan_chunk_bytes'] / gib:.0f} "
                f"GiB against a {plan['budget_bytes'] / gib:.0f} GiB budget."
            )
        elif plan["status"] == "estimation_failed":
            print(
                "WARNING: corrected-data flagdata memory could not be "
                f"estimated ({plan['reason']}); retaining whole-scan chunks."
            )

        return plan

    def postcal_checkpoint_names(self):

        telescope = str(self.telescope).lower()
        if telescope not in {"ugmrt", "meerkat"}:
            raise RuntimeError(
                "Corrected-data post-calibration flagging is implemented "
                f"only for uGMRT and MeerKAT, not {self.telescope}."
            )

        return {
            "flagging": f"{telescope}_postcal_flagging",
            "prepared": f"{telescope}_round2_prepared",
            "complete": f"{telescope}_round2_complete",
        }

    def write_atomic_json(self, path, value):

        path = Path(path)
        temporary = path.with_name(
            f".{path.name}.tmp.{os.getpid()}"
        )
        with open(temporary, "w") as handle:
            json.dump(value, handle, indent=2)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)

    def run_postcal_flagging(self):

        if not self.task_config.get("postcal_flagging", False):
            return

        names = self.postcal_checkpoint_names()
        checkpoint = names["flagging"]
        telescope_label = (
            "uGMRT" if self.telescope == "ugmrt" else "MeerKAT"
        )
        if self.checkpoint_exists(checkpoint):
            return

        state_path = self.caldir / f"{checkpoint}_state.json"
        previous = {}
        if state_path.is_file():
            with open(state_path) as handle:
                previous = json.load(handle)

        previous_version = previous.get("flag_version")
        if previous_version:
            print(
                "Restoring the saved flag state from an interrupted "
                f"or incompletely checkpointed {telescope_label} "
                "corrected-data "
                "flagging pass."
            )
            flagmanager(
                vis=self.ms,
                mode="restore",
                versionname=previous_version,
            )

        fields = list(dict.fromkeys(
            field
            for field in (
                self.get_flux_calibrators() +
                [self.get_bandpass_calibrator()] +
                self.get_phase_calibrators() +
                self.get_targets()
            )
            if field
        ))
        if not fields:
            raise RuntimeError(
                "No calibrator or target fields are available for the "
                f"{telescope_label} corrected-data flagging pass."
            )
        field_selection = ",".join(fields)

        before = flagdata(
            vis=self.ms,
            mode="summary",
            field=field_selection,
            action="calculate",
            ntime=self.postcal_flagdata_chunk_plan()["ntime"],
        )
        before_fraction = self.summary_flag_fraction(before)
        if before_fraction is None:
            raise RuntimeError(
                f"Could not measure flags before {telescope_label} "
                "post-calibration flagging."
            )

        correlation = self.absolute_parallel_hand_correlation()

        version = (
            f"before_{self.telescope}_postcal_{os.getpid()}_{time.time_ns()}"
        )
        flagmanager(
            vis=self.ms,
            mode="save",
            versionname=version,
            comment="Before corrected-data tfcrop/rflag and round-2 solves",
        )

        state = {
            "status": "running",
            "telescope": self.telescope,
            "band": self.band,
            "flag_version": version,
            "before_fraction": before_fraction,
            "correlation": correlation,
        }
        self.write_atomic_json(state_path, state)

        try:
            flagdata(
                vis=self.ms,
                mode="tfcrop",
                field=field_selection,
                datacolumn="CORRECTED",
                correlation=correlation,
                ntime=self.postcal_flagdata_chunk_plan()["ntime"],
                combinescans=False,
                timecutoff=float(self.task_config.get(
                    "postcal_tfcrop_timecutoff", 4.0
                )),
                freqcutoff=float(self.task_config.get(
                    "postcal_tfcrop_freqcutoff", 4.0
                )),
                timefit="line",
                freqfit="poly",
                flagdimension="freqtime",
                usewindowstats="none",
                extendflags=False,
                action="apply",
                flagbackup=False,
            )
            flagdata(
                vis=self.ms,
                mode="rflag",
                field=field_selection,
                datacolumn="CORRECTED",
                correlation=correlation,
                ntime=self.postcal_flagdata_chunk_plan()["ntime"],
                combinescans=False,
                timedevscale=float(self.task_config.get(
                    "postcal_rflag_timedevscale", 5.0
                )),
                freqdevscale=float(self.task_config.get(
                    "postcal_rflag_freqdevscale", 5.0
                )),
                winsize=int(self.task_config.get(
                    "postcal_rflag_winsize", 5
                )),
                extendflags=False,
                action="apply",
                flagbackup=False,
            )

            if bool(self.task_config.get("postcal_extend", True)):
                # A separate pass, so the detectors above stay scan-local and
                # non-extending: they decide what is interference, and this
                # removes the wings around what they found. combinescans
                # keeps the growth inside one scan.
                flagdata(
                    vis=self.ms,
                    mode="extend",
                    field=field_selection,
                    datacolumn="CORRECTED",
                    correlation=correlation,
                    ntime=self.postcal_flagdata_chunk_plan()["ntime"],
                    combinescans=False,
                    growtime=float(self.task_config.get(
                        "postcal_extend_growtime", 80.0
                    )),
                    growfreq=float(self.task_config.get(
                        "postcal_extend_growfreq", 80.0
                    )),
                    growaround=False,
                    flagneartime=False,
                    flagnearfreq=False,
                    extendpols=True,
                    action="apply",
                    flagbackup=False,
                )

            after = flagdata(
                vis=self.ms,
                mode="summary",
                field=field_selection,
                action="calculate",
                ntime=self.postcal_flagdata_chunk_plan()["ntime"],
            )
            after_fraction = self.summary_flag_fraction(after)
            if after_fraction is None:
                raise RuntimeError(
                    "Corrected-data flagging returned no valid summary."
                )

            new_fraction = max(0.0, after_fraction - before_fraction)
            max_new = float(self.task_config.get(
                "postcal_max_new_fraction", 0.20
            ))
            max_total = float(self.task_config.get(
                "postcal_max_total_fraction", 0.75
            ))
            max_group_new = float(self.task_config.get(
                "postcal_max_group_new_fraction", 0.30
            ))
            max_group_total = float(self.task_config.get(
                "postcal_max_group_total_fraction", 0.90
            ))
            violations = []
            if new_fraction > max_new:
                violations.append(
                    f"new fraction {new_fraction:.4f} exceeds {max_new:.4f}"
                )
            if after_fraction > max_total:
                violations.append(
                    "total fraction "
                    f"{after_fraction:.4f} exceeds {max_total:.4f}"
                )

            for category in ("field", "spw"):
                before_groups = self.summary_group_flag_fractions(
                    before,
                    category,
                )
                after_groups = self.summary_group_flag_fractions(
                    after,
                    category,
                )
                for group, group_after in after_groups.items():
                    group_before = before_groups.get(group, 0.0)
                    group_new = max(0.0, group_after - group_before)
                    if group_new > max_group_new:
                        violations.append(
                            f"{category} {group} added {group_new:.4f}, "
                            f"exceeding {max_group_new:.4f}"
                        )
                    if group_after > max_group_total:
                        violations.append(
                            f"{category} {group} total {group_after:.4f} "
                            f"exceeds {max_group_total:.4f}"
                        )

            if violations:
                raise RuntimeError(
                    f"{telescope_label} post-calibration flagging exceeded "
                    "its safety limits: " + "; ".join(violations)
                )

        except Exception as exc:
            flagmanager(
                vis=self.ms,
                mode="restore",
                versionname=version,
            )
            state["status"] = "failed_restored"
            state["failure_reason"] = str(exc)
            self.write_atomic_json(state_path, state)
            raise

        # Keep the per-field, per-scan and per-SPW breakdowns. The limits
        # above are field- and SPW-level, so a single calibrator scan can be
        # flagged entirely while the pass is still accepted; without these
        # breakdowns that loss cannot be attributed to this pass afterwards.
        breakdowns = {
            stage_name: {
                category: self.summary_group_flag_fractions(summary, category)
                for category in ("field", "scan", "spw")
            }
            for stage_name, summary in (("before", before), ("after", after))
        }

        state.update({
            "status": "completed",
            "fields": fields,
            "after_fraction": after_fraction,
            "new_fraction": new_fraction,
            "group_flag_fractions": breakdowns,
            "parameters": {
                "datacolumn": "CORRECTED",
                "correlation": correlation,
                "ntime": self.postcal_flagdata_chunk_plan()["ntime"],
                "flagdata_chunk_plan": self.postcal_flagdata_chunk_plan(),
                "combinescans": False,
                "extendflags": False,
                "extend": {
                    "applied": bool(
                        self.task_config.get("postcal_extend", True)
                    ),
                    "combinescans": False,
                    "growtime": float(self.task_config.get(
                        "postcal_extend_growtime", 80.0
                    )),
                    "growfreq": float(self.task_config.get(
                        "postcal_extend_growfreq", 80.0
                    )),
                },
                **{
                key: self.task_config.get(key)
                for key in [
                    "postcal_tfcrop_timecutoff",
                    "postcal_tfcrop_freqcutoff",
                    "postcal_rflag_timedevscale",
                    "postcal_rflag_freqdevscale",
                    "postcal_rflag_winsize",
                    "postcal_max_new_fraction",
                    "postcal_max_total_fraction",
                    "postcal_max_group_new_fraction",
                    "postcal_max_group_total_fraction",
                ]
                },
            },
        })
        self.write_atomic_json(state_path, state)
        self.write_atomic_json(
            self.caldir / f"{checkpoint}_report.json",
            state,
        )
        self.mark_checkpoint(checkpoint)

    # Backward-compatible entry point retained for external wrappers.
    def run_meerkat_postcal_flagging(self):

        return self.run_postcal_flagging()

    # ---------------------------------------------------------
    # Flags written by the provisional (round-1) solutions
    # ---------------------------------------------------------

    def round1_flag_state_path(self):

        return self.caldir / "round1_flag_state.json"

    def flag_version_table(self, version):
        """Path of the table CASA keeps a named flag version in."""

        return Path(str(self.ms) + ".flagversions") / f"flags.{version}"

    def save_flags_before_round1_applycal(self):
        """Remember the flags as they are before provisional solutions touch them.

        With corrected-data flagging enabled, calibration is solved twice:
        the first set of solutions exists only to produce CORRECTED_DATA
        that the interference detectors can read. applycal nevertheless
        writes a flag for every visibility whose provisional solution
        failed, and a solve on data that still holds interference fails
        often. Those flags have to be taken back before the second solve
        (``release_provisional_solution_flags``), which needs to know what
        the flags were before the first applycal.
        """

        if not self.task_config.get("postcal_flagging", False):
            return None
        if not bool(self.task_config.get("release_provisional_flags", True)):
            return None

        names = self.postcal_checkpoint_names()
        if (
            self.checkpoint_exists("applycal") or
            self.checkpoint_exists(names["flagging"]) or
            self.checkpoint_exists(names["prepared"]) or
            self.checkpoint_exists(names["complete"])
        ):
            # The first applycal has already run. Whatever it flagged can
            # only be taken back if this state was saved before it.
            return None

        state_path = self.round1_flag_state_path()
        if state_path.is_file():
            try:
                with open(state_path) as handle:
                    state = json.load(handle)
            except (OSError, ValueError):
                state = {}
            version = state.get("before_applycal_version")
            if version and self.flag_version_table(version).is_dir():
                # Saved by an earlier, interrupted attempt at this applycal.
                # It predates every flag that attempt wrote, so it is kept.
                return version

        version = (
            f"before_{self.telescope}_round1_applycal_"
            f"{os.getpid()}_{time.time_ns()}"
        )
        flagmanager(
            vis=self.ms,
            mode="save",
            versionname=version,
            comment="Before the provisional (round-1) applycal",
        )
        self.write_atomic_json(
            state_path,
            {"status": "saved", "before_applycal_version": version},
        )
        print(
            "Saved the flag state before the provisional applycal as "
            f"'{version}'."
        )
        return version

    def release_provisional_solution_flags(self):
        """Take back the flags that only the provisional solutions wrote.

        Three flag states are known at this point:

            A   before the round-1 applycal
            B   after it (A plus every visibility whose provisional
                solution failed), saved before the detectors ran
            C   now (B plus what tfcrop/rflag/extend found)

        The second round should start from the interference the detectors
        found and nothing else, which is ``A | (C & ~B)``. Keeping B, as
        this stage used to, made every failed provisional solution a
        permanent loss: in the run that prompted this, the first applycal
        took the flagged fraction from 32.7% to 49.2% before any final
        solution had been attempted. About half of that was two antennas
        with no solutions at all, which round 2 flags again; the rest was
        visibilities whose provisional solution failed, and the second
        solve then failed more often for want of them.

        Visibilities released here are flagged again by the round-2
        applycal if their final solution fails too. They were hidden from
        the corrected-data detectors, so they carry only the raw-data
        flagging; the residual pass in self-calibration still sees them.

        The step cannot make things worse than they were: the current
        flags are saved first and restored in full if anything fails.
        """

        if not self.task_config.get("postcal_flagging", False):
            return None
        if not bool(self.task_config.get("release_provisional_flags", True)):
            return None

        names = self.postcal_checkpoint_names()
        checkpoint = f"{str(self.telescope).lower()}_round1_flags_released"
        report_path = self.caldir / "round1_flag_release_report.json"

        if self.checkpoint_exists(checkpoint):
            return None
        if (
            self.checkpoint_exists(names["prepared"]) or
            self.checkpoint_exists(names["complete"])
        ):
            # Round 2 was already begun on the old flags; changing them now
            # would leave its tables and its data inconsistent.
            return None

        def skip(reason):
            print(
                "WARNING: the flags written by the provisional solutions "
                f"are kept, because {reason}."
            )
            self.write_atomic_json(
                report_path,
                {"status": "skipped", "reason": reason},
            )
            return None

        try:
            with open(self.round1_flag_state_path()) as handle:
                before_applycal = json.load(handle).get(
                    "before_applycal_version"
                )
        except (OSError, ValueError):
            before_applycal = None
        try:
            with open(self.caldir / f"{names['flagging']}_state.json") as handle:
                after_applycal = json.load(handle).get("flag_version")
        except (OSError, ValueError):
            after_applycal = None

        if not before_applycal:
            return skip(
                "no flag state was saved before the round-1 applycal (it "
                "ran under an earlier release)"
            )
        if not after_applycal:
            return skip("the post-applycal flag state is not recorded")

        table_a = self.flag_version_table(before_applycal)
        table_b = self.flag_version_table(after_applycal)
        for label, path in (("before", table_a), ("after", table_b)):
            if not path.is_dir():
                return skip(
                    f"the flag version saved {label} the round-1 applycal "
                    f"is missing ({path.name})"
                )

        safety_version = (
            f"before_{self.telescope}_round1_release_"
            f"{os.getpid()}_{time.time_ns()}"
        )
        flagmanager(
            vis=self.ms,
            mode="save",
            versionname=safety_version,
            comment="Before releasing provisional-solution flags",
        )

        print(
            "Releasing the flags written only by the provisional (round-1) "
            "solutions before the round-2 solve."
        )

        try:
            counts = self.combine_released_flags(table_a, table_b)
        except Exception as exc:
            flagmanager(
                vis=self.ms,
                mode="restore",
                versionname=safety_version,
            )
            print(
                "WARNING: the provisional-solution flags could not be "
                f"released ({type(exc).__name__}: {exc}). The flags were "
                "restored to what they were, and round 2 proceeds on them "
                "as it did before this step existed."
            )
            self.write_atomic_json(
                report_path,
                {
                    "status": "failed_restored",
                    "reason": f"{type(exc).__name__}: {exc}",
                    "restored_version": safety_version,
                },
            )
            return None

        total = max(1, counts["total"])
        report = {
            "status": "completed",
            "before_applycal_version": before_applycal,
            "after_applycal_version": after_applycal,
            "safety_version": safety_version,
            "total_flags": counts["total"],
            "flagged_before_release": counts["before"],
            "flagged_after_release": counts["after"],
            "released": counts["before"] - counts["after"],
            "fraction_before_release": counts["before"] / total,
            "fraction_after_release": counts["after"] / total,
            "fraction_released": (counts["before"] - counts["after"]) / total,
            "rows": counts["rows"],
            "rows_handled_singly": counts["rows_handled_singly"],
        }
        self.write_atomic_json(report_path, report)
        self.mark_checkpoint(checkpoint)
        print(
            f"Released {report['fraction_released']:.2%} of the data: "
            f"{report['fraction_before_release']:.2%} was flagged, "
            f"{report['fraction_after_release']:.2%} is now. The round-2 "
            "applycal flags again whatever its own solutions cannot "
            "calibrate."
        )
        return report

    def combine_released_flags(self, table_a, table_b, chunk_values=64000000):
        """Write ``A | (C & ~B)`` into the working MS, a block of rows at a time.

        ``table_a`` and ``table_b`` are CASA flag-version tables, which hold
        the FLAG (and FLAG_ROW) of every MAIN row. The operation is
        idempotent, so a run interrupted partway through repeats it safely.
        """

        import numpy as np

        main = table()
        saved_a = table()
        saved_b = table()
        opened = []
        counts = {
            "total": 0, "before": 0, "after": 0,
            "rows": 0, "rows_handled_singly": 0,
        }

        def released(a, b, c):
            a = np.asarray(a, dtype=bool)
            b = np.asarray(b, dtype=bool)
            c = np.asarray(c, dtype=bool)
            if a.shape != c.shape or b.shape != c.shape:
                raise RuntimeError(
                    "flag versions and the Measurement Set disagree in "
                    f"shape: {a.shape}, {b.shape}, {c.shape}"
                )
            return c, a | (c & ~b)

        def tally(current, new):
            counts["total"] += int(current.size)
            counts["before"] += int(np.count_nonzero(current))
            counts["after"] += int(np.count_nonzero(new))

        try:
            main.open(str(self.ms), nomodify=False)
            opened.append(main)
            saved_a.open(str(table_a))
            opened.append(saved_a)
            saved_b.open(str(table_b))
            opened.append(saved_b)

            nrows = int(main.nrows())
            if int(saved_a.nrows()) != nrows or int(saved_b.nrows()) != nrows:
                raise RuntimeError(
                    "flag versions do not have the Measurement Set's row "
                    f"count ({saved_a.nrows()}, {saved_b.nrows()}, {nrows})"
                )
            if nrows < 1:
                raise RuntimeError("the Measurement Set has no rows")
            for tool, label in ((saved_a, table_a), (saved_b, table_b)):
                if "FLAG" not in set(tool.colnames()):
                    raise RuntimeError(
                        f"flag version {Path(label).name} has no FLAG column"
                    )
            row_flags = all(
                "FLAG_ROW" in set(tool.colnames())
                for tool in (main, saved_a, saved_b)
            )

            probe = np.asarray(main.getcell("FLAG", 0))
            per_row = max(1, int(probe.size))
            block = max(1, int(chunk_values) // per_row)

            for start in range(0, nrows, block):
                count = min(block, nrows - start)

                if row_flags:
                    row_c, row_new = released(
                        saved_a.getcol("FLAG_ROW", startrow=start, nrow=count),
                        saved_b.getcol("FLAG_ROW", startrow=start, nrow=count),
                        main.getcol("FLAG_ROW", startrow=start, nrow=count),
                    )
                    if np.any(row_new != row_c):
                        main.putcol(
                            "FLAG_ROW", row_new, startrow=start, nrow=count
                        )

                try:
                    arrays = (
                        saved_a.getcol("FLAG", startrow=start, nrow=count),
                        saved_b.getcol("FLAG", startrow=start, nrow=count),
                        main.getcol("FLAG", startrow=start, nrow=count),
                    )
                except Exception:
                    # Rows of different shapes (spectral windows with
                    # different channel counts) cannot be read as one
                    # array, so this block is done a row at a time.
                    arrays = None

                if arrays is None:
                    for row in range(start, start + count):
                        current, new = released(
                            saved_a.getcell("FLAG", row),
                            saved_b.getcell("FLAG", row),
                            main.getcell("FLAG", row),
                        )
                        if np.any(new != current):
                            main.putcell("FLAG", row, new)
                        tally(current, new)
                    counts["rows_handled_singly"] += count
                else:
                    current, new = released(*arrays)
                    if np.any(new != current):
                        main.putcol("FLAG", new, startrow=start, nrow=count)
                    tally(current, new)
                counts["rows"] += count

            flush = getattr(main, "flush", None)
            if callable(flush):
                flush()

        finally:
            for tool in reversed(opened):
                try:
                    tool.close()
                except Exception:
                    pass

        return counts

    def prepare_second_round(self):

        names = self.postcal_checkpoint_names()
        checkpoint = names["prepared"]
        telescope_label = (
            "uGMRT" if self.telescope == "ugmrt" else "MeerKAT"
        )
        if self.checkpoint_exists(checkpoint):
            return

        round1 = self.caldir / "round1"
        round1.mkdir(exist_ok=True)

        filenames = [
            "delay.cal",
            "prebp.cal",
            "bandpass.cal",
            "phase.cal",
            "phase_scan.cal",
            "amp.cal",
            "flux.cal",
        ]

        # Validate the complete archive before moving another table. A source
        # absent here but present in round1 means a prior move was interrupted
        # and is already in the intended durable location.
        pending_moves = []
        for filename in filenames:
            source = self.caldir / filename
            destination = round1 / filename
            source_exists = source.exists()
            destination_exists = destination.exists()
            if source_exists and destination_exists:
                raise RuntimeError(
                    f"{telescope_label} round-1 archive is ambiguous because "
                    "both the "
                    f"source and destination exist for {filename}. Review "
                    f"{source} and {destination} before resuming."
                )
            if not source_exists and not destination_exists:
                raise RuntimeError(
                    f"{telescope_label} round-1 calibration table is "
                    f"missing: {source}"
                )
            if source_exists:
                pending_moves.append((source, destination))

        for source, destination in pending_moves:
            source.rename(destination)

        for stage in [
            "delay",
            "prebp",
            "bandpass",
            "phase",
            "phase_scan",
            "amp",
            "fluxscale",
            "applycal",
        ]:
            (self.caldir / f"{stage}.done").unlink(missing_ok=True)

        self.mark_checkpoint(checkpoint)

    # Backward-compatible entry point retained for external wrappers.
    def prepare_meerkat_second_round(self):

        return self.prepare_second_round()

    def run_second_round(self):

        if not self.task_config.get("postcal_flagging", False):
            self.calibration_rounds = 1
            return

        names = self.postcal_checkpoint_names()
        if self.checkpoint_exists(names["complete"]):
            self.calibration_rounds = 2
            return

        telescope_label = (
            "uGMRT" if self.telescope == "ugmrt" else "MeerKAT"
        )
        print(
            f"Running guarded {telescope_label} corrected-data flagging and "
            "fresh "
            "round-2 calibration solutions."
        )
        self.run_postcal_flagging()
        self.release_provisional_solution_flags()
        self.prepare_second_round()
        self.run_delaycal()
        self.run_prebp()
        self.run_bandpass()
        self.run_phasecal()
        self.run_phase_transfer_cal()
        self.run_ampcal()
        self.run_fluxscale()
        self.run_applycal()
        self.write_atomic_json(
            self.caldir / f"{self.telescope}_round2_report.json",
            {
                "status": "completed",
                "telescope": self.telescope,
                "band": self.band,
                "postcal_flagging_checkpoint": names["flagging"],
                "round1_archive": str(self.caldir / "round1"),
                "round2_tables": {
                    name: str(self.caldir / name)
                    for name in [
                        "delay.cal",
                        "prebp.cal",
                        "bandpass.cal",
                        "phase.cal",
                        "phase_scan.cal",
                        "amp.cal",
                        "flux.cal",
                    ]
                },
            },
        )
        self.mark_checkpoint(names["complete"])
        self.calibration_rounds = 2

    # Backward-compatible entry point retained for external wrappers.
    def run_meerkat_second_round(self):

        return self.run_second_round()

    # =========================================================
    # SPLIT TARGETS
    # =========================================================

    def statistical_weights_settings(self):
        config = self.task_config
        return {
            "enabled": bool(config.get("statwt_enabled", False)),
            "timebin": str(config.get("statwt_timebin", "120s")),
            "chanbin": int(config.get("statwt_chanbin", 16)),
            "minsamp": int(config.get("statwt_minsamp", 16)),
        }

    @staticmethod
    def statwt_marker(outms):
        return Path(str(outms) + ".statwt.json")

    def set_statistical_weights(self, outms, target):
        """Give every channel of a calibrated target the weight its noise earns.

        After calibration a channel's noise in Jy follows the system's
        equivalent flux density there (Tsys over effective area) and any
        interference left in it, not the bandpass gain: the gain scales
        signal and noise alike and divides out. Across a uGMRT Band-4 band
        the SEFD rises towards the edges and residual interference is
        patchy, so channels do differ in noise, by an amount the bandpass
        amplitude alone cannot tell. Data from the GMRT carry one weight per
        row, shared by every channel, so imaging would give the noisiest
        channels the same say as the best.
        statwt measures the scatter of the calibrated visibilities in bins of
        ``statwt_chanbin`` channels and ``statwt_timebin`` per baseline and
        correlation, and writes 1/variance (Jy^-2) to WEIGHT_SPECTRUM.

        The scatter includes the sky, which is negligible when each
        visibility's noise (a few tenths of a Jy for uGMRT) dwarfs the field's
        flux; self-calibration checks that once it has a model. Bins with
        fewer than ``statwt_minsamp`` unflagged values are flagged rather than
        weighted from a handful of numbers. A marker is written only when the
        run completes, so an interrupted run's half-reweighted MS is re-split
        rather than reused.
        """

        settings = self.statistical_weights_settings()
        if not settings["enabled"]:
            return {"status": "disabled"}
        if not CASA_COMPAT.task_available(statwt):
            raise RuntimeError(
                "statwt is not available in this CASA installation, so the "
                "calibrated targets cannot be given per-channel weights. "
                "Set statwt_enabled: false in config.yaml to proceed without "
                "them."
            )

        print(
            f"  weights → {target}: statwt over {settings['chanbin']}-channel "
            f"x {settings['timebin']} bins, at least {settings['minsamp']} "
            "samples each"
        )
        statwt(
            vis=str(outms),
            datacolumn="data",
            timebin=settings["timebin"],
            chanbin=settings["chanbin"],
            minsamp=settings["minsamp"],
            combine="",
            flagbackup=True,
        )

        record = dict(settings)
        record["status"] = "applied"
        spread = self.channel_weight_spread(outms)
        if spread:
            record.update(spread)
            print(
                "    per-channel weights span a factor of "
                f"{spread['channel_weight_ratio_p95_p5']:.0f} (5th to 95th "
                "percentile); weighting every channel alike would have left "
                "the image noise "
                f"{spread['equal_weight_noise_penalty']:.1f} times higher."
            )
        record["completed_utc"] = datetime.now(timezone.utc).isoformat()
        self.write_atomic_json(self.statwt_marker(outms), record)
        return record

    @staticmethod
    def channel_weight_spread(outms, rows_wanted=4000, block_rows=250):
        """How unequal the per-channel weights are, from a bounded sample.

        For per-channel noise sigma_c = w_c**-0.5, equal weighting gives an
        image noise sqrt(mean(sigma_c^2) * mean(1/sigma_c^2)) times the
        inverse-variance weighted one: the penalty weights avoid.
        """

        import numpy as np

        tool = table()
        opened = False
        try:
            tool.open(str(outms))
            opened = True
            if "WEIGHT_SPECTRUM" not in tool.colnames():
                return None
            nrows = int(tool.nrows())
            if nrows < 1:
                return None
            starts = list(range(0, nrows, block_rows))
            wanted = max(1, rows_wanted // block_rows)
            if len(starts) > wanted:
                step = len(starts) / float(wanted)
                starts = sorted({starts[int(i * step)] for i in range(wanted)})
            sums = None
            counts = None
            for start in starts:
                count = min(block_rows, nrows - start)
                weights = np.asarray(tool.getcol(
                    "WEIGHT_SPECTRUM", startrow=start, nrow=count
                ), dtype=float)
                flags = np.asarray(tool.getcol(
                    "FLAG", startrow=start, nrow=count
                ), dtype=bool)
                good = (~flags) & np.isfinite(weights) & (weights > 0)
                # (correlation, channel, row) -> per-channel totals
                block_sums = np.where(good, weights, 0.0).sum(axis=(0, 2))
                block_counts = good.sum(axis=(0, 2))
                sums = block_sums if sums is None else sums + block_sums
                counts = block_counts if counts is None else counts + block_counts
        except Exception:
            return None
        finally:
            if opened:
                tool.close()

        if sums is None:
            return None
        valid = counts > 0
        if np.count_nonzero(valid) < 2:
            return None
        channel_weights = sums[valid] / counts[valid]
        low, high = np.percentile(channel_weights, [5.0, 95.0])
        penalty = math.sqrt(
            float(np.mean(1.0 / channel_weights)) *
            float(np.mean(channel_weights))
        )
        return {
            "channels_with_weights": int(np.count_nonzero(valid)),
            "channel_weight_ratio_p95_p5": (
                float(high / low) if low > 0 else None
            ),
            "equal_weight_noise_penalty": penalty,
        }

    # A calibrated target is written at the channel width imaging needs, not
    # the correlator's. uGMRT Band 4 records 4096 channels of 97.7 kHz, and
    # every major cycle of every self-calibration and final image reads and
    # grids all of them: about seven minutes a cycle and five hours an image
    # for a 33-minute target. Averaging four channels (390.6 kHz) reduces the
    # peak of a source at the edge of the 5120-pixel, 1-arcsec field by 2%
    # through bandwidth smearing, and by under 1% within 30 arcmin, while a
    # quarter of the data is read and gridded. The data are Hanning-smoothed
    # first, so neighbouring channels were not independent to begin with.

    def target_channel_average_path(self):
        return self.caldir / "target_channel_average.json"

    @staticmethod
    def spectral_window_channels(vis):
        """(NUM_CHAN, widest |CHAN_WIDTH| in Hz) for every spectral window."""

        tool = table()
        opened = False
        windows = []
        try:
            tool.open(str(Path(vis) / "SPECTRAL_WINDOW"))
            opened = True
            for row in range(int(tool.nrows())):
                widths = [
                    abs(float(value))
                    for value in list(tool.getcell("CHAN_WIDTH", row))
                ]
                windows.append((
                    int(tool.getcell("NUM_CHAN", row)),
                    max(widths) if widths else 0.0,
                ))
        finally:
            if opened:
                tool.close()
        return windows

    def target_channel_average_plan(self):
        """How many native channels go into each channel of a target MS.

        The band's ``target_max_channel_width_khz`` is the widest channel
        imaging may be given. The factor is the largest that keeps every
        output channel within it and divides every window's channel count,
        so no window ends in a partial channel. A band without the setting,
        or a native channel already wider than half of it, is not averaged.
        """

        band_config = getattr(self, "band_config", None) or {}
        limit = band_config.get("target_max_channel_width_khz")
        plan = {
            "schema_version": 1,
            "max_channel_width_khz": limit,
            "width": 1,
        }
        try:
            limit_hz = float(limit) * 1.0e3
        except (TypeError, ValueError):
            return plan
        if not math.isfinite(limit_hz) or limit_hz <= 0:
            return plan

        try:
            windows = self.spectral_window_channels(self.ms)
        except Exception as error:
            # Averaging only saves time. Without the channel widths the
            # targets are written at full resolution, which is never worse.
            print(
                "WARNING: could not read the channel widths of "
                f"{self.ms} ({error}); calibrated targets keep their native "
                "channels."
            )
            plan["status"] = "channel_widths_unreadable"
            return plan

        native_hz = max((width for _, width in windows), default=0.0)
        counts = [count for count, _ in windows]
        if native_hz <= 0 or not counts or min(counts) < 1:
            plan["status"] = "channel_widths_unreadable"
            return plan

        width = int(math.floor(limit_hz / native_hz + 1.0e-9))
        while width > 1 and any(count % width for count in counts):
            width -= 1
        width = max(width, 1)
        plan.update({
            "width": width,
            "native_channel_width_khz": native_hz / 1.0e3,
            "output_channel_width_khz": native_hz * width / 1.0e3,
            "native_channels": counts,
            "output_channels": [count // width for count in counts],
        })
        return plan

    @staticmethod
    def discard_native_target(native):
        """Remove a full-resolution target copy and what statwt left with it."""

        native = Path(native)
        if not native.name.endswith(".native.ms"):
            raise ValueError(
                f"refusing to remove {native}: not a native-resolution "
                "target copy"
            )
        for path in (native, Path(str(native) + ".flagversions")):
            if path.is_dir():
                shutil.rmtree(path)
        marker = Path(str(native) + ".statwt.json")
        if marker.is_file():
            marker.unlink()

    def split_averaged_target(self, target, outms, plan):
        """Split one target at full resolution, weight it, then average it.

        statwt runs on the native channels, as it does for an unaveraged
        target, so ``statwt_chanbin`` keeps its meaning. ``split`` then
        averages with those per-channel weights and sums them into the
        averaged channel's WEIGHT_SPECTRUM; an output channel is flagged
        only when every channel in it was. The statwt marker is written
        for the averaged MS last, so an interruption anywhere before it
        leaves a target that is split again rather than reused.
        """

        outms = Path(outms)
        native = outms.with_name(f"{outms.stem}.native.ms")
        self.discard_native_target(native)

        split(
            vis=self.ms,
            field=target,
            outputvis=str(native),
            datacolumn="corrected",
            keepflags=True,
        )
        if not (native.is_dir() and (native / "table.dat").is_file()):
            raise RuntimeError(
                "split did not produce a usable target MS: "
                f"{native}"
            )

        weights = self.set_statistical_weights(native, target)

        print(
            f"  average → {target}: {plan['width']} channels of "
            f"{plan['native_channel_width_khz']:.2f} kHz into "
            f"{plan['output_channel_width_khz']:.1f} kHz"
        )
        split(
            vis=str(native),
            outputvis=str(outms),
            datacolumn="data",
            width=int(plan["width"]),
            keepflags=True,
        )
        if not (outms.is_dir() and (outms / "table.dat").is_file()):
            raise RuntimeError(
                "split did not produce a usable averaged target MS: "
                f"{outms}"
            )

        if weights.get("status") == "applied":
            record = dict(weights)
            record["channel_average"] = {
                key: plan[key]
                for key in (
                    "width",
                    "native_channel_width_khz",
                    "output_channel_width_khz",
                )
            }
            self.write_atomic_json(self.statwt_marker(outms), record)

        self.discard_native_target(native)

    def split_targets(self):

        provenance_file = self.caldir / "split_targets.provenance.json"

        if self.checkpoint_exists("split"):

            report_file = self.caldir / "split_targets.json"

            try:

                with open(report_file) as handle:

                    existing_report = json.load(handle)

            except (OSError, ValueError, TypeError):

                existing_report = {}

            if not isinstance(existing_report, dict):
                existing_report = {}

            try:
                with open(provenance_file) as handle:
                    existing_provenance = json.load(handle)
            except (OSError, ValueError, TypeError):
                existing_provenance = {}

            recorded_targets = (
                existing_provenance.get("targets", {})
                if isinstance(existing_provenance, dict)
                else {}
            )
            reusable = bool(
                existing_report and
                isinstance(recorded_targets, dict) and
                set(recorded_targets) == set(existing_report)
            )
            if reusable:
                for target, product in existing_report.items():
                    recorded = recorded_targets.get(target, {})
                    product_path = Path(product).expanduser().resolve()
                    try:
                        current_identity = (
                            IMAGING_SAFETY.measurement_set_instance_identity(
                                product_path
                            )
                        )
                    except (OSError, ValueError):
                        reusable = False
                        break
                    if (
                        not isinstance(recorded, dict) or
                        Path(recorded.get("path", "")).expanduser().resolve() !=
                        product_path or
                        recorded.get("instance_identity") != current_identity
                    ):
                        reusable = False
                        break

            if reusable:

                print(
                    "Target splitting already completed."
                )
                if self.statistical_weights_settings()["enabled"]:
                    unweighted = [
                        target for target, product in existing_report.items()
                        if not self.statwt_marker(product).is_file()
                    ]
                    if unweighted:
                        print(
                            "WARNING: these targets were split by an earlier "
                            "release, without per-channel statistical "
                            f"weights: {', '.join(unweighted)}. They are left "
                            "unchanged so their downstream products stay "
                            "consistent; run into a new --outdir to weight "
                            "them."
                        )

                return

            print(
                "WARNING: split checkpoint is stale; rebuilding the "
                "missing target products."
            )
            (self.caldir / "split.done").unlink()
            for product in existing_report.values():
                outms = Path(product)
                if not outms.exists():
                    continue
                index = 1
                while True:
                    archived = outms.with_name(
                        f"{outms.name}.incomplete.{index}"
                    )
                    if not archived.exists():
                        break
                    index += 1
                outms.rename(archived)
                print(
                    "  preserved unbound target MS as "
                    f"{archived}"
                )

        targets = (
            self.get_targets()
        )

        print(
            "Splitting calibrated targets."
        )

        channel_average = self.target_channel_average_plan()

        split_report = {}

        for target in targets:

            outms = (
                self.splitdir /
                f"{target}.ms"
            )

            print(
                f"  split → {target}"
            )

            weights_pending = (
                self.statistical_weights_settings()["enabled"] and
                not self.statwt_marker(outms).is_file()
            )

            if (
                outms.is_dir() and
                (outms / "table.dat").is_file() and
                not weights_pending
            ):

                print(
                    f"  reusing existing target MS: {outms}"
                )

            else:

                if outms.exists():

                    index = 1

                    while True:

                        archived = outms.with_name(
                            f"{outms.name}.incomplete.{index}"
                        )

                        if not archived.exists():

                            break

                        index += 1

                    outms.rename(archived)
                    print(
                        "  preserved incomplete target MS as "
                        f"{archived}"
                    )

                if channel_average["width"] > 1:

                    self.split_averaged_target(
                        target, outms, channel_average
                    )

                else:

                    split(
                        vis=self.ms,

                        field=target,

                        outputvis=str(outms),

                        datacolumn="corrected",

                        keepflags=True
                    )

                    if not (
                        outms.is_dir() and
                        (outms / "table.dat").is_file()
                    ):

                        raise RuntimeError(
                            "split did not produce a usable target MS: "
                            f"{outms}"
                        )

                    self.set_statistical_weights(outms, target)

            split_report[
                target
            ] = str(outms)

        # Save the legacy path map plus a strict, mutation-tolerant instance
        # identity for each split target. Selfcal may add MODEL_DATA to the
        # MAIN table, but it must never replace the target MS or its schema.
        split_provenance = {
            "schema_version": 1,
            "targets": {
                target: {
                    "path": str(Path(product).expanduser().resolve()),
                    "instance_identity": (
                        IMAGING_SAFETY.measurement_set_instance_identity(
                            product
                        )
                    ),
                }
                for target, product in split_report.items()
            },
        }
        # Recorded only when it changed something, so the report of a band
        # that is not averaged reads exactly as it did before.
        if channel_average["width"] > 1:
            self.write_atomic_json(
                self.target_channel_average_path(),
                channel_average,
            )
        self.write_atomic_json(
            self.caldir / "split_targets.json",
            split_report,
        )
        self.write_atomic_json(provenance_file, split_provenance)

        self.mark_checkpoint(
            "split"
        )

    # =========================================================
    # QA PLOTS
    # =========================================================

    def generate_native_qa_plots(self):
        """Render essential gain/bandpass diagnostics without CASA plotms."""

        status = {"plots": [], "warnings": []}
        plot_specs = [
            (
                "bandpass.cal",
                "bandpass_amplitude.png",
                "Bandpass amplitude (median over solutions)",
                lambda table, outfile, title: QA_PRODUCTS.plot_caltable_bandpass(
                    table, outfile, "amplitude", title
                ),
            ),
            (
                "bandpass.cal",
                "bandpass_phase.png",
                "Bandpass phase (circular mean over solutions)",
                lambda table, outfile, title: QA_PRODUCTS.plot_caltable_bandpass(
                    table, outfile, "phase", title
                ),
            ),
            (
                "phase.cal",
                "phase_gains.png",
                "Phase-calibration gains",
                lambda table, outfile, title: QA_PRODUCTS.plot_caltable_time_series(
                    table, outfile, "phase", title
                ),
            ),
            (
                "phase_scan.cal",
                "phase_scan_gains.png",
                "Scan-averaged phase gains transferred to the targets",
                lambda table, outfile, title: QA_PRODUCTS.plot_caltable_time_series(
                    table, outfile, "phase", title
                ),
            ),
            (
                "amp.cal",
                "amp_gains.png",
                "Amplitude-calibration gains",
                lambda table, outfile, title: QA_PRODUCTS.plot_caltable_time_series(
                    table, outfile, "amplitude", title
                ),
            ),
            (
                "flux.cal",
                "flux_gains.png",
                "Flux-scaled gains",
                lambda table, outfile, title: QA_PRODUCTS.plot_caltable_time_series(
                    table, outfile, "amplitude", title
                ),
            ),
        ]

        for table_name, filename, title, plotter in plot_specs:
            caltable = self.caldir / table_name

            if not caltable.exists():
                status["warnings"].append(
                    f"Calibration table is unavailable: {caltable}"
                )
                continue

            try:
                outfile = plotter(caltable, self.plotdir / filename, title)
            except Exception as exc:
                status["warnings"].append(
                    f"Could not create {filename}: {exc}"
                )
                continue

            if outfile is None:
                status["warnings"].append(
                    f"No usable data were found for {filename}."
                )
            else:
                status["plots"].append(str(outfile))

        status_path = self.caldir / "qa_plot_status.json"
        QA_PRODUCTS.write_status(status_path, status)

        if status["plots"]:
            print(
                "Calibration QA plots written: " + ", ".join(status["plots"]),
                flush=True,
            )
        else:
            print(
                "WARNING: No calibration QA plot was created. See "
                f"{status_path}.",
                flush=True,
            )

        return status

    def generate_qa_plots(self, force=False):

        status_path = self.caldir / "qa_plot_status.json"

        if self.checkpoint_exists("qa") and status_path.exists() and not force:

            print("QA plots already generated.")
            return

        print("Generating calibration QA plots.")

        # CASA 6.7 on ilifu does not expose plotms through casatasks. The
        # table-based renderer is therefore the authoritative QA path.
        self.generate_native_qa_plots()
        self.mark_checkpoint("qa")
        return

        if not CASA_COMPAT.task_available(plotms):

            print(
                "CASA plotms is unavailable in this installation; "
                "skipping optional calibration QA PNGs."
            )
            self.mark_checkpoint("qa")
            return

        qa_plots = [

            {
                "caltable": "delay.cal",
                "xaxis": "antenna",
                "yaxis": "phase",
                "plotfile": "delay_phase.png"
            },

            {
                "caltable": "bandpass.cal",
                "xaxis": "freq",
                "yaxis": "amp",
                "plotfile": "bandpass_amp.png"
            },

            {
                "caltable": "bandpass.cal",
                "xaxis": "freq",
                "yaxis": "phase",
                "plotfile": "bandpass_phase.png"
            },

            {
                "caltable": "phase.cal",
                "xaxis": "time",
                "yaxis": "phase",
                "plotfile": "phase_gains.png"
            },

            {
                "caltable": "amp.cal",
                "xaxis": "time",
                "yaxis": "amp",
                "plotfile": "amp_gains.png"
            }
        ]

        for p in qa_plots:

            caltable = (
                self.caldir /
                p["caltable"]
            )

            if not caltable.exists():

                continue

            outfile = (
                self.plotdir /
                p["plotfile"]
            )

            print(
                f"  plotms → {outfile.name}"
            )

            # plotms(vis="flux.cal", xaxis="antenna", yaxis="amp", showgui=False, highres=True, overwrite=True, plotfile=str(outfile))

            plotms(
                vis=str(caltable),

                xaxis=p["xaxis"],

                yaxis=p["yaxis"],

                showgui=False,

                highres=True,

                overwrite=True,

                plotfile=str(outfile)
            )
                #
        # Flux calibration QA
        #

        fluxcal_table = self.caldir / "flux.cal"

        if fluxcal_table.exists():

            #
            # Flux amplitude vs antenna
            #

            flux_amp = (
                self.plotdir /
                "flux_amp.png"
            )

            print(
                f"  plotms → {flux_amp.name}"
            )

            plotms(
                vis=str(fluxcal_table),

                xaxis="antenna",

                yaxis="amp",

                showgui=False,

                highres=True,

                overwrite=True,

                plotfile=str(flux_amp)
            )

            #
            # Flux amplitude vs time
            #

            flux_time = (
                self.plotdir /
                "flux_time.png"
            )

            print(
                f"  plotms → {flux_time.name}"
            )

            plotms(
                vis=str(fluxcal_table),

                xaxis="time",

                yaxis="amp",

                showgui=False,

                highres=True,

                overwrite=True,

                plotfile=str(flux_time)
            )
        '''    
        #
        # Bandpass calibration QA
        #

        bp_table = self.caldir / "bandpass.cal"

        if bp_table.exists():

            #
            # Bandpass amplitude vs channels
            #

            bp_amp = (
                self.plotdir /
                "bandpass_amp.png"
            )

            print(
                f"  plotms → {bp_amp.name}"
            )

            plotms(
                vis=str(bp_table),

                xaxis="channel",

                yaxis="amp",

                showgui=False,

                highres=True,

                overwrite=True,

                plotfile=str(bp_amp)
            )

            #
            # Bandpass phase vs channels
            #

            bp_phase = (
                self.plotdir /
                "bandpass_phase.png"
            )

            print(
                f"  plotms → {bp_phase.name}"
            )

            plotms(
                vis=str(bp_table),

                xaxis="channel",

                yaxis="phase",

                showgui=False,

                highres=True,

                overwrite=True,

                plotfile=str(bp_phase)
            )
        '''    
        
        #
        # Phasing & Amplitude calibration QA
        #

        phase_table = self.caldir / "phase.cal"

        amp_table = self.caldir / "amp.cal"

        if phase_table.exists():

            #
            # Phase gains vs time
            #

            phase_time = (
                self.plotdir /
                "phase_gains.png"
            )

            print(
                f"  plotms → {phase_time.name}"
            )

            plotms(
                vis=str(phase_table),

                xaxis="time",

                yaxis="phase",

                showgui=False,

                highres=True,

                overwrite=True,

                plotfile=str(phase_time)
            )

        if amp_table.exists():

            #
            # Amplitude gains vs time
            #

            amp_time = (
                self.plotdir /
                "amp_gains.png"
            )

            print(
                f"  plotms → {amp_time.name}"
            )

            plotms(
                vis=str(amp_table),

                xaxis="time",

                yaxis="amp",

                showgui=False,

                highres=True,

                overwrite=True,

                plotfile=str(amp_time)
            )
        
        #
        # Final corrected amplitude vs uvdist
        #

        uvdist_amp = (
            self.plotdir /
            "uvdist_amp.png"
        )

        print(
            f"  plotms → {uvdist_amp.name}"
        )

        plotms(
            vis=str(self.ms),
            xaxis="uvdist",
            yaxis="amp",
            ydatacolumn="CORRECTED",
            selectdata=False,
            showgui=False,
            highres=True,
            overwrite=True,
            plotfile=str(uvdist_amp)
        )
        
        self.mark_checkpoint("qa")

    # =========================================================
    # CALIBRATION REPORT
    # =========================================================

    def collect_caltable_qa(self):

        table_names = {
            "delay": "delay.cal",
            "prebp": "prebp.cal",
            "bandpass": "bandpass.cal",
            "phase": "phase.cal",
            "phase_scan": "phase_scan.cal",
            "amp": "amp.cal",
            "flux": "flux.cal",
        }
        collected = {}

        for role, filename in table_names.items():
            valid, reason = self.caltable_status(filename)
            sidecar = self.caltable_qa_path(filename)

            with open(sidecar) as handle:
                table_qa = json.load(handle)

            table_qa["sidecar"] = str(sidecar)
            collected[role] = table_qa

            if not valid:
                raise RuntimeError(
                    f"Calibration-table QA failed for {filename}: {reason}"
                )

        return collected

    # ---------------------------------------------------------

    def statistical_weights_record(self):
        """What statwt did to each split target, from its markers."""

        try:
            with open(self.caldir / "split_targets.json") as handle:
                targets = json.load(handle)
        except (OSError, ValueError, TypeError):
            return {}
        record = {}
        for target, product in (targets or {}).items():
            marker = self.statwt_marker(product)
            try:
                with open(marker) as handle:
                    record[target] = json.load(handle)
            except (OSError, ValueError, TypeError):
                record[target] = {
                    "status": (
                        "missing"
                        if self.statistical_weights_settings()["enabled"]
                        else "disabled"
                    )
                }
        return record

    def calibration_stability(self):
        """Summarise, in numbers, what the gain tables say about the data.

        Two things the scatter plots cannot show: how far each antenna's
        phase moves between visits to the phase calibrator, which bounds
        the interpolation error left in the target, and which antennas'
        gains sit far from the rest of the array.
        """

        try:
            summary = QA_PRODUCTS.caltable_stability_summary(
                self.caldir / "phase_scan.cal",
                self.caldir / "flux.cal",
            )
        except Exception as exc:
            return {"error": str(exc)}

        names = summary.get("field_names") or []
        for field, values in (summary.get("phase_transfer") or {}).items():
            name = names[field] if 0 <= field < len(names) else f"field {field}"
            values["field_name"] = name
            print(
                f"Phase stability on {name}: between visits a median "
                f"{values['median_gap_min']:.0f} min apart, antenna phases "
                f"moved a median of {values['median_abs_change_deg']:.0f} deg "
                f"(90th percentile {values['p90_abs_change_deg']:.0f} deg). "
                "The target's phases are interpolated across those gaps, so "
                "errors of this order stay in it until self-calibration "
                "removes them.",
                flush=True,
            )

        gains = summary.get("antenna_gains") or {}
        outliers = gains.get("outliers") or []
        if outliers:
            listed = ", ".join(
                f"{item.get('name') or item['antenna']} pol "
                f"{item['correlation']} ({item['median_gain']:.2f}, "
                f"{item['ratio_to_array']:.2f}x)"
                for item in outliers
            )
            print(
                "Antenna gains far from the array median of "
                f"{gains['array_median']:.2f}: {listed}. A gain this far "
                "off usually means an antenna or receiver chain not working "
                "like the others; statistical weights down-weight its noise, "
                "but inspect it and flag it if its data are also unstable.",
                flush=True,
            )
        return summary

    def write_report(self):

        caltable_qa = self.collect_caltable_qa()
        flux_model_provenance = self.load_flux_model_provenance()
        try:
            with open(
                self.caldir / "split_targets.provenance.json"
            ) as handle:
                split_target_identities = json.load(handle)
        except (OSError, ValueError, TypeError) as exc:
            raise RuntimeError(
                "Split-target provenance is missing or unreadable."
            ) from exc

        report = {

            "provenance_fingerprint": getattr(
                self,
                "provenance_fingerprint",
                None,
            ),

            "recipe_fingerprint": getattr(
                self,
                "recipe_fingerprint",
                None,
            ),

            "original_measurement_set":
                self.report["measurement_set"],

            "measurement_set_used":
                self.ms,

            "measurement_set_stat_sha256":
                IMAGING_SAFETY.path_stat_fingerprint(self.ms),

            "telescope": self.telescope,

            "band": self.band,

            "flux_calibrator":
                self.get_flux_calibrator(),

            "flux_calibrators":
                self.get_flux_calibrators(),

            "secondary_flux_calibrator":
                self.get_secondary_flux_calibrator(),

            "bandpass_calibrator":
                self.get_bandpass_calibrator(),

            "phase_calibrator":
                self.get_phase_calibrator(),

            "phase_calibrators":
                self.get_phase_calibrators(),

            "target_phase_map":
                self.get_target_phase_map(),

            "targets":
                self.get_targets(),

            "reference_antenna":
                self.get_reference_antenna(),

            "bad_antennas":
                self.flagging.get(
                    "bad_antennas",
                    []
                ),

            "delay_quarantine": getattr(
                self,
                "delay_quarantine",
                {"schema_version": 1, "antennas": []},
            ),

            "hanning_smoothing":
                self.band_config[
                    "hanningsmooth"
                ],

            "calibration_parameters": {
                "flux_standard": self.task_config["flux_standard"],
                "flux_model": self.flux_model,
                "flux_model_type": getattr(
                    self,
                    "flux_model_type",
                    "casa-image",
                ),
                "solution_rounds": getattr(
                    self,
                    "calibration_rounds",
                    1,
                ),
                "postcal_flagging": bool(
                    self.task_config.get("postcal_flagging", False)
                ),
                "bandpass_solnorm": self.task_config[
                    "bandpass_solnorm"
                ],
                "prebp_gaintype": self.task_config[
                    "prebp_gaintype"
                ],
                "prebp_applied_downstream": False,
                "refantmode": self.task_config.get(
                    "refantmode",
                    "strict",
                ),
                "fluxscale_fitorder": self.task_config.get(
                    "fluxscale_fitorder",
                    1,
                ),
                "target_phase_transfer_solint": self.task_config.get(
                    "target_phase_transfer_solint",
                    "inf",
                ),
                # Descriptive fields: the solve paths index these directly
                # and fail far earlier if a band profile is genuinely
                # incomplete, so the report must not abort over them.
                "phase_solint": self.band_config.get("phase_solint"),
                "amp_solint": self.band_config.get("amp_solint"),
                "calibrator_gaintype_switch": self.read_gaintype_switch(),
                "max_caltable_flagged_fraction": self.task_config[
                    "max_caltable_flagged_fraction"
                ],
                "max_caltable_group_flagged_fraction": self.task_config[
                    "max_caltable_group_flagged_fraction"
                ],
                "applymode": self.task_config["applymode"],
                "calwt": self.task_config["calwt"],
                "parang": self.task_config["parang"],
            },

            "calibration_tables": {

                "delay":
                    str(
                        self.caldir /
                        "delay.cal"
                    ),

                "prebp":
                    str(
                        self.caldir /
                        "prebp.cal"
                    ),

                "bandpass":
                    str(
                        self.caldir /
                        "bandpass.cal"
                    ),

                "phase":
                    str(
                        self.caldir /
                        "phase.cal"
                    ),

                "phase_scan":
                    str(
                        self.caldir /
                        "phase_scan.cal"
                    ),

                "amp":
                    str(
                        self.caldir /
                        "amp.cal"
                    ),

                "flux":
                    str(
                        self.caldir /
                        "flux.cal"
                    )
            },

            "calibration_table_qa": caltable_qa,

            "split_target_identities": split_target_identities,

            "statistical_weights": self.statistical_weights_record(),

            "calibration_stability": self.calibration_stability(),

            "flux_model_provenance": flux_model_provenance,

            "flux_model_provenance_path": str(
                self.flux_model_provenance_path()
            ),

            "qa_plot_status": str(
                self.caldir / "qa_plot_status.json"
            ),

            "checkpoints": {

                "setjy":
                    self.checkpoint_exists(
                        "setjy"
                    ),

                "delay":
                    self.checkpoint_exists(
                        "delay"
                    ),

                "prebp":
                    self.checkpoint_exists(
                        "prebp"
                    ),

                "bandpass":
                    self.checkpoint_exists(
                        "bandpass"
                    ),

                f"{self.telescope}_postcal_flagging":
                    self.checkpoint_exists(
                        f"{self.telescope}_postcal_flagging"
                    ),

                f"{self.telescope}_round2_complete":
                    self.checkpoint_exists(
                        f"{self.telescope}_round2_complete"
                    ),

                "phase":
                    self.checkpoint_exists(
                        "phase"
                    ),

                "phase_scan":
                    self.checkpoint_exists(
                        "phase_scan"
                    ),

                "amp":
                    self.checkpoint_exists(
                        "amp"
                    ),

                "fluxscale":
                    self.checkpoint_exists(
                        "fluxscale"
                    ),

                "applycal":
                    self.checkpoint_exists(
                        "applycal"
                    ),

                "split":
                    self.checkpoint_exists(
                        "split"
                    ),

                "qa":
                    self.checkpoint_exists(
                        "qa"
                    )
            }
        }

        try:
            with open(self.target_channel_average_path()) as handle:
                report["target_channel_average"] = json.load(handle)
        except (OSError, ValueError, TypeError):
            pass

        # What was bridged for the targets, without the per-solution list.
        try:
            with open(self.caldir / "target_transfer_bridge.json") as handle:
                bridge = json.load(handle)
            report["target_transfer_bridge"] = {
                key: (
                    {
                        role: {
                            name: value for name, value in record.items()
                            if name != "fills"
                        }
                        for role, record in bridge["tables"].items()
                    }
                    if key == "tables" else bridge[key]
                )
                for key in bridge
            }
        except (OSError, ValueError, TypeError, KeyError, AttributeError):
            pass

        outfile = (
            self.caldir /
            "calibration_report.json"
        )

        with open(outfile, "w") as f:

            json.dump(
                report,
                f,
                indent=4
            )

        print(
            f"Calibration report written to:\n"
            f"{outfile}"
        )

    # =========================================================
    # SPLIT-ONLY REPAIR
    # =========================================================

    def repair_split_targets(self):
        """Rebuild only the calibrated target MS products.

        The orchestrator selects this path when every calibration solution
        and checkpoint is still valid but a split target MS was moved,
        replaced, or left incomplete.  No solver is re-run, so the accepted
        gain/bandpass/flux solutions and their provenance are preserved.
        """

        print("\n")
        print("=" * 70)
        print("Repairing calibrated split-target products only")
        print("=" * 70)

        if not self.checkpoint_exists("applycal"):
            raise RuntimeError(
                "Split-target repair requires a completed applycal stage. "
                "Rerun calibration without --repair-split-targets-only."
            )

        opened = False
        try:
            tb.open(self.ms)
            opened = True
            columns = set(tb.colnames())
        finally:
            if opened:
                tb.close()

        if "CORRECTED_DATA" not in columns:
            raise RuntimeError(
                "Split-target repair needs CORRECTED_DATA in "
                f"{self.ms}; rerun the full calibration stage."
            )

        # collect_caltable_qa raises if any accepted table is now unusable,
        # so a repair can never publish targets from broken solutions.
        self.collect_caltable_qa()

        # The rewritten report must keep recording the solution rounds that
        # actually produced these tables. Defaulting to one round here would
        # make the orchestrator reject its own repaired stage on the next
        # resume, because two rounds are required for uGMRT and MeerKAT.
        self.calibration_rounds = 1
        if self.task_config.get("postcal_flagging", False):
            names = self.postcal_checkpoint_names()
            if not self.checkpoint_exists(names["complete"]):
                raise RuntimeError(
                    "Split-target repair requires the completed round-2 "
                    f"checkpoint {names['complete']}."
                )
            self.calibration_rounds = 2

        self.split_targets()
        self.generate_qa_plots()
        self.write_report()
        self.mark_checkpoint("calibration")

        print("\n")
        print("=" * 70)
        print("Split-target repair completed")
        print("=" * 70)

    # =========================================================
    # MAIN CALIBRATION PIPELINE
    # =========================================================

    def run(self):

        print("\n")
        print("="*70)
        print(
            f"Starting {TELESCOPE_PROFILE_MODULE.get_telescope_profile(self.telescope)['label']} "
            "Calibration Pipeline"
        )
        print("="*70)

        self.print_summary()

        #
        # Flux density model
        #
        self.run_setjy()

        #
        # Delay calibration
        #
        self.run_delaycal()

        #
        # Pre-bandpass phase calibration
        #
        self.run_prebp()

        #
        # Bandpass calibration
        #
        self.run_bandpass()

        #
        # Phase calibration
        #
        self.run_phasecal()

        #
        # Scan-averaged phase solutions transferred to the science targets
        #
        self.run_phase_transfer_cal()

        #
        # Amplitude calibration
        #
        self.run_ampcal()

        #
        # Flux scale transfer
        #
        self.run_fluxscale()

        #
        # Apply solutions
        #
        # With a second round to come, these solutions are provisional: the
        # flags they are about to write are taken back before round 2.
        self.save_flags_before_round1_applycal()
        self.run_applycal()

        # Median-based RFlag is safest after a provisional bandpass has
        # flattened the spectrum.  For configured telescopes, accepted new
        # flags are followed by a complete fresh set of calibration solutions.
        self.run_second_round()

        #
        # Split calibrated targets
        #
        self.split_targets()

        #
        # QA products
        #
        self.generate_qa_plots()

        #
        # Write calibration report
        #
        self.write_report()

        #
        # Final checkpoint
        #
        
        self.mark_checkpoint("calibration")

        print("\n")
        print("="*70)
        print(
            f"{TELESCOPE_PROFILE_MODULE.get_telescope_profile(self.telescope)['label']} "
            "Calibration Completed"
        )
        print("="*70)

        print(
            f"Outputs written to:\n"
            f"{self.caldir}"
        )

# =========================================================
# CASA ENTRY POINT
# =========================================================


# An interrupt reaches a batch stage when the job is cancelled, when it runs
# into its time limit, or when an interactive session ends. CASA would print
# a bare KeyboardInterrupt traceback from wherever it happened to be, which
# reads like a crash. Exit code 130 is the conventional one for SIGINT, and
# the orchestrator reports it as an interruption rather than a failure.
def _interrupted(stage):
    print(
        f"\n{stage} was interrupted from outside the pipeline (SIGINT): the "
        "job was cancelled, reached its time limit, or its session ended. "
        "Completed checkpoints are kept; resubmitting continues from them.",
        flush=True,
    )
    sys.exit(130)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run uGMRT or MeerKAT CASA continuum calibration."
    )

    parser.add_argument(
        "ms",
        help="Input Measurement Set.",
    )

    parser.add_argument(
        "--outdir",
        default=RUN_CONFIG.get(
            "outdir",
            "Products",
        ),
        help="Pipeline product directory.",
    )

    parser.add_argument(
        "--flux-calibrator",
        default=RUN_CONFIG.get("flux_calibrator"),
        help=(
            "Explicit primary flux-calibrator field or comma-separated "
            "list. The first field is the primary absolute reference."
        ),
    )

    parser.add_argument(
        "--secondary-flux-calibrator",
        default=RUN_CONFIG.get("secondary_flux_calibrator"),
        help=(
            "Optional additional flux-calibrator field or comma-separated "
            "list appended after --flux-calibrator."
        ),
    )

    parser.add_argument(
        "--bandpass-calibrator",
        default=RUN_CONFIG.get("bandpass_calibrator"),
        help=(
            "Explicit delay/bandpass calibrator. Defaults to the unique "
            "BANDPASS-intent field or the primary flux calibrator."
        ),
    )

    parser.add_argument(
        "--phase-calibrator",
        default=RUN_CONFIG.get("phase_calibrator"),
        help=(
            "Explicit phase-calibrator field or comma-separated list."
        ),
    )

    parser.add_argument(
        "--target-phase-map",
        default=RUN_CONFIG.get("target_phase_map"),
        metavar="TARGET:PHASE[,TARGET:PHASE...]",
        help=(
            "Explicit target-to-phase-calibrator mapping used when "
            "multiple phase calibrators are present."
        ),
    )

    parser.add_argument(
        "--targets",
        default=RUN_CONFIG.get("target_fields"),
        metavar="FIELD[,FIELD...]",
        help="Explicit comma-separated target fields.",
    )

    parser.add_argument(
        "--refant",
        default=RUN_CONFIG.get("reference_antenna"),
        help="Explicit reference antenna.",
    )

    parser.add_argument(
        "--telescope",
        choices=TELESCOPE_PROFILE_MODULE.TELESCOPE_CHOICES,
        default=RUN_CONFIG.get("telescope", "auto"),
        help="Verify the telescope identity recorded by inspection.",
    )

    parser.add_argument(
        "--flux-standard",
        default=RUN_CONFIG.get("flux_standard"),
        help="Override the telescope profile's CASA setjy flux standard.",
    )

    parser.add_argument(
        "--flux-model",
        default=RUN_CONFIG.get("flux_model"),
        help=(
            "Optional validated CASA calibrator model image. Required for "
            "MeerKAT UHF J1939-6342/J0408-6545 reductions."
        ),
    )

    parser.add_argument(
        "--flux-model-type",
        choices=["casa-image", "component-list"],
        default=RUN_CONFIG.get("flux_model_type", "casa-image"),
        help=(
            "Interpret --flux-model as a CASA image or component list for "
            "insertion with ft."
        ),
    )

    parser.add_argument(
        "--provenance-fingerprint",
        default=None,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--recipe-fingerprint",
        default=None,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--repair-split-targets-only",
        action="store_true",
        help=(
            "Rebuild only the calibrated split-target Measurement Sets "
            "from existing, validated calibration solutions."
        ),
    )

    parser.add_argument(
        "--plots-only",
        action="store_true",
        help=(
            "Regenerate calibration QA plots from existing calibration "
            "tables without recalibrating the Measurement Set."
        ),
    )

    args = parser.parse_args()

    cal = CalibrationManager(
        args.ms,
        outdir=args.outdir,
        flux_calibrator=args.flux_calibrator,
        bandpass_calibrator=args.bandpass_calibrator,
        secondary_flux_calibrator=args.secondary_flux_calibrator,
        phase_calibrator=args.phase_calibrator,
        target_phase_map=args.target_phase_map,
        target_fields=args.targets,
        reference_antenna=args.refant,
        telescope=args.telescope,
        flux_standard=args.flux_standard,
        flux_model=args.flux_model,
        flux_model_type=args.flux_model_type,
        provenance_fingerprint=args.provenance_fingerprint,
        recipe_fingerprint=args.recipe_fingerprint,
    )

    if args.plots_only:
        cal.generate_qa_plots(force=True)
    elif args.repair_split_targets_only:
        try:
            cal.repair_split_targets()
        except KeyboardInterrupt:
            _interrupted("Calibration")
    else:
        try:
            cal.run()
        except KeyboardInterrupt:
            _interrupted("Calibration")
    
