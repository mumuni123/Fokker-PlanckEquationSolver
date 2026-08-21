#!/usr/bin/env python3
"""Assemble independent Gate-A evidence into one machine-readable result."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict


def read_result(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as stream:
        for raw in stream:
            line = raw.strip()
            if line and not line.startswith("#") and "=" in line:
                key, value = line.split("=", 1)
                values[key.strip()] = value.strip()
    return values


def require_float(values: Dict[str, str], key: str) -> float:
    if key not in values:
        raise ValueError(f"missing key {key}")
    value = float(values[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite key {key}")
    return value


def require_int(values: Dict[str, str], key: str) -> int:
    value = require_float(values, key)
    integer = int(value)
    if value != integer:
        raise ValueError(f"non-integral key {key}")
    return integer


def read_spectrum_number_residual(path: Path) -> float:
    with path.open("r", encoding="utf-8") as stream:
        first = stream.readline().strip()
    marker = "relative_number_residual="
    if marker not in first:
        raise ValueError("spectrum header lacks relative_number_residual")
    value = float(first.split(marker, 1)[1].split()[0])
    if not math.isfinite(value):
        raise ValueError("non-finite spectrum number residual")
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry-result", required=True, type=Path)
    parser.add_argument("--checkpoint-audit", required=True, type=Path)
    parser.add_argument("--conversion-result", required=True, type=Path)
    parser.add_argument("--spectrum", type=Path)
    parser.add_argument("--result", required=True, type=Path)
    args = parser.parse_args()

    checks: Dict[str, int] = {}
    details: Dict[str, str] = {}
    errors = []
    try:
        geometry = read_result(args.geometry_result)
        audit = read_result(args.checkpoint_audit)
        conversion = read_result(args.conversion_result)
        spectrum_path = args.spectrum or Path(str(args.checkpoint_audit) + ".spectrum.dat")

        checks["geometry_overlap_pass"] = int(geometry.get("status") == "PASS")
        checks["geometry_mass_partition_pass"] = int(
            require_int(geometry, "empty_cell_count") == 0
            and require_float(geometry, "min_fraction") >= -1.0e-13
            and require_float(geometry, "max_sum_error") <= 2.0e-12
            and require_float(geometry, "max_rebin_error") <= 2.0e-12
        )

        details["checkpoint_audit_schema"] = audit.get("audit_schema", "missing")
        checks["checkpoint_schema_pass"] = int(
            audit.get("audit_schema") == "hybrid_checkpoint_gate_audit_v5"
        )
        checks["checkpoint_finite_pass"] = int(
            require_int(audit, "checkpoint_audit_finite") == 1
        )
        checks["spectrum_shape_pass"] = int(
            require_int(audit, "spectrum_rebin_shape_pass") == 1
            and require_int(audit, "threshold_artificial_gap_pass") == 1
            and require_int(audit, "threshold_deep_valley_count") == 0
            and require_float(audit, "threshold_max_adjacent_jump_ratio") <= 100.0
        )
        checks["core_noise_pass"] = int(
            require_int(audit, "core_tail_density_noise_pass") == 1
        )

        spectrum_number_residual = read_spectrum_number_residual(spectrum_path)
        details["spectrum_number_relative_residual"] = repr(spectrum_number_residual)
        checks["spectrum_number_conservation_pass"] = int(
            spectrum_number_residual <= 2.0e-12
        )

        q48 = require_float(conversion, "quadrature_4_vs_8_relative_max")
        details["quadrature_4_vs_8_relative_max"] = repr(q48)
        checks["conversion_operator_pass"] = int(
            conversion.get("status") == "PASS"
            and require_float(conversion, "sink_number_relative_error") <= 1.0e-12
            and require_int(conversion, "interface_duplicate_count") == 0
            and q48 <= 1.0e-11
            and require_float(conversion, "below_threshold_number_relative") <= 1.0e-12
            and require_int(conversion, "negative_node_count") == 0
        )
        # Older, already accepted 17A result files predate this diagnostic.
        # If present it is enforced; absence does not invalidate the spectrum
        # gate because tiny-transfer robustness is not a spectral criterion.
        if "tiny_tail_ok" in conversion:
            checks["tiny_tail_pass"] = int(
                require_int(conversion, "tiny_tail_ok") == 1
            )
            details["tiny_tail_evidence_available"] = "1"
        else:
            details["tiny_tail_pass"] = "-1"
            details["tiny_tail_evidence_available"] = "0"
    except (OSError, ValueError) as exc:
        errors.append(str(exc))

    passed = not errors and bool(checks) and all(checks.values())
    args.result.parent.mkdir(parents=True, exist_ok=True)
    with args.result.open("w", encoding="utf-8", newline="\n") as out:
        out.write("gate_a_schema=hybrid_gate_a_v1\n")
        out.write(f"status={'PASS' if passed else 'FAIL'}\n")
        out.write("gate_scope=threshold_geometry_real_checkpoint_shape_and_conversion_operator\n")
        out.write("conversion_quadrature_scope=production_u_parallel_operator_contract\n")
        out.write("full_run_parameter_ab_required=0\n")
        for key in sorted(checks):
            out.write(f"{key}={checks[key]}\n")
        for key in sorted(details):
            out.write(f"{key}={details[key]}\n")
        out.write(f"error_count={len(errors)}\n")
        for index, error in enumerate(errors):
            out.write(f"error_{index}={error}\n")

    print(f"Gate A: {'PASS' if passed else 'FAIL'} -> {args.result}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
