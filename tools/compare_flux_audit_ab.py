#!/usr/bin/env python3
"""Compare static-cell and read-only flux-audit runs.

The audit mode is allowed to create extra diagnostics, but it must not change
the accepted physical state.  This script deliberately refuses to claim a
flux-closure PASS unless the audit output contains an independent face ledger.
"""

from __future__ import print_function

import argparse
import hashlib
import math
import os
import sys


AUDIT_NAMES = {
    "bulk_tail_flux_accepted_steps.dat",
    "bulk_tail_flux_face_accepted_steps.dat",
    "trial_bulk_tail_flux_failures.dat",
    "bulk_tail_moment_audit_accepted_steps.dat",
    "bulk_tail_moment_audit_velocity.dat",
    "collision_flux_accepted_steps.dat",
    "collision_flux_substeps_accepted.dat",
    "conversion_source_accepted_steps.dat",
    "vpfp_step_diagnostics.dat",
    "vpfp_failure.dat",
    "run.out",
    "run.err",
    "manifest.txt",
}

TIMING_FIELDS = (
    "vlasov_s", "field_s", "beam_s", "collision_s", "mpi_collective_s",
    "wall_tail_push_s", "wall_tail_deposit_s", "wall_tail_migrate_s",
    "wall_conversion_s", "wall_diagnostics_s",
)


def accepted_run_info(root):
    """Return the accepted-state interval recorded by a production run."""
    path = os.path.join(root, "vpfp_step_diagnostics.dat")
    if not os.path.isfile(path):
        return None
    accepted = []
    with open(path, "r") as stream:
        header = None
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            if header is None:
                header = fields
                continue
            try:
                step = int(fields[0])
                time_s = float(fields[1])
                accepted_flag = int(fields[2])
            except (ValueError, IndexError):
                continue
            if accepted_flag == 1:
                accepted.append((step, time_s))
    if not accepted:
        return {
            "count": 0,
            "first_step": None,
            "first_time_s": math.nan,
            "last_step": None,
            "last_time_s": math.nan,
        }
    return {
        "count": len(accepted),
        "first_step": accepted[0][0],
        "first_time_s": accepted[0][1],
        "last_step": accepted[-1][0],
        "last_time_s": accepted[-1][1],
    }


def failure_summary(root):
    """Return whether the run recorded a rejected terminal VPFP step."""
    path = os.path.join(root, "vpfp_failure.dat")
    if not os.path.isfile(path):
        return False, ""
    with open(path, "r") as stream:
        records = [line.strip() for line in stream if line.strip()]
    if not records:
        return False, ""
    return True, records[-1]


def same_float(a, b):
    if not (math.isfinite(a) and math.isfinite(b)):
        return False
    return abs(a - b) <= 32.0 * sys.float_info.epsilon * max(
        abs(a), abs(b), 1.0e-30)


def physical_files(root):
    files = {}
    for base, _, names in os.walk(root):
        for name in names:
            if name in AUDIT_NAMES or name.startswith("trial_"):
                continue
            path = os.path.join(base, name)
            rel = os.path.relpath(path, root)
            lower = rel.lower()
            # Hash every accepted physical state artifact, not only the
            # profile files whose name happens to contain snapshot.  Audit
            # ledgers and timing files are deliberately excluded above.
            if not ("snapshot" in lower or "checkpoint" in lower or
                    os.path.basename(lower) in {
                        "state.bin", "fields.bin", "electrons.bin",
                        "beam.bin", "tail.bin"}):
                continue
            files[rel] = path
    return files


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def state_hash_equal(static_root, audit_root):
    static_files = physical_files(static_root)
    audit_files = physical_files(audit_root)
    static_names = set(static_files)
    audit_names = set(audit_files)
    common = sorted(static_names.intersection(audit_names))
    if not common or static_names != audit_names:
        return False, 0
    for rel in common:
        if sha256(static_files[rel]) != sha256(audit_files[rel]):
            return False, len(common)
    return True, len(common)


def read_flux_rows(root):
    path = os.path.join(root, "bulk_tail_flux_accepted_steps.dat")
    if not os.path.isfile(path):
        return None
    with open(path, "r") as stream:
        lines = [line.split() for line in stream if line.strip()]
    if len(lines) < 2:
        return []
    header = lines[0]
    indices = {name: i for i, name in enumerate(header)}
    rows = []
    for fields in lines[1:]:
        if len(fields) != len(header):
            continue
        row = {}
        for name, index in indices.items():
            try:
                row[name] = float(fields[index])
            except (ValueError, IndexError):
                row[name] = math.nan
        rows.append(row)
    return rows


def total_wall_seconds(root):
    path = os.path.join(root, "vpfp_step_diagnostics.dat")
    if not os.path.isfile(path):
        return math.nan
    with open(path, "r") as stream:
        lines = [line.split() for line in stream if line.strip()]
    if len(lines) < 2:
        return math.nan
    header = lines[0]
    try:
        index = header.index("wall_s")
    except ValueError:
        return math.nan
    total = 0.0
    count = 0
    for fields in lines[1:]:
        if len(fields) <= index:
            continue
        try:
            value = float(fields[index])
        except ValueError:
            continue
        if not math.isfinite(value) or value < 0.0:
            return math.nan
        total += value
        count += 1
    return total if count else math.nan


def timing_sums(root):
    """Return cumulative production timing columns without double-counting."""
    path = os.path.join(root, "vpfp_step_diagnostics.dat")
    totals = {name: math.nan for name in TIMING_FIELDS}
    if not os.path.isfile(path):
        return totals
    with open(path, "r") as stream:
        lines = [line.split() for line in stream if line.strip()]
    if len(lines) < 2:
        return totals
    indices = {name: index for index, name in enumerate(lines[0])}
    for name in TIMING_FIELDS:
        index = indices.get(name)
        if index is None:
            continue
        total = 0.0
        valid = True
        for fields in lines[1:]:
            if len(fields) <= index:
                valid = False
                break
            try:
                value = float(fields[index])
            except ValueError:
                valid = False
                break
            if not math.isfinite(value) or value < 0.0:
                valid = False
                break
            total += value
        if valid:
            totals[name] = total
    return totals


def read_face_rows(root):
    path = os.path.join(root, "bulk_tail_flux_face_accepted_steps.dat")
    if not os.path.isfile(path):
        return None
    with open(path, "r") as stream:
        lines = [line.split() for line in stream if line.strip()]
    if len(lines) < 2:
        return []
    header = lines[0]
    indices = {name: i for i, name in enumerate(header)}
    required = {
        "accepted_step", "time_s", "ix_global", "direction", "face_index",
        "transverse_index", "operator_stage", "face_export_number",
        "parcel_number",
    }
    if not required.issubset(set(indices)):
        return []
    rows = []
    for fields in lines[1:]:
        if len(fields) != len(header):
            continue
        row = {}
        for name, index in indices.items():
            try:
                row[name] = float(fields[index])
            except (ValueError, IndexError):
                row[name] = math.nan
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--static", required=True)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--result", required=True)
    parser.add_argument("--ignore-performance", action="store_true",
                        help="Do not make audit wall-time overhead a gate.")
    args = parser.parse_args()

    static_info = accepted_run_info(args.static)
    audit_info = accepted_run_info(args.audit)
    static_failed, static_failure = failure_summary(args.static)
    audit_failed, audit_failure = failure_summary(args.audit)
    static_count = None if static_info is None else static_info["count"]
    audit_count = None if audit_info is None else audit_info["count"]
    accepted_equal = (
        static_info is not None and audit_info is not None and
        static_info["count"] == audit_info["count"] and
        static_info["first_step"] == audit_info["first_step"] and
        static_info["last_step"] == audit_info["last_step"] and
        same_float(static_info["first_time_s"], audit_info["first_time_s"]) and
        same_float(static_info["last_time_s"], audit_info["last_time_s"])
    )
    comparison_reasons = []
    if static_info is None:
        comparison_reasons.append("static_missing_step_diagnostics")
    if audit_info is None:
        comparison_reasons.append("audit_missing_step_diagnostics")
    if static_info is not None and audit_info is not None and not accepted_equal:
        comparison_reasons.append("accepted_horizon_mismatch")
    if static_failed:
        comparison_reasons.append("static_vpfp_failure")
    if audit_failed:
        comparison_reasons.append("audit_vpfp_failure")
    comparison_valid = not comparison_reasons

    # Across-run hashes, flux closure and timing are meaningful only if both
    # jobs reached the identical accepted-state interval.  The per-call audit
    # gates below remain useful even when a run terminates early.
    equal, common_files = (False, 0)
    if comparison_valid:
        equal, common_files = state_hash_equal(args.static, args.audit)
    rows = read_flux_rows(args.audit)
    face_rows = read_face_rows(args.audit)
    rows_present = rows is not None and len(rows) > 0
    face_rows_present = face_rows is not None and len(face_rows) > 0
    compact_face_fields = {
        "face_audit_count", "face_audit_face_abs_sum",
        "face_audit_parcel_abs_sum", "face_audit_abs_error_sum",
        "face_audit_max_relative", "face_audit_abs_at_max_relative",
    }
    compact_face_summary = rows_present and all(
        compact_face_fields.issubset(set(row)) for row in rows or [])
    inplace_state_available = rows_present and all(
        "audit_inplace_state_bitwise_equal" in row for row in rows or [])
    inplace_rng_available = rows_present and all(
        "audit_inplace_rng_equal" in row for row in rows or [])
    inplace_ledger_available = rows_present and all(
        "audit_inplace_ledger_equal" in row for row in rows or [])
    inplace_state_ok = inplace_state_available and all(
        row.get("audit_inplace_state_bitwise_equal", math.nan) == 1.0
        for row in rows or [])
    inplace_rng_ok = inplace_rng_available and all(
        row.get("audit_inplace_rng_equal", math.nan) == 1.0
        for row in rows or [])
    inplace_ledger_ok = inplace_ledger_available and all(
        row.get("audit_inplace_ledger_equal", math.nan) == 1.0
        for row in rows or [])
    audit_valid_available = rows_present and all(
        "audit_valid" in row and "audit_failure_code" in row
        for row in rows or [])
    audit_valid_ok = audit_valid_available and all(
        row.get("audit_valid", math.nan) == 1.0 and
        row.get("audit_failure_code", math.nan) == 0.0
        for row in rows or [])

    duplicate_count = math.nan
    below_threshold = math.nan
    roundoff_discarded = math.nan
    independent_face_ledger = False
    parcel_face_max_relative = math.nan
    parcel_face_absolute_at_max_relative = math.nan
    parcel_face_relative = math.nan
    below_relative = math.nan
    face_duplicate_count = -1
    face_scale = math.nan
    parcel_scale = math.nan
    max_error_row = None
    static_wall = math.nan
    audit_wall = math.nan
    static_timings = {name: math.nan for name in TIMING_FIELDS}
    audit_timings = {name: math.nan for name in TIMING_FIELDS}
    if comparison_valid:
        duplicate_count = sum(
            row.get("duplicate_count", math.nan) for row in rows or [])
        below_threshold = sum(
            row.get("below_threshold_number", math.nan) for row in rows or [])
        roundoff_discarded = sum(
            row.get("roundoff_discarded_number", math.nan)
            for row in rows or [])
        if compact_face_summary:
            independent_face_ledger = all(
                math.isfinite(row["face_audit_face_abs_sum"]) and
                math.isfinite(row["face_audit_parcel_abs_sum"]) and
                math.isfinite(row["face_audit_abs_error_sum"]) and
                math.isfinite(row["face_audit_max_relative"])
                for row in rows or [])
            face_scale = sum(row["face_audit_face_abs_sum"]
                             for row in rows or [])
            parcel_scale = sum(row["face_audit_parcel_abs_sum"]
                               for row in rows or [])
            face_delta = sum(row["face_audit_abs_error_sum"]
                             for row in rows or [])
            if independent_face_ledger:
                parcel_face_max_relative = max(
                    row["face_audit_max_relative"] for row in rows or [])
                parcel_face_absolute_at_max_relative = max(
                    row["face_audit_abs_at_max_relative"]
                    for row in rows or [])
        else:
            independent_face_ledger = face_rows_present and all(
                math.isfinite(row.get("face_export_number", math.nan)) and
                math.isfinite(row.get("parcel_number", math.nan))
                for row in face_rows or [])
            face_scale = sum(abs(row.get("face_export_number", math.nan))
                             for row in face_rows or [])
            parcel_scale = sum(abs(row.get("parcel_number", math.nan))
                               for row in face_rows or [])
            face_delta = sum(abs(row.get("face_export_number", math.nan) -
                                 row.get("parcel_number", math.nan))
                             for row in face_rows or [])
        if independent_face_ledger and not compact_face_summary:
            face_errors = []
            for row in face_rows or []:
                face = row.get("face_export_number", math.nan)
                parcel = row.get("parcel_number", math.nan)
                if not (math.isfinite(face) and math.isfinite(parcel)):
                    face_errors = []
                    break
                absolute_error = abs(face - parcel)
                relative_error = absolute_error / max(
                    abs(face), abs(parcel), 1.0e-300)
                face_errors.append(relative_error)
                if (max_error_row is None or
                        relative_error > max_error_row["relative_error"]):
                    max_error_row = {
                        "accepted_step": int(row["accepted_step"]),
                        "time_s": row["time_s"],
                        "ix_global": int(row["ix_global"]),
                        "direction": int(row["direction"]),
                        "face_index": int(row["face_index"]),
                        "transverse_index": int(row["transverse_index"]),
                        "operator_stage": int(row["operator_stage"]),
                        "face_number": face,
                        "parcel_number": parcel,
                        "absolute_error": absolute_error,
                        "relative_error": relative_error,
                        "node_failure_reason": row.get(
                            "node_failure_reason", math.nan),
                        "reconstructed_target": row.get(
                            "reconstructed_target", math.nan),
                        "node_sum": row.get("node_sum", math.nan),
                    }
            if face_errors:
                parcel_face_max_relative = max(face_errors)
                parcel_face_absolute_at_max_relative = max_error_row["absolute_error"]
        elif compact_face_summary and face_rows_present:
            # Compact mode writes at most one worst-face record per accepted
            # step.  It is provenance only; the scalar summary above remains
            # the authoritative L1/Linf closure ledger.
            for row in face_rows or []:
                face = row.get("face_export_number", math.nan)
                parcel = row.get("parcel_number", math.nan)
                absolute_error = row.get("absolute_error", abs(face - parcel))
                relative_error = row.get(
                    "relative_error",
                    absolute_error / max(abs(face), abs(parcel), 1.0e-300))
                if not (math.isfinite(face) and math.isfinite(parcel) and
                        math.isfinite(absolute_error) and
                        math.isfinite(relative_error)):
                    continue
                if (max_error_row is None or
                        relative_error > max_error_row["relative_error"]):
                    max_error_row = {
                        "accepted_step": int(row["accepted_step"]),
                        "time_s": row["time_s"],
                        "ix_global": int(row["ix_global"]),
                        "direction": int(row["direction"]),
                        "face_index": int(row["face_index"]),
                        "transverse_index": int(row["transverse_index"]),
                        "operator_stage": int(row["operator_stage"]),
                        "face_number": face,
                        "parcel_number": parcel,
                        "absolute_error": absolute_error,
                        "relative_error": relative_error,
                        "node_failure_reason": row.get(
                            "node_failure_reason", math.nan),
                        "reconstructed_target": row.get(
                            "reconstructed_target", math.nan),
                        "node_sum": row.get("node_sum", math.nan),
                    }
        parcel_face_relative = (face_delta / face_scale
                                if independent_face_ledger and face_scale > 0.0
                                else (0.0 if independent_face_ledger else math.nan))
        below_relative = (below_threshold / max(1.0, face_scale)
                          if independent_face_ledger and math.isfinite(below_threshold)
                          else math.nan)
        static_wall = total_wall_seconds(args.static)
        audit_wall = total_wall_seconds(args.audit)
        static_timings = timing_sums(args.static)
        audit_timings = timing_sums(args.audit)
    wall_relative = (audit_wall / static_wall - 1.0
                     if math.isfinite(static_wall) and static_wall > 0.0 and
                     math.isfinite(audit_wall) else math.nan)
    wall_ok = args.ignore_performance or (
        math.isfinite(wall_relative) and wall_relative <= 0.15)

    # A negative quadrature node can be recorded by a legacy/stale audit
    # binary even when its absolute mass is many orders below the complete
    # interface ledger.  Keep the raw condition visible, but do not reject
    # physical flux closure for a zero-scale tail residue.  Structural
    # failures (non-finite nodes, non-positive support, duplicate faces) are
    # never downgraded here.
    benign_negative_node_tolerance = 1.0e-12 * max(1.0, face_scale)
    benign_negative_node_count = 0
    raw_audit_invalid_rows = []
    if audit_valid_available:
        for row in rows or []:
            if (row.get("audit_valid", math.nan) == 1.0 and
                    row.get("audit_failure_code", math.nan) == 0.0):
                continue
            raw_audit_invalid_rows.append(row)
            if (row.get("audit_failure_code", math.nan) == 1.0 and
                    row.get("audit_parcel_failure_reason", math.nan) == 1.0 and
                    math.isfinite(row.get("audit_parcel_failure_node_mass", math.nan)) and
                    abs(row["audit_parcel_failure_node_mass"]) <=
                    benign_negative_node_tolerance):
                benign_negative_node_count += 1
    audit_valid_effective = audit_valid_ok or (
        bool(raw_audit_invalid_rows) and
        benign_negative_node_count == len(raw_audit_invalid_rows))
    if comparison_valid:
        face_keys = []
        if independent_face_ledger and not compact_face_summary:
            for row in face_rows or []:
                face_keys.append((int(row["accepted_step"]),
                                  int(row["ix_global"]),
                                  int(row["direction"]),
                                  int(row["face_index"]),
                                  int(row["transverse_index"]),
                                  int(row["operator_stage"])))
        face_duplicate_count = (len(face_keys) - len(set(face_keys))
                                if not compact_face_summary else 0)
    audit_ok = (
        comparison_valid and rows_present and inplace_state_ok and
        inplace_rng_ok and inplace_ledger_ok and audit_valid_effective and
        independent_face_ledger and face_duplicate_count == 0 and
        duplicate_count == 0.0 and below_relative <= 1.0e-12 and
        parcel_face_relative <= 1.0e-12 and
        parcel_face_max_relative <= 1.0e-13 and
        wall_ok
    )
    if not comparison_valid:
        flux_ledger_status = "not_evaluated_run_interval_invalid"
    elif not independent_face_ledger:
        flux_ledger_status = "no_interface_flux"
    elif compact_face_summary:
        flux_ledger_status = "evaluated_compact_summary"
    else:
        flux_ledger_status = "evaluated"
    with open(args.result, "w") as out:
        out.write("static_accepted_step_count=%s\n" % static_count)
        out.write("audit_accepted_step_count=%s\n" % audit_count)
        out.write("accepted_step_count_equal=%d\n" % int(accepted_equal))
        for label, info in (("static", static_info), ("audit", audit_info)):
            out.write("%s_first_accepted_step=%s\n" %
                      (label, None if info is None else info["first_step"]))
            out.write("%s_first_accepted_time_s=%.17g\n" %
                      (label, math.nan if info is None else info["first_time_s"]))
            out.write("%s_last_accepted_step=%s\n" %
                      (label, None if info is None else info["last_step"]))
            out.write("%s_last_accepted_time_s=%.17g\n" %
                      (label, math.nan if info is None else info["last_time_s"]))
        out.write("static_vpfp_failure=%d\n" % int(static_failed))
        out.write("audit_vpfp_failure=%d\n" % int(audit_failed))
        out.write("comparison_valid=%d\n" % int(comparison_valid))
        out.write("comparison_failure_reason=%s\n" %
                  ("none" if comparison_valid else ",".join(comparison_reasons)))
        out.write("cross_run_metrics_evaluated=%d\n" % int(comparison_valid))
        out.write("state_hash_common_file_count=%d\n" % common_files)
        out.write("state_hash_equal=%d\n" % int(equal))
        out.write("state_hash_is_informational=1\n")
        out.write("audit_inplace_state_available=%d\n" %
                  int(inplace_state_available))
        out.write("audit_inplace_rng_available=%d\n" %
                  int(inplace_rng_available))
        out.write("audit_inplace_ledger_available=%d\n" %
                  int(inplace_ledger_available))
        out.write("audit_inplace_state_bitwise_equal=%d\n" %
                  int(inplace_state_ok))
        out.write("audit_inplace_rng_equal=%d\n" % int(inplace_rng_ok))
        out.write("audit_inplace_ledger_equal=%d\n" % int(inplace_ledger_ok))
        out.write("audit_valid_available=%d\n" % int(audit_valid_available))
        out.write("audit_valid_all_accepted_steps=%d\n" % int(audit_valid_ok))
        out.write("audit_valid_effective=%d\n" % int(audit_valid_effective))
        out.write("benign_negative_node_tolerance=%.17g\n" %
                  benign_negative_node_tolerance)
        out.write("benign_negative_node_count=%d\n" % benign_negative_node_count)
        out.write("audit_flux_rows=%d\n" % int(rows_present))
        out.write("audit_face_rows=%d\n" % int(face_rows_present))
        out.write("audit_face_row_count=%d\n" %
                  (0 if face_rows is None else len(face_rows)))
        out.write("compact_face_summary=%d\n" % int(compact_face_summary))
        out.write("audit_face_duplicate_count=%d\n" % face_duplicate_count)
        out.write("interface_duplicate_count=%.17g\n" % duplicate_count)
        out.write("below_threshold_number=%.17g\n" % below_threshold)
        out.write("roundoff_discarded_number=%.17g\n" % roundoff_discarded)
        out.write("independent_face_ledger=%d\n" % int(independent_face_ledger))
        out.write("flux_ledger_status=%s\n" % flux_ledger_status)
        out.write("face_flux_absolute_scale=%.17g\n" % face_scale)
        out.write("parcel_flux_absolute_scale=%.17g\n" % parcel_scale)
        out.write("parcel_vs_face_flux_relative_l1=%.17g\n" % parcel_face_relative)
        out.write("parcel_vs_face_flux_max_relative=%.17g\n" %
                  parcel_face_max_relative)
        out.write("parcel_vs_face_flux_absolute_at_max_relative=%.17g\n" %
                  parcel_face_absolute_at_max_relative)
        if max_error_row is None:
            out.write("parcel_vs_face_flux_max_location=not_evaluated\n")
        else:
            out.write(
                "parcel_vs_face_flux_max_location="
                "step=%d,time_s=%.17g,ix_global=%d,direction=%d,"
                "face_index=%d,transverse_index=%d,operator_stage=%d,"
                "face_number=%.17g,parcel_number=%.17g,"
                "node_failure_reason=%s,reconstructed_target=%.17g,"
                "node_sum=%.17g\n" % (
                    max_error_row["accepted_step"], max_error_row["time_s"],
                    max_error_row["ix_global"], max_error_row["direction"],
                    max_error_row["face_index"],
                    max_error_row["transverse_index"],
                    max_error_row["operator_stage"],
                    max_error_row["face_number"],
                    max_error_row["parcel_number"],
                    max_error_row["node_failure_reason"],
                    max_error_row["reconstructed_target"],
                    max_error_row["node_sum"]))
        out.write("global_parcel_vs_interface_flux_relative=%.17g\n" % parcel_face_relative)
        out.write("below_threshold_number_relative=%.17g\n" % below_relative)
        out.write("static_wall_seconds=%.17g\n" % static_wall)
        out.write("audit_wall_seconds=%.17g\n" % audit_wall)
        out.write("audit_wall_overhead_relative=%.17g\n" % wall_relative)
        out.write("audit_wall_overhead_pass=%d\n" % int(wall_ok))
        out.write("performance_gate_ignored=%d\n" % int(args.ignore_performance))
        for name in TIMING_FIELDS:
            static_value = static_timings[name]
            audit_value = audit_timings[name]
            relative = (audit_value / static_value - 1.0
                        if math.isfinite(static_value) and static_value > 0.0 and
                        math.isfinite(audit_value) else math.nan)
            out.write("static_%s=%.17g\n" % (name, static_value))
            out.write("audit_%s=%.17g\n" % (name, audit_value))
            out.write("audit_%s_overhead_relative=%.17g\n" %
                      (name, relative))
        if static_failed:
            out.write("static_vpfp_failure_record=%s\n" % static_failure)
        if audit_failed:
            out.write("audit_vpfp_failure_record=%s\n" % audit_failure)
        out.write("status=%s\n" % ("PASS" if audit_ok else "FAIL"))
    return 0 if audit_ok else 2


if __name__ == "__main__":
    sys.exit(main())
