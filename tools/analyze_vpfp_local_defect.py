#!/usr/bin/env python3
"""Summarise the Gate-F fixed-accepted-state local-defect experiment.

The program validates comparability first.  It deliberately does not promote
or reject a run on a preselected convergence order: the reported q values are
evidence for the subsequent root-cause classification.
"""

from __future__ import print_function

import argparse
import math
import os
import sys


REQUIRED_CONTRACT = (
    "checkpoint_derived_state_reconstructed",
    "same_initial_state", "same_physical_window", "collision_frozen",
    "conversion_frozen", "h10_frozen",
    "accepted_state_bitwise_equal_after_audit",
    "rng_bitwise_equal_after_audit", "ledger_bitwise_equal_after_audit",
    "source_ownership_valid", "poisson_identity_pass",
    "electrode_work_ownership_count",
    "background_boundary_energy_ownership_count",
    "beam_boundary_energy_ownership_count",
    "velocity_boundary_energy_ownership_count",
    "conversion_energy_ownership_count",
    "tail_return_energy_ownership_count",
    "collision_reservoir_ownership_count",
)
REQUIRED_INTERVALS = (1, 2, 5, 10)


def read_key_values(path):
    values = {}
    with open(path, "r") as stream:
        for line in stream:
            line = line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key] = value
    return values


def read_rows(path):
    with open(path, "r") as stream:
        header = stream.readline().split()
        if not header:
            raise ValueError("missing interval header")
        rows = []
        for line_number, line in enumerate(stream, 2):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != len(header):
                raise ValueError("line %d has invalid column count" % line_number)
            row = dict(zip(header, fields))
            for key, value in row.items():
                if key in ("dt_scale", "interval"):
                    continue
                try:
                    numeric = float(value)
                except ValueError:
                    raise ValueError("line %d invalid %s" % (line_number, key))
                if not math.isfinite(numeric):
                    raise ValueError("line %d non-finite %s" % (line_number, key))
            rows.append(row)
    return rows


def number(row, key):
    return float(row[key])


def order(a, b):
    if a <= 0.0 or b <= 0.0:
        return float("nan")
    return math.log(a / b, 2.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    parser.add_argument("--result", required=True)
    args = parser.parse_args()

    output = []
    try:
        contract = read_key_values(os.path.join(args.run, "audit_contract.result"))
        rows = read_rows(os.path.join(args.run, "local_defect_intervals.dat"))
        contract_ok = all(contract.get(key) == "1" for key in REQUIRED_CONTRACT)
        ownership_residual = float(contract.get("source_ownership_residual", "inf"))
        contract_ok = contract_ok and math.isfinite(ownership_residual)
        by_interval = {}
        for row in rows:
            scale = number(row, "dt_scale")
            interval = int(number(row, "interval"))
            by_interval.setdefault(interval, {})[scale] = row
        expected_scales = (1.0, 0.5, 0.25)
        comparable = contract_ok
        if tuple(sorted(by_interval)) != REQUIRED_INTERVALS:
            comparable = False
        for interval, levels in by_interval.items():
            if any(scale not in levels for scale in expected_scales):
                comparable = False
            for row in levels.values():
                if int(number(row, "accepted")) != 1 or \
                   int(number(row, "finite")) != 1 or \
                   int(number(row, "poisson_identity_pass")) != 1:
                    comparable = False

        q_values = []
        stable_abs_orders = []
        for interval in sorted(by_interval):
            levels = by_interval[interval]
            if not all(scale in levels for scale in expected_scales):
                continue
            signed_10 = abs(number(levels[1.0], "field_particle_pair_residual_signed"))
            signed_05 = abs(number(levels[0.5], "field_particle_pair_residual_signed"))
            signed_025 = abs(number(levels[0.25], "field_particle_pair_residual_signed"))
            abs_10 = number(levels[1.0], "field_particle_pair_residual_abs")
            abs_05 = number(levels[0.5], "field_particle_pair_residual_abs")
            abs_025 = number(levels[0.25], "field_particle_pair_residual_abs")
            q_signed = order(signed_10, signed_05)
            q_abs = order(abs_10, abs_05)
            q_signed_fine = order(signed_05, signed_025)
            q_abs_fine = order(abs_05, abs_025)
            q_values.append((q_signed, q_abs, q_signed_fine, q_abs_fine))
            # A negative order is normally a cancellation/non-asymptotic
            # prefix, not evidence of a fixed spatial incompatibility.  Use an
            # interval for root classification only when the absolute defect
            # decreases on both refinements with a reasonably stable order.
            if (math.isfinite(q_abs) and math.isfinite(q_abs_fine) and
                    q_abs > 0.0 and q_abs_fine > 0.0 and
                    abs(q_abs - q_abs_fine) <= 0.75):
                stable_abs_orders.extend((q_abs, q_abs_fine))
            output.extend([
                "q_signed_interval_%d=%.17g" % (interval, q_signed),
                "q_abs_interval_%d=%.17g" % (interval, q_abs),
                "q_signed_fine_interval_%d=%.17g" %
                    (interval, q_signed_fine),
                "q_abs_fine_interval_%d=%.17g" %
                    (interval, q_abs_fine),
            ])

        poisson_ok = contract.get("poisson_identity_pass") == "1"
        finite_q = [q for group in q_values for q in group if math.isfinite(q)]
        if contract.get("source_ownership_valid") != "1":
            classification = "source-ownership"
        elif not poisson_ok:
            classification = "poisson-space"
        elif not finite_q:
            classification = "fixed-coupling"
        elif stable_abs_orders and (
                sorted(stable_abs_orders)[len(stable_abs_orders) // 2] >= 1.0):
            classification = "time-truncation-mixed-order"
        else:
            classification = "mixed-order-unresolved"
        output.extend([
            "checkpoint_derived_state_reconstructed=%s" %
                contract.get("checkpoint_derived_state_reconstructed", "0"),
            "same_initial_state=%s" % contract.get("same_initial_state", "0"),
            "same_physical_window=%s" % contract.get("same_physical_window", "0"),
            "collision_frozen=%s" % contract.get("collision_frozen", "0"),
            "conversion_frozen=%s" % contract.get("conversion_frozen", "0"),
            "h10_frozen=%s" % contract.get("h10_frozen", "0"),
            "accepted_state_bitwise_equal_after_audit=%s" %
                contract.get("accepted_state_bitwise_equal_after_audit", "0"),
            "rng_bitwise_equal_after_audit=%s" %
                contract.get("rng_bitwise_equal_after_audit", "0"),
            "ledger_bitwise_equal_after_audit=%s" %
                contract.get("ledger_bitwise_equal_after_audit", "0"),
            "source_ownership_valid=%s" %
                contract.get("source_ownership_valid", "0"),
            "poisson_identity_pass=%d" % (1 if poisson_ok else 0),
            "classification=%s" % classification,
            "status=%s" % ("PASS" if comparable else "FAIL"),
        ])
    except (IOError, OSError, ValueError, KeyError) as error:
        output.extend([
            "source_ownership_valid=0",
            "poisson_identity_pass=0",
            "classification=stochastic-followup",
            "error=%s" % str(error).replace("\n", " "),
            "status=FAIL",
        ])
    with open(args.result, "w") as stream:
        stream.write("\n".join(output) + "\n")
    return 0 if output[-1] == "status=PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
