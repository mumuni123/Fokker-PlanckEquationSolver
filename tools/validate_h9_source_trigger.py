#!/usr/bin/env python3
"""Validate the three production-path gates required by H9 section 17.10.2.6."""

from __future__ import print_function

import argparse
import math
import os
import sys


def read_result(path):
    values = {}
    if not os.path.isfile(path):
        return values
    with open(path, "r") as stream:
        for line in stream:
            if "=" not in line:
                continue
            key, value = line.strip().split("=", 1)
            values[key] = value
    return values


def flag(values, name):
    return values.get(name) == "1"


def finite_at_most(values, name, limit):
    try:
        return math.isfinite(float(values[name])) and float(values[name]) <= limit
    except (KeyError, ValueError):
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True)
    parser.add_argument("--collision", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--result", required=True)
    args = parser.parse_args()

    source = read_result(args.source)
    collision = read_result(args.collision)
    checkpoint = read_result(args.checkpoint)
    gates = {
        "source_status": source.get("status") == "PASS",
        "source_face": flag(source, "source_face_found"),
        "source_flux": flag(source, "source_flux_seen"),
        "source_bulk_sink": flag(source, "bulk_loss_matches_face"),
        "source_tail_created": flag(source, "tail_created"),
        "source_conversion": flag(source, "conversion_complete"),
        "source_six_moment": flag(source, "six_moment_closed"),
        "source_tail_number": flag(source, "tail_number_matches_created"),
        "collision_status": collision.get("status") == "PASS",
        "collision_pairs": flag(collision, "collision_pairs_ok"),
        "collision_reaction": finite_at_most(
            collision, "collision_reaction_residual_max", 1.0e-10),
        "checkpoint_status": checkpoint.get("status") == "PASS",
        "checkpoint_roundtrip": flag(checkpoint, "checkpoint_roundtrip"),
        "checkpoint_metadata": flag(checkpoint, "interface_metadata_roundtrip"),
    }
    passed = all(gates.values())
    with open(args.result, "w") as out:
        for name in sorted(gates):
            out.write("gate_%s=%d\n" % (name, 1 if gates[name] else 0))
        out.write("status=%s\n" % ("PASS" if passed else "FAIL"))
    return 0 if passed else 2


if __name__ == "__main__":
    sys.exit(main())
