#!/usr/bin/env python3
"""Select the latest accepted checkpoint before flux conversion starts."""

from __future__ import print_function

import argparse
import os
import re
import sys


def manifest_values(path):
    values = {}
    manifest = os.path.join(path, "manifest.txt")
    if not os.path.isfile(manifest):
        return values
    with open(manifest, "r") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 2:
                values[fields[0].strip()] = fields[1].strip()
            elif fields and "=" in fields[0]:
                key, value = fields[0].split("=", 1)
                values[key.strip()] = value.strip()
    return values


def checkpoint_candidates(root):
    candidates = []
    for base, dirs, _ in os.walk(root):
        for name in dirs:
            if not name.startswith("checkpoint_"):
                continue
            path = os.path.join(base, name)
            values = manifest_values(path)
            if not values:
                continue
            try:
                time = float(values.get("time_s", values.get("time", "nan")))
            except ValueError:
                time = float("nan")
            if time != time:
                match = re.search(r"([0-9]+(?:\.[0-9]+)?)", name)
                time = float(match.group(1)) if match else float("nan")
            candidates.append((time, path, values))
    candidates.sort(key=lambda item: item[0])
    return candidates


def conversion_number_status(values):
    keys = (
        "conversion_cumulative_number",
        "conversion_number",
        "conversion_number_removed",
    )
    found = False
    for key in keys:
        if key in values:
            found = True
            try:
                value = float(values[key])
            except ValueError:
                return False, "invalid_%s=%s" % (key, values[key])
            if abs(value) > 0.0:
                return False, "nonzero_%s=%.17g" % (key, value)
    # A legacy checkpoint without an accepted conversion ledger is not proof
    # that the run was pre-conversion.  The 17E selector must not silently
    # reinterpret such a checkpoint.
    if not found:
        return False, "missing_conversion_ledger"
    return True, "zero_conversion_ledger"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--require-conversion-number-zero", action="store_true")
    parser.add_argument("--result", required=True)
    args = parser.parse_args()

    candidates = checkpoint_candidates(args.input)
    selected = None
    candidate_status = []
    for time, path, values in candidates:
        valid, reason = conversion_number_status(values)
        candidate_status.append((time, path, reason))
        if args.require_conversion_number_zero and not valid:
            continue
        selected = path
    ok = selected is not None
    with open(args.result, "w") as out:
        out.write("candidate_count=%d\n" % len(candidates))
        for index, (time, path, reason) in enumerate(candidate_status):
            out.write("candidate_%d_time_s=%.17g\n" % (index, time))
            out.write("candidate_%d_path=%s\n" % (index, path))
            out.write("candidate_%d_conversion_status=%s\n" %
                      (index, reason))
        out.write("selected_checkpoint=%s\n" % (selected or ""))
        out.write("require_conversion_number_zero=%d\n" %
                  int(args.require_conversion_number_zero))
        out.write("status=%s\n" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
