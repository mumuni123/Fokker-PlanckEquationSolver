#!/usr/bin/env python3
"""Run the read-only §15.14.9 point-4 audit.

This test intentionally does not run the coarse/fine production short window.
It only re-runs the corrected dual-face analyzer on existing output and runs
the manufactured G/G* tests.  The resulting point4.result can be copied back
with the audit output.
"""

from __future__ import print_function

import argparse
import os
import shlex
import subprocess
import sys


def read_result(path):
    values = {}
    with open(path, "r") as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def as_int(values, key, default=0):
    try:
        return int(values.get(key, default))
    except (TypeError, ValueError):
        return default


def as_float(values, key, default=float("nan")):
    try:
        return float(values.get(key, default))
    except (TypeError, ValueError):
        return default


def run_command(command):
    print("+ " + " ".join(shlex.quote(str(x)) for x in command))
    try:
        return subprocess.call(command)
    except OSError as error:
        print("command_error=%s" % error)
        return 127


def require(values, key, expected, failures, label):
    actual = values.get(key)
    if actual != str(expected):
        failures.append("%s: %s=%s, expected %s" %
                        (label, key, actual, expected))


def check_manufactured(path, label, failures):
    if not os.path.isfile(path):
        failures.append("%s result missing: %s" % (label, path))
        return
    values = read_result(path)
    require(values, "status", "PASS", failures, label)
    require(values, "manufactured_identity_pass", 1, failures, label)
    require(values, "manufactured_endpoint_weight_pass", 1, failures, label)


def parse_args():
    parser = argparse.ArgumentParser(
        description="§15.14.9 point-4 analyzer/manufactured-test audit")
    parser.add_argument("--coarse", required=True,
                        help="existing coarse production output directory")
    parser.add_argument("--fine", required=True,
                        help="existing fine production output directory")
    parser.add_argument("--dual-result", required=True,
                        help="new dual-face analyzer result path")
    parser.add_argument("--single-exe", required=True,
                        help="single-rank manufactured pairing test executable")
    parser.add_argument("--single-result", required=True,
                        help="single-rank manufactured test result path")
    parser.add_argument("--mpi-exe", default="",
                        help="MPI manufactured pairing test executable")
    parser.add_argument("--mpi-result", default="",
                        help="MPI manufactured test result path")
    parser.add_argument("--mpi-launch", default="yhrun",
                        help="MPI launcher, e.g. yhrun or mpirun")
    parser.add_argument("--mpi-ranks", type=int, default=2,
                        help="MPI ranks for the manufactured test")
    parser.add_argument("--result",
                        default="output/vpfp_pairing_gate_k1/point4.result",
                        help="point4 summary result path")
    parser.add_argument("--analyzer", default="",
                        help="corrected dual-face analyzer path")
    return parser.parse_args()


def main():
    args = parse_args()
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    analyzer = args.analyzer or os.path.join(
        root, "tools", "analyze_vpfp_k1_dual_face.py")
    failures = []

    # Point 4: analyze existing production output only.  No solver executable
    # and no production step is called here.
    rc = run_command([
        sys.executable, analyzer,
        "--coarse", args.coarse,
        "--fine", args.fine,
        "--result", args.dual_result,
    ])
    if rc != 0:
        failures.append("dual-face analyzer exit code %d" % rc)

    if not os.path.isfile(args.dual_result):
        failures.append("dual result missing: %s" % args.dual_result)
    else:
        dual = read_result(args.dual_result)
        require(dual, "status", "PASS", failures, "dual-face")
        for prefix in ("coarse", "fine"):
            require(dual, prefix + "_dual_recon_pass", 1,
                    failures, "dual-face")
            require(dual, prefix + "_missing_owner_face_count_max", 0,
                    failures, "dual-face")
            require(dual, prefix + "_owner_duplicate_face_count_max", 0,
                    failures, "dual-face")
            require(dual, prefix + "_owner_rank_error_count_max", 0,
                    failures, "dual-face")
            require(dual, prefix + "_nonfinite_face_count_max", 0,
                    failures, "dual-face")
            for key in ("max_dual_in_domain_relative",
                        "max_current_pair_reconstruction_relative"):
                value = as_float(dual, prefix + "_" + key)
                if not value == value or value > 1.0e-10:
                    failures.append("dual-face: %s_%s=%s exceeds 1e-10" %
                                    (prefix, key, value))

    # Point 4 also re-runs the manufactured tests, never the production window.
    if not os.path.isfile(args.single_exe):
        failures.append("single-rank executable missing: %s" % args.single_exe)
    elif run_command([args.single_exe, "--case", "all",
                      "--result", args.single_result]) != 0:
        failures.append("single-rank manufactured test returned nonzero")
    check_manufactured(args.single_result, "single-rank manufactured", failures)

    if args.mpi_exe:
        if not args.mpi_result:
            failures.append("--mpi-result is required with --mpi-exe")
        elif not os.path.isfile(args.mpi_exe):
            failures.append("MPI executable missing: %s" % args.mpi_exe)
        else:
            launcher = shlex.split(args.mpi_launch)
            command = launcher + ["-n", str(args.mpi_ranks), args.mpi_exe,
                                  "--case", "all", "--result", args.mpi_result]
            if run_command(command) != 0:
                failures.append("MPI manufactured test returned nonzero")
            check_manufactured(args.mpi_result, "MPI manufactured", failures)

    os.makedirs(os.path.dirname(os.path.abspath(args.result)), exist_ok=True)
    with open(args.result, "w") as stream:
        stream.write("status=%s\n" % ("PASS" if not failures else "FAIL"))
        stream.write("production_short_window_rerun=0\n")
        stream.write("dual_face_analyzer_rerun=1\n")
        stream.write("manufactured_single_rank_rerun=1\n")
        stream.write("manufactured_mpi_rerun=%d\n" %
                     (1 if args.mpi_exe else 0))
        stream.write("failure_count=%d\n" % len(failures))
        for index, failure in enumerate(failures):
            stream.write("failure_%d=%s\n" % (index, failure))

    print("point4_result=" + os.path.abspath(args.result))
    for failure in failures:
        print("FAIL: " + failure)
    print("status=" + ("PASS" if not failures else "FAIL"))
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
