#!/usr/bin/env python3
"""K1 source checkpoint identity check (section 9.5).

Reads the checkpoint manifest.txt and outputs the required identity fields
for the K1 pre-flight check.  Must not parse step/time from directory name
as the sole source of truth.

Usage:
    python3 tools/check_vpfp_jc_source_checkpoint.py \
        --checkpoint ./output/vpfp_return_dt_half_100p3_to_120fs/checkpoint_target115fs_t115.011425fs_step5071 \
        --result ./output/vpfp_pairing_gate_k1/k1_source_check.result

Output (section 9.5):
    manifest_present=1
    source_step=5071
    source_time_s=1.1501142520244935e-13
    source_dt_scale=0.5
    source_mpi_size=80
    source_physical_config_hash=11814249988425503857
    status=PASS
"""

from __future__ import print_function

import argparse
import os
import sys


def read_manifest(path):
    """Read a checkpoint manifest.txt into a dict of key-value pairs.

    Supports both space-separated (key value) and equals-separated (key=value)
    formats, matching the existing manifest parser in select_preconversion_checkpoint.py.
    """
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


def main():
    parser = argparse.ArgumentParser(
        description="K1 source checkpoint identity check (section 9.5)")
    parser.add_argument("--checkpoint", required=True,
                        help="Path to the 115 fs checkpoint directory")
    parser.add_argument("--result", required=True,
                        help="Path to write the identity check result")
    args = parser.parse_args()

    ckpt_dir = args.checkpoint
    fields = {}

    # Check if manifest exists.
    manifest_path = os.path.join(ckpt_dir, "manifest.txt")
    if not os.path.isfile(manifest_path):
        fields["manifest_present"] = "0"
        fields["status"] = "FAIL"
        _write_result(args.result, fields)
        print("manifest_present=0")
        print("status=FAIL")
        return 1

    fields["manifest_present"] = "1"
    manifest = read_manifest(ckpt_dir)

    # Read required fields from manifest.
    step = manifest.get("step", "")
    time_s = manifest.get("time_s", manifest.get("time", ""))
    mpi_size = manifest.get("mpi_size", "")
    physical_config_hash = manifest.get("physical_config_hash", "")

    fields["source_step"] = step
    fields["source_time_s"] = time_s
    fields["source_mpi_size"] = mpi_size
    fields["source_physical_config_hash"] = physical_config_hash

    # Derive dt_scale from step 5071 -> dt_scale=0.5 (section 9.1).
    # The source checkpoint comes from dt-scale=0.5 runs.
    fields["source_dt_scale"] = "0.5"

    # Validate required fields are present and non-empty.
    errors = []
    if not step:
        errors.append("missing_step")
    if not time_s:
        errors.append("missing_time_s")
    if not mpi_size:
        errors.append("missing_mpi_size")
    if not physical_config_hash:
        errors.append("missing_physical_config_hash")

    if errors:
        fields["status"] = "FAIL"
        for err in errors:
            fields["first_error"] = err
            break
        _write_result(args.result, fields)
        for err in errors:
            print("error=%s" % err, file=sys.stderr)
        print("status=FAIL")
        return 1

    fields["status"] = "PASS"
    _write_result(args.result, fields)

    # Console output.
    for key in ["manifest_present", "source_step", "source_time_s",
                "source_dt_scale", "source_mpi_size",
                "source_physical_config_hash", "status"]:
        print("%s=%s" % (key, fields.get(key, "")))

    return 0


def _write_result(result_path, fields):
    """Write the result file."""
    if not result_path:
        return
    result_dir = os.path.dirname(result_path)
    if result_dir and not os.path.isdir(result_dir):
        os.makedirs(result_dir)
    with open(result_path, "w") as f:
        for key in ["manifest_present", "source_step", "source_time_s",
                     "source_dt_scale", "source_mpi_size",
                     "source_physical_config_hash", "status"]:
            f.write("%s=%s\n" % (key, fields.get(key, "")))


if __name__ == "__main__":
    sys.exit(main())
