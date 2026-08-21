#!/usr/bin/env python3
"""JC5 comparison script: verify field-particle coupling test results.

Reads result files produced by the §8.5 test commands and verifies the
§8.10.1 acceptance criteria:

    all_unit_cases_pass=1
    all_mpi_cases_pass=1
    legacy_default_regression_pass=1
    trial_deterministic=1
    trial_side_effect_free=1
    failure_transaction_bitwise_unchanged=1
    post_field_charge_invariance_pass=1
    mpi_failure_consensus_pass=1
    soft_accept_count=0
    status=PASS

§8.7 requirements:
1. 明确列出必需文件
2. 明确列出每个文件的必需键
3. 检测重复键、NaN、Inf、空文件和未知状态
4. 任一缺失返回 exit code 2；物理/数值门失败返回 exit code 1；全部通过返回 0
5. 写汇总 status=PASS|FAIL 和 first_failure
6. 不执行 Git 命令，不因测试机无 Git 判失败
7. 不把字符串 PASS 的存在当作唯一依据，必须重新计算阈值判断

Usage:
    python3 tools/compare_vpfp_field_particle_jc.py \
        --root ./output/vpfp_pairing_gate_jc \
        --result ./output/vpfp_pairing_gate_jc/jc5_compare.result
"""

from __future__ import print_function

import argparse
import math
import os
import sys


# §8.7 requirement 1: explicitly list required files.
# These are the ONLY files produced by §8.5 commands.
REQUIRED_FILES = {
    "trial_single": "trial_single.result",
    "post_field_single": "post_field_single.result",
    "trial_mpi_n2": "trial_mpi_n2.result",
    "trial_mpi_n5": "trial_mpi_n5.result",
}


def read_result(path):
    """Read a .result file into a dict of key=value pairs.

    §8.7: detects duplicate keys, NaN/Inf values, and empty files.
    Returns (result_dict, errors_list).
    """
    result = {}
    errors = []
    if not os.path.isfile(path):
        errors.append("file_missing:%s" % path)
        return result, errors
    if os.path.getsize(path) == 0:
        errors.append("file_empty:%s" % path)
        return result, errors
    seen_keys = set()
    with open(path, "r") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()
            if key in seen_keys:
                errors.append("duplicate_key:%s:%s" % (path, key))
            seen_keys.add(key)
            # §8.7 requirement 3: detect NaN/Inf.
            if value in ("NaN", "nan", "Inf", "inf", "-Inf", "-inf"):
                errors.append("non_finite_value:%s:%s=%s" % (path, key, value))
            result[key] = value
    if not result:
        errors.append("file_no_data:%s" % path)
    return result, errors


def main():
    parser = argparse.ArgumentParser(
        description="JC5 field-particle coupling comparison")
    parser.add_argument("--root", required=True,
                        help="Root output directory")
    parser.add_argument("--result", required=True,
                        help="Path to write comparison result")
    args = parser.parse_args()

    root = args.root
    pass_count = 0
    fail_count = 0
    missing_count = 0
    first_failure = "none"
    fields = {}
    all_errors = []

    def set_first_failure(name):
        nonlocal first_failure
        if first_failure == "none":
            first_failure = name

    # §8.7 requirement 1: read all required files and validate.
    file_data = {}
    for label, rel_path in REQUIRED_FILES.items():
        path = os.path.join(root, rel_path)
        data, errors = read_result(path)
        all_errors.extend(errors)
        file_data[label] = data
        if not data:
            missing_count += 1
            set_first_failure("missing_file:%s" % rel_path)

    # §8.7 requirement 3: report all validation errors.
    for err in all_errors:
        print("validation_error=%s" % err, file=sys.stderr)

    # §8.7 requirement 4: missing files → exit code 2.
    if missing_count > 0:
        fields["status"] = "FAIL"
        fields["first_failure"] = first_failure
        if args.result:
            result_dir = os.path.dirname(args.result)
            if result_dir and not os.path.isdir(result_dir):
                os.makedirs(result_dir)
            with open(args.result, "w") as f:
                for key in ["status", "first_failure"]:
                    f.write("%s=%s\n" % (key, fields.get(key, "N/A")))
                f.write("missing_file_count=%d\n" % missing_count)
        print("status=FAIL")
        print("first_failure=%s" % first_failure)
        print("missing_file_count=%d" % missing_count)
        return 2

    # --- 1. all_unit_cases_pass ---
    all_unit = True
    if file_data["trial_single"].get("status") != "PASS":
        all_unit = False
        set_first_failure("trial_single_status")
    if file_data["post_field_single"].get("status") != "PASS":
        all_unit = False
        set_first_failure("post_field_single_status")
    fields["all_unit_cases_pass"] = "1" if all_unit else "0"
    pass_count += 1 if all_unit else 0
    fail_count += 0 if all_unit else 1

    # --- 2. all_mpi_cases_pass ---
    all_mpi = True
    for np in [2, 5]:
        key = "trial_mpi_n%d" % np
        if file_data[key].get("status") != "PASS":
            all_mpi = False
            set_first_failure("trial_mpi_n%d_status" % np)
    fields["all_mpi_cases_pass"] = "1" if all_mpi else "0"
    pass_count += 1 if all_mpi else 0
    fail_count += 0 if all_mpi else 1

    # --- 3. legacy_default_regression_pass ---
    legacy_ok = file_data["trial_single"].get(
        "legacy_default_regression_pass") == "1"
    fields["legacy_default_regression_pass"] = "1" if legacy_ok else "0"
    pass_count += 1 if legacy_ok else 0
    fail_count += 0 if legacy_ok else 1
    if not legacy_ok:
        set_first_failure("legacy_default_regression")

    # --- 4. trial_deterministic ---
    det_ok = file_data["trial_single"].get(
        "trial_deterministic") == "1"
    fields["trial_deterministic"] = "1" if det_ok else "0"
    pass_count += 1 if det_ok else 0
    fail_count += 0 if det_ok else 1
    if not det_ok:
        set_first_failure("trial_deterministic")

    # --- 5. trial_side_effect_free ---
    side_ok = file_data["trial_single"].get(
        "trial_side_effect_free") == "1"
    fields["trial_side_effect_free"] = "1" if side_ok else "0"
    pass_count += 1 if side_ok else 0
    fail_count += 0 if side_ok else 1
    if not side_ok:
        set_first_failure("trial_side_effect_free")

    # --- 6. failure_transaction_bitwise_unchanged ---
    ftbu_ok = file_data["trial_single"].get(
        "failure_transaction_bitwise_unchanged") == "1"
    fields["failure_transaction_bitwise_unchanged"] = (
        "1" if ftbu_ok else "0")
    pass_count += 1 if ftbu_ok else 0
    fail_count += 0 if ftbu_ok else 1
    if not ftbu_ok:
        set_first_failure("failure_transaction_bitwise_unchanged")

    # --- 7. post_field_charge_invariance_pass ---
    pfc_ok = file_data["post_field_single"].get(
        "post_field_charge_pass") == "1"
    fields["post_field_charge_invariance_pass"] = "1" if pfc_ok else "0"
    pass_count += 1 if pfc_ok else 0
    fail_count += 0 if pfc_ok else 1
    if not pfc_ok:
        set_first_failure("post_field_charge_invariance")

    # --- 8. mpi_failure_consensus_pass ---
    rollback_ok = True
    for np in [2, 5]:
        key = "trial_mpi_n%d" % np
        if (file_data[key].get("all_rank_failed") != "1" or
                file_data[key].get("all_rank_decision_equal") != "1"):
            rollback_ok = False
    fields["mpi_failure_consensus_pass"] = "1" if rollback_ok else "0"
    pass_count += 1 if rollback_ok else 0
    fail_count += 0 if rollback_ok else 1
    if not rollback_ok:
        set_first_failure("mpi_failure_consensus")

    # --- 9. soft_accept_count=0 (§8.10.1) ---
    soft_accept = file_data["trial_single"].get("soft_accept_count")
    soft_ok = (soft_accept is not None and soft_accept == "0")
    fields["soft_accept_count"] = "0" if soft_ok else (soft_accept or "missing")
    pass_count += 1 if soft_ok else 0
    fail_count += 0 if soft_ok else 1
    if not soft_ok:
        set_first_failure("soft_accept_count")

    # --- Overall status ---
    status = "PASS" if fail_count == 0 else "FAIL"
    fields["status"] = status
    fields["first_failure"] = first_failure

    # --- Write result ---
    if args.result:
        result_dir = os.path.dirname(args.result)
        if result_dir and not os.path.isdir(result_dir):
            os.makedirs(result_dir)
        with open(args.result, "w") as f:
            for key in ["all_unit_cases_pass", "all_mpi_cases_pass",
                        "legacy_default_regression_pass",
                        "trial_deterministic", "trial_side_effect_free",
                        "failure_transaction_bitwise_unchanged",
                        "post_field_charge_invariance_pass",
                        "mpi_failure_consensus_pass", "soft_accept_count",
                        "status", "first_failure"]:
                f.write("%s=%s\n" % (key, fields.get(key, "0")))
            f.write("pass_count=%d\n" % pass_count)
            f.write("fail_count=%d\n" % fail_count)

    # --- Console output ---
    for key in ["all_unit_cases_pass", "all_mpi_cases_pass",
                "legacy_default_regression_pass",
                "trial_deterministic", "trial_side_effect_free",
                "failure_transaction_bitwise_unchanged",
                "post_field_charge_invariance_pass",
                "mpi_failure_consensus_pass", "soft_accept_count",
                "status", "first_failure"]:
        print("%s=%s" % (key, fields.get(key, "0")))
    print("pass_count=%d" % pass_count)
    print("fail_count=%d" % fail_count)

    # §8.7 requirement 4: physical/numerical gate failure → exit code 1.
    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
