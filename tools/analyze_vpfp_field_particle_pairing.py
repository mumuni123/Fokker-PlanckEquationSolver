#!/usr/bin/env python3
"""Gate I accepted-state pairing analyzer.

Exit codes: PASS=0, physical FAIL/inconclusive=2, invalid input/schema=3.
"""

import argparse
import glob
import math
import os
import sys

PAIRING_FILE = "field_particle_power_pairing.dat"
STEP_FILE = "vpfp_step_diagnostics.dat"
THRESHOLDS = {"A": .50, "B": .50, "C": .50, "D": .50,
              "E": .20, "F": .80}
REQUIRED = {
    "step", "time_s", "dt_s", "accepted", "continuity_bulk",
    "continuity_tail", "continuity_beam", "poisson_transport_residual",
    "poisson_endpoint", "poisson_midpoint", "poisson_discrete_gradient",
    "force_work_bulk", "force_work_tail", "force_work_beam",
    "current_pair_residual", "conversion_residual", "boundary_residual",
    "full_residual", "reconstructed_residual", "reconstruction_mismatch",
    "poisson_identity_crosscheck", "poisson_identity_scale",
    "poisson_crosscheck_tolerance", "current_pair_linf",
    "roundoff_tolerance", "all_finite", "continuity_pass",
    "local_work_ledger_pass", "reconstruction_pass", "root_cause_mask",
}


def read_table(path, plain_header=False):
    if not os.path.isfile(path):
        return None, None, "missing %s" % os.path.basename(path)
    columns, rows = None, []
    try:
        with open(path, encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    if line.startswith("# columns="):
                        columns = line[len("# columns="):].split()
                    continue
                parts = line.split()
                if columns is None and plain_header:
                    columns, plain_header = parts, False
                    continue
                if columns is None:
                    return None, None, "missing columns declaration"
                if len(parts) < len(columns):
                    return None, None, "short row %d < %d" % (
                        len(parts), len(columns))
                row = {}
                for name, value in zip(columns, parts):
                    try:
                        row[name] = float(value)
                    except ValueError:
                        return None, None, "non-numeric %s=%s" % (name, value)
                rows.append(row)
    except OSError as exc:
        return None, None, str(exc)
    return (columns, rows, None) if rows else (columns, None, "no data rows")


def read_manifest(directory):
    result = {}
    for name in ("manifest.txt", "manifest.dat"):
        path = os.path.join(directory, name)
        if not os.path.isfile(path):
            continue
        with open(path, encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("=", 1) if "=" in line else line.split(None, 1)
                if len(parts) == 2:
                    result[parts[0].strip()] = parts[1].strip()
        break
    return result


def read_source_checkpoint(directory):
    """Read the *source* checkpoint ``manifest.txt`` identity (section: fix
    checkpoint identification).  The production checkpoint manifest is a
    whitespace ``key value`` table written by vpfp_checkpoint.cpp."""
    if not directory:
        return {}
    path = os.path.join(directory, "manifest.txt")
    if not os.path.isfile(path):
        return {}
    result = {}
    try:
        with open(path, encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(None, 1)
                if len(parts) == 2:
                    result[parts[0].strip()] = parts[1].strip()
    except OSError:
        return {}
    return result


def source_identity(manifest):
    """Extract the five identity fields required by the fix.  Values are
    ``None`` when missing or unparsable (never coerced to 0)."""
    def as_int(key):
        try:
            return int(manifest[key])
        except (KeyError, ValueError):
            return None

    def as_float(key):
        try:
            return float(manifest[key])
        except (KeyError, ValueError):
            return None

    return {
        "step": as_int("step"),
        "time": as_float("time"),
        "physical_config_hash": manifest.get("physical_config_hash"),
        "mpi_size": as_int("mpi_size"),
        "nx_global": as_int("nx_global"),
    }


def time_close(a, b):
    scale = max(abs(a), abs(b), 1.e-30)
    return abs(a - b) <= 256. * sys.float_info.epsilon * scale


def convergence_order(abs_coarse, abs_fine, dt_coarse, dt_fine):
    """Observed convergence order p in |R| ~ dt^p between coarse and fine.
    Near-zero order means a fixed structural incompatibility, not a
    discretization error; returns None when the estimate is undefined."""
    if not (abs_coarse > 0.0 and abs_fine > 0.0 and
            dt_coarse > 0.0 and dt_fine > 0.0 and dt_fine < dt_coarse):
        return None
    return math.log(abs_fine / abs_coarse) / math.log(dt_fine / dt_coarse)


def term_metrics(rows):
    """Per-term signed sum, absolute sum, independent exchange scale and the
    A-F explanation fractions (section I6 fix).  The fractions use the
    independent absolute exchange scale, never the strongly-cancelled
    full_residual."""
    def ssum(col):
        return sum(row[col] for row in rows)

    def asum(col):
        return sum(abs(row[col]) for row in rows)

    def diff_list(a, b):
        return [row[a] - row[b] for row in rows]

    signed = {
        "A": ssum("poisson_transport_residual"),
        "B": ssum("current_pair_residual"),
        "C": sum(diff_list("poisson_discrete_gradient", "poisson_midpoint")),
        "D": sum(row["force_work_bulk"] + row["force_work_tail"] +
                 row["force_work_beam"] for row in rows),
        "E": ssum("conversion_residual"),
        "boundary": ssum("boundary_residual"),
    }
    absolute = {
        "A": asum("poisson_transport_residual"),
        "B": asum("current_pair_residual"),
        "C": sum(abs(v) for v in diff_list("poisson_discrete_gradient",
                                           "poisson_midpoint")),
        "D": sum(abs(row["force_work_bulk"]) + abs(row["force_work_tail"]) +
                 abs(row["force_work_beam"]) for row in rows),
        "E": asum("conversion_residual"),
        "boundary": asum("boundary_residual"),
    }
    # Independent physical exchange scale: sum of |identity terms| (A + B + D
    # + E + boundary).  C is a candidate-comparison difference, not an
    # identity term, so it is excluded from the denominator.
    scale = (absolute["A"] + absolute["B"] + absolute["D"] +
             absolute["E"] + absolute["boundary"])
    denom = scale if scale > 0.0 else 1.0
    fractions = {"A": absolute["A"] / denom,
                 "B": absolute["B"] / denom,
                 "C": absolute["C"] / denom,
                 "D": absolute["D"] / denom,
                 "E": absolute["E"] / denom}
    return signed, absolute, scale, fractions


def read_profiles(directory):
    region_abs = [0., 0., 0.]
    found = False
    pattern = os.path.join(directory,
                           "field_particle_power_pairing_profile_*_rank*.dat")
    for path in glob.glob(pattern):
        columns, rows, err = read_table(path)
        if err or "pairing_residual_density" not in (columns or []) or \
                "region_id" not in (columns or []):
            continue
        found = True
        for row in rows:
            region = int(row["region_id"])
            if 0 <= region < 3:
                region_abs[region] += abs(row["pairing_residual_density"])
    total = sum(region_abs)
    return {"available": found and total > 0.,
            "left": region_abs[0] / total if total else 0.,
            "core": region_abs[1] / total if total else 0.,
            "right": region_abs[2] / total if total else 0.,
            "boundary": (region_abs[0] + region_abs[2]) / total if total else 0.}


def load_run(directory, expected_steps=None, require_no_split=False):
    columns, rows, err = read_table(os.path.join(directory, PAIRING_FILE))
    if err:
        return None, err
    missing = sorted(REQUIRED - set(columns))
    if missing:
        return None, "missing schema columns: %s" % ",".join(missing)
    for i, row in enumerate(rows):
        for name in REQUIRED:
            if not math.isfinite(row[name]):
                return None, "non-finite %s at row %d" % (name, i)
        if int(row["accepted"]) != 1:
            return None, "non-accepted pairing row %d" % i
        if i and int(row["step"]) != int(rows[i - 1]["step"]) + 1:
            return None, "step gap at row %d" % i
    if expected_steps is not None and len(rows) != expected_steps:
        return None, "accepted steps %d != expected %d" % (len(rows), expected_steps)

    _, step_rows, step_err = read_table(os.path.join(directory, STEP_FILE), True)
    if step_err:
        return None, step_err
    if len(step_rows) != len(rows):
        return None, "pairing/step row count mismatch %d/%d" % (
            len(rows), len(step_rows))
    if any(int(row.get("accepted", 0)) != 1 for row in step_rows):
        return None, "step diagnostics contains rejected row"
    if require_no_split and any(int(row.get("split", 0)) for row in step_rows):
        return None, "split step present"
    failure = os.path.join(directory, "vpfp_failure.dat")
    if os.path.isfile(failure) and os.path.getsize(failure):
        return None, "failure ledger is non-empty"

    def asum(name):
        return sum(abs(row[name]) for row in rows)

    full_abs = asum("full_residual")
    profiles = read_profiles(directory)
    signed, absolute, scale, fractions = term_metrics(rows)
    # F is a region/localization qualifier, not an identity term.
    fractions["F"] = profiles["boundary"] if profiles["available"] else 0.0
    structure = all(
        int(row[flag]) == 1
        for row in rows
        for flag in ("all_finite", "continuity_pass",
                     "local_work_ledger_pass", "reconstruction_pass"))
    structure = structure and all(
        abs(row["reconstruction_mismatch"]) <= row["roundoff_tolerance"] and
        abs(row["poisson_identity_crosscheck"]) <=
            row["poisson_crosscheck_tolerance"]
        for row in rows)
    worst = sorted(rows, key=lambda row: abs(row["full_residual"]),
                   reverse=True)[:10]
    total_dt = sum(row["dt_s"] for row in rows)
    return {"rows": rows, "accepted_steps": len(rows),
            "start": rows[0]["time_s"] - rows[0]["dt_s"],
            "end": rows[-1]["time_s"], "full_signed":
            sum(row["full_residual"] for row in rows), "full_abs": full_abs,
            "full_power": full_abs / max(total_dt, sys.float_info.min),
            "dt_mean": total_dt / max(1, len(rows)),
            "structure": structure, "fractions": fractions,
            "signed": signed, "absolute": absolute, "scale": scale,
            "profiles": profiles,
            "worst": [int(row["step"]) for row in worst],
            "manifest": read_manifest(directory)}, None


def hits(run):
    return {key for key, value in run["fractions"].items()
            if value >= THRESHOLDS[key]}


def select_root(coarse, fine):
    selected = hits(coarse) & hits(fine)
    qualifier = ""
    if "F" in selected:
        selected.remove("F")
        if coarse["profiles"]["core"] < .10 and fine["profiles"]["core"] < .10:
            qualifier = "_BOUNDARY"
    # Section I6 fix: A/B/C hitting simultaneously is inconclusive; never pick
    # the largest fraction as the root cause.
    if {"A", "B", "C"} <= selected:
        return "none", len(selected) + (1 if qualifier else 0)
    if {"A", "C"} <= selected:
        def c_signature(run):
            return all(abs(row["poisson_discrete_gradient"]) <=
                       row["roundoff_tolerance"] and
                       abs(row["poisson_midpoint"]) > row["roundoff_tolerance"]
                       for row in run["rows"])
        selected.discard("A" if c_signature(coarse) and c_signature(fine)
                         else "C")
    if {"B", "D"} <= selected:
        all_rows = coarse["rows"] + fine["rows"]
        bulk = sum(abs(row["force_work_bulk"]) for row in all_rows)
        pic = sum(abs(row["force_work_tail"]) + abs(row["force_work_beam"])
                  for row in all_rows)
        selected.discard("B" if pic > bulk else "D")
    if len(selected) == 1:
        return next(iter(selected)) + qualifier, 1
    return "none", len(selected)


def write_result(path, fields):
    if not path:
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for key in sorted(fields):
            handle.write("%s=%s\n" % (key, fields[key]))


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--run")
    parser.add_argument("--coarse")
    parser.add_argument("--fine")
    parser.add_argument("--source-checkpoint",
                        help="source checkpoint directory whose manifest.txt "
                             "provides step/time/physical_config_hash/mpi_size/"
                             "nx_global (the true experiment identity)")
    parser.add_argument("--coarse-dt-scale", type=float)
    parser.add_argument("--fine-dt-scale", type=float)
    parser.add_argument("--expected-accepted-steps", type=int)
    parser.add_argument("--require-no-split", action="store_true")
    parser.add_argument("--baseline")
    parser.add_argument("--result")
    args = parser.parse_args(argv)

    if args.run and not (args.coarse or args.fine):
        run, err = load_run(args.run, args.expected_accepted_steps,
                            args.require_no_split)
        if err:
            fields = {"status": "INVALID_INPUT", "reason": err}
            write_result(args.result, fields)
            print("status=INVALID_INPUT reason=%s" % err)
            return 3
        status = "PASS" if run["structure"] else "INVALID_AUDIT_STRUCTURE"
        fields = {"status": status, "accepted_steps": run["accepted_steps"],
                  "full_residual_signed": run["full_signed"],
                  "full_residual_abs": run["full_abs"],
                  "full_residual_power": run["full_power"],
                  "worst_steps": ",".join(map(str, run["worst"]))}
        fields.update(("fraction_" + key, value)
                      for key, value in run["fractions"].items())
        write_result(args.result, fields)
        print("status=%s" % status)
        return 0 if status == "PASS" else 3

    if not (args.coarse and args.fine):
        print("status=INVALID_INPUT reason=need --run or --coarse/--fine")
        return 3
    coarse, err_c = load_run(args.coarse, args.expected_accepted_steps,
                             args.require_no_split)
    fine, err_f = load_run(args.fine, args.expected_accepted_steps,
                           args.require_no_split)
    if err_c or err_f:
        fields = {"status": "INVALID_COMPARISON",
                  "coarse_reason": err_c or "ok", "fine_reason": err_f or "ok"}
        write_result(args.result, fields)
        print("status=INVALID_COMPARISON")
        return 3

    # Checkpoint identification (fix 1): read the *source* checkpoint
    # manifest.txt, not the coarse/fine output-directory manifests.  The two
    # runs restart from the same checkpoint, so their initial physical state
    # is identical by construction; the controlled variable is only the dt
    # scale, which must not be treated as a physical-configuration change.
    source_path = args.source_checkpoint or ""
    identity = source_identity(read_source_checkpoint(source_path))
    identity_complete = all(
        identity.get(key) is not None
        for key in ("step", "time", "physical_config_hash", "mpi_size",
                    "nx_global"))

    same_state = identity["physical_config_hash"] is not None
    same_window = False
    if identity["time"] is not None:
        source_time = identity["time"]
        same_window = (time_close(coarse["start"], source_time) and
                       time_close(fine["start"], source_time) and
                       time_close(coarse["end"], fine["end"]))

    if not identity_complete or not same_window or not same_state:
        status, root, count = "INVALID_COMPARISON", "none", 0
    elif not coarse["structure"] or not fine["structure"]:
        status, root, count = "INVALID_AUDIT_STRUCTURE", "none", 0
    else:
        root, count = select_root(coarse, fine)
        status = ("PASS_ROOT_CAUSE_UNIQUE" if count == 1 else
                  "INCONCLUSIVE_UNEXPLAINED" if count == 0 else
                  "INCONCLUSIVE_MULTIPLE_CAUSES")
    fields = {"status": status, "root_cause": root,
              "root_cause_candidate_count": count,
              "same_initial_physical_state": int(same_state),
              "same_physical_window": int(same_window),
              "source_checkpoint_path": source_path,
              "source_checkpoint_step":
                  "" if identity["step"] is None else identity["step"],
              "source_checkpoint_time":
                  "" if identity["time"] is None else identity["time"],
              "source_physical_config_hash":
                  "" if identity["physical_config_hash"] is None
                  else identity["physical_config_hash"],
              "coarse_worst_steps": ",".join(map(str, coarse["worst"])),
              "fine_worst_steps": ",".join(map(str, fine["worst"]))}
    for key in THRESHOLDS:
        fields["fraction_%s_coarse" % key] = coarse["fractions"][key]
        fields["fraction_%s_fine" % key] = fine["fractions"][key]
    # Section I6 fix: report signed sum, absolute sum, independent exchange
    # scale and the coarse/fine convergence order per candidate.
    for key in ("A", "B", "C", "D", "E"):
        fields["signed_%s_coarse" % key] = coarse["signed"][key]
        fields["signed_%s_fine" % key] = fine["signed"][key]
        fields["abs_%s_coarse" % key] = coarse["absolute"][key]
        fields["abs_%s_fine" % key] = fine["absolute"][key]
        order = convergence_order(coarse["absolute"][key],
                                  fine["absolute"][key],
                                  coarse["dt_mean"], fine["dt_mean"])
        fields["order_%s" % key] = "nan" if order is None else order
    fields["exchange_scale_coarse"] = coarse["scale"]
    fields["exchange_scale_fine"] = fine["scale"]
    full_order = convergence_order(coarse["full_abs"], fine["full_abs"],
                                   coarse["dt_mean"], fine["dt_mean"])
    fields["order_full_residual"] = "nan" if full_order is None else full_order
    write_result(args.result, fields)
    print("status=%s root_cause=%s" % (status, root))
    return 0 if status == "PASS_ROOT_CAUSE_UNIQUE" else \
           (3 if status.startswith("INVALID") else 2)


if __name__ == "__main__":
    sys.exit(main())
