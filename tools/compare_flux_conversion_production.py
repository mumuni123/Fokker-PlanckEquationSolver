#!/usr/bin/env python3
"""Compare static-cell and flux-interface accepted production runs.

The comparison is deliberately based on accepted production rows.  It does
not infer conversion quality from a proxy field or from a rank-local audit.
"""

from __future__ import print_function

import argparse
import math
import os
import sys


def read_table(root, name):
    path = os.path.join(root, name)
    if not os.path.isfile(path):
        return []
    with open(path, "r") as stream:
        lines = [line.split() for line in stream if line.strip() and
                 not line.startswith("#")]
    if len(lines) < 2:
        return []
    header = lines[0]
    rows = []
    for fields in lines[1:]:
        if len(fields) != len(header):
            continue
        row = {}
        for key, value in zip(header, fields):
            try:
                row[key] = float(value)
            except ValueError:
                row[key] = value
        rows.append(row)
    return rows


def numeric(rows, key):
    return [row[key] for row in rows
            if key in row and isinstance(row[key], (int, float)) and
            math.isfinite(row[key])]


def latest(rows, key):
    values = numeric(rows, key)
    return values[-1] if values else math.nan


def maximum(rows, key, default=math.nan):
    values = numeric(rows, key)
    return max(values) if values else default


def total(rows, key):
    values = numeric(rows, key)
    return sum(values) if values else math.nan


def accepted_rows(rows):
    return [row for row in rows if row.get("accepted", 0.0) == 1.0]


def failure_file_present(root):
    path = os.path.join(root, "vpfp_failure.dat")
    return os.path.isfile(path) and os.path.getsize(path) > 0


def finite_values(*values):
    return all(math.isfinite(value) for value in values)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--static", required=True)
    parser.add_argument("--flux", required=True)
    parser.add_argument("--result", required=True)
    args = parser.parse_args()

    static_steps = accepted_rows(read_table(
        args.static, "vpfp_step_diagnostics.dat"))
    flux_steps = accepted_rows(read_table(
        args.flux, "vpfp_step_diagnostics.dat"))
    flux_rows = read_table(args.flux, "bulk_tail_flux_accepted_steps.dat")

    exported = total(flux_rows, "exported_number")
    face_exported = total(flux_rows, "face_export_number")
    parcel_number = total(flux_rows, "parcel_number")
    duplicate = total(flux_rows, "duplicate_count")
    below = total(flux_rows, "below_threshold_number")
    static_extractor = total(flux_rows, "static_extractor_call_count")
    created = total(flux_rows, "particles_created")
    finite = finite_values(exported, face_exported, parcel_number,
                           duplicate, below, static_extractor, created)

    parcel_vs_face = (abs(face_exported - parcel_number) /
                      max(1.0, abs(face_exported))
                      if finite else math.nan)
    conversion_number = maximum(flux_rows, "conversion_number_residual")
    conversion_px = maximum(flux_rows, "conversion_px_residual")
    conversion_energy = maximum(flux_rows, "conversion_energy_residual")
    conversion_jx = maximum(flux_rows, "conversion_jx_residual")
    conversion_pixx = maximum(flux_rows, "conversion_pixx_residual")
    conversion_piperp = maximum(flux_rows, "conversion_piperp_residual")
    conversion_rho = maximum(flux_rows, "conversion_rho_l2")
    combined_number = maximum(flux_steps, "tail_number_balance_error")

    conversion_npk = max(conversion_number, conversion_px,
                          conversion_energy)
    conversion_jmom = max(conversion_jx, conversion_pixx,
                           conversion_piperp)
    static_gauss = maximum(static_steps, "gauss_linf", 0.0)
    flux_gauss = maximum(flux_steps, "gauss_linf", 0.0)
    gauss_ratio = ((flux_gauss + 1.0e-300) /
                   max(1.0e-300, static_gauss))
    static_particles = latest(static_steps, "tail_particle_count")
    flux_particles = latest(flux_steps, "tail_particle_count")
    particle_ratio = (flux_particles / max(1.0, static_particles)
                      if finite_values(static_particles, flux_particles)
                      else math.nan)
    static_wall = total(static_steps, "wall_s")
    flux_wall = total(flux_steps, "wall_s")
    wall_ratio = (flux_wall / static_wall
                  if finite_values(static_wall, flux_wall) and
                  static_wall > 0.0 else math.nan)
    below_relative = below / max(1.0, abs(exported)) if finite else math.nan
    failures = failure_file_present(args.flux)
    no_static_extractor = math.isfinite(static_extractor) and \
        static_extractor == 0.0
    accepted_step_count_equal = (bool(static_steps) and bool(flux_steps) and
                                 len(static_steps) == len(flux_steps))
    flux_all_steps_accepted = (bool(flux_steps) and
                               accepted_step_count_equal and
                               not failures)

    required_present = all(math.isfinite(value) for value in (
        combined_number, conversion_npk, conversion_jmom,
        conversion_piperp, gauss_ratio, particle_ratio, wall_ratio,
        below_relative, parcel_vs_face, conversion_rho))
    accepted = bool(static_steps) and bool(flux_steps) and bool(flux_rows)
    gates = {
        "combined_number": combined_number <= 1.0e-10,
        "conversion_npk": conversion_npk <= 1.0e-10,
        "conversion_jmom": conversion_jmom <= 1.0e-9,
        "conversion_piperp": conversion_piperp <= 1.0e-9,
        "below_threshold": below_relative <= 1.0e-12,
        "duplicate": duplicate == 0.0,
        "parcel_face": parcel_vs_face <= 1.0e-13,
        "gauss": gauss_ratio <= 1.1,
        "particles": particle_ratio <= 2.0,
        "wall": wall_ratio <= 2.0,
        "all_steps_accepted": flux_all_steps_accepted,
        "accepted_step_count_equal": accepted_step_count_equal,
        "no_static_extractor": no_static_extractor,
        "no_failure": not failures,
    }
    ok = accepted and finite and required_present and all(gates.values())
    gauss_not_regressed = gauss_ratio <= 1.1
    piperp_old_residual_removed = conversion_piperp <= 1.0e-9

    with open(args.result, "w") as out:
        out.write("static_accepted_steps=%d\n" % len(static_steps))
        out.write("flux_accepted_steps=%d\n" % len(flux_steps))
        out.write("flux_audit_rows=%d\n" % len(flux_rows))
        out.write("flux_all_steps_accepted=%d\n" %
                  int(flux_all_steps_accepted))
        out.write("accepted_step_count_equal=%d\n" %
                  int(accepted_step_count_equal))
        out.write("flux_static_extractor_call_count=%.17g\n" % static_extractor)
        out.write("flux_static_extractor_absent=%d\n" %
                  int(no_static_extractor))
        out.write("flux_exported_number=%.17g\n" % exported)
        out.write("flux_face_export_number=%.17g\n" % face_exported)
        out.write("flux_parcel_number=%.17g\n" % parcel_number)
        out.write("parcel_vs_face_flux_relative_l1=%.17g\n" % parcel_vs_face)
        out.write("flux_particles_created=%.17g\n" % created)
        out.write("flux_duplicate_count=%.17g\n" % duplicate)
        out.write("interface_duplicate_count=%.17g\n" % duplicate)
        out.write("flux_below_threshold_number=%.17g\n" % below)
        out.write("combined_number_relative_error=%.17g\n" % combined_number)
        out.write("conversion_N_Px_K_relative_l1=%.17g\n" % conversion_npk)
        out.write("conversion_Jx_Pixx_Piperp_relative_l1=%.17g\n" % conversion_jmom)
        out.write("conversion_piperp_old_residual=%.17g\n" % conversion_piperp)
        out.write("conversion_piperp_old_residual_removed=%d\n" %
                  int(piperp_old_residual_removed))
        out.write("conversion_rho_l2=%.17g\n" % conversion_rho)
        out.write("below_threshold_number_relative=%.17g\n" % below_relative)
        out.write("gauss_error_ratio_flux_over_static=%.17g\n" % gauss_ratio)
        out.write("gauss_error_not_regressed=%d\n" %
                  int(gauss_not_regressed))
        out.write("tail_particle_count_ratio_flux_over_static=%.17g\n" % particle_ratio)
        out.write("step_wall_time_ratio_flux_over_static=%.17g\n" % wall_ratio)
        out.write("required_metrics_present=%d\n" % int(required_present))
        for name, passed in sorted(gates.items()):
            out.write("gate_%s=%d\n" % (name, int(passed)))
        out.write("status=%s\n" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
