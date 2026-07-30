#include "dual_u_coupling.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

size_t face_index(int ix, int jf, int k)
{
    return (static_cast<size_t>(ix) * (Param::Nv + 1) +
            static_cast<size_t>(jf)) * Param::Nmu + k;
}

bool close_relative(double a, double b, double tolerance)
{
    return std::fabs(a - b) <=
        tolerance * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

} // namespace

int main()
{
    const int nx = 1;
    const int jf = Param::Nv / 2;
    const int k = 0;
    const double h = 0.1;
    const double acceleration = 0.5;
    const size_t face_count =
        static_cast<size_t>(nx) * (Param::Nv + 1) * Param::Nmu;
    const size_t cell_count =
        static_cast<size_t>(nx) * Param::Nvmu;

    std::vector<double> weights(Param::Nvmu, 0.0);
    weights[static_cast<size_t>(jf) * Param::Nmu + k] = 1.0;
    std::vector<double> coefficient(face_count, 0.0);
    coefficient[face_index(0, jf, k)] = 1.0;
    std::vector<double> flux(face_count, 0.0);
    flux[face_index(0, jf, k)] = acceleration;
    std::vector<double> candidate(cell_count, 1.0);
    std::vector<double> target(1, 1.25);
    std::vector<double> acceleration_by_x(1, acceleration);
    std::vector<double> upar_widths(static_cast<size_t>(Param::Nv), 1.0);
    const double mass_before =
        std::accumulate(candidate.begin(), candidate.end(), 0.0);
    const double lower_boundary_before = coefficient[face_index(0, 0, k)];
    const double upper_boundary_before =
        coefficient[face_index(0, Param::Nv, k)];

    const DualUCoupling::FinalLimitedDiagnostics exact =
        DualUCoupling::apply_final_limited_capacity_pairing(
            nx, h, target, acceleration_by_x, weights,
            upar_widths,
            candidate, coefficient, flux);
    const double mass_after =
        std::accumulate(candidate.begin(), candidate.end(), 0.0);
    const bool exact_pass =
        exact.valid &&
        exact.limited_cell_count == 0 &&
        exact.residual_after_linf <= 1.0e-13 &&
        close_relative(mass_before, mass_after, 1.0e-14) &&
        coefficient[face_index(0, 0, k)] == lower_boundary_before &&
        coefficient[face_index(0, Param::Nv, k)] == upper_boundary_before;

    std::fill(coefficient.begin(), coefficient.end(), 0.0);
    coefficient[face_index(0, jf, k)] = 1.0;
    std::fill(flux.begin(), flux.end(), 0.0);
    flux[face_index(0, jf, k)] = acceleration;
    std::fill(candidate.begin(), candidate.end(), 1.0);
    candidate[static_cast<size_t>(jf - 1) * Param::Nmu + k] = 1.0e-10;
    target[0] = 2.0;
    const double limited_mass_before =
        std::accumulate(candidate.begin(), candidate.end(), 0.0);
    const DualUCoupling::FinalLimitedDiagnostics limited =
        DualUCoupling::apply_final_limited_capacity_pairing(
            nx, h, target, acceleration_by_x, weights,
            upar_widths,
            candidate, coefficient, flux);
    const double limited_mass_after =
        std::accumulate(candidate.begin(), candidate.end(), 0.0);
    const bool limited_pass =
        limited.valid &&
        limited.limited_cell_count == 1 &&
        limited.minimum_scale < 1.0 &&
        limited.residual_after_linf > 0.0 &&
        limited.candidate_min >= 0.0 &&
        close_relative(limited_mass_before, limited_mass_after, 1.0e-14) &&
        coefficient[face_index(0, 0, k)] == lower_boundary_before &&
        coefficient[face_index(0, Param::Nv, k)] == upper_boundary_before;

    // Production-like Maxwellian support: the occupied thermal body is
    // compact in floating point and every outer-tail cell is exactly zero.
    // Tail faces carry nonzero moment weights, so a whole-plane scale would
    // be pinned to zero; the active-set correction must instead use only
    // visible donors in the thermal body.
    std::fill(weights.begin(), weights.end(), 0.0);
    for (int face = 1; face < Param::Nv; ++face)
        weights[static_cast<size_t>(face) * Param::Nmu + k] =
            0.5 + static_cast<double>(face) / Param::Nv;
    std::fill(coefficient.begin(), coefficient.end(), 0.0);
    std::fill(flux.begin(), flux.end(), 0.0);
    std::fill(candidate.begin(), candidate.end(), 0.0);
    for (int j = jf - 4; j <= jf + 4; ++j) {
        const double offset = static_cast<double>(j - jf);
        candidate[static_cast<size_t>(j) * Param::Nmu + k] =
            std::exp(-0.5 * offset * offset);
    }
    const double zero_tail_mass_before =
        std::accumulate(candidate.begin(), candidate.end(), 0.0);
    target[0] = 0.4;
    const DualUCoupling::FinalLimitedDiagnostics zero_tail =
        DualUCoupling::apply_final_limited_capacity_pairing(
            nx, h, target, acceleration_by_x, weights,
            upar_widths,
            candidate, coefficient, flux);
    const double zero_tail_mass_after =
        std::accumulate(candidate.begin(), candidate.end(), 0.0);
    const bool zero_tail_pass =
        zero_tail.valid &&
        zero_tail.corrected_cell_count == 1 &&
        zero_tail.minimum_scale > 0.0 &&
        zero_tail.correction_linf > 0.0 &&
        zero_tail.residual_after_linf < zero_tail.residual_before_linf &&
        zero_tail.candidate_min >= 0.0 &&
        close_relative(zero_tail_mass_before, zero_tail_mass_after, 1.0e-14) &&
        coefficient[face_index(0, 0, k)] == lower_boundary_before &&
        coefficient[face_index(0, Param::Nv, k)] == upper_boundary_before;

    // A tiny high-weight donor and a well-populated reserve donor must be
    // pooled by capacity.  The old whole-plane scale was pinned by the tiny
    // donor; the capacity allocation reaches the target without exhausting it
    // first and restarting an active-set solve.
    std::fill(weights.begin(), weights.end(), 0.0);
    const int small_face = jf;
    const int reserve_face = jf + 3;
    weights[static_cast<size_t>(small_face) * Param::Nmu + k] = 1.0e6;
    weights[static_cast<size_t>(reserve_face) * Param::Nmu + k] = 1.0;
    std::fill(coefficient.begin(), coefficient.end(), 0.0);
    std::fill(flux.begin(), flux.end(), 0.0);
    std::fill(candidate.begin(), candidate.end(), 0.0);
    candidate[static_cast<size_t>(small_face - 1) * Param::Nmu + k] =
        1.0e-10;
    candidate[static_cast<size_t>(reserve_face - 1) * Param::Nmu + k] =
        1.0;
    const double active_set_mass_before =
        std::accumulate(candidate.begin(), candidate.end(), 0.0);
    target[0] = 1.0;
    const DualUCoupling::FinalLimitedDiagnostics active_set =
        DualUCoupling::apply_final_limited_capacity_pairing(
            nx, h, target, acceleration_by_x, weights,
            upar_widths,
            candidate, coefficient, flux);
    const double active_set_mass_after =
        std::accumulate(candidate.begin(), candidate.end(), 0.0);
    const bool active_set_pass =
        active_set.valid &&
        active_set.corrected_cell_count == 1 &&
        active_set.limited_cell_count == 0 &&
        active_set.minimum_scale == 1.0 &&
        active_set.residual_after_linf <= 1.0e-12 &&
        active_set.candidate_min >= 0.0 &&
        close_relative(active_set_mass_before, active_set_mass_after, 1.0e-14) &&
        coefficient[face_index(0, 0, k)] == lower_boundary_before &&
        coefficient[face_index(0, Param::Nv, k)] == upper_boundary_before;

    // A numerically tiny outer-tail donor must not control the correction even
    // when its moment weight is very large.  The populated thermal donor
    // carries the correction and the tail mass remains untouched.
    std::fill(weights.begin(), weights.end(), 0.0);
    const int tail_face = 1;
    const int core_face = jf;
    weights[static_cast<size_t>(tail_face) * Param::Nmu + k] = 1.0e20;
    weights[static_cast<size_t>(core_face) * Param::Nmu + k] = 1.0;
    std::fill(coefficient.begin(), coefficient.end(), 0.0);
    std::fill(flux.begin(), flux.end(), 0.0);
    std::fill(candidate.begin(), candidate.end(), 0.0);
    candidate[k] = 1.0e-20;
    candidate[static_cast<size_t>(core_face - 1) * Param::Nmu + k] = 1.0;
    const double tiny_tail_before = candidate[k];
    const double tail_exclusion_mass_before =
        std::accumulate(candidate.begin(), candidate.end(), 0.0);
    target[0] = 0.5;
    const DualUCoupling::FinalLimitedDiagnostics tail_exclusion =
        DualUCoupling::apply_final_limited_capacity_pairing(
            nx, h, target, acceleration_by_x, weights,
            upar_widths, candidate, coefficient, flux);
    const double tail_exclusion_mass_after =
        std::accumulate(candidate.begin(), candidate.end(), 0.0);
    const bool tail_exclusion_pass =
        tail_exclusion.valid &&
        tail_exclusion.residual_after_linf <= 1.0e-12 &&
        tail_exclusion.limited_cell_count == 0 &&
        candidate[k] == tiny_tail_before &&
        close_relative(tail_exclusion_mass_before,
                       tail_exclusion_mass_after, 1.0e-14);

    // Multi-cell audit: each x cell must close independently, including mixed
    // field and residual signs.  A small global sum cannot hide a bad cell.
    const int local_nx = 4;
    const size_t local_face_count =
        static_cast<size_t>(local_nx) * (Param::Nv + 1) * Param::Nmu;
    const size_t local_cell_count =
        static_cast<size_t>(local_nx) * Param::Nvmu;
    std::fill(weights.begin(), weights.end(), 0.0);
    weights[static_cast<size_t>(core_face) * Param::Nmu + k] = 1.0;
    std::vector<double> local_coefficient(local_face_count, 0.0);
    std::vector<double> local_flux(local_face_count, 0.0);
    std::vector<double> local_candidate(local_cell_count, 0.0);
    std::vector<double> local_target;
    local_target.push_back(0.2);
    local_target.push_back(-0.3);
    local_target.push_back(0.4);
    local_target.push_back(-0.1);
    std::vector<double> local_acceleration;
    local_acceleration.push_back(0.5);
    local_acceleration.push_back(-0.4);
    local_acceleration.push_back(0.8);
    local_acceleration.push_back(-0.6);
    std::vector<double> local_mass_before(local_nx, 0.0);
    for (int ix = 0; ix < local_nx; ++ix) {
        const size_t base = static_cast<size_t>(ix) * Param::Nvmu;
        local_candidate[base +
            static_cast<size_t>(core_face - 1) * Param::Nmu + k] =
            1.0 + 0.25 * ix;
        local_candidate[base +
            static_cast<size_t>(core_face) * Param::Nmu + k] =
            0.75 + 0.1 * ix;
        local_mass_before[ix] = std::accumulate(
            local_candidate.begin() + base,
            local_candidate.begin() + base + Param::Nvmu, 0.0);
    }
    const DualUCoupling::FinalLimitedDiagnostics local_pairing =
        DualUCoupling::apply_final_limited_capacity_pairing(
            local_nx, h, local_target, local_acceleration, weights,
            upar_widths, local_candidate, local_coefficient, local_flux);
    bool local_pairing_pass =
        local_pairing.valid &&
        local_pairing.corrected_cell_count == local_nx &&
        local_pairing.residual_after_linf <= 1.0e-12;
    for (int ix = 0; ix < local_nx; ++ix) {
        const size_t base = static_cast<size_t>(ix) * Param::Nvmu;
        const double mass_after = std::accumulate(
            local_candidate.begin() + base,
            local_candidate.begin() + base + Param::Nvmu, 0.0);
        local_pairing_pass = local_pairing_pass &&
            close_relative(local_mass_before[ix], mass_after, 1.0e-14) &&
            local_coefficient[face_index(ix, 0, k)] == 0.0 &&
            local_coefficient[face_index(ix, Param::Nv, k)] == 0.0;
    }

    std::cout << "exact_pass " << (exact_pass ? 1 : 0) << "\n"
              << "exact_residual_after " << exact.residual_after_linf << "\n"
              << "limited_pass " << (limited_pass ? 1 : 0) << "\n"
              << "limited_scale " << limited.minimum_scale << "\n"
              << "limited_residual_after "
              << limited.residual_after_linf << "\n"
              << "zero_tail_pass " << (zero_tail_pass ? 1 : 0) << "\n"
              << "zero_tail_scale " << zero_tail.minimum_scale << "\n"
              << "zero_tail_residual_before "
              << zero_tail.residual_before_linf << "\n"
              << "zero_tail_residual_after "
              << zero_tail.residual_after_linf << "\n"
              << "active_set_pass " << (active_set_pass ? 1 : 0) << "\n"
              << "active_set_scale " << active_set.minimum_scale << "\n"
              << "active_set_residual_after "
              << active_set.residual_after_linf << "\n"
              << "tail_exclusion_pass " << (tail_exclusion_pass ? 1 : 0)
              << "\n"
              << "tail_exclusion_residual_after "
              << tail_exclusion.residual_after_linf << "\n"
              << "local_pairing_pass " << (local_pairing_pass ? 1 : 0)
              << "\n"
              << "local_pairing_residual_after "
              << local_pairing.residual_after_linf << "\n";
    return exact_pass && limited_pass && zero_tail_pass && active_set_pass &&
        tail_exclusion_pass && local_pairing_pass ? 0 : 1;
}
