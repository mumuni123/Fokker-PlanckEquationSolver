#ifndef GRID_H
#define GRID_H

#include "parameters.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct VelocityGrid {
    double v_min, v_max, dv;
    double mu_min, mu_max, dmu;
    double inv_dv, inv_dmu;
    double max_abs_mu;
    double max_mu_face_factor;
    double max_inv_v;
    double max_speed;

    std::vector<double> v_cells;
    std::vector<double> v_faces;
    std::vector<double> v_widths;
    std::vector<double> inv_v_widths;
    std::vector<double> inv_v_center_dist;
    std::vector<double> mu_cells;
    std::vector<double> mu_faces;
    std::vector<double> v2_cells;
    std::vector<double> inv_v_cells;
    std::vector<double> inv_v2_cells;
    std::vector<double> v2_faces;
    std::vector<double> gamma_cells;
    std::vector<double> beta_cells;
    std::vector<double> speed_cells;
    std::vector<double> chain_speed_cells;
    std::vector<double> chain_speed_faces;
    std::vector<double> mu_face_factor;
    std::vector<double> moment_weight;
    std::vector<double> inv_moment_weight;
    std::vector<double> mu_flux_scale;
    std::vector<double> mu_cfl_factor;
    std::vector<double> vx_cells;
    std::vector<double> current_weight;

    void init(double umax)
    {
        v_min = 0.0;
        v_max = umax;
        dv = (v_max - v_min) / Param::Nv;
        inv_dv = 1.0 / dv;
        mu_min = -1.0;
        mu_max = 1.0;
        dmu = (mu_max - mu_min) / Param::Nmu;
        inv_dmu = 1.0 / dmu;

        v_cells.resize(Param::Nv);
        v_faces.resize(Param::Nv + 1);
        v_widths.resize(Param::Nv);
        inv_v_widths.resize(Param::Nv);
        inv_v_center_dist.resize(Param::Nv + 1);
        mu_cells.resize(Param::Nmu);
        mu_faces.resize(Param::Nmu + 1);
        v2_cells.resize(Param::Nv);
        inv_v_cells.resize(Param::Nv);
        inv_v2_cells.resize(Param::Nv);
        v2_faces.resize(Param::Nv + 1);
        gamma_cells.resize(Param::Nv);
        beta_cells.resize(Param::Nv);
        speed_cells.resize(Param::Nv);
        chain_speed_cells.resize(Param::Nv);
        chain_speed_faces.resize(Param::Nv + 1);
        mu_face_factor.resize(Param::Nmu + 1);
        moment_weight.resize(Param::Nv);
        inv_moment_weight.resize(Param::Nv);
        mu_flux_scale.resize(Param::Nv);
        mu_cfl_factor.resize(Param::Nv);
        vx_cells.resize(Param::Nvmu);
        current_weight.resize(Param::Nvmu);

        max_abs_mu = 0.0;
        max_mu_face_factor = 0.0;
        max_inv_v = 0.0;
        max_speed = 0.0;

        const double refined_u =
            std::min(Param::momentum_refined_u, v_max);
        const int refined_cells =
            std::max(1, std::min(Param::momentum_refined_cells, Param::Nv - 1));
        for (int iv = 0; iv <= Param::Nv; ++iv) {
            if (iv <= refined_cells) {
                v_faces[iv] =
                    refined_u * static_cast<double>(iv) / refined_cells;
            } else {
                const int coarse_cells = Param::Nv - refined_cells;
                const int coarse_iv = iv - refined_cells;
                v_faces[iv] =
                    refined_u + (v_max - refined_u)
                              * static_cast<double>(coarse_iv) / coarse_cells;
            }
        }
        v_faces[0] = v_min;
        v_faces[Param::Nv] = v_max;

        for (int iv = 0; iv < Param::Nv; ++iv) {
            const double u_left = v_faces[iv];
            const double u_right = v_faces[iv + 1];
            const double width = u_right - u_left;
            const double u = 0.5 * (u_left + u_right);
            const double u_eff = std::max(u, Param::u_floor);
            const double gamma = std::sqrt(1.0 + u * u);
            const double beta = u / gamma;
            const double gamma_left = std::sqrt(1.0 + u_left * u_left);
            const double gamma_right = std::sqrt(1.0 + u_right * u_right);
            v_cells[iv] = u;
            v_widths[iv] = width;
            inv_v_widths[iv] = 1.0 / width;
            v2_cells[iv] = u_eff * u_eff;
            inv_v_cells[iv] = 1.0 / u_eff;
            inv_v2_cells[iv] = 1.0 / (u_eff * u_eff);
            gamma_cells[iv] = gamma;
            beta_cells[iv] = beta;
            speed_cells[iv] = Const::c * beta;
            chain_speed_cells[iv] =
                (width > 0.0)
                ? Const::c * (gamma_right - gamma_left) / width
                : speed_cells[iv];
            moment_weight[iv] =
                2.0 * Const::pi * dmu
                * (u_right * u_right * u_right
                   - u_left * u_left * u_left) / 3.0;
            inv_moment_weight[iv] =
                (moment_weight[iv] > 0.0) ? 1.0 / moment_weight[iv] : 0.0;
            mu_flux_scale[iv] = moment_weight[iv] * inv_dmu;
            max_inv_v = std::max(max_inv_v, inv_v_cells[iv]);
            max_speed = std::max(max_speed, speed_cells[iv]);
        }
        inv_v_center_dist[0] = inv_v_widths[0];
        inv_v_center_dist[Param::Nv] = inv_v_widths[Param::Nv - 1];
        for (int iv = 1; iv < Param::Nv; ++iv) {
            inv_v_center_dist[iv] =
                1.0 / (v_cells[iv] - v_cells[iv - 1]);
        }
        for (int iv = 0; iv <= Param::Nv; ++iv) {
            v2_faces[iv] = v_faces[iv] * v_faces[iv];
        }
        chain_speed_faces[0] = chain_speed_cells[0];
        chain_speed_faces[Param::Nv] = chain_speed_cells[Param::Nv - 1];
        for (int iv = 1; iv < Param::Nv; ++iv) {
            const double du = v_cells[iv] - v_cells[iv - 1];
            chain_speed_faces[iv] =
                (du > 0.0)
                ? Const::c * (gamma_cells[iv] - gamma_cells[iv - 1]) / du
                : 0.5 * (speed_cells[iv] + speed_cells[iv - 1]);
        }
        for (int imu = 0; imu < Param::Nmu; ++imu) {
            const double mu = mu_min + (imu + 0.5) * dmu;
            mu_cells[imu] = mu;
            max_abs_mu = std::max(max_abs_mu, std::fabs(mu));
        }
        for (int imu = 0; imu <= Param::Nmu; ++imu) {
            const double muf = mu_min + imu * dmu;
            mu_faces[imu] = muf;
            mu_face_factor[imu] = 1.0 - muf * muf;
            max_mu_face_factor = std::max(max_mu_face_factor, mu_face_factor[imu]);
        }
        for (int iv = 0; iv < Param::Nv; ++iv) {
            mu_cfl_factor[iv] =
                max_mu_face_factor * inv_v_cells[iv] * inv_dmu;
        }
        for (int iv = 0; iv < Param::Nv; ++iv) {
            for (int imu = 0; imu < Param::Nmu; ++imu) {
                const size_t k = static_cast<size_t>(iv) * Param::Nmu + imu;
                vx_cells[k] = speed_cells[iv] * mu_cells[imu];
                current_weight[k] =
                    moment_weight[iv] * chain_speed_cells[iv] * mu_cells[imu];
            }
        }
    }

    double v(int i) const { return v_cells[i]; }
    double v_face(int i) const { return v_faces[i]; }
    double mu(int j) const { return mu_cells[j]; }
    double mu_face(int j) const { return mu_faces[j]; }
};

// Cylindrical momentum grid used by the collisionless background-electron
// Vlasov--Ampere path.  The two stored indices retain Param::Nv/Param::Nmu
// extents, but represent (u_parallel, u_perp), not the legacy (u, mu) grid.
// Values in Species::f are cell-integrated masses M on this grid.
struct CylindricalVelocityGrid {
    struct NestedGridAudit {
        double core_face_identity_linf;
        double core_cell_identity_linf;
        double core_vx_identity_linf;
        double core_kinetic_energy_identity_linf;
        double symmetry_linf;
        double max_adjacent_width_ratio;
        double phase_volume_relative_error;
    } nested_audit;

    std::vector<double> upar_faces;
    std::vector<double> upar_cells;
    std::vector<double> upar_widths;
    std::vector<double> upar_center_distances;
    std::vector<double> uperp_faces;
    std::vector<double> uperp_cells;
    std::vector<double> uperp_widths;
    std::vector<double> uperp_ring_areas;
    std::vector<double> kinetic_energy;
    std::vector<double> vx;
    // P3-V.2: read-only energy-conjugate velocity representation.  Face 0
    // and face Nv are deliberately marked invalid because they are velocity
    // domain boundaries, not internal kinetic-energy exchange faces.
    std::vector<double> vx_energy_conjugate_face;
    std::vector<unsigned char> vx_energy_conjugate_face_valid;
    std::vector<double> vx_energy_conjugate_cell;
    // Moment-closure self-consistency table (section 19.3): theta points
    // and their discrete-Maxwellian second moments, built once at init.
    std::vector<double> mc_uth2_table_;
    std::vector<double> mc_moment_table_;

    void init(double umax)
    {
        init_grid(umax, Param::Nv, Param::Nmu, Param::Nv_core, Param::Nv_tail,
                  Param::momentum_upar_core_max,
                  Param::momentum_upar_extended_max, 1.25);
    }

    // Runtime-resolution grid used by the phase-2B convergence test
    // (section 16.4.2): the same smooth u_parallel mapping with nv cells and
    // no appended tail, so 96 cells reproduces the production grid bitwise.
    // The adjacent-width-ratio cap is relaxed for coarse convergence grids
    // (48 cells with the production stretch exceeds 1.25 near u=10); all
    // other validations (symmetry, phase volume, identity) stay active.
    void init(double umax, int nv, int nmu)
    {
        init_grid(umax, nv, nmu, nv, 0,
                  Param::momentum_upar_core_max,
                  Param::momentum_upar_core_max, 2.0);
    }

    void init_grid(double umax, int nv, int nmu, int nv_core, int nv_tail,
                   double core_max, double extended_max, double ratio_cap)
    {
        upar_faces.resize(static_cast<size_t>(nv) + 1);
        upar_cells.resize(static_cast<size_t>(nv));
        upar_widths.resize(static_cast<size_t>(nv));
        upar_center_distances.assign(static_cast<size_t>(nv) + 1, 0.0);
        uperp_faces.resize(static_cast<size_t>(nmu) + 1);
        uperp_cells.resize(static_cast<size_t>(nmu));
        uperp_widths.resize(static_cast<size_t>(nmu));
        uperp_ring_areas.resize(static_cast<size_t>(nmu));
        kinetic_energy.resize(static_cast<size_t>(nv) * nmu);
        vx.resize(static_cast<size_t>(nv) * nmu);

        if (nv_core <= 0 || nv_tail < 0 ||
            nv_core % 2 != 0 || nmu <= 0 ||
            core_max <= 0.0 ||
            extended_max < core_max ||
            (nv_tail == 0 && extended_max != core_max) ||
            (nv_tail > 0 && extended_max <= core_max)) {
            throw std::runtime_error("invalid nested u_parallel grid configuration");
        }
        assert(nv > 0 && nmu > 0);

        build_upar_core_faces(nv_core, nv_tail, core_max);
        append_symmetric_upar_tail_faces(nv_core, nv_tail, core_max,
                                         extended_max);

        const double sinh_uperp = std::sinh(Param::momentum_uperp_stretch);
        build_uperp_faces(umax, sinh_uperp, nmu);
        build_cell_geometry_and_moments(nv, nmu);
        build_energy_conjugate_velocity_table();
        build_moment_closure_table(
            Param::temperature_e / (Const::me * Const::c * Const::c));
        validate_nested_grid(umax, nv, nmu, nv_core, nv_tail, core_max,
                             extended_max, ratio_cap);
    }

    // Build the P3-V.2 table from the same production K[j,k] and nonuniform
    // u-grid geometry used by advect_u_parallel().  For an interior cell j,
    // w_minus = u_j-u_{j-1} and w_plus = u_{j+1}-u_j; these dimensionless
    // center distances are the local finite-volume geometry weights.  At the
    // two velocity-domain endpoints the projection is one-sided and never
    // fabricates a periodic neighbour.
    void build_energy_conjugate_velocity_table()
    {
        const size_t nv = upar_cells.size();
        const size_t nmu = uperp_cells.size();
        const double speed_tol = 4096.0 * std::numeric_limits<double>::epsilon();
        vx_energy_conjugate_face.assign(
            (nv + 1) * nmu, std::numeric_limits<double>::quiet_NaN());
        vx_energy_conjugate_face_valid.assign((nv + 1) * nmu, 0);
        vx_energy_conjugate_cell.assign(nv * nmu, 0.0);

        for (size_t j = 1; j < nv; ++j) {
            const double du = upar_cells[j] - upar_cells[j - 1];
            if (!(du > 0.0) || !std::isfinite(du))
                throw std::runtime_error("invalid u_parallel center distance for energy-conjugate velocity");
            for (size_t k = 0; k < nmu; ++k) {
                const size_t left = (j - 1) * nmu + k;
                const size_t right = j * nmu + k;
                const double value = (kinetic_energy[right] - kinetic_energy[left]) /
                    (Const::me * Const::c * du);
                const size_t face = j * nmu + k;
                if (!std::isfinite(value) ||
                    std::fabs(value) > Const::c * (1.0 + speed_tol)) {
                    throw std::runtime_error("invalid energy-conjugate u face velocity");
                }
                vx_energy_conjugate_face[face] = value;
                vx_energy_conjugate_face_valid[face] = 1;
            }
        }

        for (size_t j = 0; j < nv; ++j) {
            for (size_t k = 0; k < nmu; ++k) {
                const size_t cell = j * nmu + k;
                double value = 0.0;
                if (j == 0) {
                    value = vx_energy_conjugate_face[nmu + k];
                } else if (j + 1 == nv) {
                    value = vx_energy_conjugate_face[j * nmu + k];
                } else {
                    const double w_minus = upar_center_distances[j];
                    const double w_plus = upar_center_distances[j + 1];
                    if (!(w_minus > 0.0) || !(w_plus > 0.0) ||
                        !std::isfinite(w_minus) || !std::isfinite(w_plus)) {
                        throw std::runtime_error("invalid energy-conjugate cell projection weight");
                    }
                    value = (w_minus * vx_energy_conjugate_face[j * nmu + k] +
                             w_plus * vx_energy_conjugate_face[(j + 1) * nmu + k]) /
                            (w_minus + w_plus);
                }
                if (!std::isfinite(value) ||
                    std::fabs(value) > Const::c * (1.0 + speed_tol)) {
                    throw std::runtime_error("invalid energy-conjugate cell velocity");
                }
                vx_energy_conjugate_cell[cell] = value;
            }
        }

        // Symmetric u grids must produce an antisymmetric conjugate velocity
        // table.  This is an initialization invariant, not a runtime repair.
        double symmetry_error = 0.0;
        double low_speed_error = 0.0;
        for (size_t j = 0; j < nv; ++j) {
            const size_t mirror = nv - 1 - j;
            for (size_t k = 0; k < nmu; ++k) {
                symmetry_error = std::max(
                    symmetry_error,
                    std::fabs(vx_energy_conjugate_cell[j * nmu + k] +
                              vx_energy_conjugate_cell[mirror * nmu + k]));
                if (std::fabs(upar_cells[j]) < 0.05)
                    low_speed_error = std::max(
                        low_speed_error,
                        std::fabs(vx_energy_conjugate_cell[j * nmu + k] -
                                  vx[j * nmu + k]));
            }
        }
        if (symmetry_error > 4096.0 * std::numeric_limits<double>::epsilon() * Const::c)
            throw std::runtime_error("energy-conjugate velocity table symmetry validation failed");
        if (low_speed_error > 0.01 * Const::c)
            throw std::runtime_error("energy-conjugate velocity low-speed limit validation failed");
    }

    void build_upar_core_faces(int nv_core, int nv_tail, double core_max)
    {
        const int offset = nv_tail;
        const double sinh_upar = std::sinh(Param::momentum_upar_stretch);
        for (int j = 0; j <= nv_core; ++j) {
            // Keep this expression structurally identical to the pre-tail
            // grid so tail=0 preserves the established core bit pattern.
            const double xi = -1.0 + 2.0 * static_cast<double>(j) /
                                        nv_core;
            upar_faces[static_cast<size_t>(offset) + j] = core_max *
                std::sinh(Param::momentum_upar_stretch * xi) / sinh_upar;
        }
        upar_faces[static_cast<size_t>(offset)] = -core_max;
        upar_faces[static_cast<size_t>(offset) + nv_core / 2] = 0.0;
        upar_faces[static_cast<size_t>(offset) + nv_core] = core_max;
    }

    void append_symmetric_upar_tail_faces(int nv_core, int nv_tail,
                                          double core_max, double extended_max)
    {
        const int tail = nv_tail;
        if (tail == 0) return;
        const double width = (extended_max - core_max) /
                             static_cast<double>(tail);
        for (int t = 0; t <= tail; ++t) {
            upar_faces[static_cast<size_t>(t)] =
                -extended_max + width * static_cast<double>(t);
            upar_faces[static_cast<size_t>(tail) + nv_core + t] =
                core_max + width * static_cast<double>(t);
        }
        upar_faces.front() = -extended_max;
        upar_faces.back() = extended_max;
    }

    void build_uperp_faces(double umax, double sinh_uperp, int nmu)
    {
        for (int k = 0; k <= nmu; ++k) {
            const double eta = static_cast<double>(k) / nmu;
            uperp_faces[static_cast<size_t>(k)] = umax *
                std::sinh(Param::momentum_uperp_stretch * eta) / sinh_uperp;
        }
        uperp_faces.front() = 0.0;
        uperp_faces.back() = umax;
    }

    void build_cell_geometry_and_moments(int nv, int nmu)
    {
        for (int j = 0; j < nv; ++j) {
            upar_cells[static_cast<size_t>(j)] =
                0.5 * (upar_faces[j] + upar_faces[static_cast<size_t>(j) + 1]);
            upar_widths[static_cast<size_t>(j)] =
                upar_faces[static_cast<size_t>(j) + 1] - upar_faces[j];
        }
        for (int jf = 1; jf < nv; ++jf)
            upar_center_distances[static_cast<size_t>(jf)] =
                upar_cells[jf] - upar_cells[static_cast<size_t>(jf) - 1];
        for (int k = 0; k < nmu; ++k) {
            const double lo = uperp_faces[k];
            const double hi = uperp_faces[static_cast<size_t>(k) + 1];
            uperp_cells[static_cast<size_t>(k)] = 0.5 * (lo + hi);
            uperp_widths[static_cast<size_t>(k)] = hi - lo;
            uperp_ring_areas[static_cast<size_t>(k)] =
                Const::pi * (hi * hi - lo * lo);
        }
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                const size_t slot = static_cast<size_t>(j) * nmu + k;
                const double gamma = std::sqrt(1.0 +
                    upar_cells[j] * upar_cells[j] +
                    uperp_cells[k] * uperp_cells[k]);
                kinetic_energy[slot] = Const::me * Const::c * Const::c *
                                       (gamma - 1.0);
                vx[slot] = Const::c * upar_cells[j] / gamma;
            }
        }
    }

    void validate_nested_grid(double uperp_max, int nv, int nmu, int nv_core,
                              int nv_tail, double core_max,
                              double extended_max, double ratio_cap)
    {
        const int offset = nv_tail;
        const double sinh_upar = std::sinh(Param::momentum_upar_stretch);
        nested_audit.core_face_identity_linf = 0.0;
        nested_audit.core_cell_identity_linf = 0.0;
        nested_audit.core_vx_identity_linf = 0.0;
        nested_audit.core_kinetic_energy_identity_linf = 0.0;
        nested_audit.symmetry_linf = 0.0;
        nested_audit.max_adjacent_width_ratio = 1.0;
        for (int j = 0; j <= nv_core; ++j) {
            const double xi = -1.0 + 2.0 * static_cast<double>(j) /
                                        nv_core;
            const double expected = core_max *
                std::sinh(Param::momentum_upar_stretch * xi) / sinh_upar;
            nested_audit.core_face_identity_linf = std::max(
                nested_audit.core_face_identity_linf,
                std::fabs(upar_faces[static_cast<size_t>(offset) + j] - expected));
        }
        for (int j = 0; j < nv_core; ++j) {
            const double xi_left = -1.0 + 2.0 * static_cast<double>(j) /
                                               nv_core;
            const double xi_right = -1.0 + 2.0 * static_cast<double>(j + 1) /
                                                nv_core;
            const double face_left = core_max *
                std::sinh(Param::momentum_upar_stretch * xi_left) / sinh_upar;
            const double face_right = core_max *
                std::sinh(Param::momentum_upar_stretch * xi_right) / sinh_upar;
            const double expected = 0.5 * (face_left + face_right);
            nested_audit.core_cell_identity_linf = std::max(
                nested_audit.core_cell_identity_linf,
                std::fabs(upar_cells[static_cast<size_t>(offset) + j] - expected));
            for (int k = 0; k < nmu; ++k) {
                const double uperp = uperp_cells[k];
                const double gamma = std::sqrt(1.0 + expected * expected +
                                               uperp * uperp);
                const size_t slot =
                    (static_cast<size_t>(offset) + j) * nmu + k;
                nested_audit.core_vx_identity_linf = std::max(
                    nested_audit.core_vx_identity_linf,
                    std::fabs(vx[slot] - Const::c * expected / gamma));
                nested_audit.core_kinetic_energy_identity_linf = std::max(
                    nested_audit.core_kinetic_energy_identity_linf,
                    std::fabs(kinetic_energy[slot] - Const::me * Const::c * Const::c *
                              (gamma - 1.0)));
            }
        }
        for (int j = 0; j <= nv; ++j) {
            nested_audit.symmetry_linf = std::max(
                nested_audit.symmetry_linf,
                std::fabs(upar_faces[j] + upar_faces[nv - j]));
        }
        for (int j = 1; j < nv; ++j) {
            const double left = upar_widths[static_cast<size_t>(j) - 1];
            const double right = upar_widths[static_cast<size_t>(j)];
            nested_audit.max_adjacent_width_ratio = std::max(
                nested_audit.max_adjacent_width_ratio,
                std::max(left / right, right / left));
        }
        double ring_sum = 0.0;
        for (int k = 0; k < nmu; ++k) ring_sum += uperp_ring_areas[k];
        double upar_sum = 0.0;
        for (int j = 0; j < nv; ++j) upar_sum += upar_widths[j];
        const double expected_phase_volume =
            2.0 * extended_max * Const::pi * uperp_max * uperp_max;
        nested_audit.phase_volume_relative_error = std::fabs(
            upar_sum * ring_sum - expected_phase_volume) /
            std::max(expected_phase_volume, 1.0);

        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                                 std::max(1.0, extended_max);
        if (nested_audit.core_face_identity_linf > tolerance ||
            nested_audit.core_cell_identity_linf > tolerance ||
            nested_audit.core_vx_identity_linf > tolerance * Const::c ||
            nested_audit.core_kinetic_energy_identity_linf >
                tolerance * Const::me * Const::c * Const::c ||
            nested_audit.symmetry_linf > tolerance ||
            nested_audit.max_adjacent_width_ratio > ratio_cap + tolerance ||
            nested_audit.phase_volume_relative_error > 1.0e-13 ||
            upar_faces.front() != -extended_max ||
            upar_faces.back() != extended_max) {
            std::ostringstream message;
            message << "nested u_parallel grid validation failed: "
                    << "Nv_core=" << nv_core
                    << " Nv_tail_per_side=" << nv_tail
                    << " core_max=" << core_max
                    << " extended_max=" << extended_max
                    << " core_face_linf=" << nested_audit.core_face_identity_linf
                    << " core_cell_linf=" << nested_audit.core_cell_identity_linf
                    << " symmetry_linf=" << nested_audit.symmetry_linf
                    << " max_adjacent_width_ratio="
                    << nested_audit.max_adjacent_width_ratio
                    << " phase_volume_relative_error="
                    << nested_audit.phase_volume_relative_error
                    << " (required max_adjacent_width_ratio<=" << ratio_cap
                    << ")";
            throw std::runtime_error(message.str());
        }

#if FP_ENABLE_DEBUG_DIAGNOSTICS
        for (int j = 0; j < nv; ++j) {
            assert(std::isfinite(upar_faces[j]));
            assert(upar_faces[j + 1] > upar_faces[j]);
            assert(std::isfinite(upar_widths[j]) && upar_widths[j] > 0.0);
        }
        for (int k = 0; k < nmu; ++k) {
            assert(std::isfinite(uperp_faces[k]));
            assert(uperp_faces[k + 1] > uperp_faces[k]);
            assert(std::isfinite(uperp_widths[k]) && uperp_widths[k] > 0.0);
        }
        assert(upar_faces.front() == -extended_max &&
               upar_faces.back() == extended_max);
        assert(uperp_faces.front() == 0.0 && uperp_faces.back() == uperp_max);
#ifndef NDEBUG
        const double symmetry_tolerance = 64.0 *
            std::numeric_limits<double>::epsilon() * extended_max;
        for (int j = 0; j <= nv; ++j)
            assert(std::fabs(upar_faces[j] + upar_faces[nv - j]) <=
                   symmetry_tolerance);
#endif
        double debug_ring_sum = 0.0;
        for (int k = 0; k < nmu; ++k) debug_ring_sum += uperp_ring_areas[k];
        assert(std::fabs(debug_ring_sum - Const::pi * uperp_max * uperp_max) /
               (Const::pi * uperp_max * uperp_max) <= 1.0e-13);
        if (nv == 192)
            assert(*std::min_element(upar_widths.begin(), upar_widths.end()) >= 0.004 &&
                   *std::min_element(upar_widths.begin(), upar_widths.end()) <= 0.007);
        if (nmu == 64)
            assert(*std::min_element(uperp_widths.begin(), uperp_widths.end()) >= 0.0035 &&
                   *std::min_element(uperp_widths.begin(), uperp_widths.end()) <= 0.0070);
#endif
    }

    double cell_phase_volume(int j, int k) const
    {
        return upar_widths[j] * uperp_ring_areas[k];
    }

    // Grid-quadrature offset of the discrete Maxwellian second moment
    // (section 19.3 energy-conservation fix for the moment-closure
    // collision): the discrete Maxwellian with parameter theta has moments
    // u_th^2 = theta*(1+eps(theta)); to leading order eps ~ 1/theta, so the
    // offset delta = theta0 - u_th^2(M_theta0) is independent of theta and
    // the closure temperature that makes the equilibrium moments
    // self-consistent is theta_used = u_th^2_moments + delta.  Returns
    // delta in u^2 units for the reference temperature theta0.
    // Build the moment-closure self-consistency table (section 19.3
    // energy-conservation fix): F(theta) = u_th^2 of the discrete
    // Maxwellian with parameter theta, sampled on log-spaced points in
    // [theta_ref/64, 64*theta_ref].  The closure temperature must satisfy
    // F(theta_used) = u_th^2(current moments); the table inverts this map
    // per cell (see moment_closure_uth2_self_consistent).  theta_ref is
    // the reference temperature in u^2 units; the production grid uses
    // T_e/(m_e c^2).
    void build_moment_closure_table(double theta_ref)
    {
        const size_t nmu = uperp_cells.size();
        const int n_points = 1024;
        mc_uth2_table_.assign(static_cast<size_t>(n_points), 0.0);
        mc_moment_table_.assign(static_cast<size_t>(n_points), 0.0);
        const double log_lo = std::log(theta_ref / 64.0);
        const double log_hi = std::log(theta_ref * 64.0);
        for (int i = 0; i < n_points; ++i) {
            const double theta = std::exp(
                log_lo + (log_hi - log_lo) * static_cast<double>(i) /
                             static_cast<double>(n_points - 1));
            double n = 0.0;
            double ke = 0.0;
            const double inv2theta = 1.0 / (2.0 * theta);
            for (size_t j = 0; j < upar_cells.size(); ++j) {
                for (size_t k = 0; k < nmu; ++k) {
                    const double u2 = upar_cells[j] * upar_cells[j] +
                                      uperp_cells[k] * uperp_cells[k];
                    const double w = std::exp(-u2 * inv2theta) *
                                     cell_phase_volume(
                                         static_cast<int>(j),
                                         static_cast<int>(k));
                    n += w;
                    ke += w * kinetic_energy[j * nmu + k];
                }
            }
            mc_uth2_table_[static_cast<size_t>(i)] = theta;
            mc_moment_table_[static_cast<size_t>(i)] =
                (n > 0.0) ? (2.0 / 3.0) * (ke / n) /
                                (Const::me * Const::c * Const::c)
                          : theta;
        }
    }

    // Invert the moment-closure table: return the closure temperature
    // theta_used with F(theta_used) = uth2_mom (cubic Hermite
    // interpolation; linear at the table ends, clamped outside the range).
    double moment_closure_uth2_self_consistent(double uth2_mom) const
    {
        const size_t n = mc_moment_table_.size();
        if (n < 2) return uth2_mom;
        if (uth2_mom <= mc_moment_table_.front())
            return mc_uth2_table_.front();
        if (uth2_mom >= mc_moment_table_.back())
            return mc_uth2_table_.back();
        size_t i = 1;
        while (i + 1 < n && mc_moment_table_[i] < uth2_mom) ++i;
        const size_t il = i - 1;
        const size_t ir = i;
        const double fl = mc_moment_table_[il];
        const double fr = mc_moment_table_[ir];
        const double h = fr - fl;
        if (!(h > 0.0)) return mc_uth2_table_[il];
        const double s = (uth2_mom - fl) / h;
        if (il == 0 || ir + 1 >= n) {
            // Linear fallback at the ends.
            return mc_uth2_table_[il] + s *
                   (mc_uth2_table_[ir] - mc_uth2_table_[il]);
        }
        const double m0 = (mc_uth2_table_[ir] - mc_uth2_table_[il - 1]) /
                          (fr - mc_moment_table_[il - 1]);
        const double m1 = (mc_uth2_table_[ir + 1] - mc_uth2_table_[il]) /
                          (mc_moment_table_[ir + 1] - fl);
        const double h00 = 2.0 * s * s * s - 3.0 * s * s + 1.0;
        const double h10 = s * s * s - 2.0 * s * s + s;
        const double h01 = -2.0 * s * s * s + 3.0 * s * s;
        const double h11 = s * s * s - s * s;
        return h00 * mc_uth2_table_[il] + h10 * h * m0 +
               h01 * mc_uth2_table_[ir] + h11 * h * m1;
    }

    bool is_uniform() const
    {
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
        const double upar_reference = upar_widths.empty() ? 0.0 : upar_widths[0];
        const double uperp_reference = uperp_widths.empty() ? 0.0 : uperp_widths[0];
        for (size_t j = 1; j < upar_widths.size(); ++j)
            if (std::fabs(upar_widths[j] - upar_reference) >
                tolerance * std::max(1.0, std::fabs(upar_reference))) return false;
        for (size_t k = 1; k < uperp_widths.size(); ++k)
            if (std::fabs(uperp_widths[k] - uperp_reference) >
                tolerance * std::max(1.0, std::fabs(uperp_reference))) return false;
        return true;
    }
};

struct SpatialGrid {
    int nx_global;
    int nx_local;
    int ix_start;
    int nghost;
    int nx_total;
    double dx;
    double x_min;

    void init(int rank, int nranks) {
        init_with_domain(rank, nranks, Param::nx, Param::Lx);
    }

    // Open VPFP uses a single explicit-domain initializer; there is no
    // separate legacy periodic-domain path.
    void init_with_domain(int rank, int nranks,
                          int requested_nx_global, double length) {
        if (nranks <= 0 || rank < 0 || rank >= nranks ||
            requested_nx_global <= 0 || !std::isfinite(length) ||
            length <= 0.0) {
            throw std::runtime_error("invalid spatial domain decomposition");
        }
        if (requested_nx_global < nranks) {
            throw std::runtime_error("spatial grid has fewer cells than MPI ranks");
        }
        nx_global = requested_nx_global;
        dx = length / static_cast<double>(nx_global);
        x_min = 0.0;
        nghost = Param::Nghost;

        nx_local = nx_global / nranks;
        int remainder = nx_global % nranks;
        if (rank < remainder) {
            nx_local += 1;
            ix_start = rank * nx_local;
        } else {
            ix_start = remainder * (nx_local + 1) + (rank - remainder) * nx_local;
        }
        if (nx_local <= 0) {
            throw std::runtime_error("spatial partition leaves an empty local grid");
        }
        nx_total = nx_local + 2 * nghost;
    }

    bool owns_left_physical_boundary(int rank) const {
        return rank == 0;
    }

    bool owns_right_physical_boundary(int rank, int size) const {
        return rank == size - 1;
    }

    double length() const {
        return dx * static_cast<double>(nx_global);
    }

    double left_boundary() const {
        return 0.0;
    }

    double right_boundary() const {
        return length();
    }

    int global_cell(int local_with_ghost) const {
        return ix_start + (local_with_ghost - nghost);
    }

    double x(int i_local) const {
        int ig = ix_start + (i_local - nghost);
        return x_min + (ig + 0.5) * dx;
    }
};

// Forward declaration: the single shared (j,k) slot helper is defined with
// the other index helpers at the end of this header.
inline size_t idx2(int iv, int imu);

// Shared cell-moment formula for one cell-integrated cylindrical mass M at
// (u_parallel, u_perp) (sections 4.1 and 7.6).  This is the single production
// formula used both by Species::extract_conversion_masses and by the
// BulkTailConverter when it builds the reference loading.
inline void mass_cell_moments(double mass, double upar, double uperp,
                              double& number, double& px, double& energy,
                              double& jx_dx, double& pixx_dx,
                              double& piperp_dx)
{
    const double gamma = std::sqrt(1.0 + upar * upar + uperp * uperp);
    number = mass;                                     // m^-2
    px = Const::me * Const::c * mass * upar;           // kg m s^-1 m^-2
    energy = Const::me * Const::c * Const::c * mass *
             (gamma - 1.0);                            // J m^-2
    jx_dx = -Const::qe * Const::c * mass * upar / gamma;      // A m^-1
    pixx_dx = Const::me * Const::c * Const::c * mass *
              upar * upar / gamma;                     // J m^-2
    piperp_dx = Const::me * Const::c * Const::c * mass *
                uperp * uperp / gamma;                 // J m^-2
}

enum class VelocityFaceDirection { U_PARALLEL, U_PERP };

struct BulkTailInterfaceFace {
    VelocityFaceDirection direction;
    int face_index;
    int transverse_index;
    int bulk_iv;
    int bulk_imu;
    int tail_iv;
    int tail_imu;
    int outward_sign;
};

// Hybrid bulk/tail velocity partition (sections 5.2, 5.4 and 14.2): classifies
// every cylindrical velocity cell by its cell-centre relativistic kinetic
// energy relative to the tail conversion threshold K_out.  The arrays are
// computed once in init() and are read-only afterwards.  First-version
// classification is by the full cell at its centre K_jk (section 5.4).
struct HybridVelocityPartition {
    std::vector<double> kinetic_energy;          // J per (j,k) slot
    std::vector<unsigned char> is_conversion_cell;
    std::vector<unsigned char> is_bulk_resolved;
    std::vector<unsigned char> is_buffer_cell;
    // Explicit face-aligned topology used by the flux-interface path.  The
    // legacy is_conversion_cell mask remains unchanged for static-cell A/B.
    std::vector<unsigned char> tail_owned_cell;
    std::vector<unsigned char> bulk_owned_cell;
    std::vector<BulkTailInterfaceFace> upar_interface_faces;
    std::vector<BulkTailInterfaceFace> uperp_interface_faces;
    std::vector<int> upar_interface_lookup;
    std::vector<int> uperp_interface_lookup;
    int upar_count;
    int uperp_count;
    double convert_energy_mev;
    double buffer_width_mev;
    int upar_bins;
    int energy_bins;
    double min_conversion_energy;                // J = K_out
    double max_conversion_energy;                // J over conversion cells
    double upar_max;
    std::uint64_t config_hash;
    // Section 7.11.4 branch B: explicit threshold-aware energy edges
    // (Joules), read-only after init.  The boundary list explicitly
    // contains K_out, K_out+0.2, K_out+0.4 and K_out+0.8 MeV and then
    // widens logarithmically (doubling) up to max_conversion_energy.  The
    // THRESHOLD_AWARE_COMPRESSION loading policy bins conversion cells with
    // energy_bin_threshold_aware(); the production default policy keeps the
    // legacy uniform energy_bin() until the section 7.11 acceptance passes.
    std::vector<double> conversion_energy_edges;

    void init(const CylindricalVelocityGrid& cgrid,
              double convert_energy_mev_value,
              double buffer_width_mev_value,
              int upar_bins_value, int energy_bins_value)
    {
        convert_energy_mev = convert_energy_mev_value;
        buffer_width_mev = buffer_width_mev_value;
        upar_bins = std::max(1, upar_bins_value);
        energy_bins = std::max(1, energy_bins_value);
        if (!(convert_energy_mev > 0.0) || !(buffer_width_mev > 0.0)) {
            throw std::runtime_error(
                "invalid tail conversion partition configuration");
        }
        const size_t nslots = cgrid.kinetic_energy.size();
        if (nslots == 0 || cgrid.upar_faces.size() < 2) {
            throw std::runtime_error(
                "tail conversion partition requires a built velocity grid");
        }
        kinetic_energy = cgrid.kinetic_energy;
        upar_count = static_cast<int>(cgrid.upar_cells.size());
        uperp_count = static_cast<int>(cgrid.uperp_cells.size());
        is_conversion_cell.assign(nslots, 0);
        is_bulk_resolved.assign(nslots, 0);
        is_buffer_cell.assign(nslots, 0);
        tail_owned_cell.assign(nslots, 0);
        bulk_owned_cell.assign(nslots, 1);

        const double k_out = convert_energy_mev * 1.0e6 * Const::eV;
        const double buffer_energy = buffer_width_mev * 1.0e6 * Const::eV;
        min_conversion_energy = k_out;
        max_conversion_energy = k_out;
        upar_max = cgrid.upar_faces.back();

        double max_grid_energy = 0.0;
        for (size_t s = 0; s < nslots; ++s) {
            max_grid_energy = std::max(max_grid_energy, kinetic_energy[s]);
        }
        if (!(k_out < max_grid_energy)) {
            std::ostringstream message;
            message << "tail conversion threshold " << convert_energy_mev
                    << " MeV lies beyond the velocity grid max energy "
                    << max_grid_energy / (1.0e6 * Const::eV) << " MeV";
            throw std::runtime_error(message.str());
        }

        size_t conversion_count = 0;
        for (int j = 0; j < upar_count; ++j) {
            for (int k = 0; k < uperp_count; ++k) {
                const size_t slot = static_cast<size_t>(j) *
                                    static_cast<size_t>(uperp_count) +
                                    static_cast<size_t>(k);
                const double ke = kinetic_energy[slot];
                if (ke < k_out) {
                    is_bulk_resolved[slot] = 1;
                    continue;
                }
                is_conversion_cell[slot] = 1;
                ++conversion_count;
                max_conversion_energy =
                    std::max(max_conversion_energy, ke);
                if (ke < k_out + buffer_energy) {
                    is_buffer_cell[slot] = 1;
                }
            }
        }
        if (conversion_count == 0) {
            throw std::runtime_error(
                "tail conversion partition has no conversion cells");
        }

        // Flux-interface ownership is stricter than the legacy center-cell
        // conversion mask: every fixed positive quadrature node in a cell
        // must be above K_out.  This avoids exporting mass whose physical
        // support still intersects the bulk side of the threshold.
        static const double gl4[4] = {
            -0.8611363115940526, -0.3399810435848563,
             0.3399810435848563,  0.8611363115940526
        };
        const int nv = upar_count;
        const int nmu = uperp_count;
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                bool all_tail = true;
                const double up0 = cgrid.upar_faces[j];
                const double up1 = cgrid.upar_faces[j + 1];
                const double ut0 = cgrid.uperp_faces[k];
                const double ut1 = cgrid.uperp_faces[k + 1];
                for (int a = 0; a < 4 && all_tail; ++a) {
                    const double up = 0.5 * (up0 + up1) +
                                      0.5 * (up1 - up0) * gl4[a];
                    for (int b = 0; b < 4; ++b) {
                        const double ut = 0.5 * (ut0 + ut1) +
                                          0.5 * (ut1 - ut0) * gl4[b];
                        const double gamma = std::sqrt(1.0 + up * up + ut * ut);
                        const double node_ke = Const::me * Const::c * Const::c *
                                               (gamma - 1.0);
                        if (!(node_ke >= k_out)) {
                            all_tail = false;
                            break;
                        }
                    }
                }
                const size_t slot = static_cast<size_t>(j) *
                                    static_cast<size_t>(nmu) +
                                    static_cast<size_t>(k);
                tail_owned_cell[slot] = all_tail ? 1 : 0;
                bulk_owned_cell[slot] = all_tail ? 0 : 1;
            }
        }

        upar_interface_faces.clear();
        uperp_interface_faces.clear();
        upar_interface_lookup.assign(
            static_cast<size_t>(nv + 1) * static_cast<size_t>(nmu), -1);
        uperp_interface_lookup.assign(
            static_cast<size_t>(nv) * static_cast<size_t>(nmu + 1), -1);
        for (int face = 1; face < nv; ++face) {
            for (int k = 0; k < nmu; ++k) {
                const bool left_tail = tail_owned_cell[
                    static_cast<size_t>(face - 1) * nmu + k] != 0;
                const bool right_tail = tail_owned_cell[
                    static_cast<size_t>(face) * nmu + k] != 0;
                if (left_tail == right_tail) continue;
                BulkTailInterfaceFace item;
                item.direction = VelocityFaceDirection::U_PARALLEL;
                item.face_index = face;
                item.transverse_index = k;
                item.bulk_iv = left_tail ? face : face - 1;
                item.bulk_imu = k;
                item.tail_iv = left_tail ? face - 1 : face;
                item.tail_imu = k;
                item.outward_sign = left_tail ? -1 : 1;
                const int index = static_cast<int>(upar_interface_faces.size());
                upar_interface_faces.push_back(item);
                upar_interface_lookup[static_cast<size_t>(face) * nmu + k] = index;
            }
        }
        // The u_perp axis is never an interface.  For every other face, only
        // a bulk/tail change is represented; each shared face has one entry.
        for (int j = 0; j < nv; ++j) {
            for (int face = 1; face < nmu; ++face) {
                const bool lo_tail = tail_owned_cell[
                    static_cast<size_t>(j) * nmu + face - 1] != 0;
                const bool hi_tail = tail_owned_cell[
                    static_cast<size_t>(j) * nmu + face] != 0;
                if (lo_tail == hi_tail) continue;
                BulkTailInterfaceFace item;
                item.direction = VelocityFaceDirection::U_PERP;
                item.face_index = face;
                item.transverse_index = j;
                item.bulk_iv = j;
                item.bulk_imu = lo_tail ? face : face - 1;
                item.tail_iv = j;
                item.tail_imu = lo_tail ? face - 1 : face;
                item.outward_sign = lo_tail ? -1 : 1;
                const int index = static_cast<int>(uperp_interface_faces.size());
                uperp_interface_faces.push_back(item);
                uperp_interface_lookup[static_cast<size_t>(j) * (nmu + 1) + face] = index;
            }
        }
        // Build the explicit threshold-aware edges (section 7.11.4 branch
        // B): the 0.2 MeV structure near K_out is fixed, higher energies
        // widen logarithmically.
        conversion_energy_edges.clear();
        const double e_step = 0.2e6 * Const::eV;
        conversion_energy_edges.push_back(k_out);
        conversion_energy_edges.push_back(k_out + e_step);
        conversion_energy_edges.push_back(k_out + 2.0 * e_step);
        conversion_energy_edges.push_back(k_out + 4.0 * e_step);
        double edge = k_out + 4.0 * e_step;
        while (edge < max_conversion_energy) {
            edge *= 2.0;
            conversion_energy_edges.push_back(edge);
        }
        if (conversion_energy_edges.back() < max_conversion_energy) {
            conversion_energy_edges.push_back(max_conversion_energy);
        }

        // Deterministic configuration hash (FNV-1a over the partition inputs
        // and the velocity-grid identity).
        config_hash = 0xcbf29ce484222325ULL;
        const std::uint64_t prime = 0x100000001b3ULL;
        std::uint64_t word;
        std::memcpy(&word, &convert_energy_mev, sizeof(word));
        config_hash ^= word;
        config_hash *= prime;
        std::memcpy(&word, &buffer_width_mev, sizeof(word));
        config_hash ^= word;
        config_hash *= prime;
        config_hash ^= static_cast<std::uint64_t>(upar_bins);
        config_hash *= prime;
        config_hash ^= static_cast<std::uint64_t>(energy_bins);
        config_hash *= prime;
        config_hash ^= static_cast<std::uint64_t>(upar_count);
        config_hash *= prime;
        config_hash ^= static_cast<std::uint64_t>(uperp_count);
        config_hash *= prime;
        std::memcpy(&word, &upar_max, sizeof(word));
        config_hash ^= word;
        config_hash *= prime;
        // Section 7.11.4 branch B: fold the threshold-aware energy edges
        // into the configuration hash so a checkpoint written with a
        // different grouping cannot be silently restarted.
        for (size_t i = 0; i < conversion_energy_edges.size(); ++i) {
            std::memcpy(&word, &conversion_energy_edges[i], sizeof(word));
            config_hash ^= word;
            config_hash *= prime;
        }
        config_hash ^= static_cast<std::uint64_t>(upar_interface_faces.size());
        config_hash *= prime;
        config_hash ^= static_cast<std::uint64_t>(uperp_interface_faces.size());
        config_hash *= prime;
        for (size_t i = 0; i < upar_interface_faces.size(); ++i) {
            const BulkTailInterfaceFace& f = upar_interface_faces[i];
            config_hash ^= static_cast<std::uint64_t>(f.face_index + 1);
            config_hash *= prime;
            config_hash ^= static_cast<std::uint64_t>(f.transverse_index + 1);
            config_hash *= prime;
            config_hash ^= static_cast<std::uint64_t>(f.outward_sign + 2);
            config_hash *= prime;
        }
        for (size_t i = 0; i < uperp_interface_faces.size(); ++i) {
            const BulkTailInterfaceFace& f = uperp_interface_faces[i];
            config_hash ^= static_cast<std::uint64_t>(f.face_index + 1);
            config_hash *= prime;
            config_hash ^= static_cast<std::uint64_t>(f.transverse_index + 1);
            config_hash *= prime;
            config_hash ^= static_cast<std::uint64_t>(f.outward_sign + 2);
            config_hash *= prime;
        }
        config_hash ^= topology_mask_hash();
        config_hash *= prime;
    }

    bool is_conversion(int j, int k) const
    {
        return is_conversion_cell[slot_index(j, k)] != 0;
    }

    bool is_tail_owned(int j, int k) const
    {
        return tail_owned_cell[slot_index(j, k)] != 0;
    }

    bool is_bulk_owned(int j, int k) const
    {
        return bulk_owned_cell[slot_index(j, k)] != 0;
    }

    int upar_interface_index(int face, int imu) const
    {
        if (face < 0 || face > upar_count ||
            imu < 0 || imu >= uperp_count) return -1;
        return upar_interface_lookup[static_cast<size_t>(face) * uperp_count + imu];
    }

    int uperp_interface_index(int iv, int face) const
    {
        const int nmu = uperp_count;
        const int nv = upar_count;
        if (iv < 0 || iv >= nv || face <= 0 || face >= nmu) return -1;
        return uperp_interface_lookup[static_cast<size_t>(iv) * (nmu + 1) + face];
    }

    size_t upar_cells_count() const { return static_cast<size_t>(upar_count); }
    size_t uperp_cells_count() const { return static_cast<size_t>(uperp_count); }

    // |u_parallel| magnitude bin over [0, upar_max].
    int upar_bin(double upar) const
    {
        if (upar_bins <= 1) return 0;
        const double a = std::fabs(upar);
        const double span = std::max(upar_max, 1.0e-300);
        const int b = static_cast<int>(a / span *
                                       static_cast<double>(upar_bins));
        return std::max(0, std::min(upar_bins - 1, b));
    }

    // Energy bin over the explicit threshold-aware conversion_energy_edges
    // (section 7.11.4 branch B): upper_bound over the edges.  The bins near
    // K_out are 0.2 MeV wide (K_out/+0.2/+0.4/+0.8) and widen
    // logarithmically above, so no compression group may span a threshold
    // audit bin.  The legacy uniform binning is retained as
    // energy_bin_uniform() for the CURRENT_PRODUCTION_COMPRESSION reference
    // policy in the section 7.11.3 A/B/C tests.
    int energy_bin(double kinetic_energy_j) const
    {
        const size_t n = conversion_energy_edges.size();
        if (n < 2) return 0;
        if (kinetic_energy_j <= conversion_energy_edges.front()) return 0;
        const size_t b = static_cast<size_t>(std::upper_bound(
            conversion_energy_edges.begin(), conversion_energy_edges.end(),
            kinetic_energy_j) - conversion_energy_edges.begin());
        return std::max(0, static_cast<int>(b) - 1);
    }

    // Legacy uniform energy bin over [K_out, max_conversion_energy],
    // retained for the CURRENT_PRODUCTION_COMPRESSION reference policy in
    // the section 7.11.3 A/B/C tests (represents the pre-fix grouping).
    int energy_bin_uniform(double kinetic_energy_j) const
    {
        if (energy_bins <= 1) return 0;
        const double span =
            std::max(max_conversion_energy - min_conversion_energy, 1.0e-300);
        const int b = static_cast<int>(
            (kinetic_energy_j - min_conversion_energy) / span *
            static_cast<double>(energy_bins));
        return std::max(0, std::min(energy_bins - 1, b));
    }

    // Alias used by the THRESHOLD_AWARE_COMPRESSION loading policy; equals
    // the production energy_bin().
    int energy_bin_threshold_aware(double kinetic_energy_j) const
    {
        return energy_bin(kinetic_energy_j);
    }

    // Versioned topology components are persisted separately in the
    // checkpoint manifest.  The combined hash below remains the compact
    // restart gate used by the binary checkpoint header.
    std::uint64_t interface_topology_version() const
    {
        return 1ULL;
    }

    std::uint64_t tail_owned_mask_hash() const
    {
        std::uint64_t h = 0xcbf29ce484222325ULL;
        const std::uint64_t prime = 0x100000001b3ULL;
        for (size_t i = 0; i < tail_owned_cell.size(); ++i) {
            h ^= static_cast<std::uint64_t>(tail_owned_cell[i]);
            h *= prime;
        }
        return h;
    }

    std::uint64_t interface_face_list_hash() const
    {
        std::uint64_t h = 0xcbf29ce484222325ULL;
        const std::uint64_t prime = 0x100000001b3ULL;
        for (size_t i = 0; i < upar_interface_faces.size(); ++i) {
            const BulkTailInterfaceFace& f = upar_interface_faces[i];
            h ^= 0ULL; h *= prime;
            h ^= static_cast<std::uint64_t>(f.face_index + 1); h *= prime;
            h ^= static_cast<std::uint64_t>(f.transverse_index + 1); h *= prime;
            h ^= static_cast<std::uint64_t>(f.bulk_iv + 1); h *= prime;
            h ^= static_cast<std::uint64_t>(f.tail_iv + 1); h *= prime;
            h ^= static_cast<std::uint64_t>(f.outward_sign + 2); h *= prime;
        }
        for (size_t i = 0; i < uperp_interface_faces.size(); ++i) {
            const BulkTailInterfaceFace& f = uperp_interface_faces[i];
            h ^= 1ULL; h *= prime;
            h ^= static_cast<std::uint64_t>(f.face_index + 1); h *= prime;
            h ^= static_cast<std::uint64_t>(f.transverse_index + 1); h *= prime;
            h ^= static_cast<std::uint64_t>(f.bulk_imu + 1); h *= prime;
            h ^= static_cast<std::uint64_t>(f.tail_imu + 1); h *= prime;
            h ^= static_cast<std::uint64_t>(f.outward_sign + 2); h *= prime;
        }
        return h;
    }

    std::uint64_t topology_mask_hash() const
    {
        std::uint64_t h = 0xcbf29ce484222325ULL;
        const std::uint64_t prime = 0x100000001b3ULL;
        h ^= interface_topology_version(); h *= prime;
        std::uint64_t word = 0;
        std::memcpy(&word, &min_conversion_energy, sizeof(word));
        h ^= word; h *= prime;
        h ^= tail_owned_mask_hash(); h *= prime;
        h ^= interface_face_list_hash(); h *= prime;
        return h;
    }

    bool flux_interface_topology_valid(std::string* reason = NULL) const
    {
        const int nv = upar_count;
        const int nmu = uperp_count;
        if (nv <= 0 || nmu <= 0) {
            if (reason) *reason = "empty velocity grid";
            return false;
        }
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                if (tail_owned_cell[slot_index(j, k)] !=
                    tail_owned_cell[slot_index(nv - 1 - j, k)]) {
                    if (reason) *reason = "u_parallel mask is not mirrored";
                    return false;
                }
            }
        }
        for (size_t i = 0; i < uperp_interface_faces.size(); ++i) {
            if (uperp_interface_faces[i].face_index <= 0) {
                if (reason) *reason = "u_perp axis is an interface";
                return false;
            }
        }
        for (size_t i = 0; i < upar_interface_lookup.size(); ++i) {
            if (upar_interface_lookup[i] >=
                static_cast<int>(upar_interface_faces.size())) {
                if (reason) *reason = "invalid u_parallel lookup";
                return false;
            }
        }
        for (size_t i = 0; i < uperp_interface_lookup.size(); ++i) {
            if (uperp_interface_lookup[i] >=
                static_cast<int>(uperp_interface_faces.size())) {
                if (reason) *reason = "invalid u_perp lookup";
                return false;
            }
        }
        std::vector<unsigned char> seen(tail_owned_cell.size(), 0);
        std::vector<std::pair<int, int> > stack;
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                if (tail_owned_cell[slot_index(j, k)] &&
                    (j == 0 || j == nv - 1 || k == nmu - 1))
                    stack.push_back(std::make_pair(j, k));
            }
        }
        const int dj[4] = {-1, 1, 0, 0};
        const int dk[4] = {0, 0, -1, 1};
        while (!stack.empty()) {
            const std::pair<int, int> cell = stack.back();
            stack.pop_back();
            const size_t s = slot_index(cell.first, cell.second);
            if (seen[s]) continue;
            seen[s] = 1;
            for (int q = 0; q < 4; ++q) {
                const int jj = cell.first + dj[q];
                const int kk = cell.second + dk[q];
                if (jj >= 0 && jj < nv && kk >= 0 && kk < nmu &&
                    tail_owned_cell[slot_index(jj, kk)])
                    stack.push_back(std::make_pair(jj, kk));
            }
        }
        for (int j = 0; j < nv; ++j) {
            for (int k = 0; k < nmu; ++k) {
                if (tail_owned_cell[slot_index(j, k)] &&
                    !seen[slot_index(j, k)]) {
                    if (reason) *reason = "isolated tail component";
                    return false;
                }
            }
        }
        return true;
    }

    size_t slot_index(int j, int k) const
    {
        return static_cast<size_t>(j) * static_cast<size_t>(uperp_count) +
               static_cast<size_t>(k);
    }

    // Section 7.11.4 branch B: deterministic FNV-1a hash of the explicit
    // threshold-aware energy edges, written to the checkpoint/snapshot
    // manifests and used to reject a restart with a different grouping.
    std::uint64_t conversion_energy_edges_hash() const
    {
        std::uint64_t hash = 0xcbf29ce484222325ULL;
        const std::uint64_t prime = 0x100000001b3ULL;
        for (size_t i = 0; i < conversion_energy_edges.size(); ++i) {
            std::uint64_t word = 0;
            std::memcpy(&word, &conversion_energy_edges[i], sizeof(word));
            hash ^= word;
            hash *= prime;
        }
        return hash;
    }
};

// Beam macro-particle weight must be derived from the runtime SpatialGrid so
// that changing dx keeps the injected physical flux unchanged.  It lives
// next to SpatialGrid because parameters.h cannot include grid.h.
inline double beam_macro_weight(const SpatialGrid& grid) {
    return Param::densb * grid.dx /
           static_cast<double>(Param::beam_macro_particles_per_cell);
}

inline double gamma_from_v(double v) {
    double beta = v / Const::c;
    if (beta > 0.999999999999) beta = 0.999999999999;
    return 1.0 / std::sqrt(1.0 - beta * beta);
}

inline double gamma_from_u(double u) {
    return std::sqrt(1.0 + u * u);
}

inline double speed_from_u(double u) {
    return Const::c * u / gamma_from_u(u);
}

inline double u_from_v(double v) {
    const double beta = std::max(-0.999999999999,
                                 std::min(0.999999999999, v / Const::c));
    return beta / std::sqrt(1.0 - beta * beta);
}

inline double momentum_from_v(double v, double mass) {
    return gamma_from_v(v) * mass * v;
}

inline size_t idx3(int ix, int iv, int imu) {
    return static_cast<size_t>(ix) * Param::Nvmu
         + static_cast<size_t>(iv) * Param::Nmu
         + static_cast<size_t>(imu);
}

inline size_t idx2(int iv, int imu) {
    return static_cast<size_t>(iv) * Param::Nmu
         + static_cast<size_t>(imu);
}

#endif
