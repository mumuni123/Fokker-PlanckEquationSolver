#include "discrete_moment_operators.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

struct GridLevel {
    int nx;
    int nu;
    int nup;
    double dx;
    std::vector<double> upar_faces;
    std::vector<double> upar_cells;
    std::vector<double> dupar;
    std::vector<double> duperp;
    std::vector<double> uperp_faces;
    std::vector<double> uperp_cells;
    std::vector<double> uperp_ring_areas;
    std::vector<double> kinetic_energy;
    std::vector<double> vx;
    double upar_stretch;
    double uperp_stretch;
    bool geometry_valid;
};

struct Fluxes {
    std::vector<double> x;
    std::vector<double> cu;
    std::vector<double> u;
};

struct Result {
    int nx;
    int nu;
    int nup;
    double dt;
    double x_cfl;
    double u_cfl;
    double j_pair_relative;
    double work_pair_relative;
    double g_adjoint_relative;
    double mass_relative;
    double min_mass_mid;
    double min_dupar;
    double max_dupar;
    double min_duperp;
    double max_duperp;
    double upar_stretch;
    double uperp_stretch;
    bool geometry_valid;
};

size_t cell_index(const GridLevel& grid, int ix, int j, int k)
{
    return (static_cast<size_t>(ix) * grid.nu + j) * grid.nup + k;
}

size_t xface_index(const GridLevel& grid, int iface, int j, int k)
{
    return (static_cast<size_t>(iface) * grid.nu + j) * grid.nup + k;
}

size_t uface_index(const GridLevel& grid, int ix, int jf, int k)
{
    return (static_cast<size_t>(ix) * (grid.nu + 1) + jf) * grid.nup + k;
}

size_t vk_index(const GridLevel& grid, int j, int k)
{
    return static_cast<size_t>(j) * grid.nup + k;
}

double relative_error(double residual, double lhs, double rhs)
{
    const double scale = std::max(std::numeric_limits<double>::min(),
        std::max(std::fabs(lhs), std::fabs(rhs)));
    return std::fabs(residual) / scale;
}

GridLevel make_level(int refinement)
{
    // Fixed sinh mapping makes these nested face subsets of the production
    // nonuniform grid.  The controlled test therefore audits production
    // geometry rather than introducing an independent velocity mesh.
    CylindricalVelocityGrid production_grid;
    production_grid.init(Param::momentum_umax);
    const int stride = 4 / refinement;

    GridLevel grid;
    grid.nx = 8 * refinement;
    grid.nu = Param::Nv / stride;
    grid.nup = Param::Nmu / stride;
    grid.dx = Param::Lx / grid.nx;
    grid.upar_faces.resize(grid.nu + 1);
    grid.upar_cells.resize(grid.nu);
    grid.dupar.resize(grid.nu);
    grid.duperp.resize(grid.nup);
    grid.uperp_faces.resize(grid.nup + 1);
    grid.uperp_cells.resize(grid.nup);
    grid.uperp_ring_areas.resize(grid.nup);
    grid.kinetic_energy.resize(static_cast<size_t>(grid.nu) * grid.nup);
    grid.vx.resize(static_cast<size_t>(grid.nu) * grid.nup);
    grid.upar_stretch = Param::momentum_upar_stretch;
    grid.uperp_stretch = Param::momentum_uperp_stretch;
    grid.geometry_valid = true;

    for (int j = 0; j <= grid.nu; ++j)
        grid.upar_faces[j] = production_grid.upar_faces[j * stride];
    for (int k = 0; k <= grid.nup; ++k)
        grid.uperp_faces[k] = production_grid.uperp_faces[k * stride];
    for (int j = 0; j < grid.nu; ++j) {
        grid.dupar[j] = grid.upar_faces[j + 1] - grid.upar_faces[j];
        grid.upar_cells[j] = 0.5 * (grid.upar_faces[j] + grid.upar_faces[j + 1]);
    }
    for (int k = 0; k < grid.nup; ++k) {
        const double lo = grid.uperp_faces[k];
        const double hi = grid.uperp_faces[k + 1];
        grid.duperp[k] = hi - lo;
        grid.uperp_cells[k] = 0.5 * (lo + hi);
        grid.uperp_ring_areas[k] = Const::pi * (hi * hi - lo * lo);
    }
    const double symmetry_scale = 64.0 * std::numeric_limits<double>::epsilon() *
        Param::momentum_umax;
    for (int j = 0; j < grid.nu; ++j)
        grid.geometry_valid = grid.geometry_valid && grid.dupar[j] > 0.0 &&
            grid.upar_faces[j] < grid.upar_faces[j + 1];
    for (int k = 0; k < grid.nup; ++k)
        grid.geometry_valid = grid.geometry_valid && grid.duperp[k] > 0.0 &&
            grid.uperp_faces[k] < grid.uperp_faces[k + 1];
    for (int jf = 0; jf <= grid.nu; ++jf)
        grid.geometry_valid = grid.geometry_valid &&
            std::fabs(grid.upar_faces[jf] + grid.upar_faces[grid.nu - jf]) <=
                symmetry_scale;
    double ring_sum = 0.0;
    for (int k = 0; k < grid.nup; ++k) ring_sum += grid.uperp_ring_areas[k];
    grid.geometry_valid = grid.geometry_valid &&
        grid.uperp_faces.front() == 0.0 &&
        std::fabs(grid.uperp_faces.back() - Param::momentum_umax) <= symmetry_scale &&
        std::fabs(ring_sum - Const::pi * Param::momentum_umax * Param::momentum_umax) /
            (Const::pi * Param::momentum_umax * Param::momentum_umax) <= 1.0e-13;
    for (int j = 0; j < grid.nu; ++j) {
        for (int k = 0; k < grid.nup; ++k) {
            const size_t slot = vk_index(grid, j, k);
            const double gamma = std::sqrt(1.0 + grid.upar_cells[j] * grid.upar_cells[j] +
                grid.uperp_cells[k] * grid.uperp_cells[k]);
            grid.kinetic_energy[slot] = Const::me * Const::c * Const::c * (gamma - 1.0);
            grid.vx[slot] = Const::c * grid.upar_cells[j] / gamma;
        }
    }
    return grid;
}

Fluxes build_production_low_order_fluxes(const GridLevel& grid,
                                         const std::vector<double>& mass,
                                         const std::vector<double>& e_cell,
                                         double charge, double particle_mass)
{
    Fluxes result;
    result.x.assign(static_cast<size_t>(grid.nx) * grid.nu * grid.nup, 0.0);
    result.cu.assign(static_cast<size_t>(grid.nx) * (grid.nu + 1) * grid.nup, 0.0);
    result.u.assign(result.cu.size(), 0.0);

    // Exact production low-order x-face formula:
    // Phi_x = v_x M_donor / dx, with periodic face topology.
    for (int iface = 0; iface < grid.nx; ++iface) {
        const int left = (iface + grid.nx - 1) % grid.nx;
        const int right = iface;
        for (int j = 0; j < grid.nu; ++j) {
            for (int k = 0; k < grid.nup; ++k) {
                const double velocity = grid.vx[vk_index(grid, j, k)];
                const int donor = (velocity >= 0.0) ? left : right;
                result.x[xface_index(grid, iface, j, k)] = velocity *
                    mass[cell_index(grid, donor, j, k)] / grid.dx;
            }
        }
    }

    // Exact production low-order cylindrical u_parallel-face contract:
    // C_u = M_donor / du_donor, Phi_u = a_x C_u.
    for (int ix = 0; ix < grid.nx; ++ix) {
        const double acceleration = charge * e_cell[ix] / (particle_mass * Const::c);
        for (int jf = 1; jf < grid.nu; ++jf) {
            const int donor = (acceleration >= 0.0) ? jf - 1 : jf;
            const double donor_du = grid.dupar[donor];
            for (int k = 0; k < grid.nup; ++k) {
                const size_t face = uface_index(grid, ix, jf, k);
                result.cu[face] = Stage5::donor_cell_coefficient(
                    mass[cell_index(grid, ix, donor, k)], donor_du);
                result.u[face] = acceleration * result.cu[face];
            }
        }
        if (acceleration < 0.0) {
            for (int k = 0; k < grid.nup; ++k) {
                const size_t face = uface_index(grid, ix, 0, k);
                result.cu[face] = Stage5::donor_cell_coefficient(
                    mass[cell_index(grid, ix, 0, k)], grid.dupar[0]);
                result.u[face] = acceleration * result.cu[face];
            }
        } else if (acceleration > 0.0) {
            for (int k = 0; k < grid.nup; ++k) {
                const size_t face = uface_index(grid, ix, grid.nu, k);
                result.cu[face] = Stage5::donor_cell_coefficient(
                    mass[cell_index(grid, ix, grid.nu - 1, k)],
                    grid.dupar[grid.nu - 1]);
                result.u[face] = acceleration * result.cu[face];
            }
        }
    }
    return result;
}

std::vector<double> update_low_order(const GridLevel& grid,
                                     const std::vector<double>& mass,
                                     const Fluxes& fluxes, double h)
{
    std::vector<double> next(mass.size(), 0.0);
    for (int ix = 0; ix < grid.nx; ++ix) {
        const int right = (ix + 1) % grid.nx;
        for (int j = 0; j < grid.nu; ++j) {
            for (int k = 0; k < grid.nup; ++k) {
                next[cell_index(grid, ix, j, k)] = mass[cell_index(grid, ix, j, k)] - h * (
                    fluxes.x[xface_index(grid, right, j, k)] -
                    fluxes.x[xface_index(grid, ix, j, k)] +
                    fluxes.u[uface_index(grid, ix, j + 1, k)] -
                    fluxes.u[uface_index(grid, ix, j, k)]);
            }
        }
    }
    return next;
}

double total_mass(const std::vector<double>& mass)
{
    double result = 0.0;
    for (size_t i = 0; i < mass.size(); ++i) result += mass[i];
    return result;
}

Result run_level(int refinement)
{
    const GridLevel grid = make_level(refinement);
    const double charge = -Const::qe;
    const double particle_mass = Const::me;
    const double epsilon = 1.0e-3;
    // Resolve the reference distribution on all three nested grids.  This is
    // essential for the observed trend to measure operator convergence rather
    // than under-resolved quadrature of the initial condition.
    const double uth = 1.20;
    const double u_drift = 0.35;
    const double e_amplitude = 2.0e3;
    std::vector<double> e_face(grid.nx + 1, 0.0);
    std::vector<double> e_cell(grid.nx, 0.0);
    std::vector<double> mass(static_cast<size_t>(grid.nx) * grid.nu * grid.nup, 0.0);

    for (int iface = 0; iface < grid.nx; ++iface) {
        const double phase = 2.0 * Const::pi * iface / grid.nx;
        e_face[iface] = e_amplitude * (1.0 + 0.2 * std::sin(phase));
    }
    e_face[grid.nx] = e_face[0]; // The final face is only the periodic alias.
    for (int ix = 0; ix < grid.nx; ++ix)
        e_cell[ix] = 0.5 * (e_face[ix] + e_face[ix + 1]); // Production G.

    for (int ix = 0; ix < grid.nx; ++ix) {
        const double phase = 2.0 * Const::pi * (ix + 0.5) / grid.nx;
        const double perturbation = 1.0 + epsilon * std::cos(phase);
        for (int j = 0; j < grid.nu; ++j) {
            for (int k = 0; k < grid.nup; ++k) {
                const double upar = grid.upar_cells[j];
                const double uperp = grid.uperp_cells[k];
                const double f0 = std::exp(-0.5 * ((upar - u_drift) *
                    (upar - u_drift) + uperp * uperp) /
                    (uth * uth));
                mass[cell_index(grid, ix, j, k)] = perturbation * f0 * grid.dx *
                    grid.dupar[j] * grid.uperp_ring_areas[k];
            }
        }
    }

    double max_x_rate = 0.0;
    double max_u_rate = 0.0;
    for (int ix = 0; ix < grid.nx; ++ix) {
        const double acceleration = std::fabs(charge * e_cell[ix] /
                                               (particle_mass * Const::c));
        for (int j = 0; j < grid.nu; ++j) {
            const double u_rate = acceleration / grid.dupar[j];
            max_u_rate = std::max(max_u_rate, u_rate);
            for (int k = 0; k < grid.nup; ++k) {
                max_x_rate = std::max(max_x_rate,
                    std::fabs(grid.vx[vk_index(grid, j, k)]) / grid.dx);
            }
        }
    }
    const double dt = 0.15 / (max_x_rate + max_u_rate);

    // Use one common low-order midpoint predictor.  Both J_N and J_E below
    // are constructed from exactly this state and the same fixed E_c.
    const Fluxes predictor_flux = build_production_low_order_fluxes(
        grid, mass, e_cell, charge, particle_mass);
    const std::vector<double> mass_mid = update_low_order(grid, mass, predictor_flux, 0.5 * dt);
    const Fluxes midpoint_flux = build_production_low_order_fluxes(
        grid, mass_mid, e_cell, charge, particle_mass);

    std::vector<double> jn(grid.nx, 0.0);
    std::vector<double> je(grid.nx, 0.0);
    for (int iface = 0; iface < grid.nx; ++iface) {
        for (int j = 0; j < grid.nu; ++j)
            for (int k = 0; k < grid.nup; ++k)
                jn[iface] += charge * midpoint_flux.x[xface_index(grid, iface, j, k)];
    }
    for (int ix = 0; ix < grid.nx; ++ix) {
        for (int jf = 1; jf < grid.nu; ++jf) {
            for (int k = 0; k < grid.nup; ++k) {
                const double delta_k = grid.kinetic_energy[vk_index(grid, jf, k)] -
                    grid.kinetic_energy[vk_index(grid, jf - 1, k)];
                je[ix] += charge * delta_k * midpoint_flux.cu[
                    uface_index(grid, ix, jf, k)] /
                    (particle_mass * Const::c * grid.dx);
            }
        }
    }

    std::vector<double> gstar_je(grid.nx, 0.0);
    double current_diff_linf = 0.0;
    double current_scale_linf = 0.0;
    double face_work = 0.0;
    double cell_work = 0.0;
    double g_adjoint_face = 0.0;
    double g_adjoint_cell = 0.0;
    for (int iface = 0; iface < grid.nx; ++iface) {
        gstar_je[iface] = 0.5 * (je[(iface + grid.nx - 1) % grid.nx] + je[iface]);
        current_diff_linf = std::max(current_diff_linf,
            std::fabs(jn[iface] - gstar_je[iface]));
        current_scale_linf = std::max(current_scale_linf,
            std::max(std::fabs(jn[iface]), std::fabs(gstar_je[iface])));
        face_work += e_face[iface] * jn[iface] * grid.dx;
        g_adjoint_face += e_face[iface] * gstar_je[iface] * grid.dx;
    }
    for (int ix = 0; ix < grid.nx; ++ix) {
        cell_work += e_cell[ix] * je[ix] * grid.dx;
        g_adjoint_cell += e_cell[ix] * je[ix] * grid.dx;
    }

    Result result;
    result.nx = grid.nx;
    result.nu = grid.nu;
    result.nup = grid.nup;
    result.dt = dt;
    result.x_cfl = dt * max_x_rate;
    result.u_cfl = dt * max_u_rate;
    result.j_pair_relative = current_diff_linf /
        std::max(std::numeric_limits<double>::min(), current_scale_linf);
    result.work_pair_relative = relative_error(face_work - cell_work, face_work, cell_work);
    result.g_adjoint_relative = relative_error(g_adjoint_face - g_adjoint_cell,
                                               g_adjoint_face, g_adjoint_cell);
    result.mass_relative = relative_error(total_mass(mass_mid) - total_mass(mass),
                                          total_mass(mass_mid), total_mass(mass));
    result.min_mass_mid = *std::min_element(mass_mid.begin(), mass_mid.end());
    result.min_dupar = *std::min_element(grid.dupar.begin(), grid.dupar.end());
    result.max_dupar = *std::max_element(grid.dupar.begin(), grid.dupar.end());
    result.min_duperp = *std::min_element(grid.duperp.begin(), grid.duperp.end());
    result.max_duperp = *std::max_element(grid.duperp.begin(), grid.duperp.end());
    result.upar_stretch = grid.upar_stretch;
    result.uperp_stretch = grid.uperp_stretch;
    result.geometry_valid = grid.geometry_valid;
    return result;
}

void print_result(const Result& result)
{
    std::cout << "Nx=" << result.nx << " Nupar=" << result.nu
              << " Nuperp=" << result.nup << " dt=" << result.dt
              << " x_cfl=" << result.x_cfl << " u_cfl=" << result.u_cfl
              << " min_dupar=" << result.min_dupar
              << " max_dupar=" << result.max_dupar
              << " min_duperp=" << result.min_duperp
              << " max_duperp=" << result.max_duperp
              << " upar_stretch=" << result.upar_stretch
              << " uperp_stretch=" << result.uperp_stretch
              << " geometry_valid=" << result.geometry_valid
              << " JN_minus_GstarJE_rel=" << result.j_pair_relative
              << " work_pair_rel=" << result.work_pair_relative
              << " G_Gstar_rel=" << result.g_adjoint_relative
              << " mass_rel=" << result.mass_relative
              << " min_M_mid=" << result.min_mass_mid << "\n";
}

} // namespace

int main()
{
    // Section 7.1 joint low-order audit: Beam, PPM, CTU, FCT and Ampere
    // feedback are intentionally absent.  x is periodic and the background
    // uses the actual cylindrical (u_parallel, u_perp) production contract.
    const Result coarse = run_level(1);
    const Result medium = run_level(2);
    const Result fine = run_level(4);

    std::cout << std::scientific << std::setprecision(16)
              << "fixed_field_low_order_joint_test\n"
              << "beam_enabled=0 ppm_enabled=0 ctu_enabled=0 fct_enabled=0 "
              << "ampere_feedback_enabled=0\n"
              << "state=M0(upar-u_drift,uperp)*(1+epsilon*cos(k*x))\n"
              << "E_cell=G(E_face), JN=charge*sum(Phi_x), "
              << "JE=charge/(m*c*dx)*sum(deltaK*Cu)\n";
    print_result(coarse);
    print_result(medium);
    print_result(fine);

    const double j_order_coarse_to_medium = std::log(
        coarse.j_pair_relative / medium.j_pair_relative) / std::log(2.0);
    const double j_order_medium_to_fine = std::log(
        medium.j_pair_relative / fine.j_pair_relative) / std::log(2.0);
    const double work_order_coarse_to_medium = std::log(
        coarse.work_pair_relative / medium.work_pair_relative) / std::log(2.0);
    const double work_order_medium_to_fine = std::log(
        medium.work_pair_relative / fine.work_pair_relative) / std::log(2.0);
    std::cout << "JN_minus_GstarJE_order_coarse_to_medium="
              << j_order_coarse_to_medium
              << " JN_minus_GstarJE_order_medium_to_fine="
              << j_order_medium_to_fine
              << " work_pair_order_coarse_to_medium="
              << work_order_coarse_to_medium
              << " work_pair_order_medium_to_fine="
              << work_order_medium_to_fine << "\n";

    // The G/G* identity itself is algebraic.  The J/current-work mismatch is
    // an audit of the low-order donor discretization and must decrease under
    // simultaneous x/u_parallel/u_perp refinement with the CFL-scaled dt.
    const bool finite = std::isfinite(coarse.j_pair_relative) &&
        std::isfinite(medium.j_pair_relative) && std::isfinite(fine.j_pair_relative) &&
        std::isfinite(coarse.work_pair_relative) &&
        std::isfinite(medium.work_pair_relative) && std::isfinite(fine.work_pair_relative);
    const bool geometry = coarse.geometry_valid && medium.geometry_valid &&
        fine.geometry_valid && coarse.upar_stretch == medium.upar_stretch &&
        medium.upar_stretch == fine.upar_stretch &&
        coarse.uperp_stretch == medium.uperp_stretch &&
        medium.uperp_stretch == fine.uperp_stretch &&
        medium.min_dupar / coarse.min_dupar >= 0.45 &&
        medium.min_dupar / coarse.min_dupar <= 0.55 &&
        fine.min_dupar / medium.min_dupar >= 0.45 &&
        fine.min_dupar / medium.min_dupar <= 0.55 &&
        medium.min_duperp / coarse.min_duperp >= 0.45 &&
        medium.min_duperp / coarse.min_duperp <= 0.55 &&
        fine.min_duperp / medium.min_duperp >= 0.45 &&
        fine.min_duperp / medium.min_duperp <= 0.55;
    const bool positive = coarse.min_mass_mid >= -1.0e-14 &&
        medium.min_mass_mid >= -1.0e-14 && fine.min_mass_mid >= -1.0e-14;
    const bool conservative = coarse.mass_relative <= 1.0e-12 &&
        medium.mass_relative <= 1.0e-12 && fine.mass_relative <= 1.0e-12;
    const bool adjoint = coarse.g_adjoint_relative <= 1.0e-13 &&
        medium.g_adjoint_relative <= 1.0e-13 && fine.g_adjoint_relative <= 1.0e-13;
    const bool converges = medium.j_pair_relative < coarse.j_pair_relative &&
        fine.j_pair_relative < medium.j_pair_relative &&
        medium.work_pair_relative < coarse.work_pair_relative &&
        fine.work_pair_relative < medium.work_pair_relative;
    const bool first_order_or_roundoff =
        (std::isfinite(j_order_coarse_to_medium) &&
         std::isfinite(j_order_medium_to_fine) &&
         j_order_coarse_to_medium >= 0.7 && j_order_coarse_to_medium <= 1.3 &&
         j_order_medium_to_fine >= 0.7 && j_order_medium_to_fine <= 1.3 &&
         work_order_coarse_to_medium >= 0.7 && work_order_coarse_to_medium <= 1.3 &&
         work_order_medium_to_fine >= 0.7 && work_order_medium_to_fine <= 1.3) ||
        (fine.j_pair_relative <= 4096.0 * std::numeric_limits<double>::epsilon() &&
         fine.work_pair_relative <= 4096.0 * std::numeric_limits<double>::epsilon());
    return (finite && geometry && positive && conservative && adjoint &&
            converges && first_order_or_roundoff) ? 0 : 1;
}
