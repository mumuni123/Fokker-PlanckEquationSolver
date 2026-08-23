// A-FS-R1 manufactured-solution audit.  This test calls the production
// read-only work-identity interface and does not reproduce stable summation.
#include "grid.h"
#include "maxwell.h"
#include "open_electrostatic_solver.h"

#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {
struct Snapshot { std::vector<double> face, ex, phi, rho; };
struct Metrics {
    OpenPoissonWorkIdentity id;
    bool unchanged = false, zero_compatible = false;
    bool rho_resolvable = false, field_resolvable = false;
    double before_integral = 0.0, after_integral = 0.0;
    double integral_tol = 0.0, rho_bound = 0.0, field_bound = 0.0;
    double expected_potential_work = 0.0;
    double expected_field_change = 0.0;
};

bool parse_args(int argc, char** argv, std::string& result)
{
    result = "output/vpfp_poisson_work_identity_unit.result";
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case" && i + 1 < argc) {
            if (std::string(argv[++i]) != "all") return false;
        } else if (arg == "--result" && i + 1 < argc) {
            result = argv[++i];
        } else return false;
    }
    return true;
}

Snapshot snapshot(const EMFields& f)
{
    return {f.Ex_face, f.Ex, f.phi, f.rho};
}

bool same_bits(const Snapshot& a, const Snapshot& b)
{
    if (a.face.size() != b.face.size() || a.ex.size() != b.ex.size() ||
        a.phi.size() != b.phi.size() || a.rho.size() != b.rho.size()) return false;
    return std::memcmp(a.face.data(), b.face.data(), a.face.size()*sizeof(double)) == 0 &&
        std::memcmp(a.ex.data(), b.ex.data(), a.ex.size()*sizeof(double)) == 0 &&
        std::memcmp(a.phi.data(), b.phi.data(), a.phi.size()*sizeof(double)) == 0 &&
        std::memcmp(a.rho.data(), b.rho.data(), a.rho.size()*sizeof(double)) == 0;
}

double half_ulp(double x)
{
    const double up = std::nextafter(x, std::numeric_limits<double>::infinity());
    const double down = std::nextafter(x, -std::numeric_limits<double>::infinity());
    return 0.5 * std::max(up - x, x - down);
}

long double reduce_sum(long double value)
{
    long double result = 0.0L;
    MPI_Allreduce(&value, &result, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return result;
}

double bound(const OpenPoissonWorkIdentity& id, double factor)
{
    return factor * std::numeric_limits<double>::epsilon() * id.scale;
}

bool finite_id(const OpenPoissonWorkIdentity& id)
{
    return id.finite && id.stable_accumulation_used &&
        std::isfinite(id.scale) && id.scale >= 0.0 &&
        std::isfinite(id.field_energy_before) &&
        std::isfinite(id.field_energy_after) &&
        std::isfinite(id.field_energy_change) &&
        std::isfinite(id.electrode_work) &&
        std::isfinite(id.potential_charge_work) && std::isfinite(id.residual);
}

double ratio_to_bound(const OpenPoissonWorkIdentity& id, double factor)
{
    const double b = bound(id, factor);
    if (b > 0.0) return std::fabs(id.residual) / b;
    return std::fabs(id.residual) == 0.0
        ? 0.0 : std::numeric_limits<double>::infinity();
}

double face_x(const SpatialGrid& grid, int local_face)
{
    return grid.x_min + static_cast<double>(grid.ix_start + local_face) * grid.dx;
}

long double face_integral(const EMFields& fields, const SpatialGrid& grid, int rank)
{
    long double sum = 0.0L;
    for (int lf = 0; lf <= grid.nx_local; ++lf) {
        const int gf = grid.ix_start + lf;
        if (rank > 0 && lf == 0) continue; // shared face owned by lower rank
        const long double w = (gf == 0 || gf == grid.nx_global) ? 0.5L : 1.0L;
        sum += w * static_cast<long double>(fields.Ex_face[lf]);
    }
    return reduce_sum(sum * static_cast<long double>(grid.dx));
}

void remove_trapezoid_mean(EMFields& fields, const SpatialGrid& grid, int rank)
{
    long double local_sum = 0.0L;
    for (int lf = 0; lf <= grid.nx_local; ++lf) {
        const int gf = grid.ix_start + lf;
        if (rank > 0 && lf == 0) continue;
        const long double weight = (gf == 0 || gf == grid.nx_global)
            ? 0.5L : 1.0L;
        local_sum += weight * static_cast<long double>(fields.Ex_face[lf]);
    }
    const long double mean = reduce_sum(local_sum) /
        static_cast<long double>(grid.nx_global);
    const double mean_double = static_cast<double>(mean);
    for (double& value : fields.Ex_face) value -= mean_double;
}

void make_state(EMFields& fields, const SpatialGrid& grid,
                OpenElectrostaticSolver& solver, int rank, int size,
                double amplitude, double eta, int sign)
{
    fields.init(grid);
    const double omega = 2.0 * Const::pi / grid.length();
    for (int lf = 0; lf <= grid.nx_local; ++lf) {
        const double x = face_x(grid, lf);
        const double base = amplitude * std::cos(omega * x);
        const double shape = 0.7 * std::cos(omega * x) + 0.3 * std::sin(omega * x);
        fields.Ex_face[lf] = base + static_cast<double>(sign) * eta * amplitude * shape;
    }
    remove_trapezoid_mean(fields, grid, rank);
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        fields.Ex[grid.nghost + ix] = 0.5 * (fields.Ex_face[ix] + fields.Ex_face[ix+1]);
        const long double rho = static_cast<long double>(Const::eps0) *
            (static_cast<long double>(fields.Ex_face[ix+1]) - fields.Ex_face[ix]) /
            static_cast<long double>(grid.dx);
        fields.rho[grid.nghost + ix] = static_cast<double>(rho);
    }
    solver.reconstruct_phi(fields, rank, size);
}

std::vector<double> gauss_delta(const EMFields& before, const EMFields& after,
                                const SpatialGrid& grid)
{
    std::vector<double> delta(static_cast<size_t>(grid.nx_total), 0.0);
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const long double dl = static_cast<long double>(after.Ex_face[ix]) - before.Ex_face[ix];
        const long double dr = static_cast<long double>(after.Ex_face[ix+1]) - before.Ex_face[ix+1];
        delta[grid.nghost + ix] = static_cast<double>(
            static_cast<long double>(Const::eps0) * (dr-dl) /
            static_cast<long double>(grid.dx));
    }
    return delta;
}

Metrics evaluate(const EMFields& before, const EMFields& after,
                 const std::vector<double>& delta, const SpatialGrid& grid,
                 const OpenElectrostaticSolver& solver, int rank, int size,
                 double amplitude)
{
    Metrics m;
    const Snapshot b = snapshot(before), a = snapshot(after);
    m.id = solver.evaluate_work_identity(before, after, delta, rank, size);
    m.unchanged = same_bits(b, snapshot(before)) && same_bits(a, snapshot(after));
    m.before_integral = static_cast<double>(face_integral(before, grid, rank));
    m.after_integral = static_cast<double>(face_integral(after, grid, rank));
    m.integral_tol = 4096.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::fabs(amplitude)) * grid.length();
    m.zero_compatible = std::fabs(m.before_integral) <= m.integral_tol &&
        std::fabs(m.after_integral) <= m.integral_tol;
    const long double dx = static_cast<long double>(grid.dx);
    const long double energy_factor = static_cast<long double>(Const::eps0) * dx / 6.0L;
    long double rho_bound = 0.0L, field_bound = 0.0L;
    long double expected_potential = 0.0L, expected_field = 0.0L;
    for (int ix = 0; ix < grid.nx_local; ++ix) {
        const int l = ix, r = ix + 1;
        const long double delta_left =
            static_cast<long double>(after.Ex_face[l]) - before.Ex_face[l];
        const long double delta_right =
            static_cast<long double>(after.Ex_face[r]) - before.Ex_face[r];
        const long double expected_delta_rho =
            static_cast<long double>(Const::eps0) * (delta_right - delta_left) / dx;
        const double drho = delta[grid.nghost + ix];
        const long double phi_b = static_cast<long double>(before.phi[grid.nghost+ix]) +
            dx * (static_cast<long double>(before.Ex_face[r])-before.Ex_face[l]) / 12.0L;
        const long double phi_a = static_cast<long double>(after.phi[grid.nghost+ix]) +
            dx * (static_cast<long double>(after.Ex_face[r])-after.Ex_face[l]) / 12.0L;
        expected_potential += 0.5L * (phi_b + phi_a) * expected_delta_rho * dx;
        rho_bound += std::fabs(0.5L*(phi_b+phi_a)) * half_ulp(drho) * dx;
        const long double nl = after.Ex_face[l], nr = after.Ex_face[r];
        const long double dl = nl - before.Ex_face[l], dr = nr - before.Ex_face[r];
        const long double ol = before.Ex_face[l], oright = before.Ex_face[r];
        expected_field += energy_factor *
            (dl * (nl + ol) + dl * nr + ol * dr +
             dr * (nr + oright));
        const long double ul = half_ulp(static_cast<double>(dl));
        const long double ur = half_ulp(static_cast<double>(dr));
        field_bound += energy_factor *
            (std::fabs(2.0L*nl+nr) * ul +
             std::fabs(nl+2.0L*nr) * ur +
             ul*ul + ul*ur + ur*ur);
    }
    m.rho_bound = static_cast<double>(reduce_sum(rho_bound));
    m.field_bound = static_cast<double>(reduce_sum(field_bound));
    m.expected_potential_work = static_cast<double>(reduce_sum(expected_potential));
    m.expected_field_change = static_cast<double>(reduce_sum(expected_field));
    m.rho_resolvable = std::fabs(m.expected_potential_work) >= 128.0*m.rho_bound;
    m.field_resolvable = std::fabs(m.expected_field_change) >= 128.0*m.field_bound;
    return m;
}

bool zero_pass(const Metrics& m)
{
    return finite_id(m.id) && m.unchanged && m.zero_compatible &&
        m.rho_resolvable && m.field_resolvable &&
        std::fabs(m.id.residual) <= bound(m.id, 16384.0);
}

bool scalar_pass(const OpenPoissonWorkIdentity& id, double factor)
{
    return finite_id(id) && std::fabs(id.residual) <= bound(id, factor);
}

void write_metrics(std::ostream& out, const std::string& p, const Metrics& m)
{
    out << p << "_finite=" << (finite_id(m.id) ? 1 : 0) << "\n"
        << p << "_stable_accumulation_used=" << (m.id.stable_accumulation_used ? 1 : 0) << "\n"
        << p << "_bits_unchanged=" << (m.unchanged ? 1 : 0) << "\n"
        << p << "_before_face_integral=" << m.before_integral << "\n"
        << p << "_after_face_integral=" << m.after_integral << "\n"
        << p << "_face_integral_tolerance=" << m.integral_tol << "\n"
        << p << "_zero_endpoint_compatibility_pass=" << (m.zero_compatible ? 1 : 0) << "\n"
        << p << "_rho_delta_input_quantization_bound=" << m.rho_bound << "\n"
        << p << "_field_increment_quantization_bound=" << m.field_bound << "\n"
        << p << "_expected_potential_charge_work=" << m.expected_potential_work << "\n"
        << p << "_expected_field_energy_change=" << m.expected_field_change << "\n"
        << p << "_poisson_identity_scale=" << m.id.scale << "\n"
        << p << "_term_abs_sum_energy_before=" << m.id.term_abs_sum_energy_before << "\n"
        << p << "_term_abs_sum_energy_after=" << m.id.term_abs_sum_energy_after << "\n"
        << p << "_term_abs_sum_potential_charge=" << m.id.term_abs_sum_potential_charge << "\n"
        << p << "_rho_delta_resolvable_pass=" << (m.rho_resolvable ? 1 : 0) << "\n"
        << p << "_field_delta_resolvable_pass=" << (m.field_resolvable ? 1 : 0) << "\n"
        << p << "_residual=" << m.id.residual << "\n"
        << p << "_roundoff_bound_8192=" << bound(m.id,8192.0) << "\n"
        << p << "_roundoff_bound_16384=" << bound(m.id,16384.0) << "\n"
        << p << "_ratio_8192=" << ratio_to_bound(m.id,8192.0) << "\n"
        << p << "_ratio_16384=" << ratio_to_bound(m.id,16384.0) << "\n"
        << p << "_pass_8192=" << (scalar_pass(m.id,8192.0) ? 1 : 0) << "\n"
        << p << "_pass_16384=" << (scalar_pass(m.id,16384.0) ? 1 : 0) << "\n"
        << p << "_poisson_identity_roundoff_bound_8192=" << bound(m.id,8192.0) << "\n"
        << p << "_poisson_identity_roundoff_bound_16384=" << bound(m.id,16384.0) << "\n"
        << p << "_poisson_identity_residual_to_bound_ratio_8192=" << ratio_to_bound(m.id,8192.0) << "\n"
        << p << "_poisson_identity_residual_to_bound_ratio_16384=" << ratio_to_bound(m.id,16384.0) << "\n"
        << p << "_poisson_scalar_identity_pass_8192=" << (scalar_pass(m.id,8192.0) ? 1 : 0) << "\n"
        << p << "_poisson_scalar_identity_pass_16384=" << (scalar_pass(m.id,16384.0) ? 1 : 0) << "\n"
        << p << "_poisson_scalar_identity_pass=" << (scalar_pass(m.id,16384.0) ? 1 : 0) << "\n";
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::string result_path;
    const bool parsed = parse_args(argc, argv, result_path);
    SpatialGrid grid;
    grid.init_with_domain(rank, size, 64, 64.0e-9);
    const ElectrostaticBoundary zero_boundary = {
        ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 0.0, 0.0 };
    OpenElectrostaticSolver solver;
    solver.init(grid, zero_boundary);
    const double amplitude = 2.0e-5;

    auto make = [&](EMFields& f, double eta, int sign, OpenElectrostaticSolver& s) {
        make_state(f, grid, s, rank, size, amplitude, eta, sign);
    };
    EMFields normal_before, normal_after;
    make(normal_before,0.0,1,solver); make(normal_after,3.0e-5,1,solver);
    Metrics normal = evaluate(normal_before,normal_after,
        gauss_delta(normal_before,normal_after,grid),grid,solver,rank,size,amplitude);
    EMFields zero_before, zero_after;
    make(zero_before,0.0,1,solver); make(zero_after,0.0,1,solver);
    Metrics zero = evaluate(zero_before,zero_after,
        gauss_delta(zero_before,zero_after,grid),grid,solver,rank,size,amplitude);

    const double etas[5] = {1e-4,3e-5,1e-5,3e-6,1e-6};
    std::vector<Metrics> plus, minus;
    int selected = -1;
    for (int i=0;i<5;++i) {
        EMFields base,p,n;
        make(base,0.0,1,solver); make(p,etas[i],1,solver); make(n,etas[i],-1,solver);
        plus.push_back(evaluate(base,p,gauss_delta(base,p,grid),grid,solver,rank,size,amplitude));
        minus.push_back(evaluate(base,n,gauss_delta(base,n,grid),grid,solver,rank,size,amplitude));
        const Metrics& a=plus.back(); const Metrics& b=minus.back();
        const bool ratio_ok = a.id.field_energy_before >= 1e4*std::fabs(a.id.field_energy_change) &&
            a.id.field_energy_after >= 1e4*std::fabs(a.id.field_energy_change) &&
            b.id.field_energy_before >= 1e4*std::fabs(b.id.field_energy_change) &&
            b.id.field_energy_after >= 1e4*std::fabs(b.id.field_energy_change) &&
            a.expected_field_change > 0.0 && b.expected_field_change < 0.0 &&
            a.expected_potential_work != 0.0 && b.expected_potential_work != 0.0;
        if (selected<0 && ratio_ok && a.rho_resolvable && a.field_resolvable &&
            b.rho_resolvable && b.field_resolvable && a.zero_compatible && b.zero_compatible)
            selected=i;
    }
    const Metrics large_positive = selected>=0 ? plus[selected] : plus.back();
    const Metrics large_negative = selected>=0 ? minus[selected] : minus.back();
    const bool normal_pass=zero_pass(normal);
    const bool zero_pass_case=zero_pass(zero) && zero.id.field_energy_change_direct==0.0 &&
        zero.id.field_energy_change_from_totals==0.0 && zero.id.field_energy_change_reconstruction_error==0.0 &&
        zero.id.potential_charge_work==0.0 && zero.id.residual==0.0;
    const bool positive_pass=selected>=0 && zero_pass(large_positive) && large_positive.id.field_energy_change>0.0;
    const bool negative_pass=selected>=0 && zero_pass(large_negative) && large_negative.id.field_energy_change<0.0;
    const bool state_unchanged_pass = normal.unchanged && zero.unchanged &&
        large_positive.unchanged && large_negative.unchanged;

    const ElectrostaticBoundary endpoint_boundary = {
        ElectrostaticBoundaryType::DIRICHLET_PHI, 0.0, 100.0, -50.0 };
    OpenElectrostaticSolver endpoint_solver;
    endpoint_solver.init(grid,endpoint_boundary);
    EMFields endpoint_before,endpoint_after;
    make(endpoint_before,0.0,1,endpoint_solver); make(endpoint_after,3e-5,1,endpoint_solver);
    const std::vector<double> endpoint_delta=gauss_delta(endpoint_before,endpoint_after,grid);
    const Snapshot eb= snapshot(endpoint_before), ea=snapshot(endpoint_after);
    const OpenPoissonWorkIdentity endpoint_id=endpoint_solver.evaluate_work_identity(
        endpoint_before,endpoint_after,endpoint_delta,rank,size);
    const bool endpoint_bits=same_bits(eb,snapshot(endpoint_before)) && same_bits(ea,snapshot(endpoint_after));
    const double endpoint_ratio=ratio_to_bound(endpoint_id,16384.0);
    const bool endpoint_identity=finite_id(endpoint_id)&&endpoint_bits&&std::isfinite(endpoint_ratio)&&endpoint_ratio<=1.0;
    const bool endpoint_complete=finite_id(endpoint_id)&&endpoint_bits&&std::isfinite(endpoint_ratio);

    int local_gate=(parsed&&normal_pass&&zero_pass_case&&positive_pass&&negative_pass&&
                    state_unchanged_pass&&endpoint_complete)?1:0;
    int global_gate=0;
    MPI_Allreduce(&local_gate,&global_gate,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
    if(rank==0){
        std::ofstream out(result_path.c_str(),std::ios::trunc);
        if(!out) global_gate=0;
        out<<std::setprecision(17)
           <<"stable_accumulation_used="<<(normal.id.stable_accumulation_used?1:0)<<"\n"
           <<"selected_eta_index="<<selected<<"\n"
           <<"selected_eta="<<(selected>=0?etas[selected]:0.0)<<"\n"
           <<"case_normal_pass="<<(normal_pass?1:0)<<"\n"
           <<"case_zero_change_pass="<<(zero_pass_case?1:0)<<"\n"
           <<"case_large_baseline_small_delta_positive_pass="<<(positive_pass?1:0)<<"\n"
           <<"case_large_baseline_small_delta_negative_pass="<<(negative_pass?1:0)<<"\n"
           <<"zero_endpoint_case_normal_pass="<<(normal_pass?1:0)<<"\n"
           <<"zero_endpoint_case_zero_change_pass="<<(zero_pass_case?1:0)<<"\n"
           <<"zero_endpoint_case_large_positive_pass="<<(positive_pass?1:0)<<"\n"
           <<"zero_endpoint_case_large_negative_pass="<<(negative_pass?1:0)<<"\n"
           <<"state_unchanged_pass="<<(state_unchanged_pass?1:0)<<"\n"
           <<"zero_endpoint_production_gate_pass="<<(normal_pass&&zero_pass_case&&positive_pass&&negative_pass&&state_unchanged_pass?1:0)<<"\n"
           <<"case_endpoint_boundary_work_pass="<<(endpoint_complete?1:0)<<"\n"
           <<"poisson_identity_scale="<<large_positive.id.scale<<"\n"
           <<"poisson_identity_roundoff_bound_8192="<<bound(large_positive.id,8192.0)<<"\n"
           <<"poisson_identity_roundoff_bound_16384="<<bound(large_positive.id,16384.0)<<"\n"
           <<"poisson_identity_residual_to_bound_ratio_8192="<<ratio_to_bound(large_positive.id,8192.0)<<"\n"
           <<"poisson_identity_residual_to_bound_ratio_16384="<<ratio_to_bound(large_positive.id,16384.0)<<"\n"
           <<"poisson_scalar_identity_pass_8192="<<(scalar_pass(large_positive.id,8192.0)?1:0)<<"\n"
           <<"poisson_scalar_identity_pass_16384="<<(scalar_pass(large_positive.id,16384.0)?1:0)<<"\n"
           <<"poisson_scalar_identity_pass="<<(scalar_pass(large_positive.id,16384.0)?1:0)<<"\n"
           <<"before_face_integral="<<large_positive.before_integral<<"\n"
           <<"after_face_integral="<<large_positive.after_integral<<"\n"
           <<"before_face_integral_tolerance="<<large_positive.integral_tol<<"\n"
           <<"after_face_integral_tolerance="<<large_positive.integral_tol<<"\n"
           <<"zero_endpoint_compatibility_pass="<<(large_positive.zero_compatible?1:0)<<"\n"
           <<"rho_delta_input_quantization_bound="<<large_positive.rho_bound<<"\n"
           <<"field_increment_quantization_bound="<<large_positive.field_bound<<"\n"
           <<"expected_potential_charge_work="<<large_positive.expected_potential_work<<"\n"
           <<"expected_field_energy_change="<<large_positive.expected_field_change<<"\n"
           <<"rho_delta_resolvable_pass="<<(large_positive.rho_resolvable?1:0)<<"\n"
           <<"field_delta_resolvable_pass="<<(large_positive.field_resolvable?1:0)<<"\n"
           <<"term_abs_sum_energy_before="<<large_positive.id.term_abs_sum_energy_before<<"\n"
           <<"term_abs_sum_energy_after="<<large_positive.id.term_abs_sum_energy_after<<"\n"
           <<"term_abs_sum_potential_charge="<<large_positive.id.term_abs_sum_potential_charge<<"\n"
           <<"endpoint_nonzero_identity_pass="<<(endpoint_identity?1:0)<<"\n"
           <<"endpoint_nonzero_residual="<<endpoint_id.residual<<"\n"
           <<"endpoint_nonzero_electrode_work="<<endpoint_id.electrode_work<<"\n"
           <<"endpoint_nonzero_roundoff_bound="<<bound(endpoint_id,16384.0)<<"\n"
           <<"endpoint_nonzero_residual_to_bound_ratio="<<endpoint_ratio<<"\n"
           <<"endpoint_nonzero_known_limitation="<<(endpoint_identity?0:1)<<"\n"
           <<"endpoint_nonzero_in_production_scope=0\n"
           <<"nonzero_endpoint_diagnostic_complete="<<(endpoint_complete?1:0)<<"\n";
        write_metrics(out,"normal",normal); write_metrics(out,"zero_change",zero);
        write_metrics(out,"large_positive",large_positive); write_metrics(out,"large_negative",large_negative);
        for(int i=0;i<5;++i){
            out<<"eta_candidate_"<<i<<"="<<etas[i]<<"\n"
               <<"eta_candidate_"<<i<<"_positive_field_ratio="<<plus[i].id.field_energy_before/
                 std::max(std::numeric_limits<double>::denorm_min(),std::fabs(plus[i].id.field_energy_change))<<"\n"
               <<"eta_candidate_"<<i<<"_positive_field_after_ratio="<<plus[i].id.field_energy_after/
                 std::max(std::numeric_limits<double>::denorm_min(),std::fabs(plus[i].id.field_energy_change))<<"\n"
               <<"eta_candidate_"<<i<<"_negative_field_ratio="<<minus[i].id.field_energy_before/
                 std::max(std::numeric_limits<double>::denorm_min(),std::fabs(minus[i].id.field_energy_change))<<"\n"
               <<"eta_candidate_"<<i<<"_negative_field_after_ratio="<<minus[i].id.field_energy_after/
                 std::max(std::numeric_limits<double>::denorm_min(),std::fabs(minus[i].id.field_energy_change))<<"\n"
               <<"eta_candidate_"<<i<<"_positive_rho_resolvable="<<(plus[i].rho_resolvable?1:0)<<"\n"
               <<"eta_candidate_"<<i<<"_positive_field_resolvable="<<(plus[i].field_resolvable?1:0)<<"\n"
               <<"eta_candidate_"<<i<<"_negative_rho_resolvable="<<(minus[i].rho_resolvable?1:0)<<"\n"
               <<"eta_candidate_"<<i<<"_negative_field_resolvable="<<(minus[i].field_resolvable?1:0)<<"\n";
        }
        out<<"afs_all_cases_pass="<<global_gate<<"\n"<<"status="<<(global_gate?"PASS":"FAIL")<<"\n";
        std::cout<<"status="<<(global_gate?"PASS":"FAIL")<<"\n";
    }
    MPI_Finalize();
    return global_gate?0:1;
}
