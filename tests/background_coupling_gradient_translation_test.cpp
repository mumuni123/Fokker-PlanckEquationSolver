#include "background_coupling_test_support.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

namespace {

struct TranslationRecord {
    BackgroundCouplingTest::Norms norms;
    std::vector<double> residual_global;
    int max_face;
    int state_advanced;
    int finite;
};

TranslationRecord run_case(const SpatialGrid& sg, int rank, int size,
                           double phase)
{
    Species background;
    EMFields fields;
    const int mode = 20;
    BackgroundCouplingTest::initialize_periodic_state(
        background, fields, sg, rank, size, 0.22, 0.25, 2.0e3, mode, phase);
    const VlasovAmpereMidpointSolver::BackgroundCouplingFluxBundle bundle =
        BackgroundCouplingTest::evaluate_bundle(
            background, fields, sg, rank, size,
            BackgroundCouplingTest::stable_dt(sg));
    TranslationRecord record = {};
    record.norms = BackgroundCouplingTest::face_difference_norms(
        bundle.jn_high, bundle.gstar_je_center, sg);
    std::vector<double> local_residual(static_cast<size_t>(sg.nx_local), 0.0);
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        local_residual[static_cast<size_t>(iface)] =
            bundle.jn_high[static_cast<size_t>(iface)] -
            bundle.gstar_je_center[static_cast<size_t>(iface)];
    }
    record.residual_global = BackgroundCouplingTest::gather_owned_faces(
        local_residual, sg, rank, size);
    double local_max = -1.0;
    int local_face = -1;
    for (int iface = 0; iface < sg.nx_local; ++iface) {
        const double d = std::fabs(bundle.jn_high[static_cast<size_t>(iface)] -
                                   bundle.gstar_je_center[static_cast<size_t>(iface)]);
        if (d > local_max) {
            local_max = d;
            local_face = sg.ix_start + iface;
        }
    }
    struct { double value; int index; } local = {local_max, local_face}, global = {};
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
                  MPI_COMM_WORLD);
    record.max_face = global.index;
    record.state_advanced = bundle.state_advanced;
    record.finite = bundle.finite;
    return record;
}

double relative_difference(double a, double b)
{
    return std::fabs(a - b) / std::max(1.0e-300, std::max(std::fabs(a), std::fabs(b)));
}

double shifted_profile_relative_l2(const std::vector<double>& reference,
                                   const std::vector<double>& translated,
                                   int shift)
{
    if (reference.empty() || reference.size() != translated.size()) return 1.0;
    double mismatch = 0.0;
    double scale = 0.0;
    const int n = static_cast<int>(reference.size());
    for (int i = 0; i < n; ++i) {
        const int rotated = (i + shift % n + n) % n;
        const double d = reference[static_cast<size_t>(i)] -
            translated[static_cast<size_t>(rotated)];
        mismatch += d * d;
        scale += reference[static_cast<size_t>(i)] *
            reference[static_cast<size_t>(i)];
    }
    return std::sqrt(mismatch / std::max(1.0e-300, scale));
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    SpatialGrid sg;
    sg.init(rank, size);
    const int mode = 20;
    const double kx = BackgroundCouplingTest::wave_number(sg, mode);
    const TranslationRecord seam = run_case(sg, rank, size, 0.5 * Const::pi);
    // Move the steep feature by one quarter / one half wavelength, not by
    // L/4 or L/2 (which would be an integer number of mode-20 periods).
    const TranslationRecord quarter = run_case(
        sg, rank, size, 0.5 * Const::pi - kx * Param::Lx /
        (4.0 * static_cast<double>(mode)));
    const TranslationRecord half = run_case(
        sg, rank, size, 0.5 * Const::pi - kx * Param::Lx /
        (2.0 * static_cast<double>(mode)));
    const int quarter_shift = sg.nx_global / (4 * mode);
    const int half_shift = sg.nx_global / (2 * mode);
    double quarter_profile_relative = 0.0;
    double half_profile_relative = 0.0;
    if (rank == 0) {
        // The sign follows the chosen phase convention; report the better
        // cyclic orientation so the test checks translation, not convention.
        quarter_profile_relative = std::min(
            shifted_profile_relative_l2(seam.residual_global,
                                        quarter.residual_global, quarter_shift),
            shifted_profile_relative_l2(seam.residual_global,
                                        quarter.residual_global, -quarter_shift));
        half_profile_relative = std::min(
            shifted_profile_relative_l2(seam.residual_global,
                                        half.residual_global, half_shift),
            shifted_profile_relative_l2(seam.residual_global,
                                        half.residual_global, -half_shift));
    }
    double translation_checks[2] = {quarter_profile_relative,
                                    half_profile_relative};
    MPI_Bcast(translation_checks, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    const double translation_profile_tolerance = 1.0e-8;
    const int translation_profiles_match =
        translation_checks[0] <= translation_profile_tolerance &&
        translation_checks[1] <= translation_profile_tolerance;
    const int translation_norms_match =
        relative_difference(quarter.norms.l2, seam.norms.l2) <=
            translation_profile_tolerance &&
        relative_difference(half.norms.l2, seam.norms.l2) <=
            translation_profile_tolerance;
    int local_ok = seam.state_advanced && quarter.state_advanced &&
        half.state_advanced && seam.finite && quarter.finite && half.finite &&
        translation_profiles_match && translation_norms_match;
    MPI_Allreduce(MPI_IN_PLACE, &local_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (rank == 0) {
        std::ofstream out("output/background_coupling_gradient_translation.result");
        std::ostream& log = out ? out : std::cout;
        log << "test=background_coupling_gradient_translation\n"
            << "production_kernel=1\nbeam_enabled=0\n"
            << "mode=" << mode << "\n"
            << "seam_L1=" << seam.norms.l1 << "\n"
            << "quarter_L1=" << quarter.norms.l1 << "\n"
            << "half_L1=" << half.norms.l1 << "\n"
            << "seam_L2=" << seam.norms.l2 << "\n"
            << "quarter_L2=" << quarter.norms.l2 << "\n"
            << "half_L2=" << half.norms.l2 << "\n"
            << "seam_Linf=" << seam.norms.linf << "\n"
            << "quarter_Linf=" << quarter.norms.linf << "\n"
            << "half_Linf=" << half.norms.linf << "\n"
            << "seam_max_face=" << seam.max_face << "\n"
            << "quarter_max_face=" << quarter.max_face << "\n"
            << "half_max_face=" << half.max_face << "\n"
            << "quarter_wavelength_shift_cells=" << quarter_shift << "\n"
            << "half_wavelength_shift_cells=" << half_shift << "\n"
            << "quarter_profile_shift_relative_L2="
            << quarter_profile_relative << "\n"
            << "half_profile_shift_relative_L2="
            << half_profile_relative << "\n"
            << "translation_profile_tolerance="
            << translation_profile_tolerance << "\n"
            << "translation_profiles_match=" << translation_profiles_match << "\n"
            << "translation_norms_match=" << translation_norms_match << "\n"
            << "quarter_vs_seam_L2_relative="
            << relative_difference(quarter.norms.l2, seam.norms.l2) << "\n"
            << "half_vs_seam_L2_relative="
            << relative_difference(half.norms.l2, seam.norms.l2) << "\n"
            << "passes=" << local_ok << "\n";
        std::cout << "background_coupling_gradient_translation_test passes="
                  << local_ok << "\n";
    }
    MPI_Finalize();
    return local_ok ? 0 : 1;
}
