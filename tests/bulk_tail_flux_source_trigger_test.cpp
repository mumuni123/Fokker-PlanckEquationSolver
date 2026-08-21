// Section 17.10.2.6: controlled, material bulk-to-tail source gate.
//
// This test drives the production u_parallel PPM remap with a constant field
// chosen to move a positive bulk packet through one actual partition face.
// It then loads the emitted face parcel with the production converter.  It is
// deliberately independent of the Beam production case: no production
// threshold is changed and no artificial Beam particle is inserted.

#include "background_tail_pic.h"
#include "bulk_tail_converter.h"
#include "conservative_ppm_remap.h"
#include "grid.h"
#include "parameters.h"
#include "species.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

namespace {

struct Args {
    std::string result;
};

bool parse_args(int argc, char** argv, Args& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--result" && i + 1 < argc) {
            args.result = argv[++i];
        } else if (arg == "--case" && i + 1 < argc) {
            ++i;
        } else {
            return false;
        }
    }
    return true;
}

double physical_mass(const Species& species, const SpatialGrid& grid)
{
    double sum = 0.0;
    for (int il = 0; il < grid.nx_local; ++il) {
        const size_t base = static_cast<size_t>(grid.nghost + il) *
                            static_cast<size_t>(Param::Nvmu);
        for (int j = 0; j < Param::Nv; ++j)
            for (int k = 0; k < Param::Nmu; ++k)
                sum += species.f[base + idx2(j, k)];
    }
    return sum;
}

double tail_weight(const BackgroundTailPIC& tail)
{
    double sum = 0.0;
    for (size_t p = 0; p < tail.particles.size(); ++p)
        sum += tail.particles[p].weight;
    return sum;
}

double relative_error(double a, double b)
{
    return std::fabs(a - b) / std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

struct Metrics {
    bool source_face_found;
    bool remap_finite;
    bool source_flux_seen;
    bool bulk_loss_matches_face;
    bool tail_created;
    bool conversion_complete;
    bool six_moment_closed;
    bool tail_number_matches_created;
    double bulk_loss;
    double face_export;
    double parcel_number;
    double tail_number;
    double bulk_face_relative_error;
    double tail_created_relative_error;
    double six_moment_relative_max;

    Metrics()
        : source_face_found(false), remap_finite(false), source_flux_seen(false),
          bulk_loss_matches_face(false), tail_created(false),
          conversion_complete(false), six_moment_closed(false),
          tail_number_matches_created(false), bulk_loss(0.0), face_export(0.0),
          parcel_number(0.0), tail_number(0.0),
          bulk_face_relative_error(std::numeric_limits<double>::infinity()),
          tail_created_relative_error(std::numeric_limits<double>::infinity()),
          six_moment_relative_max(std::numeric_limits<double>::infinity())
    {}
};

bool write_result(const std::string& path, const Metrics& m, bool pass)
{
    if (path.empty()) return true;
    std::ofstream out(path.c_str(), std::ios::trunc);
    if (!out) return false;
    out << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    out << "source_face_found=" << (m.source_face_found ? 1 : 0) << "\n";
    out << "remap_finite=" << (m.remap_finite ? 1 : 0) << "\n";
    out << "source_flux_seen=" << (m.source_flux_seen ? 1 : 0) << "\n";
    out << "bulk_loss_matches_face=" << (m.bulk_loss_matches_face ? 1 : 0) << "\n";
    out << "tail_created=" << (m.tail_created ? 1 : 0) << "\n";
    out << "conversion_complete=" << (m.conversion_complete ? 1 : 0) << "\n";
    out << "six_moment_closed=" << (m.six_moment_closed ? 1 : 0) << "\n";
    out << "tail_number_matches_created="
        << (m.tail_number_matches_created ? 1 : 0) << "\n";
    out << "bulk_loss=" << m.bulk_loss << "\n";
    out << "face_export=" << m.face_export << "\n";
    out << "parcel_number=" << m.parcel_number << "\n";
    out << "tail_number=" << m.tail_number << "\n";
    out << "bulk_face_relative_error=" << m.bulk_face_relative_error << "\n";
    out << "tail_created_relative_error=" << m.tail_created_relative_error << "\n";
    out << "six_moment_relative_max=" << m.six_moment_relative_max << "\n";
    return static_cast<bool>(out);
}

Metrics run_case()
{
    Metrics m;
    SpatialGrid grid;
    grid.init_with_domain(0, 1, 8, 0.8 * Const::micro);
    Species input;
    input.init("source_bulk", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
               Const::me, Param::dens, Param::temperature_e, false, grid);
    std::fill(input.f.begin(), input.f.end(), 0.0);

    HybridVelocityPartition partition;
    partition.init(input.cgrid, 6.0, 1.0, 4, 4);
    const BulkTailInterfaceFace* source_face = NULL;
    for (size_t n = 0; n < partition.upar_interface_faces.size(); ++n) {
        const BulkTailInterfaceFace& face = partition.upar_interface_faces[n];
        // A positive outward face lets a negative Ex accelerate electrons in
        // the selected direction.  Three upstream bulk cells give PPM a
        // smooth, positive donor plateau rather than a one-cell spike.
        if (face.outward_sign > 0 && face.bulk_iv >= 3) {
            source_face = &face;
            break;
        }
    }
    if (source_face == NULL) return m;
    m.source_face_found = true;

    const int k = source_face->bulk_imu;
    const int bulk_j = source_face->bulk_iv;
    const double donor_mass = 1.0e18;
    for (int il = 0; il < grid.nx_local; ++il) {
        const size_t base = static_cast<size_t>(grid.nghost + il) *
                            static_cast<size_t>(Param::Nvmu);
        for (int j = bulk_j - 3; j <= bulk_j; ++j)
            input.f[base + idx2(j, k)] = donor_mass;
    }

    Species output = input;
    const double dt = Param::dt_multiplier / Param::omega_pe;
    const double donor_du = input.cgrid.upar_widths[static_cast<size_t>(bulk_j)];
    const double acceleration = 0.35 * donor_du / dt;
    const double electric_field = acceleration * Const::me * Const::c / input.charge;
    EMFields fields;
    fields.init(grid);
    std::fill(fields.Ex.begin(), fields.Ex.end(), electric_field);
    std::fill(fields.Ex_face.begin(), fields.Ex_face.end(), electric_field);

    ConservativePpmRemap remap;
    remap.init(grid, input.cgrid);
    BulkTailFluxBatch batch;
    batch.apply_interface_sink = true;
    const double mass_before = physical_mass(input, grid);
    const RemapDiagnostics remap_diag = remap.advect_u_parallel(
        input, output, fields, dt, 0.0, &partition, &batch, 4);
    const double mass_after = physical_mass(output, grid);
    m.remap_finite = remap_diag.finite && batch.finite && batch.nonnegative;
    m.bulk_loss = mass_before - mass_after;
    m.face_export = remap_diag.interface_face_export_number;
    m.parcel_number = remap_diag.interface_parcel_number;
    m.source_flux_seen = m.face_export > 0.0 &&
                         remap_diag.interface_parcel_count > 0 &&
                         !batch.parcels.empty();
    m.bulk_face_relative_error = relative_error(m.bulk_loss, m.face_export);
    m.bulk_loss_matches_face = m.bulk_loss > 0.0 &&
                               m.bulk_face_relative_error <= 1.0e-10;

    BackgroundTailPIC tail;
    tail.init(grid);
    BulkTailConverter converter;
    const BulkTailConversionDiagnostics conversion = converter.convert_flux_batch(
        batch, tail, grid, partition, 1, ConversionLocation::AFTER_U_SUBSTEP,
        0, 7);
    m.tail_number = tail_weight(tail);
    m.tail_created = !tail.particles.empty() && m.tail_number > 0.0;
    m.conversion_complete = conversion.complete && conversion.conservative &&
                            conversion.fidelity_ok && conversion.finite;
    m.tail_created_relative_error = relative_error(
        m.tail_number, conversion.number_created);
    m.tail_number_matches_created =
        m.tail_created_relative_error <= 1.0e-10;
    m.six_moment_relative_max = std::max(
        std::max(conversion.number_residual_rel, conversion.px_residual_rel),
        std::max(std::max(conversion.jx_residual_rel,
                          conversion.energy_residual_rel),
                 std::max(conversion.pixx_residual_rel,
                          conversion.piperp_residual_rel)));
    m.six_moment_closed = m.six_moment_relative_max <= 1.0e-9;
    return m;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Args args;
    const bool parsed = parse_args(argc, argv, args);
    Metrics m;
    if (parsed && size == 1) m = run_case();
    const bool pass = parsed && size == 1 && m.source_face_found &&
        m.remap_finite && m.source_flux_seen && m.bulk_loss_matches_face &&
        m.tail_created && m.conversion_complete && m.six_moment_closed &&
        m.tail_number_matches_created;
    const bool write_ok = write_result(args.result, m, pass);
    std::cout << "status=" << ((pass && write_ok) ? "PASS" : "FAIL")
              << " face_export=" << m.face_export
              << " tail_number=" << m.tail_number << "\n";
    MPI_Finalize();
    return pass && write_ok ? 0 : 1;
}
