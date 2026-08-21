#include "conservative_ppm_remap.h"
#include "grid.h"
#include "parameters.h"
#include "species.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string result_path;
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--result" && i + 1 < argc)
            result_path = argv[++i];
    SpatialGrid sg;
    sg.init_with_domain(rank, 1, 40, 40.0 * Const::micro);
    CylindricalVelocityGrid cg;
    cg.init(Param::momentum_umax, Param::Nv, Param::Nmu);
    HybridVelocityPartition partition;
    partition.init(cg, 6.0, 1.0, 4, 4);
    Species input, output, output8, output_null, output_no_sink, output_ghost;
    input.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
               Const::me, Param::dens, Param::temperature_e, false, sg);
    output.init("background", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                Const::me, Param::dens, Param::temperature_e, false, sg);
    output8 = output;
    output_null = output;
    output_no_sink = output;
    output_ghost = output;
    std::fill(input.f.begin(), input.f.end(), 0.0);
    if (partition.upar_interface_faces.empty()) {
        if (rank == 0) std::cerr << "no u_parallel interface in test grid\n";
        MPI_Finalize();
        return 1;
    }
    const BulkTailInterfaceFace& interface_face =
        partition.upar_interface_faces.front();
    const int ix = sg.nghost + 10;
    const int j = interface_face.bulk_iv;
    const int k = interface_face.transverse_index;
    // Keep a short non-negative plateau on the bulk side so the production
    // PPM reconstruction has no artificial edge overshoot in this contract
    // test.  The tail side remains empty and is populated only by the swept
    // interface flux.
    const int plateau_lo = interface_face.outward_sign < 0 ? j : j - 4;
    const int plateau_hi = interface_face.outward_sign < 0 ? j + 4 : j;
    for (int jj = std::max(0, plateau_lo);
         jj <= std::min(Param::Nv - 1, plateau_hi); ++jj)
        input.f[idx3(ix, jj, k)] =
            1.0e20 * cg.upar_widths[static_cast<size_t>(jj)];
    input.f_tmp = input.f;
    EMFields fields;
    fields.init(sg);
    const double local_width =
        cg.upar_faces[static_cast<size_t>(interface_face.face_index)] -
        cg.upar_faces[static_cast<size_t>(interface_face.face_index - 1)];
    const double displacement = 0.75 * local_width;
    const double acceleration =
        static_cast<double>(interface_face.outward_sign) * displacement / 1.0e-15;
    const double test_field = acceleration * Const::me * Const::c / (-Const::qe);
    std::fill(fields.Ex.begin(), fields.Ex.end(), test_field);
    std::fill(fields.Ex_face.begin(), fields.Ex_face.end(), fields.Ex[0]);
    ConservativePpmRemap remap;
    remap.init(sg, cg);
    // NULL interface regression: the interface-aware entry point with
    // exporting enabled but sink disabled must be bitwise identical to the
    // original no-interface call.  This compares the actual production
    // remap, not a test-side reconstruction.
    const RemapDiagnostics dnull = remap.advect_u_parallel(
        input, output_null, fields, 1.0e-15, 0.0, NULL, NULL, 4);
    BulkTailFluxBatch no_sink_batch;
    no_sink_batch.apply_interface_sink = false;
    const RemapDiagnostics dno_sink = remap.advect_u_parallel(
        input, output_no_sink, fields, 1.0e-15, 0.0, &partition,
        &no_sink_batch, 4);
    const bool null_interface_bitwise_equal =
        dnull.finite && dno_sink.finite && output_null.f == output_no_sink.f;
    BulkTailFluxBatch batch;
    batch.apply_interface_sink = true;
    const RemapDiagnostics d = remap.advect_u_parallel(
        input, output, fields, 1.0e-15, 0.0, &partition, &batch, 4);
    BulkTailFluxBatch batch8;
    batch8.apply_interface_sink = true;
    const RemapDiagnostics d8 = remap.advect_u_parallel(
        input, output8, fields, 1.0e-15, 0.0, &partition, &batch8, 8);

    // Guard cells are communication copies, not independent physical
    // conversion cells.  Populate one guard slice with the same donor
    // profile and verify that no flux/tail ledger changes.  This catches the
    // production failure where guard tail balances were counted without a
    // matching interior face export.
    Species ghost_input = input;
    const int ghost_ix = 0;
    for (int jj = std::max(0, plateau_lo);
         jj <= std::min(Param::Nv - 1, plateau_hi); ++jj)
        ghost_input.f[idx3(ghost_ix, jj, k)] =
            1.0e20 * cg.upar_widths[static_cast<size_t>(jj)];
    ghost_input.f_tmp = ghost_input.f;
    BulkTailFluxBatch ghost_batch;
    ghost_batch.apply_interface_sink = true;
    const RemapDiagnostics dghost = remap.advect_u_parallel(
        ghost_input, output_ghost, fields, 1.0e-15, 0.0, &partition,
        &ghost_batch, 4);
    const bool guard_ledger_invariant = dghost.finite &&
        dghost.interface_face_export_number == d.interface_face_export_number &&
        dghost.interface_parcel_number == d.interface_parcel_number &&
        dghost.tail_owned_expected_transfer_number ==
            d.tail_owned_expected_transfer_number &&
        dghost.tail_owned_roundoff_discarded_number ==
            d.tail_owned_roundoff_discarded_number &&
        dghost.tail_owned_bulk_residual == d.tail_owned_bulk_residual;

    // Regression for the interface roundoff gate.  A face transfer far below
    // the finite-volume n0*dx resolution must not create a pathological PIC
    // parcel set, and the production sink must not remove it from bulk.
    Species tiny_input = input;
    Species tiny_output = output_no_sink;
    Species tiny_output_sink = output;
    std::fill(tiny_input.f.begin(), tiny_input.f.end(), 0.0);
    const double tiny_amplitude = 1.0e-100;
    for (int jj = std::max(0, plateau_lo);
         jj <= std::min(Param::Nv - 1, plateau_hi); ++jj)
        tiny_input.f[idx3(ix, jj, k)] =
            tiny_amplitude * cg.upar_widths[static_cast<size_t>(jj)];
    tiny_input.f_tmp = tiny_input.f;
    BulkTailFluxBatch tiny_batch;
    tiny_batch.apply_interface_sink = false;
    const RemapDiagnostics dtiny = remap.advect_u_parallel(
        tiny_input, tiny_output, fields, 1.0e-15, 0.0, &partition,
        &tiny_batch, 4);
    BulkTailFluxBatch tiny_sink_batch;
    tiny_sink_batch.apply_interface_sink = true;
    const RemapDiagnostics dtiny_sink = remap.advect_u_parallel(
        tiny_input, tiny_output_sink, fields, 1.0e-15, 0.0, &partition,
        &tiny_sink_batch, 4);
    const double tiny_face = dtiny.interface_face_export_number;
    const double tiny_sink_balance = std::fabs(
        dtiny_sink.number_before - dtiny_sink.number_after -
        dtiny_sink.interface_export_number) /
        std::max(1.0, std::fabs(dtiny_sink.number_before));
    const bool tiny_tail_ok = dtiny.finite && dtiny.audit_valid &&
        dtiny.audit_parcel_failure_reason ==
            static_cast<int>(ParcelNodeFailureReason::None) &&
        tiny_face == 0.0 && dtiny.interface_parcel_number == 0.0 &&
        dtiny.interface_parcel_count == 0 &&
        dtiny.interface_roundoff_discarded_number > 0.0 &&
        dtiny_sink.finite && dtiny_sink.interface_export_number == 0.0 &&
        tiny_sink_batch.parcels.empty() &&
        dtiny_sink.interface_roundoff_discarded_number > 0.0 &&
        tiny_sink_balance <= 1.0e-12;
    const double balance = std::fabs(d.number_before - d.number_after -
                                      d.interface_export_number);
    const double scale = std::max(1.0, std::fabs(d.number_before));
    const double tail_balance_scale = std::max(
        1.0, std::fabs(d.tail_owned_expected_transfer_number));
    const double tail_balance_relative =
        d.tail_owned_bulk_residual / tail_balance_scale;
    double m4[6] = {}, m8[6] = {};
    bulk_tail_batch_moments(batch, m4);
    bulk_tail_batch_moments(batch8, m8);
    double quadrature_error = 0.0;
    for (int q = 0; q < 6; ++q)
        quadrature_error = std::max(quadrature_error,
            std::fabs(m4[q] - m8[q]) / std::max(1.0, std::fabs(m8[q])));
    const double below_relative =
        std::max(d.interface_below_threshold_number,
                 d8.interface_below_threshold_number) /
        std::max(1.0, std::max(std::fabs(d.interface_export_number),
                               std::fabs(d8.interface_export_number)));
    std::size_t negative_node_count = 0;
    const BulkTailFluxBatch* batches[2] = {&batch, &batch8};
    for (int b = 0; b < 2; ++b) {
        for (size_t p = 0; p < batches[b]->parcels.size(); ++p) {
            for (size_t q = 0; q < batches[b]->parcels[p].nodes.size(); ++q) {
                const double mass = batches[b]->parcels[p].nodes[q].mass;
                if (!std::isfinite(mass) || mass < 0.0) ++negative_node_count;
            }
        }
    }
    bool ok = null_interface_bitwise_equal && guard_ledger_invariant && tiny_tail_ok &&
              dno_sink.audit_physical_state_bitwise_equal &&
              d.finite && d8.finite &&
              d.interface_parcel_count > 0 &&
              d.interface_duplicate_count == 0 &&
              d8.interface_duplicate_count == 0 &&
              std::isfinite(d.interface_below_threshold_number) &&
              std::isfinite(d8.interface_below_threshold_number) &&
              balance / scale <= 1.0e-12 && quadrature_error <= 1.0e-11 &&
              std::fabs(d.tail_owned_expected_transfer_number -
                        d.interface_face_export_number) / tail_balance_scale <= 1.0e-12 &&
              tail_balance_relative <= 1.0e-12 &&
              below_relative <= 1.0e-12 && negative_node_count == 0;
    int global_ok = ok ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &global_ok, 1, MPI_INT, MPI_LAND,
                  MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "status=" << (global_ok ? "PASS" : "FAIL")
                  << " balance_relative=" << balance / scale
                  << " parcels=" << d.interface_parcel_count
                  << " before=" << d.number_before
                  << " after=" << d.number_after
                  << " exported=" << d.interface_export_number
                  << " tail_expected="
                  << d.tail_owned_expected_transfer_number
                  << " tail_unexplained=" << d.tail_owned_bulk_residual
                  << " tail_balance_relative=" << tail_balance_relative
                  << " guard_ledger_invariant=" << guard_ledger_invariant
                  << " below_threshold=" << d.interface_below_threshold_number
                  << " below_relative=" << below_relative
                  << " tiny_tail_ok=" << tiny_tail_ok
                  << " tiny_roundoff_discarded="
                  << dtiny.interface_roundoff_discarded_number
                  << " tiny_sink_balance_relative=" << tiny_sink_balance
                  << " negative_nodes=" << negative_node_count
                  << " quadrature_4_vs_8_relative_max=" << quadrature_error
                  << " interface_face=" << interface_face.face_index
                  << " bulk_iv=" << interface_face.bulk_iv
                  << " sign=" << interface_face.outward_sign << "\n";
    }
    if (!result_path.empty() && rank == 0) {
        std::ofstream out(result_path.c_str());
        out << "status=" << (global_ok ? "PASS" : "FAIL") << "\n"
            << "null_interface_bitwise_equal="
            << (null_interface_bitwise_equal ? 1 : 0) << "\n"
            << "audit_inplace_state_bitwise_equal="
            << (dno_sink.audit_physical_state_bitwise_equal ? 1 : 0) << "\n"
            << "guard_ledger_invariant="
            << (guard_ledger_invariant ? 1 : 0) << "\n"
            << "sink_number_relative_error=" << balance / scale << "\n"
            << "tail_owned_expected_transfer_number="
            << d.tail_owned_expected_transfer_number << "\n"
            << "tail_owned_roundoff_discarded_number="
            << d.tail_owned_roundoff_discarded_number << "\n"
            << "tail_owned_unexplained_relative="
            << tail_balance_relative << "\n"
            << "interface_duplicate_count=" << d.interface_duplicate_count << "\n"
            << "quadrature_4_vs_8_relative_max=" << quadrature_error << "\n"
            << "below_threshold_number=" << d.interface_below_threshold_number << "\n"
            << "below_threshold_number_relative=" << below_relative << "\n"
            << "tiny_tail_ok=" << (tiny_tail_ok ? 1 : 0) << "\n"
            << "tiny_tail_face_number=" << tiny_face << "\n"
            << "tiny_tail_parcel_number=" << dtiny.interface_parcel_number << "\n"
            << "tiny_tail_roundoff_discarded_number="
            << dtiny.interface_roundoff_discarded_number << "\n"
            << "tiny_tail_relative_error=0\n"
            << "tiny_sink_balance_relative=" << tiny_sink_balance << "\n"
            << "negative_node_count=" << negative_node_count << "\n";
    }
    MPI_Finalize();
    return global_ok ? 0 : 1;
}
