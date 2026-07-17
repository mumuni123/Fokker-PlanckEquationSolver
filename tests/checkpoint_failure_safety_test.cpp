#include "beam_pic.h"
#include "checkpoint.h"
#include "maxwell.h"
#include "species.h"

#include <cstdio>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <string>

namespace {

bool replace_manifest_line(const std::string& path, const std::string& key,
                           const std::string& replacement)
{
    std::ifstream input(path.c_str());
    std::ostringstream contents;
    contents << input.rdbuf();
    if (!input && !input.eof()) return false;
    std::string text = contents.str();
    const size_t begin = text.find(key);
    if (begin == std::string::npos) return false;
    const size_t end = text.find('\n', begin);
    text.replace(begin, (end == std::string::npos ? text.size() : end) - begin,
                 replacement);
    std::ofstream output(path.c_str(), std::ios::trunc);
    output << text;
    return static_cast<bool>(output);
}

bool write_fresh_checkpoint(const std::string& directory, const Species& bkg,
                            const BeamPIC& beam, const EMFields& fields,
                            const SpatialGrid& grid, int rank, int size)
{
    const CheckpointControlState control = {3, 2.0e-18, 1.0e-18, 0.0, 0, 0.0};
    std::string error;
    return write_checkpoint(directory, control, bkg, beam, fields, grid,
                            rank, size, error);
}

bool expect_collective_read_failure(const std::string& directory,
                                    const SpatialGrid& grid, int rank, int size,
                                    std::string& error)
{
    Species restored;
    restored.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
                  Const::me, Param::dens, Param::temperature_e, false, grid);
    BeamPIC beam;
    beam.init(grid);
    EMFields fields;
    fields.init(grid);
    CheckpointControlState control;
    const bool read_ok = read_checkpoint(directory, control, restored, beam,
                                         fields, grid, rank, size, error);
    const int local_failed_with_reason = (!read_ok && !error.empty()) ? 1 : 0;
    int all_failed_with_reason = 0;
    MPI_Allreduce(&local_failed_with_reason, &all_failed_with_reason, 1,
                  MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return all_failed_with_reason != 0;
}

}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    SpatialGrid grid;
    grid.init(rank, size);
    Species bkg;
    bkg.init("bkg_e", SpeciesType::BACKGROUND_ELECTRON, -Const::qe,
             Const::me, Param::dens, Param::temperature_e, false, grid);
    bkg.initialize_maxwellian();
    BeamPIC beam;
    beam.init(grid);
    EMFields fields;
    fields.init(grid);

    const char* const names[] = {"missing_rank", "bad_mpi_size", "bad_grid",
                                 "truncated_rank", "bad_checksum"};
    bool all_cases_ok = true;
    for (int which = 0; which < 5; ++which) {
        const std::string directory = std::string("checkpoint_failure_safety_tmp/") +
                                      names[which];
        bool prepared = write_fresh_checkpoint(directory, bkg, beam, fields,
                                                grid, rank, size);
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0 && prepared) {
            const std::string manifest = directory + "/manifest.dat";
            const std::string first_rank = directory + "/rank_000000.bin";
            if (which == 0) {
                prepared = std::remove(first_rank.c_str()) == 0;
            } else if (which == 1) {
                prepared = replace_manifest_line(
                    manifest, "mpi_size ", "mpi_size " + std::to_string(size + 1));
            } else if (which == 2) {
                prepared = replace_manifest_line(
                    manifest, "nx ", "nx " + std::to_string(Param::nx + 1));
            } else if (which == 3) {
                std::ofstream truncated(first_rank.c_str(),
                                        std::ios::binary | std::ios::trunc);
                truncated.write("bad", 3);
                prepared = static_cast<bool>(truncated);
            } else {
                prepared = replace_manifest_line(manifest, "rank 0 ",
                                                 "rank 0 1 1");
            }
        }
        int prepared_int = prepared ? 1 : 0;
        MPI_Bcast(&prepared_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        std::string error;
        const bool case_ok = prepared_int != 0 &&
            expect_collective_read_failure(directory, grid, rank, size, error);
        int global_case_ok = 0;
        const int local_case_ok = case_ok ? 1 : 0;
        MPI_Allreduce(&local_case_ok, &global_case_ok, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        if (rank == 0) {
            std::printf("checkpoint_failure_safety case=%s collective_failure=%d "
                        "failure_reason=%s status=%s\n", names[which],
                        global_case_ok, error.c_str(),
                        global_case_ok ? "PASS" : "FAIL");
        }
        all_cases_ok = all_cases_ok && global_case_ok != 0;
    }

    int global_ok = 0;
    const int local_ok = all_cases_ok ? 1 : 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Finalize();
    return global_ok ? 0 : 1;
}
