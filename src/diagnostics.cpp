#include "diagnostics.h"
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <mpi.h>
#include <sys/stat.h>

Diagnostics::Diagnostics()
    : snapshot_count(0)
{}

void Diagnostics::init(const std::string& dir, int mpi_rank)
{
    output_dir = dir;
    snapshot_count = 0;

    if (mpi_rank == 0) {
        mkdir(output_dir.c_str(), 0755);

        std::string fname = output_dir + "/scalars.dat";
        scalar_file.open(fname.c_str());
        scalar_file << "# step  time[fs]  N_beam  N_bkg_e  N_bkg_i  "
                    << "KE_beam  KE_bkg_e  KE_bkg_i  EM_energy  Total_energy\n";
        scalar_file << std::scientific << std::setprecision(8);
    }
}

void Diagnostics::write_scalars(double time, int step,
                                const std::vector<Species*>& species,
                                const EMFields& fields,
                                int mpi_rank, int mpi_size)
{
    // 汇总本地量
    std::vector<double> local_N(species.size());
    std::vector<double> local_KE(species.size());
    for (size_t s = 0; s < species.size(); ++s) {
        local_N[s]  = species[s]->total_particle_number();
        local_KE[s] = species[s]->total_kinetic_energy();
    }
    double local_EM = fields.total_energy();

    // MPI 归约
    std::vector<double> global_N(species.size());
    std::vector<double> global_KE(species.size());
    double global_EM = 0.0;

    MPI_Reduce(local_N.data(),  global_N.data(),  (int)species.size(),
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_KE.data(), global_KE.data(), (int)species.size(),
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_EM, &global_EM, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        double total_energy = global_EM;
        for (size_t s = 0; s < species.size(); ++s) {
            total_energy += global_KE[s];
        }

        scalar_file << step << "  "
                    << time / Const::femto << "  ";
        for (size_t s = 0; s < species.size(); ++s) {
            scalar_file << global_N[s] << "  ";
        }
        for (size_t s = 0; s < species.size(); ++s) {
            scalar_file << global_KE[s] << "  ";
        }
        scalar_file << global_EM << "  " << total_energy << "\n";
        scalar_file.flush();
    }
}

void Diagnostics::write_fields(double time,
                                const EMFields& fields,
                                const SpatialGrid& sg,
                                int mpi_rank, int mpi_size)
{
    int ng  = sg.nghost;
    int nxl = sg.nx_local;

    // 将 Ex 汇总到 rank 0
    std::vector<double> local_Ex(nxl);
    for (int i = 0; i < nxl; ++i) {
        local_Ex[i] = fields.Ex[i + ng];
    }

    // 汇总各 rank 大小
    std::vector<int> counts(mpi_size);
    std::vector<int> displs(mpi_size);
    MPI_Gather(&nxl, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        displs[0] = 0;
        for (int r = 1; r < mpi_size; ++r) {
            displs[r] = displs[r-1] + counts[r-1];
        }
    }

    int nx_global = sg.nx_global;
    std::vector<double> global_Ex(nx_global);

    MPI_Gatherv(local_Ex.data(), nxl, MPI_DOUBLE,
                global_Ex.data(), counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        std::ostringstream fname;
        fname << output_dir << "/fields_" << std::setw(5) << std::setfill('0')
              << snapshot_count << ".dat";
        std::ofstream out(fname.str().c_str());
        out << "# x[um]  Ex[V/m]\n";
        out << std::scientific << std::setprecision(8);
        for (int i = 0; i < nx_global; ++i) {
            double x = (i + 0.5) * sg.dx;
            out << x / Const::micro << "  " << global_Ex[i] << "\n";
        }
    }
}

void Diagnostics::write_px_distribution(double time,
                                        const Species& sp,
                                        int mpi_rank, int mpi_size)
{
    int ng  = sp.sgrid->nghost;
    int nxl = sp.sgrid->nx_local;
    double dp3 = sp.pgrid.dpx * sp.pgrid.dpy * sp.pgrid.dpz;

    // 对 f 在 x、py、pz 上积分得到 F(px)
    std::vector<double> local_Fpx(Param::Npx, 0.0);

    for (int ix = 0; ix < nxl; ++ix) {
        int ix_g = ix + ng;
        for (int ipx = 0; ipx < Param::Npx; ++ipx) {
            double sum = 0.0;
            for (int ipy = 0; ipy < Param::Npy; ++ipy) {
                for (int ipz = 0; ipz < Param::Npz; ++ipz) {
                    sum += sp.f[idx4(ix_g, ipx, ipy, ipz)];
                }
            }
            local_Fpx[ipx] += sum * sp.pgrid.dpy * sp.pgrid.dpz * sp.sgrid->dx;
        }
    }

    std::vector<double> global_Fpx(Param::Npx, 0.0);
    MPI_Reduce(local_Fpx.data(), global_Fpx.data(), Param::Npx,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        std::ostringstream fname;
        fname << output_dir << "/fpx_" << sp.name << "_"
              << std::setw(5) << std::setfill('0') << snapshot_count << ".dat";
        std::ofstream out(fname.str().c_str());
        out << "# px[kg*m/s]  F(px)\n";
        out << std::scientific << std::setprecision(8);
        for (int ipx = 0; ipx < Param::Npx; ++ipx) {
            out << sp.pgrid.px(ipx) << "  " << global_Fpx[ipx] << "\n";
        }
    }
}

void Diagnostics::write_density_profile(double time,
                                        const std::vector<Species*>& species,
                                        const SpatialGrid& sg,
                                        int mpi_rank, int mpi_size)
{
    int nxl = sg.nx_local;
    int nx_global = sg.nx_global;

    std::vector<int> counts(mpi_size);
    std::vector<int> displs(mpi_size);
    MPI_Gather(&nxl, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        displs[0] = 0;
        for (int r = 1; r < mpi_size; ++r) {
            displs[r] = displs[r-1] + counts[r-1];
        }
    }

    // 汇总各物种密度
    std::vector<std::vector<double> > global_n(species.size());
    for (size_t s = 0; s < species.size(); ++s) {
        global_n[s].resize(nx_global);
        MPI_Gatherv(species[s]->number_density.data(), nxl, MPI_DOUBLE,
                    global_n[s].data(), counts.data(), displs.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
    }

    if (mpi_rank == 0) {
        std::ostringstream fname;
        fname << output_dir << "/density_" << std::setw(5) << std::setfill('0')
              << snapshot_count << ".dat";
        std::ofstream out(fname.str().c_str());
        out << "# x[um]";
        for (size_t s = 0; s < species.size(); ++s) {
            out << "  n_" << species[s]->name;
        }
        out << "\n";
        out << std::scientific << std::setprecision(8);
        for (int i = 0; i < nx_global; ++i) {
            out << (i + 0.5) * sg.dx / Const::micro;
            for (size_t s = 0; s < species.size(); ++s) {
                out << "  " << global_n[s][i];
            }
            out << "\n";
        }
    }

    ++snapshot_count;
}
