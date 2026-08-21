#include "beam_pic.h"
#include "tail_bulk_return.h"
#include "vpfp_checkpoint.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct Options {
    std::string checkpoint, result;
    std::vector<double> thresholds;
    std::vector<int> radii;
    double tolerance;
    double minimum_feasible_fraction;
    Options() : tolerance(1.0e-12), minimum_feasible_fraction(0.90) {
        thresholds = {5.25, 5.5, 5.75};
        radii = {1, 2, 3};
    }
};
template <typename T> bool parse_list(const std::string& text,
                                      std::vector<T>& values) {
    std::stringstream in(text); std::string token; values.clear();
    while (std::getline(in, token, ',')) {
        std::stringstream item(token); T value; char extra = 0;
        if (!(item >> value) || (item >> extra)) return false;
        values.push_back(value);
    }
    return !values.empty();
}
bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a(argv[i]);
        if (i + 1 >= argc) return false;
        if (a == "--checkpoint") o.checkpoint = argv[++i];
        else if (a == "--result") o.result = argv[++i];
        else if (a == "--thresholds-mev" || a == "--return-thresholds-mev") {
            if (!parse_list(argv[++i], o.thresholds)) return false;
        } else if (a == "--stencil-radii") {
            if (!parse_list(argv[++i], o.radii)) return false;
        } else if (a == "--moment-tolerance") {
            char* end=NULL; o.tolerance=std::strtod(argv[++i],&end);
            if(end==argv[i]||*end!='\0'||!(o.tolerance>0.0))return false;
        } else if (a == "--minimum-feasible-fraction") {
            char* end=NULL;
            o.minimum_feasible_fraction=std::strtod(argv[++i],&end);
            if(end==argv[i]||*end!='\0'||
               !(o.minimum_feasible_fraction>=0.0) ||
               !(o.minimum_feasible_fraction<=1.0)) return false;
        } else return false;
    }
    return !o.checkpoint.empty() && !o.result.empty();
}
std::uint64_t hash_bytes(const void* data, size_t bytes,
                         std::uint64_t seed = UINT64_C(1469598103934665603)) {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (size_t i = 0; i < bytes; ++i) { seed ^= p[i]; seed *= UINT64_C(1099511628211); }
    return seed;
}
std::uint64_t state_hash(const Species& bulk, const BackgroundTailPIC& tail) {
    std::uint64_t h = hash_bytes(bulk.f.data(), bulk.f.size() * sizeof(double));
    for (size_t i = 0; i < tail.particles.size(); ++i) {
        const BackgroundTailParticle& p = tail.particles[i];
        h = hash_bytes(&p.x, sizeof(p.x), h); h = hash_bytes(&p.ux, sizeof(p.ux), h);
        h = hash_bytes(&p.uy, sizeof(p.uy), h); h = hash_bytes(&p.uz, sizeof(p.uz), h);
        h = hash_bytes(&p.weight, sizeof(p.weight), h); h = hash_bytes(&p.id, sizeof(p.id), h);
        h = hash_bytes(&p.return_residence_steps, sizeof(p.return_residence_steps), h);
    }
    return h;
}
double particle_ke(const BackgroundTailParticle& p) {
    return Const::me * Const::c * Const::c *
        (std::sqrt(1.0 + p.ux*p.ux + p.uy*p.uy + p.uz*p.uz) - 1.0);
}
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank=0,size=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
    Options o; int valid=parse(argc,argv,o)?1:0;
    MPI_Allreduce(MPI_IN_PLACE,&valid,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
    if(!valid){if(rank==0)std::cerr<<"invalid tail return audit arguments\n";MPI_Finalize();return 2;}
    SpatialGrid grid; grid.init_with_domain(rank,size,Param::nx,Param::Lx);
    Species bulk; bulk.init("background_electrons",SpeciesType::BACKGROUND_ELECTRON,
        -Const::qe,Const::me,Param::dens,Param::temperature_e,false,grid);
    BeamPIC beam; beam.init(grid); EMFields fields; fields.init(grid);
    VpfpCheckpointControl control={}; VpfpCheckpointTailState state; std::string error;
    int read_ok=read_vpfp_checkpoint(o.checkpoint,control,bulk,beam,fields,grid,
        &state,rank,size,error)&&state.present;
    MPI_Allreduce(MPI_IN_PLACE,&read_ok,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
    if(!read_ok){if(rank==0)std::cerr<<"checkpoint read failed: "<<error<<"\n";MPI_Finalize();return 3;}
    HybridVelocityPartition partition;
    partition.init(bulk.cgrid,state.config.convert_energy_mev,
        state.config.buffer_width_mev,state.config.upar_bins,state.config.energy_bins);
    const std::uint64_t before=state_hash(bulk,state.tail);
    std::ostringstream rows; rows<<std::setprecision(17);
    bool pass=true; bool finite_all=true;
    double central_fraction=0.0; bool central_seen=false;
    for(size_t ti=0;ti<o.thresholds.size();++ti){for(size_t ri=0;ri<o.radii.size();++ri){
        const double threshold=o.thresholds[ti]; const int radius=o.radii[ri];
        Species trial_bulk=bulk; BackgroundTailPIC trial_tail=state.tail;
        double local_candidate_n=0.0,local_candidate_e=0.0;
        unsigned long long local_candidate_macro=0;
        for(size_t p=0;p<trial_tail.particles.size();++p){
            const double ke=particle_ke(trial_tail.particles[p]);
            if(ke<threshold*1.0e6*Const::eV){local_candidate_n+=trial_tail.particles[p].weight;
                local_candidate_e+=trial_tail.particles[p].weight*ke;++local_candidate_macro;}
        }
        double candidate[2]={local_candidate_n,local_candidate_e};
        MPI_Allreduce(MPI_IN_PLACE,candidate,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        unsigned long long candidate_macro=local_candidate_macro;
        MPI_Allreduce(MPI_IN_PLACE,&candidate_macro,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
        TailBulkReturnConfig cfg; cfg.enabled=true; cfg.return_energy_mev=threshold;
        cfg.residence_steps=1; cfg.max_stencil_radius=radius; cfg.moment_tolerance=o.tolerance;
        TailBulkReturnDiagnostics d;
        const bool ok=TailBulkReturn(cfg).apply(trial_bulk,trial_tail,grid,partition,
            control.step+1,rank,size,d);
        unsigned long long group_counts[8] = {
            static_cast<unsigned long long>(d.attempted_groups),
            static_cast<unsigned long long>(d.committed_groups),
            static_cast<unsigned long long>(d.deferred_infeasible_groups),
            static_cast<unsigned long long>(d.particles_removed),
            static_cast<unsigned long long>(d.projection_invalid_input_cells),
            static_cast<unsigned long long>(d.projection_insufficient_support_cells),
            static_cast<unsigned long long>(d.projection_infeasible_invariant_cells),
            static_cast<unsigned long long>(
                d.projection_representation_incompatible_cells) };
        MPI_Allreduce(MPI_IN_PLACE,group_counts,8,MPI_UNSIGNED_LONG_LONG,
                      MPI_SUM,MPI_COMM_WORLD);
        double local_core_deferred=0.0,local_boundary_deferred=0.0;
        const int edge=grid.nx_global/10;
        for(size_t p=0;p<trial_tail.particles.size();++p){
            if(particle_ke(trial_tail.particles[p])>=threshold*1.0e6*Const::eV)continue;
            const int ix=static_cast<int>(std::floor(trial_tail.particles[p].x/grid.dx));
            if(ix>=edge&&ix<grid.nx_global-edge)local_core_deferred+=trial_tail.particles[p].weight;
            else local_boundary_deferred+=trial_tail.particles[p].weight;
        }
        double deferred[2]={local_core_deferred,local_boundary_deferred};
        MPI_Allreduce(MPI_IN_PLACE,deferred,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        const double feasible_fraction=d.number/std::max(candidate[0],1.0e-300);
        if(std::fabs(threshold-5.5)<=1.0e-12 && radius==3){
            central_fraction=feasible_fraction; central_seen=true;
        }
        rows<<"threshold_"<<ti<<"_radius_"<<radius<<"_candidate_number="<<candidate[0]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_candidate_energy="<<candidate[1]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_candidate_macro_particles="<<candidate_macro<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_returned_number="<<d.number<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_feasible_number_fraction="<<feasible_fraction<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_attempted_groups="<<group_counts[0]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_committed_groups="<<group_counts[1]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_deferred_groups="<<group_counts[2]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_particles_removed="<<group_counts[3]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_projection_invalid_input_cells="<<group_counts[4]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_projection_insufficient_support_cells="<<group_counts[5]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_projection_infeasible_invariant_cells="<<group_counts[6]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_projection_representation_incompatible_cells="<<group_counts[7]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_core_deferred_number="<<deferred[0]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_boundary_deferred_number="<<deferred[1]<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_invariant_residual_max="
            <<std::max(d.number_residual,
                       std::max(d.px_residual,d.energy_residual))<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_representation_residual_max="
            <<std::max(d.jx_residual,
                       std::max(d.pixx_residual,d.piperp_residual))<<'\n'
            <<"threshold_"<<ti<<"_radius_"<<radius<<"_moment_residual_max="
            <<std::max(d.number_residual,std::max(d.px_residual,std::max(d.jx_residual,
              std::max(d.energy_residual,std::max(d.pixx_residual,d.piperp_residual)))))<<'\n';
        finite_all=finite_all&&ok&&d.finite;
        pass=pass&&ok&&d.finite;
    }}
    const std::uint64_t after=state_hash(bulk,state.tail);
    int unchanged=before==after?1:0; MPI_Allreduce(MPI_IN_PLACE,&unchanged,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
    pass=pass&&unchanged;
    const bool feasibility_gate_required =
        control.time >= 99.0 * Const::femto;
    const bool feasibility_gate_pass = central_seen &&
        central_fraction >= o.minimum_feasible_fraction;
    if (feasibility_gate_required) pass=pass&&feasibility_gate_pass;
    if(rank==0){std::ofstream out(o.result.c_str());out<<std::setprecision(17)
        <<"schema=tail_bulk_return_checkpoint_audit_v2\ncheckpoint_step="<<control.step
        <<"\ncheckpoint_time_s="<<control.time<<"\nmpi_size="<<size
        <<"\nprojection_schema="<<tail_bulk_return_projection_schema()
        <<"\nfinite="<<(finite_all?1:0)
        <<"\nfeasible_candidate_number_fraction="
        <<(central_seen?central_fraction:0.0)
        <<"\nminimum_feasible_fraction="<<o.minimum_feasible_fraction
        <<"\nfeasibility_gate_required="<<(feasibility_gate_required?1:0)
        <<"\nfeasibility_gate_pass="<<(feasibility_gate_pass?1:0)
        <<"\naudit_read_only_state_unchanged="<<unchanged<<'\n'<<rows.str()
        <<"status="<<(pass?"PASS":"FAIL")<<'\n';
        std::cout<<"status="<<(pass?"PASS":"FAIL")<<'\n';}
    MPI_Finalize(); return pass?0:4;
}
