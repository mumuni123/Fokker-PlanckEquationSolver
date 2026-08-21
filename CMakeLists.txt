cmake_minimum_required(VERSION 3.16)
project(OpenBeamVPFPSolver LANGUAGES CXX)

# Build workflow: configure once with
#   cmake .. -DCMAKE_CXX_COMPILER=mpicxx [-DCMAKE_BUILD_TYPE=Release] [-DBUILD_TESTING=ON]
# then build exactly what is needed with
#   make -j4 fp_solver
#   make -j4 grid_contract_test beam_weight_grid_test open_electrostatic_solver_test
# BUILD_TESTING defaults to ON (include(CTest)), so passing it explicitly is
# optional; test executables are ordinary make targets when testing is enabled.
include(CTest)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED)

set(FP_UPAR_TAIL_CELLS_PER_SIDE 0 CACHE STRING
    "Number of appended cylindrical u_parallel cells on each side")
set(FP_UPAR_EXTENDED_MAX 20.0 CACHE STRING
    "Outer u_parallel bound; must equal the core half-width when no tail cells are appended")

set(VPFP_SOURCES
    src/main_vpfp.cpp
    src/grid.h
    src/parameters.h
    src/species.cpp
    src/beam_pic.cpp
    src/background_tail_pic.cpp
    src/bulk_tail_converter.cpp
    src/tail_subcell_quadrature.cpp
    src/bulk_tail_moment_audit.cpp
    src/bulk_tail_moment_audit_io.cpp
    src/tail_moment_constraint.cpp
    src/tail_population_controller.cpp
    src/tail_bulk_return.cpp
    src/background_tail_collision.cpp
    src/background_tail_nanbu_perez.cpp
    src/background_tail_collision_sde.cpp
    src/hybrid_collision_step.cpp
    src/maxwell.cpp
    src/open_boundary.cpp
    src/open_electrostatic_solver.cpp
    src/joint_phase_space_midpoint.cpp
    src/conservative_ppm_remap.cpp
    src/field_particle_power_audit.cpp
    src/bulk_tail_flux_parcel.cpp
    src/vlasov_split_step.cpp
    src/vpfp_integrator.cpp
    src/vpfp_checkpoint.cpp
    src/collision_coefficients.cpp
    src/cylindrical_fp_collision.cpp
    src/vpfp_diagnostics.cpp)

# fp_solver is now the production open-boundary VPFP-PIC executable.
add_executable(fp_solver ${VPFP_SOURCES})
target_include_directories(fp_solver PRIVATE src)
target_compile_definitions(fp_solver PRIVATE
    FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
    FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
    FP_ENABLE_DEBUG_DIAGNOSTICS=0)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fp_solver PRIVATE -O3 -Wall -Wextra)
elseif(MSVC)
    target_compile_options(fp_solver PRIVATE /O2 /W4)
endif()
target_link_libraries(fp_solver PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX)

# Section 17.11.1: read-only Gate-A audit.  It reuses the production
# checkpoint/grid/tail implementations but deliberately omits main_vpfp.cpp;
# the audit never advances a physical state or writes a checkpoint.
set(VPFP_CORE_SOURCES ${VPFP_SOURCES})
list(REMOVE_ITEM VPFP_CORE_SOURCES src/main_vpfp.cpp)
add_executable(hybrid_checkpoint_gate_audit
    src/hybrid_checkpoint_gate_audit.cpp
    ${VPFP_CORE_SOURCES})
target_include_directories(hybrid_checkpoint_gate_audit PRIVATE src)
target_compile_definitions(hybrid_checkpoint_gate_audit PRIVATE
    FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
    FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
    FP_ENABLE_DEBUG_DIAGNOSTICS=0)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(hybrid_checkpoint_gate_audit PRIVATE -O3 -Wall -Wextra)
elseif(MSVC)
    target_compile_options(hybrid_checkpoint_gate_audit PRIVATE /O2 /W4)
endif()
target_link_libraries(hybrid_checkpoint_gate_audit PRIVATE
    MPI::MPI_CXX OpenMP::OpenMP_CXX)

# Read-only accepted-checkpoint audit used after the 120 fs controller-off
# run. It reports the spatial distribution of Tail particles below selected
# kinetic-energy thresholds and never advances or rewrites the checkpoint.
add_executable(low_energy_tail_checkpoint_audit
    src/low_energy_tail_checkpoint_audit.cpp
    src/low_energy_tail_audit.cpp
    ${VPFP_CORE_SOURCES})
target_include_directories(low_energy_tail_checkpoint_audit PRIVATE src)
target_compile_definitions(low_energy_tail_checkpoint_audit PRIVATE
    FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
    FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
    FP_ENABLE_DEBUG_DIAGNOSTICS=0)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(low_energy_tail_checkpoint_audit PRIVATE
        -O3 -Wall -Wextra)
elseif(MSVC)
    target_compile_options(low_energy_tail_checkpoint_audit PRIVATE /O2 /W4)
endif()
target_link_libraries(low_energy_tail_checkpoint_audit PRIVATE
    MPI::MPI_CXX OpenMP::OpenMP_CXX)

# Phase-1 acceptance tests (section 13.14 and 16.2).  MPI multi-rank cases
# (2 and 5 ranks) are executed explicitly with yhrun on the cluster, not
# registered in plain ctest.
if(BUILD_TESTING)
    add_executable(vpfp_time_control_test tests/vpfp_time_control_test.cpp)
    target_include_directories(vpfp_time_control_test PRIVATE src)
    add_test(NAME vpfp_time_control_test COMMAND vpfp_time_control_test)

    # Section 17.15.7.3: validate the read-only stage-energy audit schema and
    # analyzer without linking or advancing the production solver.
    add_executable(vpfp_stage_energy_audit_test
        tests/vpfp_stage_energy_audit_test.cpp)
    target_link_libraries(vpfp_stage_energy_audit_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_stage_energy_audit_test PRIVATE -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_stage_energy_audit_test PRIVATE /O2 /W4)
    endif()
    add_test(NAME vpfp_stage_energy_audit_test
             COMMAND vpfp_stage_energy_audit_test --case all
                     --result output/vpfp_stage_energy_audit_unit.result)
    set_tests_properties(vpfp_stage_energy_audit_test PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    # Section 10.5.3: direct production-Poisson discrete-work identity.
    add_executable(vpfp_poisson_work_identity_test
        tests/vpfp_poisson_work_identity_test.cpp
        src/open_electrostatic_solver.cpp
        src/maxwell.cpp)
    target_include_directories(vpfp_poisson_work_identity_test PRIVATE src)
    target_link_libraries(vpfp_poisson_work_identity_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(vpfp_poisson_work_identity_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_poisson_work_identity_test PRIVATE -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_poisson_work_identity_test PRIVATE /O2 /W4)
    endif()
    add_test(NAME vpfp_poisson_work_identity_test
             COMMAND vpfp_poisson_work_identity_test --case all
                     --result output/vpfp_poisson_work_identity_unit.result)
    set_tests_properties(vpfp_poisson_work_identity_test PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    # Gate F 10.5.2: checkpoint-backed, read-only dt/dt/2/dt/4 audit.  It
    # links the production core but deliberately excludes main_vpfp.cpp, so
    # no CLI production state, checkpoint, or output path can be mutated.
    add_executable(vpfp_local_defect_checkpoint_audit
        tests/vpfp_local_defect_checkpoint_audit.cpp
        ${VPFP_CORE_SOURCES})
    target_include_directories(vpfp_local_defect_checkpoint_audit PRIVATE src)
    target_compile_definitions(vpfp_local_defect_checkpoint_audit PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_local_defect_checkpoint_audit PRIVATE
            -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_local_defect_checkpoint_audit PRIVATE /O2 /W4)
    endif()
    target_link_libraries(vpfp_local_defect_checkpoint_audit PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)

    # Gate C (section 7.8): read-only discrete-force-work audit.  Drives the
    # production ConservativePpmRemap::advect_u_parallel, BackgroundTailPIC::kick
    # and BeamPIC::last_field_work directly; no production flux is re-derived.
    add_executable(vpfp_force_work_audit_test
        tests/vpfp_force_work_audit_test.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/species.cpp
        src/maxwell.cpp
        src/background_tail_pic.cpp
        src/beam_pic.cpp)
    target_include_directories(vpfp_force_work_audit_test PRIVATE src)
    target_link_libraries(vpfp_force_work_audit_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(vpfp_force_work_audit_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_force_work_audit_test PRIVATE -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_force_work_audit_test PRIVATE /O2 /W4)
    endif()
    add_test(NAME vpfp_force_work_audit_test
             COMMAND vpfp_force_work_audit_test --case all
                     --result output/vpfp_force_work_audit_unit.result)
    set_tests_properties(vpfp_force_work_audit_test PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    # Gate I (section 4.7.2): read-only x-face swept-number exposure test.
    # Links the production ConservativePpmRemap and FieldParticlePowerAudit
    # (no production flux re-derived); MPI multi-rank is not required.
    add_executable(vpfp_x_transport_flux_audit_test
        tests/vpfp_x_transport_flux_audit_test.cpp
        src/conservative_ppm_remap.cpp
        src/field_particle_power_audit.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(vpfp_x_transport_flux_audit_test PRIVATE src)
    target_link_libraries(vpfp_x_transport_flux_audit_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(vpfp_x_transport_flux_audit_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_x_transport_flux_audit_test PRIVATE -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_x_transport_flux_audit_test PRIVATE /O2 /W4)
    endif()
    add_test(NAME vpfp_x_transport_flux_audit_test
             COMMAND vpfp_x_transport_flux_audit_test --case all
                     --result output/vpfp_pairing_gate_i/x_transport_unit.result)
    set_tests_properties(vpfp_x_transport_flux_audit_test PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    # P3.0: drives production x/u remaps and reports their electric-work
    # pairing without changing the production state or re-deriving fluxes.
    add_executable(vpfp_x_u_power_pairing_test
        tests/vpfp_x_u_power_pairing_test.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(vpfp_x_u_power_pairing_test PRIVATE src)
    target_link_libraries(vpfp_x_u_power_pairing_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(vpfp_x_u_power_pairing_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_x_u_power_pairing_test PRIVATE -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_x_u_power_pairing_test PRIVATE /O2 /W4)
    endif()
    add_test(NAME vpfp_x_u_power_pairing_test
             COMMAND vpfp_x_u_power_pairing_test --case all
                     --result output/p3_x_u_pairing_unit.result)
    set_tests_properties(vpfp_x_u_power_pairing_test PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    # P3-V.2.7: independent energy-conjugate velocity table contract test.
    add_executable(vpfp_energy_conjugate_velocity_test
        tests/vpfp_energy_conjugate_velocity_test.cpp)
    target_include_directories(vpfp_energy_conjugate_velocity_test PRIVATE src)
    target_compile_definitions(vpfp_energy_conjugate_velocity_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_energy_conjugate_velocity_test PRIVATE -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_energy_conjugate_velocity_test PRIVATE /O2 /W4)
    endif()
    target_link_libraries(vpfp_energy_conjugate_velocity_test PRIVATE MPI::MPI_CXX)
    add_test(NAME vpfp_energy_conjugate_velocity_test
             COMMAND vpfp_energy_conjugate_velocity_test
                     --result output/p3v2_energy_conjugate_velocity_unit.result)
    set_tests_properties(vpfp_energy_conjugate_velocity_test PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    # J0: pure joint phase-space center-flux and production G/G* audit.  This
    # target does not link the production integrator and cannot alter the
    # legacy Strang/PPM/FCT path.
    add_executable(joint_phase_space_midpoint_unit_test
        tests/joint_phase_space_midpoint_unit_test.cpp
        src/joint_phase_space_midpoint.cpp
        src/open_electrostatic_solver.cpp
        src/maxwell.cpp)
    target_include_directories(joint_phase_space_midpoint_unit_test PRIVATE src)
    target_compile_definitions(joint_phase_space_midpoint_unit_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(joint_phase_space_midpoint_unit_test PRIVATE
            -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(joint_phase_space_midpoint_unit_test PRIVATE /O2 /W4)
    endif()
    target_link_libraries(joint_phase_space_midpoint_unit_test PRIVATE MPI::MPI_CXX)
    add_test(NAME joint_phase_space_midpoint_unit_test
             COMMAND joint_phase_space_midpoint_unit_test --case all
                     --result output/joint_midpoint_j0.result)
    set_tests_properties(joint_phase_space_midpoint_unit_test PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    # J1: background-only joint midpoint smoke.  It links the production core
    # but is explicitly constrained by the executable to beam/tail/collision
    # off; no J1 path is enabled in the default fp_solver mode.
    add_executable(joint_phase_space_midpoint_energy_test
        tests/joint_phase_space_midpoint_energy_test.cpp
        ${VPFP_CORE_SOURCES})
    target_include_directories(joint_phase_space_midpoint_energy_test PRIVATE src)
    target_compile_definitions(joint_phase_space_midpoint_energy_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(joint_phase_space_midpoint_energy_test PRIVATE
            -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(joint_phase_space_midpoint_energy_test PRIVATE /O2 /W4)
    endif()
    target_link_libraries(joint_phase_space_midpoint_energy_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    add_test(NAME joint_phase_space_midpoint_energy_test
             COMMAND joint_phase_space_midpoint_energy_test
                     --case smooth-background
                     --result output/joint_midpoint_j1_unit.result)
    set_tests_properties(joint_phase_space_midpoint_energy_test PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    # Gate I (section 4.7.3): single-rank four-identity pairing test.  Links
    # the production core (minus main) so the production VlasovSplitStep,
    # OpenElectrostaticSolver, BeamPIC, BackgroundTailPIC and the new
    # FieldParticlePowerAudit are driven directly, never re-implemented.
    add_executable(vpfp_field_particle_pairing_test
        tests/vpfp_field_particle_pairing_test.cpp
        ${VPFP_CORE_SOURCES})
    target_include_directories(vpfp_field_particle_pairing_test PRIVATE src)
    target_compile_definitions(vpfp_field_particle_pairing_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_field_particle_pairing_test PRIVATE -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_field_particle_pairing_test PRIVATE /O2 /W4)
    endif()
    target_link_libraries(vpfp_field_particle_pairing_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    add_test(NAME vpfp_field_particle_pairing_test
             COMMAND vpfp_field_particle_pairing_test --case all
                     --result output/vpfp_pairing_gate_i/pairing_unit.result)
    set_tests_properties(vpfp_field_particle_pairing_test PROPERTIES
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME vpfp_field_particle_pairing_injected_faults_test
             COMMAND vpfp_field_particle_pairing_test --case injected-faults
                     --result output/vpfp_pairing_gate_i/injected_faults.result)
    set_tests_properties(vpfp_field_particle_pairing_injected_faults_test
        PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    # Gate I (section 4.7.4): MPI 2/5-rank pairing test.  Built as a plain
    # target; not registered in ctest (executed explicitly with yhrun).
    add_executable(vpfp_field_particle_pairing_mpi_test
        tests/vpfp_field_particle_pairing_mpi_test.cpp
        ${VPFP_CORE_SOURCES})
    target_include_directories(vpfp_field_particle_pairing_mpi_test PRIVATE src)
    target_compile_definitions(vpfp_field_particle_pairing_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_field_particle_pairing_mpi_test PRIVATE -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_field_particle_pairing_mpi_test PRIVATE /O2 /W4)
    endif()
    target_link_libraries(vpfp_field_particle_pairing_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)

    # JC0 (section 3): post-field charge-invariance precheck gate.  Drives the
    # production VpfpIntegrator post-field transaction (C2, collision-face
    # conversion and H10 return) through the read-only JC0 hook; no production
    # collision/conversion/return formula is re-implemented.  Built as a plain
    # target; the single-rank and 2/5-rank MPI runs are executed explicitly
    # with yhrun, not registered in ctest.
    add_executable(vpfp_post_field_charge_invariance_test
        tests/vpfp_post_field_charge_invariance_test.cpp
        ${VPFP_CORE_SOURCES})
    target_include_directories(vpfp_post_field_charge_invariance_test PRIVATE src)
    target_compile_definitions(vpfp_post_field_charge_invariance_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_post_field_charge_invariance_test PRIVATE
            -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_post_field_charge_invariance_test PRIVATE /O2 /W4)
    endif()
    target_link_libraries(vpfp_post_field_charge_invariance_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)

    # JC1 (section 4.11): field-particle trial ownership and candidate-step
    # tests.  Built as a plain target; single-rank runs executed explicitly.
    add_executable(vpfp_field_particle_trial_test
        tests/vpfp_field_particle_trial_test.cpp
        ${VPFP_CORE_SOURCES})
    target_include_directories(vpfp_field_particle_trial_test PRIVATE src)
    target_compile_definitions(vpfp_field_particle_trial_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_field_particle_trial_test PRIVATE
            -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_field_particle_trial_test PRIVATE /O2 /W4)
    endif()
    target_link_libraries(vpfp_field_particle_trial_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)

    # JC2 (section 5.10): MPI face-map test for 2/5-rank runs.
    add_executable(vpfp_field_particle_trial_mpi_test
        tests/vpfp_field_particle_trial_mpi_test.cpp
        ${VPFP_CORE_SOURCES})
    target_include_directories(vpfp_field_particle_trial_mpi_test PRIVATE src)
    target_compile_definitions(vpfp_field_particle_trial_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_field_particle_trial_mpi_test PRIVATE
            -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_field_particle_trial_mpi_test PRIVATE /O2 /W4)
    endif()
    target_link_libraries(vpfp_field_particle_trial_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)

    # JC3 (section 6.11): post-field integration gate.  Combines the
    # converged Picard trial with the one-shot C2/conversion/return path and
    # the atomic commit, validating the JC0 combined-charge gate.
    add_executable(vpfp_field_particle_post_field_test
        tests/vpfp_field_particle_post_field_test.cpp
        ${VPFP_CORE_SOURCES})
    target_include_directories(vpfp_field_particle_post_field_test PRIVATE src)
    target_compile_definitions(vpfp_field_particle_post_field_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(vpfp_field_particle_post_field_test PRIVATE
            -O3 -Wall -Wextra)
    elseif(MSVC)
        target_compile_options(vpfp_field_particle_post_field_test PRIVATE /O2 /W4)
    endif()
    target_link_libraries(vpfp_field_particle_post_field_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)

    # Gate I (section 4.7.7 item 5): the Python analyzer must byte-compile
    # before the cluster runs.  Registered as a plain ctest; the 2/5-rank
    # explicit MPI execution remains outside ctest (section 4.7.7 item 2).
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME analyze_vpfp_field_particle_pairing_py_compile
                 COMMAND ${Python3_EXECUTABLE} -m py_compile
                         ${CMAKE_SOURCE_DIR}/tools/analyze_vpfp_field_particle_pairing.py)
        set_tests_properties(analyze_vpfp_field_particle_pairing_py_compile
            PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()

    add_executable(low_energy_tail_audit_test
        tests/low_energy_tail_audit_test.cpp
        src/low_energy_tail_audit.cpp)
    target_include_directories(low_energy_tail_audit_test PRIVATE src)
    target_link_libraries(low_energy_tail_audit_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    add_test(NAME low_energy_tail_audit_test
             COMMAND low_energy_tail_audit_test)

    add_executable(grid_contract_test tests/grid_contract_test.cpp)
    target_include_directories(grid_contract_test PRIVATE src)
    target_link_libraries(grid_contract_test PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(grid_contract_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME grid_contract_test
             COMMAND grid_contract_test
                     --expected-nx 8000
                     --expected-length-um 40
                     --expected-dx-um 0.005)

    add_executable(beam_weight_grid_test
        tests/beam_weight_grid_test.cpp
        src/beam_pic.cpp
        src/maxwell.cpp)
    target_include_directories(beam_weight_grid_test PRIVATE src)
    target_link_libraries(beam_weight_grid_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(beam_weight_grid_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME beam_weight_grid_test
             COMMAND beam_weight_grid_test
                     --dx-um 0.005
                     --particles-per-cell 1000)

    add_executable(open_electrostatic_solver_test
        tests/open_electrostatic_solver_test.cpp
        src/open_electrostatic_solver.cpp
        src/maxwell.cpp)
    target_include_directories(open_electrostatic_solver_test PRIVATE src)
    target_link_libraries(open_electrostatic_solver_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(open_electrostatic_solver_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME open_electrostatic_solver_test
             COMMAND open_electrostatic_solver_test
                     --case all
                     --field-boundary dirichlet-phi)

    add_executable(conservative_x_remap_test
        tests/conservative_x_remap_test.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(conservative_x_remap_test PRIVATE src)
    target_link_libraries(conservative_x_remap_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(conservative_x_remap_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME conservative_x_remap_test
             COMMAND conservative_x_remap_test --case all)

    add_executable(open_boundary_remap_test
        tests/open_boundary_remap_test.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(open_boundary_remap_test PRIVATE src)
    target_link_libraries(open_boundary_remap_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(open_boundary_remap_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME open_boundary_remap_test
             COMMAND open_boundary_remap_test
                     --case reservoir-and-absorbing)

    add_executable(conservative_x_remap_mpi_test
        tests/conservative_x_remap_mpi_test.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(conservative_x_remap_mpi_test PRIVATE src)
    target_link_libraries(conservative_x_remap_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(conservative_x_remap_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME conservative_x_remap_mpi_test
             COMMAND conservative_x_remap_mpi_test
                     --case crossing-and-multicell)

    add_executable(conservative_upar_remap_test
        tests/conservative_upar_remap_test.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(conservative_upar_remap_test PRIVATE src)
    target_link_libraries(conservative_upar_remap_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(conservative_upar_remap_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME conservative_upar_remap_test
             COMMAND conservative_upar_remap_test --case all)

    add_executable(vlasov_poisson_order_test
        tests/vlasov_poisson_order_test.cpp
        src/vlasov_split_step.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/open_electrostatic_solver.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/bulk_tail_converter.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/tail_population_controller.cpp
        src/background_tail_collision.cpp
        src/background_tail_nanbu_perez.cpp
        src/background_tail_collision_sde.cpp
        src/hybrid_collision_step.cpp
        src/maxwell.cpp
        src/collision_coefficients.cpp
        src/cylindrical_fp_collision.cpp
        src/vpfp_integrator.cpp
        src/field_particle_power_audit.cpp
        src/tail_bulk_return.cpp)
    target_include_directories(vlasov_poisson_order_test PRIVATE src)
    target_link_libraries(vlasov_poisson_order_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(vlasov_poisson_order_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME vlasov_poisson_order_test
             COMMAND vlasov_poisson_order_test --case free-stream-reservoir)

    add_executable(beam_leapfrog_test
        tests/beam_leapfrog_test.cpp
        src/beam_pic.cpp
        src/maxwell.cpp)
    target_include_directories(beam_leapfrog_test PRIVATE src)
    target_link_libraries(beam_leapfrog_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(beam_leapfrog_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME beam_leapfrog_test
             COMMAND beam_leapfrog_test --case zero-and-constant-field)

    add_executable(beam_injection_event_test
        tests/beam_injection_event_test.cpp
        src/beam_pic.cpp
        src/maxwell.cpp)
    target_include_directories(beam_injection_event_test PRIVATE src)
    target_link_libraries(beam_injection_event_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(beam_injection_event_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME beam_injection_event_test
             COMMAND beam_injection_event_test
                     --case before-and-after-midpoint)

    add_executable(beam_open_boundary_mpi_test
        tests/beam_open_boundary_mpi_test.cpp
        src/vlasov_split_step.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/open_electrostatic_solver.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/bulk_tail_converter.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/tail_population_controller.cpp
        src/background_tail_collision.cpp
        src/background_tail_nanbu_perez.cpp
        src/background_tail_collision_sde.cpp
        src/hybrid_collision_step.cpp
        src/maxwell.cpp
        src/collision_coefficients.cpp
        src/cylindrical_fp_collision.cpp
        src/vpfp_integrator.cpp
        src/field_particle_power_audit.cpp
        src/tail_bulk_return.cpp)
    target_include_directories(beam_open_boundary_mpi_test PRIVATE src)
    target_link_libraries(beam_open_boundary_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(beam_open_boundary_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME beam_open_boundary_mpi_test
             COMMAND beam_open_boundary_mpi_test
                     --case migration-left-right-poisson)

    # Stage H1 targets (section 15 H1, section 17.3): tail PIC
    # infrastructure.  The single-rank tests run plainly; the MPI test is
    # executed explicitly with yhrun -n 1 / -n 2 / -n 5 on the cluster.
    add_executable(background_tail_pusher_test
        tests/background_tail_pusher_test.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(background_tail_pusher_test PRIVATE src)
    target_link_libraries(background_tail_pusher_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(background_tail_pusher_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME background_tail_pusher_test
             COMMAND background_tail_pusher_test
                     --case zero-and-constant-field)

    add_executable(background_tail_deposition_test
        tests/background_tail_deposition_test.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(background_tail_deposition_test PRIVATE src)
    target_link_libraries(background_tail_deposition_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(background_tail_deposition_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME background_tail_deposition_test
             COMMAND background_tail_deposition_test --case all)

    add_executable(background_tail_shape_difference_test
        tests/background_tail_shape_difference_test.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(background_tail_shape_difference_test PRIVATE src)
    target_link_libraries(background_tail_shape_difference_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(background_tail_shape_difference_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME background_tail_shape_difference_test
             COMMAND background_tail_shape_difference_test)

    add_executable(background_tail_stagger_test
        tests/background_tail_stagger_test.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(background_tail_stagger_test PRIVATE src)
    target_link_libraries(background_tail_stagger_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(background_tail_stagger_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME background_tail_stagger_test
             COMMAND background_tail_stagger_test)

    add_executable(background_tail_open_boundary_mpi_test
        tests/background_tail_open_boundary_mpi_test.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(background_tail_open_boundary_mpi_test PRIVATE src)
    target_link_libraries(background_tail_open_boundary_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(background_tail_open_boundary_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME background_tail_open_boundary_mpi_test
             COMMAND background_tail_open_boundary_mpi_test
                     --case migration-consistency)

    # Stage H2 targets (sections 15 H2 and 17.4): standalone conservative
    # bulk-to-tail converter.  All run plainly; the multibin test is also
    # executed explicitly with yhrun -n 2 / -n 5 on the cluster for the
    # rank-consistency acceptance.
    add_executable(bulk_tail_single_cell_test
        tests/bulk_tail_single_cell_test.cpp
        src/bulk_tail_converter.cpp
        src/bulk_tail_flux_parcel.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_single_cell_test PRIVATE src)
    target_link_libraries(bulk_tail_single_cell_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_single_cell_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_single_cell_test
             COMMAND bulk_tail_single_cell_test)

    add_executable(bulk_tail_multibin_test
        tests/bulk_tail_multibin_test.cpp
        src/bulk_tail_converter.cpp
        src/bulk_tail_flux_parcel.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_multibin_test PRIVATE src)
    target_link_libraries(bulk_tail_multibin_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_multibin_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_multibin_test
             COMMAND bulk_tail_multibin_test)

    add_executable(bulk_tail_poisson_invariance_test
        tests/bulk_tail_poisson_invariance_test.cpp
        src/bulk_tail_converter.cpp
        src/bulk_tail_flux_parcel.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/open_electrostatic_solver.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_poisson_invariance_test PRIVATE src)
    target_link_libraries(bulk_tail_poisson_invariance_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_poisson_invariance_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_poisson_invariance_test
             COMMAND bulk_tail_poisson_invariance_test)

    add_executable(bulk_tail_transaction_test
        tests/bulk_tail_transaction_test.cpp
        src/bulk_tail_converter.cpp
        src/bulk_tail_flux_parcel.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_transaction_test PRIVATE src)
    target_link_libraries(bulk_tail_transaction_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_transaction_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_transaction_test
             COMMAND bulk_tail_transaction_test)

    # Section 7.11.2/7.11.3 threshold-interface audit targets (runbook
    # 17.10.1).  The support audit is read-only; the interface tests select
    # the loading policy through BulkTailConverter::set_loading_policy (no
    # compile macros or globals).
    add_executable(hybrid_threshold_support_audit_test
        tests/hybrid_threshold_support_audit_test.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(hybrid_threshold_support_audit_test PRIVATE src)
    target_link_libraries(hybrid_threshold_support_audit_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(hybrid_threshold_support_audit_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME hybrid_threshold_support_audit_test
             COMMAND hybrid_threshold_support_audit_test
                     --convert-energy-mev 6.0)

    add_executable(bulk_tail_threshold_interface_test
        tests/bulk_tail_threshold_interface_test.cpp
        src/bulk_tail_converter.cpp
        src/bulk_tail_flux_parcel.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_threshold_interface_test PRIVATE src)
    target_link_libraries(bulk_tail_threshold_interface_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_threshold_interface_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_threshold_interface_test
             COMMAND bulk_tail_threshold_interface_test
                     --case all --policy all --bin-width-mev 0.1)

    add_executable(bulk_tail_threshold_interface_mpi_test
        tests/bulk_tail_threshold_interface_mpi_test.cpp
        src/bulk_tail_converter.cpp
        src/bulk_tail_flux_parcel.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_threshold_interface_mpi_test PRIVATE src)
    target_link_libraries(bulk_tail_threshold_interface_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_threshold_interface_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_threshold_interface_mpi_test
             COMMAND bulk_tail_threshold_interface_mpi_test
                     --case all --policy threshold-aware --bin-width-mev 0.1)

    # Section 7.11.9: read-only finite-volume spectrum reference.  This is
    # deliberately not linked into the converter or integrator; it audits
    # physical velocity-cell support before any subcell loading is allowed.
    add_executable(bulk_tail_cell_volume_spectrum_test
        tests/bulk_tail_cell_volume_spectrum_test.cpp
        src/tail_subcell_quadrature.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_cell_volume_spectrum_test PRIVATE src)
    target_link_libraries(bulk_tail_cell_volume_spectrum_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_cell_volume_spectrum_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_cell_volume_spectrum_test
             COMMAND bulk_tail_cell_volume_spectrum_test
                     --convert-energy-mev 6.0
                     --uperp-rings 0,1,2,3,4,5,6,7
                     --bin-widths-mev 0.05,0.1,0.2)

    # Gate-A v5 geometry regression: exact cylindrical cell/energy-shell
    # overlap used by the checkpoint spectrum audit.
    add_executable(tail_energy_shell_overlap_test
        tests/tail_energy_shell_overlap_test.cpp
        src/tail_subcell_quadrature.cpp)
    target_include_directories(tail_energy_shell_overlap_test PRIVATE src)
    target_compile_definitions(tail_energy_shell_overlap_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME tail_energy_shell_overlap_test
             COMMAND tail_energy_shell_overlap_test)

    # Section 7.11.10 candidate.  It explicitly enables the conservative
    # subcell loader only inside this test; the production converter remains
    # gated until the 7.11.9 reference has been reviewed.
    add_executable(bulk_tail_subcell_loading_test
        tests/bulk_tail_subcell_loading_test.cpp
        src/bulk_tail_converter.cpp
        src/bulk_tail_flux_parcel.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_subcell_loading_test PRIVATE src)
    target_link_libraries(bulk_tail_subcell_loading_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_subcell_loading_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)

    add_executable(bulk_tail_real_moment_audit_test
        tests/bulk_tail_real_moment_audit_test.cpp
        src/bulk_tail_moment_audit.cpp src/bulk_tail_moment_audit_io.cpp
        src/tail_subcell_quadrature.cpp
        src/tail_moment_constraint.cpp src/species.cpp src/maxwell.cpp)
    target_include_directories(bulk_tail_real_moment_audit_test PRIVATE src)
    target_link_libraries(bulk_tail_real_moment_audit_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    add_test(NAME bulk_tail_real_moment_audit_test
             COMMAND bulk_tail_real_moment_audit_test)

    add_executable(bulk_tail_real_moment_audit_mpi_test
        tests/bulk_tail_real_moment_audit_mpi_test.cpp
        src/bulk_tail_moment_audit.cpp src/bulk_tail_moment_audit_io.cpp
        src/tail_subcell_quadrature.cpp
        src/tail_moment_constraint.cpp src/species.cpp src/maxwell.cpp)
    target_include_directories(bulk_tail_real_moment_audit_mpi_test PRIVATE src)
    target_link_libraries(bulk_tail_real_moment_audit_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    add_test(NAME bulk_tail_real_moment_audit_mpi_test
             COMMAND bulk_tail_real_moment_audit_mpi_test)

    # Section 7.11.16C / 16B: fixed-(192,64) local velocity-grid offline
    # prototype.  The module is independent of the production grid default:
    # tail_interface_grid_design.cpp never calls CylindricalVelocityGrid::init
    # for the candidate (it reuses the public geometry/moment builders) and
    # the tests only link the production module plus its shared quadrature /
    # moment-constraint dependencies.
    add_executable(tail_interface_grid_design_test
        tests/tail_interface_grid_design_test.cpp
        src/tail_interface_grid_design.cpp
        src/tail_subcell_quadrature.cpp
        src/tail_moment_constraint.cpp)
    target_include_directories(tail_interface_grid_design_test PRIVATE src)
    target_link_libraries(tail_interface_grid_design_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    add_test(NAME tail_interface_grid_design_test
             COMMAND tail_interface_grid_design_test --case all)

    add_executable(tail_interface_grid_replay_test
        tests/tail_interface_grid_replay_test.cpp
        src/tail_interface_grid_design.cpp
        src/tail_subcell_quadrature.cpp
        src/tail_moment_constraint.cpp)
    target_include_directories(tail_interface_grid_replay_test PRIVATE src)
    target_link_libraries(tail_interface_grid_replay_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    add_test(NAME tail_interface_grid_replay_test
             COMMAND tail_interface_grid_replay_test
                     --input
                     ${CMAKE_CURRENT_SOURCE_DIR}/output/hybrid_h9_real_moment_audit_12fs/bulk_tail_moment_audit_velocity_histogram.dat
                     --profiles G0,Gx,Gp,G2
                     --ax-values 0.5,1.0,2.0
                     --aperp-values 0.5,1.0,2.0
                     --sigma-x-cells 1,2,4
                     --sigma-perp-cells 1,2,4)

    add_executable(tail_interface_grid_replay_mpi_test
        tests/tail_interface_grid_replay_mpi_test.cpp
        src/tail_interface_grid_design.cpp
        src/tail_subcell_quadrature.cpp
        src/tail_moment_constraint.cpp)
    target_include_directories(tail_interface_grid_replay_mpi_test PRIVATE src)
    target_link_libraries(tail_interface_grid_replay_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    add_test(NAME tail_interface_grid_replay_mpi_test
             COMMAND tail_interface_grid_replay_mpi_test
                     --input
                     ${CMAKE_CURRENT_SOURCE_DIR}/output/hybrid_h9_real_moment_audit_12fs/bulk_tail_moment_audit_velocity_histogram.dat
                     --profiles G0,Gx,Gp,G2
                     --ax-values 0.5,1.0,2.0
                     --aperp-values 0.5,1.0,2.0
                     --sigma-x-cells 1,2,4
                     --sigma-perp-cells 1,2,4)

    # Section 7.11.10 feasibility audit: separates solver correctness from
    # incompatibility between a cell-centre six-moment target and a spread
    # nonnegative subcell representation.
    add_executable(bulk_tail_subcell_feasibility_test
        tests/bulk_tail_subcell_feasibility_test.cpp
        src/tail_subcell_quadrature.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_subcell_feasibility_test PRIVATE src)
    target_link_libraries(bulk_tail_subcell_feasibility_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_subcell_feasibility_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_subcell_feasibility_test
             COMMAND bulk_tail_subcell_feasibility_test --case all)

    # Stage H3 target (sections 15 H3 and 17.5): no-beam hybrid
    # Vlasov-Poisson integration driving the production VpfpIntegrator with
    # the tail PIC + converter wired in.  Single rank.
    add_executable(hybrid_no_beam_test
        tests/hybrid_no_beam_test.cpp
        src/vpfp_integrator.cpp
        src/field_particle_power_audit.cpp
        src/tail_bulk_return.cpp
        src/vlasov_split_step.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/open_electrostatic_solver.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/bulk_tail_converter.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/tail_population_controller.cpp
        src/background_tail_collision.cpp
        src/background_tail_nanbu_perez.cpp
        src/background_tail_collision_sde.cpp
        src/hybrid_collision_step.cpp
        src/maxwell.cpp
        src/collision_coefficients.cpp
        src/cylindrical_fp_collision.cpp)
    target_include_directories(hybrid_no_beam_test PRIVATE src)
    target_link_libraries(hybrid_no_beam_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(hybrid_no_beam_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME hybrid_no_beam_test
             COMMAND hybrid_no_beam_test --case all)

    # Stage H4 target (sections 15 H4 and 17.6): Beam + collisionless hybrid
    # short integration.  Single rank.
    add_executable(hybrid_beam_short_test
        tests/hybrid_beam_short_test.cpp
        src/vpfp_integrator.cpp
        src/field_particle_power_audit.cpp
        src/tail_bulk_return.cpp
        src/vlasov_split_step.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/open_electrostatic_solver.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/bulk_tail_converter.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/tail_population_controller.cpp
        src/background_tail_collision.cpp
        src/background_tail_nanbu_perez.cpp
        src/background_tail_collision_sde.cpp
        src/hybrid_collision_step.cpp
        src/maxwell.cpp
        src/collision_coefficients.cpp
        src/cylindrical_fp_collision.cpp)
    target_include_directories(hybrid_beam_short_test PRIVATE src)
    target_link_libraries(hybrid_beam_short_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(hybrid_beam_short_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME hybrid_beam_short_test
             COMMAND hybrid_beam_short_test --case all)

    # H9 multi-rank control-flow gate: a collision failure on only one rank
    # must be reduced to a common rejected step before any rank returns.
    add_executable(hybrid_collision_failure_consensus_mpi_test
        tests/hybrid_collision_failure_consensus_mpi_test.cpp
        src/vpfp_integrator.cpp
        src/field_particle_power_audit.cpp
        src/tail_bulk_return.cpp
        src/vlasov_split_step.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/open_electrostatic_solver.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/bulk_tail_converter.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/tail_population_controller.cpp
        src/background_tail_collision.cpp
        src/background_tail_nanbu_perez.cpp
        src/background_tail_collision_sde.cpp
        src/hybrid_collision_step.cpp
        src/maxwell.cpp
        src/collision_coefficients.cpp
        src/cylindrical_fp_collision.cpp)
    target_include_directories(hybrid_collision_failure_consensus_mpi_test
        PRIVATE src)
    target_link_libraries(hybrid_collision_failure_consensus_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(hybrid_collision_failure_consensus_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)

    # Stage H5 targets (sections 15 H5 and 17.7): TailPopulationController
    # (section 7.10) plus the shared TailMomentConstraint module.  The
    # single-rank test runs plainly; the MPI test is executed explicitly
    # with yhrun -n 1 / -n 2 / -n 5 on the cluster.
    add_executable(tail_population_controller_test
        tests/tail_population_controller_test.cpp
        src/tail_population_controller.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(tail_population_controller_test PRIVATE src)
    target_link_libraries(tail_population_controller_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(tail_population_controller_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME tail_population_controller_test
             COMMAND tail_population_controller_test --case all)

    add_executable(tail_population_controller_mpi_test
        tests/tail_population_controller_mpi_test.cpp
        src/tail_population_controller.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(tail_population_controller_mpi_test PRIVATE src)
    target_link_libraries(tail_population_controller_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(tail_population_controller_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME tail_population_controller_mpi_test
             COMMAND tail_population_controller_mpi_test)

    # Stage H10 (section 16): tail-to-bulk hysteresis, conservative
    # projection, transaction, disabled-mode and MPI seam acceptance gates.
    set(TAIL_BULK_RETURN_TEST_SOURCES
        src/tail_bulk_return.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    foreach(test_name
            tail_bulk_return_none_regression_test
            tail_bulk_return_hysteresis_test
            tail_bulk_return_parameter_switch_test
            tail_bulk_return_projection_test
            tail_bulk_return_transaction_test
            tail_bulk_return_mpi_test)
        add_executable(${test_name}
            tests/${test_name}.cpp
            ${TAIL_BULK_RETURN_TEST_SOURCES})
        target_include_directories(${test_name} PRIVATE src tests)
        target_link_libraries(${test_name} PRIVATE
            MPI::MPI_CXX OpenMP::OpenMP_CXX)
        target_compile_definitions(${test_name} PRIVATE
            FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
            FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
            FP_ENABLE_DEBUG_DIAGNOSTICS=0)
        add_test(NAME ${test_name} COMMAND ${test_name})
    endforeach()

    add_executable(tail_bulk_return_checkpoint_test
        tests/tail_bulk_return_checkpoint_test.cpp
        src/vpfp_checkpoint.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(tail_bulk_return_checkpoint_test PRIVATE src tests)
    target_link_libraries(tail_bulk_return_checkpoint_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(tail_bulk_return_checkpoint_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME tail_bulk_return_checkpoint_test
             COMMAND tail_bulk_return_checkpoint_test)

    add_executable(tail_bulk_return_checkpoint_audit
        src/tail_bulk_return_checkpoint_audit.cpp
        src/tail_bulk_return.cpp
        src/tail_moment_constraint.cpp
        src/vpfp_checkpoint.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(tail_bulk_return_checkpoint_audit PRIVATE src)
    target_link_libraries(tail_bulk_return_checkpoint_audit PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(tail_bulk_return_checkpoint_audit PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)

    # Stage H6 targets (sections 12.3 and 15 H6): checkpoint schema-v2
    # round trip and the one-step checkpoint/restart equivalence gate.
    add_executable(checkpoint_roundtrip_test
        tests/checkpoint_roundtrip_test.cpp
        src/vpfp_checkpoint.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(checkpoint_roundtrip_test PRIVATE src)
    target_link_libraries(checkpoint_roundtrip_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(checkpoint_roundtrip_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME checkpoint_roundtrip_test
             COMMAND checkpoint_roundtrip_test)

    add_executable(checkpoint_restart_equivalence_test
        tests/checkpoint_restart_equivalence_test.cpp
        src/vpfp_checkpoint.cpp
        src/vpfp_integrator.cpp
        src/field_particle_power_audit.cpp
        src/tail_bulk_return.cpp
        src/vlasov_split_step.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/open_electrostatic_solver.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/bulk_tail_converter.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/tail_population_controller.cpp
        src/background_tail_collision.cpp
        src/background_tail_nanbu_perez.cpp
        src/background_tail_collision_sde.cpp
        src/hybrid_collision_step.cpp
        src/maxwell.cpp
        src/collision_coefficients.cpp
        src/cylindrical_fp_collision.cpp)
    target_include_directories(checkpoint_restart_equivalence_test PRIVATE src)
    target_link_libraries(checkpoint_restart_equivalence_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(checkpoint_restart_equivalence_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME checkpoint_restart_equivalence_test
             COMMAND checkpoint_restart_equivalence_test)

    # Stage H7 target (sections 15 H7 and 17.9): bulk cylindrical
    # Fokker-Planck collision with velocity-dependent coefficients,
    # non-diagonal diffusion support and collision-reservoir accounting.
    add_executable(cylindrical_fp_collision_test
        tests/cylindrical_fp_collision_test.cpp
        src/cylindrical_fp_collision.cpp
        src/collision_coefficients.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(cylindrical_fp_collision_test PRIVATE src)
    target_link_libraries(cylindrical_fp_collision_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(cylindrical_fp_collision_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME cylindrical_fp_collision_test
             COMMAND cylindrical_fp_collision_test --case all)

    # Stage H8 targets (sections 15 H8 and 17.9): tail collision backends
    # (Nanbu-Perez, SDE), the unified tail collision entry and the hybrid
    # bulk+tail half-step.
    set(H8_TAIL_SOURCES
        src/background_tail_collision.cpp
        src/background_tail_nanbu_perez.cpp
        src/background_tail_collision_sde.cpp
        src/collision_coefficients.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    add_executable(background_tail_nanbu_perez_test
        tests/background_tail_nanbu_perez_test.cpp ${H8_TAIL_SOURCES})
    target_include_directories(background_tail_nanbu_perez_test PRIVATE src)
    target_link_libraries(background_tail_nanbu_perez_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(background_tail_nanbu_perez_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME background_tail_nanbu_perez_test
             COMMAND background_tail_nanbu_perez_test)

    add_executable(background_tail_collision_moments_test
        tests/background_tail_collision_moments_test.cpp ${H8_TAIL_SOURCES})
    target_include_directories(background_tail_collision_moments_test PRIVATE src)
    target_link_libraries(background_tail_collision_moments_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(background_tail_collision_moments_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME background_tail_collision_moments_test
             COMMAND background_tail_collision_moments_test)

    add_executable(background_tail_collision_equilibrium_test
        tests/background_tail_collision_equilibrium_test.cpp ${H8_TAIL_SOURCES})
    target_include_directories(background_tail_collision_equilibrium_test
        PRIVATE src)
    target_link_libraries(background_tail_collision_equilibrium_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(background_tail_collision_equilibrium_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME background_tail_collision_equilibrium_test
             COMMAND background_tail_collision_equilibrium_test)

    add_executable(background_tail_collision_weight_test
        tests/background_tail_collision_weight_test.cpp ${H8_TAIL_SOURCES})
    target_include_directories(background_tail_collision_weight_test PRIVATE src)
    target_link_libraries(background_tail_collision_weight_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(background_tail_collision_weight_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME background_tail_collision_weight_test
             COMMAND background_tail_collision_weight_test)

    add_executable(tail_collision_isotropisation_test
        tests/tail_collision_isotropisation_test.cpp ${H8_TAIL_SOURCES})
    target_include_directories(tail_collision_isotropisation_test PRIVATE src)
    target_link_libraries(tail_collision_isotropisation_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(tail_collision_isotropisation_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME tail_collision_isotropisation_test
             COMMAND tail_collision_isotropisation_test)

    add_executable(tail_collision_two_population_test
        tests/tail_collision_two_population_test.cpp ${H8_TAIL_SOURCES})
    target_include_directories(tail_collision_two_population_test PRIVATE src)
    target_link_libraries(tail_collision_two_population_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(tail_collision_two_population_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME tail_collision_two_population_test
             COMMAND tail_collision_two_population_test)

    set(H8_HYBRID_SOURCES
        src/hybrid_collision_step.cpp
        src/cylindrical_fp_collision.cpp
        src/background_tail_nanbu_perez.cpp
        src/background_tail_collision_sde.cpp
        src/background_tail_collision.cpp
        src/collision_coefficients.cpp
        src/background_tail_pic.cpp
        src/species.cpp
        src/maxwell.cpp)
    add_executable(hybrid_collision_pair_balance_test
        tests/hybrid_collision_pair_balance_test.cpp ${H8_HYBRID_SOURCES})
    target_include_directories(hybrid_collision_pair_balance_test PRIVATE src)
    target_link_libraries(hybrid_collision_pair_balance_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(hybrid_collision_pair_balance_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME hybrid_collision_pair_balance_test
             COMMAND hybrid_collision_pair_balance_test)

    add_executable(hybrid_collision_transaction_test
        tests/hybrid_collision_transaction_test.cpp ${H8_HYBRID_SOURCES})
    target_include_directories(hybrid_collision_transaction_test PRIVATE src)
    target_link_libraries(hybrid_collision_transaction_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(hybrid_collision_transaction_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME hybrid_collision_transaction_test
             COMMAND hybrid_collision_transaction_test)

    # Section 7.11.17 flux-interface contracts.  These tests are deliberately
    # small and independent of the production main: topology, parcel moment
    # construction, and the one-way u_parallel sink are checked separately.
    add_executable(bulk_tail_flux_interface_test
        tests/bulk_tail_flux_interface_test.cpp
        src/tail_subcell_quadrature.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_flux_interface_test PRIVATE src)
    target_link_libraries(bulk_tail_flux_interface_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_flux_interface_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_flux_interface_test
             COMMAND bulk_tail_flux_interface_test)

    add_executable(upar_flux_parcel_test
        tests/upar_flux_parcel_test.cpp
        src/bulk_tail_flux_parcel.cpp)
    target_include_directories(upar_flux_parcel_test PRIVATE src)
    target_link_libraries(upar_flux_parcel_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(upar_flux_parcel_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME upar_flux_parcel_test COMMAND upar_flux_parcel_test)

    add_executable(upar_flux_sink_test
        tests/upar_flux_sink_test.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/species.cpp
        src/maxwell.cpp)
    target_include_directories(upar_flux_sink_test PRIVATE src)
    target_link_libraries(upar_flux_sink_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(upar_flux_sink_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME upar_flux_sink_test COMMAND upar_flux_sink_test)

    # Section 7.11.17.5: parcel-to-tail loader and transaction contracts.
    # These targets deliberately share the production converter, moment
    # compressor and tail state implementation; the tests do not duplicate
    # the loader math.
    set(BULK_TAIL_FLUX_LOADER_SOURCES
        src/bulk_tail_converter.cpp
        src/bulk_tail_flux_parcel.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/species.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)

    add_executable(bulk_tail_flux_loader_test
        tests/bulk_tail_flux_loader_test.cpp
        ${BULK_TAIL_FLUX_LOADER_SOURCES})
    target_include_directories(bulk_tail_flux_loader_test PRIVATE src)
    target_link_libraries(bulk_tail_flux_loader_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_flux_loader_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_flux_loader_test COMMAND bulk_tail_flux_loader_test)

    # Regression for the production six-moment compressor.  A large local
    # flux group must reduce in bounded work instead of rebuilding an RREF
    # over the complete active set for every eliminated support.
    add_executable(tail_moment_constraint_large_group_test
        tests/tail_moment_constraint_large_group_test.cpp
        src/tail_moment_constraint.cpp)
    target_include_directories(tail_moment_constraint_large_group_test PRIVATE src)
    target_link_libraries(tail_moment_constraint_large_group_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(tail_moment_constraint_large_group_test PRIVATE
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME tail_moment_constraint_large_group_test
             COMMAND tail_moment_constraint_large_group_test)

    # Section 17.10.2.6: a material source must cross a real u_parallel
    # interface through the production remap before the production converter
    # is allowed to create tail particles.
    add_executable(bulk_tail_flux_source_trigger_test
        tests/bulk_tail_flux_source_trigger_test.cpp
        src/conservative_ppm_remap.cpp
        src/open_boundary.cpp
        ${BULK_TAIL_FLUX_LOADER_SOURCES})
    target_include_directories(bulk_tail_flux_source_trigger_test PRIVATE src)
    target_link_libraries(bulk_tail_flux_source_trigger_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_flux_source_trigger_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_flux_source_trigger_test
             COMMAND bulk_tail_flux_source_trigger_test)

    add_executable(bulk_tail_flux_transaction_test
        tests/bulk_tail_flux_transaction_test.cpp
        ${BULK_TAIL_FLUX_LOADER_SOURCES})
    target_include_directories(bulk_tail_flux_transaction_test PRIVATE src)
    target_link_libraries(bulk_tail_flux_transaction_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_flux_transaction_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_flux_transaction_test
             COMMAND bulk_tail_flux_transaction_test)

    add_executable(bulk_tail_flux_loader_mpi_test
        tests/bulk_tail_flux_loader_mpi_test.cpp
        ${BULK_TAIL_FLUX_LOADER_SOURCES})
    target_include_directories(bulk_tail_flux_loader_mpi_test PRIVATE src)
    target_link_libraries(bulk_tail_flux_loader_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_flux_loader_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_flux_loader_mpi_test
             COMMAND bulk_tail_flux_loader_mpi_test)

    # Section 7.11.17.8 / 17D: checkpoint metadata, particle IDs and
    # cumulative ledgers must survive both a single-rank round trip and a
    # repartitioned MPI restart.
    add_executable(bulk_tail_flux_checkpoint_test
        tests/bulk_tail_flux_checkpoint_test.cpp
        src/vpfp_checkpoint.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_flux_checkpoint_test PRIVATE src)
    target_link_libraries(bulk_tail_flux_checkpoint_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_flux_checkpoint_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_flux_checkpoint_test
             COMMAND bulk_tail_flux_checkpoint_test)

    add_executable(bulk_tail_flux_restart_mpi_test
        tests/checkpoint_restart_equivalence_test.cpp
        src/vpfp_checkpoint.cpp
        src/vpfp_integrator.cpp
        src/field_particle_power_audit.cpp
        src/tail_bulk_return.cpp
        src/vlasov_split_step.cpp
        src/conservative_ppm_remap.cpp
        src/bulk_tail_flux_parcel.cpp
        src/open_boundary.cpp
        src/open_electrostatic_solver.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/bulk_tail_converter.cpp
        src/tail_subcell_quadrature.cpp
        src/bulk_tail_moment_audit.cpp
        src/tail_moment_constraint.cpp
        src/tail_population_controller.cpp
        src/hybrid_collision_step.cpp
        src/background_tail_collision_sde.cpp
        src/background_tail_nanbu_perez.cpp
        src/maxwell.cpp
        src/collision_coefficients.cpp
        src/cylindrical_fp_collision.cpp)
    target_include_directories(bulk_tail_flux_restart_mpi_test PRIVATE src)
    target_link_libraries(bulk_tail_flux_restart_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_flux_restart_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_flux_restart_mpi_test
             COMMAND bulk_tail_flux_restart_mpi_test)

    # Metadata/ID/ledger round-trip at a fixed global state.  This is kept as
    # a separate target from the full production one-step equivalence above:
    # the latter exercises the actual Vlasov/field continuation, while this
    # small test isolates repartition-independent checkpoint metadata and is
    # inexpensive to run at 1/2/5 ranks.
    add_executable(bulk_tail_flux_checkpoint_mpi_test
        tests/bulk_tail_flux_restart_mpi_test.cpp
        src/vpfp_checkpoint.cpp
        src/species.cpp
        src/beam_pic.cpp
        src/background_tail_pic.cpp
        src/maxwell.cpp)
    target_include_directories(bulk_tail_flux_checkpoint_mpi_test PRIVATE src)
    target_link_libraries(bulk_tail_flux_checkpoint_mpi_test PRIVATE
        MPI::MPI_CXX OpenMP::OpenMP_CXX)
    target_compile_definitions(bulk_tail_flux_checkpoint_mpi_test PRIVATE
        FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
        FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
        FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    add_test(NAME bulk_tail_flux_checkpoint_mpi_test
             COMMAND bulk_tail_flux_checkpoint_mpi_test)

    # Section 7.11.17.7 / 17F: the production cylindrical collision operator
    # must expose final face fluxes and preserve the trial transaction.  The
    # same source is compiled into focused cases so failures identify the
    # contract that broke rather than collapsing into one aggregate result.
    set(COLLISION_FLUX_TEST_SOURCES
        src/cylindrical_fp_collision.cpp
        src/collision_coefficients.cpp
        src/species.cpp
        src/maxwell.cpp)
    foreach(_case_name IN ITEMS
            cylindrical_collision_flux_test
            cylindrical_collision_flux_interface_test
            cylindrical_collision_flux_cross_test
            cylindrical_collision_flux_transaction_test)
        add_executable(${_case_name}
            tests/cylindrical_collision_flux_test.cpp
            ${COLLISION_FLUX_TEST_SOURCES})
        target_include_directories(${_case_name} PRIVATE src)
        target_link_libraries(${_case_name} PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX)
        target_compile_definitions(${_case_name} PRIVATE
            FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
            FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
            FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    endforeach()
    target_compile_definitions(cylindrical_collision_flux_interface_test PRIVATE
        VPFP_COLLISION_FLUX_CASE=2)
    target_compile_definitions(cylindrical_collision_flux_cross_test PRIVATE
        VPFP_COLLISION_FLUX_CASE=3)
    target_compile_definitions(cylindrical_collision_flux_transaction_test PRIVATE
        VPFP_COLLISION_FLUX_CASE=4)
    add_test(NAME cylindrical_collision_flux_test
             COMMAND cylindrical_collision_flux_test)
    add_test(NAME cylindrical_collision_flux_interface_test
             COMMAND cylindrical_collision_flux_interface_test)
    add_test(NAME cylindrical_collision_flux_cross_test
             COMMAND cylindrical_collision_flux_cross_test)
    add_test(NAME cylindrical_collision_flux_transaction_test
             COMMAND cylindrical_collision_flux_transaction_test)

    # Documented section-17F target names.  These wrappers include the same
    # production collision-face test implementation, so the target and test
    # names match the execution document without duplicating solver formulas.
    set(DOCUMENTED_COLLISION_FLUX_TESTS
        collision_flux_interface_test
        collision_flux_interface_cross_test
        collision_flux_conversion_transaction_test
        collision_flux_interface_mpi_test)
    foreach(_case_name IN ITEMS ${DOCUMENTED_COLLISION_FLUX_TESTS})
        add_executable(${_case_name}
            tests/${_case_name}.cpp
            ${COLLISION_FLUX_TEST_SOURCES})
        target_include_directories(${_case_name} PRIVATE src)
        target_link_libraries(${_case_name} PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX)
        target_compile_definitions(${_case_name} PRIVATE
            FP_VELOCITY_GRID_UPAR_TAIL_CELLS_PER_SIDE=${FP_UPAR_TAIL_CELLS_PER_SIDE}
            FP_VELOCITY_GRID_UPAR_EXTENDED_MAX=${FP_UPAR_EXTENDED_MAX}
            FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    endforeach()
    add_test(NAME collision_flux_interface_test
             COMMAND collision_flux_interface_test)
    add_test(NAME collision_flux_interface_cross_test
             COMMAND collision_flux_interface_cross_test)
    add_test(NAME collision_flux_conversion_transaction_test
             COMMAND collision_flux_conversion_transaction_test)
    add_test(NAME collision_flux_interface_mpi_test
             COMMAND collision_flux_interface_mpi_test)
endif()
