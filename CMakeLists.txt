cmake_minimum_required(VERSION 3.5)
project(FokkerPlanckSolver CXX)

# 为兼容 GCC 4.8.5 使用 C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -Wall")
option(FP_ENABLE_DEBUG_DIAGNOSTICS "Compile debug diagnostics support" OFF)
option(FP_BUILD_PPM_CTU_TESTS "Build standalone PPM reconstruction checks" OFF)
option(FP_BUILD_NONUNIFORM_RECONSTRUCTION_TESTS
       "Build standalone nonuniform MUSCL reconstruction checks" OFF)
option(FP_BUILD_FIXED_FIELD_LOW_ORDER_TEST
       "Build the section-7.1 standalone low-order audits" OFF)
option(FP_BUILD_PRODUCTION_AMPERE_TEST
       "Build the section-7.2 MPI production-solver Ampere audit" OFF)
option(FP_BUILD_LANGMUIR_WAVE_TEST
       "Build the section-7.3 MPI production Langmuir-wave tests" OFF)
option(FP_BUILD_CHECKPOINT_TESTS
       "Build checkpoint/restart infrastructure tests" OFF)
option(FP_BUILD_PERIODIC_STAGGERED_TESTS
       "Build periodic G/G* adjoint and x-face topology tests" OFF)
option(FP_BUILD_BACKGROUND_COUPLING_TESTS
       "Build section-11 production-kernel background coupling audits" OFF)
option(FP_BUILD_MIDPOINT_PREDICTOR_TEST
       "Build the standalone protected midpoint field predictor test" OFF)
set(FP_VELOCITY_GRID_UPAR_STRETCH "6.2" CACHE STRING
    "Smooth cylindrical u_parallel sinh mapping parameter")
set(FP_VELOCITY_GRID_UPERP_STRETCH "5.9" CACHE STRING
    "Smooth cylindrical u_perp sinh mapping parameter")
set(FP_AMPERE_TEST_SUITE_FINEST_NV "192" CACHE STRING
    "Finest Nv used to select the common physical 7.2 test timestep")
# Reconfiguration defines one immutable result namespace for the complete
# coarse/base/fine batch.  Individual executables only overwrite their own
# resolution file inside that namespace.
string(TIMESTAMP FP_TEST_RUN_ID "%Y%m%dT%H%M%SZ" UTC)
set(FP_TEST_RESULT_DIR "output/production_ampere_7_2_${FP_TEST_RUN_ID}")
set(FP_LANGMUIR_RESULT_DIR "output/production_langmuir_7_3_${FP_TEST_RUN_ID}")
add_definitions(
    -DFP_VELOCITY_GRID_UPAR_STRETCH=${FP_VELOCITY_GRID_UPAR_STRETCH}
    -DFP_VELOCITY_GRID_UPERP_STRETCH=${FP_VELOCITY_GRID_UPERP_STRETCH})

if(FP_BUILD_MIDPOINT_PREDICTOR_TEST)
    add_executable(midpoint_field_predictor_test
        tests/midpoint_field_predictor_test.cpp)
    target_include_directories(midpoint_field_predictor_test PRIVATE src)
endif()

# MPI: 通过 mpicxx 封装器自动处理（在 PBS 脚本中通过
# -DCMAKE_CXX_COMPILER=mpicxx 指定）。
# 如果用户直接用 g++ 编译，则回退查找 MPI。
find_package(MPI QUIET)
if(MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

if(FP_BUILD_CHECKPOINT_TESTS)
    if(NOT MPI_CXX_FOUND)
        message(FATAL_ERROR "FP_BUILD_CHECKPOINT_TESTS requires MPI")
    endif()
    add_executable(checkpoint_roundtrip_test
        tests/checkpoint_roundtrip_test.cpp src/checkpoint.cpp src/species.cpp
        src/maxwell.cpp src/periodic_staggered_operators.cpp src/beam_pic.cpp)
    target_include_directories(checkpoint_roundtrip_test PRIVATE src)
    target_link_libraries(checkpoint_roundtrip_test ${MPI_CXX_LIBRARIES})

    add_executable(restart_one_step_equivalence_test
        tests/restart_one_step_equivalence_test.cpp
        src/checkpoint.cpp src/species.cpp src/maxwell.cpp src/periodic_staggered_operators.cpp src/beam_pic.cpp
        src/vlasov_ampere_midpoint.cpp src/cylindrical_vlasov_midpoint.cpp
        src/nonuniform_reconstruction.cpp)
    target_include_directories(restart_one_step_equivalence_test PRIVATE src)
    target_link_libraries(restart_one_step_equivalence_test ${MPI_CXX_LIBRARIES})

    add_executable(beam_rng_injection_continuity_test
        tests/beam_rng_injection_continuity_test.cpp
        src/checkpoint.cpp src/species.cpp src/maxwell.cpp src/periodic_staggered_operators.cpp src/beam_pic.cpp)
    target_include_directories(beam_rng_injection_continuity_test PRIVATE src)
    target_link_libraries(beam_rng_injection_continuity_test ${MPI_CXX_LIBRARIES})

    add_executable(checkpoint_failure_safety_test
        tests/checkpoint_failure_safety_test.cpp
        src/checkpoint.cpp src/species.cpp src/maxwell.cpp src/periodic_staggered_operators.cpp src/beam_pic.cpp)
    target_include_directories(checkpoint_failure_safety_test PRIVATE src)
    target_link_libraries(checkpoint_failure_safety_test ${MPI_CXX_LIBRARIES})

    add_executable(beam_dt_half_time_consistency_test
        tests/beam_dt_half_time_consistency_test.cpp
        src/species.cpp src/maxwell.cpp src/periodic_staggered_operators.cpp src/beam_pic.cpp)
    target_include_directories(beam_dt_half_time_consistency_test PRIVATE src)
    target_link_libraries(beam_dt_half_time_consistency_test ${MPI_CXX_LIBRARIES})

    add_executable(fixed_midpoint_operator_audit_test
        tests/fixed_midpoint_operator_audit_test.cpp
        src/checkpoint.cpp src/species.cpp src/maxwell.cpp src/periodic_staggered_operators.cpp src/beam_pic.cpp
        src/vlasov_ampere_midpoint.cpp src/cylindrical_vlasov_midpoint.cpp
        src/nonuniform_reconstruction.cpp)
    target_include_directories(fixed_midpoint_operator_audit_test PRIVATE src)
    target_link_libraries(fixed_midpoint_operator_audit_test ${MPI_CXX_LIBRARIES})
endif()

if(FP_BUILD_LANGMUIR_WAVE_TEST)
    if(NOT MPI_CXX_FOUND)
        message(FATAL_ERROR "FP_BUILD_LANGMUIR_WAVE_TEST requires MPI")
    endif()
    add_executable(production_langmuir_wave_test
        tests/production_langmuir_wave_test.cpp
        src/species.cpp
        src/vlasov_ampere_midpoint.cpp
        src/cylindrical_vlasov_midpoint.cpp
        src/nonuniform_reconstruction.cpp
        src/maxwell.cpp
        src/periodic_staggered_operators.cpp
        src/beam_pic.cpp)
    target_include_directories(production_langmuir_wave_test PRIVATE src)
    target_compile_definitions(production_langmuir_wave_test PRIVATE
        FP_TEST_RUN_ID="${FP_TEST_RUN_ID}"
        FP_LANGMUIR_RESULT_DIR="${FP_LANGMUIR_RESULT_DIR}")
    if(FP_ENABLE_DEBUG_DIAGNOSTICS)
        target_compile_definitions(production_langmuir_wave_test
                                   PRIVATE FP_ENABLE_DEBUG_DIAGNOSTICS=1)
    else()
        target_compile_definitions(production_langmuir_wave_test
                                   PRIVATE FP_ENABLE_DEBUG_DIAGNOSTICS=0)
    endif()
    target_link_libraries(production_langmuir_wave_test ${MPI_CXX_LIBRARIES})
    if(MPI_CXX_COMPILE_FLAGS)
        set_target_properties(production_langmuir_wave_test PROPERTIES
                              COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    endif()
    if(MPI_CXX_LINK_FLAGS)
        set_target_properties(production_langmuir_wave_test PROPERTIES
                              LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
    endif()

    add_executable(production_langmuir_wave_audit
        tests/production_langmuir_wave_audit.cpp)
    target_include_directories(production_langmuir_wave_audit PRIVATE src)
endif()

if(FP_BUILD_PERIODIC_STAGGERED_TESTS)
    if(NOT MPI_CXX_FOUND)
        message(FATAL_ERROR "FP_BUILD_PERIODIC_STAGGERED_TESTS requires MPI")
    endif()
    add_executable(periodic_staggered_adjoint_test
        tests/periodic_staggered_adjoint_test.cpp
        src/periodic_staggered_operators.cpp)
    target_include_directories(periodic_staggered_adjoint_test PRIVATE src)
    target_link_libraries(periodic_staggered_adjoint_test ${MPI_CXX_LIBRARIES})

    add_executable(periodic_x_flux_topology_test
        tests/periodic_x_flux_topology_test.cpp
        src/species.cpp src/maxwell.cpp src/periodic_staggered_operators.cpp)
    target_include_directories(periodic_x_flux_topology_test PRIVATE src)
    target_link_libraries(periodic_x_flux_topology_test ${MPI_CXX_LIBRARIES})
endif()

if(FP_BUILD_BACKGROUND_COUPLING_TESTS)
    if(NOT MPI_CXX_FOUND)
        message(FATAL_ERROR "FP_BUILD_BACKGROUND_COUPLING_TESTS requires MPI")
    endif()
    function(add_background_coupling_kernel_test target_name source_file nx nv nmu)
        add_executable(${target_name}
            ${source_file}
            src/species.cpp
            src/vlasov_ampere_midpoint.cpp
            src/cylindrical_vlasov_midpoint.cpp
            src/nonuniform_reconstruction.cpp
            src/checkpoint.cpp
            src/maxwell.cpp
            src/periodic_staggered_operators.cpp
            src/beam_pic.cpp)
        target_include_directories(${target_name} PRIVATE src tests)
        target_compile_definitions(${target_name} PRIVATE
            FP_SPATIAL_GRID_NX=${nx}
            FP_VELOCITY_GRID_NV=${nv}
            FP_VELOCITY_GRID_NMU=${nmu})
        target_link_libraries(${target_name} ${MPI_CXX_LIBRARIES})
        if(MPI_CXX_COMPILE_FLAGS)
            set_target_properties(${target_name} PROPERTIES
                                  COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
        endif()
        if(MPI_CXX_LINK_FLAGS)
            set_target_properties(${target_name} PROPERTIES
                                  LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
        endif()
    endfunction()

    # 11.2 uses a compact spatial/velocity refinement sequence.  The same
    # source calls the production bundle; only immutable grid macros differ.
    add_background_coupling_kernel_test(background_coupling_smooth_periodic_test_coarse
        tests/background_coupling_smooth_periodic_test.cpp 1000 48 32)
    add_background_coupling_kernel_test(background_coupling_smooth_periodic_test
        tests/background_coupling_smooth_periodic_test.cpp 2000 96 64)
    add_background_coupling_kernel_test(background_coupling_smooth_periodic_test_fine
        tests/background_coupling_smooth_periodic_test.cpp 4000 192 128)
    add_executable(background_coupling_smooth_convergence_audit
        tests/background_coupling_smooth_convergence_audit.cpp)

    # 11.3--11.5 use the established production resolution and differ only
    # in what they audit from the same fixed production operator call.
    add_background_coupling_kernel_test(background_coupling_gradient_translation_test
        tests/background_coupling_gradient_translation_test.cpp 4000 96 64)
    add_background_coupling_kernel_test(background_coupling_velocity_consistency_test
        tests/background_coupling_velocity_consistency_test.cpp 4000 96 64)
    add_background_coupling_kernel_test(background_coupling_midpoint_state_test
        tests/background_coupling_midpoint_state_test.cpp 4000 96 64)

    # Direction-separated convergence audits for x, u_parallel and u_perp.
    # FCT is disabled inside the shared test source for every target so all
    # three resolutions use the same raw high-order operator branch.
    function(add_directional_grid_triplet prefix kind
             nx0 nu0 nmu0 nx1 nu1 nmu1 nx2 nu2 nmu2)
        add_background_coupling_kernel_test(${prefix}_coarse
            tests/background_coupling_directional_grid_test.cpp ${nx0} ${nu0} ${nmu0})
        add_background_coupling_kernel_test(${prefix}
            tests/background_coupling_directional_grid_test.cpp ${nx1} ${nu1} ${nmu1})
        add_background_coupling_kernel_test(${prefix}_fine
            tests/background_coupling_directional_grid_test.cpp ${nx2} ${nu2} ${nmu2})
        target_compile_definitions(${prefix}_coarse PRIVATE FP_DIRECTIONAL_GRID_KIND=${kind})
        target_compile_definitions(${prefix} PRIVATE FP_DIRECTIONAL_GRID_KIND=${kind})
        target_compile_definitions(${prefix}_fine PRIVATE FP_DIRECTIONAL_GRID_KIND=${kind})
    endfunction()
    add_directional_grid_triplet(background_coupling_x_refinement 1
        1000 96 64 2000 96 64 4000 96 64)
    add_directional_grid_triplet(background_coupling_u_parallel_refinement 2
        2000 48 64 2000 96 64 2000 192 64)
    add_directional_grid_triplet(background_coupling_u_perp_refinement 3
        2000 96 32 2000 96 64 2000 96 128)
    add_executable(background_coupling_directional_grid_audit
        tests/background_coupling_directional_grid_audit.cpp)

    # Test-only common-speed A/B.  The alternative speed is selected while
    # the production x high-order face flux is reconstructed, never by
    # rescaling an already constructed flux in the test.
    add_background_coupling_kernel_test(background_coupling_common_speed_flux_test
        tests/background_coupling_common_speed_flux_test.cpp 2000 96 64)

    # Section 12: frozen production u_parallel reconstruction/operator
    # audits.  These are small-grid diagnostics only; they never select a
    # different production transport or Ampere path.
    add_background_coupling_kernel_test(background_coupling_upar_dual_audit
        tests/background_coupling_upar_dual_audit.cpp 800 96 64)
    target_compile_definitions(background_coupling_upar_dual_audit PRIVATE
        FP_UPAR_AUDIT_DETAILED=1)
    # Uses the production checkpoint grid.  It is separate from the compact
    # manufactured-state scan because midpoint checkpoints are grid-exact.
    add_background_coupling_kernel_test(background_coupling_upar_checkpoint_audit
        tests/background_coupling_upar_dual_audit.cpp 4000 96 64)
    target_compile_definitions(background_coupling_upar_checkpoint_audit PRIVATE
        FP_UPAR_AUDIT_DETAILED=1)
    add_background_coupling_kernel_test(background_coupling_upar_transpose_audit_coarse
        tests/background_coupling_upar_dual_audit.cpp 800 48 64)
    add_background_coupling_kernel_test(background_coupling_upar_transpose_audit
        tests/background_coupling_upar_dual_audit.cpp 800 96 64)
    add_background_coupling_kernel_test(background_coupling_upar_transpose_audit_fine
        tests/background_coupling_upar_dual_audit.cpp 800 192 64)
    add_executable(background_coupling_upar_transpose_convergence_audit
        tests/background_coupling_upar_transpose_convergence_audit.cpp)
    # Section 15.2 stage-1: fixed midpoint, no Beam, no FCT.  It compares the
    # untouched legacy path with the opt-in local dual-u rank-one operator.
    add_background_coupling_kernel_test(background_coupling_dual_u_center_test
        tests/background_coupling_dual_u_center_test.cpp 800 96 64)

    # Section 15.3 stage-2: fixed production midpoint operator, no Beam and
    # no FCT.  Every producer compares the untouched legacy coupling with the
    # opt-in dual-u operator and writes the same conservation/pairing ledger.
    set(FP_DUAL_U_STAGE2_GIT_REVISION "unknown")
    set(FP_DUAL_U_STAGE2_GIT_DIRTY -1)
    execute_process(COMMAND git rev-parse --short=12 HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        RESULT_VARIABLE FP_DUAL_U_STAGE2_GIT_RESULT
        OUTPUT_VARIABLE FP_DUAL_U_STAGE2_GIT_REVISION_RAW
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(FP_DUAL_U_STAGE2_GIT_RESULT EQUAL 0)
        set(FP_DUAL_U_STAGE2_GIT_REVISION
            "${FP_DUAL_U_STAGE2_GIT_REVISION_RAW}")
        execute_process(COMMAND git diff --quiet --ignore-submodules HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            RESULT_VARIABLE FP_DUAL_U_STAGE2_DIRTY_RESULT
            ERROR_QUIET)
        if(FP_DUAL_U_STAGE2_DIRTY_RESULT EQUAL 0)
            set(FP_DUAL_U_STAGE2_GIT_DIRTY 0)
        else()
            set(FP_DUAL_U_STAGE2_GIT_DIRTY 1)
        endif()
    endif()
    function(add_dual_u_stage2_test target_name case_id nx nv nmu)
        add_background_coupling_kernel_test(${target_name}
            tests/background_coupling_dual_u_stage2_test.cpp ${nx} ${nv} ${nmu})
        target_compile_definitions(${target_name} PRIVATE
            FP_DUAL_U_STAGE2_CASE=${case_id}
            FP_DUAL_U_STAGE2_GIT_REVISION="${FP_DUAL_U_STAGE2_GIT_REVISION}"
            FP_DUAL_U_STAGE2_GIT_DIRTY=${FP_DUAL_U_STAGE2_GIT_DIRTY})
    endfunction()
    add_dual_u_stage2_test(background_coupling_dual_u_zero_field_test
        1 800 96 64)
    add_dual_u_stage2_test(background_coupling_dual_u_field_symmetry_test
        2 800 96 64)
    add_dual_u_stage2_test(background_coupling_dual_u_phase_test
        3 800 96 64)
    add_dual_u_stage2_test(background_coupling_dual_u_nonuniform_grid_test
        4 800 96 64)
    add_dual_u_stage2_test(background_coupling_dual_u_manufactured_coarse
        5 800 48 64)
    add_dual_u_stage2_test(background_coupling_dual_u_manufactured_base
        5 800 96 64)
    add_dual_u_stage2_test(background_coupling_dual_u_manufactured_fine
        5 800 192 64)
    add_dual_u_stage2_test(background_coupling_dual_u_endpoint_flux_test
        6 800 96 64)
    add_dual_u_stage2_test(background_coupling_dual_u_mpi_consistency_test
        7 800 96 64)
    add_executable(background_coupling_dual_u_convergence_audit
        tests/background_coupling_dual_u_convergence_audit.cpp)
    add_executable(background_coupling_dual_u_mpi_consistency_audit
        tests/background_coupling_dual_u_mpi_consistency_audit.cpp)

    # Section 15.4 stage-3: reconnect the opt-in dual-u high candidate to the
    # unchanged donor-cell base and production shared-face FCT.  One fixed
    # operator test audits all six required low/high/final layers.
    add_background_coupling_kernel_test(background_coupling_dual_u_muscl_fct_test
        tests/background_coupling_dual_u_stage3_test.cpp 800 96 64)
    # Section 13.3: unit-level conservation, boundary and positivity-budget
    # checks for the final limited dual-u pairing correction.
    add_executable(final_dual_u_limited_pairing_test
        tests/final_dual_u_limited_pairing_test.cpp)
    target_include_directories(final_dual_u_limited_pairing_test PRIVATE src)
    target_compile_definitions(final_dual_u_limited_pairing_test PRIVATE
        FP_SPATIAL_GRID_NX=8
        FP_VELOCITY_GRID_NV=16
        FP_VELOCITY_GRID_NMU=8)
    add_executable(face_pairing_regularized_unit_test
        tests/face_pairing_regularized_unit_test.cpp)
    target_include_directories(face_pairing_regularized_unit_test PRIVATE src)
    add_background_coupling_kernel_test(face_pairing_checkpoint_audit
        tests/face_pairing_checkpoint_audit.cpp 4000 96 64)
    # Section 15.5 stage-4: compare accepted legacy/dual production midpoint
    # dumps.  It reads dumps only and never changes the production solver.
    add_background_coupling_kernel_test(background_coupling_stage4_ab_compare
        tests/background_coupling_stage4_ab_compare.cpp 4000 96 64)
    # Run the three producers before the aggregate reader.  This target is
    # intentionally serial so the convergence audit can never observe a
    # partially written fine-grid result.
    add_custom_target(background_coupling_upar_transpose_audit_chain
        COMMAND ${CMAKE_COMMAND} -E rm -f
            output/background_coupling_upar_transpose_Nx800_Nu48_Nuperp64.result
            output/background_coupling_upar_transpose_Nx800_Nu96_Nuperp64.result
            output/background_coupling_upar_transpose_Nx800_Nu192_Nuperp64.result
            output/background_coupling_upar_transpose_convergence_audit.result
        COMMAND $<TARGET_FILE:background_coupling_upar_transpose_audit_coarse>
        COMMAND $<TARGET_FILE:background_coupling_upar_transpose_audit>
        COMMAND $<TARGET_FILE:background_coupling_upar_transpose_audit_fine>
        COMMAND $<TARGET_FILE:background_coupling_upar_transpose_convergence_audit>
            output/background_coupling_upar_transpose_Nx800_Nu48_Nuperp64.result
            output/background_coupling_upar_transpose_Nx800_Nu96_Nuperp64.result
            output/background_coupling_upar_transpose_Nx800_Nu192_Nuperp64.result
        DEPENDS
            background_coupling_upar_transpose_audit_coarse
            background_coupling_upar_transpose_audit
            background_coupling_upar_transpose_audit_fine
            background_coupling_upar_transpose_convergence_audit
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        USES_TERMINAL)
endif()

if(FP_BUILD_PPM_CTU_TESTS)
    add_executable(ppm_ctu_reconstruction_test tests/ppm_ctu_reconstruction_test.cpp)
    target_include_directories(ppm_ctu_reconstruction_test PRIVATE src)
endif()

if(FP_BUILD_NONUNIFORM_RECONSTRUCTION_TESTS)
    add_executable(nonuniform_reconstruction_test
        tests/nonuniform_reconstruction_test.cpp
        src/nonuniform_reconstruction.cpp)
    target_include_directories(nonuniform_reconstruction_test PRIVATE src)
endif()

if(FP_BUILD_FIXED_FIELD_LOW_ORDER_TEST)
    add_executable(fixed_field_low_order_test
        tests/fixed_field_low_order_test.cpp)
    target_include_directories(fixed_field_low_order_test PRIVATE src)

    add_executable(fixed_field_low_order_joint_test
        tests/fixed_field_low_order_joint_test.cpp)
    target_include_directories(fixed_field_low_order_joint_test PRIVATE src)

    add_executable(cylindrical_velocity_grid_test
        tests/cylindrical_velocity_grid_test.cpp)
    target_include_directories(cylindrical_velocity_grid_test PRIVATE src)

    # Fixed Maxwellian, signed small-field audit of the direct C_u -> J_E
    # construction.  This is standalone and never alters production modes.
    add_executable(fixed_maxwellian_energy_current_scaling_test
        tests/fixed_maxwellian_energy_current_scaling_test.cpp
        src/nonuniform_reconstruction.cpp)
    target_include_directories(fixed_maxwellian_energy_current_scaling_test
                               PRIVATE src)

endif()

if(FP_BUILD_PRODUCTION_AMPERE_TEST)
    if(NOT MPI_CXX_FOUND)
        message(FATAL_ERROR "FP_BUILD_PRODUCTION_AMPERE_TEST requires MPI")
    endif()
    function(add_production_ampere_velocity_test target_name nv nmu)
        add_executable(${target_name}
            tests/production_self_consistent_ampere_test.cpp
            src/species.cpp
            src/vlasov_ampere_midpoint.cpp
            src/cylindrical_vlasov_midpoint.cpp
            src/nonuniform_reconstruction.cpp
            src/maxwell.cpp
            src/periodic_staggered_operators.cpp
            src/beam_pic.cpp)
        target_include_directories(${target_name} PRIVATE src)
        target_compile_definitions(${target_name} PRIVATE
            FP_VELOCITY_GRID_NV=${nv}
            FP_VELOCITY_GRID_NMU=${nmu}
            FP_AMPERE_TEST_SUITE_FINEST_NV=${FP_AMPERE_TEST_SUITE_FINEST_NV}
            FP_TEST_RUN_ID=\"${FP_TEST_RUN_ID}\"
            FP_TEST_RESULT_DIR=\"${FP_TEST_RESULT_DIR}\")
        if(FP_ENABLE_DEBUG_DIAGNOSTICS)
            target_compile_definitions(${target_name}
                                       PRIVATE FP_ENABLE_DEBUG_DIAGNOSTICS=1)
        else()
            target_compile_definitions(${target_name}
                                       PRIVATE FP_ENABLE_DEBUG_DIAGNOSTICS=0)
        endif()
        target_link_libraries(${target_name} ${MPI_CXX_LIBRARIES})
        if(MPI_CXX_COMPILE_FLAGS)
            set_target_properties(${target_name} PROPERTIES
                                  COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
        endif()
        if(MPI_CXX_LINK_FLAGS)
            set_target_properties(${target_name} PROPERTIES
                                  LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
        endif()
    endfunction()

    # The unqualified target retains the production 96 x 64 grid.  The other
    # two targets provide a factor-two velocity-grid convergence sequence.
    add_production_ampere_velocity_test(production_self_consistent_ampere_test
                                        96 64)
    add_production_ampere_velocity_test(production_self_consistent_ampere_test_coarse
                                        48 32)
    add_production_ampere_velocity_test(production_self_consistent_ampere_test_fine
                                        192 128)

    # Manual controlled-FCT regression: a strictly positive but steep
    # x/u_parallel state must trigger the production limiter.  It is separate
    # from the section-7.2 convergence targets and never changes fp_solver.
    add_production_ampere_velocity_test(controlled_fct_activation_test 96 64)
    target_compile_definitions(controlled_fct_activation_test PRIVATE
                               FP_FCT_ACTIVATION_TEST=1)

    # One-substep signed-field scaling audit of the production Vlasov-Ampere
    # path.  It is an explicit test target and never changes fp_solver.
    add_executable(production_energy_current_scaling_test
        tests/production_energy_current_scaling_test.cpp
        src/species.cpp
        src/vlasov_ampere_midpoint.cpp
        src/cylindrical_vlasov_midpoint.cpp
        src/nonuniform_reconstruction.cpp
        src/maxwell.cpp
        src/periodic_staggered_operators.cpp
        src/beam_pic.cpp)
    target_include_directories(production_energy_current_scaling_test PRIVATE src)
    target_compile_definitions(production_energy_current_scaling_test PRIVATE
        FP_TEST_RUN_ID=\"${FP_TEST_RUN_ID}\")
    target_link_libraries(production_energy_current_scaling_test ${MPI_CXX_LIBRARIES})
    if(MPI_CXX_COMPILE_FLAGS)
        set_target_properties(production_energy_current_scaling_test PROPERTIES
                              COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    endif()
    if(MPI_CXX_LINK_FLAGS)
        set_target_properties(production_energy_current_scaling_test PROPERTIES
                              LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
    endif()

    add_executable(production_self_consistent_ampere_convergence_audit
        tests/production_self_consistent_ampere_convergence_audit.cpp)

    # This target is intentionally manual: it is not attached to ALL and does
    # not submit a scheduler job.  Invoke it inside an already allocated MPI
    # job to run the complete 6.10 evidence chain with isolated result files.
    if(FP_BUILD_FIXED_FIELD_LOW_ORDER_TEST)
        add_custom_target(low_order_audit_suite
            COMMAND ${CMAKE_COMMAND} -E make_directory ${FP_TEST_RESULT_DIR}
            COMMAND $<TARGET_FILE:cylindrical_velocity_grid_test>
            COMMAND $<TARGET_FILE:fixed_field_low_order_test>
            COMMAND $<TARGET_FILE:fixed_field_low_order_joint_test>
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:production_self_consistent_ampere_test>
                ${MPIEXEC_POSTFLAGS}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:production_self_consistent_ampere_test_coarse>
                ${MPIEXEC_POSTFLAGS}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:production_self_consistent_ampere_test_fine>
                ${MPIEXEC_POSTFLAGS}
            COMMAND $<TARGET_FILE:production_self_consistent_ampere_convergence_audit>
                ${FP_TEST_RESULT_DIR}/production_ampere_7_2_Nv48_Nmu32.result
                ${FP_TEST_RESULT_DIR}/production_ampere_7_2_Nv96_Nmu64.result
                ${FP_TEST_RESULT_DIR}/production_ampere_7_2_Nv192_Nmu128.result
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS cylindrical_velocity_grid_test fixed_field_low_order_test
                fixed_field_low_order_joint_test
                production_self_consistent_ampere_test
                production_self_consistent_ampere_test_coarse
                production_self_consistent_ampere_test_fine
                production_self_consistent_ampere_convergence_audit)
    endif()
endif()

# 查找 OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# 源文件列表
set(SOURCES
    src/main.cpp
    src/config.cpp
    src/runtime_options.cpp
    src/checkpoint.cpp
    src/species.cpp
    src/vlasov_ampere_midpoint.cpp
    src/cylindrical_vlasov_midpoint.cpp
    src/nonuniform_reconstruction.cpp
    src/maxwell.cpp
    src/periodic_staggered_operators.cpp
    src/collision.cpp
    src/diagnostics.cpp
    src/beam_pic.cpp
)

add_executable(fp_solver ${SOURCES})
target_include_directories(fp_solver PRIVATE src)
if(FP_ENABLE_DEBUG_DIAGNOSTICS)
    target_compile_definitions(fp_solver PRIVATE FP_ENABLE_DEBUG_DIAGNOSTICS=1)
else()
    target_compile_definitions(fp_solver PRIVATE FP_ENABLE_DEBUG_DIAGNOSTICS=0)
endif()

# 链接库
if(MPI_CXX_FOUND)
    target_link_libraries(fp_solver ${MPI_CXX_LIBRARIES})
    if(MPI_CXX_COMPILE_FLAGS)
        set_target_properties(fp_solver PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    endif()
    if(MPI_CXX_LINK_FLAGS)
        set_target_properties(fp_solver PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
    endif()
endif()
