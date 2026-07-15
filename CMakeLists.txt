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

# MPI: 通过 mpicxx 封装器自动处理（在 PBS 脚本中通过
# -DCMAKE_CXX_COMPILER=mpicxx 指定）。
# 如果用户直接用 g++ 编译，则回退查找 MPI。
find_package(MPI QUIET)
if(MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
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
    src/species.cpp
    src/vlasov_ampere_midpoint.cpp
    src/cylindrical_vlasov_midpoint.cpp
    src/nonuniform_reconstruction.cpp
    src/maxwell.cpp
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
