cmake_minimum_required(VERSION 2.8.12)
project(FokkerPlanckSolver CXX)

# 为兼容 GCC 4.8.5 使用 C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -Wall")

# MPI: 通过 mpicxx 封装器自动处理（在 PBS 脚本中通过
# -DCMAKE_CXX_COMPILER=mpicxx 指定）。
# 如果用户直接用 g++ 编译，则回退查找 MPI。
find_package(MPI QUIET)
if(MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

# 查找 OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# 源文件列表
set(SOURCES
    src/main.cpp
    src/species.cpp
    src/vlasov.cpp
    src/maxwell.cpp
    src/collision.cpp
    src/diagnostics.cpp
)

add_executable(fp_solver ${SOURCES})
target_include_directories(fp_solver PRIVATE src)

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
