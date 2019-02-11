# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/binning_directory/CMakeLists.txt
#
# CMakeLists.txt file for the binning_directory module tests.
# ===========================================================================

cmake_minimum_required (VERSION 3.2)
project (sims CXX)
message (STATUS "Configuring tests/binning_directory")

# ----------------------------------------------------------------------------
# Define paths
# ----------------------------------------------------------------------------

# CMAKE_SOURCE_DIR is seqan/tests/binning_directory/benchmark (the path of this CMakeLists.txt)
set(SEQAN_ROOT "${CMAKE_SOURCE_DIR}/../../..")
message(STATUS ${SEQAN_ROOT})

# SeqAn include
set (CMAKE_INCLUDE_PATH ${SEQAN_ROOT}/include ${CMAKE_INCLUDE_PATH})

# SeqAn cmake module
set (CMAKE_PREFIX_PATH ${SEQAN_ROOT}/util/cmake ${CMAKE_PREFIX_PATH})
set (CMAKE_MODULE_PATH ${SEQAN_ROOT}/util/cmake ${CMAKE_MODULE_PATH})

include_directories (${CMAKE_INCLUDE_PATH})

# ----------------------------------------------------------------------------
# Load SeqAn module
# ----------------------------------------------------------------------------

find_package (SeqAn REQUIRED)

# add include dir of test directory for some global test utility functions.
set (SEQAN_INCLUDE_DIRS "${SEQAN_INCLUDE_DIRS};${CMAKE_SOURCE_DIR}/include")

# import package's information (libraries are added with targets below)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})

# add warning levels
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Werror -Wall -Wextra")

# CMAKE_SOURCE_DIR is seqan/include/seqan/kmer/test (the path of this CMakeLists.txt)
set(SEQAN_ROOT "${CMAKE_SOURCE_DIR}/../../../..")
set(BENCHMARK_ROOT "${PROJECT_BINARY_DIR}/vendor/benchmark")

# ----------------------------------------------------------------------------
# Use CCache if possible
# ----------------------------------------------------------------------------

find_program(CCACHE_PROGRAM ccache)
# TODO(marehr): Once we require cmake>=3.4 use CMAKE_CXX_COMPILER_LAUNCHER instead
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif(CCACHE_PROGRAM)

# ----------------------------------------------------------------------------
# Add google/benchmark
# ----------------------------------------------------------------------------

include(ExternalProject)
# TODO(marehr): Once we require cmake>=3.4 add "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER};"
# to CMAKE_ARGS
set(
    google_benchmark_args
    "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}"
    "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    "-DBENCHMARK_ENABLE_TESTING=false"
    # "-DBENCHMARK_ENABLE_LTO=true"
)
ExternalProject_Add(
    google_benchmark
    PREFIX google_benchmark
    GIT_REPOSITORY "https://github.com/google/benchmark.git"
    SOURCE_DIR "${BENCHMARK_ROOT}"
    CMAKE_ARGS "${google_benchmark_args}"
    UPDATE_DISCONNECTED yes
)

add_library (benchmark STATIC IMPORTED)
add_dependencies(benchmark google_benchmark)
set_target_properties(benchmark PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libbenchmark.a")

macro(seqan3_benchmark target_cpp)
    # $target_cpp = "pod_tuple_benchmark.cpp"
    #   * will register the global TARGET name "pod_tuple_benchmark" and
    #   * will register the benchmark name "core_pod_tuple" if
    #     pod_tuple_benchmark.cpp is in benchmark/core/
    #
    # NOTE(marehr): ".+/benchmark/" REGEX is greedy, that means
    # /benchmark/benchmark/benchmark/hello_benchmark.cpp will result in an empty `benchmark_path`
    string(REGEX REPLACE "_benchmark.cpp$" "" target_name ${target_cpp})
    string(REGEX REPLACE ".+/performance/" "" benchmark_path ${CMAKE_CURRENT_LIST_DIR})
    string(REGEX REPLACE ".+/unit/" "" benchmark_path ${CMAKE_CURRENT_LIST_DIR})
    string(REGEX REPLACE "/" "_" benchmark_path ${benchmark_path})

    add_executable(${target_name}_benchmark ${target_cpp})
    target_link_libraries(${target_name}_benchmark "pthread" "benchmark")
    add_dependencies(${target_name}_benchmark benchmark)
endmacro()

# define all include folders
include_directories("${BENCHMARK_ROOT}/include/")
include_directories("${SEQAN_ROOT}/include/")
include_directories (${SDSL_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
# Add the benchmarks
# ----------------------------------------------------------------------------

add_executable(bench g_benchmark.cpp)
target_link_libraries(bench "pthread" "benchmark")
add_dependencies(bench benchmark)

add_executable(sim sim.cpp)

add_executable(count count.cpp
                     helper.h)
target_link_libraries(count "pthread")

#add_subdirectories()
