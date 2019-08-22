# Copyright 2017-2018 Scality
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# Source files.
set(BENCH_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/benchmark.cpp

  CACHE
  INTERNAL
  ""
  FORCE
)


if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
  message(WARNING "Benchmark results with an un-optimised build may be misleading")
endif()

############
# Executable
############

# Dependencies.
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(BENCH_DRIVER ${PROJECT_NAME}_bench)

add_executable(${BENCH_DRIVER}
  ${BENCH_SRC}
)
add_coverage(${BENCH_DRIVER})

# Link librairies.
target_link_libraries(${BENCH_DRIVER}
  ${STATIC_LIB}
  Threads::Threads
)

if (NOT APPLE)
    # Workaround a bug on some version of Ubuntu
    # See https://bugs.launchpad.net/ubuntu/+source/gcc-defaults/+bug/1228201
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
endif()

###########
# Execution
###########

add_custom_target(benchmark
    COMMAND ${SCRIPT_DIR}/benchmark.sh "$<TARGET_FILE:${BENCH_DRIVER}>"
  COMMENT "run the benchmark"
)
add_dependencies(benchmark ${BENCH_DRIVER})

#########################################################
# Google Benchmark
#########################################################

# Download and install GoogleBenchmark

include(ExternalProject)

# NOTE: we do not enable link-time optimisation as it require GCC_AR and
# GCC_RANLIB for gcc or LLVMAR_EXECUTABLE, LLVMNM_EXECUTABLE and
# LLVMRANLIB_EXECUTABLE for clang.
# If these packages are available, we can enable this option by adding the
# following flags: "-DBENCHMARK_ENABLE_LTO=true"
set(CMAKE_ARGS_LIST
    "-DBENCHMARK_ENABLE_GTEST_TESTS=OFF"
    "-DCMAKE_BUILD_TYPE=Release"
    "-DCMAKE_CXX_FLAGS=-fPIC"
    "-DCMAKE_AR=${CMAKE_AR}"
)

if(CMAKE_CXX_COMPILER MATCHES "^(\/arm64*|\/ARM64*)")
  list(APPEND CMAKE_ARGS_LIST "-DCMAKE_CROSSCOMPILING=1" "-DRUN_HAVE_POSIX_REGEX=0")
endif()


ExternalProject_Add(gbench
  PREFIX          ${CMAKE_CURRENT_BINARY_DIR}/gbench
  GIT_REPOSITORY  https://github.com/google/benchmark.git
  GIT_TAG         v1.5.0
  CMAKE_ARGS      ${CMAKE_ARGS_LIST}
  INSTALL_COMMAND ""
)

# Get paths of the installed GoogleBenchmark.
ExternalProject_Get_Property(gbench source_dir binary_dir)
set(GBENCH_SOURCE_DIR ${source_dir})

# GoogleTest target (to be used as a dependency by our test driver).
add_library(libgbench IMPORTED STATIC GLOBAL)
add_dependencies(libgbench gbench)
set_target_properties(libgbench PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/src/libbenchmark.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

##############################
# Components performance tests
##############################

set(COMP_PERF  comp_perf)

# Source files.
set(COMP_PERF_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/components/main.cpp
)

add_executable(${COMP_PERF}
  ${COMP_PERF_SRC}
)

target_link_libraries(${COMP_PERF}
  libgbench
  ${STATIC_LIB}
)
# Use SYSTEM so that our strict compilers settings are not applied on this code.
target_include_directories(${COMP_PERF} SYSTEM
    PUBLIC "${GBENCH_SOURCE_DIR}/include"
)

###########
# Execution
###########

add_custom_target(perf
  COMMAND "$<TARGET_FILE:${COMP_PERF}>" --benchmark_counters_tabular=true
  COMMENT "run the component performance"
)
add_dependencies(perf ${COMP_PERF})
