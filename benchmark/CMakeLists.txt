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

############
# ISA-L
############
include_directories(SYSTEM ./libs/usr/include)
link_directories(./libs/usr/lib)

link_libraries(isal)
link_libraries(boost_program_options)


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
