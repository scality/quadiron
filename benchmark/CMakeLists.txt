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
  ${GMPXX_LIBRARIES}
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
