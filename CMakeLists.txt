cmake_minimum_required(VERSION 3.5)

######################
# Project informations
######################

project(ntl CXX)

enable_testing()

set(CMAKE_CXX_STANDARD          11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)
set(VERSION_PATCH 0)
execute_process(
  COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE VERSION_REV
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(VERSION_REV)
  set(VERSION_REV "-${VERSION_REV}")
endif()

set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}${VERSION_REV})
set(PACKAGE ${CMAKE_PROJECT_NAME})

################
# Useful defines
################

# Paths
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SCRIPT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
set(GENERATE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin/${CMAKE_BUILD_TYPE})

# Targets
set(OBJECT_LIB  ${CMAKE_PROJECT_NAME}_object)
set(SHARED_LIB  ${CMAKE_PROJECT_NAME}_shared)
set(STATIC_LIB  ${CMAKE_PROJECT_NAME}_static)
set(TEST_DRIVER ${CMAKE_PROJECT_NAME}_test)
set(NTL_TOOLS   ${CMAKE_PROJECT_NAME}_tools)

include_directories(${SOURCE_DIR})

file(MAKE_DIRECTORY ${GENERATE_DIR})
# Use SYSTEM so that our strict compilers settings are not applied to generated
# code.
include_directories(SYSTEM
    ${GENERATE_DIR}
)

####################
# Default build type
####################

# From https://blog.kitware.com/cmake-and-the-default-build-type/

# Set a default build type if none was specified.
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}"
    CACHE
    STRING
    "Choose the type of build."
    FORCE
  )
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo"
  )
endif()

################
# Compiler flags
################

include(CheckCXXCompilerFlag)

set(COMMON_CXX_FLAGS
  -pipe
  -Wall
)

# Manually add -Werror, for some reasons I can't make it works in the foreach…
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  add_compile_options(-Werror)
endif()

foreach(flag ${COMMON_CXX_FLAGS})
  check_cxx_compiler_flag(${flag} has_flag_${flag})
  if(has_flag_${flag})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
  endif(has_flag_${flag})
endforeach(flag)

# Debug flags
set(DEBUG_CXX_FLAGS
  -O0 -g3 -ggdb3 -fno-limit-debug-info
  -ftrapv
)
foreach(flag ${DEBUG_CXX_FLAGS})
  check_cxx_compiler_flag(${flag} has_flag_${flag})
  if(has_flag_${flag})
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${flag}")
  endif(has_flag_${flag})
endforeach(flag)

# Release flags
foreach(flag -O3 -flto -DNDEBUG)
  check_cxx_compiler_flag(${flag} has_flag_${flag})
  if(has_flag_${flag})
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${flag}")
  endif(has_flag_${flag})
endforeach(flag)

###########
# NTL build
###########

find_package(codecov)

# Dependencies.
find_package(GMPXX REQUIRED)
include_directories(SYSTEM ${GMPXX_INCLUDE_DIRS})

# NTL
add_subdirectory(src)
# Documentation
add_subdirectory(documentation)
# Tools
add_subdirectory(tools)
# Tests
add_subdirectory(test)
# Benchmark
add_subdirectory(benchmark)

coverage_evaluate()

# Be sure to recompile the test driver if needed before running the test suite.
add_custom_target(check
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS ${TEST_DRIVER}
  DEPENDS ${NTL_TOOLS}
)

#####################
# Formatting & Linter
#####################

file(GLOB_RECURSE ALL_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/benchmark/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/test/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/tools/*.h
)

set(ALL_SRCS
  ${BENCH_SRC}
  ${LIB_SRC}
  ${TEST_SRC}
  ${TOOLS_SRC}
  ${ALL_HEADERS}
)

add_custom_target(format
  COMMAND clang-format -i ${ALL_SRCS}
  COMMENT "fix code formatting"
)
add_dependencies(format ${OBJECT_LIB})

set(CLANG_LINTS_LIST
  android-*
  boost-*
  clang-analyzer-*
  mpi-*
)
string(REPLACE ";" "," CLANG_LINTS "${CLANG_LINTS_LIST}")

add_custom_target(lint
  COMMAND ${SCRIPT_DIR}/run-clang-tidy.py -checks=${CLANG_LINTS} -header-filter=.*
  COMMENT "run the linter"
)
add_dependencies(lint ${OBJECT_LIB})

add_custom_target(fix-lint
  COMMAND ${SCRIPT_DIR}/run-clang-tidy.py -fix -checks=${CLANG_LINTS} -header-filter=.*
  COMMENT "run the linter and apply proposed fixes"
)
add_dependencies(fix-lint ${OBJECT_LIB})

####################
# Targets for the CI
####################

# clang-tidy returns 0, even if there are lint errors, so we have to be hackish.
# Should be no longer needed when https://reviews.llvm.org/D39105 is merged…
set(LINT_ERROR_PATTERN "'(error|warning):'")
add_custom_target(check-lint
  COMMAND test `${SCRIPT_DIR}/run-clang-tidy.py -checks=${CLANG_LINTS} -header-filter=.* | grep -Pc ${LINT_ERROR_PATTERN}` -eq 0
  COMMENT "check for lint errors"
)
add_dependencies(check-lint ${OBJECT_LIB})

add_custom_target(check-format
  COMMAND test `clang-format -output-replacements-xml ${ALL_SRCS} ${TEST_SRC} | grep -c offset` -eq 0
  COMMENT "check code formatting"
)
add_dependencies(check-format ${OBJECT_LIB})
