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
cmake_minimum_required(VERSION 3.5)

######################
# Project informations
######################

project(quadiron CXX)

enable_testing()

set(CMAKE_CXX_STANDARD          14)
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
string(TOUPPER "${CMAKE_PROJECT_NAME}" UPCASE_PROJECT_NAME)

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
set(OBJECT_LIB  object)
set(SHARED_LIB  shared)
set(STATIC_LIB  static)
set(UNIT_TESTS  unit_tests)
set(EC_DRIVER   ec_driver)

file(MAKE_DIRECTORY ${GENERATE_DIR})

##################
# Setting for SIMD
##################
set(USE_SIMD "OFF" CACHE STRING "SIMD vectorization")
set_property(CACHE USE_SIMD PROPERTY STRINGS OFF ON SSE AVX)

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
  -Wall -Wextra -Wdeprecated -Wpedantic
  -Wdocumentation -Wdocumentation-pedantic
  -Wno-documentation-unknown-command  # Doesn't recognize @cite :(
  -Wunused-exception-parameter -Wunused-macros -Wunused-member-function -Wunused-template
  -Wgcc-compat -Wgnu
  -Warray-bounds-pointer-arithmetic
  -Wassign-enum
  -Wbad-function-cast -Wcast-qual -Wcast-align -Wold-style-cast -Wuseless-cast
  -Wconditional-uninitialized
  -Wduplicate-method-arg -Wduplicate-method-match
  -Wfloat-equal -Wdouble-promotion
  -Wformat-pedantic
  -Winit-self
  -Wimplicit-fallthrough=4
  -Wlogical-op
  -Wmissing-noreturn -Wmissing-prototypes -Wmissing-include-dirs -Wmissing-declarations
  -Wnested-anon-types
  -Wnon-virtual-dtor
  -Woverlength-strings -Wwrite-strings
  -Woverloaded-virtual
  -Wredundant-decls
  -Wshift-sign-overflow
  -Wstring-compare -Wstring-conversion -Wstring-plus-char
  -Wstrict-null-sentinel
  -Wsometimes-uninitialized
  -Wswitch-enum
  -Wunreachable-code-aggressive
)

# Option for enabling/disabling SIMD flags is for both of debug and release
if (USE_SIMD STREQUAL "ON")
  list(APPEND COMMON_CXX_FLAGS "-march=native")
  add_definitions(-DQUADIRON_USE_SIMD)
elseif (USE_SIMD STREQUAL "SSE")
  list(APPEND COMMON_CXX_FLAGS "-msse4.1")
  add_definitions(-DQUADIRON_USE_SIMD)
elseif (USE_SIMD STREQUAL "AVX")
  list(APPEND COMMON_CXX_FLAGS "-mavx2")
  add_definitions(-DQUADIRON_USE_SIMD)
endif()

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
)
foreach(flag ${DEBUG_CXX_FLAGS})
  check_cxx_compiler_flag(${flag} has_flag_${flag})
  if(has_flag_${flag})
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${flag}")
  endif(has_flag_${flag})
endforeach(flag)

# Seems buggy on Clang: see https://bugs.llvm.org/show_bug.cgi?id=16404
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ftrapv")
endif()

# Release flags
set(RELEASE_CXX_FLAGS
  -O3 -DNDEBUG
)
foreach(flag ${RELEASE_CXX_FLAGS})
  check_cxx_compiler_flag(${flag} has_flag_${flag})
  if(has_flag_${flag})
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${flag}")
  endif(has_flag_${flag})
endforeach(flag)

################
# QuadIron build
################

find_package(codecov)

# QuadIron
add_subdirectory(src)
# Documentation
add_subdirectory(documentation)
# Tests
add_subdirectory(test)
# Benchmark
add_subdirectory(benchmark)

coverage_evaluate()

# Be sure to recompile the test driver if needed before running the test suite.
add_custom_target(check
  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
  DEPENDS ${UNIT_TESTS}
  DEPENDS ${EC_DRIVER}
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
  bugprone-*
  cert-*
  -cert-dcl58-cpp
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

##################
# Uninstall target
##################

# From https://cmake.org/Wiki/CMake_FAQ#Can_I_do_.22make_uninstall.22_with_CMake.3F
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE
    @ONLY
  )

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
  )
endif()
