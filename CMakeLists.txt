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
set(GENERATE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin/${CMAKE_BUILD_TYPE})

# Targets
set(OBJECT_LIB ${CMAKE_PROJECT_NAME}_object)
set(SHARED_LIB ${CMAKE_PROJECT_NAME}_shared)
set(STATIC_LIB ${CMAKE_PROJECT_NAME}_static)

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
# Tools
add_subdirectory(tools)
# Tests
add_subdirectory(test)
# Benchmark
add_subdirectory(benchmark)

coverage_evaluate()
