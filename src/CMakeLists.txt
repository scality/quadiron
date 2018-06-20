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
include(GNUInstallDirs)

# Source files.
set(LIB_SRC
  ${SOURCE_DIR}/core.cpp
  ${SOURCE_DIR}/fec_vectorisation.cpp
  ${SOURCE_DIR}/misc.cpp
  ${SOURCE_DIR}/gf_nf4.cpp
  ${SOURCE_DIR}/gf_ring.cpp
  ${SOURCE_DIR}/property.cpp
  ${SOURCE_DIR}/vec_vector.cpp

  CACHE
  INTERNAL
  ""
  FORCE
)

# Generate build_info.h (with compile-time information).
configure_file(${SOURCE_DIR}/build_info.in ${GENERATE_DIR}/build_info.h @ONLY)

###########
# Libraries
###########

set(INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}")

set(OBJECT_SYS_INCLUDES
  $<BUILD_INTERFACE:${GENERATE_DIR}>
  $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)

set(OBJECT_INCLUDES
  $<BUILD_INTERFACE:${SOURCE_DIR}>
  $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)


# Build an Object Library (can be reused for both static and dynamic libs).
add_library(${OBJECT_LIB} OBJECT ${LIB_SRC})
add_coverage(${OBJECT_LIB})
set_property(TARGET ${OBJECT_LIB} PROPERTY POSITION_INDEPENDENT_CODE 1)
target_include_directories(${OBJECT_LIB}        PUBLIC ${OBJECT_INCLUDES})
target_include_directories(${OBJECT_LIB} SYSTEM PUBLIC ${OBJECT_SYS_INCLUDES})

# Dynamic library.
add_library(${SHARED_LIB} SHARED $<TARGET_OBJECTS:${OBJECT_LIB}>)
# Static library.
add_library(${STATIC_LIB} STATIC $<TARGET_OBJECTS:${OBJECT_LIB}>)

# Set properties/add dependencies.
foreach(lib ${SHARED_LIB} ${STATIC_LIB})
  set_target_properties(${lib} PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME})
  target_include_directories(${lib}        PUBLIC ${OBJECT_INCLUDES})
  target_include_directories(${lib} SYSTEM PUBLIC ${OBJECT_SYS_INCLUDES})
endforeach()

##############
# Installation
##############

install(TARGETS ${SHARED_LIB} ${STATIC_LIB}
  EXPORT ${UPCASE_PROJECT_NAME}
  COMPONENT libraries
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${INSTALL_INCLUDE_DIR}
)
install(DIRECTORY .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}
  COMPONENT headers
  FILES_MATCHING PATTERN "*.h"
)
install(FILES ${GENERATE_DIR}/build_info.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}
  COMPONENT headers
)


include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_NAME          "${CMAKE_PROJECT_NAME}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
  "${UPCASE_PROJECT_NAME} — NTT Based Erasure Coding Library"
)
set(CPACK_RESOURCE_FILE_LICENSE    "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README     "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_VERSION_MAJOR    "${VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR    "${VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH    "${VERSION_PATCH}")

# Source distributions.
set(CPACK_SOURCE_GENERATOR         TBZ2)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${VERSION}")
set(CPACK_SOURCE_IGNORE_FILES
  ".swp$"
  "/.git/"
  "/build/"
  "/.gitignore"
  "${CPACK_SOURCE_IGNORE_FILES}"
)

include(CPack)

cpack_add_component_group(runtime
  DISPLAY_NAME Runtime
  DESCRIPTION  "Components required to run softwares using ${UPCASE_PROJECT_NAME}"
  EXPANDED
)
cpack_add_component_group(development
  DISPLAY_NAME "Development"
  DESCRIPTION  "Components required to develop softwares using ${UPCASE_PROJECT_NAME}"
  EXPANDED
)

cpack_add_component(libraries
  DISPLAY_NAME "${UPCASE_PROJECT_NAME} libraries"
  DESCRIPTION  "Libraries required to run programs using ${UPCASE_PROJECT_NAME}"
  REQUIRED
  GROUP        runtime
)
cpack_add_component(headers
  DISPLAY_NAME "${UPCASE_PROJECT_NAME} development headers"
  DESCRIPTION  "C++ header files required to build programs using ${UPCASE_PROJECT_NAME}"
  DISABLED
  GROUP        development
  DEPENDS      libraries
)

###############################
# Export QuadIron configuration
###############################
include(CMakePackageConfigHelpers)

install(EXPORT ${UPCASE_PROJECT_NAME}
  NAMESPACE ${UPCASE_PROJECT_NAME}::
  FILE ${UPCASE_PROJECT_NAME}Config.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${UPCASE_PROJECT_NAME}
)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${UPCASE_PROJECT_NAME}ConfigVersion.cmake
  VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}"
  COMPATIBILITY SameMajorVersion
)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${UPCASE_PROJECT_NAME}ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${UPCASE_PROJECT_NAME}
  COMPONENT headers
)
