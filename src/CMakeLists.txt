# Source files.
set(LIB_SRC
  ${SOURCE_DIR}/misc.cpp
  ${SOURCE_DIR}/property.cpp
  ${SOURCE_DIR}/vec.cpp

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

# Build an Object Library (can be reused for both static and dynamic libs).
add_library(${OBJECT_LIB} OBJECT ${LIB_SRC})
add_coverage(${OBJECT_LIB})
set_property(TARGET ${OBJECT_LIB} PROPERTY POSITION_INDEPENDENT_CODE 1)

# Dynamic library.
add_library(${SHARED_LIB} SHARED $<TARGET_OBJECTS:${OBJECT_LIB}>)
set_target_properties(${SHARED_LIB} PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME})
# Static library.
add_library(${STATIC_LIB} STATIC $<TARGET_OBJECTS:${OBJECT_LIB}>)
set_target_properties(${STATIC_LIB} PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME})
