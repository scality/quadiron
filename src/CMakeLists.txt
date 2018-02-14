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

# Dependencies.
find_package(GMPXX REQUIRED)

set(OBJECT_SYS_INCLUDES
  ${GENERATE_DIR}
  ${GMPXX_INCLUDE_DIRS}
)

set(OBJECT_INCLUDES
  ${SOURCE_DIR}
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

  target_link_libraries(${lib}
    ${GMPXX_LIBRARIES}
  )
endforeach()
