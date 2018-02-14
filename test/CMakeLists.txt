#############
# Build tests
#############

# Source files.
set(TEST_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/arith_utest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fec_utest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fft_utest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gf_utest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mat_utest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/poly_utest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/rs_utest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/vecp_utest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/vec_utest.cpp

  CACHE
  INTERNAL
  ""
  FORCE
)

add_executable(${TEST_DRIVER}
  ${TEST_SRC}
)
add_coverage(${TEST_DRIVER})

target_link_libraries(${TEST_DRIVER}
  ${STATIC_LIB}
)

# Don't disable assert when compiling tests…
add_definitions(-UNDEBUG)

add_test(
  NAME unit_test
  COMMAND ${TEST_DRIVER}
)
