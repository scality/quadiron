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
  ${CMAKE_CURRENT_SOURCE_DIR}/buffers_utest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/vector_utest.cpp

  CACHE
  INTERNAL
  ""
  FORCE
)

add_executable(${UNIT_TESTS}
  ${TEST_SRC}
)
add_coverage(${UNIT_TESTS})

target_link_libraries(${UNIT_TESTS}
  ${STATIC_LIB}
)

# Don't disable assert when compiling tests…
add_definitions(-UNDEBUG)

add_test(
  NAME unit_test
  COMMAND ${UNIT_TESTS}
)
