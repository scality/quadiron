# Source files.
#TODO: add find_pminus1_compo.cpp, is_prime.cpp and find_primes.cpp
set(TOOLS_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/ec.cpp

  CACHE
  INTERNAL
  ""
  FORCE
)

################
# Erasure Coding
################

set(EC_BIN_NAME ec)

add_executable(${EC_BIN_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}/ec.cpp
)
add_coverage(${EC_BIN_NAME})

# Link librairies.
target_link_libraries(${EC_BIN_NAME}
  ${STATIC_LIB}
)

#TODO: add targets for find_pminus1_compo, is_prime and find_primes

#############
# Root target
#############

add_custom_target(${NTTEC_TOOLS}
  COMMENT "build the NTTEC tools"
)

#TODO: add find_pminus1_compo, is_prime and find_primes targets as dependencies.
add_dependencies(${NTTEC_TOOLS}
  ${EC_BIN_NAME}
)

#######
# Tests
#######

add_test(
  NAME ec_test
  COMMAND ${SCRIPT_DIR}/test_ec.sh "$<TARGET_FILE:${EC_BIN_NAME}>"
)
