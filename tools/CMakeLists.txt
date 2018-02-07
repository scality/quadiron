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

# Link librairies.
target_link_libraries(${EC_BIN_NAME}
  ${STATIC_LIB}
  ${GMPXX_LIBRARIES}
)

#TODO: add targets for find_pminus1_compo, is_prime and find_primes

#############
# Root target
#############

add_custom_target(tools
  COMMENT "build the NTL tools"
)

#TODO: add find_pminus1_compo, is_prime and find_primes targets as dependencies.
add_dependencies(tools
  ${EC_BIN_NAME}
)
