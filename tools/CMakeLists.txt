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
