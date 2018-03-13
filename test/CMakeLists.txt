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

######################
# Erasure Coding Tests
######################

add_executable(${EC_DRIVER}
  ${CMAKE_CURRENT_SOURCE_DIR}/ec_driver.cpp
)
add_coverage(${EC_DRIVER})

# Link librairies.
target_link_libraries(${EC_DRIVER}
  ${STATIC_LIB}
)

# Tests
add_test(
  NAME ec_test
  COMMAND ${SCRIPT_DIR}/test_ec.sh "$<TARGET_FILE:${EC_DRIVER}>"
)
