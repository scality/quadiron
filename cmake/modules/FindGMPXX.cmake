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

# - Try to find C++ binding for GMP (GNU Multiple Precision Arithmetic Library)
#
# Once done this will define
#  GMPXX_FOUND        - System has C++ binding for GMP.
#  GMPXX_INCLUDE_DIRS - The include directories for the C++ binding of GMP.
#  GMPXX_LIBRARIES    - The libraries needed to use the C++ binding of GMP.

include(LibFindMacros)

# Dependencies.
libfind_package(GMPXX GMP)

# Use pkg-config to get hints about paths.
libfind_pkg_check_modules(GMPXX_PKGCONF gmpxx)

# Include directory.
find_path(GMPXX_INCLUDE_DIR
  NAMES gmpxx.h
  PATHS ${GMPXX_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself.
find_library(GMPXX_LIBRARY
  NAMES gmpxx
  PATHS ${GMPXX_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir/libraries variables and let libfind_process do the rest.
libfind_process(GMPXX)
