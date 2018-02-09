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
