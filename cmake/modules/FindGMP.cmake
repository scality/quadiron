# - Try to find GMP (GNU Multiple Precision Arithmetic Library)
#
# Once done this will define
#  GMP_FOUND        - System has GMP.
#  GMP_INCLUDE_DIRS - The GMP include directories.
#  GMP_LIBRARIES    - The libraries needed to use GMP.

include(LibFindMacros)

# Use pkg-config to get hints about paths.
libfind_pkg_check_modules(GMP_PKGCONF gmp)

# Include directory.
find_path(GMP_INCLUDE_DIR
  NAMES gmp.h
  PATHS ${GMP_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself.
find_library(GMP_LIBRARY
  NAMES gmp
  PATHS ${GMP_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir/libraries variables and let libfind_process do the rest.
libfind_process(GMP)
