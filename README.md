# QuadIron

[![CircleCI][badgepub]](https://circleci.com/gh/scality/quadiron)

The QuadIron library is a C++ library, written in C++11, that offers a
streaming API to use the different flavors of NTT-based erasure codes.

The library focuses primarily on high fragmentation i.e. n >> k.

It includes general modular arithmetic routines, algorithms for
manipulating rings of integers modulo n, finite fields (including
binary, prime and (non binary) extension fields), polynomial
operations, implementations of different flavors of discrete Fourier
transforms and forward error correction (FEC) codes algorithms.

The library also includes an abstraction for writing systematic and
non-systematic codes, although for applications requiring high
fragmentation systematic codes are not necessarily useful.

# Build from source

## Build

```sh
mkdir build
cd build
cmake -G 'Unix Makefiles' ..
make
```

**WARNING**: If you are compiling the code with Clang, don't use a 4.X version:
there is [a bug](https://bugs.llvm.org/show_bug.cgi?id=36723) that affects SIMD
code in QuadIron.

### Targets

- `doc`: build the documentation using Doxygen
- `bench`: build the QuadIron benchmark (build mode "Release" is recommended)
- `shared`: build the QuadIron shared library
- `static`: build the QuadIron static library
- `unit_tests`: build the unit tests.
- `check`: run the test suite
- `benchmark`: run the QuadIron benchmark (build mode "Release" is recommended)
- `package`: generate a binary installer
- `package_source`: generate a source installer (a tarball with the sources)
- `install`: install the library in `CMAKE_INSTALL_PREFIX`.
- `uninstall`: uninstall the library.

### Code coverage

By default the code coverage is not enabled.

To generate the code coverage reports:
1. set the option `ENABLE_COVERAGE` to `ON`.
2. recompile with `make`.
3. run the tests with `make check`.
4. extract coverage data with `make gcov`.
5. generate a report with `make lcov`.
6. open the report in `${CMAKE_BINARY_DIR}/lcov/html/all_targets`.

Note that, even though code coverage is supported by both Clang and GCC, result
with GCC seems more reliable (not surprising as we are using gcov).

### SIMD vectorization

QuadIron can be accelerated using SIMD vectorization. This is controled by the
`USE_SIMD` parameter, that can have one of the following values:
- **OFF** (default value): no SIMD vectorisation (except the one done by the
  compiler)
- **ON**: select the best SIMD instructions set supported QuadIron and the
  machine
- **SSE**: use SSE4.1 SIMD instructions
- **AVX**: use AVX2 SIMD instructions

[badgepub]: https://circleci.com/gh/scality/quadiron.svg?style=svg
