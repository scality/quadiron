# NTTEC

[![CircleCI][badgepub]](https://circleci.com/gh/vrancurel/nttec)

stuff

# Build from source

## Build

```sh
mkdir build
cd build
cmake -G 'Unix Makefiles' ..
make
```

### Targets

- `doc`: build the documentation using Doxygen
- `bench`: build the NTTEC benchmark (build mode "Release" is recommended)
- `shared`: build the NTTEC shared library
- `static`: build the NTTEC static library
- `unit_tests`: build the unit tests.
- `check`: run the test suite
- `tools`: build the NTTEC tools
- `benchmark`: run the NTTEC benchmark (build mode "Release" is recommended)
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

[badgepub]: https://circleci.com/gh/vrancurel/nttec.svg?style=shield&circle-token=:circle-token
