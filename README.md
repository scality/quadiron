# NTTEC

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
- `nttec_bench`: build the NTTEC benchmark (build mode "Release" is recommended)
- `nttec_shared`: build the NTTEC shared library
- `nttec_static`: build the NTTEC static library
- `nttec_test`: build the test driver.
- `check`: run the test suite
- `tools`: build the NTTEC tools
- `benchmark`: run the NTTEC benchmark (build mode "Release" is recommended)

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
