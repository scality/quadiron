# NTL

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
- `ntl_bench`: build the NTL benchmark (build mode "Release" is recommended)
- `ntl_shared`: build the NTL shared library
- `ntl_static`: build the NTL static library
- `ntl_test`: build the test driver.
- `check`: run the test suite
- `tools`: build the NTL tools
- `benchmark`: run the NTL benchmark (build mode "Release" is recommended)

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
