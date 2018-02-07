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

- `ntl_bench`: build the NTL benchmark (build mode "Release" is recommended)
- `ntl_shared`: build the NTL shared library
- `ntl_static`: build the NTL static library
- `ntl_test`: build the test driver.
- `check`: run the test suite
- `tools`: build the NTL tools
- `benchmark`: run the NTL benchmark (build mode "Release" is recommended)
