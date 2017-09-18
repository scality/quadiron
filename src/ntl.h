/* -*- mode: c++ -*- */
#pragma once

#include <iostream>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <gmpxx.h>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdint>
#include <vector>
#include <map>
#include <utility>

#include "big_int.h"

template<typename Type> struct Double {};
template<>           struct Double<uint32_t> {typedef uint64_t T;};
template<>           struct Double<uint64_t> {typedef __uint128_t T;};
template<>           struct Double<__uint128_t> {typedef uint256_t T;};
template<>           struct Double<mpz_class> {typedef mpz_class T;};

template<typename Type> struct SignedDouble {};
template<>           struct SignedDouble<uint32_t> {typedef int64_t T;};
template<>           struct SignedDouble<uint64_t> {typedef __int128_t T;};
template<>           struct SignedDouble<__uint128_t> {typedef int256_t T;};
template<>           struct SignedDouble<mpz_class> {typedef mpz_class T;};

typedef enum
  {
    NTL_EX_NOT_FOUND,
    NTL_EX_MAT_NOT_INVERTIBLE,
    NTL_EX_INVAL,
    NTL_EX_OVERFLOW,
    NTL_EX_IO,
    NTL_EX_DIV_BY_ZERO
  } NtlException;

#include "gf.h"
#include "gfp.h"
#include "gf2n.h"
#include "misc.h"
#include "config.h"
#include "vec.h"
#include "vvec.h"
#include "v2vec.h"
#include "vmvec.h"
#include "mat.h"
#include "fft.h"
#include "fft2.h"
#include "fftn.h"
#include "fftln.h"
#include "fft2k.h"
#include "fftpf.h"
#include "poly.h"
#include "fec.h"
#include "fecgf2nrs.h"
#include "fecfntrs.h"
#include "fecgf2nfftrs.h"
