
#ifndef __NTL_H__
#define __NTL_H__ 1

#include <iostream>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <gmpxx.h>

template<typename Type> struct Double {};
template<>           struct Double<uint32_t> {typedef uint64_t T;};
template<>           struct Double<uint64_t> {typedef unsigned __int128 T;};
template<>           struct Double<mpz_class> {typedef mpz_class T;};

template<typename Type> struct SignedDouble {};
template<>           struct SignedDouble<uint32_t> {typedef int64_t T;};
template<>           struct SignedDouble<uint64_t> {typedef __int128 T;};
template<>           struct SignedDouble<mpz_class> {typedef mpz_class T;};

typedef enum
  {
    NTL_EX_NOT_FOUND,
    NTL_EX_MAT_NOT_INVERTIBLE,
    NTL_EX_INVAL,
  } NtlException;

#include "gf.h"
#include "gfp.h"
#include "gf2n.h"
#include "vec.h"
#include "mat.h"
#include "misc.h"

#endif
