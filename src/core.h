#ifndef __NTTEC_CORE_H__
#define __NTTEC_CORE_H__

#include <cstdint>

#include <gmpxx.h>

#include "big_int.h"

namespace nttec {

template <typename Type>
struct Double {
};
template <>
struct Double<uint32_t> {
    typedef uint64_t T;
};
template <>
struct Double<uint64_t> {
    typedef __uint128_t T;
};
template <>
struct Double<__uint128_t> {
    typedef uint256_t T;
};
template <>
struct Double<mpz_class> {
    typedef mpz_class T;
};

template <typename Type>
struct SignedDouble {
};
template <>
struct SignedDouble<uint32_t> {
    typedef int64_t T;
};
template <>
struct SignedDouble<uint64_t> {
    typedef __int128_t T;
};
template <>
struct SignedDouble<__uint128_t> {
    typedef int256_t T;
};
template <>
struct SignedDouble<mpz_class> {
    typedef mpz_class T;
};

template <typename T>
struct compT {
    T val;
    uint32_t flag;
};

typedef enum {
    NTTEC_EX_NOT_FOUND,
    NTTEC_EX_MAT_NOT_INVERTIBLE,
    NTTEC_EX_INVAL,
    NTTEC_EX_OVERFLOW,
    NTTEC_EX_IO,
    NTTEC_EX_DIV_BY_ZERO,
} NttecException;

} // namespace nttec

#endif
