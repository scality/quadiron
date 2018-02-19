#ifndef __NTTEC_CORE_H__
#define __NTTEC_CORE_H__

#include <cstdint>

#include <gmpxx.h>

#include "big_int.h"

namespace nttec {

template <typename Type>
struct DoubleSize {
};
template <>
struct DoubleSize<uint32_t> {
    typedef uint64_t T;
};
template <>
struct DoubleSize<uint64_t> {
    typedef __uint128_t T;
};
template <>
struct DoubleSize<__uint128_t> {
    typedef UInt256 T;
};
template <>
struct DoubleSize<mpz_class> {
    typedef mpz_class T;
};

template <typename Type>
struct SignedDoubleSize {
};
template <>
struct SignedDoubleSize<uint32_t> {
    typedef int64_t T;
};
template <>
struct SignedDoubleSize<uint64_t> {
    typedef __int128_t T;
};
template <>
struct SignedDoubleSize<__uint128_t> {
    typedef Int256 T;
};
template <>
struct SignedDoubleSize<mpz_class> {
    typedef mpz_class T;
};

/** A group of values stored as one.
 *
 * This allows faster processing, as the values can be processed as one.
 */
template <typename T>
struct GroupedValues {
    // A group of several values.
    T values;

    /** Per-value flags.
     *
     * For now, only the first n bits (n being the number of values stored) are
     * used.
     * When the bit is set, the corresponding value should be 0 and that means
     * that the real value is Fn-1.
     */
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

/** Return the version string of NTTEC.
 *
 * The version string has the form MAJOR.MINOR.PATCH-REVISION, where '-REVISION'
 * is optional (only present for development version).
 *
 * @return the version string.
 */
const char* get_version();

} // namespace nttec

#endif
