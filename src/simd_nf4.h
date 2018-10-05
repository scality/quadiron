/*
 * Copyright 2017-2018 Scality
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __QUAD_SIMD_NF4_H__
#define __QUAD_SIMD_NF4_H__

#include <x86intrin.h>

#include <simd/simd.h>

namespace quadiron {
namespace simd {

typedef uint32_t aint32 __attribute__((aligned(ALIGNMENT)));
typedef __uint128_t NF4Type;

/** Return NF4Type integer from a _m128i register */
static inline NF4Type m128i_to_uint128(__m128i v)
{
    NF4Type i;
    _mm_store_si128((__m128i*)&i, v);

    return i; // NOLINT(clang-analyzer-core.uninitialized.UndefReturn)
}

inline NF4Type expand16(uint16_t* arr, int n)
{
    // since n <= 4
    uint16_t _arr[4] __attribute__((aligned(ALIGNMENT))) = {0, 0, 0, 0};
    std::copy_n(arr, n, _arr);

    __m128i b = _mm_set_epi16(0, 0, 0, 0, _arr[3], _arr[2], _arr[1], _arr[0]);

    return m128i_to_uint128(b);
}

inline NF4Type expand32(uint32_t* arr, int n)
{
    // since n <= 4
    uint32_t _arr[4] __attribute__((aligned(simd::ALIGNMENT))) = {0, 0, 0, 0};
    std::copy_n(arr, n, _arr);

    __m128i b = _mm_set_epi32(_arr[3], _arr[2], _arr[1], _arr[0]);

    return m128i_to_uint128(b);
}

inline GroupedValues<__uint128_t> unpack(__uint128_t a, int n)
{
    uint16_t ai[8];
    NF4Type values;

    __m128i _a = _mm_loadu_si128((__m128i*)&a);
    ai[0] = _mm_extract_epi16(_a, 0);
    ai[1] = _mm_extract_epi16(_a, 1);
    ai[2] = _mm_extract_epi16(_a, 2);
    ai[3] = _mm_extract_epi16(_a, 3);
    ai[4] = _mm_extract_epi16(_a, 4);
    ai[5] = _mm_extract_epi16(_a, 5);
    ai[6] = _mm_extract_epi16(_a, 6);
    ai[7] = _mm_extract_epi16(_a, 7);

    const uint32_t flag =
        ai[1] | (!!ai[3] << 1u) | (!!ai[5] << 2u) | (!!ai[7] << 3u);

    __m128i val = _mm_set_epi16(0, 0, 0, 0, ai[6], ai[4], ai[2], ai[0]);
    _mm_store_si128((__m128i*)&values, val);

    GroupedValues<__uint128_t> b = {values, flag};

    return b;
}

inline void unpack(__uint128_t a, GroupedValues<__uint128_t>& b, int n)
{
    uint16_t ai[8];
    NF4Type values;

    __m128i _a = _mm_loadu_si128((__m128i*)&a);
    ai[0] = _mm_extract_epi16(_a, 0);
    ai[1] = _mm_extract_epi16(_a, 1);
    ai[2] = _mm_extract_epi16(_a, 2);
    ai[3] = _mm_extract_epi16(_a, 3);
    ai[4] = _mm_extract_epi16(_a, 4);
    ai[5] = _mm_extract_epi16(_a, 5);
    ai[6] = _mm_extract_epi16(_a, 6);
    ai[7] = _mm_extract_epi16(_a, 7);

    const uint32_t flag =
        ai[1] | (!!ai[3] << 1u) | (!!ai[5] << 2u) | (!!ai[7] << 3u);

    __m128i val = _mm_set_epi16(0, 0, 0, 0, ai[6], ai[4], ai[2], ai[0]);
    _mm_store_si128((__m128i*)&values, val);

    b.flag = flag;
    b.values = values; // NOLINT(clang-analyzer-core.uninitialized.Assign)
}

inline NF4Type pack(__uint128_t a)
{
    __m128i _a = _mm_loadu_si128((__m128i*)&a);
    __m128i b = _mm_set_epi32(
        _mm_extract_epi16(_a, 3),
        _mm_extract_epi16(_a, 2),
        _mm_extract_epi16(_a, 1),
        _mm_extract_epi16(_a, 0));

    return m128i_to_uint128(b);
}

inline NF4Type pack(__uint128_t a, uint32_t flag)
{
    aint32 b0, b1, b2, b3;
    __m128i _a = _mm_loadu_si128((__m128i*)&a);

    if (flag & 1)
        b0 = 65536;
    else
        b0 = _mm_extract_epi16(_a, 0);
    flag >>= 1;
    if (flag & 1)
        b1 = 65536;
    else
        b1 = _mm_extract_epi16(_a, 1);
    flag >>= 1;
    if (flag & 1)
        b2 = 65536;
    else
        b2 = _mm_extract_epi16(_a, 2);
    flag >>= 1;
    if (flag & 1)
        b3 = 65536;
    else
        b3 = _mm_extract_epi16(_a, 3);

    __m128i b = _mm_set_epi32(b3, b2, b1, b0);

    return m128i_to_uint128(b);
}

/* ================= Basic operations for NF4 ================= */

#if defined(__AVX2__)

inline VecType CAST_TO_DOUBLE(HalfVecType x)
{
    return _mm256_castsi128_si256(x);
}

inline void STORE_LOW(HalfVecType* address, VecType reg)
{
    _mm_store_si128(address, _mm256_castsi256_si128(reg));
}

inline NF4Type add(NF4Type a, NF4Type b)
{
    HalfVecType res;
    VecType _a = CAST_TO_DOUBLE((HalfVecType)a);
    VecType _b = CAST_TO_DOUBLE((HalfVecType)b);
    STORE_LOW(&res, ADD_MOD(_a, _b, F4));
    return (NF4Type)res;
}

inline NF4Type sub(NF4Type a, NF4Type b)
{
    HalfVecType res;
    VecType _a = CAST_TO_DOUBLE((HalfVecType)a);
    VecType _b = CAST_TO_DOUBLE((HalfVecType)b);
    STORE_LOW(&res, SUB_MOD(_a, _b, F4));
    return (NF4Type)res;
}

inline NF4Type mul(NF4Type a, NF4Type b)
{
    HalfVecType res;
    VecType _a = CAST_TO_DOUBLE((HalfVecType)a);
    VecType _b = CAST_TO_DOUBLE((HalfVecType)b);
    STORE_LOW(&res, MULFULL_MOD(_a, _b, F4));
    return (NF4Type)res;
}

inline void
add_buf_to_two_bufs_rem(unsigned n, NF4Type* x, NF4Type* x_half, NF4Type* y)
{
    // add last _y[] to x and x_next
    HalfVecType* _x = reinterpret_cast<HalfVecType*>(x);
    HalfVecType* _x_half = reinterpret_cast<HalfVecType*>(x_half);
    HalfVecType* _y = reinterpret_cast<HalfVecType*>(y);
    for (unsigned i = 0; i < n; ++i) {
        VecType _x_p = CAST_TO_DOUBLE(_x[i]);
        VecType _x_next_p = CAST_TO_DOUBLE(_x_half[i]);
        VecType _y_p = CAST_TO_DOUBLE(_y[i]);

        STORE_LOW(_x + i, ADD_MOD(_x_p, _y_p, F4));
        STORE_LOW(_x_half + i, ADD_MOD(_x_next_p, _y_p, F4));
    }
}

inline void hadamard_mul_rem(unsigned n, NF4Type* x, NF4Type* y)
{
    HalfVecType* _x = reinterpret_cast<HalfVecType*>(x);
    HalfVecType* _y = reinterpret_cast<HalfVecType*>(y);
    for (unsigned i = 0; i < n; ++i) {
        VecType _x_p = CAST_TO_DOUBLE(_x[i]);
        VecType _y_p = CAST_TO_DOUBLE(_y[i]);

        STORE_LOW(_x + i, MULFULL_MOD(_x_p, _y_p, F4));
    }
}

inline void
hadamard_mul_doubled_rem(unsigned n, NF4Type* x, NF4Type* x_half, NF4Type* y)
{
    HalfVecType* _x = reinterpret_cast<HalfVecType*>(x);
    HalfVecType* _x_half = reinterpret_cast<HalfVecType*>(x_half);
    HalfVecType* _y = reinterpret_cast<HalfVecType*>(y);
    for (unsigned i = 0; i < n; ++i) {
        VecType _x_p = CAST_TO_DOUBLE(_x[i]);
        VecType _x_next_p = CAST_TO_DOUBLE(_x_half[i]);
        VecType _y_p = CAST_TO_DOUBLE(_y[i]);

        STORE_LOW(_x + i, MULFULL_MOD(_x_p, _y_p, F4));
        STORE_LOW(_x_half + i, MULFULL_MOD(_x_next_p, _y_p, F4));
    }
}

#elif defined(__SSE4_1__)

inline NF4Type add(NF4Type a, NF4Type b)
{
    VecType res;
    STORE(&res, ADD_MOD((VecType)a, (VecType)b, F4));
    return (NF4Type)res;
}

inline NF4Type sub(NF4Type a, NF4Type b)
{
    VecType res;
    STORE(&res, SUB_MOD((VecType)a, (VecType)b, F4));
    return (NF4Type)res;
}

inline NF4Type mul(NF4Type a, NF4Type b)
{
    VecType res;
    STORE(&res, MULFULL_MOD((VecType)a, (VecType)b, F4));
    return (NF4Type)res;
}

inline void
add_buf_to_two_bufs_rem(unsigned n, NF4Type* x, NF4Type* x_half, NF4Type* y)
{
    // do nothing
}

inline void hadamard_mul_rem(unsigned n, NF4Type* x, NF4Type* y)
{
    // do nothing
}

inline void
hadamard_mul_doubled_rem(unsigned n, NF4Type* x, NF4Type* x_half, NF4Type* y)
{
    // do nothing
}

#endif

/* ==================== Operations for NF4 =================== */

/** Add buffer `y` to two halves of `x`. `x` is of length `n` */
inline void add_buf_to_two_bufs(unsigned n, NF4Type* _x, NF4Type* _y)
{
    unsigned i;
    VecType* x = reinterpret_cast<VecType*>(_x);
    VecType* y = reinterpret_cast<VecType*>(_y);

    const unsigned ratio = sizeof(*x) / sizeof(*_x);
    const unsigned half_len = n / 2;
    const unsigned vec_len = half_len / ratio;
    const unsigned num_len = vec_len * ratio;
    const unsigned rem_len = half_len - num_len;

    NF4Type* x_half = _x + half_len;
    VecType* x_next = reinterpret_cast<VecType*>(x_half);

    // add y to the first half of `x`
    for (i = 0; i < vec_len; ++i) {
        x[i] = ADD_MOD(x[i], y[i], F4);
    }

    // add y to the second half of `x`
    for (i = 0; i < vec_len; ++i) {
        x_next[i] = ADD_MOD(x_next[i], y[i], F4);
    }

    if (rem_len > 0) {
        add_buf_to_two_bufs_rem(
            rem_len, _x + num_len, x_half + num_len, _y + num_len);
    }
}

inline void hadamard_mul(unsigned n, NF4Type* _x, NF4Type* _y)
{
    unsigned i;
    VecType* x = reinterpret_cast<VecType*>(_x);
    VecType* y = reinterpret_cast<VecType*>(_y);

    const unsigned ratio = sizeof(*x) / sizeof(*_x);
    const unsigned vec_len = n / ratio;
    const unsigned num_len = vec_len * ratio;
    const unsigned rem_len = n - num_len;

    // multiply y to the first half of `x`
    for (i = 0; i < vec_len; ++i) {
        x[i] = MULFULL_MOD(x[i], y[i], F4);
    }

    if (rem_len > 0) {
        // add last _y[] to x
        hadamard_mul_rem(rem_len, _x + num_len, _y + num_len);
    }
}

} // namespace simd
} // namespace quadiron

#endif
