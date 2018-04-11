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

#ifndef __NTTEC_SIMD_128_U32_H__
#define __NTTEC_SIMD_128_U32_H__

#include <x86intrin.h>

namespace nttec {

/** The namespace simd contains functions that are accelerated by using SIMD
 *  operations over 128bits
 *
 *  It supports operations on 32-bit numbers
 */
namespace simd {
namespace u32 {

typedef __m128i am128i __attribute__((aligned(ALIGN_SIZE)));

/* ==================== Essential Operations =================== */

const aint32 F4 = 65537;
const aint32 F3 = 257;

// Disable `cert-err58-cpp` on these: AFAIK they cannot throw.
// (probably a false positive present in Clang 5 and fixed in Clang 6).
const am128i F4_m128i = _mm_set1_epi32(65537);       // NOLINT(cert-err58-cpp)
const am128i F4minus1_m128i = _mm_set1_epi32(65536); // NOLINT(cert-err58-cpp)
const am128i F3_m128i = _mm_set1_epi32(257);         // NOLINT(cert-err58-cpp)
const am128i F3minus1_m128i = _mm_set1_epi32(256);   // NOLINT(cert-err58-cpp)

/** Return aint128 integer from a _m128i register */
static inline aint128 m128i_to_uint128(am128i v)
{
    aint128 i;
    _mm_store_si128((am128i*)&i, v);

    return i; // NOLINT(clang-analyzer-core.uninitialized.UndefReturn)
}

/** Perform a%card where a is a addition of two numbers whose elements are
 *  symbols of GF(card) */
inline aint128 mod_after_add(am128i a, aint32 card)
{
    const am128i _card = _mm_set1_epi32(card);
    const am128i _card_minus_1 = _mm_set1_epi32(card - 1);

    am128i cmp = _mm_cmpgt_epi32(a, _card_minus_1);
    am128i b = _mm_sub_epi32(a, _mm_and_si128(_card, cmp));

    return m128i_to_uint128(b);
}

/** Perform addition of two numbers a, b whose elements are of GF(card) */
inline aint128 add(aint128 a, aint128 b, aint32 card = F4)
{
    am128i _a = _mm_load_si128((am128i*)&a);
    am128i _b = _mm_load_si128((am128i*)&b);
    am128i c = _mm_add_epi32(_a, _b);

    // Modulo
    return mod_after_add(c, card);
}

/** Perform subtraction of a by b where a, b whose elements are symbols of
 *  GF(card)
 * sub(a, b) = a - b if a >= b, or
 *             card + a - b, otherwise
 */
inline aint128 sub(aint128 a, aint128 b, aint32 card = F4)
{
    const am128i _card = _mm_set1_epi32(card);

    am128i _a = _mm_load_si128((am128i*)&a);
    am128i _b = _mm_load_si128((am128i*)&b);

    am128i cmp = _mm_cmpgt_epi32(_b, _a);
    am128i _a1 = _mm_add_epi32(_a, _mm_and_si128(_card, cmp));

    am128i c = _mm_sub_epi32(_a1, _b);

    return m128i_to_uint128(c);
}

/** Perform a%card where a is a multiplication of two numbers whose elements are
 *  symbols of GF(F4)
 *
 *  We find v in a = u * card + v
 *  a is expressed also as: a = hi * (card-1) + lo
 *  where hi and lo is 16-bit for F4 (or 8-bit for F3) high and low parts of a
 *  hence, v = (lo - hi) % F4
 *      v = lo - hi, if lo >= hi
 *          or
 *          F4 + lo - hi, otherwise
 */
inline aint128 mod_after_multiply_f4(am128i a)
{
    const am128i mask = _mm_set1_epi32(F4 - 2);

    am128i lo = _mm_and_si128(a, mask);

    am128i a_shift = _mm_srli_si128(a, 2);
    am128i hi = _mm_and_si128(a_shift, mask);

    am128i cmp = _mm_cmpgt_epi32(hi, lo);
    am128i _lo = _mm_add_epi32(lo, _mm_and_si128(F4_m128i, cmp));

    am128i b = _mm_sub_epi32(_lo, hi);

    return m128i_to_uint128(b);
}

inline aint128 mod_after_multiply_f3(am128i a)
{
    const am128i mask = _mm_set1_epi32(F3 - 2);

    am128i lo = _mm_and_si128(a, mask);

    am128i a_shift = _mm_srli_si128(a, 1);
    am128i hi = _mm_and_si128(a_shift, mask);

    am128i cmp = _mm_cmpgt_epi32(hi, lo);
    am128i _lo = _mm_add_epi32(lo, _mm_and_si128(F3_m128i, cmp));

    am128i b = _mm_sub_epi32(_lo, hi);

    return m128i_to_uint128(b);
}

inline aint128 mul_f4(aint128 a, aint128 b)
{
    am128i _a = _mm_load_si128((am128i*)&a);
    am128i _b = _mm_load_si128((am128i*)&b);

    am128i c = _mm_mullo_epi32(_a, _b);

    // filter elements of both of a & b = card-1
    am128i cmp = _mm_and_si128(
        _mm_cmpeq_epi32(_a, F4minus1_m128i),
        _mm_cmpeq_epi32(_b, F4minus1_m128i));

    const am128i one = _mm_set1_epi32(1);
    c = _mm_add_epi32(c, _mm_and_si128(one, cmp));

    // Modulo
    return mod_after_multiply_f4(c);
}

inline aint128 mul_f3(aint128 a, aint128 b)
{
    am128i _a = _mm_load_si128((am128i*)&a);
    am128i _b = _mm_load_si128((am128i*)&b);

    am128i c = _mm_mullo_epi32(_a, _b);

    // filter elements of both of a & b = card-1
    am128i cmp = _mm_and_si128(
        _mm_cmpeq_epi32(_a, F3minus1_m128i),
        _mm_cmpeq_epi32(_b, F3minus1_m128i));

    c = _mm_xor_si128(c, _mm_and_si128(F4_m128i, cmp));

    // Modulo
    return mod_after_multiply_f3(c);
}

/** Perform multiplication of two numbers a, b whose elements are of GF(card)
 *  where `card` is a prime Fermat number, i.e. card = Fx with x < 5
 *  Currently, it supports only for F3 and F4
 */
inline aint128 mul(aint128 a, aint128 b, aint32 card = F4)
{
    assert(card == F4 || card == F3);
    if (card == F4)
        return mul_f4(a, b);
    return mul_f3(a, b);
}

/** Perform a multiplication of a coefficient `a` to each element of `src` and
 *  add result to correspondent element of `dest`
 */
inline void mul_coef_to_buf(
    const aint32 a,
    aint32* src,
    aint32* dest,
    size_t len,
    aint32 card = F4)
{
    const am128i _coef = _mm_set1_epi32(a);
    aint128 coef = m128i_to_uint128(_coef);

    aint128* _src = static_cast<aint128*>(static_cast<void*>(src));
    aint128* _dest = static_cast<aint128*>(static_cast<void*>(dest));
    size_t _len = len / 4;
    size_t _last_len = len - _len * 4;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform multiplication
        _dest[i] = mul(coef, _src[i], card);
    }
    if (_last_len > 0) {
        uint64_t coef_64 = (uint64_t)a;
        for (i = _len * 4; i < len; i++) {
            // perform multiplication
            dest[i] = (aint32)((coef_64 * src[i]) % card);
        }
    }
}

inline void
add_two_bufs(aint32* src, aint32* dest, size_t len, aint32 card = F4)
{
    aint128* _src = reinterpret_cast<aint128*>(src);
    aint128* _dest = reinterpret_cast<aint128*>(dest);
    size_t _len = len / 4;
    size_t _last_len = len - _len * 4;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform addition
        _dest[i] = add(_src[i], _dest[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * 4; i < len; i++) {
            // perform addition
            aint32 tmp = src[i] + dest[i];
            dest[i] = (tmp >= card) ? (tmp - card) : tmp;
        }
    }
}

inline void sub_two_bufs(
    aint32* bufa,
    aint32* bufb,
    aint32* res,
    size_t len,
    aint32 card = F4)
{
    aint128* _bufa = static_cast<aint128*>(static_cast<void*>(bufa));
    aint128* _bufb = static_cast<aint128*>(static_cast<void*>(bufb));
    aint128* _res = static_cast<aint128*>(static_cast<void*>(res));
    size_t _len = len / 4;
    size_t _last_len = len - _len * 4;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform subtraction
        _res[i] = sub(_bufa[i], _bufb[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * 4; i < len; i++) {
            // perform subtraction
            if (bufa[i] >= bufb[i])
                res[i] = bufa[i] - bufb[i];
            else
                res[i] = card - (bufb[i] - bufa[i]);
        }
    }
}

inline void
mul_two_bufs(aint32* src, aint32* dest, size_t len, aint32 card = F4)
{
    aint128* _src = reinterpret_cast<aint128*>(src);
    aint128* _dest = reinterpret_cast<aint128*>(dest);
    size_t _len = len / 4;
    size_t _last_len = len - _len * 4;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform multiplicaton
        _dest[i] = mul(_src[i], _dest[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * 4; i < len; i++) {
            // perform multiplicaton
            dest[i] = uint32_t((uint64_t(src[i]) * dest[i]) % card);
        }
    }
}

/* ==================== Operations for NF4 =================== */

/** Add buffer `y` to two halves of `x`. `x` is of length `n` */
inline void add_buf_to_two_bufs(int n, aint128* x, aint128* y, aint32 card = F4)
{
    int i;
    const int half = n / 2;
    aint128* x_next = x + half;

    // add y to the first half of `x`
    for (i = 0; i < half; i++) {
        x[i] = add(x[i], y[i], card);
    }

    // add y to the second half of `x`
    for (i = 0; i < half; i++) {
        x_next[i] = add(x_next[i], y[i], card);
    }
}

inline void hadamard_mul(int n, aint128* x, aint128* y)
{
    int i;
    const int half = n / 2;
    aint128* x_next = x + half;

    // multiply y to the first half of `x`
    for (i = 0; i < half; i++) {
        x[i] = mul(x[i], y[i]);
    }

    // multiply y to the second half of `x`
    for (i = 0; i < half; i++) {
        x_next[i] = mul(x_next[i], y[i]);
    }
}

inline aint128 expand16(aint16* arr, int n)
{
    // since n <= 4
    uint16_t _arr[4] __attribute__((aligned(ALIGN_SIZE))) = {0, 0, 0, 0};
    std::copy_n(arr, n, _arr);

    am128i b = _mm_set_epi64(
        _mm_setzero_si64(), _mm_set_pi16(_arr[3], _arr[2], _arr[1], _arr[0]));

    return m128i_to_uint128(b);
}

inline aint128 expand32(aint32* arr, int n)
{
    // since n <= 4
    uint32_t _arr[4] __attribute__((aligned(ALIGN_SIZE))) = {0, 0, 0, 0};
    std::copy_n(arr, n, _arr);

    am128i b = _mm_set_epi32(_arr[3], _arr[2], _arr[1], _arr[0]);

    return m128i_to_uint128(b);
}

inline GroupedValues<aint128> unpack(aint128 a, int n)
{
    aint32 flag = 0;
    uint32_t ai[4] __attribute__((aligned(ALIGN_SIZE)));
    uint32_t bi[4] __attribute__((aligned(ALIGN_SIZE))) = {0, 0, 0, 0};
    aint128 values;
    int i;

    am128i _a = _mm_load_si128((am128i*)&a);
    ai[0] = _mm_extract_epi32(_a, 0);
    ai[1] = _mm_extract_epi32(_a, 1);
    ai[2] = _mm_extract_epi32(_a, 2);
    ai[3] = _mm_extract_epi32(_a, 3);
    for (i = 0; i < n; i++) {
        if (ai[i] == 65536)
            flag |= (1 << i);
        else
            bi[i] = (aint16)ai[i];
    }
    am128i val = _mm_set_epi64(
        _mm_setzero_si64(), _mm_set_pi16(bi[3], bi[2], bi[1], bi[0]));
    _mm_store_si128((am128i*)&values, val);

    GroupedValues<aint128> b = {values, flag};

    return b;
}

inline aint128 pack(aint128 a)
{
    am128i _a = _mm_load_si128((am128i*)&a);
    am128i b = _mm_set_epi32(
        _mm_extract_epi16(_a, 3),
        _mm_extract_epi16(_a, 2),
        _mm_extract_epi16(_a, 1),
        _mm_extract_epi16(_a, 0));

    return m128i_to_uint128(b);
}

inline aint128 pack(aint128 a, aint32 flag)
{
    aint32 b0, b1, b2, b3;
    am128i _a = _mm_load_si128((am128i*)&a);

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

    am128i b = _mm_set_epi32(b3, b2, b1, b0);

    return m128i_to_uint128(b);
}

} // namespace u32
} // namespace simd
} // namespace nttec

#endif
