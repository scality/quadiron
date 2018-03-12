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

#ifndef __NTTEC_SIMD_128_H__
#define __NTTEC_SIMD_128_H__

#include <x86intrin.h>

namespace nttec {

/** The namespace simd contains functions that are accelerated by using SIMD
 *  operations over 128bits
 *
 *  Currently, it supports operations on 32-bit numbers
 */
namespace simd {

/* ==================== Essential Operations =================== */

const uint32_t F4 = 65537;
const uint32_t F3 = 257;
const __m128i F4_m128i = _mm_set1_epi32(65537);
const __m128i F4minus1_m128i = _mm_set1_epi32(65536);
const __m128i F3_m128i = _mm_set1_epi32(257);
const __m128i F3minus1_m128i = _mm_set1_epi32(256);

/** Return __uint128_t integer from a _m128i register */
static inline __uint128_t m128i_to_uint128(__m128i v)
{
    __uint128_t i;
    _mm_storeu_si128((__m128i*)&i, v);

    return i; // NOLINT(clang-analyzer-core.uninitialized.UndefReturn)
}

/** Perform a%card where a is a addition of two numbers whose elements are
 *  symbols of GF(card) */
inline __uint128_t mod_after_add(__m128i a, uint32_t card)
{
    const __m128i _card = _mm_set1_epi32(card);
    const __m128i _card_minus_1 = _mm_set1_epi32(card - 1);

    __m128i cmp = _mm_cmpgt_epi32(a, _card_minus_1);
    __m128i b = _mm_sub_epi32(a, _mm_and_si128(_card, cmp));

    return m128i_to_uint128(b);
}

/** Perform addition of two numbers a, b whose elements are of GF(card) */
inline __uint128_t add(__uint128_t a, __uint128_t b, uint32_t card = F4)
{
    __m128i _a = _mm_loadu_si128((__m128i*)&a);
    __m128i _b = _mm_loadu_si128((__m128i*)&b);
    __m128i c = _mm_add_epi32(_a, _b);

    // Modulo
    return mod_after_add(c, card);
}

/** Perform subtraction of a by b where a, b whose elements are symbols of
 *  GF(card)
 * sub(a, b) = a - b if a >= b, or
 *             card + a - b, otherwise
 */
inline __uint128_t sub(__uint128_t a, __uint128_t b, uint32_t card = F4)
{
    const __m128i _card = _mm_set1_epi32(card);

    __m128i _a = _mm_loadu_si128((__m128i*)&a);
    __m128i _b = _mm_loadu_si128((__m128i*)&b);

    __m128i cmp = _mm_cmpgt_epi32(_b, _a);
    __m128i _a1 = _mm_add_epi32(_a, _mm_and_si128(_card, cmp));

    __m128i c = _mm_sub_epi32(_a1, _b);

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
inline __uint128_t mod_after_multiply_f4(__m128i a)
{
    const __m128i mask = _mm_set1_epi32(F4 - 2);

    __m128i lo = _mm_and_si128(a, mask);

    __m128i a_shift = _mm_srli_si128(a, 2);
    __m128i hi = _mm_and_si128(a_shift, mask);

    __m128i cmp = _mm_cmpgt_epi32(hi, lo);
    __m128i _lo = _mm_add_epi32(lo, _mm_and_si128(F4_m128i, cmp));

    __m128i b = _mm_sub_epi32(_lo, hi);

    return m128i_to_uint128(b);
}

inline __uint128_t mod_after_multiply_f3(__m128i a)
{
    const __m128i mask = _mm_set1_epi32(F3 - 2);

    __m128i lo = _mm_and_si128(a, mask);

    __m128i a_shift = _mm_srli_si128(a, 1);
    __m128i hi = _mm_and_si128(a_shift, mask);

    __m128i cmp = _mm_cmpgt_epi32(hi, lo);
    __m128i _lo = _mm_add_epi32(lo, _mm_and_si128(F3_m128i, cmp));

    __m128i b = _mm_sub_epi32(_lo, hi);

    return m128i_to_uint128(b);
}

inline __uint128_t mul_f4(__uint128_t a, __uint128_t b)
{
    __m128i _a = _mm_loadu_si128((__m128i*)&a);
    __m128i _b = _mm_loadu_si128((__m128i*)&b);

    __m128i c = _mm_mullo_epi32(_a, _b);

    // filter elements of both of a & b = card-1
    __m128i cmp = _mm_and_si128(
        _mm_cmpeq_epi32(_a, F4minus1_m128i),
        _mm_cmpeq_epi32(_b, F4minus1_m128i));

    const __m128i one = _mm_set1_epi32(1);
    c = _mm_add_epi32(c, _mm_and_si128(one, cmp));

    // Modulo
    return mod_after_multiply_f4(c);
}

inline __uint128_t mul_f3(__uint128_t a, __uint128_t b)
{
    __m128i _a = _mm_loadu_si128((__m128i*)&a);
    __m128i _b = _mm_loadu_si128((__m128i*)&b);

    __m128i c = _mm_mullo_epi32(_a, _b);

    // filter elements of both of a & b = card-1
    __m128i cmp = _mm_and_si128(
        _mm_cmpeq_epi32(_a, F3minus1_m128i),
        _mm_cmpeq_epi32(_b, F3minus1_m128i));

    c = _mm_xor_si128(c, _mm_and_si128(F3_m128i, cmp));

    // Modulo
    return mod_after_multiply_f3(c);
}

/** Perform multiplication of two numbers a, b whose elements are of GF(card)
 *  where `card` is a prime Fermat number, i.e. card = Fx with x < 5
 *  Currently, it supports only for F3 and F4
 */
inline __uint128_t mul(__uint128_t a, __uint128_t b, uint32_t card = F4)
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
    const uint32_t a,
    uint32_t* src,
    uint32_t* dest,
    size_t len,
    uint32_t card = F4)
{
    const __m128i _coef = _mm_set1_epi32(a);
    __uint128_t coef = m128i_to_uint128(_coef);

    __uint128_t* _src = static_cast<__uint128_t*>(static_cast<void*>(src));
    __uint128_t* _dest = static_cast<__uint128_t*>(static_cast<void*>(dest));
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
            dest[i] = (uint32_t)((coef_64 * src[i]) % card);
        }
    }
}

inline void
add_two_bufs(uint32_t* src, uint32_t* dest, size_t len, uint32_t card = F4)
{
    __uint128_t* _src = static_cast<__uint128_t*>(static_cast<void*>(src));
    __uint128_t* _dest = static_cast<__uint128_t*>(static_cast<void*>(dest));
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
            uint32_t tmp = src[i] + dest[i];
            dest[i] = (tmp >= card) ? (tmp - card) : tmp;
        }
    }
}

inline void sub_two_bufs(
    uint32_t* bufa,
    uint32_t* bufb,
    uint32_t* res,
    size_t len,
    uint32_t card = F4)
{
    __uint128_t* _bufa = static_cast<__uint128_t*>(static_cast<void*>(bufa));
    __uint128_t* _bufb = static_cast<__uint128_t*>(static_cast<void*>(bufb));
    __uint128_t* _res = static_cast<__uint128_t*>(static_cast<void*>(res));
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

/* ==================== Operations for NF4 =================== */

inline __uint128_t expand16(uint16_t* arr, int n)
{
    // since n <= 4
    uint16_t _arr[4] = {0, 0, 0, 0};
    std::copy_n(arr, n, _arr);

    __m128i b = _mm_set_epi64(
        _mm_setzero_si64(), _mm_set_pi16(_arr[3], _arr[2], _arr[1], _arr[0]));

    return m128i_to_uint128(b);
}

inline __uint128_t expand32(uint32_t* arr, int n)
{
    // since n <= 4
    uint32_t _arr[4] = {0, 0, 0, 0};
    std::copy_n(arr, n, _arr);

    __m128i b = _mm_set_epi32(_arr[3], _arr[2], _arr[1], _arr[0]);

    return m128i_to_uint128(b);
}

/** Add buffer `y` to two halves of `x`. `x` is of length `n` */
inline void add_buf_to_two_bufs(int n, __uint128_t* x, __uint128_t* y)
{
    int i;
    int half = n / 2;
    __uint128_t* x_next = x + half;

    // add y to the first half of `x`
    for (i = 0; i < half; i++) {
        x[i] = add(x[i], y[i]);
    }

    // add y to the second half of `x`
    for (i = 0; i < half; i++) {
        x_next[i] = add(x_next[i], y[i]);
    }
}

inline void hadamard_mul(int n, __uint128_t* x, __uint128_t* y)
{
    int i;
    int half = n / 2;
    __uint128_t* x_next = x + half;

    // multiply y to the first half of `x`
    for (i = 0; i < half; i++) {
        x[i] = mul(x[i], y[i]);
    }

    // multiply y to the second half of `x`
    for (i = 0; i < half; i++) {
        x_next[i] = mul(x_next[i], y[i]);
    }
}

inline GroupedValues<__uint128_t> unpack(__uint128_t a, int n)
{
    uint32_t flag = 0;
    uint32_t ai[4];
    uint16_t bi[4] = {0, 0, 0, 0};
    __uint128_t values;
    int i;

    __m128i _a = _mm_loadu_si128((__m128i*)&a);
    ai[0] = _mm_extract_epi32(_a, 0);
    ai[1] = _mm_extract_epi32(_a, 1);
    ai[2] = _mm_extract_epi32(_a, 2);
    ai[3] = _mm_extract_epi32(_a, 3);
    for (i = 0; i < n; i++) {
        if (ai[i] == 65536)
            flag |= (1 << i);
        else
            bi[i] = (uint16_t)ai[i];
    }
    __m128i val = _mm_set_epi64(
        _mm_setzero_si64(), _mm_set_pi16(bi[3], bi[2], bi[1], bi[0]));
    _mm_storeu_si128((__m128i*)&values, val);

    GroupedValues<__uint128_t> b = {values, flag};

    return b;
}

inline __uint128_t pack(__uint128_t a)
{
    __m128i _a = _mm_loadu_si128((__m128i*)&a);
    __m128i b = _mm_set_epi32(
        _mm_extract_epi16(_a, 3),
        _mm_extract_epi16(_a, 2),
        _mm_extract_epi16(_a, 1),
        _mm_extract_epi16(_a, 0));

    return m128i_to_uint128(b);
}

inline __uint128_t pack(__uint128_t a, uint32_t flag)
{
    uint32_t b0, b1, b2, b3;
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

} // namespace simd
} // namespace nttec

#endif
