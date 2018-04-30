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

#ifndef __NTTEC_SIMD_256_U32_H__
#define __NTTEC_SIMD_256_U32_H__

#include <x86intrin.h>

namespace nttec {

/** The namespace simd contains functions that are accelerated by using SIMD
 *  operations over 256bits
 *
 *  It supports operations on 32-bit numbers
 */
namespace simd {
namespace u32 {

typedef __m256i m256i;

/* ==================== Essential Operations =================== */

// Disable `cert-err58-cpp` on these: AFAIK they cannot throw.
// (probably a false positive present in Clang 5 and fixed in Clang 6).
const m256i F4_m256i = _mm256_set1_epi32(65537);       // NOLINT(cert-err58-cpp)
const m256i F4minus1_m256i = _mm256_set1_epi32(65536); // NOLINT(cert-err58-cpp)
const m256i F3_m256i = _mm256_set1_epi32(257);         // NOLINT(cert-err58-cpp)
const m256i F3minus1_m256i = _mm256_set1_epi32(256);   // NOLINT(cert-err58-cpp)

/** Perform a%card where a is a addition of two numbers whose elements are
 *  symbols of GF(card) */
inline m256i mod_after_add(m256i a, aint32 card)
{
    const m256i _card = _mm256_set1_epi32(card);
    const m256i _card_minus_1 = _mm256_set1_epi32(card - 1);

    m256i cmp = _mm256_cmpgt_epi32(a, _card_minus_1);
    m256i b = _mm256_sub_epi32(a, _mm256_and_si256(_card, cmp));

    return b;
}

/** Perform addition of two numbers a, b whose elements are of GF(card) */
inline m256i add(m256i a, m256i b, aint32 card = F4)
{
    m256i _a = _mm256_load_si256(&a);
    m256i _b = _mm256_load_si256(&b);
    m256i c = _mm256_add_epi32(_a, _b);

    // Modulo
    return mod_after_add(c, card);
}

/** Perform subtraction of a by b where a, b whose elements are symbols of
 *  GF(card)
 * sub(a, b) = a - b if a >= b, or
 *             card + a - b, otherwise
 */
inline m256i sub(m256i a, m256i b, aint32 card = F4)
{
    const m256i _card = _mm256_set1_epi32(card);

    m256i _a = _mm256_load_si256(&a);
    m256i _b = _mm256_load_si256(&b);

    m256i cmp = _mm256_cmpgt_epi32(_b, _a);
    m256i _a1 = _mm256_add_epi32(_a, _mm256_and_si256(_card, cmp));

    return _mm256_sub_epi32(_a1, _b);
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
inline m256i mod_after_multiply_f4(m256i a)
{
    const m256i mask = _mm256_set1_epi32(F4 - 2);

    m256i lo = _mm256_and_si256(a, mask);

    m256i a_shift = _mm256_srli_si256(a, 2);
    m256i hi = _mm256_and_si256(a_shift, mask);

    m256i cmp = _mm256_cmpgt_epi32(hi, lo);
    m256i _lo = _mm256_add_epi32(lo, _mm256_and_si256(F4_m256i, cmp));

    return _mm256_sub_epi32(_lo, hi);
}

inline m256i mod_after_multiply_f3(m256i a)
{
    const m256i mask = _mm256_set1_epi32(F3 - 2);

    m256i lo = _mm256_and_si256(a, mask);

    m256i a_shift = _mm256_srli_si256(a, 1);
    m256i hi = _mm256_and_si256(a_shift, mask);

    m256i cmp = _mm256_cmpgt_epi32(hi, lo);
    m256i _lo = _mm256_add_epi32(lo, _mm256_and_si256(F3_m256i, cmp));

    return _mm256_sub_epi32(_lo, hi);
}

inline m256i mul_f4(m256i a, m256i b)
{
    m256i _a = _mm256_load_si256(&a);
    m256i _b = _mm256_load_si256(&b);

    m256i c = _mm256_mullo_epi32(_a, _b);

    // filter elements of both of a & b = card-1
    m256i cmp = _mm256_and_si256(
        _mm256_cmpeq_epi32(_a, F4minus1_m256i),
        _mm256_cmpeq_epi32(_b, F4minus1_m256i));

    const m256i one = _mm256_set1_epi32(1);
    c = _mm256_add_epi32(c, _mm256_and_si256(one, cmp));

    // Modulo
    return mod_after_multiply_f4(c);
}

inline m256i mul_f3(m256i a, m256i b)
{
    m256i _a = _mm256_load_si256(&a);
    m256i _b = _mm256_load_si256(&b);

    m256i c = _mm256_mullo_epi32(_a, _b);

    // filter elements of both of a & b = card-1
    m256i cmp = _mm256_and_si256(
        _mm256_cmpeq_epi32(_a, F3minus1_m256i),
        _mm256_cmpeq_epi32(_b, F3minus1_m256i));

    c = _mm256_xor_si256(c, _mm256_and_si256(F4_m256i, cmp));

    // Modulo
    return mod_after_multiply_f3(c);
}

/** Perform multiplication of two numbers a, b whose elements are of GF(card)
 *  where `card` is a prime Fermat number, i.e. card = Fx with x < 5
 *  Currently, it supports only for F3 and F4
 */
inline m256i mul(m256i a, m256i b, aint32 card = F4)
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
    const m256i coef = _mm256_set1_epi32(a);

    m256i* _src = reinterpret_cast<m256i*>(src);
    m256i* _dest = reinterpret_cast<m256i*>(dest);
    unsigned ratio = sizeof(*_src) / sizeof(*src);
    size_t _len = len / ratio;
    size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform multiplication
        _dest[i] = mul(coef, _src[i], card);
    }
    if (_last_len > 0) {
        uint64_t coef_64 = (uint64_t)a;
        for (i = _len * ratio; i < len; i++) {
            // perform multiplication
            dest[i] = (aint32)((coef_64 * src[i]) % card);
        }
    }
}

inline void
add_two_bufs(aint32* src, aint32* dest, size_t len, aint32 card = F4)
{
    m256i* _src = reinterpret_cast<m256i*>(src);
    m256i* _dest = reinterpret_cast<m256i*>(dest);
    unsigned ratio = sizeof(*_src) / sizeof(*src);
    size_t _len = len / ratio;
    size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform addition
        _dest[i] = add(_src[i], _dest[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
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
    m256i* _bufa = reinterpret_cast<m256i*>(bufa);
    m256i* _bufb = reinterpret_cast<m256i*>(bufb);
    m256i* _res = reinterpret_cast<m256i*>(res);
    unsigned ratio = sizeof(*_bufa) / sizeof(*bufa);
    size_t _len = len / ratio;
    size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform subtraction
        _res[i] = sub(_bufa[i], _bufb[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
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
    m256i* _src = reinterpret_cast<m256i*>(src);
    m256i* _dest = reinterpret_cast<m256i*>(dest);
    unsigned ratio = sizeof(*_src) / sizeof(*src);
    size_t _len = len / ratio;
    size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform multiplicaton
        _dest[i] = mul(_src[i], _dest[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
            // perform multiplicaton
            dest[i] = uint32_t((uint64_t(src[i]) * dest[i]) % card);
        }
    }
}

/* ==================== Operations for NF4 =================== */
typedef __m128i m128i;

/** Return aint128 integer from a _m128i register */
inline aint128 m256i_to_uint128(m256i v)
{
    aint128 hi, lo;
    _mm256_storeu2_m128i((m128i*)&hi, (m128i*)&lo, v);
    return lo;
}

inline __uint128_t add(__uint128_t a, __uint128_t b)
{
    m256i _a = _mm256_castsi128_si256((m128i)a);
    m256i _b = _mm256_castsi128_si256((m128i)b);
    m256i res = add(_a, _b, F4);
    return m256i_to_uint128(res);
}

inline __uint128_t sub(__uint128_t a, __uint128_t b)
{
    m256i _a = _mm256_castsi128_si256((m128i)a);
    m256i _b = _mm256_castsi128_si256((m128i)b);
    m256i res = sub(_a, _b, F4);
    return m256i_to_uint128(res);
}

inline __uint128_t mul(__uint128_t a, __uint128_t b)
{
    m256i _a = _mm256_castsi128_si256((m128i)a);
    m256i _b = _mm256_castsi128_si256((m128i)b);
    m256i res = mul(_a, _b, F4);
    return m256i_to_uint128(res);
}

/** Add buffer `y` to two halves of `x`. `x` is of length `n` */
inline void
add_buf_to_two_bufs(int n, aint128* _x, aint128* _y, aint32 card = F4)
{
    int i;
    m256i* x = reinterpret_cast<m256i*>(_x);
    m256i* y = reinterpret_cast<m256i*>(_y);

    const unsigned ratio = sizeof(*x) / sizeof(*_x);
    const int len_y = n / 2;
    const int len_y_256 = len_y / ratio;
    const int last_len_y = len_y - len_y_256 * ratio;

    aint128* x_half = _x + len_y;
    m256i* x_next = reinterpret_cast<m256i*>(x_half);

    // add y to the first half of `x`
    for (i = 0; i < len_y_256; i++) {
        x[i] = add(x[i], y[i], card);
    }

    // add y to the second half of `x`
    for (i = 0; i < len_y_256; i++) {
        x_next[i] = add(x_next[i], y[i], card);
    }

    if (last_len_y > 0) {
        // add last _y[] to x and x_next
        for (i = len_y_256 * ratio; i < len_y; i++) {
            m256i _x_p = _mm256_castsi128_si256((m128i)_x[i]);
            m256i _x_next_p = _mm256_castsi128_si256((m128i)x_half[i]);
            m256i _y_p = _mm256_castsi128_si256((m128i)_y[i]);

            _x_p = add(_x_p, _y_p, card);
            _x_next_p = add(_x_next_p, _y_p, card);
        }
    }
}

inline void hadamard_mul(int n, aint128* _x, aint128* _y)
{
    int i;
    m256i* x = reinterpret_cast<m256i*>(_x);
    m256i* y = reinterpret_cast<m256i*>(_y);

    const unsigned ratio = sizeof(*x) / sizeof(*_x);
    const int len_y = n / 2;
    const int len_y_256 = len_y / ratio;
    const int last_len_y = len_y - len_y_256 * ratio;

    aint128* x_half = _x + len_y;
    m256i* x_next = reinterpret_cast<m256i*>(x_half);

    // multiply y to the first half of `x`
    for (i = 0; i < len_y_256; i++) {
        x[i] = mul(x[i], y[i]);
    }

    // multiply y to the second half of `x`
    for (i = 0; i < len_y_256; i++) {
        x_next[i] = mul(x_next[i], y[i]);
    }

    if (last_len_y > 0) {
        // add last _y[] to x and x_next
        for (i = len_y_256 * ratio; i < len_y; i++) {
            m256i _x_p = _mm256_castsi128_si256((m128i)_x[i]);
            m256i _x_next_p = _mm256_castsi128_si256((m128i)x_half[i]);
            m256i _y_p = _mm256_castsi128_si256((m128i)_y[i]);

            _x_p = mul(_x_p, _y_p);
            _x_next_p = mul(_x_next_p, _y_p);
        }
    }
}

} // namespace u32
} // namespace simd
} // namespace nttec

#endif
