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

#ifndef __QUAD_SIMD_128_U32_H__
#define __QUAD_SIMD_128_U32_H__

#include <x86intrin.h>

namespace quad {
namespace simd {

/* ==================== Essential Operations =================== */

/** Perform a%card where a is a addition of two numbers whose elements are
 *  symbols of GF(card) */
inline m128i mod_after_add(m128i a, aint32 card)
{
    const m128i _card = _mm_set1_epi32(card);
    const m128i _card_minus_1 = _mm_set1_epi32(card - 1);

    m128i cmp = _mm_cmpgt_epi32(a, _card_minus_1);
    m128i b = _mm_sub_epi32(a, _mm_and_si128(_card, cmp));

    return b;
}

/** Perform addition of two numbers a, b whose elements are of GF(card) */
inline m128i add(m128i a, m128i b, aint32 card)
{
    m128i _a = _mm_load_si128(&a);
    m128i _b = _mm_load_si128(&b);
    m128i c = _mm_add_epi32(_a, _b);

    // Modulo
    return mod_after_add(c, card);
}

/** Perform subtraction of a by b where a, b whose elements are symbols of
 *  GF(card)
 * sub(a, b) = a - b if a >= b, or
 *             card + a - b, otherwise
 */
inline m128i sub(m128i a, m128i b, aint32 card)
{
    const m128i _card = _mm_set1_epi32(card);

    m128i _a = _mm_load_si128(&a);
    m128i _b = _mm_load_si128(&b);

    m128i cmp = _mm_cmpgt_epi32(_b, _a);
    m128i _a1 = _mm_add_epi32(_a, _mm_and_si128(_card, cmp));

    return _mm_sub_epi32(_a1, _b);
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
inline m128i mod_after_multiply_f4(m128i a)
{
    const m128i mask = _mm_set1_epi32(F4 - 2);

    m128i lo = _mm_and_si128(a, mask);

    m128i a_shift = _mm_srli_si128(a, 2);
    m128i hi = _mm_and_si128(a_shift, mask);

    m128i cmp = _mm_cmpgt_epi32(hi, lo);
    m128i _lo = _mm_add_epi32(lo, _mm_and_si128(F4_m128i, cmp));

    return _mm_sub_epi32(_lo, hi);
}

inline m128i mod_after_multiply_f3(m128i a)
{
    const m128i mask = _mm_set1_epi32(F3 - 2);

    m128i lo = _mm_and_si128(a, mask);

    m128i a_shift = _mm_srli_si128(a, 1);
    m128i hi = _mm_and_si128(a_shift, mask);

    m128i cmp = _mm_cmpgt_epi32(hi, lo);
    m128i _lo = _mm_add_epi32(lo, _mm_and_si128(F3_m128i, cmp));

    return _mm_sub_epi32(_lo, hi);
}

inline m128i mul_f4(m128i a, m128i b)
{
    m128i _a = _mm_load_si128(&a);
    m128i _b = _mm_load_si128(&b);

    m128i c = _mm_mullo_epi32(_a, _b);

    // filter elements of both of a & b = card-1
    m128i cmp = _mm_and_si128(
        _mm_cmpeq_epi32(_a, F4minus1_m128i),
        _mm_cmpeq_epi32(_b, F4minus1_m128i));

    const m128i one = _mm_set1_epi32(1);
    c = _mm_add_epi32(c, _mm_and_si128(one, cmp));

    // Modulo
    return mod_after_multiply_f4(c);
}

inline m128i mul_f3(m128i a, m128i b)
{
    m128i _a = _mm_load_si128(&a);
    m128i _b = _mm_load_si128(&b);

    m128i c = _mm_mullo_epi32(_a, _b);

    // filter elements of both of a & b = card-1
    m128i cmp = _mm_and_si128(
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
inline m128i mul(m128i a, m128i b, aint32 card)
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
    const m128i coef = _mm_set1_epi32(a);

    m128i* _src = reinterpret_cast<m128i*>(src);
    m128i* _dest = reinterpret_cast<m128i*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

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
    m128i* _src = reinterpret_cast<m128i*>(src);
    m128i* _dest = reinterpret_cast<m128i*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

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
    m128i* _bufa = reinterpret_cast<m128i*>(bufa);
    m128i* _bufb = reinterpret_cast<m128i*>(bufb);
    m128i* _res = reinterpret_cast<m128i*>(res);
    const unsigned ratio = sizeof(*_bufa) / sizeof(*bufa);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

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
    m128i* _src = reinterpret_cast<m128i*>(src);
    m128i* _dest = reinterpret_cast<m128i*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

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

/** Return aint128 integer from a _m128i register */
static inline aint128 m128i_to_uint128(m128i v)
{
    aint128 i;
    _mm_store_si128((m128i*)&i, v);

    return i; // NOLINT(clang-analyzer-core.uninitialized.UndefReturn)
}

inline __uint128_t add(__uint128_t a, __uint128_t b)
{
    m128i res = add((m128i)a, (m128i)b, F4);
    return m128i_to_uint128(res);
}

inline __uint128_t sub(__uint128_t a, __uint128_t b)
{
    m128i res = sub((m128i)a, (m128i)b, F4);
    return m128i_to_uint128(res);
}

inline __uint128_t mul(__uint128_t a, __uint128_t b)
{
    m128i res = mul((m128i)a, (m128i)b, F4);
    return m128i_to_uint128(res);
}

/** Add buffer `y` to two halves of `x`. `x` is of length `n` */
inline void
add_buf_to_two_bufs(int n, aint128* _x, aint128* _y, aint32 card = F4)
{
    int i;
    const int half = n / 2;
    m128i* x = reinterpret_cast<m128i*>(_x);
    m128i* y = reinterpret_cast<m128i*>(_y);
    m128i* x_next = x + half;

    // add y to the first half of `x`
    for (i = 0; i < half; i++) {
        x[i] = add(x[i], y[i], card);
    }

    // add y to the second half of `x`
    for (i = 0; i < half; i++) {
        x_next[i] = add(x_next[i], y[i], card);
    }
}

inline void hadamard_mul(int n, aint128* _x, aint128* _y)
{
    int i;
    m128i* x = reinterpret_cast<m128i*>(_x);
    m128i* y = reinterpret_cast<m128i*>(_y);

    // multiply y to `x`
    for (i = 0; i < n; i++) {
        x[i] = mul(x[i], y[i], F4);
    }
}

inline void hadamard_mul_doubled(int n, aint128* _x, aint128* _y)
{
    int i;
    const int half = n / 2;
    m128i* x = reinterpret_cast<m128i*>(_x);
    m128i* y = reinterpret_cast<m128i*>(_y);
    m128i* x_next = x + half;

    // multiply y to the first half of `x`
    for (i = 0; i < half; i++) {
        x[i] = mul(x[i], y[i], F4);
    }

    // multiply y to the second half of `x`
    for (i = 0; i < half; i++) {
        x_next[i] = mul(x_next[i], y[i], F4);
    }
}

} // namespace simd
} // namespace quad

#endif
