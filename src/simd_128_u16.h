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

#ifndef __NTTEC_SIMD_128_U16_H__
#define __NTTEC_SIMD_128_U16_H__

#include <x86intrin.h>

namespace nttec {

/** The namespace simd contains functions that are accelerated by using SIMD
 *  operations over 128bits
 *
 *  It supports operations on 16-bit numbers
 */
namespace simd {
namespace u16 {

typedef __m128i am128i __attribute__((aligned(ALIGN_SIZE)));

/* ==================== Essential Operations =================== */

const aint16 F3 = 257;

// Disable `cert-err58-cpp` on these: AFAIK they cannot throw.
// (probably a false positive present in Clang 5 and fixed in Clang 6).
const am128i F3_m128i = _mm_set1_epi16(257);       // NOLINT(cert-err58-cpp)
const am128i F3minus1_m128i = _mm_set1_epi16(256); // NOLINT(cert-err58-cpp)

/** Return aint128 integer from a _m128i register */
static inline aint128 m128i_to_uint128(am128i v)
{
    aint128 i;
    _mm_store_si128((am128i*)&i, v);

    return i; // NOLINT(clang-analyzer-core.uninitialized.UndefReturn)
}

/** Perform a%card where a is a addition of two numbers whose elements are
 *  symbols of GF(card) */
inline aint128 mod_after_add(am128i a, aint16 card)
{
    const am128i _card = _mm_set1_epi16(card);
    const am128i _card_minus_1 = _mm_set1_epi16(card - 1);

    am128i cmp = _mm_cmpgt_epi16(a, _card_minus_1);
    am128i b = _mm_sub_epi16(a, _mm_and_si128(_card, cmp));

    return m128i_to_uint128(b);
}

/** Perform addition of two numbers a, b whose elements are of GF(card) */
inline aint128 add(aint128 a, aint128 b, aint16 card = F3)
{
    am128i _a = _mm_load_si128((am128i*)&a);
    am128i _b = _mm_load_si128((am128i*)&b);
    am128i c = _mm_add_epi16(_a, _b);

    // Modulo
    return mod_after_add(c, card);
}

/** Perform subtraction of a by b where a, b whose elements are symbols of
 *  GF(card)
 * sub(a, b) = a - b if a >= b, or
 *             card + a - b, otherwise
 */
inline aint128 sub(aint128 a, aint128 b, aint16 card = F3)
{
    const am128i _card = _mm_set1_epi16(card);

    am128i _a = _mm_load_si128((am128i*)&a);
    am128i _b = _mm_load_si128((am128i*)&b);

    am128i cmp = _mm_cmpgt_epi16(_b, _a);
    am128i _a1 = _mm_add_epi16(_a, _mm_and_si128(_card, cmp));

    am128i c = _mm_sub_epi16(_a1, _b);

    return m128i_to_uint128(c);
}

inline aint128 mod_after_multiply(am128i a)
{
    const am128i mask = _mm_set1_epi16(F3 - 2);

    am128i lo = _mm_and_si128(a, mask);

    am128i a_shift = _mm_srli_si128(a, 1);
    am128i hi = _mm_and_si128(a_shift, mask);

    am128i cmp = _mm_cmpgt_epi16(hi, lo);
    am128i _lo = _mm_add_epi16(lo, _mm_and_si128(F3_m128i, cmp));

    am128i b = _mm_sub_epi16(_lo, hi);

    return m128i_to_uint128(b);
}

inline aint128 mul(aint128 a, aint128 b)
{
    am128i _a = _mm_load_si128((am128i*)&a);
    am128i _b = _mm_load_si128((am128i*)&b);

    am128i c = _mm_mullo_epi16(_a, _b);

    // filter elements of both of a & b = card-1
    am128i cmp = _mm_and_si128(
        _mm_cmpeq_epi16(_a, F3minus1_m128i),
        _mm_cmpeq_epi16(_b, F3minus1_m128i));

    const am128i one = _mm_set1_epi16(1);
    c = _mm_add_epi16(c, _mm_and_si128(one, cmp));

    // Modulo
    return mod_after_multiply(c);
}

/** Perform multiplication of two numbers a, b whose elements are of GF(card)
 *  where `card` is a prime Fermat number, i.e. card = Fx with x < 5
 *  Currently, it supports only for F3
 */
inline aint128 mul(aint128 a, aint128 b, aint16 card)
{
    // FIXME: generalize card
    assert(card == F3);
    return mul(a, b);
}

/** Perform a multiplication of a coefficient `a` to each element of `src` and
 *  add result to correspondent element of `dest`
 */
inline void mul_coef_to_buf(
    const aint16 a,
    aint16* src,
    aint16* dest,
    size_t len,
    aint16 card = F3)
{
    const am128i _coef = _mm_set1_epi16(a);
    aint128 coef = m128i_to_uint128(_coef);

    aint128* _src = reinterpret_cast<aint128*>(src);
    aint128* _dest = reinterpret_cast<aint128*>(dest);
    unsigned ratio = sizeof(*_src) / sizeof(*src);
    size_t _len = len / ratio;
    size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform multiplication
        _dest[i] = mul(coef, _src[i], card);
    }
    if (_last_len > 0) {
        uint32_t coef_doubled = (uint32_t)a;
        for (i = _len * ratio; i < len; i++) {
            // perform multiplication
            dest[i] = (aint16)((coef_doubled * src[i]) % card);
        }
    }
}

inline void
add_two_bufs(aint16* src, aint16* dest, size_t len, aint16 card = F3)
{
    aint128* _src = reinterpret_cast<aint128*>(src);
    aint128* _dest = reinterpret_cast<aint128*>(dest);
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
            aint16 tmp = src[i] + dest[i];
            dest[i] = (tmp >= card) ? (tmp - card) : tmp;
        }
    }
}

inline void sub_two_bufs(
    aint16* bufa,
    aint16* bufb,
    aint16* res,
    size_t len,
    aint16 card = F3)
{
    aint128* _bufa = reinterpret_cast<aint128*>(bufa);
    aint128* _bufb = reinterpret_cast<aint128*>(bufb);
    aint128* _res = reinterpret_cast<aint128*>(res);
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
mul_two_bufs(aint16* src, aint16* dest, size_t len, aint16 card = F3)
{
    aint128* _src = reinterpret_cast<aint128*>(src);
    aint128* _dest = reinterpret_cast<aint128*>(dest);
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
            dest[i] = uint32_t(src[i]) * uint32_t(dest[i]) % card;
        }
    }
}

} // namespace u16
} // namespace simd
} // namespace nttec

#endif
