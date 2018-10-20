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

#ifndef __QUAD_SIMD_128_U16_H__
#define __QUAD_SIMD_128_U16_H__

#include <x86intrin.h>

#include "property.h"
#include "simd/simd.h"

namespace quadiron {
namespace simd {

/* ==================== Essential Operations =================== */

/** Perform a%card where a is a addition of two numbers whose elements are
 *  symbols of GF(card) */
inline m128i mod_after_add(m128i a, aint16 card)
{
    const m128i _card = _mm_set1_epi16(card);
    const m128i _card_minus_1 = _mm_set1_epi16(card - 1);

    m128i cmp = _mm_cmpgt_epi16(a, _card_minus_1);
    m128i b = _mm_sub_epi16(a, _mm_and_si128(_card, cmp));

    return b;
}

/** Perform addition of two numbers a, b whose elements are of GF(card) */
inline m128i add(m128i a, m128i b, aint16 card)
{
    m128i _a = _mm_load_si128(&a);
    m128i _b = _mm_load_si128(&b);
    m128i c = _mm_add_epi16(_a, _b);

    // Modulo
    return mod_after_add(c, card);
}

/** Perform subtraction of a by b where a, b whose elements are symbols of
 *  GF(card)
 * sub(a, b) = a - b if a >= b, or
 *             card + a - b, otherwise
 */
inline m128i sub(m128i a, m128i b, aint16 card)
{
    const m128i _card = _mm_set1_epi16(card);

    m128i _a = _mm_load_si128(&a);
    m128i _b = _mm_load_si128(&b);

    m128i cmp = _mm_cmpgt_epi16(_b, _a);
    m128i _a1 = _mm_add_epi16(_a, _mm_and_si128(_card, cmp));

    return _mm_sub_epi16(_a1, _b);
}

/** Negate `a`
 * @return 0 if (a == 0), else card - a
 */
inline m128i neg(m128i a, aint16 card = F3)
{
    const m128i _card = _mm_set1_epi16(card);
    m128i _a = _mm_load_si128(&a);
    m128i _b = _mm_setzero_si128();

    m128i cmp = _mm_cmpgt_epi16(_a, _b);

    return _mm_sub_epi16(_mm_and_si128(cmp, _card), _a);
}

inline m128i mod_after_multiply(m128i a)
{
    const m128i mask = _mm_set1_epi16(F3 - 2);

    m128i lo = _mm_and_si128(a, mask);

    m128i a_shift = _mm_srli_si128(a, 1);
    m128i hi = _mm_and_si128(a_shift, mask);

    m128i cmp = _mm_cmpgt_epi16(hi, lo);
    m128i _lo = _mm_add_epi16(lo, _mm_and_si128(F3_m128i_u16, cmp));

    return _mm_sub_epi16(_lo, hi);
}

inline m128i mul(m128i a, m128i b)
{
    m128i _a = _mm_load_si128(&a);
    m128i _b = _mm_load_si128(&b);

    m128i c = _mm_mullo_epi16(_a, _b);

    // filter elements of both of a & b = card-1
    m128i cmp = _mm_and_si128(
        _mm_cmpeq_epi16(_a, F3minus1_m128i_u16),
        _mm_cmpeq_epi16(_b, F3minus1_m128i_u16));

    const m128i one = _mm_set1_epi16(1);
    c = _mm_add_epi16(c, _mm_and_si128(one, cmp));

    // Modulo
    return mod_after_multiply(c);
}

/** Perform multiplication of two numbers a, b whose elements are of GF(card)
 *  where `card` is a prime Fermat number, i.e. card = Fx with x < 5
 *  Currently, it supports only for F3
 */
inline m128i mul(m128i a, m128i b, aint16 card)
{
    // FIXME: generalize card
    assert(card == F3);
    return mul(a, b);
}

/** Apply an element-wise negation to a buffer
 */
inline void neg(size_t len, aint16* buf, aint16 card = F3)
{
    m128i* _buf = reinterpret_cast<m128i*>(buf);
    unsigned ratio = sizeof(*_buf) / sizeof(*buf);
    size_t _len = len / ratio;
    size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        _buf[i] = neg(_buf[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
            if (buf[i])
                buf[i] = card - buf[i];
        }
    }
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
    const m128i coef = _mm_set1_epi16(a);

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
        uint32_t coef_doubled = a;
        for (i = _len * ratio; i < len; i++) {
            // perform multiplication
            dest[i] = narrow_cast<aint16>((coef_doubled * src[i]) % card);
        }
    }
}

inline void
add_two_bufs(aint16* src, aint16* dest, size_t len, aint16 card = F3)
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
mul_two_bufs(aint16* src, aint16* dest, size_t len, aint16 card = F3)
{
    m128i* _src = reinterpret_cast<m128i*>(src);
    m128i* _dest = reinterpret_cast<m128i*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform multiplicaton
        _dest[i] = mul(_src[i], _dest[i], F3);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
            // perform multiplicaton
            // dest[i] = uint32_t(src[i]) * uint32_t(dest[i]) % card;
            dest[i] = uint16_t((uint32_t(src[i]) * dest[i]) % card);
        }
    }
}

/*
 * buf1[i] = buf1[i] + coef * buf2[i]
 * buf2[i] = buf1[i] - coef * buf2[i]
 */
inline void butterfly_ct(
    uint16_t coef,
    aint16* buf1,
    aint16* buf2,
    size_t len,
    uint32_t card = F3)
{
    const m128i _coef = _mm_set1_epi16(coef);
    m128i* _buf1 = reinterpret_cast<m128i*>(buf1);
    m128i* _buf2 = reinterpret_cast<m128i*>(buf2);

    for (size_t i = 0; i < len; ++i) {
        m128i a = mul(_coef, _buf2[i], card);
        _buf2[i] = sub(_buf1[i], a, card);
        _buf1[i] = add(_buf1[i], a, card);
    }
}

/*
 * buf1[i] = buf1[i] + buf2[i]
 * buf2[i] = coef * (buf1[i] - buf2[i])
 */
inline void butterfly_gs(
    uint16_t coef,
    aint16* buf1,
    aint16* buf2,
    size_t len,
    uint16_t card = F3)
{
    const m128i _coef = _mm_set1_epi16(coef);
    m128i* _buf1 = reinterpret_cast<m128i*>(buf1);
    m128i* _buf2 = reinterpret_cast<m128i*>(buf2);

    for (size_t i = 0; i < len; ++i) {
        m128i a = _buf1[i];
        m128i b = _buf2[i];
        m128i c = sub(a, b, card);
        _buf1[i] = add(a, b, card);
        _buf2[i] = mul(_coef, c, card);
    }
}

inline void encode_post_process(
    vec::Buffers<uint16_t>& output,
    std::vector<Properties>& props,
    off_t offset,
    unsigned code_len,
    uint16_t threshold,
    size_t vecs_nb)
{
    const unsigned vec_size = simd::countof<uint16_t>();

    const m128i _threshold = _mm_set1_epi16(threshold);
    uint16_t max = 1 << (sizeof(uint16_t) * 8 - 1);
    const m128i mask_hi = _mm_set1_epi16(max);
    const unsigned element_size = sizeof(uint16_t);

    for (unsigned frag_id = 0; frag_id < code_len; ++frag_id) {
        uint16_t* chunk = output.get(frag_id);
        m128i* buf = reinterpret_cast<m128i*>(chunk);
        for (unsigned vec_id = 0; vec_id < vecs_nb; ++vec_id) {
            const m128i a = _mm_load_si128(&(buf[vec_id]));
            const m128i b = _mm_cmpeq_epi16(_threshold, a);
            const m128i c = _mm_and_si128(mask_hi, b);
            uint16_t d = _mm_movemask_epi8(c);

            while (d > 0) {
                unsigned byte_idx = __builtin_ctz(d);
                unsigned element_idx = byte_idx / element_size;
                off_t _offset = offset + vec_id * vec_size + element_idx;
                props[frag_id].add(_offset, 1);
                d ^= 1 << byte_idx;
            }
        }
    }
}

} // namespace simd
} // namespace quadiron

#endif
