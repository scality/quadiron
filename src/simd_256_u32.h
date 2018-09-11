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

#ifndef __QUAD_SIMD_256_U32_H__
#define __QUAD_SIMD_256_U32_H__

#include <x86intrin.h>

namespace quadiron {
namespace simd {

/* ==================== Essential Operations =================== */

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
inline m256i add(m256i a, m256i b, aint32 card)
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
inline m256i sub(m256i a, m256i b, aint32 card)
{
    const m256i _card = _mm256_set1_epi32(card);

    m256i _a = _mm256_load_si256(&a);
    m256i _b = _mm256_load_si256(&b);

    m256i cmp = _mm256_cmpgt_epi32(_b, _a);
    m256i _a1 = _mm256_add_epi32(_a, _mm256_and_si256(_card, cmp));

    return _mm256_sub_epi32(_a1, _b);
}

/** Negate `a`
 * @return 0 if (a == 0), else card - a
 */
inline m256i neg(m256i a, aint32 card = F4)
{
    const m256i _card = _mm256_set1_epi32(card);
    m256i _a = _mm256_load_si256(&a);
    m256i _b = _mm256_setzero_si256();

    m256i cmp = _mm256_cmpgt_epi32(_a, _b);

    return _mm256_sub_epi32(_mm256_and_si256(cmp, _card), _a);
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

inline m256i mul_f4_simple(m256i a, m256i b)
{
    m256i _a = _mm256_load_si256(&a);
    m256i _b = _mm256_load_si256(&b);

    m256i c = _mm256_mullo_epi32(_a, _b);

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

inline m256i mul_f3_simple(m256i a, m256i b)
{
    m256i _a = _mm256_load_si256(&a);
    m256i _b = _mm256_load_si256(&b);

    m256i c = _mm256_mullo_epi32(_a, _b);

    // Modulo
    return mod_after_multiply_f3(c);
}

/** Perform multiplication of two numbers a, b whose elements are of GF(card)
 *  where `card` is a prime Fermat number, i.e. card = Fx with x < 5
 *  Currently, it supports only for F3 and F4
 */
inline m256i mul(m256i a, m256i b, aint32 card)
{
    assert(card == F4 || card == F3);
    if (card == F4)
        return mul_f4(a, b);
    return mul_f3(a, b);
}

inline m256i mul_simple(m256i a, m256i b, aint32 card)
{
    assert(card == F4 || card == F3);
    if (card == F4)
        return mul_f4_simple(a, b);
    return mul_f3_simple(a, b);
}

/** Apply an element-wise negation to a buffer
 */
inline void neg(size_t len, aint32* buf, aint32 card = F4)
{
    m256i* _buf = reinterpret_cast<m256i*>(buf);
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
    const aint32 a,
    aint32* src,
    aint32* dest,
    size_t len,
    aint32 card = F4)
{
    const m256i coef = _mm256_set1_epi32(a);

    m256i* _src = reinterpret_cast<m256i*>(src);
    m256i* _dest = reinterpret_cast<m256i*>(dest);
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
    m256i* _src = reinterpret_cast<m256i*>(src);
    m256i* _dest = reinterpret_cast<m256i*>(dest);
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
    m256i* _bufa = reinterpret_cast<m256i*>(bufa);
    m256i* _bufb = reinterpret_cast<m256i*>(bufb);
    m256i* _res = reinterpret_cast<m256i*>(res);
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
    m256i* _src = reinterpret_cast<m256i*>(src);
    m256i* _dest = reinterpret_cast<m256i*>(dest);
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

/*
 * buf1[i] = buf1[i] + coef * buf2[i]
 * buf2[i] = buf1[i] - coef * buf2[i]
 */
inline void butterfly_ct(
    uint32_t coef,
    aint32* buf1,
    aint32* buf2,
    size_t len,
    uint32_t card = F4)
{
    const m256i _coef = _mm256_set1_epi32(coef);
    m256i* _buf1 = reinterpret_cast<m256i*>(buf1);
    m256i* _buf2 = reinterpret_cast<m256i*>(buf2);

    for (size_t i = 0; i < len; ++i) {
        m256i a = mul(_coef, _buf2[i], card);
        _buf2[i] = sub(_buf1[i], a, card);
        _buf1[i] = add(_buf1[i], a, card);
    }
}

// outputA = inputA + inputB
// outputB = inputA - inputB
inline void butterfly_step(
    m256i* inputA,
    m256i* inputB,
    m256i* outputA,
    m256i* outputB,
    uint32_t _card)
{
    const m256i card = (_card == F3) ? F3_m256i : F4_m256i;
    const m256i card_1 = (_card == F3) ? F3minus1_m256i : F4minus1_m256i;

    // --------------------------------------
    // outputB = inputA - inputB
    // --------------------------------------
    m256i a = _mm256_load_si256(inputA);
    m256i b = _mm256_load_si256(inputB);
    m256i cmp_1 = _mm256_cmpgt_epi32(b, a);
    m256i res_1 = _mm256_add_epi32(a, _mm256_and_si256(card, cmp_1));

    _mm256_store_si256(outputB, _mm256_sub_epi32(res_1, b));

    // --------------------------------------
    // outputA = symbA + symbB
    // --------------------------------------
    m256i res_2 = _mm256_add_epi32(a, b);
    // modulo
    m256i cmp_2 = _mm256_cmpgt_epi32(res_2, card_1);
    m256i c = _mm256_sub_epi32(res_2, _mm256_and_si256(card, cmp_2));

    _mm256_store_si256(outputA, c);
}

// for each pair (P, Q) = (buf[i], buf[i + m]):
// P = P + Q
// Q = P - Q
inline void butterfly_ct_1(
    vec::Buffers<uint32_t>& buf,
    unsigned start,
    unsigned m,
    unsigned step,
    size_t len,
    uint32_t card = F4)
{
    for (unsigned i = start; i < buf.get_n(); i += step) {
        uint32_t* a = buf.get(i);
        uint32_t* b = buf.get(i + m);
        m256i* _a = reinterpret_cast<m256i*>(a);
        m256i* _b = reinterpret_cast<m256i*>(b);
        // perform butterfly operation for Cooley-Tukey FFT algorithm
        for (size_t j = 0; j < len; ++j) {
            butterfly_step(&(_a[j]), &(_b[j]), &(_a[j]), &(_b[j]), card);
        }
    }
}

// for each pair (P, Q) = (buf[i], buf[i + m]):
// P = P - Q
// Q = P + Q
inline void butterfly_ct_2(
    vec::Buffers<uint32_t>& buf,
    unsigned start,
    unsigned m,
    unsigned step,
    size_t len,
    uint32_t card = F4)
{
    for (unsigned i = start; i < buf.get_n(); i += step) {
        uint32_t* a = buf.get(i);
        uint32_t* b = buf.get(i + m);
        m256i* _a = reinterpret_cast<m256i*>(a);
        m256i* _b = reinterpret_cast<m256i*>(b);
        // perform butterfly operation for Cooley-Tukey FFT algorithm
        for (size_t j = 0; j < len; ++j) {
            butterfly_step(&(_a[j]), &(_b[j]), &(_b[j]), &(_a[j]), card);
        }
    }
}

// output = coef * input
inline void
butterfly_mul(m256i* coef, m256i* input, m256i* output, uint32_t _card)
{
    const m256i card = (_card == F3) ? F3_m256i : F4_m256i;
    const m256i card_2 = (_card == F3) ? F3minus2_m256i : F4minus2_m256i;

    // --------------------------------------
    // compute coef * symbB
    // --------------------------------------
    m256i _coef = _mm256_load_si256(coef);
    m256i b = _mm256_load_si256(input);
    m256i res = _mm256_mullo_epi32(_coef, b);
    // modulo
    m256i lo = _mm256_and_si256(res, card_2);
    m256i res_shift =
        (_card == F3) ? _mm256_srli_si256(res, 1) : _mm256_srli_si256(res, 2);
    m256i hi = _mm256_and_si256(res_shift, card_2);

    m256i cmp_1 = _mm256_cmpgt_epi32(hi, lo);
    m256i _lo = _mm256_add_epi32(lo, _mm256_and_si256(card, cmp_1));

    m256i res_2 = _mm256_sub_epi32(_lo, hi);

    _mm256_store_si256(output, res_2);
}

// symbA = symbA + coef * symbB
// symbB = symbA - coef * symbB
inline void
butterfly_ct_3_step(m256i* coef, m256i* symbA, m256i* symbB, uint32_t card)
{
    // --------------------------------------
    // compute coef * symbB
    // --------------------------------------
    m256i coef_x_symbB;
    butterfly_mul(coef, symbB, &coef_x_symbB, card);

    // --------------------------------------
    // symbA = symbA + coef_x_symbB
    // symbB = symbA - coef_x_symbB
    // --------------------------------------
    butterfly_step(symbA, &coef_x_symbB, symbA, symbB, card);
}

// for each pair (P, Q) = (buf[i], buf[i + m]):
// P = P + c * Q
// Q = P - c * Q
inline void butterfly_ct_3(
    uint32_t coef,
    vec::Buffers<uint32_t>& buf,
    unsigned start,
    unsigned m,
    unsigned step,
    size_t len,
    uint32_t card = F4)
{
    m256i _coef = _mm256_set1_epi32(coef);
    for (unsigned i = start; i < buf.get_n(); i += step) {
        uint32_t* a = buf.get(i);
        uint32_t* b = buf.get(i + m);
        m256i* _a = reinterpret_cast<m256i*>(a);
        m256i* _b = reinterpret_cast<m256i*>(b);
        // perform butterfly operation for Cooley-Tukey FFT algorithm
        for (size_t j = 0; j < len; ++j) {
            butterfly_ct_3_step(&_coef, &(_a[j]), &(_b[j]), card);
        }
    }
}

/*
 * buf1[i] = buf1[i] + buf2[i]
 * buf2[i] = coef * (buf1[i] - buf2[i])
 */
inline void butterfly_gs(
    uint32_t coef,
    aint32* buf1,
    aint32* buf2,
    size_t len,
    uint32_t card = F4)
{
    const m256i _coef = _mm256_set1_epi32(coef);
    m256i* _buf1 = reinterpret_cast<m256i*>(buf1);
    m256i* _buf2 = reinterpret_cast<m256i*>(buf2);

    for (size_t i = 0; i < len; ++i) {
        m256i a = add(_buf1[i], _buf2[i], card);
        _buf2[i] = mul(_coef, sub(_buf1[i], _buf2[i], card), card);
        _buf1[i] = a;
    }
}

inline void encode_post_process(
    vec::Buffers<uint32_t>& output,
    std::vector<Properties>& props,
    off_t offset,
    unsigned code_len,
    uint32_t threshold,
    size_t vecs_nb)
{
    const unsigned vec_size = ALIGN_SIZE / sizeof(uint32_t);

    const m256i _threshold = _mm256_set1_epi32(threshold);
    const uint32_t max = 1 << (sizeof(uint32_t) * 8 - 1);
    const m256i mask_hi = _mm256_set1_epi32(max);
    const unsigned element_size = sizeof(uint32_t);

    for (unsigned frag_id = 0; frag_id < code_len; ++frag_id) {
        uint32_t* chunk = output.get(frag_id);
        m256i* buf = reinterpret_cast<m256i*>(chunk);
        for (unsigned vec_id = 0; vec_id < vecs_nb; ++vec_id) {
            const m256i a = _mm256_load_si256(&(buf[vec_id]));
            const m256i b = _mm256_cmpeq_epi32(_threshold, a);
            const m256i c = _mm256_and_si256(mask_hi, b);
            uint32_t d = _mm256_movemask_epi8(c);

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

/* ==================== Operations for NF4 =================== */
typedef __m128i m128i;

/** Return aint128 integer from a _m128i register */
inline aint128 m256i_to_uint128(m256i v)
{
    aint128 hi, lo;
    _mm256_storeu2_m128i((m128i*)&hi, (m128i*)&lo, v);
    return lo; // NOLINT(clang-analyzer-core.uninitialized.UndefReturn)
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
    const int len_256 = n / ratio;
    const int last_len = n - len_256 * ratio;

    // multiply y to the first half of `x`
    for (i = 0; i < len_256; i++) {
        x[i] = mul(x[i], y[i], F4);
    }

    if (last_len > 0) {
        // add last _y[] to x
        for (i = len_256 * ratio; i < n; i++) {
            m256i _x_p = _mm256_castsi128_si256((m128i)_x[i]);
            m256i _y_p = _mm256_castsi128_si256((m128i)_y[i]);

            _x_p = mul(_x_p, _y_p, F4);
        }
    }
}

inline void hadamard_mul_doubled(int n, aint128* _x, aint128* _y)
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
        x[i] = mul(x[i], y[i], F4);
    }

    // multiply y to the second half of `x`
    for (i = 0; i < len_y_256; i++) {
        x_next[i] = mul(x_next[i], y[i], F4);
    }

    if (last_len_y > 0) {
        // add last _y[] to x and x_next
        for (i = len_y_256 * ratio; i < len_y; i++) {
            m256i _x_p = _mm256_castsi128_si256((m128i)_x[i]);
            m256i _x_next_p = _mm256_castsi128_si256((m128i)x_half[i]);
            m256i _y_p = _mm256_castsi128_si256((m128i)_y[i]);

            _x_p = mul(_x_p, _y_p, F4);
            _x_next_p = mul(_x_next_p, _y_p, F4);
        }
    }
}

} // namespace simd
} // namespace quadiron

#endif
