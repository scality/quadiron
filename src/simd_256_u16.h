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

#ifndef __QUAD_SIMD_256_U16_H__
#define __QUAD_SIMD_256_U16_H__

#include <x86intrin.h>

namespace quadiron {
namespace simd {

#define F3_u16 _mm256_set1_epi16(257)
#define F3m1_u16 _mm256_set1_epi16(256)

/* ==================== Essential Operations =================== */
// Following functions are used for AVX2 w/ u16 only

inline m256i SET1(uint16_t val)
{
    return _mm256_set1_epi16(val);
}
inline m256i ADD16(m256i x, m256i y)
{
    return _mm256_add_epi16(x, y);
}
inline m256i SUB16(m256i x, m256i y)
{
    return _mm256_sub_epi16(x, y);
}
inline m256i MUL16(m256i x, m256i y)
{
    return _mm256_mullo_epi16(x, y);
}

inline m256i CMPEQ16(m256i x, m256i y)
{
    return _mm256_cmpeq_epi16(x, y);
}
inline m256i CMPGT16(m256i x, m256i y)
{
    return _mm256_cmpgt_epi16(x, y);
}
inline m256i MINU16(m256i x, m256i y)
{
    return _mm256_min_epu16(x, y);
}

#define MASK8_LO (_mm256_set1_epi16(0x80))
#define BLEND8(x, y, mask) (_mm256_blendv_epi8(x, y, mask))

// z = x + y mod q
// Input are loaded to registers
// Output is register
inline m256i ADD_MOD(m256i x, m256i y, uint16_t q)
{
    m256i res = ADD16(x, y);
    return MINU16(res, SUB16(res, F3_u16));
}

// z = x - y mod q => z = q + x - y mod q
// Input are loaded to registers
// Output is register
inline m256i SUB_MOD(m256i x, m256i y, uint16_t q)
{
    m256i res = SUB16(x, y);
    return MINU16(res, SUB16(ADD16(x, F3_u16), y));
}

// y = 0 - x mod q => y = q - x mod q
// Input are loaded to registers
// Output is register
inline m256i NEG_MOD(m256i x, uint16_t q)
{
    m256i res = SUB16(F3_u16, x);
    return MINU16(res, SUB16(res, F3_u16));
}

// z = x * y mod q
// Input are loaded to registers
// Output is register
// Note: we assume that at least `x` or `y` is less than `q-1` so it's
// not necessary to verify overflow on multiplying elements
inline m256i MUL_MOD(m256i x, m256i y, uint16_t q)
{
    m256i res = MUL16(x, y);
    m256i lo = BLEND8(ZERO, res, MASK8_LO);
    m256i hi = BLEND8(ZERO, SHIFTR_1(res), MASK8_LO);
    return SUB_MOD(lo, hi, q);
}

// z = x * y mod q
// Input are loaded to registers
// Output is register
inline m256i MULFULL_MOD(m256i x, m256i y, uint16_t q)
{
    m256i res = MUL16(x, y);

    // filter elements of both of a & b = card-1
    m256i cmp = AND(CMPEQ16(x, F3m1_u16), CMPEQ16(y, F3m1_u16));
    res = ADD16(res, AND(ONE, cmp));

    m256i lo = BLEND8(ZERO, res, MASK8_LO);
    m256i hi = BLEND8(ZERO, SHIFTR_1(res), MASK8_LO);
    return SUB_MOD(lo, hi, q);
}

// butterfly CT with r == 1
inline void BUTTERFLY_1(m256i* x, m256i* y, uint16_t q)
{
    m256i add = ADD_MOD(*x, *y, q);
    *y = SUB_MOD(*x, *y, q);
    *x = add;
}

// butterfly CT with r == q - 1
inline void BUTTERFLY_2(m256i* x, m256i* y, uint16_t q)
{
    m256i add = ADD_MOD(*x, *y, q);
    *x = SUB_MOD(*x, *y, q);
    *y = add;
}

// butterfly CT with 1 < r < q - 1
inline void BUTTERFLY_3(m256i c, m256i* x, m256i* y, uint16_t q)
{
    m256i z = MUL_MOD(c, *y, q);
    *y = SUB_MOD(*x, z, q);
    *x = ADD_MOD(*x, z, q);
}

// butterfly GS w/ r = q - 1
inline void BUTTERFLY_4(m256i* x, m256i* y, uint16_t q)
{
    m256i add = ADD_MOD(*x, *y, q);
    *y = SUB_MOD(*y, *x, q);
    *x = add;
}

// butterfly GS w/ 1 < r < q - 1
// x = x + y mod q
// y = z * (x - y) mod q
inline void BUTTERFLY_5(m256i c, m256i* x, m256i* y, uint16_t q)
{
    m256i sub = SUB_MOD(*x, *y, q);
    *x = ADD_MOD(*x, *y, q);
    *y = MUL_MOD(c, sub, q);
}

/**
 * Vectorized butterly CT step
 *
 * For each pair (P, Q) = (buf[i], buf[i + m]) for step = 2 * m and coef `r`
 *      P = P + r * Q
 *      Q = P - r * Q
 *
 * @param buf - working buffers
 * @param r - coefficient
 * @param start - index of buffer among `m` ones
 * @param m - current group size
 * @param len - number of vectors per buffer
 * @param card - modulo cardinal
 */
inline void butterfly_ct_step(
    vec::Buffers<uint16_t>& buf,
    uint16_t r,
    unsigned start,
    unsigned m,
    size_t len,
    uint16_t card)
{
    const unsigned step = m << 1;
    m256i c = SET1(r);

#define BUTTERFLY_CT(x, y)                                                     \
    (EITHER(                                                                   \
        r == 1,                                                                \
        BUTTERFLY_1(x, y, card),                                               \
        EITHER(                                                                \
            r < card - 1,                                                      \
            BUTTERFLY_3(c, x, y, card),                                        \
            BUTTERFLY_2(x, y, card))));

    const size_t end = len - 1;
    const unsigned bufs_nb = buf.get_n();
    // #pragma omp parallel for
    // #pragma unroll
    const std::vector<uint16_t*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        m256i x1, y1;
        m256i x2, y2;
        m256i* __restrict p = reinterpret_cast<m256i*>(mem[i]);
        m256i* __restrict q = reinterpret_cast<m256i*>(mem[i + m]);

        // #pragma omp parallel for
        size_t j = 0;
        // #pragma unroll
        for (; j < end; j += 2) {
            x1 = LOAD(p + j);
            y1 = LOAD(q + j);
            x2 = LOAD(p + j + 1);
            y2 = LOAD(q + j + 1);

            BUTTERFLY_CT(&x1, &y1);
            BUTTERFLY_CT(&x2, &y2);

            // Store back to memory
            STORE(p + j, x1);
            STORE(p + j + 1, x2);
            STORE(q + j, y1);
            STORE(q + j + 1, y2);
        }
        for (; j < len; ++j) {
            x1 = LOAD(p + j);
            y1 = LOAD(q + j);

            BUTTERFLY_CT(&x1, &y1);

            // Store back to memory
            STORE(p + j, x1);
            STORE(q + j, y1);
        }
    }
}

/**
 * Vectorized butterly CT on two-layers at a time
 *
 * For each quadruple
 * (P, Q, R, S) = (buf[i], buf[i + m], buf[i + 2 * m], buf[i + 3 * m])
 * First layer: butterfly on (P, Q) and (R, S) for step = 2 * m
 *      coef r1 = W[start * n / (2 * m)]
 *      P = P + r1 * Q
 *      Q = P - r1 * Q
 *      R = R + r1 * S
 *      S = R - r1 * S
 * Second layer: butterfly on (P, R) and (Q, S) for step = 4 * m
 *      coef r2 = W[start * n / (4 * m)]
 *      coef r3 = W[(start + m) * n / (4 * m)]
 *      P = P + r2 * R
 *      R = P - r2 * R
 *      Q = Q + r3 * S
 *      S = Q - r3 * S
 *
 * @param buf - working buffers
 * @param r1 - coefficient for the 1st layer
 * @param r2 - 1st coefficient for the 2nd layer
 * @param r3 - 2nd coefficient for the 2nd layer
 * @param start - index of buffer among `m` ones
 * @param m - current group size
 * @param len - number of vectors per buffer
 * @param card - modulo cardinal
 */
inline void butterfly_ct_two_layers_step(
    vec::Buffers<uint16_t>& buf,
    unsigned r1,
    unsigned r2,
    unsigned r3,
    unsigned start,
    unsigned m,
    size_t len,
    uint16_t card)
{
    const unsigned step = m << 2;
    m256i c1 = SET1(r1);
    m256i c2 = SET1(r2);
    m256i c3 = SET1(r3);

#define BUTTERFLY_R1(c, x, y)                                                  \
    (EITHER(                                                                   \
        r1 == 1,                                                               \
        BUTTERFLY_1(x, y, card),                                               \
        EITHER(                                                                \
            r1 < card - 1,                                                     \
            BUTTERFLY_3(c, x, y, card),                                        \
            BUTTERFLY_2(x, y, card))));
#define BUTTERFLY_R2(c, x, y)                                                  \
    (EITHER(                                                                   \
        r2 == 1,                                                               \
        BUTTERFLY_1(x, y, card),                                               \
        EITHER(                                                                \
            r2 < card - 1,                                                     \
            BUTTERFLY_3(c, x, y, card),                                        \
            BUTTERFLY_2(x, y, card))));
#define BUTTERFLY_R3(c, x, y)                                                  \
    (EITHER(                                                                   \
        r3 == 1,                                                               \
        BUTTERFLY_1(x, y, card),                                               \
        EITHER(                                                                \
            r3 < card - 1,                                                     \
            BUTTERFLY_3(c, x, y, card),                                        \
            BUTTERFLY_2(x, y, card))));

    const size_t end = len - 1;
    const unsigned bufs_nb = buf.get_n();
    // #pragma omp parallel for
    // #pragma unroll
    const std::vector<uint16_t*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        m256i x1, y1, u1, v1;
        m256i x2, y2, u2, v2;
        m256i* __restrict p = reinterpret_cast<m256i*>(mem[i]);
        m256i* __restrict q = reinterpret_cast<m256i*>(mem[i + m]);
        m256i* __restrict r = reinterpret_cast<m256i*>(mem[i + 2 * m]);
        m256i* __restrict s = reinterpret_cast<m256i*>(mem[i + 3 * m]);

        // #pragma omp parallel for
        size_t j = 0;
        // #pragma unroll
        for (; j < end; j += 2) {
            // First layer (c1, x, y) & (c1, u, v)
            x1 = LOAD(p + j);
            y1 = LOAD(q + j);
            x2 = LOAD(p + j + 1);
            y2 = LOAD(q + j + 1);

            u1 = LOAD(r + j);
            v1 = LOAD(s + j);
            u2 = LOAD(r + j + 1);
            v2 = LOAD(s + j + 1);

            BUTTERFLY_R1(c1, &x1, &y1);
            BUTTERFLY_R1(c1, &x2, &y2);

            BUTTERFLY_R1(c1, &u1, &v1);
            BUTTERFLY_R1(c1, &u2, &v2);

            // Second layer (c2, x, u) & (c3, y, v)
            BUTTERFLY_R2(c2, &x1, &u1);
            BUTTERFLY_R2(c2, &x2, &u2);

            BUTTERFLY_R3(c3, &y1, &v1);
            BUTTERFLY_R3(c3, &y2, &v2);

            // Store back to memory
            STORE(p + j, x1);
            STORE(p + j + 1, x2);
            STORE(q + j, y1);
            STORE(q + j + 1, y2);

            STORE(r + j, u1);
            STORE(r + j + 1, u2);
            STORE(s + j, v1);
            STORE(s + j + 1, v2);
        }
        for (; j < len; ++j) {
            // First layer (c1, x, y) & (c1, u, v)
            x1 = LOAD(p + j);
            y1 = LOAD(q + j);
            u1 = LOAD(r + j);
            v1 = LOAD(s + j);

            BUTTERFLY_R1(c1, &x1, &y1);
            BUTTERFLY_R1(c1, &u1, &v1);
            // Second layer (c2, x, u) & (c3, y, v)
            BUTTERFLY_R2(c2, &x1, &u1);
            BUTTERFLY_R3(c3, &y1, &v1);
            // Store back to memory
            STORE(p + j, x1);
            STORE(q + j, y1);
            STORE(r + j, u1);
            STORE(s + j, v1);
        }
    }
}

/**
 * Vectorized butterly GS step
 *
 * For each pair (P, Q) = (buf[i], buf[i + m]) for step = 2 * m and coef `r`
 *      P = P + Q
 *      Q = r * (P - Q)
 *
 * @param buf - working buffers
 * @param r - coefficient
 * @param start - index of buffer among `m` ones
 * @param m - current group size
 * @param len - number of vectors per buffer
 * @param card - modulo cardinal
 */
inline void butterfly_gs_step(
    vec::Buffers<uint16_t>& buf,
    uint16_t r,
    unsigned start,
    unsigned m,
    size_t len,
    uint16_t card)
{
    const unsigned step = m << 1;
    m256i c = SET1(r);

#define BUTTERFLY_GS(x, y)                                                     \
    (EITHER(                                                                   \
        r == 1,                                                                \
        BUTTERFLY_1(x, y, card),                                               \
        EITHER(                                                                \
            r < card - 1,                                                      \
            BUTTERFLY_5(c, x, y, card),                                        \
            BUTTERFLY_4(x, y, card))));

    const size_t end = len - 1;
    const unsigned bufs_nb = buf.get_n();
    // #pragma omp parallel for
    // #pragma unroll
    const std::vector<uint16_t*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        m256i x1, y1;
        m256i x2, y2;
        m256i* __restrict p = reinterpret_cast<m256i*>(mem[i]);
        m256i* __restrict q = reinterpret_cast<m256i*>(mem[i + m]);

        // #pragma omp parallel for
        size_t j = 0;
        // #pragma unroll
        for (; j < end; j += 2) {
            x1 = LOAD(p + j);
            y1 = LOAD(q + j);
            x2 = LOAD(p + j + 1);
            y2 = LOAD(q + j + 1);

            BUTTERFLY_GS(&x1, &y1);
            BUTTERFLY_GS(&x2, &y2);

            // Store back to memory
            STORE(p + j, x1);
            STORE(p + j + 1, x2);
            STORE(q + j, y1);
            STORE(q + j + 1, y2);
        }
        for (; j < len; ++j) {
            x1 = LOAD(p + j);
            y1 = LOAD(q + j);

            BUTTERFLY_GS(&x1, &y1);

            // Store back to memory
            STORE(p + j, x1);
            STORE(q + j, y1);
        }
    }
}

inline void add_props_16(
    Properties& props,
    m256i threshold,
    m256i mask,
    m256i symb,
    off_t offset)
{
    const m256i b = CMPEQ16(threshold, symb);
    const m256i c = AND(mask, b);
    uint32_t d = MVMSK8(c);
    const unsigned element_size = sizeof(uint16_t);
    while (d > 0) {
        unsigned byte_idx = __builtin_ctz(d);
        off_t _offset = offset + byte_idx / element_size;
        props.add(_offset, OOR_MARK);
        d ^= 1 << byte_idx;
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
    const unsigned element_size = sizeof(uint16_t);
    const unsigned vec_size = ALIGN_SIZE / element_size;
    const uint16_t max = 1 << (element_size * 8 - 1);
    const m256i _threshold = SET1(threshold);
    const m256i mask_hi = SET1(max);

    // #pragma unroll
    const std::vector<uint16_t*>& mem = output.get_mem();
    for (unsigned frag_id = 0; frag_id < code_len; ++frag_id) {
        m256i* __restrict buf = reinterpret_cast<m256i*>(mem[frag_id]);

        size_t vec_id = 0;
        size_t end = vecs_nb - 3;
        // #pragma unroll
        for (; vec_id < end; vec_id += 4) {
            m256i a1 = LOAD(buf + vec_id);
            m256i a2 = LOAD(buf + vec_id + 1);
            m256i a3 = LOAD(buf + vec_id + 2);
            m256i a4 = LOAD(buf + vec_id + 3);

            if (TESTZ(a1, _threshold) == 0) {
                const off_t curr_offset = offset + vec_id * vec_size;
                add_props_16(
                    props[frag_id], _threshold, mask_hi, a1, curr_offset);
            }
            if (TESTZ(a2, _threshold) == 0) {
                const off_t curr_offset = offset + (vec_id + 1) * vec_size;
                add_props_16(
                    props[frag_id], _threshold, mask_hi, a2, curr_offset);
            }
            if (TESTZ(a3, _threshold) == 0) {
                const off_t curr_offset = offset + (vec_id + 2) * vec_size;
                add_props_16(
                    props[frag_id], _threshold, mask_hi, a3, curr_offset);
            }
            if (TESTZ(a4, _threshold) == 0) {
                const off_t curr_offset = offset + (vec_id + 3) * vec_size;
                add_props_16(
                    props[frag_id], _threshold, mask_hi, a4, curr_offset);
            }
        }
        for (; vec_id < vecs_nb; ++vec_id) {
            m256i a = LOAD(buf + vec_id);
            uint16_t c = TESTZ(a, _threshold);
            if (c == 0) {
                const off_t curr_offset = offset + vec_id * vec_size;
                add_props_16(
                    props[frag_id], _threshold, mask_hi, a, curr_offset);
            }
        }
    }
}

/* ==================== Operations =================== */
/** Perform a multiplication of a coefficient `a` to each element of `src` and
 *  add result to correspondent element of `dest`
 *
 * @note: 1 < `a` < card - 1
 */
inline void mul_coef_to_buf(
    const uint16_t a,
    aint16* src,
    aint16* dest,
    size_t len,
    uint16_t card)
{
    const m256i coef = SET1(a);

    m256i* __restrict _src = reinterpret_cast<m256i*>(src);
    m256i* __restrict _dest = reinterpret_cast<m256i*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform multiplication
        _dest[i] = MUL_MOD(coef, _src[i], card);
    }
    if (_last_len > 0) {
        uint32_t coef_32 = (uint32_t)a;
        for (i = _len * ratio; i < len; i++) {
            // perform multiplication
            dest[i] = (aint16)((coef_32 * src[i]) % card);
        }
    }
}

inline void add_two_bufs(aint16* src, aint16* dest, size_t len, aint16 card)
{
    m256i* __restrict _src = reinterpret_cast<m256i*>(src);
    m256i* __restrict _dest = reinterpret_cast<m256i*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform addition
        _dest[i] = ADD_MOD(_src[i], _dest[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
            // perform addition
            aint16 tmp = src[i] + dest[i];
            dest[i] = (tmp >= card) ? (tmp - card) : tmp;
        }
    }
}

inline void
sub_two_bufs(aint16* bufa, aint16* bufb, aint16* res, size_t len, aint16 card)
{
    m256i* __restrict _bufa = reinterpret_cast<m256i*>(bufa);
    m256i* __restrict _bufb = reinterpret_cast<m256i*>(bufb);
    m256i* __restrict _res = reinterpret_cast<m256i*>(res);
    const unsigned ratio = sizeof(*_bufa) / sizeof(*bufa);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform subtraction
        _res[i] = SUB_MOD(_bufa[i], _bufb[i], card);
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

inline void mul_two_bufs(aint16* src, aint16* dest, size_t len, aint16 card)
{
    m256i* __restrict _src = reinterpret_cast<m256i*>(src);
    m256i* __restrict _dest = reinterpret_cast<m256i*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform multiplicaton
        _dest[i] = MULFULL_MOD(_src[i], _dest[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
            // perform multiplicaton
            dest[i] = uint16_t((uint64_t(src[i]) * dest[i]) % card);
        }
    }
}

/** Apply an element-wise negation to a buffer
 */
inline void neg(size_t len, aint16* buf, aint16 card)
{
    m256i* _buf = reinterpret_cast<m256i*>(buf);
    unsigned ratio = sizeof(*_buf) / sizeof(*buf);
    size_t _len = len / ratio;
    size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        _buf[i] = NEG_MOD(_buf[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
            if (buf[i])
                buf[i] = card - buf[i];
        }
    }
}

} // namespace simd
} // namespace quadiron

#endif
