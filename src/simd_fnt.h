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

#ifndef __QUAD_SIMD_FNT_H__
#define __QUAD_SIMD_FNT_H__

#include <x86intrin.h>

namespace quadiron {
namespace simd {

/* ================= Vectorized Operations ================= */

/**
 * Butterfly Cooley-Tukey operation
 *
 * x <- x + r * y
 * y <- x - r * y
 *
 * @param rp1 coefficient `r` plus one
 * @param c a register stores coefficient `r`
 * @param x working register
 * @param y working register
 * @param q modular
 */
template <typename T>
inline void BUTTERFLY_CT(T rp1, VecType c, VecType* x, VecType* y, T q)
{
    VecType z = (rp1 == 2) ? *y : MUL_MOD(c, *y, q);
    if (rp1 < q) {
        *y = SUB_MOD(*x, z, q);
        *x = ADD_MOD(*x, z, q);
    } else { // i.e. r == q - 1
        *y = ADD_MOD(*x, z, q);
        *x = SUB_MOD(*x, z, q);
    }
}

/**
 * Butterfly Genteleman-Sande operation
 *
 * x <- x + y
 * y <- r * (x - y)
 *
 * @param rp1 coefficient `r` plus one
 * @param c a register stores coefficient `r`
 * @param x working register
 * @param y working register
 * @param q modular
 */
template <typename T>
inline void BUTTERFLY_GS(T rp1, VecType c, VecType* x, VecType* y, T q)
{
    VecType add = ADD_MOD(*x, *y, q);
    if (rp1 == 2) {
        *y = SUB_MOD(*x, *y, q);
    } else if (rp1 < q) {
        VecType sub = SUB_MOD(*x, *y, q);
        *y = MUL_MOD(c, sub, q);
    } else { // i.e. r == q - 1
        *y = SUB_MOD(*y, *x, q);
    }
    *x = add;
}

/**
 * Butterfly Genteleman-Sande simple operation where y = 0
 *
 * x <- x, i.e. no operation
 * y <- r * x
 *
 * @param rp1 coefficient `r` plus one
 * @param c a register stores coefficient `r`
 * @param x working register
 * @param q modular
 * @return r * x
 */
template <typename T>
inline VecType BUTTERFLY_GS_SIMPLE(T rp1, VecType c, VecType x, T q)
{
    if (rp1 == 2) {
        return x;
    } else if (rp1 < q) {
        return MUL_MOD(c, x, q);
    } else {
        return NEG_MOD(x, q);
    }
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
 * @param step - next loop
 * @param len - number of vectors per buffer
 * @param card - modulo cardinal
 */
template <typename T>
inline void butterfly_ct_step(
    vec::Buffers<T>& buf,
    T r,
    unsigned start,
    unsigned m,
    unsigned step,
    size_t len,
    T card)
{
    if (len == 0) {
        return;
    }
    const T rp1 = r + 1;
    VecType c = SET1(r);

    const size_t end = (len > 1) ? len - 1 : 0;
    const unsigned bufs_nb = buf.get_n();
    // #pragma omp parallel for
    // #pragma unroll
    const std::vector<T*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        VecType x1, y1;
        VecType x2, y2;
        VecType* __restrict p = reinterpret_cast<VecType*>(mem[i]);
        VecType* __restrict q = reinterpret_cast<VecType*>(mem[i + m]);

        // #pragma omp parallel for
        size_t j = 0;
        // #pragma unroll
        for (; j < end; j += 2) {
            x1 = LOAD(p + j);
            y1 = LOAD(q + j);

            BUTTERFLY_CT(rp1, c, &x1, &y1, card);

            x2 = LOAD(p + j + 1);
            y2 = LOAD(q + j + 1);

            BUTTERFLY_CT(rp1, c, &x2, &y2, card);

            // Store back to memory
            STORE(p + j, x1);
            STORE(p + j + 1, x2);
            STORE(q + j, y1);
            STORE(q + j + 1, y2);
        }
        for (; j < len; ++j) {
            x1 = LOAD(p + j);
            y1 = LOAD(q + j);

            BUTTERFLY_CT(rp1, c, &x1, &y1, card);

            // Store back to memory
            STORE(p + j, x1);
            STORE(q + j, y1);
        }
    }
}

template <typename T>
inline static void do_butterfly_ct_2_layers(
    const std::vector<T*>& mem,
    T r1,
    T r2,
    T r3,
    unsigned start,
    unsigned m,
    size_t len,
    T card)
{
    const T r1p1 = r1 + 1;
    const T r2p1 = r2 + 1;
    const T r3p1 = r3 + 1;

    VecType c1 = SET1(r1);
    VecType c2 = SET1(r2);
    VecType c3 = SET1(r3);

    VecType* __restrict p = reinterpret_cast<VecType*>(mem[start]);
    VecType* __restrict q = reinterpret_cast<VecType*>(mem[start + m]);
    VecType* __restrict r = reinterpret_cast<VecType*>(mem[start + 2 * m]);
    VecType* __restrict s = reinterpret_cast<VecType*>(mem[start + 3 * m]);

    // #pragma omp parallel for
    size_t j = 0;
    const size_t end = (len > 1) ? len - 1 : 0;
    // #pragma unroll
    while (j < end) {
        // First layer (c1, x, y) & (c1, u, v)
        VecType x1 = LOAD(p);
        VecType x2 = LOAD(p + 1);
        VecType y1 = LOAD(q);
        VecType y2 = LOAD(q + 1);

        BUTTERFLY_CT(r1p1, c1, &x1, &y1, card);
        BUTTERFLY_CT(r1p1, c1, &x2, &y2, card);

        VecType u1 = LOAD(r);
        VecType u2 = LOAD(r + 1);
        VecType v1 = LOAD(s);
        VecType v2 = LOAD(s + 1);

        BUTTERFLY_CT(r1p1, c1, &u1, &v1, card);
        BUTTERFLY_CT(r1p1, c1, &u2, &v2, card);

        // Second layer (c2, x, u) & (c3, y, v)
        BUTTERFLY_CT(r2p1, c2, &x1, &u1, card);
        BUTTERFLY_CT(r2p1, c2, &x2, &u2, card);

        BUTTERFLY_CT(r3p1, c3, &y1, &v1, card);
        BUTTERFLY_CT(r3p1, c3, &y2, &v2, card);

        // Store back to memory
        STORE(p, x1);
        STORE(p + 1, x2);
        STORE(q, y1);
        STORE(q + 1, y2);

        STORE(r, u1);
        STORE(r + 1, u2);
        STORE(s, v1);
        STORE(s + 1, v2);
        p = p + 2;
        q = q + 2;
        r = r + 2;
        s = s + 2;
        j = j + 2;
    };

    for (; j < len; ++j) {
        // First layer (c1, x, y) & (c1, u, v)
        VecType x1 = LOAD(p + j);
        VecType y1 = LOAD(q + j);
        VecType u1 = LOAD(r + j);
        VecType v1 = LOAD(s + j);

        // BUTTERFLY_3_test(c1, &x1, &y1, &u1, &v1, card);
        BUTTERFLY_CT(r1p1, c1, &x1, &y1, card);
        BUTTERFLY_CT(r1p1, c1, &u1, &v1, card);
        BUTTERFLY_CT(r2p1, c2, &x1, &u1, card);
        BUTTERFLY_CT(r3p1, c3, &y1, &v1, card);

        // Store back to memory
        STORE(p + j, x1);
        STORE(q + j, y1);
        STORE(r + j, u1);
        STORE(s + j, v1);
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
template <typename T>
inline void butterfly_ct_two_layers_step(
    vec::Buffers<T>& buf,
    T r1,
    T r2,
    T r3,
    unsigned start,
    unsigned m,
    size_t len,
    T card)
{
    if (len == 0) {
        return;
    }
    const unsigned step = m << 2;
    const unsigned bufs_nb = buf.get_n();

    // #pragma omp parallel for
    // #pragma unroll
    const std::vector<T*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        do_butterfly_ct_2_layers(mem, r1, r2, r3, i, m, len, card);
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
template <typename T>
inline void butterfly_gs_step(
    vec::Buffers<T>& buf,
    T r,
    unsigned start,
    unsigned m,
    size_t len,
    T card)
{
    if (len == 0) {
        return;
    }
    const unsigned step = m << 1;
    const T rp1 = r + 1;
    VecType c = SET1(r);

    const size_t end = (len > 3) ? len - 3 : 0;
    const unsigned bufs_nb = buf.get_n();
    // #pragma omp parallel for
    // #pragma unroll
    const std::vector<T*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        VecType x1, x2, x3, x4;
        VecType y1, y2, y3, y4;
        VecType* __restrict p = reinterpret_cast<VecType*>(mem[i]);
        VecType* __restrict q = reinterpret_cast<VecType*>(mem[i + m]);

        // #pragma omp parallel for
        size_t j = 0;
        // #pragma unroll
        for (; j < end; j += 4) {
            x1 = LOAD(p + j);
            x2 = LOAD(p + j + 1);
            x3 = LOAD(p + j + 2);
            x4 = LOAD(p + j + 3);
            y1 = LOAD(q + j);
            y2 = LOAD(q + j + 1);
            y3 = LOAD(q + j + 2);
            y4 = LOAD(q + j + 3);

            BUTTERFLY_GS(rp1, c, &x1, &y1, card);
            BUTTERFLY_GS(rp1, c, &x2, &y2, card);
            BUTTERFLY_GS(rp1, c, &x3, &y3, card);
            BUTTERFLY_GS(rp1, c, &x4, &y4, card);

            // Store back to memory
            STORE(p + j, x1);
            STORE(p + j + 1, x2);
            STORE(p + j + 2, x3);
            STORE(p + j + 3, x4);
            STORE(q + j, y1);
            STORE(q + j + 1, y2);
            STORE(q + j + 2, y3);
            STORE(q + j + 3, y4);
        }
        for (; j < len; ++j) {
            x1 = LOAD(p + j);
            y1 = LOAD(q + j);

            BUTTERFLY_GS(rp1, c, &x1, &y1, card);

            // Store back to memory
            STORE(p + j, x1);
            STORE(q + j, y1);
        }
    }
}

/**
 * Vectorized butterly GS step
 *
 * For each pair (P, Q) = (buf[i], buf[i + m]) for step = 2 * m and coef `r`
 *      Q = r * P
 *
 * @param buf - working buffers
 * @param r - coefficient
 * @param start - index of buffer among `m` ones
 * @param m - current group size
 * @param len - number of vectors per buffer
 * @param card - modulo cardinal
 */
template <typename T>
inline void butterfly_gs_step_simple(
    vec::Buffers<T>& buf,
    T r,
    unsigned start,
    unsigned m,
    size_t len,
    T card)
{
    if (len == 0) {
        return;
    }
    const unsigned step = m << 1;
    const T rp1 = r + 1;
    VecType c = SET1(r);

    const size_t end = (len > 1) ? len - 1 : 0;
    const unsigned bufs_nb = buf.get_n();
    // #pragma omp parallel for
    // #pragma unroll
    const std::vector<T*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        VecType x1, y1;
        VecType x2, y2;
        VecType* __restrict p = reinterpret_cast<VecType*>(mem[i]);
        VecType* __restrict q = reinterpret_cast<VecType*>(mem[i + m]);

        // #pragma omp parallel for
        size_t j = 0;
        // #pragma unroll
        for (; j < end; j += 2) {
            x1 = LOAD(p + j);
            x2 = LOAD(p + j + 1);

            y1 = BUTTERFLY_GS_SIMPLE(rp1, c, x1, card);
            y2 = BUTTERFLY_GS_SIMPLE(rp1, c, x2, card);

            // Store back to memory
            STORE(q + j, y1);
            STORE(q + j + 1, y2);
        }
        for (; j < len; ++j) {
            x1 = LOAD(p + j);

            y1 = BUTTERFLY_GS_SIMPLE(rp1, c, x1, card);

            // Store back to memory
            STORE(q + j, y1);
        }
    }
}

template <typename T>
inline void encode_post_process(
    vec::Buffers<T>& output,
    std::vector<Properties>& props,
    off_t offset,
    unsigned code_len,
    T threshold,
    size_t vecs_nb)
{
    const unsigned element_size = sizeof(T);
    const unsigned vec_size = countof<T>();
    const T max = 1 << (element_size * 8 - 1);
    const VecType _threshold = SET1(threshold);
    const VecType mask_hi = SET1(max);

    // #pragma unroll
    const std::vector<T*>& mem = output.get_mem();
    for (unsigned frag_id = 0; frag_id < code_len; ++frag_id) {
        VecType* __restrict buf = reinterpret_cast<VecType*>(mem[frag_id]);

        size_t vec_id = 0;
        size_t end = (vecs_nb > 3) ? vecs_nb - 3 : 0;
        // #pragma unroll
        for (; vec_id < end; vec_id += 4) {
            VecType a1 = LOAD(buf + vec_id);
            VecType a2 = LOAD(buf + vec_id + 1);
            VecType a3 = LOAD(buf + vec_id + 2);
            VecType a4 = LOAD(buf + vec_id + 3);

            if (TESTZ(a1, _threshold) == 0) {
                const off_t curr_offset = offset + vec_id * vec_size;
                ADD_PROPS(
                    props[frag_id], _threshold, mask_hi, a1, curr_offset, max);
            }
            if (TESTZ(a2, _threshold) == 0) {
                const off_t curr_offset = offset + (vec_id + 1) * vec_size;
                ADD_PROPS(
                    props[frag_id], _threshold, mask_hi, a2, curr_offset, max);
            }
            if (TESTZ(a3, _threshold) == 0) {
                const off_t curr_offset = offset + (vec_id + 2) * vec_size;
                ADD_PROPS(
                    props[frag_id], _threshold, mask_hi, a3, curr_offset, max);
            }
            if (TESTZ(a4, _threshold) == 0) {
                const off_t curr_offset = offset + (vec_id + 3) * vec_size;
                ADD_PROPS(
                    props[frag_id], _threshold, mask_hi, a4, curr_offset, max);
            }
        }
        for (; vec_id < vecs_nb; ++vec_id) {
            VecType a = LOAD(buf + vec_id);
            uint32_t c = TESTZ(a, _threshold);
            if (c == 0) {
                const off_t curr_offset = offset + vec_id * vec_size;
                ADD_PROPS(
                    props[frag_id], _threshold, mask_hi, a, curr_offset, max);
            }
        }
    }
}

} // namespace simd
} // namespace quadiron

#endif
