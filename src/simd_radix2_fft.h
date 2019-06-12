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

#ifndef __QUAD_SIMD_RADIX2_FFT_H__
#define __QUAD_SIMD_RADIX2_FFT_H__

#include <x86intrin.h>

namespace quadiron {
namespace simd {

enum class CtGsCase {
    SIMPLE,
    NORMAL,
    EXTREME,
};

template <typename T>
constexpr CtGsCase get_case(T r, T q)
{
    if (r == 1) {
        return CtGsCase::SIMPLE;
    } else if (r < q - 1) {
        return CtGsCase::NORMAL;
    } else {
        return CtGsCase::EXTREME;
    }
}

/* ================= Vectorized Operations ================= */

/**
 * Butterfly Cooley-Tukey operation
 *
 * x <- x + r * y
 * y <- x - r * y
 *
 * @param ct_case coefficient case
 * @param c a register stores coefficient `r`
 * @param x working register
 * @param y working register
 */
template <typename T>
inline void
butterfly_ct(CtGsCase ct_case, const VecType& c, VecType& x, VecType& y)
{
    VecType z = y;
    switch (ct_case) {
    case CtGsCase::SIMPLE:
        y = mod_sub<T>(x, z);
        x = mod_add<T>(x, z);
        break;
    case CtGsCase::EXTREME:
        y = mod_add<T>(x, z);
        x = mod_sub<T>(x, z);
        break;
    case CtGsCase::NORMAL:
        z = mod_mul<T>(c, y);
        y = mod_sub<T>(x, z);
        x = mod_add<T>(x, z);
        break;
    }
}

/**
 * Butterfly Genteleman-Sande operation
 *
 * x <- x + y
 * y <- r * (x - y)
 *
 * @param gs_case coefficient case
 * @param c a register stores coefficient `r`
 * @param x working register
 * @param y working register
 */
template <typename T>
inline void
butterfly_gs(CtGsCase gs_case, const VecType& c, VecType& x, VecType& y)
{
    VecType add = mod_add<T>(x, y);
    switch (gs_case) {
    case CtGsCase::SIMPLE:
        y = mod_sub<T>(x, y);
        break;
    case CtGsCase::EXTREME:
        y = mod_sub<T>(y, x);
        break;
    case CtGsCase::NORMAL:
        VecType sub = mod_sub<T>(x, y);
        y = mod_mul<T>(c, sub);
        break;
    }
    x = add;
}

/**
 * Butterfly Genteleman-Sande simple operation where y = 0
 *
 * x <- x, i.e. no operation
 * y <- r * x
 *
 * @param gs_case coefficient case
 * @param c a register stores coefficient `r`
 * @param x working register
 */
template <typename T>
inline void butterfly_simple_gs(CtGsCase gs_case, const VecType& c, VecType& x)
{
    switch (gs_case) {
    case CtGsCase::EXTREME:
        x = mod_neg<T>(x);
        break;
    case CtGsCase::NORMAL:
        x = mod_mul<T>(c, x);
        break;
    case CtGsCase::SIMPLE:
        break;
    }
}
/**
 * Vectorized butterfly CT step
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
    const CtGsCase ct_case = get_case<T>(r, card);
    const VecType c = set_one(r);

    const size_t end = (len > 1) ? len - 1 : 0;
    const unsigned bufs_nb = buf.get_n();
    const std::vector<T*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        VecType* p = reinterpret_cast<VecType*>(mem[i]);
        VecType* q = reinterpret_cast<VecType*>(mem[i + m]);

        size_t j = 0;
        for (; j < end; j += 2) {
            VecType x1 = load_to_reg(p);
            VecType y1 = load_to_reg(q);

            butterfly_ct<T>(ct_case, c, x1, y1);

            VecType x2 = load_to_reg(p + 1);
            VecType y2 = load_to_reg(q + 1);

            butterfly_ct<T>(ct_case, c, x2, y2);

            // Store back to memory
            store_to_mem(p++, x1);
            store_to_mem(p++, x2);
            store_to_mem(q++, y1);
            store_to_mem(q++, y2);
        }
        for (; j < len; ++j) {
            VecType x1 = load_to_reg(p);
            VecType y1 = load_to_reg(q);

            butterfly_ct<T>(ct_case, c, x1, y1);

            // Store back to memory
            store_to_mem(p++, x1);
            store_to_mem(q++, y1);
        }
    }
}

template <typename T>
inline void do_butterfly_ct_2_layers(
    const std::vector<T*>& mem,
    T r1,
    T r2,
    T r3,
    unsigned start,
    unsigned m,
    size_t len,
    T card)
{
    const CtGsCase case1 = get_case<T>(r1, card);
    const CtGsCase case2 = get_case<T>(r2, card);
    const CtGsCase case3 = get_case<T>(r3, card);

    VecType c1 = set_one(r1);
    VecType c2 = set_one(r2);
    VecType c3 = set_one(r3);

    VecType* p = reinterpret_cast<VecType*>(mem[start]);
    VecType* q = reinterpret_cast<VecType*>(mem[start + m]);
    VecType* r = reinterpret_cast<VecType*>(mem[start + 2 * m]);
    VecType* s = reinterpret_cast<VecType*>(mem[start + 3 * m]);

    size_t j = 0;
    const size_t end = (len > 1) ? len - 1 : 0;
    while (j < end) {
        VecType x1 = load_to_reg(p);
        VecType y1 = load_to_reg(q);
        VecType u1 = load_to_reg(r);
        VecType v1 = load_to_reg(s);

        butterfly_ct<T>(case1, c1, x1, y1);
        butterfly_ct<T>(case1, c1, u1, v1);
        butterfly_ct<T>(case2, c2, x1, u1);
        butterfly_ct<T>(case3, c3, y1, v1);

        VecType x2 = load_to_reg(p + 1);
        VecType y2 = load_to_reg(q + 1);
        VecType u2 = load_to_reg(r + 1);
        VecType v2 = load_to_reg(s + 1);

        butterfly_ct<T>(case1, c1, x2, y2);
        butterfly_ct<T>(case1, c1, u2, v2);
        butterfly_ct<T>(case2, c2, x2, u2);
        butterfly_ct<T>(case3, c3, y2, v2);

        store_to_mem(p++, x1);
        store_to_mem(p++, x2);
        store_to_mem(q++, y1);
        store_to_mem(q++, y2);
        store_to_mem(r++, u1);
        store_to_mem(r++, u2);
        store_to_mem(s++, v1);
        store_to_mem(s++, v2);

        j += 2;
    }
    for (; j < len; ++j) {
        VecType x1 = load_to_reg(p);
        VecType y1 = load_to_reg(q);
        VecType u1 = load_to_reg(r);
        VecType v1 = load_to_reg(s);

        butterfly_ct<T>(case1, c1, x1, y1);
        butterfly_ct<T>(case1, c1, u1, v1);
        butterfly_ct<T>(case2, c2, x1, u1);
        butterfly_ct<T>(case3, c3, y1, v1);

        store_to_mem(p++, x1);
        store_to_mem(q++, y1);
        store_to_mem(r++, u1);
        store_to_mem(s++, v1);
    }
}

/**
 * Vectorized butterfly CT on two-layers at a time
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

    const std::vector<T*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        do_butterfly_ct_2_layers(mem, r1, r2, r3, i, m, len, card);
    }
}

/**
 * Vectorized butterfly GS step
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
    const CtGsCase gs_case = get_case<T>(r, card);
    VecType c = set_one(r);

    const unsigned bufs_nb = buf.get_n();
    const std::vector<T*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        VecType* p = reinterpret_cast<VecType*>(mem[i]);
        VecType* q = reinterpret_cast<VecType*>(mem[i + m]);

        for (size_t j = 0; j < len; ++j) {
            VecType x1 = load_to_reg(p);
            VecType y1 = load_to_reg(q);

            butterfly_gs<T>(gs_case, c, x1, y1);

            store_to_mem(p++, x1);
            store_to_mem(q++, y1);
        }
    }
}

/**
 * Vectorized butterfly GS step
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
    const CtGsCase gs_case = get_case<T>(r, card);
    VecType c = set_one(r);

    const unsigned bufs_nb = buf.get_n();
    const std::vector<T*>& mem = buf.get_mem();
    for (unsigned i = start; i < bufs_nb; i += step) {
        VecType* p = reinterpret_cast<VecType*>(mem[i]);
        VecType* q = reinterpret_cast<VecType*>(mem[i + m]);

        for (size_t j = 0; j < len; ++j) {
            VecType x = load_to_reg(p++);

            butterfly_simple_gs<T>(gs_case, c, x);

            store_to_mem(q++, x);
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
    const unsigned vec_size = countof<T>();
    const T max = 1U << (sizeof(T) * CHAR_BIT - 1);
    const VecType _threshold = set_one(threshold);
    const VecType mask_hi = set_one(max);

    const std::vector<T*>& mem = output.get_mem();
    for (unsigned frag_id = 0; frag_id < code_len; ++frag_id) {
        VecType* buf = reinterpret_cast<VecType*>(mem[frag_id]);

        size_t vec_id = 0;
        size_t end = (vecs_nb > 3) ? vecs_nb - 3 : 0;
        for (; vec_id < end; vec_id += 4) {
            VecType a1 = load_to_reg(buf + vec_id);
            VecType a2 = load_to_reg(buf + vec_id + 1);
            VecType a3 = load_to_reg(buf + vec_id + 2);
            VecType a4 = load_to_reg(buf + vec_id + 3);

            if (!and_is_zero(a1, _threshold)) {
                const off_t curr_offset = offset + vec_id * vec_size;
                add_props(
                    props[frag_id], _threshold, mask_hi, a1, curr_offset, max);
            }
            if (!and_is_zero(a2, _threshold)) {
                const off_t curr_offset = offset + (vec_id + 1) * vec_size;
                add_props(
                    props[frag_id], _threshold, mask_hi, a2, curr_offset, max);
            }
            if (!and_is_zero(a3, _threshold)) {
                const off_t curr_offset = offset + (vec_id + 2) * vec_size;
                add_props(
                    props[frag_id], _threshold, mask_hi, a3, curr_offset, max);
            }
            if (!and_is_zero(a4, _threshold)) {
                const off_t curr_offset = offset + (vec_id + 3) * vec_size;
                add_props(
                    props[frag_id], _threshold, mask_hi, a4, curr_offset, max);
            }
        }
        for (; vec_id < vecs_nb; ++vec_id) {
            VecType a = load_to_reg(buf + vec_id);
            if (!and_is_zero(a, _threshold)) {
                const off_t curr_offset = offset + vec_id * vec_size;
                add_props(
                    props[frag_id], _threshold, mask_hi, a, curr_offset, max);
            }
        }
    }
}

} // namespace simd
} // namespace quadiron

#endif
