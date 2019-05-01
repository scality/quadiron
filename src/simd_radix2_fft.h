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
inline CtGsCase get_case(T r, T q)
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
    default:
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
    default:
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
    default:
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

    const unsigned bufs_nb = buf.get_n();
    const std::vector<T*>& mem = buf.get_mem();
    const std::vector<uint8_t*>& meta = buf.get_meta();
    for (unsigned i = start; i < bufs_nb; i += step) {
        VecType* p = reinterpret_cast<VecType*>(mem[i]);
        VecType* q = reinterpret_cast<VecType*>(mem[i + m]);
        MetaType* m_p = reinterpret_cast<MetaType*>(meta[i]);
        MetaType* m_q = reinterpret_cast<MetaType*>(meta[i + m]);

        for (size_t j = 0; j < len; ++j) {
            VecType x1 = load_to_reg(p);
            VecType y1 = load_to_reg(q);

            VecType x1_lo, x1_hi;
            VecType y1_lo, y1_hi;

            unpack<T>(m_p[j], x1, x1_hi, x1_lo);
            unpack<T>(m_q[j], y1, y1_hi, y1_lo);

            butterfly_ct<T>(ct_case, c, x1_lo, y1_lo);
            butterfly_ct<T>(ct_case, c, x1_hi, y1_hi);

            pack<T>(x1_lo, x1_hi, x1, m_p[j]);
            pack<T>(y1_lo, y1_hi, y1, m_q[j]);

            // Store back to memory
            store_to_mem(p++, x1);
            store_to_mem(q++, y1);
        }
    }
}

template <typename T>
inline void do_butterfly_ct_2_layers(
    const std::vector<T*>& mem,
    const std::vector<uint8_t*>& meta,
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

    MetaType* m_p = reinterpret_cast<MetaType*>(meta[start]);
    MetaType* m_q = reinterpret_cast<MetaType*>(meta[start + m]);
    MetaType* m_r = reinterpret_cast<MetaType*>(meta[start + 2 * m]);
    MetaType* m_s = reinterpret_cast<MetaType*>(meta[start + 3 * m]);

    for (size_t j = 0; j < len; ++j) {
        VecType x1 = load_to_reg(p);
        VecType y1 = load_to_reg(q);
        VecType u1 = load_to_reg(r);
        VecType v1 = load_to_reg(s);

        VecType x1_lo, x1_hi;
        VecType y1_lo, y1_hi;
        VecType u1_lo, u1_hi;
        VecType v1_lo, v1_hi;

        unpack<T>(m_p[j], x1, x1_hi, x1_lo);
        unpack<T>(m_q[j], y1, y1_hi, y1_lo);
        unpack<T>(m_r[j], u1, u1_hi, u1_lo);
        unpack<T>(m_s[j], v1, v1_hi, v1_lo);

        butterfly_ct<T>(case1, c1, x1_lo, y1_lo);
        butterfly_ct<T>(case1, c1, x1_hi, y1_hi);
        butterfly_ct<T>(case1, c1, u1_lo, v1_lo);
        butterfly_ct<T>(case1, c1, u1_hi, v1_hi);

        butterfly_ct<T>(case2, c2, x1_lo, u1_lo);
        butterfly_ct<T>(case2, c2, x1_hi, u1_hi);

        butterfly_ct<T>(case3, c3, y1_lo, v1_lo);
        butterfly_ct<T>(case3, c3, y1_hi, v1_hi);

        pack<T>(x1_lo, x1_hi, x1, m_p[j]);
        pack<T>(y1_lo, y1_hi, y1, m_q[j]);
        pack<T>(u1_lo, u1_hi, u1, m_r[j]);
        pack<T>(v1_lo, v1_hi, v1, m_s[j]);

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
    const std::vector<uint8_t*>& meta = buf.get_meta();
    for (unsigned i = start; i < bufs_nb; i += step) {
        do_butterfly_ct_2_layers(mem, meta, r1, r2, r3, i, m, len, card);
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
    const std::vector<uint8_t*>& meta = buf.get_meta();
    for (unsigned i = start; i < bufs_nb; i += step) {
        VecType* p = reinterpret_cast<VecType*>(mem[i]);
        VecType* q = reinterpret_cast<VecType*>(mem[i + m]);
        MetaType* m_p = reinterpret_cast<MetaType*>(meta[i]);
        MetaType* m_q = reinterpret_cast<MetaType*>(meta[i + m]);

        for (size_t j = 0; j < len; ++j) {
            VecType x1 = load_to_reg(p);
            VecType y1 = load_to_reg(q);

            VecType x1_lo, x1_hi;
            VecType y1_lo, y1_hi;

            unpack<T>(m_p[j], x1, x1_hi, x1_lo);
            unpack<T>(m_q[j], y1, y1_hi, y1_lo);

            butterfly_gs<T>(gs_case, c, x1_lo, y1_lo);
            butterfly_gs<T>(gs_case, c, x1_hi, y1_hi);

            pack<T>(x1_lo, x1_hi, x1, m_p[j]);
            pack<T>(y1_lo, y1_hi, y1, m_q[j]);

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
    const std::vector<uint8_t*>& meta = buf.get_meta();
    for (unsigned i = start; i < bufs_nb; i += step) {
        VecType* p = reinterpret_cast<VecType*>(mem[i]);
        VecType* q = reinterpret_cast<VecType*>(mem[i + m]);
        MetaType* m_p = reinterpret_cast<MetaType*>(meta[i]);
        MetaType* m_q = reinterpret_cast<MetaType*>(meta[i + m]);

        for (size_t j = 0; j < len; ++j) {
            VecType x = load_to_reg(p++);
            VecType x_lo, x_hi;

            unpack<T>(m_p[j], x, x_hi, x_lo);

            butterfly_simple_gs<T>(gs_case, c, x_lo);
            butterfly_simple_gs<T>(gs_case, c, x_hi);

            pack<T>(x_lo, x_hi, x, m_q[j]);

            store_to_mem(q++, x);
        }
    }
}

} // namespace simd
} // namespace quadiron

#endif
