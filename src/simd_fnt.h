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

template <typename T>
inline VecType card(T q);
template <>
inline VecType card<uint16_t>(uint16_t)
{
    return F3_U16;
}
template <>
inline VecType card<uint32_t>(uint32_t q)
{
    return (q == F3) ? F3_U32 : F4_U32;
}

template <typename T>
inline VecType card_minus_one(T q);
template <>
inline VecType card_minus_one<uint16_t>(uint16_t)
{
    return F3_MINUS_ONE_U16;
}
template <>
inline VecType card_minus_one<uint32_t>(uint32_t q)
{
    return (q == F3) ? F3_MINUS_ONE_U32 : F4_MINUS_ONE_U32;
}

template <typename T>
inline VecType get_low_half(VecType x, T q)
{
    return (q == F3) ? BLEND8(ZERO, x, MASK8_LO) : BLEND16(ZERO, x, 0x55);
}

template <typename T>
inline VecType get_high_half(VecType x, T q)
{
    return (q == F3) ? BLEND8(ZERO, SHIFTR(x, 1), MASK8_LO)
                     : BLEND16(ZERO, SHIFTR(x, 2), 0x55);
}

/* ================= Basic Operations ================= */

/**
 * Modular addition
 *
 * @param x input register
 * @param y input register
 * @param q modulo
 * @return (x + y) mod q
 */
template <typename T>
inline VecType mod_add(VecType x, VecType y, T q)
{
    const VecType res = add<T>(x, y);
    return min<T>(res, sub<T>(res, card(q)));
}

/**
 * Modular subtraction for packed unsigned 32-bit integers
 *
 * @param x input register
 * @param y input register
 * @param q modulo
 * @return (x - y) mod q
 */
template <typename T>
inline VecType mod_sub(VecType x, VecType y, T q)
{
    const VecType res = sub<T>(x, y);
    return min<T>(res, add<T>(res, card(q)));
}

/**
 * Modular negation for packed unsigned 32-bit integers
 *
 * @param x input register
 * @param q modulo
 * @return (-x) mod q
 */
template <typename T>
inline VecType mod_neg(VecType x, T q)
{
    const VecType res = sub<T>(card(q), x);
    return min<T>(res, sub<T>(res, card(q)));
}

/**
 * Modular multiplication for packed unsigned 32-bit integers
 *
 * @note We assume that at least `x` or `y` is less than `q-1` so it's
 * not necessary to verify overflow on multiplying elements
 *
 * @param x input register
 * @param y input register
 * @param q modulo
 * @return (x * y) mod q
 */
template <typename T>
inline VecType mod_mul(VecType x, VecType y, T q)
{
    const VecType res = mul<T>(x, y);
    const VecType lo = get_low_half(res, q);
    const VecType hi = get_high_half(res, q);
    return mod_sub(lo, hi, q);
}

/**
 * Modular general multiplication for packed unsigned 32-bit integers
 *
 * @note It's necessary to verify overflow on multiplying elements
 *
 * @param x input register
 * @param y input register
 * @param q modulo
 * @return (x * y) mod q
 */
template <typename T>
inline VecType mod_mul_safe(VecType x, VecType y, T q)
{
    const VecType res = mod_mul(x, y, q);

    // filter elements of both of a & b = card-1
    const VecType cmp = bit_and(
        compare_eq<T>(x, card_minus_one(q)),
        compare_eq<T>(y, card_minus_one(q)));

    if (is_zero(cmp)) {
        return res;
    }
    return (q == F3) ? bit_xor(res, bit_and(F4_U32, cmp))
                     : add<T>(res, bit_and(ONE_U32, cmp));
}

/**
 * Update property for a given register for packed unsigned 32-bit integers
 *
 * @param props properties bound to fragments
 * @param threshold register storing max value in its elements
 * @param mask a specific mask
 * @param symb input register
 * @param offset offset in the data fragments
 */
template <typename T>
inline void add_props(
    Properties& props,
    VecType threshold,
    VecType mask,
    VecType symb,
    off_t offset,
    T)
{
    const VecType b = compare_eq<T>(threshold, symb);
    const VecType c = bit_and(mask, b);
    auto d = msb8_mask(c);
    const unsigned element_size = sizeof(T);
    while (d > 0) {
        const unsigned byte_idx = __builtin_ctz(d);
        const size_t _offset = offset + byte_idx / element_size;
        props.add(_offset, OOR_MARK);
        d ^= 1 << byte_idx;
    }
}

} // namespace simd
} // namespace quadiron

#endif
