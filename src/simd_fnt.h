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
inline VecType card();
template <>
inline VecType card<uint16_t>()
{
    return set_one<uint16_t>(257);
}
template <>
inline VecType card<uint32_t>()
{
    return set_one<uint32_t>(65537);
}

template <typename T>
inline VecType card_minus_one();
template <>
inline VecType card_minus_one<uint16_t>()
{
    return set_one<uint16_t>(256);
}
template <>
inline VecType card_minus_one<uint32_t>()
{
    return set_one<uint32_t>(65536);
}

const int I_MASK8_LO = 0b01010101;

template <typename T>
inline VecType get_low_half(const VecType& x);
template <>
inline VecType get_low_half<uint16_t>(const VecType& x)
{
    return BLEND8(zero(), x, MASK8_LO);
}
template <>
inline VecType get_low_half<uint32_t>(const VecType& x)
{
    return BLEND16(zero(), x, I_MASK8_LO);
}

template <typename T>
inline VecType get_high_half(const VecType& x);
template <>
inline VecType get_high_half<uint16_t>(const VecType& x)
{
    return BLEND8(zero(), SHIFTR(x, 1), MASK8_LO);
}
template <>
inline VecType get_high_half<uint32_t>(const VecType& x)
{
    return BLEND16(zero(), SHIFTR(x, 2), I_MASK8_LO);
}

/* ================= Basic Operations ================= */

/**
 * Modular addition
 *
 * @param x input register
 * @param y input register
 * @return (x + y) mod q
 */
template <typename T>
inline VecType mod_add(const VecType& x, const VecType& y)
{
    const VecType res = add<T>(x, y);
    return min<T>(res, sub<T>(res, card<T>()));
}

/**
 * Modular subtraction for packed unsigned 32-bit integers
 *
 * @param x input register
 * @param y input register
 * @return (x - y) mod q
 */
template <typename T>
inline VecType mod_sub(const VecType& x, const VecType& y)
{
    const VecType res = sub<T>(x, y);
    return min<T>(res, add<T>(res, card<T>()));
}

/**
 * Modular negation for packed unsigned 32-bit integers
 *
 * @param x input register
 * @return (-x) mod q
 */
template <typename T>
inline VecType mod_neg(const VecType& x)
{
    const VecType res = sub<T>(card<T>(), x);
    return min<T>(res, sub<T>(res, card<T>()));
}

/**
 * Modular multiplication for packed unsigned 32-bit integers
 *
 * @note We assume that at least `x` or `y` is less than `q-1` so it's
 * not necessary to verify overflow on multiplying elements
 *
 * @param x input register
 * @param y input register
 * @return (x * y) mod q
 */
template <typename T>
inline VecType mod_mul(const VecType& x, const VecType& y)
{
    const VecType res = mul<T>(x, y);
    const VecType lo = get_low_half<T>(res);
    const VecType hi = get_high_half<T>(res);
    return mod_sub<T>(lo, hi);
}

/**
 * Modular general multiplication for packed unsigned 32-bit integers
 *
 * @note It's necessary to verify overflow on multiplying elements
 *
 * @param x input register
 * @param y input register
 * @return (x * y) mod q
 */
template <typename T>
inline VecType mod_mul_safe(const VecType& x, const VecType& y)
{
    const VecType res = mod_mul<T>(x, y);

    // filter elements of both of a & b = card-1
    const VecType cmp = bit_and(
        compare_eq<T>(x, card_minus_one<T>()),
        compare_eq<T>(y, card_minus_one<T>()));

    if (is_zero(cmp)) {
        return res;
    }
    return add<T>(res, bit_and(one<T>(), cmp));
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
