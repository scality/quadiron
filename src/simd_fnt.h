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

#include "vec_buffers.h"

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
    return BLEND16(zero(), x, 0x55);
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
    return BLEND16(zero(), SHIFTR(x, 2), 0x55);
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
 * Update property given an output Buffers
 *
 * @param output output Buffers
 * @param props properties bound to fragments
 * @param offset offset in the data fragments
 * @param code_len erasure codes' length
 * @param vecs_nb number of vectors corresponding to the data
 */
template <typename T>
inline void encode_post_process(
    vec::Buffers<T>& output,
    std::vector<Properties>& props,
    off_t offset,
    unsigned code_len,
    size_t vecs_nb)
{
    // nb of elements per vector
    const unsigned vec_size = countof<T>();
    // size of meta element in bits
    const unsigned ele_size_in_bits = sizeof(MetaType) * CHAR_BIT / vec_size;
    // mask to get meta
    const T mask = ((static_cast<T>(1) << ele_size_in_bits) - 1);

    const std::vector<uint8_t*>& meta = output.get_meta();
    for (unsigned frag_id = 0; frag_id < code_len; ++frag_id) {
        const MetaType* meta_frag = reinterpret_cast<MetaType*>(meta[frag_id]);
        for (size_t vec_id = 0; vec_id < vecs_nb; ++vec_id) {
            if (meta_frag[vec_id]) {
                const off_t curr_offset = offset + vec_id * vec_size;
                MetaType val = meta_frag[vec_id];
                unsigned idx = 0;
                while (val) {
                    const T m_val = val & mask;
                    if (m_val) {
                        const size_t _offset = curr_offset + idx;
                        props[frag_id].add(_offset, m_val);
                    }

                    val >>= ele_size_in_bits;
                    idx++;
                }
            }
        }
    }
}

} // namespace simd
} // namespace quadiron

#endif
