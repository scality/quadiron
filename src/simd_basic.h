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

#ifndef __QUAD_SIMD_BASIC_H__
#define __QUAD_SIMD_BASIC_H__

#include <x86intrin.h>

namespace quadiron {
namespace simd {

template <typename T>
inline VecType card(T q);
template <>
inline VecType card<uint16_t>(uint16_t q)
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
inline VecType card_minus_one<uint16_t>(uint16_t q)
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
 * @param max a dummy variable
 */
template <typename T>
inline void add_props(
    Properties& props,
    VecType threshold,
    VecType mask,
    VecType symb,
    off_t offset,
    T max)
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

/* ==================== Operations for RingModN =================== */
/** Perform a multiplication of a coefficient `a` to each element of `src` and
 *  add result to correspondent element of `dest`
 *
 * @note: 1 < `a` < card - 1
 */
template <typename T>
inline void mul_coef_to_buf(const T a, T* src, T* dest, size_t len, T card)
{
    const VecType coef = set_one(a);

    VecType* __restrict _src = reinterpret_cast<VecType*>(src);
    VecType* __restrict _dest = reinterpret_cast<VecType*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i = 0;
    const size_t end = (_len > 3) ? _len - 3 : 0;
    for (; i < end; i += 4) {
        _dest[i] = mod_mul(coef, _src[i], card);
        _dest[i + 1] = mod_mul(coef, _src[i + 1], card);
        _dest[i + 2] = mod_mul(coef, _src[i + 2], card);
        _dest[i + 3] = mod_mul(coef, _src[i + 3], card);
    }
    for (; i < _len; ++i) {
        _dest[i] = mod_mul(coef, _src[i], card);
    }

    if (_last_len > 0) {
        const DoubleSizeVal<T> coef_double = DoubleSizeVal<T>(a);
        for (size_t i = _len * ratio; i < len; i++) {
            dest[i] = (T)((coef_double * src[i]) % card);
        }
    }
}

template <typename T>
inline void add_two_bufs(T* src, T* dest, size_t len, T card)
{
    VecType* __restrict _src = reinterpret_cast<VecType*>(src);
    VecType* __restrict _dest = reinterpret_cast<VecType*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        _dest[i] = mod_add(_src[i], _dest[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
            const T tmp = src[i] + dest[i];
            dest[i] = (tmp >= card) ? (tmp - card) : tmp;
        }
    }
}

template <typename T>
inline void sub_two_bufs(T* bufa, T* bufb, T* res, size_t len, T card)
{
    VecType* __restrict _bufa = reinterpret_cast<VecType*>(bufa);
    VecType* __restrict _bufb = reinterpret_cast<VecType*>(bufb);
    VecType* __restrict _res = reinterpret_cast<VecType*>(res);
    const unsigned ratio = sizeof(*_bufa) / sizeof(*bufa);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform subtraction
        _res[i] = mod_sub(_bufa[i], _bufb[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
            // perform subtraction
            if (bufa[i] >= bufb[i]) {
                res[i] = bufa[i] - bufb[i];
            } else {
                res[i] = card - (bufb[i] - bufa[i]);
            }
        }
    }
}

template <typename T>
inline void mul_two_bufs(T* src, T* dest, size_t len, T card)
{
    VecType* __restrict _src = reinterpret_cast<VecType*>(src);
    VecType* __restrict _dest = reinterpret_cast<VecType*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform multiplicaton
        _dest[i] = mod_mul_safe(_src[i], _dest[i], card);
    }
    if (_last_len > 0) {
        for (i = _len * ratio; i < len; i++) {
            // perform multiplicaton
            dest[i] = T((DoubleSizeVal<T>(src[i]) * dest[i]) % card);
        }
    }
}

/** Apply an element-wise negation to a buffer
 */
template <typename T>
inline void neg(size_t len, T* buf, T card)
{
    VecType* _buf = reinterpret_cast<VecType*>(buf);
    const unsigned ratio = sizeof(*_buf) / sizeof(*buf);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        _buf[i] = mod_neg(_buf[i], card);
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
