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
inline VecType CARD(T q)
{
    return (q == F3) ? F3_u32 : F4_u32;
}

template <typename T>
inline VecType CARD_M_1(T q)
{
    return (q == F3) ? F3m1_u32 : F4m1_u32;
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
inline VecType ADD_MOD(VecType x, VecType y, T q)
{
    const VecType res = Add<T>(x, y);
    return MIN<T>(res, Sub<T>(res, CARD(q)));
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
inline VecType SUB_MOD(VecType x, VecType y, T q)
{
    const VecType res = Sub<T>(x, y);
    return MIN<T>(res, Add<T>(res, CARD(q)));
}

/**
 * Modular negation for packed unsigned 32-bit integers
 *
 * @param x input register
 * @param q modulo
 * @return (-x) mod q
 */
template <typename T>
inline VecType NEG_MOD(VecType x, T q)
{
    const VecType res = Sub<T>(CARD(q), x);
    return MIN<T>(res, Sub<T>(res, CARD(q)));
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
inline VecType MUL_MOD(VecType x, VecType y, T q)
{
    const VecType res = Mul<T>(x, y);
    const VecType lo =
        (q == F3) ? BLEND8(ZERO, res, MASK8_LO) : BLEND16(ZERO, res, 0x55);
    const VecType hi = (q == F3) ? BLEND8(ZERO, SHIFTR(res, 1), MASK8_LO)
                                 : BLEND16(ZERO, SHIFTR(res, 2), 0x55);
    return SUB_MOD(lo, hi, q);
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
inline VecType MULFULL_MOD(VecType x, VecType y, T q)
{
    const VecType res = MUL_MOD(x, y, q);

    // filter elements of both of a & b = card-1
    const VecType cmp =
        And(CompareEq<T>(x, CARD_M_1(q)), CompareEq<T>(y, CARD_M_1(q)));

    if (IsZero(cmp) == 1) {
        return res;
    }
    return (q == F3) ? Xor(res, And(F4_u32, cmp))
                     : Add<T>(res, And(ONE32, cmp));
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
inline void ADD_PROPS(
    Properties& props,
    VecType threshold,
    VecType mask,
    VecType symb,
    off_t offset,
    T max)
{
    const VecType b = CompareEq<T>(threshold, symb);
    const VecType c = And(mask, b);
    auto d = Msb8Mask(c);
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
    const VecType coef = SetOne(a);

    VecType* __restrict _src = reinterpret_cast<VecType*>(src);
    VecType* __restrict _dest = reinterpret_cast<VecType*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i = 0;
    const size_t end = (_len > 3) ? _len - 3 : 0;
    for (; i < end; i += 4) {
        _dest[i] = MUL_MOD(coef, _src[i], card);
        _dest[i + 1] = MUL_MOD(coef, _src[i + 1], card);
        _dest[i + 2] = MUL_MOD(coef, _src[i + 2], card);
        _dest[i + 3] = MUL_MOD(coef, _src[i + 3], card);
    }
    for (; i < _len; ++i) {
        _dest[i] = MUL_MOD(coef, _src[i], card);
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
        _dest[i] = ADD_MOD(_src[i], _dest[i], card);
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
        _res[i] = SUB_MOD(_bufa[i], _bufb[i], card);
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
        _dest[i] = MULFULL_MOD(_src[i], _dest[i], card);
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
