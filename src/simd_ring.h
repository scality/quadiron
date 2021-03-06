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

#ifndef __QUAD_SIMD_RING_H__
#define __QUAD_SIMD_RING_H__

#include <x86intrin.h>

namespace quadiron {
namespace simd {

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

    VecType* _src = reinterpret_cast<VecType*>(src);
    VecType* _dest = reinterpret_cast<VecType*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i = 0;
    const size_t end = (_len > 3) ? _len - 3 : 0;
    for (; i < end; i += 4) {
        _dest[i] = mod_mul<T>(coef, _src[i]);
        _dest[i + 1] = mod_mul<T>(coef, _src[i + 1]);
        _dest[i + 2] = mod_mul<T>(coef, _src[i + 2]);
        _dest[i + 3] = mod_mul<T>(coef, _src[i + 3]);
    }
    for (; i < _len; ++i) {
        _dest[i] = mod_mul<T>(coef, _src[i]);
    }

    if (_last_len > 0) {
        const DoubleSizeVal<T> coef_double = DoubleSizeVal<T>(a);
        for (size_t i = _len * ratio; i < len; i++) {
            dest[i] = static_cast<T>((coef_double * src[i]) % card);
        }
    }
}

template <typename T>
inline void add_two_bufs(T* src, T* dest, size_t len, T card)
{
    VecType* _src = reinterpret_cast<VecType*>(src);
    VecType* _dest = reinterpret_cast<VecType*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        _dest[i] = mod_add<T>(_src[i], _dest[i]);
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
    VecType* _bufa = reinterpret_cast<VecType*>(bufa);
    VecType* _bufb = reinterpret_cast<VecType*>(bufb);
    VecType* _res = reinterpret_cast<VecType*>(res);
    const unsigned ratio = sizeof(*_bufa) / sizeof(*bufa);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform subtraction
        _res[i] = mod_sub<T>(_bufa[i], _bufb[i]);
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
    VecType* _src = reinterpret_cast<VecType*>(src);
    VecType* _dest = reinterpret_cast<VecType*>(dest);
    const unsigned ratio = sizeof(*_src) / sizeof(*src);
    const size_t _len = len / ratio;
    const size_t _last_len = len - _len * ratio;

    size_t i;
    for (i = 0; i < _len; i++) {
        // perform multiplicaton
        _dest[i] = mod_mul_safe<T>(_src[i], _dest[i]);
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
        _buf[i] = mod_neg<T>(_buf[i]);
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
