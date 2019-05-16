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

#include "gf_ring.h"

namespace quadiron {
namespace simd {

/* ==================== Operations for RingModN =================== */
/** Perform a multiplication of a coefficient `c` to each element of `buf_id`th
 * buffers of `src` and store result to correspondent element of `dest`
 *
 * @note: Buffers `src` and `dest` have meta
 * @note: 1 < `a` < card - 1
 */
template <typename T>
inline void mul_coef_to_buf(
    const gf::RingModN<T>& gf,
    T c,
    vec::Buffers<T>& src,
    vec::Buffers<T>& dest,
    size_t buf_id)
{
    const size_t size = src.get_size();
    const unsigned ratio = simd::countof<T>();
    const size_t simd_vec_len = size / ratio;
    const size_t simd_trailing_len = size - simd_vec_len * ratio;

    const VecType coef = set_one(c);

    VecType* s_data = reinterpret_cast<VecType*>(src.get(buf_id));
    VecType* d_data = reinterpret_cast<VecType*>(dest.get(buf_id));

    MetaType* s_meta = reinterpret_cast<MetaType*>(src.get_meta(buf_id));
    MetaType* d_meta = reinterpret_cast<MetaType*>(dest.get_meta(buf_id));

    for (size_t i = 0; i < simd_vec_len; ++i) {
        VecType lo, hi;
        VecType x = load_to_reg(&s_data[i]);

        unpack<T>(s_meta[i], x, hi, lo);
        hi = mod_mul<T>(coef, hi);
        lo = mod_mul<T>(coef, lo);
        pack<T>(lo, hi, x, d_meta[i]);

        store_to_mem(&d_data[i], x);
    }

    if (simd_trailing_len) {
        const size_t simd_offset = simd_vec_len * ratio;
        for (size_t i = simd_offset; i < size; ++i) {
            T hi, lo;
            src.get(buf_id, i, hi, lo);
            dest.set(buf_id, i, gf.mul(c, hi), gf.mul(c, lo));
        }
    }
}

/** Apply an element-wise negation to a buffer
 * @note: Buffers `src` and `dest` have meta
 */
template <typename T>
inline void neg(const gf::RingModN<T>& gf, vec::Buffers<T>& buf, size_t buf_id)
{
    const size_t size = buf.get_size();
    const unsigned ratio = simd::countof<T>();
    const size_t simd_vec_len = size / ratio;
    const size_t simd_trailing_len = size - simd_vec_len * ratio;

    VecType* vec_data = reinterpret_cast<VecType*>(buf.get(buf_id));
    MetaType* vec_meta = reinterpret_cast<MetaType*>(buf.get_meta(buf_id));

    for (size_t i = 0; i < simd_vec_len; ++i) {
        VecType lo, hi;
        VecType x = load_to_reg(&vec_data[i]);

        unpack<T>(vec_meta[i], x, hi, lo);
        hi = mod_neg<T>(hi);
        lo = mod_neg<T>(lo);
        pack<T>(lo, hi, x, vec_meta[i]);

        store_to_mem(&vec_data[i], x);
    }

    if (simd_trailing_len) {
        const size_t simd_offset = simd_vec_len * ratio;
        for (size_t i = simd_offset; i < size; ++i) {
            T hi, lo;
            buf.get(buf_id, i, hi, lo);
            buf.set(buf_id, i, gf.neg(hi), gf.neg(lo));
        }
    }
}

} // namespace simd
} // namespace quadiron

#endif
