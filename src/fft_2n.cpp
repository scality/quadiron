/* -*- mode: c++ -*- */
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

#include "fft_2n.h"

/*
 * The file includes vectorized operations used by Radix2 classes
 */

#ifdef QUADIRON_USE_SIMD

#include "simd.h"

namespace quadiron {
namespace fft {

template <>
void Radix2<uint16_t>::butterfly_ct_two_layers_step(
    vec::Buffers<uint16_t>& buf,
    unsigned start,
    unsigned m)
{
    const unsigned ratio = simd::countof<uint16_t>();
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;
    const unsigned coefIndex = start * this->n / m / 2;
    const uint16_t r1 = vec_W[coefIndex];
    const uint16_t r2 = vec_W[coefIndex / 2];
    const uint16_t r3 = vec_W[coefIndex / 2 + this->n / 4];

    // perform vector operations
    simd::butterfly_ct_two_layers_step(
        buf, r1, r2, r3, start, m, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        const unsigned step = m << 2;
        size_t offset = vec_len * ratio;
        //  ---------
        // First layer
        //  ---------
        const uint16_t r1 = W->get(start * this->n / m / 2);
        // first pair
        butterfly_ct_step_slow(buf, r1, start, m, step, offset);
        // second pair
        butterfly_ct_step_slow(buf, r1, start + 2 * m, m, step, offset);
        //  ---------
        // Second layer
        //  ---------
        // first pair
        const uint16_t r2 = W->get(start * this->n / m / 4);
        butterfly_ct_step_slow(buf, r2, start, 2 * m, step, offset);
        // second pair
        const uint16_t r3 = W->get((start + m) * this->n / m / 4);
        butterfly_ct_step_slow(buf, r3, start + m, 2 * m, step, offset);
    }
}

template <>
void Radix2<uint16_t>::butterfly_ct_step(
    vec::Buffers<uint16_t>& buf,
    uint16_t r,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = simd::countof<uint16_t>();
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;

    // perform vector operations
    simd::butterfly_ct_step(buf, r, start, m, step, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        size_t offset = vec_len * ratio;
        butterfly_ct_step_slow(buf, r, start, m, step, offset);
    }
}

template <>
void Radix2<uint16_t>::butterfly_gs_step(
    vec::Buffers<uint16_t>& buf,
    uint16_t coef,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = simd::countof<uint16_t>();
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;

    // perform vector operations
    simd::butterfly_gs_step(buf, coef, start, m, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        size_t offset = vec_len * ratio;
        butterfly_gs_step_slow(buf, coef, start, m, step, offset);
    }
}

template <>
void Radix2<uint16_t>::butterfly_gs_step_simple(
    vec::Buffers<uint16_t>& buf,
    uint16_t coef,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = simd::countof<uint16_t>();
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;

    // perform vector operations
    simd::butterfly_gs_step_simple(buf, coef, start, m, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        size_t offset = vec_len * ratio;
        butterfly_gs_step_simple_slow(buf, coef, start, m, step, offset);
    }
}

template <>
void Radix2<uint32_t>::butterfly_ct_two_layers_step(
    vec::Buffers<uint32_t>& buf,
    unsigned start,
    unsigned m)
{
    const unsigned ratio = simd::countof<uint32_t>();
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;
    const unsigned coefIndex = start * this->n / m / 2;
    const uint32_t r1 = vec_W[coefIndex];
    const uint32_t r2 = vec_W[coefIndex / 2];
    const uint32_t r3 = vec_W[coefIndex / 2 + this->n / 4];

    // perform vector operations
    simd::butterfly_ct_two_layers_step(
        buf, r1, r2, r3, start, m, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        const unsigned step = m << 2;
        size_t offset = vec_len * ratio;
        //  ---------
        // First layer
        //  ---------
        const uint32_t r1 = W->get(start * this->n / m / 2);
        // first pair
        butterfly_ct_step_slow(buf, r1, start, m, step, offset);
        // second pair
        butterfly_ct_step_slow(buf, r1, start + 2 * m, m, step, offset);
        //  ---------
        // Second layer
        //  ---------
        // first pair
        const uint32_t r2 = W->get(start * this->n / m / 4);
        butterfly_ct_step_slow(buf, r2, start, 2 * m, step, offset);
        // second pair
        const uint32_t r3 = W->get((start + m) * this->n / m / 4);
        butterfly_ct_step_slow(buf, r3, start + m, 2 * m, step, offset);
    }
}

template <>
void Radix2<uint32_t>::butterfly_ct_step(
    vec::Buffers<uint32_t>& buf,
    uint32_t r,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = simd::countof<uint32_t>();
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;

    // perform vector operations
    simd::butterfly_ct_step(buf, r, start, m, step, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        size_t offset = vec_len * ratio;
        butterfly_ct_step_slow(buf, r, start, m, step, offset);
    }
}

template <>
void Radix2<uint32_t>::butterfly_gs_step(
    vec::Buffers<uint32_t>& buf,
    uint32_t coef,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = simd::countof<uint32_t>();
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;

    // perform vector operations
    simd::butterfly_gs_step(buf, coef, start, m, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        size_t offset = vec_len * ratio;
        butterfly_gs_step_slow(buf, coef, start, m, step, offset);
    }
}

template <>
void Radix2<uint32_t>::butterfly_gs_step_simple(
    vec::Buffers<uint32_t>& buf,
    uint32_t coef,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = simd::countof<uint32_t>();
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;

    // perform vector operations
    simd::butterfly_gs_step_simple(buf, coef, start, m, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        size_t offset = vec_len * ratio;
        butterfly_gs_step_simple_slow(buf, coef, start, m, step, offset);
    }
}

} // namespace fft
} // namespace quadiron

#endif // #ifdef QUADIRON_USE_SIMD
