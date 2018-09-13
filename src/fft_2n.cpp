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
void Radix2<uint32_t>::butterfly_ct_two_layers_step(
    vec::Buffers<uint32_t>& buf,
    unsigned start,
    unsigned m)
{
    const unsigned ratio = ALIGN_SIZE / sizeof(uint32_t);
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;
    const uint32_t card = this->gf->card();
    const uint32_t coefIndex = start * this->n / m / 2;
    const uint32_t r1 = W->get(coefIndex);
    const uint32_t r2 = W->get(coefIndex / 2);
    const uint32_t r3 = W->get(coefIndex / 2 + this->n / 4);

    // perform vector operations
    simd::butterfly_ct_two_layers_step(buf, r1, r2, r3, start, m, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        const unsigned step = m << 2;
        size_t offset = vec_len * ratio;
        //  ---------
        // First layer
        //  ---------
        const uint32_t r1 = W->get(start * this->n / m / 2);
        // first pair
        butterfly_ct_3_offset(r1, buf, start, m, step, offset);
        // second pair
        butterfly_ct_3_offset(r1, buf, start + 2 * m, m, step, offset);
        //  ---------
        // Second layer
        //  ---------
        // first pair
        const uint32_t r2 = W->get(start * this->n / m / 4);
        butterfly_ct_3_offset(r2, buf, start, 2 * m, step, offset);
        // second pair
        const uint32_t r3 = W->get((start + m) * this->n / m / 4);
        butterfly_ct_3_offset(r3, buf, start + m, 2 * m, step, offset);
    }
}


template <>
void Radix2<uint32_t>::butterfly_ct_1(
    vec::Buffers<uint32_t>& buf,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = ALIGN_SIZE / sizeof(uint32_t);
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;
    const uint32_t card = this->gf->card();

    // perform vector operations
    simd::butterfly_ct_1(buf, start, m, step, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        for (int i = start; i < this->n; i += step) {
            uint32_t* a = buf.get(i);
            uint32_t* b = buf.get(i + m);
            // perform butterfly operation for Cooley-Tukey FFT algorithm
            for (size_t j = vec_len * ratio; j < len; ++j) {
                uint32_t x = this->gf->add(a[j], b[j]);
                b[j] = this->gf->sub(a[j], b[j]);
                a[j] = x;
            }
        }
    }
}

template <>
void Radix2<uint32_t>::butterfly_ct_2(
    vec::Buffers<uint32_t>& buf,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = ALIGN_SIZE / sizeof(uint32_t);
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;
    const uint32_t card = this->gf->card();

    // perform vector operations
    simd::butterfly_ct_2(buf, start, m, step, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        for (int i = start; i < this->n; i += step) {
            uint32_t* a = buf.get(i);
            uint32_t* b = buf.get(i + m);
            // perform butterfly operation for Cooley-Tukey FFT algorithm
            for (size_t j = vec_len * ratio; j < len; ++j) {
                uint32_t x = this->gf->sub(a[j], b[j]);
                b[j] = this->gf->add(a[j], b[j]);
                a[j] = x;
            }
        }
    }
}

template <>
void Radix2<uint32_t>::butterfly_ct_3(
    uint32_t coef,
    vec::Buffers<uint32_t>& buf,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = ALIGN_SIZE / sizeof(uint32_t);
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;
    const uint32_t card = this->gf->card();

    // perform vector operations
    simd::butterfly_ct_3(coef, buf, start, m, step, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        butterfly_ct_3_offset(coef, buf, start, m, step, vec_len * ratio);
    }
}

template <>
void Radix2<uint32_t>::butterfly_gs_2(
    vec::Buffers<uint32_t>& buf,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = ALIGN_SIZE / sizeof(uint32_t);
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;
    const uint32_t card = this->gf->card();

    // perform vector operations
    simd::butterfly_gs_2(buf, start, m, step, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        for (int i = start; i < this->n; i += step) {
            uint32_t* a = buf.get(i);
            uint32_t* b = buf.get(i + m);
            // perform butterfly operation for Cooley-Tukey FFT algorithm
            for (size_t j = vec_len * ratio; j < len; ++j) {
                uint32_t x = this->gf->add(a[j], b[j]);
                b[j] = this->gf->sub(b[j], a[j]);
                a[j] = x;
            }
        }
    }
}

template <>
void Radix2<uint32_t>::butterfly_gs_3(
    uint32_t coef,
    vec::Buffers<uint32_t>& buf,
    unsigned start,
    unsigned m,
    unsigned step)
{
    const unsigned ratio = ALIGN_SIZE / sizeof(uint32_t);
    const size_t len = this->pkt_size;
    const size_t vec_len = len / ratio;
    const size_t last_len = len - vec_len * ratio;
    const uint32_t card = this->gf->card();

    // perform vector operations
    simd::butterfly_gs_3(coef, buf, start, m, step, vec_len, card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        for (int i = start; i < this->n; i += step) {
            uint32_t* a = buf.get(i);
            uint32_t* b = buf.get(i + m);
            // perform butterfly operation for Cooley-Tukey FFT algorithm
            for (size_t j = vec_len * ratio; j < len; ++j) {
                uint32_t x = this->gf->sub(a[j], b[j]);
                a[j] = this->gf->add(a[j], b[j]);
                b[j] = this->gf->mul(coef, x);
            }
        }
    }
}

} // namespace fft
} // namespace quadiron

#endif // #ifdef QUADIRON_USE_SIMD
