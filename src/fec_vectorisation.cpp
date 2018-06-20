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

#include "fec_rs_fnt.h"

/*
 * The file includes vectorized operations used by FEC classes
 */

#ifdef NTTEC_USE_SIMD

#include "simd.h"

namespace nttec {
namespace fec {

template <>
void RsFnt<uint16_t>::encode_post_process(
    vec::Buffers<uint16_t>& output,
    std::vector<Properties>& props,
    off_t offset)
{
    size_t size = this->pkt_size;
    uint16_t threshold = this->gf->card_minus_one();
    unsigned code_len = this->code_len;

    // number of elements per vector register
    unsigned vec_size = ALIGN_SIZE / sizeof(uint16_t);
    // number of vector registers per fragment packet
    size_t vecs_nb = size / vec_size;
    // odd number of elements not vectorized
    size_t last_len = size - vecs_nb * vec_size;

    simd::encode_post_process(
        output, props, offset, code_len, threshold, vecs_nb);

    if (last_len > 0) {
        for (unsigned i = 0; i < code_len; ++i) {
            uint16_t* chunk = output.get(i);
            for (int j = vecs_nb * vec_size; j < size; ++j) {
                if (chunk[j] == threshold) {
                    const ValueLocation loc(offset + j, i);
                    props[i].add(loc, "@");
                }
            }
        }
    }
}

template <>
void RsFnt<uint32_t>::encode_post_process(
    vec::Buffers<uint32_t>& output,
    std::vector<Properties>& props,
    off_t offset)
{
    size_t size = this->pkt_size;
    uint32_t threshold = this->gf->card_minus_one();
    unsigned code_len = this->code_len;

    // number of elements per vector register
    unsigned vec_size = ALIGN_SIZE / sizeof(uint32_t);
    // number of vector registers per fragment packet
    size_t vecs_nb = size / vec_size;
    // odd number of elements not vectorized
    size_t last_len = size - vecs_nb * vec_size;

    simd::encode_post_process(
        output, props, offset, code_len, threshold, vecs_nb);

    if (last_len > 0) {
        for (unsigned i = 0; i < code_len; ++i) {
            uint32_t* chunk = output.get(i);
            for (int j = vecs_nb * vec_size; j < size; ++j) {
                if (chunk[j] == threshold) {
                    const ValueLocation loc(offset + j, i);
                    props[i].add(loc, "@");
                }
            }
        }
    }
}

} // namespace fec
} // namespace nttec

#endif // #ifdef NTTEC_USE_SIMD
