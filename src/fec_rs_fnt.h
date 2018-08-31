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
#ifndef __QUAD_FEC_RS_FNT_H__
#define __QUAD_FEC_RS_FNT_H__

#include "arith.h"
#include "fec_base.h"
#include "fft_2n.h"
#include "gf_prime.h"
#include "vec_buf_zero_ext.h"
#include "vec_buffers.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

namespace quadiron {
namespace fec {

/** Reed-Solomon (RS) erasure code based on Fermat Number Transform (FNT).
 *
 * This class implements a Forward Error Correction (FEC) code based on FNT
 * (original idea from @cite fnt-rs).
 *
 * It works on \f$GF(2^{2^k} + 1)\f$
 */
template <typename T>
class RsFnt : public FecCode<T> {
  public:
    RsFnt(
        unsigned word_size,
        unsigned n_data,
        unsigned n_parities,
        size_t pkt_size = 8)
        : FecCode<T>(
              FecType::NON_SYSTEMATIC,
              word_size,
              n_data,
              n_parities,
              pkt_size)
    {
        this->fec_init();
    }

    inline void check_params() override
    {
        assert(this->word_size < 4);
    }

    inline void init_gf() override
    {
        // warning all fermat numbers >= to F_5 (2^32+1) are composite!!!
        T gf_p = (1ULL << (8 * this->word_size)) + 1;
        this->gf = gf::alloc<gf::Field<T>, gf::Prime<T>>(gf_p);

        assert(
            arith::jacobi<T>(this->gf->get_primitive_root(), this->gf->card())
            == -1);
    }

    inline void init_fft() override
    {
        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        this->n =
            this->gf->get_code_len_high_compo(this->n_parities + this->n_data);

        // compute root of order n-1 such as r^(n-1) mod q == 1
        this->r = this->gf->get_nth_root(this->n);

        int m = arith::get_smallest_power_of_2<int>(this->n_data);
        this->fft = std::unique_ptr<fft::Radix2<T>>(
            new fft::Radix2<T>(*(this->gf), this->n, m, this->pkt_size));

        unsigned len_2k = this->gf->get_code_len_high_compo(2 * this->n_data);
        this->fft_2k = std::unique_ptr<fft::Radix2<T>>(
            new fft::Radix2<T>(*(this->gf), len_2k, len_2k, this->pkt_size));
    }

    inline void init_others() override
    {
        // vector stores r^{-i} for i = 0, ... , k
        T inv_r = this->gf->inv(this->r);
        this->inv_r_powers = std::unique_ptr<vec::Vector<T>>(
            new vec::Vector<T>(*(this->gf), this->n_data + 1));
        for (unsigned i = 0; i <= this->n_data; i++)
            this->inv_r_powers->set(i, this->gf->exp(inv_r, i));

        // vector stores r^{i} for i = 0, ... , n-1
        this->r_powers = std::unique_ptr<vec::Vector<T>>(
            new vec::Vector<T>(*(this->gf), this->n));
        for (unsigned i = 0; i < this->n; i++) {
            this->r_powers->set(i, this->gf->exp(this->r, i));
        }
    }

    int get_n_outputs() override
    {
        return this->n;
    }

    /**
     * Encode vector
     *
     * @param output must be n
     * @param props must be exactly n
     * @param offset used to locate special values
     * @param words must be n_data
     */
    void encode(
        vec::Vector<T>& output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>& words) override
    {
        this->fft->fft(output, words);
        encode_post_process(output, props, offset);
    }

    void encode_post_process(
        vec::Vector<T>& output,
        std::vector<Properties>& props,
        off_t offset) override
    {
        // max_value = 2^x
        T thres = this->gf->card() - 1;
        // check for out of range value in output
        for (unsigned i = 0; i < this->code_len; i++) {
            if (output.get(i) & thres) {
                props[i].add(ValueLocation(offset, i), "@");
                output.set(i, 0);
            }
        }
    }

    void encode(
        vec::Buffers<T>& output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>& words) override
    {
        this->fft->fft(output, words);
        encode_post_process(output, props, offset);
    }

    void encode_post_process(
        vec::Buffers<T>& output,
        std::vector<Properties>& props,
        off_t offset) override
    {
        // check for out of range value in output
        unsigned size = output.get_size();
        T thres = (this->gf->card() - 1);
        for (unsigned i = 0; i < this->code_len; ++i) {
            T* chunk = output.get(i);
            for (unsigned j = 0; j < size; ++j) {
                if (chunk[j] & thres) {
                    const ValueLocation loc(offset + j, i);

                    props[i].add(loc, "@");
                }
            }
        }
    }

    void decode_add_data(int fragment_index, int row) override
    {
        // not applicable
        assert(false);
    }

    void decode_add_parities(int fragment_index, int row) override
    {
        // we can't anticipate here
    }

    void decode_build() override
    {
        // nothing to do
    }
};

#ifdef QUADIRON_USE_SIMD

/* Operations are vectorized by SIMD */

template <>
void RsFnt<uint16_t>::encode_post_process(
    vec::Buffers<uint16_t>& output,
    std::vector<Properties>& props,
    off_t offset);

template <>
void RsFnt<uint32_t>::encode_post_process(
    vec::Buffers<uint32_t>& output,
    std::vector<Properties>& props,
    off_t offset);

#endif // #ifdef QUADIRON_USE_SIMD

} // namespace fec
} // namespace quadiron

#endif
