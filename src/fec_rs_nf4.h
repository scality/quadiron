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
#ifndef __QUAD_FEC_RS_NF4_H__
#define __QUAD_FEC_RS_NF4_H__

#include <string>

#include "fec_base.h"
#include "fft_2n.h"
#include "gf_base.h"
#include "gf_nf4.h"
#include "polynomial.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

namespace quadiron {
namespace fec {

/** Reed-Solomon (RS) Erasure code over `n` GF(F<sub>4</sub>). */
template <typename T>
class RsNf4 : public FecCode<T> {
  public:
    RsNf4(unsigned word_size, unsigned n_data, unsigned n_parities)
        : FecCode<T>(FecType::NON_SYSTEMATIC, word_size, n_data, n_parities)
    {
        this->fec_init();
    }

    inline void check_params() override
    {
        assert(this->word_size >= 2);
        assert(this->word_size <= 8);
    }

    inline void init_gf() override
    {
        gf_n = this->word_size / 2;

        // NOTE: ngff4 is wrapped in this->gf that will release it
        // Hence don't release ngff4 manually
        ngff4 = new gf::NF4<T>(gf_n);
        this->gf = std::unique_ptr<gf::Field<T>>(ngff4);

        sub_field = &(ngff4->get_sub_field());
    }

    inline void init_fft() override
    {
        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        this->n =
            sub_field->get_code_len_high_compo(this->n_parities + this->n_data);

        // compute root of order n-1 such as r^(n-1) mod q == (1, ..,1)
        this->r = ngff4->get_nth_root(this->n);

        int m = arith::get_smallest_power_of_2<int>(this->n_data);
        this->fft = std::unique_ptr<fft::Radix2<T>>(
            new fft::Radix2<T>(*ngff4, this->n, m));

        this->fft_full = std::unique_ptr<fft::Radix2<T>>(
            new fft::Radix2<T>(*ngff4, this->n));

        unsigned len_2k = this->gf->get_code_len_high_compo(2 * this->n_data);
        this->fft_2k = std::unique_ptr<fft::Radix2<T>>(
            new fft::Radix2<T>(*ngff4, len_2k, len_2k));
    }

    inline void init_others() override
    {
        // vector stores r^{-i} for i = 0, ... , k
        const T inv_r = ngff4->inv(this->r);
        this->inv_r_powers = std::unique_ptr<vec::Vector<T>>(
            new vec::Vector<T>(*ngff4, this->n_data + 1));
        for (unsigned i = 0; i <= this->n_data; i++)
            this->inv_r_powers->set(i, ngff4->exp(inv_r, i));

        // vector stores r^{i} for i = 0, ... , k
        this->r_powers = std::unique_ptr<vec::Vector<T>>(
            new vec::Vector<T>(*ngff4, this->n));
        for (unsigned i = 0; i < this->n; i++) {
            this->r_powers->set(i, ngff4->exp(this->r, i));
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
        vec::Vector<T>* output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* words) override
    {
        // std::cout << "words:"; words->dump();
        for (unsigned i = 0; i < this->n_data; i++) {
            words->set(i, ngff4->pack(words->get(i)));
        }
        // std::cout << "pack words:"; words->dump();
        vec::ZeroExtended<T> vwords(words, this->n);
        this->fft->fft(output, &vwords);
        // std::cout << "encoded:"; output->dump();
        for (unsigned i = 0; i < this->code_len; i++) {
            T val = output->get(i);
            GroupedValues<T> true_val = ngff4->unpack(val);
            if (true_val.flag > 0) {
                props[i].add(
                    ValueLocation(offset, i), std::to_string(true_val.flag));
                // std::cout << "\ni:" << true_val.flag << " at buf" << buf <<
                // std::endl; std::cout << "encode: val:" << val << " <- " <<
                // true_val.val << std::endl;
            }
            output->set(i, true_val.values);
        }
        // std::cout << "unpacked:"; output->dump();
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

  private:
    const gf::Field<uint32_t>* sub_field;
    gf::NF4<T>* ngff4;
    int gf_n;

  protected:
    std::unique_ptr<DecodeContext<T>>
    init_context_dec(vec::Vector<T>& fragments_ids, size_t size) override
    {
        if (this->inv_r_powers == nullptr) {
            throw LogicError("FEC base: vector (inv_r)^i must be initialized");
        }
        if (this->r_powers == nullptr) {
            throw LogicError("FEC base: vector r^i must be initialized");
        }
        if (this->fft == nullptr) {
            throw LogicError("FEC base: FFT must be initialized");
        }
        if (this->fft_full == nullptr) {
            throw LogicError("FEC base: FFT full must be initialized");
        }

        int k = this->n_data; // number of fragments received
        // vector x=(x_0, x_1, ..., x_k-1)
        vec::Vector<T> vx(*(this->gf), k);
        for (int i = 0; i < k; ++i) {
            vx.set(
                i,
                this->gf->exp(this->r, ngff4->replicate(fragments_ids.get(i))));
        }

        std::unique_ptr<DecodeContext<T>> context =
            std::unique_ptr<DecodeContext<T>>(new DecodeContext<T>(
                *(this->gf),
                *(this->fft),
                *(this->fft_2k),
                fragments_ids,
                vx,
                k,
                this->n,
                -1,
                size));

        return context;
    }

    void decode_prepare(
        const DecodeContext<T>& context,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* words) override
    {
        T true_val;
        const vec::Vector<T>& fragments_ids = context.get_fragments_id();
        // std::cout << "fragments_ids:"; fragments_ids->dump();
        int k = this->n_data; // number of fragments received
        for (int i = 0; i < k; ++i) {
            const int j = fragments_ids.get(i);
            auto data = props[j].get(ValueLocation(offset, j));

            if (data) {
                uint32_t flag = std::stoul(*data);
                true_val = ngff4->pack(words->get(i), flag);
            } else {
                true_val = ngff4->pack(words->get(i));
            }
            words->set(i, true_val);
        }
    }

    void decode_apply(
        const DecodeContext<T>& context,
        vec::Vector<T>* output,
        vec::Vector<T>* words) override
    {
        // decode_apply: do the same thing as in fec_base
        FecCode<T>::decode_apply(context, output, words);
        // unpack decoded symbols
        for (unsigned i = 0; i < this->n_data; ++i) {
            output->set(i, ngff4->unpack(output->get(i)).values);
        }
    }
};

} // namespace fec
} // namespace quadiron

#endif
