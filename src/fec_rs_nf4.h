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
#include "vec_vector.h"

namespace quadiron {
namespace fec {

/** Reed-Solomon (RS) Erasure code over `n` GF(F<sub>4</sub>). */
template <typename T>
class RsNf4 : public FecCode<T> {
  public:
    RsNf4(
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
        assert(this->word_size >= 2);
        assert(this->word_size <= 8);
    }

    inline void init_gf() override
    {
        gf_n = this->word_size / 2;
        this->gf = gf::alloc<gf::Field<T>, gf::NF4<T>>(gf_n);
        ngff4 = static_cast<gf::NF4<T>*>(this->gf.get());
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

        int m = arith::ceil2<int>(this->n_data);
        this->fft = std::unique_ptr<fft::Radix2<T>>(
            new fft::Radix2<T>(*ngff4, this->n, m, this->pkt_size));

        unsigned len_2k = this->gf->get_code_len_high_compo(2 * this->n_data);
        this->fft_2k = std::unique_ptr<fft::Radix2<T>>(
            new fft::Radix2<T>(*ngff4, len_2k, len_2k, this->pkt_size));
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

        work_buf =
            std::make_unique<vec::Buffers<T>>(this->n_data, this->pkt_size);
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
        // std::cout << "words:"; words.dump();
        for (unsigned i = 0; i < this->n_data; i++) {
            words.set(i, ngff4->pack(words.get(i)));
        }
        // std::cout << "pack words:"; words.dump();
        this->fft->fft(output, words);
        encode_post_process(output, props, offset);
    }

    void encode_post_process(
        vec::Vector<T>& output,
        std::vector<Properties>& props,
        off_t offset) override
    {
        // std::cout << "encoded:"; output.dump();
        GroupedValues<T> true_val;
        for (unsigned i = 0; i < this->code_len; i++) {
            T val = output.get(i);
            ngff4->unpack(val, true_val);
            if (true_val.flag > 0) {
                props[i].add(offset, true_val.flag);
                // std::cout << "\ni:" << true_val.flag << " at buf" << buf <<
                // std::endl; std::cout << "encode: val:" << val << " <- " <<
                // true_val.val << std::endl;
            }
            output.set(i, true_val.values);
        }
        // std::cout << "unpacked:"; output.dump();
    }

    void decode_add_data(int, int) override
    {
        // not applicable
        assert(false);
    }

  private:
    std::unique_ptr<vec::Buffers<T>> work_buf;
    const gf::Field<uint32_t>* sub_field;
    gf::NF4<T>* ngff4;
    int gf_n;

  protected:
    std::unique_ptr<DecodeContext<T>> init_context_dec(
        vec::Vector<T>& fragments_ids,
        std::vector<Properties>& input_props,
        size_t size,
        vec::Buffers<T>* output) override
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
                input_props,
                vx,
                k,
                this->n,
                -1,
                size,
                output));

        return context;
    }

    void decode_prepare(
        DecodeContext<T>& context,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>& words) override
    {
        T true_val;
        const vec::Vector<T>& fragments_ids = context.get_fragments_id();
        int k = this->n_data; // number of fragments received
        for (int i = 0; i < k; ++i) {
            const int j = fragments_ids.get(i);
            if (props[j].is_marked(context.props_indices[j], offset)) {
                true_val = ngff4->pack(
                    words.get(i), props[j].marker(context.props_indices[j]));
                context.props_indices.at(j)++;
            } else {
                true_val = ngff4->pack(words.get(i));
            }
            words.set(i, true_val);
        }
    }

    void decode_apply(
        DecodeContext<T>& context,
        vec::Vector<T>& output,
        vec::Vector<T>& words) override
    {
        // decode_apply: do the same thing as in fec_base
        FecCode<T>::decode_apply(context, output, words);
        // unpack decoded symbols
        for (unsigned i = 0; i < this->n_data; ++i) {
            output.set(i, ngff4->unpack(output.get(i)).values);
        }
    }

    /********** Encoding & Decoding using Buffers **********/

    void encode(
        vec::Buffers<T>& output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>& words) override
    {
        // as Buffers has not meta, `words` contains only `buf_size = pkt_size *
        // word_size` bytes from source. It is stored in first `buf_size /
        // sizeof(T)` elements of `words`

        const unsigned nb_words_per_element = sizeof(T) / this->word_size;
        for (unsigned i = 0; i < this->n_data; ++i) {
            T* chunk = words.get(i);
            T* work = work_buf->get(i);

            size_t u = 0;
            T element = chunk[u];
            for (size_t j = 0, u = 0; j < this->pkt_size; ++j) {
                work[j] = ngff4->pack(element);

                (j + 1) % nb_words_per_element == 0
                    ? element = chunk[++u]
                    : element = static_cast<T>(element)
                                >> (CHAR_BIT * this->word_size);
            }
        }

        this->fft->fft(output, *work_buf);
        encode_post_process(output, props, offset);
    }

    void encode_post_process(
        vec::Buffers<T>& output,
        std::vector<Properties>& props,
        off_t offset) override
    {
        // as Buffers has not meta, output write only `buf_size = pkt_size *
        // word_size` bytes to destination
        // This data should be stored in first `buf_size / sizeof(T)` elements.

        const unsigned nb_words_per_element = sizeof(T) / this->word_size;
        const size_t size = output.get_size();
        GroupedValues<T> true_val;
        for (unsigned frag_id = 0; frag_id < this->code_len; ++frag_id) {
            T* chunk = output.get(frag_id);

            size_t out_symb_id = 0;
            unsigned symb_offset = 0;
            T element = 0;
            for (size_t symb_id = 0; symb_id < size; ++symb_id) {
                ngff4->unpack(chunk[symb_id], true_val);
                if (true_val.flag > 0) {
                    const off_t loc = offset + symb_id;
                    props[frag_id].add(loc, true_val.flag);
                }

                // for test
                chunk[symb_id] = 0;

                element |= (static_cast<T>(true_val.values) << symb_offset);
                if ((symb_id + 1) % nb_words_per_element == 0) {
                    chunk[out_symb_id] = element;
                    out_symb_id++;
                    symb_offset = 0;
                    element = 0;
                } else {
                    symb_offset += (CHAR_BIT * this->word_size);
                }
            }
            // for the last no-full element
            if (symb_offset > 0) {
                chunk[out_symb_id] = element;
            }
        }
    }

    void decode_prepare(
        DecodeContext<T>& context,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>& words) override
    {
        const vec::Vector<T>& fragments_ids = context.get_fragments_id();
        // as Buffers has not meta, `words` contains only `buf_size = pkt_size *
        // word_size` bytes from source. It is stored in first `buf_size /
        // sizeof(T)` elements of `words`

        const unsigned nb_words_per_element = sizeof(T) / this->word_size;
        for (unsigned i = 0; i < this->n_data; ++i) {
            const int frag_id = fragments_ids.get(i);
            T* chunk = words.get(i);
            T* work = work_buf->get(i);

            size_t u = 0;
            T element = chunk[u];
            for (size_t j = 0, u = 0; j < this->pkt_size; ++j) {
                if (props[frag_id].is_marked(
                        context.props_indices[frag_id], offset + j)) {
                    // pack marked symbol
                    work[j] = ngff4->pack(
                        element,
                        props[frag_id].marker(context.props_indices[frag_id]));
                    context.props_indices.at(frag_id)++;
                } else {
                    // pack un-marked symbol
                    work[j] = ngff4->pack(element);
                }

                (j + 1) % nb_words_per_element == 0
                    ? element = chunk[++u]
                    : element = static_cast<T>(element)
                                >> (CHAR_BIT * this->word_size);
            }
        }
    }

    void decode_apply(
        DecodeContext<T>& context,
        vec::Buffers<T>& output,
        vec::Buffers<T>&) override
    {
        // as Buffers has not meta, output write only `buf_size = pkt_size *
        // word_size` bytes to destination
        // This data should be stored in first `buf_size / sizeof(T)` elements.

        // decode_apply: do the same thing as in fec_base
        FecCode<T>::decode_apply(context, output, *work_buf);

        const unsigned nb_words_per_element = sizeof(T) / this->word_size;
        GroupedValues<T> true_val;
        // unpack decoded symbols
        for (unsigned frag_id = 0; frag_id < this->n_data; ++frag_id) {
            T* chunk = output.get(frag_id);

            size_t out_symb_id = 0;
            unsigned symb_offset = 0;
            T element = 0;
            for (size_t symb_id = 0; symb_id < this->pkt_size; ++symb_id) {
                ngff4->unpack(chunk[symb_id], true_val);
                // for test
                chunk[symb_id] = 0;

                element |= (static_cast<T>(true_val.values) << symb_offset);
                if ((symb_id + 1) % nb_words_per_element == 0) {
                    chunk[out_symb_id] = element;
                    out_symb_id++;
                    symb_offset = 0;
                    element = 0;
                } else {
                    symb_offset += (CHAR_BIT * this->word_size);
                }
            }
            // for the last no-full element
            if (symb_offset > 0) {
                chunk[out_symb_id] = element;
            }
        }
    }
};

} // namespace fec
} // namespace quadiron

#endif
