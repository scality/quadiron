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
#ifndef __NTTEC_FEC_RS_LEO_H__
#define __NTTEC_FEC_RS_LEO_H__

// use <leopard.h> instead when linking against installed
#include "leopard.h"

#include "fec_base.h"
#include "gf_bin_ext.h"
#include "vec_matrix.h"
#include "vec_vector.h"

namespace nttec {
namespace fec {

typedef unsigned char u8;
/** Reed-Solomon (RS) Erasure code over GF(2<sup>n</sup>) (Cauchy or
 *  Vandermonde) using the Leopard library
 */
template <typename T>
class RsLeo : public FecCode<T> {
  public:
    RsMatrixType mat_type;

    RsLeo(
        unsigned word_size,
        unsigned n_data,
        unsigned n_parities,
        RsMatrixType type,
        size_t pkt_size = 8 * 1024)
        : FecCode<T>(
              FecType::SYSTEMATIC,
              word_size,
              n_data,
              n_parities,
              pkt_size)
    {
        mat_type = type;
        this->fec_init();
    }

    ~RsLeo() = default;

    inline void check_params() override
    {
        if (this->word_size > 2) {
            assert(false); // not support yet
            exit(1);
        }
        if (this->code_len > 65536) {
            assert(false); // not support yet
            exit(1);
        }
    }

    inline void init_gf() override {}

    inline void init_fft() override {}

    inline void init_others() override
    {
        if (0 != leo_init()) {
            std::cout << "Failed to initialize" << std::endl;
            exit(1);
        }
        encode_work_count =
            leo_encode_work_count(this->n_data, this->n_parities);
        decode_work_count =
            leo_decode_work_count(this->n_data, this->n_parities);

        unsigned enc_data_len = encode_work_count - this->n_parities;
        unsigned dec_data_len = decode_work_count - this->n_data;

        enc_data = std::unique_ptr<vec::Buffers<T>>(
            new vec::Buffers<T>(enc_data_len, this->pkt_size));
        dec_data = std::unique_ptr<vec::Buffers<T>>(
            new vec::Buffers<T>(dec_data_len, this->pkt_size));

        original_data = std::unique_ptr<std::vector<T*>>(
            new std::vector<T*>(this->n_data, nullptr));
        encode_work_data = std::unique_ptr<std::vector<T*>>(
            new std::vector<T*>(encode_work_count, nullptr));
    }

    int get_n_outputs() override
    {
        return this->n_parities;
    }

    void encode(
        vec::Vector<T>& output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>& words) override
    {
    }

    void encode(
        vec::Buffers<T>& output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>& words) override
    {
        vec::Buffers<T> work_data(output, *enc_data);
        std::vector<T*>* vec_work = work_data.get_mem();

        std::vector<T*>* vec_data = words.get_mem();

        LeopardResult encodeResult = leo_encode(
            this->buf_size,
            this->n_data,
            this->n_parities,
            encode_work_count,
            (void**)&vec_data->at(0),
            (void**)&vec_work->at(0));

        if (encodeResult != Leopard_Success) {
            if (encodeResult == Leopard_TooMuchData) {
                std::cout << "Skipping this test: Parameters are unsupported "
                             "by the codec"
                          << std::endl;
                return;
            }
            std::cout << "Error: Leopard encode failed with result="
                      << encodeResult << ": " << leo_result_string(encodeResult)
                      << std::endl;
        }
    }

    void reset_for_new_dec()
    {
        std::fill(original_data->begin(), original_data->end(), nullptr);
        std::fill(encode_work_data->begin(), encode_work_data->end(), nullptr);
    }

    void
    prepare_for_new_dec(vec::Buffers<T>& words)
    {
        for (unsigned i = 0; i < this->n_data; ++i) {
            unsigned frag_id = fragments_ids->get(i);
            if (frag_id < this->n_data) {
                original_data->at(frag_id) = words.get(i);
            } else {
                unsigned parity_id = frag_id - this->n_data;
                encode_work_data->at(parity_id) = words.get(i);
            }
        }
    }

    void decode_add_data(int fragment_index, int row) override
    {
        if (fragment_index == 0) {
            reset_for_new_dec();
        }
    }

    void decode_add_parities(int fragment_index, int row) override
    {
        if (fragment_index == 0) {
            reset_for_new_dec();
        }
    }

    void decode_build() override {}

    std::unique_ptr<DecodeContext<T>>
    init_context_dec(vec::Vector<T>& fragments_ids, size_t size) override
    {
        this->fragments_ids = &fragments_ids;
        std::unique_ptr<DecodeContext<T>> context;
        return context;
    }

    void decode(
        const DecodeContext<T>& context,
        vec::Vector<T>& output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>& words) override
    {
    }

    void decode(
        const DecodeContext<T>& context,
        vec::Buffers<T>& output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>& words) override
    {
        prepare_for_new_dec(words);

        vec::Buffers<T> work_data(output, *dec_data);

        std::vector<T*>* vec_work = work_data.get_mem();

        LeopardResult decodeResult = leo_decode(
            this->buf_size,
            this->n_data,
            this->n_parities,
            decode_work_count,
            (void**)&original_data->at(0),
            (void**)&encode_work_data->at(0),
            (void**)&vec_work->at(0));

        if (decodeResult != Leopard_Success) {
            std::cout << "Error: Leopard decode failed with result="
                      << decodeResult << ": " << leo_result_string(decodeResult)
                      << std::endl;
        }

        // copy received data fragments from words
        for (unsigned i = 0; i < this->n_data; ++i) {
            unsigned frag_id = fragments_ids->get(i);
            if (frag_id < this->n_data) {
                output.copy(frag_id, words.get(i));
            }
        }
    }

  private:
    unsigned encode_work_count;
    unsigned decode_work_count;
    std::unique_ptr<vec::Buffers<T>> enc_data = nullptr;
    std::unique_ptr<vec::Buffers<T>> dec_data = nullptr;

    const vec::Vector<T>* fragments_ids;
    std::unique_ptr<std::vector<T*>> original_data = nullptr;
    std::unique_ptr<std::vector<T*>> encode_work_data = nullptr;
};

} // namespace fec
} // namespace nttec

#endif
