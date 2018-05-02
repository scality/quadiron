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
#ifndef __NTTEC_FEC_RS_ISAL_H__
#define __NTTEC_FEC_RS_ISAL_H__

// use <isa-l.h> instead when linking against installed
#include "isa-l/include/erasure_code.h"

#include "fec_base.h"
#include "gf_bin_ext.h"
#include "vec_matrix.h"
#include "vec_vector.h"

namespace nttec {
namespace fec {

// FIXME: Does ISA-L support only erasure codes on GF(2^8)?
#define MMAX 255
#define KMAX 255

typedef unsigned char u8;
/** Reed-Solomon (RS) Erasure code over GF(2<sup>n</sup>) (Cauchy or
 *  Vandermonde) using the ISA-L library
 */
template <typename T>
class RsIsal : public FecCode<T> {
  public:
    RsMatrixType mat_type;

    RsIsal(
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

    ~RsIsal() = default;

    inline void check_params() override
    {
        assert(
            mat_type == RsMatrixType::VANDERMONDE
            || mat_type == RsMatrixType::CAUCHY);

        // FIXME: Does ISA-L support only erasure codes on GF(2^8)?
        if (this->word_size > 1 || sizeof(T) > 1) {
            assert(false); // not support yet
            exit(1);
        }
    }

    inline void init_gf() override {}

    inline void init_fft() override {}

    inline void init_others() override
    {
        encode_mat = std::unique_ptr<u8>(new u8[this->code_len * this->n_data]);
        decode_mat_tmp =
            std::unique_ptr<u8>(new u8[this->n_data * this->n_data]);
        decode_mat = std::unique_ptr<u8>(new u8[this->n_data * this->n_data]);
        enc_g_tbls =
            std::unique_ptr<u8>(new u8[this->n_data * this->n_parities * 32]);
        dec_g_tbls =
            std::unique_ptr<u8>(new u8[this->n_data * this->n_data * 32]);

        // Pick an encode matrix.
        if (mat_type == RsMatrixType::CAUCHY) {
            gf_gen_cauchy1_matrix(
                encode_mat.get(), this->code_len, this->n_data);
        } else if (mat_type == RsMatrixType::VANDERMONDE) {
            gf_gen_rs_matrix(encode_mat.get(), this->code_len, this->n_data);
        }
        // Initialize enc_g_tbls from encode matrix
        ec_init_tables(
            this->n_data,
            this->n_parities,
            encode_mat.get() + this->n_data * this->n_data,
            enc_g_tbls.get());
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

    void do_multiply(vec::Buffers<T>& output, vec::Buffers<T>& words, u8* tbls)
    {
        const unsigned input_len = words.get_n();
        const unsigned output_len = output.get_n();
        std::vector<T*>* vec_data = words.get_mem();
        std::vector<T*>* vec_coding = output.get_mem();

        u8** data = reinterpret_cast<u8**>(vec_data->data());
        u8** coding = reinterpret_cast<u8**>(vec_coding->data());
        // Generate EC parity blocks from sources

#if NTTEC_USE_AVX2
        ec_encode_data_avx2(
            this->pkt_size, input_len, output_len, tbls, data, coding);
#elif NTTEC_USE_SSE4
        ec_encode_data_sse(
            this->pkt_size, input_len, output_len, tbls, data, coding);
#else
        ec_encode_data(
            this->pkt_size, input_len, output_len, tbls, data, coding);
#endif
    }

    void encode(
        vec::Buffers<T>& output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>& words) override
    {
        do_multiply(output, words, enc_g_tbls.get());
    }

    void decode_add_data(int fragment_index, int row) override
    {
        std::copy_n(
            encode_mat.get() + row * this->n_data,
            this->n_data,
            decode_mat_tmp.get() + fragment_index * this->n_data);
    }

    void decode_add_parities(int fragment_index, int row) override
    {
        int _row = row + this->n_data;
        std::copy_n(
            encode_mat.get() + _row * this->n_data,
            this->n_data,
            decode_mat_tmp.get() + fragment_index * this->n_data);
    }

    void decode_build() override
    {
        gf_invert_matrix(decode_mat_tmp.get(), decode_mat.get(), this->n_data);
        // Initialize dec_g_tbls from decode matrix
        ec_init_tables(
            this->n_data, this->n_data, decode_mat.get(), dec_g_tbls.get());
    }

    std::unique_ptr<DecodeContext<T>>
    init_context_dec(vec::Vector<T>& fragments_ids, size_t size) override
    {
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
        do_multiply(output, words, dec_g_tbls.get());
    }

  private:
    std::unique_ptr<u8> encode_mat = nullptr;
    std::unique_ptr<u8> decode_mat = nullptr;
    std::unique_ptr<u8> decode_mat_tmp = nullptr;
    std::unique_ptr<u8> enc_g_tbls = nullptr;
    std::unique_ptr<u8> dec_g_tbls = nullptr;
};

} // namespace fec
} // namespace nttec

#endif
