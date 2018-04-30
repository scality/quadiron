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

    ~RsIsal()
    {
        if (mat)
            delete[] mat;
        if (decode_mat)
            delete[] decode_mat;
        if (g_tbls)
            delete[] g_tbls;
    }

    inline void check_params()
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

    inline void init_gf() {}

    inline void init_fft() {}

    inline void init_others()
    {
        mat = new u8[this->code_len * this->n_data];
        decode_mat = new u8[this->n_data * this->n_data];
        g_tbls = new u8[this->n_data * this->n_parities * 32];

        // Pick an encode matrix.
        if (mat_type == RsMatrixType::CAUCHY) {
            gf_gen_cauchy1_matrix(mat, this->code_len, this->n_data);
        } else if (mat_type == RsMatrixType::VANDERMONDE) {
            gf_gen_rs_matrix(mat, this->code_len, this->n_data);
        }
        // Initialize g_tbls from encode matrix
        ec_init_tables(
            this->n_data,
            this->n_parities,
            &mat[this->n_data * this->n_data],
            g_tbls);
    }

    int get_n_outputs()
    {
        return this->n_parities;
    }

    void encode(
        vec::Vector<T>* output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* words)
    {
    }

    void encode(
        vec::Buffers<T>* output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>* words)
    {
        // printf(
        //     " encode (m,k,p)=(%d,%d,%d) len=%zu\n",
        //     this->code_len,
        //     this->n_data,
        //     this->n_parities,
        //     this->pkt_size);

        std::vector<T*>* vec_data = words->get_mem();
        std::vector<T*>* vec_coding = output->get_mem();

        u8** data = reinterpret_cast<u8**>(vec_data->data());
        u8** coding = reinterpret_cast<u8**>(vec_coding->data());
        // Generate EC parity blocks from sources

#if NTTEC_USE_SIMD == AVX2
        ec_encode_data_avx2(
            this->pkt_size,
            this->n_data,
            this->n_parities,
            g_tbls,
            data,
            coding);
#elif NTTEC_USE_SIMD == SSE4
        ec_encode_data_sse(
            this->pkt_size,
            this->n_data,
            this->n_parities,
            g_tbls,
            data,
            coding);
#else
        ec_encode_data(
            this->pkt_size,
            this->n_data,
            this->n_parities,
            g_tbls,
            data,
            coding);
#endif
    }

    void decode_add_data(int fragment_index, int row)
    {
        // for each data available generate the corresponding identity
        for (int j = 0; j < this->n_data; j++) {
            if (row == j)
                // decode_mat->set(fragment_index, j, 1);
                decode_mat[fragment_index * this->code_len + j] = 1;
            else
                // decode_mat->set(fragment_index, j, 0);
                decode_mat[fragment_index * this->code_len + j] = 0;
        }
    }

    void decode_add_parities(int fragment_index, int row)
    {
        // copy corresponding row in vandermonde matrix
        for (int j = 0; j < this->n_data; j++) {
            // decode_mat->set(fragment_index, j, mat->get(row, j));
            decode_mat[fragment_index * this->code_len + j] =
                mat[row * this->code_len + j];
        }
    }

    void decode_build()
    {
        gf_invert_matrix(decode_mat, decode_mat, this->n_data);
        // ec_init_tables(, nerrs, decode_mat, g_tbls);
    }

    void decode(
        ContextDec<T>* context,
        vec::Vector<T>* output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* words)
    {
        // decode_mat->mul(output, words);
    }

    ContextDec<T>* init_context_dec(vec::Vector<T>* fragments_ids)
    {
        return nullptr;
    }

  private:
    u8* mat = nullptr;
    u8* decode_mat = nullptr;
    u8* g_tbls = nullptr;
};

} // namespace fec
} // namespace nttec

#endif
