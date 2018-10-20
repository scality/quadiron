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
#ifndef __QUAD_FEC_RS_GF2N_H__
#define __QUAD_FEC_RS_GF2N_H__

#include "fec_base.h"
#include "gf_bin_ext.h"
#include "vec_matrix.h"
#include "vec_vector.h"

namespace quadiron {
namespace fec {

enum class RsMatrixType { VANDERMONDE, CAUCHY };

/** Reed-Solomon (RS) Erasure code over GF(2<sup>n</sup>) (Cauchy or
 *  Vandermonde).
 */
template <typename T>
class RsGf2n : public FecCode<T> {
  public:
    RsMatrixType mat_type;

    RsGf2n(
        unsigned word_size,
        unsigned n_data,
        unsigned n_parities,
        RsMatrixType type)
        : FecCode<T>(FecType::SYSTEMATIC, word_size, n_data, n_parities)
    {
        mat_type = type;
        this->fec_init();
    }

    inline void check_params() override
    {
        assert(
            mat_type == RsMatrixType::VANDERMONDE
            || mat_type == RsMatrixType::CAUCHY);

        if (this->word_size > 16)
            assert(false); // not support yet
    }

    inline void init_gf() override
    {
        unsigned gf_n = 8 * this->word_size;
        this->gf = gf::alloc<gf::Field<T>, gf::BinExtension<T>>(gf_n);
    }

    inline void init_fft() override {}

    inline void init_others() override
    {
        this->mat = std::unique_ptr<vec::Matrix<T>>(
            new vec::Matrix<T>(*(this->gf), this->n_parities, this->n_data));
        if (mat_type == RsMatrixType::CAUCHY) {
            mat->cauchy();
        } else if (mat_type == RsMatrixType::VANDERMONDE) {
            mat->vandermonde_suitable_for_ec();
        }

        // has to be a n_data*n_data invertible square matrix
        decode_mat = std::unique_ptr<vec::Matrix<T>>(new vec::Matrix<T>(
            *(this->gf), mat->get_n_cols(), mat->get_n_cols()));
    }

    int get_n_outputs() override
    {
        return this->n_parities;
    }

    void encode(
        vec::Vector<T>& output,
        std::vector<Properties>&,
        off_t,
        vec::Vector<T>& words) override
    {
        mat->mul(&output, &words);
    }

    void decode_add_data(int fragment_index, int row) override
    {
        // for each data available generate the corresponding identity
        for (int j = 0; j < mat->get_n_cols(); j++) {
            if (row == j)
                decode_mat->set(fragment_index, j, 1);
            else
                decode_mat->set(fragment_index, j, 0);
        }
    }

    void decode_add_parities(int fragment_index, int row) override
    {
        // copy corresponding row in vandermonde matrix
        for (int j = 0; j < mat->get_n_cols(); j++) {
            decode_mat->set(fragment_index, j, mat->get(row, j));
        }
    }

    void decode_build() override
    {
        decode_mat->inv();
    }

    void decode(
        const DecodeContext<T>&,
        vec::Vector<T>& output,
        const std::vector<Properties>&,
        off_t,
        vec::Vector<T>& words) override
    {
        decode_mat->mul(&output, &words);
    }

    std::unique_ptr<DecodeContext<T>>
    init_context_dec(vec::Vector<T>&, size_t, vec::Buffers<T>*) override
    {
        std::unique_ptr<DecodeContext<T>> context;
        return context;
    }

  private:
    std::unique_ptr<vec::Matrix<T>> mat = nullptr;
    std::unique_ptr<vec::Matrix<T>> decode_mat = nullptr;
};

} // namespace fec
} // namespace quadiron

#endif
