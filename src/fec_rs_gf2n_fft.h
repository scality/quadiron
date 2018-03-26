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
#ifndef __NTTEC_FEC_RS_GF2N_FFT_H__
#define __NTTEC_FEC_RS_GF2N_FFT_H__

#include "fec_base.h"
#include "fft_ct.h"
#include "gf_bin_ext.h"
#include "polynomial.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

namespace nttec {
namespace fec {

/** Reed-Solomon (RS) Erasure code over GF(2<sup>n</sup>)using FFT. */
template <typename T>
class RsGf2nFft : public FecCode<T> {
  public:
    // NOTE: only NON_SYSTEMATIC is supported now
    RsGf2nFft(unsigned word_size, unsigned n_data, unsigned n_parities)
        : FecCode<T>(FecType::NON_SYSTEMATIC, word_size, n_data, n_parities)
    {
        if (word_size > 16)
            assert(false); // not support yet
        unsigned gf_n = 8 * word_size;
        this->gf = new gf::BinExtension<T>(gf_n);

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        this->n = this->gf->get_code_len_high_compo(n_parities + n_data);

        // compute root of order n such as r^n == 1
        this->r = this->gf->get_nth_root(this->n);

        // std::cerr << "n_parities=" << n_parities << "\n";
        // std::cerr << "n_data=" << n_data << "\n";
        // std::cerr << "n=" << n << "\n";
        // std::cerr << "r=" << r << "\n";

        this->fft = new fft::CooleyTukey<T>(this->gf, this->n);

        this->fft_full = std::unique_ptr<fft::CooleyTukey<T>>(
            new fft::CooleyTukey<T>(this->gf, this->n));

        // vector stores r^{-i} for i = 0, ... , k
        T inv_r = this->gf->inv(this->r);
        this->inv_r_powers = std::unique_ptr<vec::Vector<T>>(
            new vec::Vector<T>(this->gf, this->n_data + 1));
        for (int i = 0; i <= this->n_data; i++)
            this->inv_r_powers->set(i, this->gf->exp(inv_r, i));
    }

    ~RsGf2nFft()
    {
        delete this->fft;
        delete this->gf;
    }

    int get_n_outputs()
    {
        return this->n;
    }

    /**
     * Encode vector
     *
     * @param output must be n
     * @param props special values dictionary must be exactly n_data
     * @param offset used to locate special values
     * @param words must be n_data
     */
    void encode(
        vec::Vector<T>* output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* words)
    {
        vec::ZeroExtended<T> vwords(words, this->n);
        this->fft->fft(output, &vwords);
    }

    void decode_add_data(int fragment_index, int row)
    {
        // not applicable
        assert(false);
    }

    void decode_add_parities(int fragment_index, int row)
    {
        // we can't anticipate here
    }

    void decode_build()
    {
        // nothing to do
    }
};

} // namespace fec
} // namespace nttec

#endif
