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
#ifndef __NTTEC_FEC_RS_GFP_FFT_H__
#define __NTTEC_FEC_RS_GFP_FFT_H__

#include "arith.h"
#include "fec_base.h"
#include "fft_2n.h"
#include "fft_base.h"
#include "fft_ct.h"
#include "gf_prime.h"
#include "polynomial.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

namespace nttec {
namespace fec {

/** Reed-Solomon (RS) Erasure code over prime Galois Fields and FFT.
 *
 * @note Necessary condition: p ≥ 2<sup>8*word_size</sup> to cover all
 * possible values of symbols, as each symbol is 8*word_size bits.
 *
 * A length-k vector X is encoded into a length-n vector Y whose elements are
 * GF(p) elements. Hence, elements of Y could be greater than
 * 2<sup>8*word_size</sup>
 *
 * We can deal with this as below.
 *
 * @note To facilitate our implementation, we choose
 *   p < 2 * 2<sup>8*word_size</sup>.
 *
 * Supposing that an element Y[i] > 2<sup>8*word_size</sup>, we will store its
 * value modulo 2<sup>8*word_size</sup>, i.e.
 *
 *   Y<sub>p</sub>[i] = Y[i] % 2<sup>8*word_size</sup>
 *
 * A flag will be used to mark the position `i`.
 *
 * To decode X, we read Y[i] from Y<sub>p</sub>[i] as:
 *
 *    Y[i] = Y<sub>p</sub>[i] + 2<sup>8*word_size</sup>
 *
 * Because p < 2 * 2<sup>8*word_size</sup>, a single bool is enough as flag.
 */
template <typename T>
class RsGfpFft : public FecCode<T> {
  public:
    RsGfpFft(unsigned word_size, unsigned n_data, unsigned n_parities)
        : FecCode<T>(FecType::NON_SYSTEMATIC, word_size, n_data, n_parities)
    {
        // warning all fermat numbers >= to F_5 (2^32+1) are composite!!!
        T gf_p = 0;
        if (word_size < 4) {
            gf_p = (1ULL << (8 * word_size)) + 1;
            this->limit_value = (1ULL << (8 * word_size));
        } else if (word_size == 4) {
            gf_p = (T)4294991873ULL;              // p-1=2^13 29^1 101^1 179^1
            this->limit_value = (T)4294967296ULL; // 2^32
        } else {
            assert(false); // not support yet
            exit(1);
        }

        assert(gf_p >= this->limit_value);
        // we choose gf_p for a simple implementation
        assert(gf_p / 2 < this->limit_value);

        this->gf = new gf::Prime<T>(gf_p);
        assert(
            arith::jacobi<T>(this->gf->get_primitive_root(), this->gf->card())
            == -1);

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        this->n = this->gf->get_code_len_high_compo(n_parities + n_data);

        // compute root of order n-1 such as r^(n-1) mod q == 1
        this->r = this->gf->get_nth_root(this->n);

        // std::cerr << "limit_value=" << limit_value << "\n";
        // std::cerr << "gf_p=" << gf_p << "\n";
        // std::cerr << "n=" << n << "\n";
        // std::cerr << "r=" << r << "\n";

        if (arith::is_power_of_2<T>(this->n)) {
            this->fft = new fft::Radix2<T>(this->gf, this->n);
            this->fft_full = std::unique_ptr<fft::Radix2<T>>(
                new fft::Radix2<T>(this->gf, this->n));
        } else {
            this->fft = new fft::CooleyTukey<T>(this->gf, this->n);
            this->fft_full = std::unique_ptr<fft::CooleyTukey<T>>(
                new fft::CooleyTukey<T>(this->gf, this->n));
        }

        // vector stores r^{-i} for i = 0, ... , k
        T inv_r = this->gf->inv(this->r);
        this->inv_r_powers = std::unique_ptr<vec::Vector<T>>(
            new vec::Vector<T>(this->gf, this->n_data + 1));
        for (int i = 0; i <= this->n_data; i++)
            this->inv_r_powers->set(i, this->gf->exp(inv_r, i));
    }

    ~RsGfpFft()
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
     * @param props must be exactly n
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
        // check for out of range value in output
        for (unsigned i = 0; i < this->code_len; i++) {
            if (output->get(i) >= this->limit_value) {
                props[i].add(ValueLocation(offset, i), "@");
                output->set(i, output->get(i) % this->limit_value);
            }
        }
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

  private:
    // fft::FourierTransform<T>* fft = nullptr;
    T limit_value;

  protected:
    void decode_prepare(
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words,
        vec::Vector<T>* vx,
        int* vx_zero)
    {
        int k = this->n_data; // number of fragments received
        // vector x=(x_0, x_1, ..., x_k-1)
        for (int i = 0; i < k; i++) {
            vx->set(i, this->gf->exp(this->r, fragments_ids->get(i)));
        }

        for (int i = 0; i < k; i++) {
            const int j = fragments_ids->get(i);
            auto data = props[j].get(ValueLocation(offset, j));

            // Check if the symbol is a special case whick is marked by "@".
            // In encoded data, its value was subtracted by the predefined
            // limite_value. This operation restore its value.
            if (data && *data == "@") {
                words->set(i, words->get(i) + limit_value);
            }
        }
    }
};

} // namespace fec
} // namespace nttec

#endif
