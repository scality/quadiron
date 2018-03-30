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
#ifndef __NTTEC_FEC_RS_GF2N_FFT_ADD_H__
#define __NTTEC_FEC_RS_GF2N_FFT_ADD_H__

#include "arith.h"
#include "fec_base.h"
#include "fft_add.h"
#include "gf_bin_ext.h"
#include "polynomial.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

namespace nttec {
namespace fec {

/** Reed-Solomon (RS) Erasure code over GF(2<sup>n</sup>) using additive FFT. */
template <typename T>
class RsGf2nFftAdd : public FecCode<T> {
  public:
    // NOTE: only NON_SYSTEMATIC is supported now
    RsGf2nFftAdd(unsigned word_size, unsigned n_data, unsigned n_parities)
        : FecCode<T>(FecType::NON_SYSTEMATIC, word_size, n_data, n_parities)
    {
        this->fec_init();
    }

    ~RsGf2nFftAdd()
    {
        if (this->gf)
            delete this->gf;
        if (betas)
            delete betas;
    }

    inline void check_params()
    {
        if (this->word_size > 16)
            assert(false); // not support yet
    }

    inline void init_gf()
    {
        unsigned gf_n = 8 * this->word_size;
        this->gf = new gf::BinExtension<T>(gf_n);
    }

    inline void init_fft()
    {
        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = smallest power of 2 and at least (n_parities + n_data)
        this->n =
            arith::get_smallest_power_of_2<T>(this->n_data + this->n_parities);

        T m = arith::log2<T>(this->n);

        this->fft = std::unique_ptr<fft::Additive<T>>(
            new fft::Additive<T>(this->gf, m));
    }

    inline void init_others()
    {
        // subspace spanned by <beta_i>
        this->betas = new vec::Vector<T>(this->gf, this->n);
        this->fft->compute_B(this->betas);
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

  private:
    vec::Vector<T>* betas = nullptr;

  protected:
    void decode_prepare(
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words,
        vec::Vector<T>* vx,
        int* vx_zero)
    {
        int _vx_zero = -1;
        // vector x=(x_0, x_1, ..., x_k-1)
        for (int i = 0; i < this->n_data; i++) {
            int _vx = this->betas->get(fragments_ids->get(i));
            vx->set(i, _vx);
            if (_vx == 0)
                _vx_zero = i;
        }
        *vx_zero = _vx_zero;
    }

    // Lagrange interpolation
    void decode_lagrange(
        vec::Vector<T>* output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words,
        vec::Vector<T>* vx,
        int vx_zero)
    {
        int k = this->n_data; // number of fragments received
        Polynomial<T> A(this->gf), _A(this->gf);

        // compute A(x) = prod_j(x-x_j)
        A.set(0, 1);
        for (int i = 0; i < k; i++) {
            A.mul_to_x_plus_coef(this->gf->sub(0, vx->get(i)));
        }
        // std::cout << "A(x)="; A.dump();

        // compute A'(x) since A_i(x_i) = A'_i(x_i)
        _A.copy(&A);
        _A.derivative();
        // std::cout << "A'(x)="; _A.dump();

        // evaluate n_i=v_i/A'_i(x_i)
        vec::Vector<T> _n(this->gf, k);
        for (int i = 0; i < k; i++) {
            _n.set(i, this->gf->div(words->get(i), _A.eval(vx->get(i))));
        }

        // We have to find the numerator of the following expression:
        // P(x)/A(x) = S(x) + R(x)
        //  where S(x) = sum_{0 <= i <= k-1, i != vx_zero}(n_i/(x-x_i)) mod x^n
        //        R(x) = _n[vx_zero] / x
        // using Taylor series we rewrite the expression into
        // S(x) = sum_i=0_k-1(sum_j=0_n-1(n_i*x_i^(-j-1)*x^j))
        Polynomial<T> S(this->gf);
        for (int j = 0; j <= k - 1; j++) {
            T val = 0;
            for (int i = 0; i <= k - 1; i++) {
                if (i == vx_zero)
                    continue;
                // perform Taylor series at 0
                T xi_j_1 = this->gf->inv(this->gf->exp(vx->get(i), j + 1));
                val = this->gf->add(val, this->gf->mul(_n.get(i), xi_j_1));
            }
            S.set(j, val);
        }
        // std::cout << "S:"; S.dump();
        S.mul(&A, k - 1);
        // std::cout << "S x A:"; S.dump();
        if (vx_zero > -1) {
            assert(A.get(0) == 0);
            // P(x) = A(x)*S(x) + _n[vx_zero] * A(x) / x
            //  as S(x) does not include the term of vx_zero
            // Note: A(0) = 0 since vx_zero exists
            int deg = A.degree();
            T val = _n.get(vx_zero);
            for (int i = 1; i <= deg; i++)
                S.set(
                    i - 1,
                    this->gf->add(S.get(i - 1), this->gf->mul(val, A.get(i))));
            // std::cout << "S x A + vx_zero:"; S.dump();
        }

        // output is n_data length
        for (unsigned i = 0; i < this->n_data; i++)
            output->set(i, S.get(i));
    }

    // Lagrange interpolation w/ vector
    void decode_vec_lagrange(
        vec::Vector<T>* output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words,
        vec::Vector<T>* vx,
        int vx_zero)
    {
        // std::cout << "WORDS:"; words->dump();
        int k = this->n_data; // number of fragments received

        vec::Poly<T> A(this->gf, this->n);
        vec::Poly<T> _A_fft(this->gf, this->n);
        vec::Poly<T> S(this->gf, k);

        // compute A(x) = prod_j(x-x_j)
        A.zero();
        A.set(0, 1);
        for (int i = 0; i < k; i++) {
            A.mul_to_x_plus_coef(this->gf->sub(0, vx->get(i)));
        }
        // std::cout << "A(x)="; A.dump();

        // compute A'(x) since A_i(x_i) = A'_i(x_i)
        vec::Poly<T> _A(A);
        _A.derivative();
        this->fft->fft(&_A_fft, &_A);
        // std::cout << "A'(x)="; _A.dump();

        // evaluate n_i=v_i/A'_i(x_i)
        vec::Vector<T> _n(this->gf, k);
        for (int i = 0; i < k; i++) {
            _n.set(
                i,
                this->gf->div(
                    words->get(i), _A_fft.get(fragments_ids->get(i))));
        }

        // We have to find the numerator of the following expression:
        // P(x)/A(x) = S(x) + R(x)
        //  where S(x) = sum_{0 <= i <= k-1, i != vx_zero}(n_i/(x-x_i)) mod x^n
        //        R(x) = _n[vx_zero] / x
        // using Taylor series we rewrite the expression into
        // S(x) = sum_i=0_k-1(sum_j=0_n-1(n_i*x_i^(-j-1)*x^j))
        for (int j = 0; j <= k - 1; j++) {
            T val = 0;
            for (int i = 0; i <= k - 1; i++) {
                if (i == vx_zero)
                    continue;
                // perform Taylor series at 0
                T xi_j_1 = this->gf->inv(this->gf->exp(vx->get(i), j + 1));
                val = this->gf->add(val, this->gf->mul(_n.get(i), xi_j_1));
            }
            S.set(j, val);
        }
        // std::cout << "S:"; S.dump();
        S.mul(&A, k - 1);
        // std::cout << "S x A:"; S.dump();
        if (vx_zero > -1) {
            assert(A.get(0) == 0);
            // P(x) = A(x)*S(x) + _n[vx_zero] * A(x) / x
            //  as S(x) does not include the term of vx_zero
            // Note: A(0) = 0 since vx_zero exists
            int deg = A.get_deg();
            T val = _n.get(vx_zero);
            for (int i = 1; i <= deg; i++)
                S.set(
                    i - 1,
                    this->gf->add(S.get(i - 1), this->gf->mul(val, A.get(i))));
            // std::cout << "S x A + vx_zero:"; S.dump();
        }

        // output is n_data length
        for (unsigned i = 0; i < this->n_data; i++)
            output->set(i, S.get(i));
    }

  private:
    // this overwrite multiplicative FFT
    std::unique_ptr<fft::Additive<T>> fft = nullptr;
};

} // namespace fec
} // namespace nttec

#endif
