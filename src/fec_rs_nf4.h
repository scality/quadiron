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
#ifndef __NTTEC_FEC_RS_NF4_H__
#define __NTTEC_FEC_RS_NF4_H__

#include <string>

#include "fec_base.h"
#include "fft_2n.h"
#include "gf_base.h"
#include "gf_nf4.h"
#include "polynomial.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

namespace nttec {
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

    ~RsNf4()
    {
        if (ngff4)
            delete ngff4;
    }

    inline void check_params()
    {
        assert(this->word_size >= 2);
        assert(this->word_size <= 8);
    }

    inline void init_gf()
    {
        gf_n = this->word_size / 2;

        ngff4 = new gf::NF4<T>(gf_n);
        this->gf = ngff4;

        sub_field = ngff4->get_sub_field();
    }

    inline void init_fft()
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
            new fft::Radix2<T>(ngff4, this->n, m));

        this->fft_full =
            std::unique_ptr<fft::Radix2<T>>(new fft::Radix2<T>(ngff4, this->n));
    }

    inline void init_others()
    {
        // vector stores r^{-i} for i = 0, ... , k
        T inv_r = ngff4->inv(this->r);
        this->inv_r_powers = std::unique_ptr<vec::Vector<T>>(
            new vec::Vector<T>(ngff4, this->n_data + 1));
        for (int i = 0; i <= this->n_data; i++)
            this->inv_r_powers->set(i, ngff4->exp(inv_r, i));
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
    gf::Field<uint32_t>* sub_field;
    gf::NF4<T>* ngff4;
    int gf_n;

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
            vx->set(
                i,
                this->gf->exp(
                    this->r, ngff4->replicate(fragments_ids->get(i))));
        }
        // std::cout << "vx"; vx->dump();

        T true_val;
        for (int i = 0; i < k; i++) {
            const int j = fragments_ids->get(i);
            auto data = props[j].get(ValueLocation(offset, j));

            if (data) {
                uint32_t flag = std::stoul(*data);
                // std::cout << "\nflag at " << buf << ":" << flag << std::endl;
                true_val = ngff4->pack(words->get(i), flag);
                // std::cout << "word:" << words->get(i) << " -> " << true_val
                // << std::endl;
            } else {
                true_val = ngff4->pack(words->get(i));
            }
            words->set(i, true_val);
        }
        // std::cout << "words packed"; words->dump();
    }

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
        Polynomial<T> A(ngff4), _A(ngff4);

        // compute A(x) = prod_j(x-x_j)
        T one = ngff4->get_unit();
        A.set(0, one);
        for (int i = 0; i < k; i++) {
            A.mul_to_x_plus_coef(this->gf->sub(0, vx->get(i)));
        }
        // std::cout << "A(x)="; A.dump();

        // compute A'(x) since A_i(x_i) = A'_i(x_i)
        // _A.derivative();
        for (int i = 1; i <= A.degree(); i++)
            _A.set(i - 1, ngff4->mul(A.get(i), ngff4->replicate(i)));

        // std::cout << "A'(x)="; _A.dump();

        // evaluate n_i=v_i/A'_i(x_i)
        vec::Vector<T> _n(ngff4, k);
        for (int i = 0; i < k; i++) {
            _n.set(i, ngff4->div(words->get(i), _A.eval(vx->get(i))));
        }
        // std::cout << "_n="; _n.dump();

        // compute N'(x) = sum_i{n_i * x^z_i}
        Polynomial<T> N_p(ngff4);
        for (int i = 0; i <= k - 1; i++) {
            N_p.set(fragments_ids->get(i), _n.get(i));
        }

        // //std::cout << "N_p="; N_p.dump();
        // We have to find the numerator of the following expression:
        // P(x)/A(x) = sum_i=0_k-1(n_i/(x-x_i)) mod x^n
        // using Taylor series we rewrite the expression into
        // P(x)/A(x) = -sum_i=0_k-1(sum_j=0_n-1(n_i*x_i^(-j-1)*x^j))
        Polynomial<T> S(ngff4);
        for (T i = 0; i <= k - 1; i++) {
            S.set(i, N_p.eval(this->inv_r_powers->get(i + 1)));
        }
        // std::cout << "S="; S.dump();
        S.neg();
        S.mul(&A, k - 1);
        // std::cout << "S="; S.dump();
        // No need to mod x^n since only last n_data coefs are obtained
        // output is n_data length
        for (unsigned i = 0; i < this->n_data; i++)
            output->set(i, ngff4->unpack(S.get(i)).values);
        // std::cout << "decoded"; output->dump();
    }

    void decode_vec_lagrange(
        vec::Vector<T>* output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words,
        vec::Vector<T>* vx,
        int vx_zero)
    {
        int k = this->n_data; // number of fragments received
        vec::Poly<T> A(this->gf, k + 1);
        vec::Poly<T> _A(this->gf, this->n);
        vec::Poly<T> _A_fft(this->gf, this->n);
        vec::Poly<T> N_p(this->gf, this->n);
        vec::Poly<T> N_p_ifft(this->gf, this->n);
        vec::Poly<T> S(this->gf, k);

        // compute A(x) = prod_j(x-x_j)
        T one = ngff4->get_unit();
        A.set(0, one);
        for (int i = 0; i < k; i++) {
            A.mul_to_x_plus_coef(this->gf->sub(0, vx->get(i)));
        }
        // std::cout << "A(x)="; A.dump();

        // compute A'(x) since A_i(x_i) = A'_i(x_i)
        // _A.derivative();
        _A.zero();
        for (int i = 1; i <= A.get_deg(); i++)
            _A.set(i - 1, ngff4->mul(A.get(i), ngff4->replicate(i)));

        this->fft->fft(&_A_fft, &_A);
        // std::cout << "A'(x)="; _A.dump();

        // evaluate n_i=v_i/A'_i(x_i)
        vec::Vector<T> _n(ngff4, k);
        for (int i = 0; i < k; i++) {
            // _n.set(i, ngff4->div(words->get(i), _A.eval(vx->get(i))));
            _n.set(
                i,
                ngff4->div(words->get(i), _A_fft.get(fragments_ids->get(i))));
        }
        // std::cout << "_n="; _n.dump();

        // compute N'(x) = sum_i{n_i * x^z_i}
        N_p.zero();
        for (int i = 0; i <= k - 1; i++) {
            N_p.set(fragments_ids->get(i), _n.get(i));
        }
        // std::cout << "N_p="; N_p.dump();
        this->fft_full->fft_inv(&N_p_ifft, &N_p);

        // We have to find the numerator of the following expression:
        // P(x)/A(x) = sum_i=0_k-1(n_i/(x-x_i)) mod x^n
        // using Taylor series we rewrite the expression into
        // P(x)/A(x) = -sum_i=0_k-1(sum_j=0_n-1(n_i*x_i^(-j-1)*x^j))
        for (int i = 0; i <= k - 1; i++) {
            S.set(i, N_p_ifft.get(i + 1));
        }
        // std::cout << "S="; S.dump();
        S.neg();
        S.mul(&A, k - 1);
        // std::cout << "S="; S.dump();
        // No need to mod x^n since only last n_data coefs are obtained
        // output is n_data length
        for (unsigned i = 0; i < this->n_data; i++)
            output->set(i, ngff4->unpack(S.get(i)).values);
        // std::cout << "decoded"; output->dump();
    }
};

} // namespace fec
} // namespace nttec

#endif
