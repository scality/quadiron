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
#ifndef __QUAD_FEC_RS_GF2N_FFT_ADD_H__
#define __QUAD_FEC_RS_GF2N_FFT_ADD_H__

#include "arith.h"
#include "fec_base.h"
#include "fft_add.h"
#include "gf_bin_ext.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

namespace quadiron {
namespace fec {

/** Reed-Solomon (RS) Erasure code over GF(2<sup>n</sup>) using additive FFT. */
template <typename T>
class RsGf2nFftAdd : public FecCode<T> {
  public:
    using FecCode<T>::decode_apply;
    using FecCode<T>::decode_prepare;
    using FecCode<T>::encode;

    // NOTE: only NON_SYSTEMATIC is supported now
    RsGf2nFftAdd(unsigned word_size, unsigned n_data, unsigned n_parities)
        : FecCode<T>(FecType::NON_SYSTEMATIC, word_size, n_data, n_parities)
    {
        this->fec_init();
    }

    inline void check_params() override
    {
        if (this->word_size > 16)
            assert(false); // not support yet
    }

    inline void init_gf() override
    {
        unsigned gf_n = 8 * this->word_size;
        this->gf = gf::alloc<gf::Field<T>, gf::BinExtension<T>>(gf_n);
    }

    inline void init_fft() override
    {
        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = smallest power of 2 and at least (n_parities + n_data)
        this->n =
            arith::get_smallest_power_of_2<T>(this->n_data + this->n_parities);

        T m = arith::log2<T>(this->n);

        this->fft = std::unique_ptr<fft::Additive<T>>(
            new fft::Additive<T>(*(this->gf), m));
    }

    inline void init_others() override
    {
        // subspace spanned by <beta_i>
        this->betas = std::unique_ptr<vec::Vector<T>>(
            new vec::Vector<T>(*(this->gf), this->n));
        this->fft->compute_B(*betas);
    }

    int get_n_outputs() override
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
        vec::Vector<T>& output,
        std::vector<Properties>&,
        off_t,
        vec::Vector<T>& words) override
    {
        vec::ZeroExtended<T> vwords(words, this->n);
        this->fft->fft(output, vwords);
    }

    void decode_add_data(int, int) override
    {
        // not applicable
        assert(false);
    }

  private:
    std::unique_ptr<vec::Vector<T>> betas = nullptr;

  protected:
    std::unique_ptr<DecodeContext<T>> init_context_dec(
        vec::Vector<T>& fragments_ids,
        size_t,
        vec::Buffers<T>*) override
    {
        if (this->betas == nullptr) {
            throw LogicError("FEC FFT ADD: vector 'betas' must be initialized");
        }
        if (this->fft == nullptr) {
            throw LogicError("FEC FFT ADD: FFT must be initialized");
        }

        int k = this->n_data; // number of fragments received
        // vector x=(x_0, x_1, ..., x_k-1)
        vec::Vector<T> vx(*(this->gf), k);

        int vx_zero = -1;
        for (unsigned i = 0; i < this->n_data; ++i) {
            T val = betas->get(fragments_ids.get(i));
            vx.set(i, val);
            if (val == 0) {
                vx_zero = i;
            }
        }

        std::unique_ptr<DecodeContext<T>> context =
            std::unique_ptr<DecodeContext<T>>(new DecodeContext<T>(
                *(this->gf),
                *fft,
                *(this->fft_2k),
                fragments_ids,
                vx,
                this->n_data,
                this->n,
                vx_zero));

        return context;
    }

    void decode_prepare(
        const DecodeContext<T>&,
        const std::vector<Properties>&,
        off_t,
        vec::Vector<T>&) override
    {
        // nothing to do
    }

    void decode_apply(
        const DecodeContext<T>& context,
        vec::Vector<T>& output,
        vec::Vector<T>& words) override
    {
        const vec::Vector<T>& fragments_ids = context.get_fragments_id();
        vec::Poly<T>& A = context.get_poly(CtxPoly::A);
        vec::Vector<T>& inv_A_i = context.get_vector(CtxVec::INV_A_I);
        vec::Vector<T>& vec1_n = context.get_vector(CtxVec::N1);
        vec::Poly<T>& S = context.get_poly(CtxPoly::S);

        int k = this->n_data; // number of fragments received
        int vx_zero = context.vx_zero;

        // FIXME: split this step in decode_init as multiplicative FFT
        vec::Vector<T> vx(*(this->gf), k);
        for (unsigned i = 0; i < this->n_data; ++i) {
            vx.set(i, this->betas->get(fragments_ids.get(i)));
        }

        // compute N'(x) = sum_i{n_i * x^z_i}
        // where n_i=v_i/A'_i(x_i)
        vec1_n.zero_fill();
        for (int i = 0; i <= k - 1; ++i) {
            vec1_n.set(i, this->gf->mul(words.get(i), inv_A_i.get(i)));
        }

        // We have to find the numerator of the following expression:
        // P(x)/A(x) = S(x) + R(x)
        //  where S(x) = sum_{0 <= i <= k-1, i != vx_zero}(n_i/(x-x_i)) mod x^n
        //        R(x) = _n[vx_zero] / x
        // using Taylor series we rewrite the expression into
        // S(x) = sum_i=0_k-1(sum_j=0_n-1(n_i*x_i^(-j-1)*x^j))
        for (int j = 0; j <= k - 1; j++) {
            T val = 0;
            for (int i = 0; i <= k - 1; ++i) {
                if (i == vx_zero)
                    continue;
                // perform Taylor series at 0
                T xi_j_1 = this->gf->inv(this->gf->exp(vx.get(i), j));
                val = this->gf->add(val, this->gf->mul(vec1_n.get(i), xi_j_1));
            }
            S.set(j, val);
        }
        S.mul(&A, k - 1);
        if (vx_zero > -1) {
            assert(A.get(0) == 0);
            // P(x) = A(x)*S(x) + _n[vx_zero] * A(x) / x
            //  as S(x) does not include the term of vx_zero
            // Note: A(0) = 0 since vx_zero exists
            int deg = A.get_deg();
            T val = vec1_n.get(vx_zero);
            for (int i = 1; i <= deg; ++i)
                S.set(
                    i - 1,
                    this->gf->add(S.get(i - 1), this->gf->mul(val, A.get(i))));
        }

        // No need to mod x^n since only last n_data coefs are obtained
        // output is n_data length
        for (unsigned i = 0; i < this->n_data; ++i)
            output.set(i, S.get(i));
    }

  private:
    // this overwrite multiplicative FFT
    std::unique_ptr<fft::Additive<T>> fft = nullptr;
};

} // namespace fec
} // namespace quadiron

#endif
