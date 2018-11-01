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
#ifndef __QUAD_FFT_CT_H__
#define __QUAD_FFT_CT_H__

#include "arith.h"
#include "fft_2.h"
#include "fft_base.h"
#include "fft_naive.h"
#include "gf_base.h"
#include "vec_vector.h"
#include "vec_view.h"

namespace quadiron {
namespace fft {

/** FFT implementation using the Cooley–Tukey algorithm
 *
 * Implementation based on @cite fft-ct-mitra
 *
 * Suppose \f$n = n_1 \times n_2\f$, the DFT
 *
 * \f[
 *   X[k] = \sum_{i=0}^{n-1} x_i \times {w}^{ik}
 * \f]
 *
 * is split into two smaller DFTs by using the index mapping
 *
 * \f[
 *   i = i_1 + n_1 \times i_2
 * \f]
 * \f[
 *   k = n_2 \times k_1 + k_2
 * \f]
 * where
 *   \f$1 \leq i_1, k_1 \leq n_1-1\f$
 * and
 *   \f$1 \leq i_2, k_2 \leq n_2-1\f$.
 *
 * Then we have
 * \f[
 *   X[n_2 k_1 + k_2] = \sum_{i_1}
 *     \left\{
 *       \sum_{i_2}{x[i_1 + n_1 i_2] \times {w_2}^{i_2  k_2}}
 *     \right\}
 *     \times w^{i_1 k_2} \times {w_1}^{i_1 k_1}
 * \f]
 *
 * where
 * \f$w\f$ is \f$n^\text{th}\f$ root,
 * \f$w_1 = w^{n_2}, \text{ and } w_2 = w^{n_1}\f$
 *
 * - Step1: calculate the inner DFT, i.e. \f$\sum_{i_2}\f$
 * - Step2: multiply to twiddle factors \f$w^{i_1 k_2}\f$
 * - Step3: calculate outer DFT, i.e. \f$\sum_{i_1}\f$
 */
template <typename T>
class CooleyTukey : public FourierTransform<T> {
  public:
    using FourierTransform<T>::fft;
    using FourierTransform<T>::ifft;
    using FourierTransform<T>::fft_inv;

    CooleyTukey(
        const gf::Field<T>& gf,
        T n,
        int id = 0,
        std::vector<T>* factors = nullptr,
        T _w = 0);
    ~CooleyTukey();
    void fft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void ifft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void fft_inv(vec::Vector<T>& output, vec::Vector<T>& input) override;

  private:
    void _fft(vec::Vector<T>& output, vec::Vector<T>& input, bool inv);

    bool loop;
    bool first_layer_fft;
    T n1;
    T n2;
    T w, w1, w2, inv_w;
    vec::Vector<T>* G = nullptr;
    vec::View<T>* Y = nullptr;
    vec::View<T>* X = nullptr;
    FourierTransform<T>* dft_outer = nullptr;
    FourierTransform<T>* dft_inner = nullptr;
    std::vector<T> prime_factors;
    void mul_twiddle_factors(bool inv);
};

/** Initialize the FFT.
 *
 * n-th root will be constructed with primitive root
 *
 * @param id index in the list of factors of n
 */
template <typename T>
CooleyTukey<T>::CooleyTukey(
    const gf::Field<T>& gf,
    T n,
    int id,
    std::vector<T>* factors,
    T _w)
    : FourierTransform<T>(gf, n)
{
    if (factors == nullptr) {
        first_layer_fft = true;
        this->prime_factors = arith::get_prime_factors<T>(n);
        // w is of order-n
        w = gf.get_nth_root(n);
    } else {
        this->prime_factors = *factors;
        first_layer_fft = false;
        w = _w;
    }
    inv_w = gf.inv(w);
    n1 = prime_factors[id];
    n2 = n / n1;

    w1 = gf.exp(w, n2); // order of w1 = n1
    if (n1 == 2) {
        this->dft_outer = new fft::Size2<T>(gf);
    } else {
        this->dft_outer = new fft::Naive<T>(gf, n1, w1);
    }

    if (n2 > 1) {
        loop = true;
        w2 = gf.exp(w, n1); // order of w2 = n2
        T _n2 = n / n1;
        // TODO: fix that fft::Radix2 does not work here
        // if (_is_power_of_2<T>(_n2))
        //   this->dft_inner = new fft::Radix2<T>(gf, _n2);
        // else
        this->dft_inner =
            new CooleyTukey<T>(gf, _n2, id + 1, &this->prime_factors, w2);
        this->G = new vec::Vector<T>(gf, this->n);
        this->Y = new vec::View<T>(this->G);
        this->X = new vec::View<T>(this->G);
    } else
        loop = false;
}

template <typename T>
CooleyTukey<T>::~CooleyTukey()
{
    if (dft_outer)
        delete dft_outer;
    if (dft_inner)
        delete dft_inner;
    if (X)
        delete X;
    if (Y)
        delete Y;
    if (G)
        delete G;
}

template <typename T>
void CooleyTukey<T>::mul_twiddle_factors(bool inv)
{
    T factor;
    T _w;
    if (inv)
        _w = inv_w;
    else
        _w = w;
    T base = 1;
    for (T i1 = 1; i1 < n1; i1++) {
        base = this->gf->mul(base, _w); // base = _w^i1
        factor = base;                  // init factor = base^1
        for (T k2 = 1; k2 < n2; k2++) {
            T loc = i1 + n1 * k2;
            ;
            G->set(loc, this->gf->mul(G->get(loc), factor));
            // next factor = base^(k2+1)
            factor = this->gf->mul(factor, base);
        }
    }
}

template <typename T>
void CooleyTukey<T>::_fft(
    vec::Vector<T>& output,
    vec::Vector<T>& input,
    bool inv)
{
    X->set_vec(&input);
    X->set_len(n2);
    Y->set_len(n2);
    for (T i1 = 0; i1 < n1; i1++) {
        Y->set_map(i1, n1);
        X->set_map(i1, n1);
        if (inv)
            this->dft_inner->fft_inv(*Y, *X);
        else
            this->dft_inner->fft(*Y, *X);
        // std::cout << "X:"; X->dump();
        // std::cout << "Y:"; Y->dump();
    }

    // multiply to twiddle factors
    mul_twiddle_factors(inv);

    X->set_vec(&output);
    X->set_len(n1);
    Y->set_len(n1);
    for (T k2 = 0; k2 < n2; k2++) {
        Y->set_map(k2 * n1, 1);
        X->set_map(k2, n2);
        if (inv)
            this->dft_outer->fft_inv(*X, *Y);
        else
            this->dft_outer->fft(*X, *Y);
        // std::cout << "Y:"; Y->dump();
        // std::cout << "X:"; X->dump();
    }
}

template <typename T>
void CooleyTukey<T>::fft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    if (!loop)
        return dft_outer->fft(output, input);
    else
        return _fft(output, input, false);
}

template <typename T>
void CooleyTukey<T>::fft_inv(vec::Vector<T>& output, vec::Vector<T>& input)
{
    if (!loop)
        dft_outer->fft_inv(output, input);
    else
        _fft(output, input, true);
}

template <typename T>
void CooleyTukey<T>::ifft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    fft_inv(output, input);

    /*
     * We need to divide output to `N` for the inverse formular
     */

    if (this->first_layer_fft && (this->inv_n_mod_p > 1)) {
        output.mul_scalar(this->inv_n_mod_p);
    }
}

} // namespace fft
} // namespace quadiron

#endif
