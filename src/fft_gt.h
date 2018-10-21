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
#ifndef __QUAD_FFT_GT_H__
#define __QUAD_FFT_GT_H__

#include "arith.h"
#include "fft_2.h"
#include "fft_2n.h"
#include "fft_base.h"
#include "fft_ct.h"
#include "fft_naive.h"
#include "gf_base.h"
#include "vec_vector.h"
#include "vec_view.h"

namespace quadiron {
namespace fft {

/** FFT implementation using the Good-Thomas algorithm.
 *
 * Suppose \f$n = n_1 \times n_2\f$ where \f$n_1\f$ and \f$n_2\f$ are relatively
 * prime, i.e. \f$\gcd(n_1, n_2) = 1\f$
 *
 * \f[
 *   X[k] = \sum_{i=0}^{n-1} x_i \times {w}^{ik}
 * \f]
 *
 * is plit into two smaller DFTs by using the index mapping
 *
 * \f[
 *   i = a \times i_1 + b \times i_2
 * \f]
 * \f[
 *   k = c \times k_1 + d \times k_2
 * \f]
 *
 * with
 *   \f$a = n_2\f$,
 *   \f$b = n_1\f$,
 *   \f$c = n_2*inv\_mod(n_2, n_1)\f$,
 *   \f$d = n_1*inv\_mod(n_1, n_2)\f$
 * where
 *   \f$1 \leq i_1, k_1 \leq n_1-1\f$
 * and
 *   \f$1 \leq i_2, k_2 \leq n_2-1\f$.
 * and
 *   inv_mod(u, v) = multiplicative inverse of u reduced modulo v,
 *   i.e. \f$(u \times inv\_mod(u, v)) \mod v = 1\f$
 *
 * \f[
 *   X[c \times k1 + d \times k2] = \sum_{i_1}
 *     \left\{
 *       \sum_{i_2}{x[a \times i_1 + b \times i_2] \times {w_2}^{i_2  k_2}}
 *     \right\}
 *     \times {w_1}^{i_1 k_1}
 * \f]
 *
 * where
 * \f$w\f$ is \f$n^\text{th}\f$ root,
 * \f$w_1 = w^{n_2}, \text{ and } w_2 = w^{n_1}\f$
 *
 * Note, the index mapping leads to
 *  \f$w^{i k} = w_{1}^{i_1 k_1} \times w_{2}^{i_2 k_2}\f$
 *
 * - Step1: calculate DFT of the inner parenthese, i.e. \f$\sum_{i_2}\f$
 * - Step2: calculate DFT of the outer parenthese, i.e. \f$\sum_{i_1}\f$
 *
 * @see <a
 * href="https://en.wikipedia.org/wiki/Prime-factor_FFT_algorithm">Prime-factor
 * FFT algorithm</a>
 */
template <typename T>
class GoodThomas : public FourierTransform<T> {
  public:
    using FourierTransform<T>::fft;
    using FourierTransform<T>::ifft;
    using FourierTransform<T>::fft_inv;

    GoodThomas(
        const gf::Field<T>& gf,
        T n,
        int id = 0,
        std::vector<T>* factors = nullptr,
        T _w = 0);
    ~GoodThomas();
    void fft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void ifft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void fft_inv(vec::Vector<T>& output, vec::Vector<T>& input) override;

  private:
    void _fft(vec::Vector<T>& output, vec::Vector<T>& input, bool inv);
    T _inverse_mod(T nb, T mod);

    bool loop;
    bool first_layer_fft;
    T n1;
    T n2;
    T w, w1, w2;
    T a, b, c, d;
    vec::Vector<T>* G = nullptr;
    vec::View<T>* Y = nullptr;
    vec::View<T>* X = nullptr;
    FourierTransform<T>* dft_outer = nullptr;
    FourierTransform<T>* dft_inner = nullptr;
    std::vector<T> prime_factors;
};

/**
 * n-th root will be constructed with primitive root
 *
 * @param gf
 * @param n
 * @param id: index in the list of factors of n
 *
 * @return
 */
template <typename T>
GoodThomas<T>::GoodThomas(
    const gf::Field<T>& gf,
    T n,
    int id,
    std::vector<T>* factors,
    T _w)
    : FourierTransform<T>(gf, n)
{
    if (factors == nullptr) {
        first_layer_fft = true;
        this->prime_factors = arith::get_coprime_factors<T>(n);
        // w is of order-n
        w = gf.get_nth_root(n);
    } else {
        this->prime_factors = *factors;
        first_layer_fft = false;
        w = _w;
    }
    n1 = prime_factors[id];
    n2 = n / n1;

    a = n2;
    b = n1;
    c = n2 * (this->_inverse_mod(n2, n1));
    d = n1 * (this->_inverse_mod(n1, n2));

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
        if (arith::is_power_of_2<T>(_n2)) {
            this->dft_inner = new fft::Radix2<T>(gf, _n2);
        } else {
            this->dft_inner = new fft::CooleyTukey<T>(
                gf, _n2, id + 1, &this->prime_factors, w2);
        }
        this->G = new vec::Vector<T>(gf, this->n);
        this->Y = new vec::View<T>(this->G);
        this->X = new vec::View<T>(this->G);
    } else
        loop = false;
}

template <typename T>
GoodThomas<T>::~GoodThomas()
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

/*
 * It calculates multiplicative inverse of a number modulo mod
 *  NOTE: it should be optimized
 */
template <typename T>
T GoodThomas<T>::_inverse_mod(T nb, T mod)
{
    for (T n = 1; n < mod; n++) {
        if ((n * nb) % mod == 1)
            return n;
    }
    return 0;
}

template <typename T>
void GoodThomas<T>::_fft(
    vec::Vector<T>& output,
    vec::Vector<T>& input,
    bool inv)
{
    X->set_vec(&input);
    X->set_len(n2);
    Y->set_len(n2);
    for (T i1 = 0; i1 < n1; i1++) {
        Y->set_map(i1, n1);
        X->set_map((a * i1) % this->n, b);
        if (inv)
            this->dft_inner->fft_inv(*Y, *X);
        else
            this->dft_inner->fft(*Y, *X);
        // std::cout << "X:"; X->dump();
        // std::cout << "Y:"; Y->dump();
    }

    X->set_vec(&output);
    X->set_len(n1);
    Y->set_len(n1);
    for (T k2 = 0; k2 < n2; k2++) {
        Y->set_map(k2 * n1, 1);
        X->set_map((d * k2) % this->n, c);
        if (inv)
            this->dft_outer->fft_inv(*X, *Y);
        else
            this->dft_outer->fft(*X, *Y);
        // std::cout << "Y:"; Y->dump();
        // std::cout << "X:"; X->dump();
    }
}

template <typename T>
void GoodThomas<T>::fft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    if (!loop)
        return dft_outer->fft(output, input);
    else
        return _fft(output, input, false);
}

template <typename T>
void GoodThomas<T>::fft_inv(vec::Vector<T>& output, vec::Vector<T>& input)
{
    if (!loop)
        dft_outer->fft_inv(output, input);
    else
        _fft(output, input, true);
}

template <typename T>
void GoodThomas<T>::ifft(vec::Vector<T>& output, vec::Vector<T>& input)
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
