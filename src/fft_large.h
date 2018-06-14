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
#ifndef __NTTEC_FFT_LARGE_H__
#define __NTTEC_FFT_LARGE_H__

#include "arith.h"
#include "fft_base.h"
#include "gf_base.h"
#include "vec_matrix.h"
#include "vec_vector.h"

namespace nttec {
namespace fft {

/** FFT implementation specialized for large FFT size.
 *
 * The implementation comes from @cite meunier.
 */
template <typename T>
class Large : public FourierTransform<T> {
  public:
    Large(const gf::Field<T>& gf, int l, T w);
    ~Large();
    void fft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void ifft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void fft_inv(vec::Vector<T>& output, vec::Vector<T>& input) override;

  private:
    void
    _fft(vec::Vector<T>& output, vec::Vector<T>& input, vec::Vector<T>& _W);

    int l;
    T w;
    T inv_w;
    vec::Vector<T>* W;
    vec::Vector<T>* inv_W;
    vec::Matrix<T>* tmp2;
    vec::Matrix<T>* tmp3;
    vec::Matrix<T>* tmp4;
    vec::Vector<T>& tmp5;
    int _get_p(int i, int j);
    int _get_p0(int i, int j, int s);
    int _get_p1(int i, int j, int s);
    void _pre_compute_consts();
};

template <typename T>
Large<T>::Large(const gf::Field<T>& gf, int l, T w)
    : FourierTransform<T>(gf, arith::exp2<T>(l))
{
    this->l = l;
    this->w = w;
    this->inv_w = gf->inv(w);

    this->W = new vec::Vector<T>(gf, this->n + 1);
    this->inv_W = new vec::Vector<T>(gf, this->n + 1);
    gf->compute_omegas_cached(W, this->n, w);
    gf->compute_omegas_cached(inv_W, this->n, inv_w);

    _pre_compute_consts();
}

template <typename T>
Large<T>::~Large()
{
    delete this->inv_W;
    delete this->W;
    delete this->tmp2;
    delete this->tmp3;
    delete this->tmp4;
    delete this->tmp5;
}

/**
 * Representing powers of root as a vector is more practical than a
 * matrix (e.g. (n+k=2^15) => would be 1 billion entries !)
 *
 * @param _W vector of powers of roots
 * @param _w Nth root of unity
 */
/**
 * Simulate the bit matrix
 *
 * [0][0]=undef   [0][1]=0   [0][2]=0   ... [0][n]=0    <= 0
 * [1][0]=undef   [1][1]=1   [1][2]=0   ... [1][n]=0    <= 1
 * [2][0]=undef   [2][1]=0   [2][2]=1   ... [2][n]=0    <= 2
 * ...
 * [N-1][0]=undef [N-1][1]=1 [N-1][2]=1 ... [N-1][n]=1  <= N-1
 *
 * where the numbers from 0 .. N-1 are encoded in reverse order
 *
 * @param i 0 <= i <= N-1
 * @param j 1 <= j <= n
 *
 * @return
 */
template <typename T>
int Large<T>::_get_p(int i, int j)
{
    assert(i >= 0 && i <= this->n - 1);
    assert(j >= 1 && j <= this->l);

    int x = i;
    int y = 1;

    do {
        if (y == j) {
            if (x & 1)
                return 1;
            else
                return 0;
        }
        y++;
    } while (x >>= 1);

    return 0;
}

/**
 * @return _get_p(i, j) except when j==s it returns 0
 */
template <typename T>
int Large<T>::_get_p0(int i, int j, int s)
{
    assert(i >= 0 && i <= this->n - 1);
    assert(j >= 1 && j <= this->l);

    return (j == s) ? 0 : _get_p(i, j);
}

/**
 * @return _get_p(i, j) except when j==s it returns 1
 */
template <typename T>
int Large<T>::_get_p1(int i, int j, int s)
{
    assert(i >= 0 && i <= this->n - 1);
    assert(j >= 1 && j <= this->l);

    return (j == s) ? 1 : _get_p(i, j);
}

template <typename T>
void Large<T>::_pre_compute_consts()
{
    tmp2 = new vec::Matrix<T>(this->gf, this->l + 1, this->n);
    tmp3 = new vec::Matrix<T>(this->gf, this->l + 1, this->n);
    tmp4 = new vec::Matrix<T>(this->gf, this->l + 1, this->n);
    tmp5 = new vec::Vector<T>(this->gf, this->n);

    tmp2->zero_fill();
    tmp3->zero_fill();
    tmp4->zero_fill();
    tmp5->zero_fill();

    for (int i = 1; i <= this->l; i++) {
        for (int j = 0; j <= this->n - 1; j++) {
            T _tmp1 = 0;
            for (int k = 1; k <= i; k++)
                _tmp1 += _get_p(j, this->l - i + k) * arith::exp2<T>(i - k);

            T _tmp2 = 0;
            for (int k = 1; k <= this->l; k++)
                _tmp2 += _get_p0(j, k, this->l - i + 1) * arith::exp2<T>(k - 1);
            tmp2->set(i, j, _tmp2);

            T _tmp3 = 0;
            for (int k = 1; k <= this->l; k++)
                _tmp3 += _get_p1(j, k, this->l - i + 1) * arith::exp2<T>(k - 1);
            tmp3->set(i, j, _tmp3);

            T _tmp4 = _tmp1 * arith::exp2<T>(this->l - i);
            tmp4->set(i, j, _tmp4);
        }
    }

    for (int i = 0; i <= this->n - 1; i++) {
        T _tmp5 = 0;
        for (int k = 1; k <= this->l; k++)
            _tmp5 += _get_p(i, this->l - k + 1) * arith::exp2<T>(k - 1);
        tmp5->set(i, _tmp5);
    }
}

template <typename T>
void Large<T>::_fft(
    vec::Vector<T>& output,
    vec::Vector<T>& input,
    vec::Vector<T>& _W)
{
    vec::Matrix<T> phi(this->gf, this->l + 1, this->n);

    // compute phi[0][i]
    for (int i = 0; i <= this->n - 1; i++)
        phi.set(0, i, input->get(i));

    for (int i = 1; i <= this->l; i++) {
        for (int j = 0; j <= this->n - 1; j++) {
            DoubleSizeVal<T> val = DoubleSizeVal<T>(_W->get(tmp4->get(i, j)))
                                       * phi.get(i - 1, tmp3->get(i, j))
                                   + phi.get(i - 1, tmp2->get(i, j));

            phi.set(i, j, val % this->gf->card());
        }
    }

    // compute FFT
    for (int i = 0; i <= this->n - 1; i++)
        output.set(i, phi.get(this->l, tmp5->get(i)));
}

template <typename T>
void Large<T>::fft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    _fft(output, input, W);
}

template <typename T>
void Large<T>::fft_inv(vec::Vector<T>& output, vec::Vector<T>& input)
{
    _fft(output, input, inv_W);
}

template <typename T>
void Large<T>::ifft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    fft_inv(output, input);
    if (this->inv_n_mod_p > 1)
        output.mul_scalar(this->inv_n_mod_p);
}

} // namespace fft
} // namespace nttec

#endif
