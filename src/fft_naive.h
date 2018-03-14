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
#ifndef __NTTEC_FFT_NAIVE_H__
#define __NTTEC_FFT_NAIVE_H__

#include "fft_base.h"
#include "gf_base.h"
#include "vec_matrix.h"
#include "vec_vector.h"

namespace nttec {
namespace fft {

/** Naïve implementation of the Discrete Fourier Transform (DFT).
 *
 * This algorithm can handle any kind of size (even if it's neither a power of
 * two nor a prime number), but is limited to very small inputs because of its
 * algorithmic complexity: O(n²).
 */
template <typename T>
class Naive : public FourierTransform<T> {
  public:
    Naive(gf::Field<T>* gf, int n, T w);
    ~Naive();
    void fft(vec::Vector<T>* output, vec::Vector<T>* input);
    void ifft(vec::Vector<T>* output, vec::Vector<T>* input);
    void fft_inv(vec::Vector<T>* output, vec::Vector<T>* input);

  private:
    T w;
    T inv_w;
    vec::Matrix<T>* W;
    vec::Matrix<T>* inv_W;
    void compute_W(vec::Matrix<T>* _W, T _w);
    void
    _fft(vec::Vector<T>* output, vec::Vector<T>* input, vec::Matrix<T>* _W);
};

template <typename T>
Naive<T>::Naive(gf::Field<T>* gf, int n, T w) : FourierTransform<T>(gf, n)
{
    this->w = w;
    this->inv_w = gf->inv(w);
    this->W = new vec::Matrix<T>(gf, this->n, this->n);
    this->inv_W = new vec::Matrix<T>(gf, this->n, this->n);

    compute_W(W, w);
    compute_W(inv_W, this->inv_w);
}

template <typename T>
Naive<T>::~Naive()
{
    delete this->inv_W;
    delete this->W;
}

/**
 * Compute matrix
 *
 * @param _W matrix of powers of roots
 * @param _w nth root of unity
 */
template <typename T>
void Naive<T>::compute_W(vec::Matrix<T>* _W, T _w)
{
    for (int i = 0; i <= this->n - 1; i++) {
        for (int j = 0; j <= this->n - 1; j++) {
            _W->set(i, j, this->gf->exp(_w, (i * j) % this->n));
        }
    }
    // _W->dump();
}

template <typename T>
void Naive<T>::_fft(
    vec::Vector<T>* output,
    vec::Vector<T>* input,
    vec::Matrix<T>* _W)
{
    _W->mul(output, input);
}

template <typename T>
void Naive<T>::fft(vec::Vector<T>* output, vec::Vector<T>* input)
{
    _fft(output, input, W);
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void Naive<T>::fft_inv(vec::Vector<T>* output, vec::Vector<T>* input)
{
    _fft(output, input, inv_W);
}

template <typename T>
void Naive<T>::ifft(vec::Vector<T>* output, vec::Vector<T>* input)
{
    _fft(output, input, inv_W);
    if (this->inv_n_mod_p > 1)
        output->mul_scalar(this->inv_n_mod_p);
}

} // namespace fft
} // namespace nttec

#endif
