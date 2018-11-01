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
#ifndef __QUAD_FFT_BASE_H__
#define __QUAD_FFT_BASE_H__

#include "gf_base.h"
#include "vec_buffers.h"
#include "vec_vector.h"

namespace quadiron {

/** Various Fast Fourier Transform (FFT) implementations. */
namespace fft {

/** Base class for Fourier Transform on Galois Fields.
 *
 *  The Fourier Transform is applied on vector of size `n` with ω as
 *  n<sup>th</sup> root of unity.
 */
template <typename T>
class FourierTransform {
  public:
    virtual ~FourierTransform();
    int get_n();
    const gf::Field<T>& get_gf();
    /** Compute the Fourier Transform. */
    virtual void fft(vec::Vector<T>& output, vec::Vector<T>& input) = 0;
    virtual void
    fft(vec::Buffers<T>& /* output */, vec::Buffers<T>& /* input */){};
    /** Compute the Inverse Fourier Transform. */
    virtual void ifft(vec::Vector<T>& output, vec::Vector<T>& input) = 0;
    virtual void
    ifft(vec::Buffers<T>& /* output */, vec::Buffers<T>& /* input */){};
    /** Compute the summation for the inverse FFT formula. */
    virtual void fft_inv(vec::Vector<T>& output, vec::Vector<T>& input) = 0;
    virtual void
    fft_inv(vec::Buffers<T>& /* output */, vec::Buffers<T>& /* input */){};

  protected:
    const gf::Field<T>* gf;
    int n;
    T inv_n_mod_p;
    vec::Vector<T>* vec_inv_n = nullptr;
    FourierTransform(const gf::Field<T>& gf, int n, bool additive = false);
};

template <typename T>
FourierTransform<T>::FourierTransform(
    const gf::Field<T>& gf,
    int n,
    bool additive)
{
    this->gf = &gf;
    this->n = n;

    if (!additive) {
        this->inv_n_mod_p = gf.get_inv_n_mod_p(n);

        this->vec_inv_n = new vec::Vector<T>(gf, n);
        for (int i = 0; i < n; i++) {
            this->vec_inv_n->set(i, this->inv_n_mod_p);
        }
    }
}

template <typename T>
FourierTransform<T>::~FourierTransform()
{
    if (vec_inv_n != nullptr)
        delete vec_inv_n;
}

template <typename T>
int FourierTransform<T>::get_n()
{
    return n;
}

template <typename T>
const gf::Field<T>& FourierTransform<T>::get_gf()
{
    return *gf;
}

} // namespace fft
} // namespace quadiron

#endif
