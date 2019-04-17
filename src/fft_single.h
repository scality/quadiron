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
#ifndef __QUAD_FFT_SINGLE_H__
#define __QUAD_FFT_SINGLE_H__

#include "fft_base.h"
#include "gf_base.h"
#include "vec_buffers.h"
#include "vec_vector.h"

namespace quadiron {
namespace fft {

/** FFT specialization for when all but the first element are zeros.
 *
 * In such case, all elements of output are equals the first element of input.
 */
template <typename T>
class Single : public FourierTransform<T> {
  public:
    explicit Single(const gf::Field<T>& gf, int n);
    ~Single() = default;
    void fft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void ifft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void fft_inv(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void fft(vec::Buffers<T>& output, vec::Buffers<T>& input) override;
    void ifft(vec::Buffers<T>& output, vec::Buffers<T>& input) override;
    void fft_inv(vec::Buffers<T>& output, vec::Buffers<T>& input) override;
};

template <typename T>
Single<T>::Single(const gf::Field<T>& gf, int n) : FourierTransform<T>(gf, n)
{
}

template <typename T>
void Single<T>::fft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    T val = input.get(0);
    output.fill(val);
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void Single<T>::fft_inv(vec::Vector<T>& output, vec::Vector<T>& input)
{
    fft(output, input);
}

template <typename T>
void Single<T>::ifft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    output.zero_fill();
    output.set(0, input.get(0));
}

template <typename T>
void Single<T>::fft(vec::Buffers<T>& output, vec::Buffers<T>& input)
{
    for (int i = 0; i < this->n; i++) {
        output.copy(input, 0, i);
    }
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void Single<T>::fft_inv(vec::Buffers<T>& output, vec::Buffers<T>& input)
{
    fft(output, input);
}

template <typename T>
void Single<T>::ifft(vec::Buffers<T>& output, vec::Buffers<T>& input)
{
    output.zero_fill();
    output.copy(input, 0, 0);
}

} // namespace fft
} // namespace quadiron

#endif
