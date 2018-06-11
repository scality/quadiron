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
#ifndef __NTTEC_FFT_2_H__
#define __NTTEC_FFT_2_H__

#include "fft_base.h"
#include "gf_base.h"
#include "vec_buffers.h"
#include "vec_vector.h"

namespace nttec {
namespace fft {

/** FFT implementation specialized for FFT size=2.
 *
 *  Since w<sup>2</sup> = 1 we have either w = -1 or w = 1 (we exclude)
 *
 *    And w<sup>-1</sup> = w = -1
 */
template <typename T>
class Size2 : public FourierTransform<T> {
  public:
    explicit Size2(const gf::Field<T>& gf);
    ~Size2() = default;
    void fft(vec::Vector<T>* output, vec::Vector<T>* input) override;
    void ifft(vec::Vector<T>* output, vec::Vector<T>* input) override;
    void fft_inv(vec::Vector<T>* output, vec::Vector<T>* input) override;
    void fft(vec::Buffers<T>* output, vec::Buffers<T>* input) override;
    void ifft(vec::Buffers<T>* output, vec::Buffers<T>* input) override;
    void fft_inv(vec::Buffers<T>* output, vec::Buffers<T>* input) override;
};

template <typename T>
Size2<T>::Size2(const gf::Field<T>& gf) : FourierTransform<T>(gf, 2)
{
}

template <typename T>
void Size2<T>::fft(vec::Vector<T>* output, vec::Vector<T>* input)
{
    T a = input->get(0);
    T b = input->get(1);
    output->set(0, this->gf->add(a, b));
    output->set(1, this->gf->sub(a, b));
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void Size2<T>::fft_inv(vec::Vector<T>* output, vec::Vector<T>* input)
{
    fft(output, input);
}

template <typename T>
void Size2<T>::ifft(vec::Vector<T>* output, vec::Vector<T>* input)
{
    fft_inv(output, input);
    if (this->inv_n_mod_p > 1)
        output->mul_scalar(this->inv_n_mod_p);
}

template <typename T>
void Size2<T>::fft(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    output->copy(input);
    size_t buf_len = input->get_size();
    this->gf->add_two_bufs(input->get(1), output->get(0), buf_len);
    this->gf->sub_two_bufs(
        input->get(0), output->get(1), output->get(1), buf_len);
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void Size2<T>::fft_inv(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    fft(output, input);
}

template <typename T>
void Size2<T>::ifft(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    fft_inv(output, input);
    if (this->inv_n_mod_p > 1)
        this->gf->mul_vec_to_vecp(this->vec_inv_n, output, output);
}

} // namespace fft
} // namespace nttec

#endif
