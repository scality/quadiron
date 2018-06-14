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
#ifndef __QUAD_FFT_2N_H__
#define __QUAD_FFT_2N_H__

#include "fft_2.h"
#include "fft_base.h"
#include "fft_single.h"
#include "gf_base.h"
#include "vec_doubled.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

namespace quadiron {
namespace fft {

/** Implementation of the radix-2 decimation-in-time (DIT) FFT
 *
 * if \f$n = 2k\f$, then
 *
 * \f[
 *   F_{2k}(A) =
 *
 *   \begin{bmatrix}
 *     F_k(A_{\text{even}}) \\
 *     F_k(A_{\text{even}})
 *   \end{bmatrix}
 *
 *   +
 *
 *   \begin{bmatrix}
 *     1 & w & w^2 \dots w^{2k-1}
 *   \end{bmatrix}
 *
 *   \times
 *
 *   \begin{bmatrix}
 *     F_k(A_{\text{odd}}) \\
 *     F_k(A_{\text{odd}})
 *   \end{bmatrix}
 * \f]
 */
template <typename T>
class Radix2 : public FourierTransform<T> {
  public:
    Radix2(
        const gf::Field<T>& gf,
        int n,
        int m = 0,
        size_t pkt_size = 0,
        int N = 0);
    ~Radix2() = default;
    void fft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void ifft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void fft_inv(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void fft(vec::Buffers<T>& output, vec::Buffers<T>& input) override;
    void ifft(vec::Buffers<T>& output, vec::Buffers<T>& input) override;
    void fft_inv(vec::Buffers<T>& output, vec::Buffers<T>& input) override;

  private:
    void _fftp(vec::Buffers<T>& output, vec::Buffers<T>& input, bool inv);
    void _fft(vec::Vector<T>& output, vec::Vector<T>& input, bool inv);

    bool bypass;
    int k;
    int m; // number of real input elements
    int N;
    T w;
    T inv_w;
    size_t pkt_size;
    std::unique_ptr<vec::Vector<T>> W = nullptr;
    std::unique_ptr<vec::Vector<T>> inv_W = nullptr;
    std::unique_ptr<vec::Vector<T>> W_half = nullptr;
    std::unique_ptr<vec::Vector<T>> inv_W_half = nullptr;

    std::unique_ptr<fft::FourierTransform<T>> fft_trivial = nullptr;
    std::unique_ptr<fft::Radix2<T>> fftk = nullptr;

    std::unique_ptr<vec::Vector<T>> even = nullptr;
    std::unique_ptr<vec::Vector<T>> _even = nullptr;
    std::unique_ptr<vec::Vector<T>> odd = nullptr;
    std::unique_ptr<vec::Vector<T>> _odd = nullptr;
    std::unique_ptr<vec::Doubled<T>> veven = nullptr;
    std::unique_ptr<vec::Doubled<T>> vodd = nullptr;

    std::unique_ptr<vec::Buffers<T>> tmp_buf = nullptr;
};

/**
 * n-th root will be constructed with primitive root
 *
 * @param gf
 * @param n FFT length, for now must be a power of 2
 * @param m length of input vector without zero padding. It allows shorterning
 *  operation cycles
 * @param pkt_size size of packet, i.e. number of symbols per chunk will be
 *  received and processed at a time
 * @param N original length
 *
 * @return
 */
template <typename T>
Radix2<T>::Radix2(const gf::Field<T>& gf, int n, int m, size_t pkt_size, int N)
    : FourierTransform<T>(gf, n)
{
    w = gf.get_nth_root(n);
    inv_w = gf.inv(w);
    this->pkt_size = pkt_size;
    this->N = N;
    if (this->N == 0) {
        this->N = n;
    }
    this->m = m > 0 ? m : n;

    // Note: we need k to calculate ifft where k is checked == N/2
    // It's necessary for FFT of length N = 2.
    k = n / 2;

    if (this->m == 1) {
        bypass = true;
        this->fft_trivial =
            std::unique_ptr<fft::Single<T>>(new fft::Single<T>(gf, this->n));
    } else if (this->n <= 2) {
        bypass = true;
        this->fft_trivial =
            std::unique_ptr<fft::Size2<T>>(new fft::Size2<T>(gf));
    } else { // (this->m > 1 && this->n > 2)
        bypass = false;

        W = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, n));
        inv_W = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, n));
        gf.compute_omegas(W.get(), n, w);
        gf.compute_omegas(inv_W.get(), n, inv_w);

        W_half = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, k));
        inv_W_half = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, k));
        for (int i = 0; i < k; i++) {
            W_half->set(i, W->get(i));
            inv_W_half->set(i, inv_W->get(i));
        }

        int next_m = this->m / 2;
        this->fftk = std::unique_ptr<Radix2<T>>(
            new Radix2<T>(gf, k, next_m, pkt_size, this->N));

        this->even = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, k));
        this->_even =
            std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, k));
        this->veven =
            std::unique_ptr<vec::Doubled<T>>(new vec::Doubled<T>(_even.get()));
        this->odd = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, k));
        this->_odd = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, k));
        this->vodd =
            std::unique_ptr<vec::Doubled<T>>(new vec::Doubled<T>(_odd.get()));

        // if (this->pkt_size > 0)
        this->tmp_buf =
            std::unique_ptr<vec::Buffers<T>>(new vec::Buffers<T>(k, pkt_size));
    }
}

template <typename T>
void Radix2<T>::_fft(vec::Vector<T>& output, vec::Vector<T>& input, bool inv)
{
    for (int i = 0; i < this->n; i++) {
        if (i % 2 == 0)
            even->set(i / 2, input.get(i));
        else
            odd->set(i / 2, input.get(i));
    }
    // this->even->dump();
    // this->odd->dump();
    if (inv) {
        fftk->fft_inv(*_even, *even);
        fftk->fft_inv(*_odd, *odd);
    } else {
        fftk->fft(*_even, *even);
        fftk->fft(*_odd, *odd);
    }

    if (inv)
        output.copy(inv_W.get(), this->n);
    else
        output.copy(W.get(), this->n);
    output.hadamard_mul(vodd.get());
    output.add(veven.get());
}

template <typename T>
void Radix2<T>::fft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    // input.dump();

    if (bypass)
        return fft_trivial->fft(output, input);
    else
        return _fft(output, input, false);
}

template <typename T>
void Radix2<T>::fft_inv(vec::Vector<T>& output, vec::Vector<T>& input)
{
    if (bypass)
        return fft_trivial->fft_inv(output, input);
    else
        return _fft(output, input, true);
}

template <typename T>
void Radix2<T>::ifft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    fft_inv(output, input);

    /*
     * We need to divide output to `N` for the inverse formular
     */
    if ((this->k == this->N / 2) && (this->inv_n_mod_p > 1))
        output.mul_scalar(this->inv_n_mod_p);
}

template <typename T>
void Radix2<T>::_fftp(vec::Buffers<T>& output, vec::Buffers<T>& input, bool inv)
{
    int half = this->n / 2;
    size_t size = input.get_size();
    std::vector<T*> even_mem(half, nullptr);
    std::vector<T*> odd_mem(half, nullptr);
    vec::Buffers<T> i_even(half, size, &even_mem);
    vec::Buffers<T> i_odd(half, size, &odd_mem);

    // separate even and odd elements of input
    input.separate_even_odd(&i_even, &i_odd);

    vec::Buffers<T> o_even(&output, 0, half);
    vec::Buffers<T> o_odd(&output, half, this->n);

    if (inv) {
        fftk->fft_inv(o_even, i_even);
        fftk->fft_inv(o_odd, i_odd);
    } else {
        fftk->fft(o_even, i_even);
        fftk->fft(o_odd, i_odd);
    }

    /*
     * output[i] = even[i] + w * odd[i] for 0 <= i < n/2
     * output[i] = even[i] - w * odd[i] otherwise
     */

    // set tmp_buf = w * o_odd
    if (inv)
        this->gf->mul_vec_to_vecp(*inv_W_half, o_odd, *tmp_buf);
    else
        this->gf->mul_vec_to_vecp(*W_half, o_odd, *tmp_buf);

    // substract o_even by tmp_buf and store in o_dd: o_even - w * o_odd
    this->gf->sub_vecp_to_vecp(&o_even, tmp_buf.get(), &o_odd);
    // add tmp_buf to o_even to get: o_even + w * o_odd
    this->gf->add_vecp_to_vecp(tmp_buf.get(), &o_even);
}

template <typename T>
void Radix2<T>::fft(vec::Buffers<T>& output, vec::Buffers<T>& input)
{
    if (bypass)
        return fft_trivial->fft(output, input);
    else
        return _fftp(output, input, false);
}

template <typename T>
void Radix2<T>::fft_inv(vec::Buffers<T>& output, vec::Buffers<T>& input)
{
    if (bypass)
        return fft_trivial->fft_inv(output, input);
    else
        return _fftp(output, input, true);
}

template <typename T>
void Radix2<T>::ifft(vec::Buffers<T>& output, vec::Buffers<T>& input)
{
    fft_inv(output, input);

    /*
     * We need to divide output to `N` for the inverse formular
     */
    if ((this->k == this->N / 2) && (this->inv_n_mod_p > 1))
        this->gf->mul_vec_to_vecp(*(this->vec_inv_n), output, output);
}

} // namespace fft
} // namespace quadiron

#endif
