/* -*- mode: c++ -*- */
#ifndef __NTTEC_FFT2K_H__
#define __NTTEC_FFT2K_H__

#include "dft.h"
#include "fft1.h"
#include "fft2.h"

namespace nttec {
namespace fft {

/**
 * Cooley-Tukey Algorithm
 *
 * if n = 2k, then
 *
 *           | F_k(A_even) |   |     1    |   | F_k(A_odd) |
 * F_2k(A) = | F_k(A_even) | + |     w    | . | F_k(A_odd) |
 *                             |    w^2   |
 *                             |    ...   |
 *                             | w^(2k-1) |
 *
 * @param output
 * @param input
 */
template <typename T>
class FFT2K : public DFT<T> {
  private:
    bool bypass;
    int k;
    int m; // number of real input elements
    int N;
    T w;
    T inv_w;
    size_t pkt_size;
    vec::Vec<T>* W = nullptr;
    vec::Vec<T>* inv_W = nullptr;
    vec::Vec<T>* W_half = nullptr;
    vec::Vec<T>* inv_W_half = nullptr;

    DFT<T>* fft_trivial = nullptr;
    FFT2K<T>* fftk = nullptr;

    vec::Vec<T>* even = nullptr;
    vec::Vec<T>* _even = nullptr;
    vec::Vec<T>* odd = nullptr;
    vec::Vec<T>* _odd = nullptr;
    vec::V2Vec<T>* veven = nullptr;
    vec::V2Vec<T>* vodd = nullptr;

    vec::Vecp<T>* tmp_buf = nullptr;

  public:
    FFT2K(gf::GF<T>* gf, int n, int m = 0, size_t pkt_size = 0, int N = 0);
    ~FFT2K();
    void fft(vec::Vec<T>* output, vec::Vec<T>* input);
    void ifft(vec::Vec<T>* output, vec::Vec<T>* input);
    void fft_inv(vec::Vec<T>* output, vec::Vec<T>* input);
    void fft(vec::Vecp<T>* output, vec::Vecp<T>* input);
    void ifft(vec::Vecp<T>* output, vec::Vecp<T>* input);
    void fft_inv(vec::Vecp<T>* output, vec::Vecp<T>* input);

  private:
    void _fftp(vec::Vecp<T>* output, vec::Vecp<T>* input, bool inv);
    void _fft(vec::Vec<T>* output, vec::Vec<T>* input, bool inv);
};

/**
 * n-th root will be constructed with primitive root
 *
 * @param gf
 * @param n for now must be a power of 2
 * @param pkt_size size of packet, i.e. number of symbols per chunk will be
 *  received and processed at a time
 * @param N original length
 *
 * @return
 */
template <typename T>
FFT2K<T>::FFT2K(gf::GF<T>* gf, int n, int m, size_t pkt_size, int N)
    : DFT<T>(gf, n)
{
    w = gf->get_nth_root(n);
    inv_w = gf->inv(w);
    this->pkt_size = pkt_size;
    this->N = N;
    if (this->N == 0)
        this->N = n;
    if (m > 0)
        this->m = m;
    else
        this->m = n;

    if (this->m == 1) {
        bypass = true;
        this->fft_trivial = new FFT1<T>(gf, this->n);
    } else if (this->n <= 2) {
        bypass = true;
        this->fft_trivial = new FFT2<T>(gf);
    } else { // (this->m > 1 && this->n > 2)
        bypass = false;
        k = n / 2;

        W = new vec::Vec<T>(gf, n);
        inv_W = new vec::Vec<T>(gf, n);
        gf->compute_omegas(W, n, w);
        gf->compute_omegas(inv_W, n, inv_w);

        W_half = new vec::Vec<T>(gf, k);
        inv_W_half = new vec::Vec<T>(gf, k);
        for (int i = 0; i < k; i++) {
            W_half->set(i, W->get(i));
            inv_W_half->set(i, inv_W->get(i));
        }

        int next_m = this->m / 2;
        this->fftk = new FFT2K<T>(gf, k, next_m, pkt_size, this->N);

        this->even = new vec::Vec<T>(this->gf, k);
        this->_even = new vec::Vec<T>(this->gf, k);
        this->veven = new vec::V2Vec<T>(this->_even);
        this->odd = new vec::Vec<T>(this->gf, k);
        this->_odd = new vec::Vec<T>(this->gf, k);
        this->vodd = new vec::V2Vec<T>(this->_odd);

        if (this->pkt_size > 0)
            this->tmp_buf = new vec::Vecp<T>(k, pkt_size);
    }
}

template <typename T>
FFT2K<T>::~FFT2K()
{
    if (bypass) {
        if (fft_trivial != nullptr)
            delete fft_trivial;
    } else {
        if (fftk != nullptr)
            delete fftk;
        if (W != nullptr)
            delete W;
        if (inv_W != nullptr)
            delete inv_W;
        if (W_half != nullptr)
            delete W_half;
        if (inv_W_half != nullptr)
            delete inv_W_half;
        if (even != nullptr)
            delete even;
        if (_even != nullptr)
            delete _even;
        if (odd != nullptr)
            delete odd;
        if (_odd != nullptr)
            delete _odd;
        if (veven != nullptr)
            delete veven;
        if (vodd != nullptr)
            delete vodd;
        if (pkt_size > 0 && tmp_buf != nullptr)
            delete tmp_buf;
    }
}

template <typename T>
void FFT2K<T>::_fft(vec::Vec<T>* output, vec::Vec<T>* input, bool inv)
{
    for (int i = 0; i < this->n; i++) {
        if (i % 2 == 0)
            this->even->set(i / 2, input->get(i));
        else
            this->odd->set(i / 2, input->get(i));
    }
    // this->even->dump();
    // this->odd->dump();
    if (inv) {
        fftk->fft_inv(this->_even, this->even);
        fftk->fft_inv(this->_odd, this->odd);
    } else {
        fftk->fft(this->_even, this->even);
        fftk->fft(this->_odd, this->odd);
    }

    if (inv)
        output->copy(inv_W, this->n);
    else
        output->copy(W, this->n);
    output->hadamard_mul(this->vodd);
    output->add(this->veven);
}

template <typename T>
void FFT2K<T>::fft(vec::Vec<T>* output, vec::Vec<T>* input)
{
    // input->dump();

    if (bypass)
        return fft_trivial->fft(output, input);
    else
        return _fft(output, input, false);
}

template <typename T>
void FFT2K<T>::fft_inv(vec::Vec<T>* output, vec::Vec<T>* input)
{
    if (bypass)
        return fft_trivial->fft_inv(output, input);
    else
        return _fft(output, input, true);
}

template <typename T>
void FFT2K<T>::ifft(vec::Vec<T>* output, vec::Vec<T>* input)
{
    fft_inv(output, input);

    /*
     * We need to divide output to `N` for the inverse formular
     */
    if ((this->k == this->N / 2) && (this->inv_n_mod_p > 1))
        output->mul_scalar(this->inv_n_mod_p);
}

template <typename T>
void FFT2K<T>::_fftp(vec::Vecp<T>* output, vec::Vecp<T>* input, bool inv)
{
    int half = this->n / 2;
    size_t size = input->get_size();
    std::vector<T*> even_mem(half, nullptr);
    std::vector<T*> odd_mem(half, nullptr);
    vec::Vecp<T> i_even(half, size, &even_mem);
    vec::Vecp<T> i_odd(half, size, &odd_mem);

    // separate even and odd elements of input
    input->separate_even_odd(&i_even, &i_odd);

    vec::Vecp<T> o_even(output, 0, half);
    vec::Vecp<T> o_odd(output, half, this->n);

    if (inv) {
        fftk->fft_inv(&o_even, &i_even);
        fftk->fft_inv(&o_odd, &i_odd);
    } else {
        fftk->fft(&o_even, &i_even);
        fftk->fft(&o_odd, &i_odd);
    }

    /*
     * output[i] = even[i] + w * odd[i] for 0 <= i < n/2
     * output[i] = even[i] - w * odd[i] otherwise
     */
    // prepare tm_buf
    vec::Vecp<T>* _o_odd = nullptr;
    if (this->pkt_size == 0) {
        _o_odd = new vec::Vecp<T>(half, size);
        tmp_buf = _o_odd;
    }

    // set tmp_buf = w * o_odd
    if (inv)
        this->gf->mul_vec_to_vecp(inv_W_half, &o_odd, tmp_buf);
    else
        this->gf->mul_vec_to_vecp(W_half, &o_odd, tmp_buf);

    // substract o_even by tmp_buf and store in o_dd: o_even - w * o_odd
    this->gf->sub_vecp_to_vecp(&o_even, tmp_buf, &o_odd);
    // add tmp_buf to o_even to get: o_even + w * o_odd
    this->gf->add_vecp_to_vecp(tmp_buf, &o_even);

    if (_o_odd != nullptr)
        delete _o_odd;
}

template <typename T>
void FFT2K<T>::fft(vec::Vecp<T>* output, vec::Vecp<T>* input)
{
    if (bypass)
        return fft_trivial->fft(output, input);
    else
        return _fftp(output, input, false);
}

template <typename T>
void FFT2K<T>::fft_inv(vec::Vecp<T>* output, vec::Vecp<T>* input)
{
    if (bypass)
        return fft_trivial->fft_inv(output, input);
    else
        return _fftp(output, input, true);
}

template <typename T>
void FFT2K<T>::ifft(vec::Vecp<T>* output, vec::Vecp<T>* input)
{
    fft_inv(output, input);

    /*
     * We need to divide output to `N` for the inverse formular
     */
    if ((this->k == this->N / 2) && (this->inv_n_mod_p > 1))
        this->gf->mul_vec_to_vecp(this->vec_inv_n, output, output);
}

} // namespace fft
} // namespace nttec

#endif
