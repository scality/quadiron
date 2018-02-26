/* -*- mode: c++ -*- */
#ifndef __NTTEC_FFT_2N_H__
#define __NTTEC_FFT_2N_H__

#include "fft_2.h"
#include "fft_base.h"
#include "fft_single.h"
#include "gf_base.h"
#include "vec_doubled.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

namespace nttec {
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
    Radix2(gf::Field<T>* gf, int n, int m = 0, size_t pkt_size = 0, int N = 0);
    ~Radix2();
    void fft(vec::Vector<T>* output, vec::Vector<T>* input);
    void ifft(vec::Vector<T>* output, vec::Vector<T>* input);
    void fft_inv(vec::Vector<T>* output, vec::Vector<T>* input);
    void fft(vec::Buffers<T>* output, vec::Buffers<T>* input);
    void ifft(vec::Buffers<T>* output, vec::Buffers<T>* input);
    void fft_inv(vec::Buffers<T>* output, vec::Buffers<T>* input);

  private:
    void _fftp(vec::Buffers<T>* output, vec::Buffers<T>* input, bool inv);
    void _fft(vec::Vector<T>* output, vec::Vector<T>* input, bool inv);

    bool bypass;
    int k;
    int m; // number of real input elements
    int N;
    T w;
    T inv_w;
    size_t pkt_size;
    vec::Vector<T>* W = nullptr;
    vec::Vector<T>* inv_W = nullptr;
    vec::Vector<T>* W_half = nullptr;
    vec::Vector<T>* inv_W_half = nullptr;

    fft::FourierTransform<T>* fft_trivial = nullptr;
    fft::Radix2<T>* fftk = nullptr;

    vec::Vector<T>* even = nullptr;
    vec::Vector<T>* _even = nullptr;
    vec::Vector<T>* odd = nullptr;
    vec::Vector<T>* _odd = nullptr;
    vec::Doubled<T>* veven = nullptr;
    vec::Doubled<T>* vodd = nullptr;

    vec::Buffers<T>* tmp_buf = nullptr;
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
Radix2<T>::Radix2(gf::Field<T>* gf, int n, int m, size_t pkt_size, int N)
    : FourierTransform<T>(gf, n)
{
    w = gf->get_nth_root(n);
    inv_w = gf->inv(w);
    this->pkt_size = pkt_size;
    this->N = N;
    if (this->N == 0) {
        this->N = n;
    }
    this->m = m > 0 ? m : n;

    if (this->m == 1) {
        bypass = true;
        this->fft_trivial = new fft::Single<T>(gf, this->n);
    } else if (this->n <= 2) {
        bypass = true;
        this->fft_trivial = new fft::Size2<T>(gf);
    } else { // (this->m > 1 && this->n > 2)
        bypass = false;
        k = n / 2;

        W = new vec::Vector<T>(gf, n);
        inv_W = new vec::Vector<T>(gf, n);
        gf->compute_omegas(W, n, w);
        gf->compute_omegas(inv_W, n, inv_w);

        W_half = new vec::Vector<T>(gf, k);
        inv_W_half = new vec::Vector<T>(gf, k);
        for (int i = 0; i < k; i++) {
            W_half->set(i, W->get(i));
            inv_W_half->set(i, inv_W->get(i));
        }

        int next_m = this->m / 2;
        this->fftk = new Radix2<T>(gf, k, next_m, pkt_size, this->N);

        this->even = new vec::Vector<T>(this->gf, k);
        this->_even = new vec::Vector<T>(this->gf, k);
        this->veven = new vec::Doubled<T>(this->_even);
        this->odd = new vec::Vector<T>(this->gf, k);
        this->_odd = new vec::Vector<T>(this->gf, k);
        this->vodd = new vec::Doubled<T>(this->_odd);

        if (this->pkt_size > 0)
            this->tmp_buf = new vec::Buffers<T>(k, pkt_size);
    }
}

template <typename T>
Radix2<T>::~Radix2()
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
void Radix2<T>::_fft(vec::Vector<T>* output, vec::Vector<T>* input, bool inv)
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
void Radix2<T>::fft(vec::Vector<T>* output, vec::Vector<T>* input)
{
    // input->dump();

    if (bypass)
        return fft_trivial->fft(output, input);
    else
        return _fft(output, input, false);
}

template <typename T>
void Radix2<T>::fft_inv(vec::Vector<T>* output, vec::Vector<T>* input)
{
    if (bypass)
        return fft_trivial->fft_inv(output, input);
    else
        return _fft(output, input, true);
}

template <typename T>
void Radix2<T>::ifft(vec::Vector<T>* output, vec::Vector<T>* input)
{
    fft_inv(output, input);

    /*
     * We need to divide output to `N` for the inverse formular
     */
    if ((this->k == this->N / 2) && (this->inv_n_mod_p > 1))
        output->mul_scalar(this->inv_n_mod_p);
}

template <typename T>
void Radix2<T>::_fftp(vec::Buffers<T>* output, vec::Buffers<T>* input, bool inv)
{
    int half = this->n / 2;
    size_t size = input->get_size();
    std::vector<T*> even_mem(half, nullptr);
    std::vector<T*> odd_mem(half, nullptr);
    vec::Buffers<T> i_even(half, size, &even_mem);
    vec::Buffers<T> i_odd(half, size, &odd_mem);

    // separate even and odd elements of input
    input->separate_even_odd(&i_even, &i_odd);

    vec::Buffers<T> o_even(output, 0, half);
    vec::Buffers<T> o_odd(output, half, this->n);

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
    vec::Buffers<T>* _o_odd = nullptr;
    if (this->pkt_size == 0) {
        _o_odd = new vec::Buffers<T>(half, size);
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
void Radix2<T>::fft(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    if (bypass)
        return fft_trivial->fft(output, input);
    else
        return _fftp(output, input, false);
}

template <typename T>
void Radix2<T>::fft_inv(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    if (bypass)
        return fft_trivial->fft_inv(output, input);
    else
        return _fftp(output, input, true);
}

template <typename T>
void Radix2<T>::ifft(vec::Buffers<T>* output, vec::Buffers<T>* input)
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
