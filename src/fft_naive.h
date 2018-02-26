/* -*- mode: c++ -*- */
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
