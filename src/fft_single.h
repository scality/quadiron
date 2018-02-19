/* -*- mode: c++ -*- */
#ifndef __NTTEC_FFT_SINGLE_H__
#define __NTTEC_FFT_SINGLE_H__

#include "fft_base.h"
#include "gf_base.h"
#include "vec_buffers.h"
#include "vec_vector.h"

namespace nttec {
namespace fft {

/** FFT specialization for when all but the first element are zeros.
 *
 * In such case, all elements of output are equals the first element of input.
 */
template <typename T>
class Single : public FourierTransform<T> {
  public:
    explicit Single(gf::Field<T>* gf, int n);
    ~Single();
    void fft(vec::Vector<T>* output, vec::Vector<T>* input);
    void ifft(vec::Vector<T>* output, vec::Vector<T>* input);
    void fft_inv(vec::Vector<T>* output, vec::Vector<T>* input);
    void fft(vec::Buffers<T>* output, vec::Buffers<T>* input);
    void ifft(vec::Buffers<T>* output, vec::Buffers<T>* input);
    void fft_inv(vec::Buffers<T>* output, vec::Buffers<T>* input);
};

template <typename T>
Single<T>::Single(gf::Field<T>* gf, int n) : FourierTransform<T>(gf, n)
{
}

template <typename T>
Single<T>::~Single()
{
}

template <typename T>
void Single<T>::fft(vec::Vector<T>* output, vec::Vector<T>* input)
{
    T val = input->get(0);
    output->fill(val);
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void Single<T>::fft_inv(vec::Vector<T>* output, vec::Vector<T>* input)
{
    fft(output, input);
}

template <typename T>
void Single<T>::ifft(vec::Vector<T>* output, vec::Vector<T>* input)
{
    output->zero_fill();
    output->set(0, input->get(0));
}

template <typename T>
void Single<T>::fft(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    T* buf = input->get(0);
    for (int i = 0; i < this->n; i++)
        output->copy(i, buf);
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void Single<T>::fft_inv(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    fft(output, input);
}

template <typename T>
void Single<T>::ifft(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    output->zero_fill();
    output->copy(0, input->get(0));
}

} // namespace fft
} // namespace nttec

#endif
