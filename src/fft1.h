/* -*- mode: c++ -*- */
#ifndef __NTTEC_FFT1_H__
#define __NTTEC_FFT1_H__

#include "dft.h"
#include "vec.h"

namespace nttec {
namespace fft {

/**
 * Algorithm for FFT(n) where input's elements except the first one are all zero
 * All elements of output equals the first element of input
 */
template <typename T>
class FFT1 : public DFT<T> {
  public:
    explicit FFT1(gf::GF<T>* gf, int n);
    ~FFT1();
    void fft(vec::Vec<T>* output, vec::Vec<T>* input);
    void ifft(vec::Vec<T>* output, vec::Vec<T>* input);
    void fft_inv(vec::Vec<T>* output, vec::Vec<T>* input);
    void fft(vec::Vecp<T>* output, vec::Vecp<T>* input);
    void ifft(vec::Vecp<T>* output, vec::Vecp<T>* input);
    void fft_inv(vec::Vecp<T>* output, vec::Vecp<T>* input);
};

template <typename T>
FFT1<T>::FFT1(gf::GF<T>* gf, int n) : DFT<T>(gf, n)
{
}

template <typename T>
FFT1<T>::~FFT1()
{
}

template <typename T>
void FFT1<T>::fft(vec::Vec<T>* output, vec::Vec<T>* input)
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
void FFT1<T>::fft_inv(vec::Vec<T>* output, vec::Vec<T>* input)
{
    fft(output, input);
}

template <typename T>
void FFT1<T>::ifft(vec::Vec<T>* output, vec::Vec<T>* input)
{
    output->zero_fill();
    output->set(0, input->get(0));
}

template <typename T>
void FFT1<T>::fft(vec::Vecp<T>* output, vec::Vecp<T>* input)
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
void FFT1<T>::fft_inv(vec::Vecp<T>* output, vec::Vecp<T>* input)
{
    fft(output, input);
}

template <typename T>
void FFT1<T>::ifft(vec::Vecp<T>* output, vec::Vecp<T>* input)
{
    output->zero_fill();
    output->copy(0, input->get(0));
}

} // namespace fft
} // namespace nttec

#endif
