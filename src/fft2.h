/* -*- mode: c++ -*- */
#ifndef __NTTEC_FFT2_H__
#define __NTTEC_FFT2_H__

#include "dft.h"
#include "vec.h"

namespace nttec {
namespace fft {

/**
 * Algorithm for FFT(n=2)
 *  Since w^2 = 1 => w = -1 or w = 1 (we exclude)
 *    And w^(-1) = w = -1;
 */
template <typename T>
class FFT2 : public DFT<T> {
  public:
    explicit FFT2(gf::GF<T>* gf);
    ~FFT2();
    void fft(vec::Vec<T>* output, vec::Vec<T>* input);
    void ifft(vec::Vec<T>* output, vec::Vec<T>* input);
    void fft_inv(vec::Vec<T>* output, vec::Vec<T>* input);
    void fft(vec::Vecp<T>* output, vec::Vecp<T>* input);
    void ifft(vec::Vecp<T>* output, vec::Vecp<T>* input);
    void fft_inv(vec::Vecp<T>* output, vec::Vecp<T>* input);
};

template <typename T>
FFT2<T>::FFT2(gf::GF<T>* gf) : DFT<T>(gf, 2)
{
}

template <typename T>
FFT2<T>::~FFT2()
{
}

template <typename T>
void FFT2<T>::fft(vec::Vec<T>* output, vec::Vec<T>* input)
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
void FFT2<T>::fft_inv(vec::Vec<T>* output, vec::Vec<T>* input)
{
    fft(output, input);
}

template <typename T>
void FFT2<T>::ifft(vec::Vec<T>* output, vec::Vec<T>* input)
{
    fft_inv(output, input);
    if (this->inv_n_mod_p > 1)
        output->mul_scalar(this->inv_n_mod_p);
}

template <typename T>
void FFT2<T>::fft(vec::Vecp<T>* output, vec::Vecp<T>* input)
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
void FFT2<T>::fft_inv(vec::Vecp<T>* output, vec::Vecp<T>* input)
{
    fft(output, input);
}

template <typename T>
void FFT2<T>::ifft(vec::Vecp<T>* output, vec::Vecp<T>* input)
{
    fft_inv(output, input);
    if (this->inv_n_mod_p > 1)
        this->gf->mul_vec_to_vecp(this->vec_inv_n, output, output);
}

} // namespace fft
} // namespace nttec

#endif
