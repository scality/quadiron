/* -*- mode: c++ -*- */
#ifndef __NTL_FFT2_H__
#define __NTL_FFT2_H__

#include "dft.h"
#include "vec.h"

/**
 * Algorithm for FFT(n=2)
 *  Since w^2 = 1 => w = -1 or w = 1 (we exclude)
 *    And w^(-1) = w = -1;
 */
template<typename T>
class FFT2 : public DFT<T>
{
 public:
  explicit FFT2(GF<T> *gf);
  ~FFT2();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
  void fft_inv(Vec<T> *output, Vec<T> *input);
  void fft(Vecp<T> *output, Vecp<T> *input);
  void ifft(Vecp<T> *output, Vecp<T> *input);
  void fft_inv(Vecp<T> *output, Vecp<T> *input);
};

template <typename T>
FFT2<T>::FFT2(GF<T> *gf) : DFT<T>(gf, 2)
{}

template <typename T>
FFT2<T>::~FFT2()
{}

template <typename T>
void FFT2<T>::fft(Vec<T> *output, Vec<T> *input)
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
void FFT2<T>::fft_inv(Vec<T> *output, Vec<T> *input)
{
  fft(output, input);
}

template <typename T>
void FFT2<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  fft_inv(output, input);
  if (this->inv_n_mod_p > 1)
    output->mul_scalar(this->inv_n_mod_p);
}

template <typename T>
void FFT2<T>::fft(Vecp<T> *output, Vecp<T> *input)
{
  output->copy(input);
  size_t buf_len = input->get_size();
  this->gf->add_two_bufs(input->get(1), output->get(0) , buf_len);
  this->gf->sub_two_bufs(input->get(0), output->get(1), output->get(1),
    buf_len);
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void FFT2<T>::fft_inv(Vecp<T> *output, Vecp<T> *input)
{
  fft(output, input);
}

template <typename T>
void FFT2<T>::ifft(Vecp<T> *output, Vecp<T> *input)
{
  fft_inv(output, input);
  if (this->inv_n_mod_p > 1)
    this->gf->mul_vec_to_vecp(this->vec_inv_n, output, output);
}

#endif
