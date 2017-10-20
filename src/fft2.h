/* -*- mode: c++ -*- */

#pragma once

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
  output->set(0, this->gf->add(input->get(0), input->get(1)));
  output->set(1, this->gf->sub(input->get(0), input->get(1)));
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
