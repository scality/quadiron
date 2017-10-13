/* -*- mode: c++ -*- */

#pragma once

/*
 * FFT over the field gf and over vectors of size n with w as n-th
 * root of unity
 */
template<typename T>
class FFT
{
 protected:
  GF<T> *gf;
  int n;
  T inv_n_mod_p;
  FFT(GF<T> *gf, int n);
 public:
  virtual ~FFT();
  int get_n();
  GF<T> *get_gf();
  virtual void fft(Vec<T> *output, Vec<T> *input) = 0;
  virtual void ifft(Vec<T> *output, Vec<T> *input) = 0;
  virtual void fft_inv(Vec<T> *output, Vec<T> *input) = 0;
};

template <typename T>
FFT<T>::FFT(GF<T> *gf, int n)
{
  this->gf = gf;
  this->n = n;
  this->inv_n_mod_p = gf->inv(n) % gf->get_sub_field()->card();
}

template <typename T>
FFT<T>::~FFT()
{
}

template <typename T>
int FFT<T>::get_n()
{
  return n;
}

template <typename T>
GF<T> *FFT<T>::get_gf()
{
  return gf;
}

