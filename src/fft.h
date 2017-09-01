/* -*- mode: c++ -*- */

#pragma once

/*
 * FFT over the field gf and over vectors of size n with w as n-th
 * root of unity
 */
template<typename T>
class FFT
{
 public:
  GF<T> *gf;
  int n;
  T inv_n_mod_p;
 protected:
  FFT(GF<T> *gf, int n);
 public:
  virtual void fft(Vec<T> *output, Vec<T> *input) = 0;
  virtual void ifft(Vec<T> *output, Vec<T> *input) = 0;
};

template <typename T>
FFT<T>::FFT(GF<T> *gf, int n)
{
  this->gf = gf;
  this->n = n;
  this->inv_n_mod_p = (gf->inv(n) % gf->p);
}
