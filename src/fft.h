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
  T inv_n;
  T w;
  T inv_w;
protected:
  FFT(GF<T> *gf, int n, T w);
public:
  virtual void fft(Vec<T> *output, Vec<T> *input) = 0;
  virtual void ifft(Vec<T> *output, Vec<T> *input) = 0;
};

template <typename T>
FFT<T>::FFT(GF<T> *gf, int n, T w)
{
  this->gf = gf;
  this->n = n;
  this->inv_n = gf->inv(n);
  this->w = w;
  this->inv_w = gf->inv(w);
}
