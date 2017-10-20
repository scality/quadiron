/* -*- mode: c++ -*- */

#pragma once

/*
 * DFT over the field gf and over vectors of size n with w as n-th
 * root of unity
 */
template<typename T>
class DFT
{
 protected:
  GF<T> *gf;
  int n;
  T inv_n_mod_p;
  DFT(GF<T> *gf, int n);
 public:
  virtual ~DFT();
  int get_n();
  GF<T> *get_gf();
  virtual void fft(Vec<T> *output, Vec<T> *input) = 0;
  virtual void ifft(Vec<T> *output, Vec<T> *input) = 0;
  virtual void fft_inv(Vec<T> *output, Vec<T> *input) = 0;
};

template <typename T>
DFT<T>::DFT(GF<T> *gf, int n)
{
  this->gf = gf;
  this->n = n;
  this->inv_n_mod_p = gf->inv(n) % gf->get_sub_field()->card();
}

template <typename T>
DFT<T>::~DFT()
{
}

template <typename T>
int DFT<T>::get_n()
{
  return n;
}

template <typename T>
GF<T> *DFT<T>::get_gf()
{
  return gf;
}
