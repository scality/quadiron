/* -*- mode: c++ -*- */
#ifndef __NTL_DFTN_H__
#define __NTL_DFTN_H__

#include "dft.h"
#include "mat.h"

/**
 * Algorithm for very small n
 */
template<typename T>
class DFTN : public DFT<T>
{
 private:
  T w;
  T inv_w;
  Mat<T> *W;
  Mat<T> *inv_W;
  void compute_W(Mat<T> *_W, T _w);
  void _fft(Vec<T> *output, Vec<T> *input, Mat<T> *_W);
 public:
  DFTN(GF<T> *gf, int n, T w);
  ~DFTN();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
  void fft_inv(Vec<T> *output, Vec<T> *input);
};

template <typename T>
DFTN<T>::DFTN(GF<T> *gf, int n, T w) : DFT<T>(gf, n)
{
  this->w = w;
  this->inv_w = gf->inv(w);
  this->W = new Mat<T>(gf, this->n, this->n);
  this->inv_W = new Mat<T>(gf, this->n, this->n);

  compute_W(W, w);
  compute_W(inv_W, this->inv_w);
}

template <typename T>
DFTN<T>::~DFTN()
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
void DFTN<T>::compute_W(Mat<T> *_W, T _w)
{
  for (int i = 0; i <= this->n-1; i++) {
    for (int j = 0; j <= this->n-1; j++) {
      _W->set(i, j, this->gf->exp(_w, (i*j) % this->n));
    }
  }
  // _W->dump();
}

template <typename T>
void DFTN<T>::_fft(Vec<T> *output, Vec<T> *input, Mat<T> *_W)
{
  _W->mul(output, input);
}

template <typename T>
void DFTN<T>::fft(Vec<T> *output, Vec<T> *input)
{
  _fft(output, input, W);
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void DFTN<T>::fft_inv(Vec<T> *output, Vec<T> *input)
{
  _fft(output, input, inv_W);
}

template <typename T>
void DFTN<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  _fft(output, input, inv_W);
  if (this->inv_n_mod_p > 1)
    output->mul_scalar(this->inv_n_mod_p);
}

#endif
