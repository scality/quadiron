
/* -*- mode: c++ -*- */

#pragma once

/** 
 * Cooley-Tukey Algorithm
 *
 * if n = 2k, then
 *
 *           | F_k(A_even) |   |     1    |   | F_k(A_odd) |
 * F_2k(A) = | F_k(A_even) | + |     w    | . | F_k(A_odd) |
 *                             |    w^2   |
 *                             |    ...   |
 *                             | w^(2k-1) |
 *
 * @param output 
 * @param input 
 */
template<typename T>
class FFT2K : public FFT<T>
{
 private:

 public:
  FFT2K(GF<T> *gf, int n, T w);
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
};

template <typename T>
FFT2K<T>::FFT2K(GF<T> *gf, int n, T w) : FFT<T>(gf, n, w)
{
}

template <typename T>
void FFT2K<T>::fft(Vec<T> *output, Vec<T> *input)
{
}

template <typename T>
void FFT2K<T>::ifft(Vec<T> *output, Vec<T> *input)
{
}
