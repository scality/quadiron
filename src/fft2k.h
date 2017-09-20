
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
  FFT2<T> *fft2 = NULL;
  bool bypass;
  int k;
  int N;
  T w;
  T inv_w;
  Vec<T> *W = NULL;
  Vec<T> *inv_W = NULL;
  FFT2K<T> *fftk = NULL;
  Vec<T> *even = NULL;
  Vec<T> *_even = NULL;
  Vec<T> *odd = NULL;
  Vec<T> *_odd = NULL;
  V2Vec<T> *veven = NULL;
  V2Vec<T> *vodd = NULL;
 public:
  FFT2K(GF<T> *gf, int n, int N = 0);
  ~FFT2K();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
  void fft_inv(Vec<T> *output, Vec<T> *input);
 private:
  void _fft(Vec<T> *output, Vec<T> *input, bool inv);
};

/**
 * n-th root will be constructed with primitive root
 *
 * @param gf
 * @param n for now must be a power of 2
 * @param N original length
 *
 * @return
 */
template <typename T>
FFT2K<T>::FFT2K(GF<T> *gf, int n, int N) : FFT<T>(gf, n)
{
  w = gf->get_nth_root(n);
  inv_w = gf->inv(w);
  this->N = N;
  if (this->N == 0)
    this->N = n;

  if (n != 2 && n % 2 == 0) {
    bypass = false;

    W = new Vec<T>(gf, n + 1);
    inv_W = new Vec<T>(gf, n + 1);
    gf->compute_omegas(W, n, w);
    gf->compute_omegas(inv_W, n, inv_w);

    k = n / 2;
    this->fftk = new FFT2K<T>(gf, k, this->N);

    this->even = new Vec<T>(this->gf, k);
    this->_even = new Vec<T>(this->gf, k);
    this->veven = new V2Vec<T>(this->_even);
    this->odd = new Vec<T>(this->gf, k);
    this->_odd = new Vec<T>(this->gf, k);
    this->vodd = new V2Vec<T>(this->_odd);
  } else {
    bypass = true;

    W = nullptr;
    inv_W = nullptr;

    this->fft2 = new FFT2<T>(gf);
  }
}

template <typename T>
FFT2K<T>::~FFT2K()
{
  if (bypass)
    if (fft2) delete fft2;
  else {
    if (fftk) delete fftk;
    if (W) delete W;
    if (inv_W) delete inv_W;
    if (even) delete even;
    if (_even) delete _even;
    if (odd) delete odd;
    if (_odd) delete _odd;
    if (veven) delete veven;
    if (vodd) delete vodd;
  }
}

template <typename T>
void FFT2K<T>::_fft(Vec<T> *output, Vec<T> *input, bool inv)
{
  for (int i = 0; i < this->n; i++) {
    if (i % 2 == 0)
      this->even->set(i/2, input->get(i));
    else
      this->odd->set(i/2, input->get(i));
  }
  // this->even->dump();
  // this->odd->dump();
  if (inv) {
    fftk->fft_inv(this->_even, this->even);
    fftk->fft_inv(this->_odd, this->odd);
  } else {
    fftk->fft(this->_even, this->even);
    fftk->fft(this->_odd, this->odd);
  }

  if (inv)
    output->copy(inv_W, this->n);
  else
    output->copy(W, this->n);
  output->hadamard_mul(this->vodd);
  output->add(this->veven);
}

template <typename T>
void FFT2K<T>::fft(Vec<T> *output, Vec<T> *input)
{
  // input->dump();

  if (bypass)
    return fft2->fft(output, input);
  else
    return _fft(output, input, false);
}

template <typename T>
void FFT2K<T>::fft_inv(Vec<T> *output, Vec<T> *input)
{
  if (bypass)
    return fft2->fft_inv(output, input);
  else
    return _fft(output, input, true);
}

template <typename T>
void FFT2K<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  fft_inv(output, input);

  /*
   * We need to divide output to `N` for the inverse formular
   */
  if ((this->k == this->N / 2) && (this->inv_n_mod_p > 1))
      output->mul_scalar(this->inv_n_mod_p);
}
