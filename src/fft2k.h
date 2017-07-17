
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
  FFTN<T> *fftn;
  bool bypass;
  int k;
  int N;
  T w;
  T inv_w;
  Vec<T> *W;
  Vec<T> *inv_W;
  FFT2K<T> *fftk;
 public:
  FFT2K(GF<T> *gf, int n, int R, int N = 0);
  ~FFT2K();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
 private:
  void _fft(Vec<T> *output, Vec<T> *input, bool inv);
};

/**
 * n-th root will be constructed with primitive root
 *
 * @param gf
 * @param n for now must be a power of 2
 * @param R primitive root
 * @param N original length
 *
 * @return
 */
template <typename T>
FFT2K<T>::FFT2K(GF<T> *gf, int n, int R, int N) : FFT<T>(gf, n)
{
  w = gf->get_nth_root(n, R);
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
    this->fftk = new FFT2K<T>(gf, k, R, this->N);
  } else {
    bypass = true;

    W = nullptr;
    inv_W = nullptr;

    this->fftn = new FFTN<T>(gf, n, w);
  }
}

template <typename T>
FFT2K<T>::~FFT2K()
{
  if (bypass)
    delete fftn;
  else
    delete fftk;
}

/**
 * Virtual (v->get_n()*2) x 1 vertical vector for the need of
 * Cooley-Tukey algorithm
 *
 * | v_0 |
 * | v_1 |
 * | ...
 * | v_b |
 * | v_0 |
 * | v_1 |
 * | ...
 * | v_b |
 */
template<typename T>
class V2Vec : public Vec<T>
{
 private:
  Vec<T> *vec;
 public:
  explicit V2Vec(Vec<T> *vec);
  int get_n(void);
  T get(int i);
};

template <typename T>
V2Vec<T>::V2Vec(Vec<T> *vec) : Vec<T>(NULL, 0)
{
  this->vec = vec;
}

template <typename T>
int V2Vec<T>::get_n(void)
{
  return vec->n * 2;
}

template <typename T>
T V2Vec<T>::get(int i)
{
  assert(i >= 0 && i < 2*vec->n);

  if (i < vec->n)
    return vec->get(i);
  else
    return vec->get(i - vec->n);
}

template <typename T>
void FFT2K<T>::_fft(Vec<T> *output, Vec<T> *input, bool inv)
{
  Vec<T> even(this->gf, k);
  Vec<T> odd(this->gf, k);
  for (int i = 0; i < this->n; i++) {
    if (i % 2 == 0)
      even.set(i/2, input->get(i));
    else
      odd.set(i/2, input->get(i));
  }
  // even.dump();
  // odd.dump();
  Vec<T> _even(this->gf, k);
  Vec<T> _odd(this->gf, k);
  if (inv) {
    fftk->ifft(&_even, &even);
    fftk->ifft(&_odd, &odd);
  } else {
    fftk->fft(&_even, &even);
    fftk->fft(&_odd, &odd);
  }

  V2Vec<T> veven(&_even);
  V2Vec<T> vodd(&_odd);

  if (inv)
    output->copy(inv_W, this->n);
  else
    output->copy(W, this->n);
  output->hadamard_mul(&vodd);
  output->add(&veven);
  /*
   * We need to divide output to `N` for the inverse formular
   * We need to multiply output to `2` since they are computed by fftn for n=2
   */
  if (inv && this->k == this->N / 2)
    output->mul_scalar(this->gf->div(2, this->N));
}

template <typename T>
void FFT2K<T>::fft(Vec<T> *output, Vec<T> *input)
{
  // input->dump();

  if (bypass)
    return fftn->fft(output, input);
  else
    return _fft(output, input, false);
}

template <typename T>
void FFT2K<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  if (bypass)
    return fftn->ifft(output, input);
  else
    return _fft(output, input, true);
}
