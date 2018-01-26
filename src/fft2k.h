
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
class FFT2K : public DFT<T>
{
 private:
  bool bypass;
  int k;
  int N;
  T w;
  T inv_w;
  Vec<T> *W = nullptr;
  Vec<T> *W_half = nullptr;
  Vec<T> *inv_W = nullptr;
  Vec<T> *inv_W_half = nullptr;

  FFT2<T> *fft2 = nullptr;
  FFT2K<T> *fftk = nullptr;

  Vec<T> *even = nullptr;
  Vec<T> *_even = nullptr;
  Vec<T> *odd = nullptr;
  Vec<T> *_odd = nullptr;
  V2Vec<T> *veven = nullptr;
  V2Vec<T> *vodd = nullptr;
 public:
  FFT2K(GF<T> *gf, int n, int N = 0);
  ~FFT2K();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
  void fft_inv(Vec<T> *output, Vec<T> *input);
  void fft(Vecp<T> *output, Vecp<T> *input);
  void ifft(Vecp<T> *output, Vecp<T> *input);
  void fft_inv(Vecp<T> *output, Vecp<T> *input);
 private:
  void _fftp(Vecp<T> *output, Vecp<T> *input, bool inv);
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
FFT2K<T>::FFT2K(GF<T> *gf, int n, int N) : DFT<T>(gf, n)
{
  w = gf->get_nth_root(n);
  inv_w = gf->inv(w);
  this->N = N;
  if (this->N == 0)
    this->N = n;

  if (n != 2 && n % 2 == 0) {
    bypass = false;
    k = n / 2;

    W = new Vec<T>(gf, n);
    inv_W = new Vec<T>(gf, n);
    gf->compute_omegas(W, n, w);
    gf->compute_omegas(inv_W, n, inv_w);

    W_half = new Vec<T>(gf, k);
    inv_W_half = new Vec<T>(gf, k);
    for (int i = 0; i < k; i++) {
      W_half->set(i, W->get(i));
      inv_W_half->set(i, inv_W->get(i));
    }

    this->fftk = new FFT2K<T>(gf, k, this->N);

    this->even = new Vec<T>(this->gf, k);
    this->_even = new Vec<T>(this->gf, k);
    this->veven = new V2Vec<T>(this->_even);
    this->odd = new Vec<T>(this->gf, k);
    this->_odd = new Vec<T>(this->gf, k);
    this->vodd = new V2Vec<T>(this->_odd);
  } else {
    bypass = true;
    this->fft2 = new FFT2<T>(gf);
  }
}

template <typename T>
FFT2K<T>::~FFT2K()
{
  if (bypass) {
    if (fft2 != nullptr) delete fft2;
  } else {
    if (fftk != nullptr) delete fftk;
    if (W != nullptr) delete W;
    if (inv_W != nullptr) delete inv_W;
    if (W_half != nullptr) delete W_half;
    if (inv_W_half != nullptr) delete inv_W_half;
    if (even != nullptr) delete even;
    if (_even != nullptr) delete _even;
    if (odd != nullptr) delete odd;
    if (_odd != nullptr) delete _odd;
    if (veven != nullptr) delete veven;
    if (vodd != nullptr) delete vodd;
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

template <typename T>
void FFT2K<T>::_fftp(Vecp<T> *output, Vecp<T> *input, bool inv)
{
  int half = this->n / 2;
  size_t size = input->get_size();
  std::vector<T*> even_mem(half, nullptr);
  std::vector<T*> odd_mem(half, nullptr);
  Vecp<T> i_even(half, size, &even_mem);
  Vecp<T> i_odd(half, size, &odd_mem);

  // separate even and odd elements of input
  input->separate_even_odd(&i_even, &i_odd);

  Vecp<T> *o_even = output->slice(0, half);
  Vecp<T> *o_odd = output->slice(half, this->n);

  if (inv) {
    fftk->fft_inv(o_even, &i_even);
    fftk->fft_inv(o_odd, &i_odd);
  } else {
    fftk->fft(o_even, &i_even);
    fftk->fft(o_odd, &i_odd);
  }

  /*
   * output[i] = even[i] + w * odd[i] for 0 <= i < n/2
   * output[i] = even[i] - w * odd[i] otherwise
   */
  // tmp vec to store o_odd
  Vecp<T> _o_odd(o_odd);

  // multiply _o_odd by w or inv_w
  if (inv)
    this->gf->mul_vec_to_vecp(inv_W_half, &_o_odd);
  else
    this->gf->mul_vec_to_vecp(W_half, &_o_odd);
  // set o_odd = o_even to set two halves of output = o_even
  o_odd->copy(o_even);
  // add _o_odd to o_even to get: even + w * odd
  this->gf->add_vecp_to_vecp(&_o_odd, o_even);
  // substract o_odd by _o_odd to get: even - w * odd
  this->gf->sub_vecp_to_vecp(o_odd, &_o_odd, o_odd);

  delete o_even;
  delete o_odd;
}

template <typename T>
void FFT2K<T>::fft(Vecp<T> *output, Vecp<T> *input)
{
  if (bypass)
    return fft2->fft(output, input);
  else
    return _fftp(output, input, false);
}

template <typename T>
void FFT2K<T>::fft_inv(Vecp<T> *output, Vecp<T> *input)
{
  if (bypass)
    return fft2->fft_inv(output, input);
  else
    return _fftp(output, input, true);
}

template <typename T>
void FFT2K<T>::ifft(Vecp<T> *output, Vecp<T> *input)
{
  fft_inv(output, input);

  /*
   * We need to divide output to `N` for the inverse formular
   */
  if ((this->k == this->N / 2) && (this->inv_n_mod_p > 1))
      this->gf->mul_vec_to_vecp(this->vec_inv_n, output);
}
