
/* -*- mode: c++ -*- */

#pragma once

/**
 * Cooley-Tukey FFT algorithm
 *  Implementation based on S.K. Matra's slides
 *    http://sip.cua.edu/res/docs/courses/ee515/chapter08/ch8-2.pdf
 *
 * Suppose n = n1*n2
 *  X[k] = sum_i x_i * w_n^ik is plit into two smaller DFTs
 * by using the index mapping
 *  i = i1 + n1*i2
 *  k = n2*k1 + k2
 *  where 1 <= i1, k1 <= n1-1, 1 <= i2, k2 <= n2-1
 *
 *  X[n2*k1 + k2] =
 *    sum_i1 [
 *      (
 *        sum_i2 x[i1 + n1*i2] * w_2^(i2*k2)
 *      ) * w^(i1*k2)
 *    ] w_1^(i1*k1)
 *
 *  where
 *    w is nth root, w_1 = w^n2, w_2 = w^n1
 *
 *  Step1: calculate DFT of the inner parenthese, i.e. sum_i2...
 *  Step2: multiply to twiddle factors w^(i1*k2)
 *  Step3: calculate DFT of the outer parenthese, i.e. sum_i1
 *
 * @param output
 * @param input
 */
template<typename T>
class FFTCT : public DFT<T>
{
 private:
  bool loop;
  bool first_layer_fft;
  T n1;
  T n2;
  T w, w1, w2, inv_w;
  Vec<T> *G = NULL;
  VcVec<T> *Y = NULL;
  VcVec<T> *X = NULL;
  DFT<T> *dft_outer = NULL;
  DFT<T> *dft_inner = NULL;
  std::vector<T>* prime_factors = NULL;
  void mul_twiddle_factors(bool inv);
 public:
  FFTCT(GF<T> *gf, T n, int id = 0, std::vector<T>* factors = NULL, T _w = 0);
  ~FFTCT();
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
 * @param n
 * @param id: index in the list of factors of n
 *
 * @return
 */
template <typename T>
FFTCT<T>::FFTCT(GF<T> *gf, T n, int id, std::vector<T>* factors, T _w) :
  DFT<T>(gf, n)
{
  if (factors == NULL) {
    this->prime_factors = new std::vector<T>();
    first_layer_fft = true;
    _get_prime_factors<T>(n, this->prime_factors);
    // w is of order-n
    w = gf->get_nth_root(n);
  } else {
    this->prime_factors = factors;
    first_layer_fft = false;
    w = _w;
  }
  inv_w = gf->inv(w);
  n1 = prime_factors->at(id);
  n2 = n/n1;

  w1 = gf->exp(w, n2);  // order of w1 = n1
  if (n1 == 2)
    this->dft_outer = new FFT2<T>(gf);
  else
    this->dft_outer = new DFTN<T>(gf, n1, w1);

  if (n2 > 1) {
    loop = true;
    w2 = gf->exp(w, n1);  // order of w2 = n2
    T _n2 = n/n1;
    // TODO: fix that FFT2K does not work here
    // if (_is_power_of_2<T>(_n2))
    //   this->dft_inner = new FFT2K<T>(gf, _n2);
    // else
    this->dft_inner = new FFTCT<T>(gf, _n2, id+1, this->prime_factors, w2);
    this->G = new Vec<T>(this->gf, this->n);
    this->Y = new VcVec<T>(this->G);
    this->X = new VcVec<T>(this->G);
  } else
    loop = false;
}

template <typename T>
FFTCT<T>::~FFTCT()
{
  if (prime_factors && first_layer_fft) delete prime_factors;
  if (dft_outer) delete dft_outer;
  if (dft_inner) delete dft_inner;
  if (X) delete X;
  if (Y) delete Y;
  if (G) delete G;
}

template <typename T>
void FFTCT<T>::mul_twiddle_factors(bool inv)
{
  T factor;
  T _w;
  if (inv)
    _w = inv_w;
  else
    _w = w;
  T base = 1;
  for (T i1 = 1; i1 < n1; i1++) {
    base = this->gf->mul(base, _w);  // base = _w^i1
    factor = base;  // init factor = base^1
    for (T k2 = 1; k2 < n2; k2++) {
      T loc = i1+n1*k2;;
      G->set(loc, this->gf->mul(G->get(loc), factor));
      // next factor = base^(k2+1)
      factor = this->gf->mul(factor, base);
    }
  }
}

template <typename T>
void FFTCT<T>::_fft(Vec<T> *output, Vec<T> *input, bool inv)
{
  X->set_vec(input);
  X->set_len(n2);
  Y->set_len(n2);
  for (T i1 = 0; i1 < n1; i1++) {
    Y->set_map(i1, n1);
    X->set_map(i1, n1);
    if (inv)
      this->dft_inner->fft_inv(Y, X);
    else
      this->dft_inner->fft(Y, X);
    // std::cout << "X:"; X->dump();
    // std::cout << "Y:"; Y->dump();
  }

  // multiply to twiddle factors
  mul_twiddle_factors(inv);

  X->set_vec(output);
  X->set_len(n1);
  Y->set_len(n1);
  for (T k2 = 0; k2 < n2; k2++) {
    Y->set_map(k2*n1, 1);
    X->set_map(k2, n2);
    if (inv)
      this->dft_outer->fft_inv(X, Y);
    else
      this->dft_outer->fft(X, Y);
    // std::cout << "Y:"; Y->dump();
    // std::cout << "X:"; X->dump();
  }
}

template <typename T>
void FFTCT<T>::fft(Vec<T> *output, Vec<T> *input)
{
  if (!loop)
    return dft_outer->fft(output, input);
  else
    return _fft(output, input, false);
}

template <typename T>
void FFTCT<T>::fft_inv(Vec<T> *output, Vec<T> *input)
{
  if (!loop)
    dft_outer->fft_inv(output, input);
  else
    _fft(output, input, true);
}

template <typename T>
void FFTCT<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  fft_inv(output, input);

  /*
   * We need to divide output to `N` for the inverse formular
   */

  if (this->first_layer_fft && (this->inv_n_mod_p > 1)) {
      output->mul_scalar(this->inv_n_mod_p);
  }
}
