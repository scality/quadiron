
/* -*- mode: c++ -*- */

#pragma once

/**
 * Prime-factor FFT algorithm
 *  https://en.wikipedia.org/wiki/Prime-factor_FFT_algorithm
 *
 * Suppose n = n1*n2 where n1 and n2 are relatively prime, i.e. gcd(n1, n2) = 1
 *  X[k] = sum_i x_i * w_n^ik is plit into two smaller DFTs
 *  by using the index mapping
 *    i = a*i1 + b*i2
 *    k = c*k1 + d*k2
 *    with a = n2, b = n1, c = n2*inv_mod(n2, n1), d = n1*inv_mod(n1, n2)
 *      where 1 <= i1, k1 <= n1-1, 1 <= i2, k2 <= n2-1
 *        and inv_mod(u, v) := multiplicative inverse of u reduced modulo v,
 *          i.e. (u * inv_mod(u, v)) % v == 1
 *
 *  X[c*k1 + d*k2] =
 *    sum_i1 [
 *        sum_i2 x[a*i1 + b*i2] * w_2^(i2*k2)
 *    ] w_1^(i1*k1)
 *
 *  where
 *    w is nth root, w_1 = w^n2, w_2 = w^n1
 *
 *  Note, the index mapping leads to w^(i*k) = w_1^(i1*k1) * w_2^(i2*k2)
 *
 *  Step1: calculate DFT of the inner parenthese, i.e. sum_i2...
 *  Step2: calculate DFT of the outer parenthese, i.e. sum_i1
 *
 * @param output
 * @param input
 */
template<typename T>
class FFTPF : public FFT<T>
{
 private:
  bool loop;
  T N;
  T n1;
  T n2;
  T w, w1, w2;
  T a, b, c, d;
  std::vector<T>* prime_factors;
  Vec<T> *G = NULL;
  VmVec<T> *Y = NULL;
  VmVec<T> *X = NULL;
  FFTN<T> *dft;
  FFTPF<T> *fftpf = NULL;
 public:
  FFTPF(GF<T> *gf, T n, int id = 0, std::vector<T>* prime_factors = NULL,
    T _w = 0);
  ~FFTPF();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
 private:
  void _fft(Vec<T> *output, Vec<T> *input, bool inv);
  T _inverse_mod(T nb, T mod);
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
FFTPF<T>::FFTPF(GF<T> *gf, T n, int id, std::vector<T>* prime_factors, T _w) :
  FFT<T>(gf, n)
{
  if (id == 0) {
    this->N = n;
    prime_factors = gf->_get_coprime_factors(n);
    // w is of order-n
    w = gf->get_nth_root(n);
  } else {
    this->prime_factors = prime_factors;
    w = _w;
  }
  n1 = prime_factors->at(id);
  n2 = n/n1;

  a = n2;
  b = n1;
  c = n2*(this->_inverse_mod(n2, n1));
  d = n1*(this->_inverse_mod(n1, n2));

  w1 = gf->exp(w, n2);  // order of w1 = n1
  this->dft = new FFTN<T>(gf, n1, w1);

  if (n2 > 1) {
    loop = true;
    w2 = gf->exp(w, n1);  // order of w2 = n2
    this->fftpf = new FFTPF<T>(gf, n/n1, id+1, prime_factors, w2);
    this->G = new Vec<T>(this->gf, this->n);
    this->Y = new VmVec<T>(this->G);
    this->X = new VmVec<T>(this->G);
  } else
    loop = false;
}

template <typename T>
FFTPF<T>::~FFTPF()
{
  delete dft;
  if (fftpf) delete fftpf;
  if (X) delete X;
  if (Y) delete Y;
  if (G) delete G;
}

/*
 * It calculates multiplicative inverse of a number modulo mod
 *  NOTE: it should be optimized
 */
template <typename T>
T FFTPF<T>::_inverse_mod(T nb, T mod)
{
  for (T n = 1; n < mod; n++) {
    if ((n * nb) % mod == 1) return n;
  }
  return 0;
}

template <typename T>
void FFTPF<T>::_fft(Vec<T> *output, Vec<T> *input, bool inv)
{
  X->set_vec(input);
  X->set_len(n2);
  Y->set_len(n2);
  for (T i1 = 0; i1 < n1; i1++) {
    Y->set_map(i1, n1);
    X->set_map((a*i1) % this->n, b);
    if (inv)
      this->fftpf->ifft(Y, X);
    else
      this->fftpf->fft(Y, X);
      // std::cout << "X:"; X->dump();
      // std::cout << "Y:"; Y->dump();
  }

  X->set_vec(output);
  X->set_len(n1);
  Y->set_len(n1);
  for (T k2 = 0; k2 < n2; k2++) {
    Y->set_map(k2*n1, 1);
    X->set_map((d*k2) % this->n, c);
    if (inv)
      this->dft->fft_inv(X, Y);
    else
      this->dft->fft(X, Y);
      // std::cout << "Y:"; Y->dump();
      // std::cout << "X:"; X->dump();
  }
}

template <typename T>
void FFTPF<T>::fft(Vec<T> *output, Vec<T> *input)
{
  if (!loop)
    return dft->fft(output, input);
  else
    return _fft(output, input, false);
}

template <typename T>
void FFTPF<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  if (!loop)
    dft->fft_inv(output, input);
  else
    _fft(output, input, true);

  /*
   * We need to divide output to `N` for the inverse formular
   */
  if (this->N == this->n) {
    T coef = this->gf->inv(this->N) % this->gf->p;
    if (coef > 0)
      output->mul_scalar(coef);
  }
}
