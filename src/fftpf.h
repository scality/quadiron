
/* -*- mode: c++ -*- */

#pragma once

/**
 * Prime-factor FFT algorithm
 *  https://en.wikipedia.org/wiki/Prime-factor_FFT_algorithm
 *
 * Suppose n = n1*n2 where n1 and n2 are relatively prime
 *  X[k] = sum_i x_i * w_n^ik is plit into two smaller ones
 *
 * 1. Calculate a, b, c, d as
 *    a = n2, b = n1, c = n1+1, d = n2+1
 * 2. Index mapping
 *  k = c*k1 + d*k2
 *  i = a*i1 + b*i2
 *  where 1 <= i1, k1 <= n1-1, 1 <= i2, k2 <= n2-1
 * 3. DFT formular
 *    Y = DFT(X, N), i.e.
 *  X, Y are length-N vector and
 *    Y[k] = sum_i X_i * w_N^ik,
 *  where w is of order N, i.e. w^N=1
 * 4. Hence,
 *    X[c*k1 + d*k2] = DFT(G(i2), n2),
 *  where G(i2) is a length-n1 vector:
 *    G(i2) = DFT(y, n1),
 *  where y[i1] = x[a*i1 + b*i2]
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
  FFTN<T> *dft;
  FFTPF<T> *fft_ct = NULL;
 public:
  FFTPF(GF<T> *gf, int n, int id = 0, std::vector<T>* prime_factors = NULL,
    int _w = 0);
  ~FFTPF();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
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
FFTPF<T>::FFTPF(GF<T> *gf, int n, int id, std::vector<T>* prime_factors,
  int _w) : FFT<T>(gf, n)
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
  c = n2*(gf->inv(n2) % n1);
  d = n1*(gf->inv(n1)% n2);

  w1 = gf->exp(w, n2);  // order of w1 = n1
  assert(gf->get_nth_root(n1) == w1);
  this->dft = new FFTN<T>(gf, n1, w1);

  if (id < prime_factors->size() - 1) {
    loop = true;
    w2 = gf->exp(w, n1);  // order of w2 = n2
    this->fft_ct = new FFTPF<T>(gf, n/n1, id+1, prime_factors, w2);
  } else
    loop = false;
}

template <typename T>
FFTPF<T>::~FFTPF()
{
  delete dft;
  if (fft_ct) delete fft_ct;
}

template <typename T>
void FFTPF<T>::_fft(Vec<T> *output, Vec<T> *input, bool inv)
{
  // std::cout << "n1:" << n1 << " n2:" << n2 << std::endl;
  // std::cout << "a:" << a << " b:" << b << " c:" << c << " d:" << d << std::endl;
  Vec<T> G(this->gf, this->n);
  G.zero_fill();

  for (T i1 = 0; i1 < n1; i1++) {
    T offset = (a*i1) % this->n;
    VmVec<T> Gcol(&G, n2, offset, b);
    Vec<T> _Gcol(this->gf, n2);
    VmVec<T> x(input, n2, offset, b);
    if (inv)
      this->fft_ct->ifft(&_Gcol, x.toVec());
    else
      this->fft_ct->fft(&_Gcol, x.toVec());
    Gcol.import(&_Gcol);
    // std::cout << "x:"; x.dump();
    // std::cout << "_Gcol:"; _Gcol.dump();
    // std::cout << "Gcol:"; Gcol.dump();
    // std::cout << "G:"; G.dump();
  }

  for (T k2 = 0; k2 < n2; k2++) {
    T offset = (d*k2) % this->n;
    VmVec<T> Grow(&G, n1, offset, c);
    VmVec<T> X(output, n1, offset, c);
    Vec<T> _X(this->gf, n1);
    if (inv)
      this->dft->fft_inv(&_X, Grow.toVec());
    else
      this->dft->fft(&_X, Grow.toVec());
    X.import(&_X);
    // std::cout << "Grow:"; Grow.dump();
    // std::cout << "_X:"; _X.dump();
    // std::cout << "X:"; X.dump();
    // std::cout << "output:"; output->dump();
  }

  /*
   * We need to divide output to `N` for the inverse formular
   * We need to multiply output to `2` since they are computed by fftn for n=2
   */
  if (inv && loop)
    output->mul_scalar(this->gf->inv(this->N) % this->gf->p);
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
    return dft->ifft(output, input);
  else
    return _fft(output, input, true);
}
