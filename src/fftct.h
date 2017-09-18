
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
class FFTCT : public FFT<T>
{
 private:
  bool loop;
  T N;
  T n1;
  T n2;
  T w, w1, w2, inv_w;
  std::vector<T>* prime_factors;
  FFTN<T> *dft;
  FFTCT<T> *fftct = NULL;
 public:
  FFTCT(GF<T> *gf, T n, int id = 0, std::vector<T>* prime_factors = NULL,
    T _w = 0);
  ~FFTCT();
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
FFTCT<T>::FFTCT(GF<T> *gf, T n, int id, std::vector<T>* prime_factors, T _w) :
  FFT<T>(gf, n)
{
  if (id == 0) {
    this->N = n;
    prime_factors = gf->_get_prime_factors(n);
    // w is of order-n
    w = gf->get_nth_root(n);
  } else {
    this->prime_factors = prime_factors;
    w = _w;
  }
  inv_w = gf->inv(w);
  n1 = prime_factors->at(id);
  n2 = n/n1;

  w1 = gf->exp(w, n2);  // order of w1 = n1
  this->dft = new FFTN<T>(gf, n1, w1);

  if (n2 > 1) {
    loop = true;
    w2 = gf->exp(w, n1);  // order of w2 = n2
    this->fftct = new FFTCT<T>(gf, n/n1, id+1, prime_factors, w2);
  } else
    loop = false;
}

template <typename T>
FFTCT<T>::~FFTCT()
{
  delete dft;
  if (fftct) delete fftct;
}

template <typename T>
void FFTCT<T>::_fft(Vec<T> *output, Vec<T> *input, bool inv)
{
  Vec<T> G(this->gf, this->n);
  G.zero_fill();

  for (T i1 = 0; i1 < n1; i1++) {
    VmVec<T> Gcol(&G, n2, i1, n1);
    VmVec<T> x(input, n2, i1, n1);
    if (inv)
      this->fftct->ifft(&Gcol, &x);
    else
      this->fftct->fft(&Gcol, &x);
    // std::cout << "x:"; x.dump();
    // std::cout << "Gcol:"; Gcol.dump();
    // std::cout << "G:"; G.dump();
  }

  // multiply to twiddle factors
  T factor;
  T loc;
  for (T i1 = 1; i1 < n1; i1++) {
    for (T k2 = 1; k2 < n2; k2++) {
        loc = i1+n1*k2;
      if (inv)
        factor = this->gf->exp(inv_w, i1*k2);
      else
        factor = this->gf->exp(w, i1*k2);

      G.set(loc, this->gf->mul(G.get(loc), factor));
    }
  }

  for (T k2 = 0; k2 < n2; k2++) {
    VmVec<T> Grow(&G, n1, k2*n1, 1);
    VmVec<T> X(output, n1, k2, n2);
    if (inv)
      this->dft->fft_inv(&X, &Grow);
    else
      this->dft->fft(&X, &Grow);
    // std::cout << "Grow:"; Grow.dump();
    // std::cout << "X:"; X.dump();
    // std::cout << "output:"; output->dump();
  }
}

template <typename T>
void FFTCT<T>::fft(Vec<T> *output, Vec<T> *input)
{
  if (!loop)
    return dft->fft(output, input);
  else
    return _fft(output, input, false);
}

template <typename T>
void FFTCT<T>::ifft(Vec<T> *output, Vec<T> *input)
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
