/* -*- mode: c++ -*- */

#pragma once

/**
 * Algorithm for very large n=2^l
 */
template<typename T>
class FFTLN : public DFT<T>
{
 private:
  int l;
  T w;
  T inv_w;
  Vec<T> *W;
  Vec<T> *inv_W;
  Mat<T> *tmp2;
  Mat<T> *tmp3;
  Mat<T> *tmp4;
  Vec<T> *tmp5;
  int _get_p(int i, int j);
  int _get_p0(int i, int j, int s);
  int _get_p1(int i, int j, int s);
  void _pre_compute_consts();
  void _fft(Vec<T> *output, Vec<T> *input, Vec<T> *_W);
 public:
  FFTLN(GF<T> *gf, int l, T w);
  ~FFTLN();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
  void fft_inv(Vec<T> *output, Vec<T> *input);
};

template <typename T>
FFTLN<T>::FFTLN(GF<T> *gf, int l, T w) : DFT<T>(gf, _exp2<T>(l))
{
  this->l = l;
  this->w = w;
  this->inv_w = gf->inv(w);

  this->W = new Vec<T>(gf, this->n + 1);
  this->inv_W = new Vec<T>(gf, this->n + 1);
  gf->compute_omegas_cached(W, this->n, w);
  gf->compute_omegas_cached(inv_W, this->n, inv_w);

  _pre_compute_consts();
}

template <typename T>
FFTLN<T>::~FFTLN()
{
  delete this->inv_W;
  delete this->W;
  delete this->tmp2;
  delete this->tmp3;
  delete this->tmp4;
  delete this->tmp5;
}

/**
 * Representing powers of root as a vector is more practical than a
 * matrix (e.g. (n+k=2^15) => would be 1 billion entries !)
 *
 * @param _W vector of powers of roots
 * @param _w Nth root of unity
 */
/**
 * Simulate the bit matrix
 *
 * [0][0]=undef   [0][1]=0   [0][2]=0   ... [0][n]=0    <= 0
 * [1][0]=undef   [1][1]=1   [1][2]=0   ... [1][n]=0    <= 1
 * [2][0]=undef   [2][1]=0   [2][2]=1   ... [2][n]=0    <= 2
 * ...
 * [N-1][0]=undef [N-1][1]=1 [N-1][2]=1 ... [N-1][n]=1  <= N-1
 *
 * where the numbers from 0 .. N-1 are encoded in reverse order
 *
 * @param i 0 <= i <= N-1
 * @param j 1 <= j <= n
 *
 * @return
 */
template <typename T>
int FFTLN<T>::_get_p(int i, int j)
{
  assert(i >= 0 && i <= this->n-1);
  assert(j >= 1 && j <= this->l);

  int x = i;
  int y = 1;

  do {
    if (y == j) {
      if (x & 1)
        return 1;
      else
        return 0;
    }
    y++;
  } while (x >>= 1);

  return 0;
}

/**
 * @return _get_p(i, j) except when j==s it returns 0
 */
template <typename T>
int FFTLN<T>::_get_p0(int i, int j, int s)
{
  assert(i >=0 && i <= this->n-1);
  assert(j >=1 && j <= this->l);

  return (j == s) ? 0 : _get_p(i, j);
}

/**
 * @return _get_p(i, j) except when j==s it returns 1
 */
template <typename T>
int FFTLN<T>::_get_p1(int i, int j, int s)
{
  assert(i >=0 && i <= this->n-1);
  assert(j >=1 && j <= this->l);

  return (j == s) ? 1 : _get_p(i, j);
}

template <typename T>
void FFTLN<T>::_pre_compute_consts()
{
  tmp2 = new Mat<T>(this->gf, this->l+1, this->n);
  tmp3 = new Mat<T>(this->gf, this->l+1, this->n);
  tmp4 = new Mat<T>(this->gf, this->l+1, this->n);
  tmp5 = new Vec<T>(this->gf, this->n);

  tmp2->zero_fill();
  tmp3->zero_fill();
  tmp4->zero_fill();
  tmp5->zero_fill();

  for (int i = 1; i <= this->l; i++) {
    for (int j = 0; j <= this->n-1; j++) {
      T _tmp1 = 0;
      for (int k = 1; k <= i; k++)
        _tmp1 += _get_p(j, this->l-i+k) * _exp2<T>(i-k);

      T _tmp2 = 0;
      for (int k = 1; k <= this->l; k++)
        _tmp2 += _get_p0(j, k, this->l-i+1) * _exp2<T>(k-1);
      tmp2->set(i, j, _tmp2);

      T _tmp3 = 0;
      for (int k = 1; k <= this->l; k++)
        _tmp3 += _get_p1(j, k, this->l-i+1) * _exp2<T>(k-1);
      tmp3->set(i, j, _tmp3);

      T _tmp4 = _tmp1 * _exp2<T>(this->l-i);
      tmp4->set(i, j, _tmp4);
    }
  }

  for (int i = 0; i <= this->n-1; i++) {
    T _tmp5 = 0;
    for (int k = 1; k <= this->l; k++)
      _tmp5 += _get_p(i, this->l-k+1) * _exp2<T>(k-1);
    tmp5->set(i, _tmp5);
  }
}

template <typename T>
void FFTLN<T>::_fft(Vec<T> *output, Vec<T> *input, Vec<T> *_W)
{
  Mat<T> phi(this->gf, this->l+1, this->n);

  // compute phi[0][i]
  for (int i = 0; i <= this->n-1; i++)
    phi.set(0, i, input->get(i));

  for (int i = 1; i <= this->l; i++) {
    for (int j = 0; j <= this->n-1; j++) {
      DoubleT<T>val =
        DoubleT<T>(_W->get(tmp4->get(i, j))) *
        phi.get(i-1, tmp3->get(i, j)) +
        phi.get(i-1, tmp2->get(i, j));

      phi.set(i, j, val % this->gf->card());
    }
  }

  // compute FFT
  for (int i = 0; i <= this->n-1; i++)
    output->set(i, phi.get(this->l, tmp5->get(i)));
}

template <typename T>
void FFTLN<T>::fft(Vec<T> *output, Vec<T> *input)
{
  _fft(output, input, W);
}

template <typename T>
void FFTLN<T>::fft_inv(Vec<T> *output, Vec<T> *input)
{
  _fft(output, input, inv_W);
}

template <typename T>
void FFTLN<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  fft_inv(output, input);
  if (this->inv_n_mod_p > 1)
    output->mul_scalar(this->inv_n_mod_p);
}
