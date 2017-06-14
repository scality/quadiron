/* -*- mode: c++ -*- */

#pragma once

template<typename T>
class FFT
{
public:
  GF<T> *gf;
  
private:
  int n;
  int N;
  T inv_N;
  T w;
  T inv_w;
  Vec<T> *W;
  Vec<T> *inv_W;
  void compute_W(Vec<T> *_W, T _w);
  int _get_p(int i, int j);
  int _get_p0(int i, int j, int s);
  int _get_p1(int i, int j, int s);
  void _fft(Vec<T> *output, Vec<T> *input, Vec<T> *_W);
 public:
  FFT(GF<T> *gf, int n, T w);
  ~FFT();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
};

template <typename T>
FFT<T>::FFT(GF<T> *gf, int n, T w)
{
  this->gf = gf;
  this->n = n;
  this->N = __gf64._exp(2, n);
  this->inv_N = gf->inv(N);
  this->w = w;
  this->inv_w = gf->inv(w);
  this->W = new Vec<T>(gf, this->N + 1);
  this->inv_W = new Vec<T>(gf, this->N + 1);

  compute_W(W, w);
  compute_W(inv_W, inv_w);
}

template <typename T>
FFT<T>::~FFT()
{
  delete this->inv_W;
  delete this->W;
}

template <typename T>
void FFT<T>::compute_W(Vec<T> *_W, T _w)
{
  std::ostringstream filename;  

  filename << "W" << _w << ".cache";

  if (-1 == access(filename.str().c_str(), F_OK)) {
    std::ofstream file;
    file.open(filename.str().c_str(), std::ios::out);
    for (int i = 0;i <= N;i++) {
      _W->set(i, gf->exp(_w, i));
      file << _W->get(i) << "\n";
    }
  } else {
    std::ifstream file;
    int i = 0;
    file.open(filename.str().c_str(), std::ios::in);
    T tmp;
    while (file >> tmp) {
      _W->set(i, tmp);
      i++;
    }
    assert(i == N + 1);
  }
}

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
int FFT<T>::_get_p(int i, int j)
{
  assert(i >= 0 && i <= N-1);
  assert(j >= 1 && j <= n);
  
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
int FFT<T>::_get_p0(int i, int j, int s)
{
  assert(i >=0 && i <= N-1);
  assert(j >=1 && j <= n);
  
  return (j == s) ? 0 : _get_p(i, j);
}

/** 
 * @return _get_p(i, j) except when j==s it returns 1
 */
template <typename T>
int FFT<T>::_get_p1(int i, int j, int s)
{
  assert(i >=0 && i <= N-1);
  assert(j >=1 && j <= n);
  
  return (j == s) ? 1 : _get_p(i, j);
}

template <typename T>
void FFT<T>::_fft(Vec<T> *output, Vec<T> *input, Vec<T> *_W)
{
  Mat<T> *phi = new Mat<T>(gf, n+1, N);

  //compute phi[0][i]
  for (int i = 0;i <= N-1;i++)
    phi->set(0, i, input->get(i));
  
  //compute phi[1][i]
  for (int i = 0;i <= N-1;i++) {

    T tmp1 = 0;
    for (int k = 1;k <= n;k++)
      tmp1 += _get_p0(i, k, n) * __gf64._exp(2, k-1);

    T tmp2 = 0;
    for (int k = 1;k <= n;k++)
      tmp2 += _get_p1(i, k, n) * __gf64._exp(2, k-1);

    DoubleT<T> val = 
      DoubleT<T>(_W->get(__gf64._exp(2, n-1) * _get_p(i, n))) *
      phi->get(0, tmp2) +
      phi->get(0, tmp1);

    phi->set(1, i, val % gf->card());
  }

  //compute phi[2,3] to phi[n,N-1]
  for (int mm = 2; mm <= n;mm++) {

    for (int i = 0;i <= N-1;i++) {
      
      T tmp1 = 0;
      for (int k = 1;k <= mm;k++)
        tmp1 += _get_p(i, n-mm+k) * __gf64._exp(2, mm-k);

      T tmp2 = 0;
      for (int k = 1;k <= n;k++)
        tmp2 += _get_p0(i, k, n-mm+1) * __gf64._exp(2, k-1);

      T tmp3 = 0;
      for (int k = 1;k <= n;k++)
        tmp3 += _get_p1(i, k, n-mm+1) * __gf64._exp(2, k-1);

      DoubleT<T>val = 
        DoubleT<T>(_W->get(__gf64._exp(2, n-mm) * tmp1)) *
        phi->get(mm-1, tmp3) +
        phi->get(mm-1, tmp2);

      phi->set(mm, i, val % gf->card());
    }
  }

  //compute FFT
  for (int i = 0;i <= N-1;i++) {
    
    T tmp = 0;
    for (int k = 1;k <= n;k++)
      tmp += _get_p(i, n-k+1) * __gf64._exp(2, k-1);
    
    output->set(i, phi->get(n, tmp));
  }

  delete phi;
}

template <typename T>
void FFT<T>::fft(Vec<T> *output, Vec<T> *input)
{
  _fft(output, input, W);
}

template <typename T>
void FFT<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  _fft(output, input, inv_W);
  output->mul_scalar(inv_N);
}
