
#include "ntl.h"

template <typename T>
FFT<T>::FFT(GF<T> *gf, T omega, int q)
{
  this->gf = gf;
  this->omega = omega;
  this->q = q;
  this->N = gf->__exp(2, q);
  this->W = new Vec<T>(gf, this->N);

  compute_omegas();
}

template <typename T>
FFT<T>::~FFT()
{
  delete this->W;
}

template <typename T>
void FFT<T>::compute_omegas(void)
{
  std::ostringstream filename;  

  filename << "W" << omega << ".cache";

  if (-1 == access(filename.str().c_str(), F_OK)) {
    std::ofstream file;
    file.open(filename.str().c_str(), std::ios::out);
    for (int i = 0;i < N;i++) {
      VEC_ITEM(W, i) = gf->exp(omega, i);
      file << VEC_ITEM(W, i) << "\n";
    }
  } else {
    std::ifstream file;
    int i = 0;
    file.open(filename.str().c_str(), std::ios::in);
    while (!file.eof()) {
      file >> VEC_ITEM(W, i);
      i++;
    }
    assert(i == N + 1);
  }
}

/** 
 * Simulate the bit matrix
 *
 * [0][0]=0   [0][1]=0   ... [0][N-1]=0    <= 0
 * [1][0]=1   [1][1]=0   ... [1][N-1]=0    <= 1
 * [2][0]=0   [2][1]=1   ... [2][N-1]=0    <= 2
 * ...
 * [N-1][0]=1 [N-1][1]=1 ... [N-1][N-1]=0  <= N-1
 *
 * where the numbers from 0 .. N-1 are encoded in reverse order
 * 
 * @param i 
 * @param j 
 * 
 * @return 
 */
template <typename T>
int FFT<T>::_get_p(int i, int j)
{
  assert(i >=0 && i < N - 1);
  assert(j >=0 && j < N - 1);
  
  int x = i;
  int y = 0;
  
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

template <typename T>
int FFT<T>::_get_p0(int i, int j, int x)
{
  assert(i >=0 && i < N);
  assert(j >=1 && j <= q);
  
  return (j == x) ? 0 : _get_p(i, j);
}

template <typename T>
int FFT<T>::_get_p1(int i, int j, int x)
{
  assert(i >=0 && i < N);
  assert(j >=1 && j <= q);
  
  return (j == x) ? 1 : _get_p(i, j);
}

template <typename T>
void FFT<T>::fft(Vec<T> *output, Vec<T> *input)
{
  Mat<T> *phi = new Mat<T>(gf, q+1, N-1);

  //compute phi[0][i]
  for (int i = 0;i < N-1;i++)
    MAT_ITEM(phi, 0, i) = VEC_ITEM(input, i);
  
  //compute phi[1][i]
  for (int i = 0;i < N-1;i++) {

    T tmp0 = 0;
    for (int k = 1;k < q;k++)
      tmp0 += _get_p0(i, k, q) * gf->__exp(2, k-1);

    T tmp1 = 0;
    for (int k = 1;k < q;k++)
      tmp1 += _get_p1(i, k, q) * gf->__exp(2, k-1);

    //XXX redefine with GF operations
    MAT_ITEM(phi, 1, i) = 
      (MAT_ITEM(phi, 0, tmp0) +
       VEC_ITEM(W, gf->__exp(2, q-1) * _get_p(i, q)) *
       MAT_ITEM(phi, 0, tmp1)) % gf->card();
  }

  Vec<T> *t = new Vec<T>(gf, q);

  //compute phi[2,3] to phi[q,i]
  for (int mm = 2; mm < q;mm++) {
    for (int i = 0;i < N-1;i++) {

      T tmp1 = 0;
      for (int k = 1;k < mm;k++)
        tmp1 += _get_p(i, q-mm+k) * gf->__exp(2, mm-k);

      VEC_ITEM(t, mm) = gf->__exp(2, q-mm) * tmp1;

      T tmp2 = 0;
      for (int k = 1;k < mm;k++)
        tmp2 += _get_p0(i, k, q-mm+1) * gf->__exp(2, k-1);

      T tmp3 = 0;
      for (int k = 1;k < q;k++)
        tmp3 += _get_p1(i, k, q-mm+1) * gf->__exp(2, k-1);

      MAT_ITEM(phi, mm, i) = 
        (MAT_ITEM(phi, mm-1, tmp2) +
         VEC_ITEM(W, 
                  VEC_ITEM(t, mm) *
                  MAT_ITEM(phi, mm-1, tmp3))) % gf->card();
    }

    //compute FFT
    for (int i = 0;i < N-1;i++) {

      T tmp = 0;
      for (int k = 1;k < q;k++)
        tmp += _get_p(i, q-k+1) * gf->__exp(2, k-1);

      VEC_ITEM(output, i) =
        MAT_ITEM(phi, q, tmp);
    }
  }


  delete t;
  delete phi;
}
