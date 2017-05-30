
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
      W->set(i, gf->exp(omega, i));
      std::cerr << i << " " << N << "\n";
      file << W->get(i) << "\n";
    }
  } else {
    std::ifstream file;
    int i = 0;
    file.open(filename.str().c_str(), std::ios::in);
    T tmp;
    while (file >> tmp) {
      W->set(i, tmp);
      i++;
    }
    assert(i == N);
  }
}

/** 
 * Simulate the bit matrix
 *
 * [0][0]=undef   [0][1]=0   [0][2]=0   ... [0][q]=0    <= 0
 * [1][0]=undef   [1][1]=1   [1][2]=0   ... [1][q]=0    <= 1
 * [2][0]=undef   [2][1]=0   [2][2]=1   ... [2][q]=0    <= 2
 * ...
 * [N-1][0]=undef [N-1][1]=1 [N-1][2]=1 ... [N-1][q]=1  <= N-1
 *
 * where the numbers from 0 .. N-1 are encoded in reverse order
 * 
 * @param i 0 <= i <= N-1
 * @param j 1 <= j <= q
 * 
 * @return
 */
template <typename T>
int FFT<T>::_get_p(int i, int j)
{
  assert(i >= 0 && i <= N-1);
  assert(j >= 1 && j <= q);
  
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
  assert(j >=1 && j <= q);
  
  return (j == s) ? 0 : _get_p(i, j);
}

/** 
 * @return _get_p(i, j) except when j==s it returns 1
 */
template <typename T>
int FFT<T>::_get_p1(int i, int j, int s)
{
  assert(i >=0 && i <= N-1);
  assert(j >=1 && j <= q);
  
  return (j == s) ? 1 : _get_p(i, j);
}

template <typename T>
void FFT<T>::fft(Vec<T> *output, Vec<T> *input)
{
  Mat<T> *phi = new Mat<T>(gf, q+1, N);

  //compute phi[0][i]
  for (int i = 0;i <= N-1;i++)
    phi->set(0, i, input->get(i));
  
  //compute phi[1][i]
  for (int i = 0;i <= N-1;i++) {

    T tmp0 = 0;
    for (int k = 1;k <= q;k++)
      tmp0 += _get_p0(i, k, q) * gf->__exp(2, k-1);

    T tmp1 = 0;
    for (int k = 1;k <= q;k++)
      tmp1 += _get_p1(i, k, q) * gf->__exp(2, k-1);

    phi->set(1, i,
            (phi->get(0, tmp0) +
             W->get(gf->__exp(2, q-1) * _get_p(i, q)) *
             phi->get(0, tmp1)) % gf->card());
  }

  //compute phi[2,3] to phi[q,N-1]
  for (int mm = 2; mm <= q;mm++) {

    for (int i = 0;i <= N-1;i++) {
      
      T tmp1 = 0;
      for (int k = 1;k <= mm;k++)
        tmp1 += _get_p(i, q-mm+k) * gf->__exp(2, mm-k);

      T tmp2 = 0;
      for (int k = 1;k <= q;k++)
        tmp2 += _get_p0(i, k, q-mm+1) * gf->__exp(2, k-1);

      T tmp3 = 0;
      for (int k = 1;k <= q;k++)
        tmp3 += _get_p1(i, k, q-mm+1) * gf->__exp(2, k-1);

      phi->set(mm, i,
              (phi->get(mm-1, tmp2) +
               W->get(gf->__exp(2, q-mm) * tmp1) *
               phi->get(mm-1, tmp3)) % gf->card());
    }

    //compute FFT
    for (int i = 0;i <= N-1;i++) {

      T tmp = 0;
      for (int k = 1;k <= q;k++)
        tmp += _get_p(i, q-k+1) * gf->__exp(2, k-1);
      
      output->set(i, phi->get(q, tmp));
    }
  }

  delete phi;
}
