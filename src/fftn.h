/* -*- mode: c++ -*- */

#pragma once

/** 
 * Algorithm for very small n
 */
template<typename T>
class FFTN : public FFT<T>
{
private:
  Mat<T> *W;
  Mat<T> *inv_W;
  void compute_W(Mat<T> *_W, T _w);
  void _fft(Vec<T> *output, Vec<T> *input, Mat<T> *_W);
 public:
  FFTN(GF<T> *gf, int n, T w);
  ~FFTN();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
};

template <typename T>
FFTN<T>::FFTN(GF<T> *gf, int n, T w) : FFT<T>(gf, n, w)
{
  this->W = new Mat<T>(gf, this->n, this->n);
  this->inv_W = new Mat<T>(gf, this->n, this->n);

  compute_W(W, w);
  compute_W(inv_W, this->inv_w);
}

template <typename T>
FFTN<T>::~FFTN()
{
  delete this->inv_W;
  delete this->W;
}

/** 
 * Compute matrix 
 * 
 * @param _W matrix of powers of roots
 * @param _w nth root of unity
 */
template <typename T>
void FFTN<T>::compute_W(Mat<T> *_W, T _w)
{
  for (int i = 0;i <= this->n-1;i++) {
    for (int j = 0;j <= this->n-1;j++) {
      _W->set(i, j, this->gf->exp(_w, i*j));
    }
  }
  //_W->dump();
}

template <typename T>
void FFTN<T>::_fft(Vec<T> *output, Vec<T> *input, Mat<T> *_W)
{
  _W->mul(output, input);
}

template <typename T>
void FFTN<T>::fft(Vec<T> *output, Vec<T> *input)
{
  _fft(output, input, W);
}

template <typename T>
void FFTN<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  _fft(output, input, inv_W);
  output->mul_scalar(this->inv_n);
}
