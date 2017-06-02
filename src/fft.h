
#ifndef __FFT_H__
#define __FFT_H__ 1

template<typename T>
class FFT
{
 private:
  GF<T> *gf;
  int n;
  int N;
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
  void _debug(void);
};

#endif
