
#ifndef __FFT_H__
#define __FFT_H__ 1

template<typename T>
class FFT
{
 private:
  GF<T> *gf;
  T w;
  T invw;
  int q;
  int N;
  Vec<T> *W;
  Vec<T> *invW;
  void compute_W(Vec<T> *_W, T _w);
  int _get_p(int i, int j);
  int _get_p0(int i, int j, int s);
  int _get_p1(int i, int j, int s);
  void _fft(Vec<T> *output, Vec<T> *input, Vec<T> *_W);
 public:
  FFT(GF<T> *gf, T w, int q);
  ~FFT();
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
};

#endif
