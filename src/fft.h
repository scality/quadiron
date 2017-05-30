
#ifndef __FFT_H__
#define __FFT_H__ 1

template<typename T>
class FFT
{
 private:
  GF<T> *gf;
  T omega;
  int q;
  int N;
  Vec<T> *W;
  void compute_omegas(void);
 public:
  FFT(GF<T> *gf, T omega, int q);
  ~FFT();
  int _get_p(int i, int j);
  int _get_p0(int i, int j, int s);
  int _get_p1(int i, int j, int s);
  void fft(Vec<T> *output, Vec<T> *input);
};

#endif
