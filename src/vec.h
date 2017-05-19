
#ifndef __VEC_H__
#define __VEC_H__ 1

template<typename T>
class Vec
{
 public:
  GF<T> *gf;
  int n;
  T n_cols;
  T *mem;
#define VEC_ITEM(vec, i) ((vec)->mem[(i)])
  Vec(GF<T> *gf, int n);
  void zero(void);
};

#endif
