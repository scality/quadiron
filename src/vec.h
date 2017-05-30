
#ifndef __VEC_H__
#define __VEC_H__ 1

template<typename T>
class Vec
{
 public:
  GF<T> *gf;
  int n;
  T *mem;
  Vec(GF<T> *gf, int n);
  ~Vec();
  void zero_fill(void);
  void set(int i, T val);
  T get(int i);
  void dump(void);
};

#endif
