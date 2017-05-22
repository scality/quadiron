
#ifndef __GF2N_H__
#define __GF2N_H__

#include "gf.h"

template<typename T>
class GF2N : public GF<T>
{
 private:
  T my_card;
  T prim_poly;
  T *gflog;
  T *gfilog;
  void setup_tables(void);

 public:
  GF2N(T n);
  T card(void);
  bool check(T a);
  T neg(T a);
  T add(T a, T b);
  T sub(T a, T b);
  T mul(T a, T b);
  T div(T a, T b);
  T pow(T a);
  T log(T a);
  T inv(T a);
  T weak_rand(void);
};

#endif
