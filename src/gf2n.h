
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
  T zero(void);
  T one(void);
  T check(T a);
  bool eq(T a, T b);
  T add(T a, T b);
  T sub(T a, T b);
  T mul(T a, T b);
  T div(T a, T b);
  T pow(T a, T b);
  T log(T a, T b);
  T inv(T a);
  T weak_rand(void);
};

#endif
