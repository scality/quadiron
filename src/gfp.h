
#ifndef __GFP_H__
#define __GFP_H__

#include "gf.h"

template<typename T>
class GFP : public GF<T>
{
public:
  GFP(T p);
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
};

#endif
