
#ifndef __GFP_H__
#define __GFP_H__

#include "gf.h"

template<typename T>
class GFP : public GF<T>
{
public:
  GFP(T p);
  T card(void);
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
