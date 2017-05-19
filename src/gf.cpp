
#include "ntl.h"

template <typename T>
GF<T>::GF(T p, T n)
{
  this->p = p;
  this->n = n;
}

template <typename T>
static T ipow(T a, T b)
{
  T r = 1;

  while (b)
    {
      if (b & 1)
        r *= a;
      b >>= 1;
      a *= a;
    }

  return r;
}

template <typename T>
T GF<T>::generic_card(GF *gf)
{
  return ipow(p, n);
}

template <typename T>
T GF<T>::generic_pow(GF *gf, T a, T b)
{
  T r;
  int i;

  if (0 == b)
    return 1;

  if (1 == b)
    return a;

  r = a;
  for (i = 1; i < b;i++) {
    r = gf->mul(r, a);
  }

  return r;
}

template <typename T>
T GF<T>::generic_trial_mult_log(GF *gf, T a, T b)
{
  T r;

  for (r = 1;r < card();r++) {
    if (gf->pow(a, r) == b)
      return r;
  }

  //not found
  throw NTL_EX_NOT_FOUND;
}
