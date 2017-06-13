/* -*- mode: c++ -*- */
#pragma once
#include "gf.h"

template<typename T>
class GFP : public GF<T>
{
public:
  GFP(T p);
  T card(void);
  bool check(T a);
  T max(T a, T b);
  T min(T a, T b);
  T neg(T a);
  T add(T a, T b);
  T sub(T a, T b);
  T mul(T a, T b);
  T div(T a, T b);
  T inv(T a);
  T invExp(T a);
  T invBezout(T a);
  T exp(T a, T b);
  T log(T a, T b);
};

template <typename T>
GFP<T>::GFP(T p) : GF<T>(p, 1)
{
  //XXX shall check that p is prime
  assert(p != 2);
}

template <typename T>
T GFP<T>::card(void)
{
  return this->p;
}

template <typename T>
bool GFP<T>::check(T a)
{
  return (a >= 0 && a < this->p);
}

template <typename T>
T GFP<T>::max(T a, T b)
{
  assert(check(a));
  assert(check(b));

  return (a >= b) ? a : b;
}

template <typename T>
T GFP<T>::min(T a, T b)
{
  assert(check(a));
  assert(check(b));

  return (a < b) ? a : b;
}

template <typename T>
T GFP<T>::neg(T a)
{
  assert(check(a));

  return sub(0, a);
}

template <typename T>
T GFP<T>::add(T a, T b)
{
  assert(check(a));
  assert(check(b));

  return (a + b) % this->p;
}

template <typename T>
T GFP<T>::sub(T a, T b)
{
  assert(check(a));
  assert(check(b));

  if (a >= b)
    return a - b;
  else
    return this->p - (b - a);
}

template <typename T>
T GFP<T>::mul(T a, T b)
{
  assert(check(a));
  assert(check(b));

  return T((DoubleT<T>(a) * b) % this->p);
}

template <typename T>
T GFP<T>::div(T a, T b)
{
  assert(check(a));
  assert(check(b));

  T inv_b = inv(b);

  return mul(a, inv_b);
}

/** 
 * Inverse by exponentiation
 * 
 * @param a 
 * 
 * @return 
 */
template <typename T>
T GFP<T>::invExp(T a)
{
  assert(check(a));

  return GF<T>::_mod_exp(a, this->p - 2, this->p);
}

/** 
 * Inverse by Bezout
 * 
 * @param a 
 * 
 * @return 
 */
template <typename T>
T GFP<T>::invBezout(T a)
{
  assert(check(a));

  SignedDoubleT<T> x = a;
  SignedDoubleT<T> n = this->p;
  SignedDoubleT<T> bezout[2];

  GF<T>::_extended_gcd(x, n, bezout, NULL);
  if (bezout[0] < 0)
    bezout[0] = this->p + bezout[0];
  return bezout[0];
}

template <typename T>
T GFP<T>::inv(T a)
{
  return invBezout(a);
}

template <typename T>
T GFP<T>::exp(T a, T b)
{
  assert(check(a));
  assert(check(b));

  return GF<T>::_mod_exp(a, b, this->p);
}

template <typename T>
T GFP<T>::log(T a, T b)
{
  assert(check(a));

  return GF<T>::logNaive(a, b);
}
