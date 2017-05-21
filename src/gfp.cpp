
#include "ntl.h"

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
T GFP<T>::zero(void)
{
  return 0;
}

template <typename T>
T GFP<T>::one(void)
{
  return 1;
}

template <typename T>
T GFP<T>::check(T a)
{
  assert(a >= 0 && a < this->p);
  return a;
}

template <typename T>
bool GFP<T>::eq(T a, T b)
{
  return a == b;
}

template <typename T>
T GFP<T>::add(T a, T b)
{
  return sub(a, this->p - b);
}

template <typename T>
T GFP<T>::sub(T a, T b)
{
  if (a >= b)
    return a - b;
  else
    return this->p - (b - a);
}

template <typename T>
T GFP<T>::mul(T a, T b)
{
  using DoubleT = typename Double<T>::T;
  return T((DoubleT(a) * b) % this->p);
}

template <typename T>
T GFP<T>::div(T a, T b)
{
  T inv_b = inv(b);
  return mul(a, inv_b);
}

template <typename T>
T GFP<T>::inv(T a)
{
  return pow(a, this->p - 2);
}

template <typename T>
T GFP<T>::pow(T a, T b)
{
  return GF<T>::generic_pow(this, a, b);
}

template <typename T>
T GFP<T>::log(T a, T b)
{
  return GF<T>::generic_trial_mult_log(this, a, b);
}

template <typename T>
T GFP<T>::weak_rand(void)
{
  T r;

 retry:
  r = rand() % this->p;
  if (zero() == r)
    goto retry;
  return r;
}
