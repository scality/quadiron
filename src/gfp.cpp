
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
bool GFP<T>::check(T a)
{
  return (a >= 0 && a < this->p);
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
  assert(check(a) && check(b));

  return (a + b) % this->p;
}

template <typename T>
T GFP<T>::sub(T a, T b)
{
  assert(check(a) && check(b));

  if (a >= b)
    return a - b;
  else
    return this->p - (b - a);
}

template <typename T>
T GFP<T>::mul(T a, T b)
{
  assert(check(a) && check(b));

  using DoubleT = typename Double<T>::T;

  return T((DoubleT(a) * b) % this->p);
}

template <typename T>
T GFP<T>::div(T a, T b)
{
  assert(check(a) && check(b));

  T inv_b = inv(b);

  return mul(a, inv_b);
}

template <typename T>
T GFP<T>::inv(T a)
{
  assert(check(a));

  return GF<T>::generic_mod_exp(this, a, this->p - 2, this->p);
}

template <typename T>
T GFP<T>::pow(T a)
{
  assert(check(a));

  return GF<T>::generic_mod_exp(this, this->p, a, this->p);
}

template <typename T>
T GFP<T>::log(T a)
{
  assert(check(a));

  return GF<T>::generic_trial_mult_log(this, this->p, a, this->p);
}

template <typename T>
T GFP<T>::weak_rand(void)
{
  T r;

 retry:
  r = rand() % this->p;
  if (0 == r)
    goto retry;
  return r;
}
