
#include "main.h"

template <typename T>
GF2N<T>::GF2N(T n) : GF<T>(2, n)
{
  if (n == 4)
    this->prim_poly = 0x13;
  if (n == 8)
    this->prim_poly = 0x11d;
  else if (n == 16)
    this->prim_poly = 0x1100b;
  else
    assert(false); //XXX generate polynomial

  T one = 1;
  this->my_card = one << n;
  setup_tables();
}

template <typename T>
void GF2N<T>::setup_tables(void)
{
  T b, log;

  gflog = new T[my_card];
  gfilog = new T[my_card];
  b = 1;
  for (log = 0; log < my_card - 1; log++) {
    gflog[b] = log;
    gfilog[log] = b;
    b = b << 1;
    if (b & my_card) 
      b = b ^ prim_poly;
  }
}

template <typename T>
T GF2N<T>::card(void)
{
  return my_card;
}

template <typename T>
T GF2N<T>::add(T a, T b)
{
  return a ^ b;
}

template <typename T>
T GF2N<T>::sub(T a, T b)
{
  return a ^ b;
}

template <typename T>
T GF2N<T>::mul(T a, T b)
{
  int sum_log;

  if (a == 0 || b == 0)
    return 0;
  sum_log = gflog[a] + gflog[b];
  if (sum_log >= my_card - 1) 
    sum_log -= my_card - 1;

  return gfilog[sum_log];
}

template <typename T>
T GF2N<T>::div(T a, T b)
{
  int diff_log;

  if (a == 0) 
    return 0;
  if (b == 0) 
    return -1;
  diff_log = gflog[a] - gflog[b];
  if (diff_log < 0) 
    diff_log += my_card - 1;

  return gfilog[diff_log];
}

template <typename T>
T GF2N<T>::inv(T a)
{
  return div(1, a);
}

template <typename T>
T GF2N<T>::pow(T a, T b)
{
  return GF<T>::generic_pow(this, a, b);
}

template <typename T>
T GF2N<T>::log(T a, T b)
{
  return GF<T>::generic_trial_mult_log(this, a, b);
}

template <typename T>
T GF2N<T>::weak_rand(void)
{
  return rand() % my_card;
}

template class GF2N<uint32_t>;
