
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
      if (1 == (b & 1))
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

/** 
 * Naive exponentiation
 * Note: honor Galois field multiplication 
 *
 * @param gf 
 * @param base 
 * @param exponent 
 * 
 * @return 
 */
template <typename T>
T GF<T>::generic_naive_exp(GF *gf, T base, T exponent)
{
  T result;
  T i;

  if (0 == exponent)
    return 1;

  if (1 == exponent)
    return base;

  result = base;
  for (i = 1; i < exponent;i++) {
    result = gf->mul(result, base);
  }

  return result;
}

/** 
 * Modular exponentation taken from Applied Cryptography by Bruce Schneier.
 * Note: does not respect Galois field multiplication
 * 
 * @param gf 
 * @param base 
 * @param exponent 
 * @param modulus 
 * 
 * @return
 */
template <typename T>
T GF<T>::generic_mod_exp(GF *gf, T base, T exponent, T modulus)
{
  if (1 == modulus)
    return 0;
  //XXX assert (modulus - 1) * (modulus - 1) does not overflow base
  T result = 1;
  base = base % modulus;
  while (exponent > 0) {
    if (exponent % 2 == 1)
      result = (result * base) % modulus;
    exponent = exponent >> 1;
    base = (base * base) % modulus;
  }

  return result;
}

/** 
 * Naive brute force algorithm
 * 
 * @param gf 
 * @param base 
 * @param exponent 
 * @param modulus 
 * 
 * @throw NTL_EX_NOT_FOUND if result is not found
 * return
 */
template <typename T>
T GF<T>::generic_trial_mult_log(GF *gf, T base, T exponent, T modulus)
{
  T result;

  for (result = 1;result < card();result++) {
    if (generic_mod_exp(gf, base, result, modulus) == exponent)
      return result;
  }

  //not found
  throw NTL_EX_NOT_FOUND;
}
