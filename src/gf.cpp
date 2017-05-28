/**
 * @file   gf.cpp
 * @author vr <vr@vr-VirtualBox>
 * @date   Thu May 25 08:43:45 2017
 * 
 * @brief  Generic arithmetic routines
 * 
 * @note   Routines starting with underscore are not using GF functions
 */

#include "ntl.h"

template <typename T>
GF<T>::GF(T p, T n)
{
  this->p = p;
  this->n = n;
}

template <typename T>
T GF<T>::_card()
{
  return _exp(p, n);
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
T GF<T>::exp(T base, T exponent)
{
  T result;
  T i;

  if (0 == exponent)
    return 1;

  if (1 == exponent)
    return base;

  result = base;
  for (i = 1; i < exponent;i++) {
    result = this->mul(result, base);
  }

  return result;
}

/** 
 * Regular exponentation 
 * 
 * @param gf 
 * @param base 
 * @param exponent 
 * @param modulus 
 * 
 * @return
 */
template <typename T>
T GF<T>::_exp(T base, T exponent)
{
  T result = 1;
  while (exponent > 0) {
    if (exponent % 2 == 1)
      result *= base;
    exponent >>= 1;
    base *= base;
  }

  return result;
}

/** 
 * Modular exponentation taken from Applied Cryptography by Bruce Schneier.
 * 
 * @param gf 
 * @param base 
 * @param exponent 
 * @param modulus 
 * 
 * @return
 */
template <typename T>
T GF<T>::_mod_exp(T base, T exponent, T modulus)
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
T GF<T>::_trial_mult_log(T base, T exponent, T modulus)
{
  T result;

  for (result = 1;result < card();result++) {
    if (_mod_exp(base, result, modulus) == exponent)
      return result;
  }

  //not found
  throw NTL_EX_NOT_FOUND;
}

template <typename T>
T GF<T>::_log2(T exponent)
{
  T result = 0;

  if (exponent == 0) 
    throw NTL_EX_INVAL;

  if (exponent == 1) 
    return 0;

  while (exponent > 1) {
    exponent >>= 1;
    result++;
  }

  return result;
}

/** 
 * Extended Euclidean algorithm (from Wikipedia)
 * Perform GCD of a and b
 *
 * @param a
 * @param b
 * @param bezout_coef[output] computed bezout coefficients, might be NULL
 * @param quotient_gcd[output] computed quotient by the GCD, might be NULL
 * 
 * XXX take care of the signs of input
 *
 * @return the GCD
 */
template <typename T>
SignedDoubleT<T> GF<T>::_extended_gcd(SignedDoubleT<T> a, SignedDoubleT<T> b, SignedDoubleT<T> bezout_coef[2], SignedDoubleT<T> quotient_gcd[2])
{
  SignedDoubleT<T> s = 0;
  SignedDoubleT<T> old_s = 1;
  SignedDoubleT<T> t = 1;
  SignedDoubleT<T> old_t = 0;
  SignedDoubleT<T> r = b;
  SignedDoubleT<T> old_r = a;
  SignedDoubleT<T> quotient;
  SignedDoubleT<T> tmp;
  
  while (0 != r) {
    quotient = old_r / r;

    tmp = r;
    r = old_r - quotient * r;
    old_r = tmp;

    tmp = s;
    s = old_s - quotient * s;
    old_s = tmp;

    tmp = t;
    t = old_t - quotient * t;
    old_t = tmp;
  }
  
  if (bezout_coef) {
    bezout_coef[0] = old_s;
    bezout_coef[1] = old_t;
  }

  if (quotient_gcd) {
    quotient_gcd[0] = t;
    quotient_gcd[1] = s;
  }

  return old_r;
}

/** 
 * Chinese remainder theorem (from Wikipedia)
 * 
 * @param n_mod number of moduli
 * @param a the a's (integers)
 * @param n the n's (moduli)
 * 
 * XXX check if there is a solution
 *
 * @return the solution
 * @throw NTL_EX_NO_SOLUTION
 */
template <typename T>
T GF<T>::_chinese_remainder(int n_mod, T a[], T n[])
{
  int i;
  SignedDoubleT<T> _N, x;
  SignedDoubleT<T> *N = new SignedDoubleT<T>[n_mod];
  SignedDoubleT<T> *M = new SignedDoubleT<T>[n_mod];

  _N = 1;
  for (i = 0;i < n_mod;i++)
    _N *= n[i];

  for (i = 0;i < n_mod;i++) {
    SignedDoubleT<T> bezout[2];

    N[i] = _N / n[i];
    _extended_gcd(N[i], n[i], bezout, NULL);
    M[i] = bezout[0]; //XXX % n[i];
  }
  
  x = 0;
  for (i = 0;i < n_mod;i++)
    x += N[i] * M[i] * a[i];

  delete[] N;
  delete[] M;

  x = x % _N;

  if (x < 0)
    x = _N + x;
  
  return x;
}

/** 
 * check if q is a quadractic residue in GF_p^n
 * 
 * @param q 
 * 
 * @return boolean
 */
template <typename T>
bool GF<T>::is_quadratic_residue(T q)
{
  T i;

  for (i = 0; i < card();i++) {
    if (exp(i, 2) == q)
      return true;
  }

  return false;
}

/** 
 * Compute the jacobi symbol of the 2 numbers
 * https://groups.google.com/forum/#!topic/sci.crypt/v9_cNF06XjU
 * 
 * @param n
 * @param m 
 * 
 * @return the jacobi symbol
 * @throw NT_EX_INVAL if b is even or b is negative
 */
template <typename T>
int GF<T>::_jacobi(SignedDoubleT<T> n, SignedDoubleT<T> m)
{
  SignedDoubleT<T> t;
  int jac;

  //m must be odd
  if ((m % 2) == 0) 
    throw NTL_EX_INVAL;

  //m must be positive
  if (m < 0)
    throw NTL_EX_INVAL;
  
  jac = 1;
  while (m != 1) {
    //if the gcd of $n$ and $m$ is $>1$ Jacobi returns $0$
    if ((n == 0) || ((n % 2 == 0) && (m % 2 == 0))) {
      jac = 0;
      m = 1;
    }
    //$J(n,2*m) = J(n,m) * J(n,2) = J(n,m) * (-1)^{(n^2-1)/8}$
    else if (m % 2 == 0) {
      if ((n % 8 == 3) || (n % 8 == 5)) 
        jac = -jac;
      m = m / 2;
    }
    //$J(2*n,m) = J(n,m) * J(2,m) = J(n,m) * (-1)^{(m^2-1)/8}$
    else if (n % 2 == 0) {
      if ((m % 8 == 3) || (m % 8 == 5))
        jac = -jac;
      n = n / 2;
    }
    //$J(-n,m) = J(n,m) * J(-1,m) = J(n,m) * (-1)^{(m-1)/2}$
    else if (n < 0) {
      if (m % 4 == 3)
        jac = -jac;
      n = -n;
    }
    //$J(n,m) = J(m,n) * (-1)^{(n-1)*(m-1)/4}$ (quadratic reciprocity)
    else {
      if ((n % 4 == 3) && (m % 4 == 3))
        jac = -jac;
      t = n;
      n = m % n;
      m = t;
    }
  }

  return jac;
}

template <typename T>
bool GF<T>::_solovay_strassen1(T a, T n)
{
  T _n = (n - 1) / 2;
  T _j = _mod_exp(a, _n, n);
  int j = _jacobi(a, n);
  if (j < 0)
    _j = _j - n;
  return j == _j;
}

/** 
 * Perform the Solvay-Strassen primality test
 * 
 * @param gf 
 * @param n check if n is prime
 * 
 * @return true if n is prime else false
 */
template <typename T>
bool GF<T>::_solovay_strassen(T n)
{
  int ok = 0;
  for (int i = 0;i < 100;i++) {
    T a = _weak_rand(n);
    if (!_solovay_strassen1(a, n))
      return false;
  }
  return true;
}

/** 
 * Returns a number n such as 0 < n < max
 * 
 * @param max 
 * 
 * @return 
 */
template <typename T>
T GF<T>::_weak_rand(T max)
{
  T r;

 retry:
  r = rand() % max;
  if (0 == r)
    goto retry;
  return r;
}

template <typename T>
T GF<T>::weak_rand(void)
{
  return _weak_rand(card());
}
