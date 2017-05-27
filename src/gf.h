
#ifndef __GF_H__
#define __GF_H__ 1

template<typename T>
using DoubleT = typename Double<T>::T;

template<typename T>
using SignedDoubleT = typename SignedDouble<T>::T;

/** 
 * Generic class for Galois Fields of p^n
 * 
 * @note Functions starting with _ are convenience function that does
 *       not call GF specific operations
 *
 * @param p primer number
 * @param n exponent
 */
template<typename T>
class GF
{
 protected:
  T p;
  T n;
  GF(T p, T n);

 public:
  virtual T card(void) = 0;
  virtual bool check(T a) = 0;
  virtual T neg(T a) = 0;
  virtual T add(T a, T b) = 0;
  virtual T sub(T a, T b) = 0;
  virtual T mul(T a, T b) = 0;
  virtual T div(T a, T b) = 0;
  virtual T pow(T a) = 0;
  virtual T log(T a) = 0;
  virtual T inv(T a) = 0;
  T _card();
  T exp(T base, T exponent);
  T _exp(T base, T exponent);
  T _mod_exp(T base, T exponent, T modulus);
  T _trial_mult_log(T base, T exponent, T modulus);
  T _log2(T exponent);
  SignedDoubleT<T> _extended_gcd(SignedDoubleT<T> a, SignedDoubleT<T> b, SignedDoubleT<T> bezout_coef[2], SignedDoubleT<T> quotient_gcd[2]);
  T _chinese_remainder(int n_mod, T a[], T n[]);
  bool is_quadratic_residue(T q);
  int _jacobi(SignedDoubleT<T> n, SignedDoubleT<T> m);
  T weak_rand(void);
};

#endif
