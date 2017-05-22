
#ifndef __GF_H__
#define __GF_H__ 1

template<typename T>
class GF
{
 protected:
  T p;
  T n;
  GF(T p, T n);

 public:
  T generic_card(GF *gf);
  T generic_naive_exp(GF *gf, T base, T exponent);
  T generic_mod_exp(GF *gf, T base, T exponent, T modulus);
  T generic_trial_mult_log(GF *gf, T base, T exponent, T modulus);
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
  virtual T weak_rand(void) = 0;
};

extern void gf_utest();

#endif
