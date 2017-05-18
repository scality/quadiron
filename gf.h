
#ifndef __GF_H__
#define __GF_H__ 1

template<typename T>
class GF
{
 protected:
  T p;
  T n;
  GF(T p, T n);
  T generic_card(GF *gf);
  T generic_pow(GF *gf, T a, T b);
  T generic_trial_mult_log(GF *gf, T a, T b);

 public:
  virtual T card(void) {};
  virtual T add(T a, T b) {};
  virtual T sub(T a, T b) {};
  virtual T mul(T a, T b) {};
  virtual T div(T a, T b) {};
  virtual T pow(T a, T b) {};
  virtual T log(T a, T b) {};
  virtual T inv(T a) {};
  virtual T weak_rand(void) {};
};

extern void gf_utest();

#endif
