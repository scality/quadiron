
#include "main.h"

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
  throw NTLEX_NOT_FOUND;
}

template class GF<uint32_t>;

void gf_utest1(GF<uint32_t> *gf)
{
  int i;
  
  for (i = 0; i < 100;i++) {
    //pick a non-zero random number from the field
  retry:
    uint32_t x, y, z, t;
    x = gf->weak_rand();
    if (0 == x)
      goto retry;
    //std::cout << "x=" << x << "\n";
    y = gf->inv(x);
    //std::cout << "inv(x)=" << y << "\n";
    assert(gf->mul(x, y) == 1);
    //pick another non-zero random number from the field
  retry2:
    y = gf->weak_rand();
    if (0 == y)
      goto retry2;
    //std::cout << "y=" << y << "\n";
    try {
      z = gf->log(x, y);
    } catch (...) {
      //std::cout << "not found\n";
      continue ;
    }
    //std::cout << "log" << x << "(" << y << ")=" << z << "\n";
    t = gf->pow(x, z);
    //std::cout << x << "^" << z << "=" << t << "\n";
    assert(t == y);
  }
}

void gf_utest()
{
  srand(time(0));
  GFP<uint32_t> gf5(5);
  GF2N<uint32_t> gf256(8);
  std::cout << "gf5\n";
  gf_utest1(&gf5);
  std::cout << "gf256\n";
  gf_utest1(&gf256);
}
