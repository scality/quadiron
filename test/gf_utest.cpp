
#include "ntl.h"
#include "gf.cpp"
#include "gfp.cpp"
#include "gf2n.cpp"

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;

void gf_utest1(GF<uint32_t> *gf)
{
  int i;
  
  for (i = 0; i < 100;i++) {
    uint32_t x, y;
    x = gf->weak_rand();
    //std::cout << "x=" << x << "\n";
    y = gf->inv(x);
    //std::cout << "inv(x)=" << y << "\n";
    assert(gf->mul(x, y) == 1);
  }
}

void gf_utest2(GF<uint32_t> *gf)
{
  int i;
  
  for (i = 0; i < 100;i++) {
    uint32_t x, y, z, t;
    x = gf->weak_rand();
    //std::cout << "x=" << x << "\n";
    y = gf->weak_rand();
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
  //std::cout << "gf5\n";
  gf_utest1(&gf5);
  gf_utest2(&gf5);
  //std::cout << "gf256\n";
  gf_utest1(&gf256);
  gf_utest2(&gf256);
}
