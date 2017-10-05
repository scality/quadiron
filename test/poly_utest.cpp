
#include "ntl.h"

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;
template class Mat<uint32_t>;
template class Vec<uint32_t>;

void poly_utest1()
{
  GFP<uint32_t> gfp(1000);

  Poly<uint32_t> p0(&gfp);
  assert(0 == p0.degree());
}

/**
 * Examples taken from https://www.doc.ic.ac.uk/~mrh/330tutor/ch04s02.html
 */
void poly_utest2()
{
  std::cout << "poly_utest2\n";

  GFP<uint32_t> gfp(11);
  Poly<uint32_t> p1(&gfp);
  p1.set(5, 1);
  p1.set(3, 3);
  p1.set(0, 4);
  // p1.dump();
  Poly<uint32_t> p2(&gfp);
  p2.set(6, 6);
  p2.set(3, 4);
  // p2.dump();
  Poly<uint32_t> p3(&gfp);
  p1._add(&p3, &p1, &p2);
  // p3.dump();
  assert(p3.degree() == 6);
  assert(p3.get(6) == 6);
  assert(p3.get(5) == 1);
  assert(p3.get(3) == 7);
  assert(p3.get(0) == 4);
}

void poly_utest3()
{
  std::cout << "poly_utest3\n";

  GFP<uint32_t> gfp(11);
  Poly<uint32_t> p1(&gfp);
  p1.set(5, 1);
  p1.set(3, 3);
  p1.set(0, 4);
  // p1.dump();
  Poly<uint32_t> p2(&gfp);
  p2.set(6, 6);
  p2.set(3, 4);
  // p2.dump();
  Poly<uint32_t> p3(&gfp);
  p1._sub(&p3, &p1, &p2);
  assert(p3.degree() == 6);
  assert(p3.get(6) == 5);
  assert(p3.get(5) == 1);
  assert(p3.get(3) == 10);
  assert(p3.get(0) == 4);
}

void poly_utest4()
{
  std::cout << "poly_utest4\n";

  GFP<uint32_t> gfp(11);
  Poly<uint32_t> p1(&gfp);
  p1.set(5, 1);
  p1.set(3, 3);
  p1.set(0, 4);
  // p1.dump();
  Poly<uint32_t> p2(&gfp);
  p2.set(6, 6);
  p2.set(3, 4);
  // p2.dump();
  Poly<uint32_t> p3(&gfp);
  p1._mul(&p3, &p1, &p2);
  assert(p3.degree() == 11);
  assert(p3.get(11) == 6);
  assert(p3.get(9) == 7);
  assert(p3.get(8) == 4);
  assert(p3.get(6) == 3);
  assert(p3.get(3) == 5);
}

void poly_utest5()
{
  std::cout << "poly_utest5\n";
  GFP<uint32_t> gfp(11);
  Poly<uint32_t> p1(&gfp);
  p1.set(6, 3);
  p1.set(4, 7);
  p1.set(3, 4);
  p1.set(0, 5);
  // p1.dump();
  Poly<uint32_t> p2(&gfp);
  p2.set(4, 1);
  p2.set(3, 3);
  p2.set(0, 4);
  // p2.dump();
  Poly<uint32_t> p3(&gfp);
  Poly<uint32_t> p4(&gfp);
  p1._div(&p3, &p4, &p1, &p2);
  // p3.dump();
  // p4.dump();
  assert(p3.degree() == 2);
  assert(p3.get(2) == 3);
  assert(p3.get(1) == 2);
  assert(p3.get(0) == 1);
  assert(p4.degree() == 3);
  assert(p4.get(3) == 1);
  assert(p4.get(2) == 10);
  assert(p4.get(1) == 3);
  assert(p4.get(0) == 1);
}

void poly_utest6()
{
  std::cout << "poly_utest6\n";
  GFP<uint32_t> gfp(11);
  Poly<uint32_t> p1(&gfp);
  p1.set(4, 4);
  p1.set(1, 1);
  p1.set(0, 8);
  Poly<uint32_t> p2(&gfp);
  p1._derivative(&p2, &p1);
  assert(p2.degree() == 3);
  assert(p2.get(3) == 5);
  assert(p2.get(0) == 1);
}

void poly_utest7()
{
  std::cout << "poly_utest7\n";
  GFP<uint32_t> gfp(11);
  Poly<uint32_t> p1(&gfp);
  p1.set(4, 4);
  p1.set(1, 1);
  p1.set(0, 8);
  uint32_t result = p1.eval(3);
  assert(result == 5);
}

void taylor_expand(GF2N<uint32_t> *gf) {
  Poly<uint32_t> p1(gf);
  std::vector<Poly<uint32_t>> p2;
  uint32_t n = (gf->weak_rand() % 30) + 6;
  uint32_t t = (gf->weak_rand() % 6) + 2;
  for (int i = 0; i < n; i++)
    p1.set(i, gf->weak_rand());
  // p1.dump();
  p1.taylor_expand(&p2, n, t);
  Poly<uint32_t> _p1(gf);
  _p1.inv_taylor_expand(&p2, t);
  // _p1.dump();
  assert(_p1.equal(&p1));
}

// taylor expansion on (x^t - x)
void poly_utest8()
{
  std::cout << "poly_utest8\n";
  GF2N<uint32_t> gf(8);
  for (int i = 0; i < 1000; i++)
    taylor_expand(&gf);
}

void taylor_expand_t2(GF2N<uint32_t> *gf) {
  Poly<uint32_t> p1(gf);
  uint32_t n = (gf->weak_rand() % 30) + 6;
  uint32_t t = 2;
  uint32_t m = n/2;
  if (m * 2 < n) m++;
  Vec<uint32_t> G0(gf, m);
  Vec<uint32_t> G1(gf, m);
  for (int i = 0; i < n; i++)
    p1.set(i, gf->weak_rand());
  // p1.dump();
  p1.taylor_expand_t2(&G0, &G1, n);
  Poly<uint32_t> _p1(gf);
  _p1.inv_taylor_expand_t2(&G0, &G1);
  // _p1.dump();
  assert(_p1.equal(&p1));
}

// taylor expansion on (x^2 - x)
void poly_utest9()
{
  std::cout << "poly_utest9\n";
  GF2N<uint32_t> gf(8);
  for (int i = 0; i < 1000; i++)
    taylor_expand_t2(&gf);
}

void poly_utest()
{
  std::cout << "poly_utest\n";

  poly_utest1();
  poly_utest2();
  poly_utest3();
  poly_utest4();
  poly_utest5();
  poly_utest6();
  poly_utest7();
  poly_utest8();
  poly_utest9();
}
