
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

template <typename T>
bool check_taylor_expansion(Poly<T> *f, T n, T t, std::vector<Poly<T>> *res) {
  Poly<T> y(f->gf), z(f->gf), g(f->gf);

  y.set(1, 1);
  y.set(t, 1);

  z.copy(&y);
  g.copy(&(res->at(0)));
  for (std::vector<int>::size_type i = 1; i < res->size(); i++) {
    Poly<T> tmp(f->gf);
    tmp.copy(&(res->at(i)));
    tmp.mul(&z);
    g.add(&tmp);
    z.mul(&y);
  }
  return f->equal(&g);
}

// taylor expansion on (x^t - x)
void poly_utest8()
{
  std::cout << "poly_utest8\n";
  GF2N<uint32_t> gf(8);
  Poly<uint32_t> p1(&gf);
  std::vector<Poly<uint32_t>> p2;
  uint32_t n = 8;
  uint32_t t = 2;
  p1.set(6, 2);
  p1.set(4, 4);
  p1.set(1, 1);
  p1.set(0, 8);
  // p1.dump();
  p1.taylor_expand(&p2, n, t);
  // std::cout << "taylor_expand done" << std::endl;
  // for(std::vector<int>::size_type i = 0; i != p2.size(); i++) {
  //   std::cout << i << ":";
  //   p2.at(i).dump();
  // }
  bool ok = check_taylor_expansion(&p1, n, t, &p2);
  assert(ok);
}

// taylor expansion on (x^2 - x)
void poly_utest9()
{
  std::cout << "poly_utest9\n";
  GF2N<uint32_t> gf(8);
  Poly<uint32_t> p1(&gf);
  Poly<uint32_t> G0(&gf);
  Poly<uint32_t> G1(&gf);
  uint32_t n = 8;
  uint32_t t = 2;
  p1.set(6, 2);
  p1.set(4, 4);
  p1.set(1, 1);
  p1.set(0, 8);
  // p1.dump();
  p1.taylor_expand_t2(&G0, &G1, n);
  std::vector<Poly<uint32_t>> p2;
  int deg = G0.degree() > G1.degree() ? G0.degree() : G1.degree();
  for (int i = 0; i <= deg; i++) {
    Poly<uint32_t> tmp(&gf);
    tmp.set(0, G0.get(i));
    tmp.set(1, G1.get(i));
    p2.push_back(tmp);
  }
  bool ok = check_taylor_expansion(&p1, n, t, &p2);
  assert(ok);
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
