
#include "ntl.h"
#include "gf.cpp"
#include "gfp.cpp"
#include "gf2n.cpp"
#include "mat.cpp"
#include "vec.cpp"

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;
template class Mat<uint32_t>;
template class Vec<uint32_t>;

/** 
 * see http://web.eecs.utk.edu/~plank/plank/papers/CS-03-504.html
 */
void rs_utest1()
{
  std::cout << "rs_utest1\n";

  GF2N<uint32_t> gf16(4);
  Mat<uint32_t> mat(&gf16, 3, 3);

  mat.vandermonde_suitable_for_ec();
  //mat.dump();

  assert(MAT_ITEM(&mat, 0, 0) == 1);
  assert(MAT_ITEM(&mat, 0, 1) == 1);
  assert(MAT_ITEM(&mat, 0, 2) == 1);
  assert(MAT_ITEM(&mat, 1, 0) == 15);
  assert(MAT_ITEM(&mat, 1, 1) == 8);
  assert(MAT_ITEM(&mat, 1, 2) == 6);
  assert(MAT_ITEM(&mat, 2, 0) == 14);
  assert(MAT_ITEM(&mat, 2, 1) == 9);
  assert(MAT_ITEM(&mat, 2, 2) == 6);
}

void rs_utest2()
{
  std::cout << "rs_utest2\n";

  GF2N<uint32_t> gf16(4);
  Mat<uint32_t> mat(&gf16, 3, 3);

  mat.cauchy();
}

void rs_utest3()
{
  std::cout << "rs_utest3\n";

  GF2N<uint32_t> gf256(8);

  assert(gf256.mul(3, 7) == 9);
  assert(gf256.mul(13, 10) == 114);
  assert(gf256.div(13, 10) == 40);
  assert(gf256.div(3, 7) == 211);
}

void rs_utest()
{
  std::cout << "rs_utest\n";

  rs_utest1();
  rs_utest2();
  rs_utest3();
}
