
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
 * see examples from http://web.eecs.utk.edu/~plank/plank/papers/CS-96-332.html
 * 
 */
void rs_utest1()
{
  GF2N<uint32_t> gf16(4);

  assert(gf16.mul(3, 7) == 9);
  assert(gf16.mul(13, 10) == 11);  
  assert(gf16.div(13, 10) == 3);
  assert(gf16.div(3, 7) == 10);

  Mat<uint32_t> mat(&gf16, 3, 3);
  Vec<uint32_t> vec(&gf16, 3);
  Vec<uint32_t> output(&gf16, 3);

  mat.vandermonde();
  VEC_ITEM(&vec, 0) = 3;
  VEC_ITEM(&vec, 1) = 13;
  VEC_ITEM(&vec, 2) = 9;
  mat.mult(&output, &mat, &vec);
  assert(VEC_ITEM(&output, 0) == 7);
  assert(VEC_ITEM(&output, 1) == 2);
  assert(VEC_ITEM(&output, 2) == 9);

  VEC_ITEM(&vec, 0) = 3;
  VEC_ITEM(&vec, 1) = 1;
  VEC_ITEM(&vec, 2) = 9;
  mat.mult(&output, &mat, &vec);
  assert(VEC_ITEM(&output, 0) == 11);
  assert(VEC_ITEM(&output, 1) == 9);
  assert(VEC_ITEM(&output, 2) == 12);

  MAT_ITEM(&mat, 0, 0) = 1;
  MAT_ITEM(&mat, 0, 1) = 0;
  MAT_ITEM(&mat, 0, 2) = 0;
  MAT_ITEM(&mat, 1, 0) = 1;
  MAT_ITEM(&mat, 1, 1) = 1;
  MAT_ITEM(&mat, 1, 2) = 1;
  MAT_ITEM(&mat, 2, 0) = 1;
  MAT_ITEM(&mat, 2, 1) = 2;
  MAT_ITEM(&mat, 2, 2) = 3;
  VEC_ITEM(&vec, 0) = 3;
  VEC_ITEM(&vec, 1) = 11;
  VEC_ITEM(&vec, 2) = 9;
  mat.inv();
  mat.mult(&output, &mat, &vec);
  assert(VEC_ITEM(&output, 1) == 1);
  assert(VEC_ITEM(&output, 2) == 9);
}

/** 
 * see http://web.eecs.utk.edu/~plank/plank/papers/CS-03-504.html
 */
void rs_utest2()
{
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

void rs_utest3()
{
  GF2N<uint32_t> gf16(4);
  Mat<uint32_t> mat(&gf16, 3, 3);

  mat.cauchy();
}

void rs_utest4()
{
  GF2N<uint32_t> gf256(8);

  assert(gf256.mul(3, 7) == 9);
  assert(gf256.mul(13, 10) == 114);
  assert(gf256.div(13, 10) == 40);
  assert(gf256.div(3, 7) == 211);
}

void rs_utest()
{
  rs_utest1();
  rs_utest2();
  rs_utest3();
  rs_utest4();
}
