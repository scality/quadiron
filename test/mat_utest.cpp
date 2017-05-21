
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

void mat_utest1()
{
  GFP<uint32_t> gf11(11);
  Mat<uint32_t> mat(&gf11, 3, 3);

  MAT_ITEM(&mat, 0, 0) = 2;
  MAT_ITEM(&mat, 0, 1) = 1;
  MAT_ITEM(&mat, 0, 2) = 2;
  MAT_ITEM(&mat, 1, 0) = 1;
  MAT_ITEM(&mat, 1, 1) = 2;
  MAT_ITEM(&mat, 1, 2) = 9;
  MAT_ITEM(&mat, 2, 0) = 1;
  MAT_ITEM(&mat, 2, 1) = 2;
  MAT_ITEM(&mat, 2, 2) = 7;
  
  mat.inv();

  assert(MAT_ITEM(&mat, 0, 0) == 8);
  assert(MAT_ITEM(&mat, 0, 1) == 6);
  assert(MAT_ITEM(&mat, 0, 2) == 1);
  assert(MAT_ITEM(&mat, 1, 0) == 7);
  assert(MAT_ITEM(&mat, 1, 1) == 9);
  assert(MAT_ITEM(&mat, 1, 2) == 10);
  assert(MAT_ITEM(&mat, 2, 0) == 0);
  assert(MAT_ITEM(&mat, 2, 1) == 6);
  assert(MAT_ITEM(&mat, 2, 2) == 5);
}

void mat_utest2()
{
  GFP<uint32_t> gf29(29);
  Mat<uint32_t> mat(&gf29, 3, 3);

  MAT_ITEM(&mat, 0, 0) = 22;
  MAT_ITEM(&mat, 0, 1) = 27;
  MAT_ITEM(&mat, 0, 2) = 18;
  MAT_ITEM(&mat, 1, 0) = 18;
  MAT_ITEM(&mat, 1, 1) = 28;
  MAT_ITEM(&mat, 1, 2) = 5;
  MAT_ITEM(&mat, 2, 0) = 4;
  MAT_ITEM(&mat, 2, 1) = 17;
  MAT_ITEM(&mat, 2, 2) = 1;
  
  mat.inv();

  assert(MAT_ITEM(&mat, 0, 0) == 1);
  assert(MAT_ITEM(&mat, 0, 1) == 18);
  assert(MAT_ITEM(&mat, 0, 2) == 8);
  assert(MAT_ITEM(&mat, 1, 0) == 2);
  assert(MAT_ITEM(&mat, 1, 1) == 8);
  assert(MAT_ITEM(&mat, 1, 2) == 11);
  assert(MAT_ITEM(&mat, 2, 0) == 20);
  assert(MAT_ITEM(&mat, 2, 1) == 24);
  assert(MAT_ITEM(&mat, 2, 2) == 14);
}

void mat_utest()
{
  mat_utest1();
  mat_utest2();
}
