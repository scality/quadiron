
#include "ntl.h"

template <typename T>
Vec<T>::Vec(GF<T> *gf, int n)
{
  this->gf = gf;
  this->n = n;
  this->mem = new T[n];
}

template <typename T>
void Vec<T>::zero(void)
{
  int i;

  for (i = 0;i < n;i++)
    VEC_ITEM(this, i) = 0;
}
