
#include "ntl.h"

template <typename T>
Vec<T>::Vec(GF<T> *gf, int n)
{
  this->gf = gf;
  this->n;
  this->mem = new T[n];
}
