
#include "ntl.h"

template <typename T>
Vec<T>::Vec(GF<T> *gf, int n)
{
  this->gf = gf;
  this->n = n;
  this->mem = new T[n];
}

template <typename T>
Vec<T>::~Vec()
{
  delete[] this->mem;
}

template <typename T>
void Vec<T>::zero_fill(void)
{
  int i;

  for (i = 0;i < n;i++)
    VEC_ITEM(this, i) = 0;
}

template <typename T>
void Vec<T>::dump(void)
{
  int i, j;
  
  std::cout << "--\n";
  for (i = 0;i < n;i++)
    std::cout << " " << VEC_ITEM(this, i);
  std::cout << "\n";
}
