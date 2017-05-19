
#include "ntl.h"

template <typename T>
Mat<T>::Mat(GF<T> *gf, int n_rows, int n_cols)
{
  this->gf = gf;
  this->n_rows = n_rows;
  this->n_cols = n_cols;
  this->mem = new T[n_rows * n_cols];
}
