/* -*- mode: c++ -*- */
#ifndef __NTL_VVEC_H__
#define __NTL_VVEC_H__

#include "vec.h"

/**
 * Virtual vector that returns 0 beyond managed vector length
 */
template <typename T>
class VVec : public Vec<T> {
  private:
    Vec<T>* vec;
    int n;

  public:
    VVec(Vec<T>* vec, int n);
    int get_n(void);
    T get(int i);
};

template <typename T>
VVec<T>::VVec(Vec<T>* vec, int n)
    : Vec<T>(vec->rn, n, vec->get_mem(), vec->get_mem_len())
{
    this->vec = vec;
    this->n = n;
}

template <typename T>
int VVec<T>::get_n(void)
{
    return n;
}

template <typename T>
T VVec<T>::get(int i)
{
    assert(i >= 0 && i < n);

    return (i < vec->get_n()) ? vec->get(i) : 0;
}

#endif
