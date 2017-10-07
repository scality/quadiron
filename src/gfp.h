/* -*- mode: c++ -*- */
#pragma once
#include "gf.h"

template<typename T>
class GFP : public GF<T>
{
private:
  T p;

 public:
  explicit GFP(T p);
  T inv_exp(T a);
};

template <typename T>
GFP<T>::GFP(T p) : GF<T>(p, 1)
{
  this->p = p;
}

/**
 * Inverse by exponentiation
 *
 * @param a
 *
 * @return
 */
template <typename T>
T GFP<T>::inv_exp(T a)
{
  assert(this->check(a));

  return this->exp(a, this->p - 2);
}
