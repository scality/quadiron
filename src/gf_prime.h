/* -*- mode: c++ -*- */
#ifndef __NTTEC_GF_PRIME_H__
#define __NTTEC_GF_PRIME_H__

#include "gf_base.h"

namespace nttec {
namespace gf {

/** A Galois Field whose order is a prime number. */
template <typename T>
class Prime : public gf::Field<T> {
  public:
    explicit Prime(T p);
    T inv_exp(T a);
};

template <typename T>
Prime<T>::Prime(T p) : gf::Field<T>(p, 1)
{
}

/**
 * Inverse by exponentiation
 *
 * @param a
 *
 * @return
 */
template <typename T>
T Prime<T>::inv_exp(T a)
{
    assert(this->check(a));

    return this->exp(a, this->p - 2);
}

} // namespace gf
} // namespace nttec

#endif
