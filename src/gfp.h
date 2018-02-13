/* -*- mode: c++ -*- */
#ifndef __NTTEC_GFP_H__
#define __NTTEC_GFP_H__

#include "gf.h"

namespace nttec {
namespace gf {

template <typename T>
class GFP : public GF<T> {
  public:
    explicit GFP(T p);
    T inv_exp(T a);
};

template <typename T>
GFP<T>::GFP(T p) : GF<T>(p, 1)
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
T GFP<T>::inv_exp(T a)
{
    assert(this->check(a));

    return this->exp(a, this->p - 2);
}

} // namespace gf
} // namespace nttec

#endif
