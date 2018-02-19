/* -*- mode: c++ -*- */
#ifndef __NTTEC_GF_BASE_H__
#define __NTTEC_GF_BASE_H__

#include "arith.h"
#include "core.h"
#include "rn.h"

namespace nttec {

/** Galois Fields handling. */
namespace gf {

template <typename T>
class Prime;

/** Generic class for Galois Fields of q=p<sup>n/<sup>. */
template <typename T>
class Field : public RingModN<T> {
  public:
    Field(T p, int n);
    virtual ~Field();
    Field<T>* get_sub_field();
    T get_p();
    int get_n();

  protected:
    T p;
    int n;
    Prime<T>* sub_field;
};

template <typename T>
Field<T>::Field(T p, int n) : RingModN<T>(arith::exp<T>(p, n))
{
    // XXX shall check that p is prime
    this->p = p;
    this->n = n;
    if (n == 1)
        this->sub_field = nullptr;
    else
        this->sub_field = new Prime<T>(p);
}

template <typename T>
Field<T>::~Field()
{
    if (sub_field)
        delete sub_field;
}

/**
 * return the field in which is based the extension field (or the field
 * itself if n == 1)
 */
template <typename T>
Field<T>* Field<T>::get_sub_field()
{
    if (this->sub_field)
        return this->sub_field;
    else
        return this;
}

template <typename T>
T Field<T>::get_p()
{
    return p;
}

template <typename T>
int Field<T>::get_n()
{
    return n;
}

} // namespace gf
} // namespace nttec

#endif
