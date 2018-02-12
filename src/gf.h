/* -*- mode: c++ -*- */
#ifndef __NTTEC_GF_H__
#define __NTTEC_GF_H__

#include <vector>

#include "arith.h"
#include "core.h"
#include "rn.h"

template <typename T>
class GFP;

/**
 * Generic class for Galois Fields of q=p^n
 *
 * @param p primer number
 * @param n exponent
 */
template <typename T>
class GF : public RN<T> {
  protected:
    T p;
    int n;
    GFP<T>* sub_field;

  public:
    GF(T p, int n);
    virtual ~GF();
    GF<T>* get_sub_field();
    T get_p();
    int get_n();
};

template <typename T>
GF<T>::GF(T p, int n) : RN<T>(_exp<T>(p, n))
{
    // XXX shall check that p is prime
    this->p = p;
    this->n = n;
    if (n == 1)
        this->sub_field = nullptr;
    else
        this->sub_field = new GFP<T>(p);
}

template <typename T>
GF<T>::~GF()
{
    if (sub_field)
        delete sub_field;
}

/**
 * return the field in which is based the extension field (or the field
 * itself if n == 1)
 */
template <typename T>
GF<T>* GF<T>::get_sub_field()
{
    if (this->sub_field)
        return this->sub_field;
    else
        return this;
}

template <typename T>
T GF<T>::get_p()
{
    return p;
}

template <typename T>
int GF<T>::get_n()
{
    return n;
}

#endif
