/* -*- mode: c++ -*- */
#pragma once

template<typename T>
class GFP;

/**
 * Generic class for Galois Fields of q=p^n
 *
 * @param p primer number
 * @param n exponent
 */
template<typename T>
class GF : public CG<T>
{
private:
  T p;
  T n;
  GFP<T> *sub_field;

public:
  GF(T p, T n);
  virtual ~GF();
  GF<T> *get_sub_field();
};

template <typename T>
GF<T>::GF(T p, T n) : CG<T>(this->arith->exp(p, n))
{
  // XXX shall check that p is prime
  this->p = p;
  this->n = n;
  if (n == 1)
    this->sub_field = NULL;
  else
    this->sub_field = new GFP<T>(p);
}

template <typename T>
GF<T>::~GF()
{
  if (sub_field) delete sub_field;
}

/** 
 * return the field in which is based the extension field (or the field
 * itself if n == 1)
 */
template <typename T>
GF<T> *GF<T>::get_sub_field()
{
  if (this->sub_field)
    return this->sub_field;
  else
    return this;
}
