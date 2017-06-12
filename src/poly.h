/* -*- mode: c++ -*- */

#pragma once

template<typename T>
class Poly
{
public:

  struct Term: std::map <T, T>
  {
    bool is_key(T &s) const
    {
      return count(s) != 0;
    }
  };
  
  GF<T> *gf;

  Term terms;
  
private:

 public:
  Poly(GF<T> *gf);
  void clear();
  void copy(Poly<T> *src);
  T degree();
  T lead();
  bool is_zero();
  T get(T exponent);
  void set(T exponent, T coef);
  void add(Poly<T> *result, Poly<T> *a, Poly<T> *b);
  void sub(Poly<T> *result, Poly<T> *a, Poly<T> *b);
  void mul(Poly<T> *result, Poly<T> *a, Poly<T> *b);
  void div(Poly<T> *q, Poly<T> *r, Poly<T> *n, Poly<T> *d);
  void derivative(Poly<T> *result);
  void dump();
};

template <typename T>
Poly<T>::Poly(GF<T> *gf)
{
  this->gf = gf;
}

template <typename T>
void Poly<T>::clear()
{
  terms.clear();
}

template <typename T>
void Poly<T>::copy(Poly<T> *src)
{
  clear();

  for (int i = src->degree(); i >= 0; i--)
    set(i, src->get(i));
}

template <typename T>
T Poly<T>::degree()
{
  return (terms.rbegin() == terms.rend()) ? 0 : terms.rbegin()->first;
}

template <typename T>
T Poly<T>::lead()
{
  return get(degree());
}

template <typename T>
bool Poly<T>::is_zero()
{
  return degree() == 0 && get(0) == 0;
}

template <typename T>
T Poly<T>::get(T exponent)
{
  typename Term::const_iterator it = terms.find(exponent);
  
  return (it == terms.end()) ? 0 : it->second;
}

template <typename T>
void Poly<T>::set(T exponent, T coef)
{
  typename Term::const_iterator it = terms.find(exponent);

  if (it == terms.end() && coef == 0)
    return ;
  terms[exponent] = coef;
}

template <typename T>
void Poly<T>::add(Poly<T> *result, Poly<T> *a, Poly<T> *b)
{
  result->clear();

  T max = gf->max(a->degree(), b->degree());

  for (int i = max; i >= 0; i--)
    result->set(i, 
                gf->add(a->get(i), b->get(i)));
}

template <typename T>
void Poly<T>::sub(Poly<T> *result, Poly<T> *a, Poly<T> *b)
{
  result->clear();

  T max = gf->max(a->degree(), b->degree());

  for (int i = max; i >= 0; i--)
    result->set(i, 
                gf->sub(a->get(i), b->get(i)));
}

template <typename T>
void Poly<T>::mul(Poly<T> *result, Poly<T> *a, Poly<T> *b)
{
  result->clear();

  for (int i = a->degree(); i >= 0; i--)
    for (int j = b->degree(); j >= 0; j--)
      result->set(i + j,
                  gf->add(result->get(i + j),
                          gf->mul(a->get(i), b->get(j))));
}

/** 
 * long division algorithm (source Wikipedia)
 * 
 * @param q quotient
 * @param r remainder
 * @param n dividend
 * @param d divisor
 */
template <typename T>
void Poly<T>::div(Poly<T> *q, Poly<T> *r, Poly<T> *n, Poly<T> *d)
{
  if (d->is_zero())
    throw NTL_EX_DIV_BY_ZERO;

  q->clear();
  r->copy(n);
  while (!r->is_zero() && (r->degree() >= d->degree())) {
    Poly<T> t(gf), _q(gf), _r(gf), _m(gf);
    t.set(r->degree() - d->degree(), 
          gf->div(r->lead(), d->lead()));
    add(&_q, q, &t);
    q->copy(&_q);
    mul(&_m, &t, d);
    sub(&_r, r, &_m);
    r->copy(&_r);
  }
}

template <typename T>
void Poly<T>::derivative(Poly<T> *result)
{
  result->clear();

  for (int i = degree(); i > 0; i--)
    result->set(i - 1, gf->mul(get(i), i));
}

template <typename T>
void Poly<T>::dump()
{
  typename Term::const_reverse_iterator it = terms.rbegin();
  
  for (;it != terms.rend(); it++)
    std::cout << " " << it->second << "x^" << it->first;
  std::cout << "\n";
}
