/* -*- mode: c++ -*- */

#pragma once

template<typename T>
class Poly
{
 public:
  struct Term: std::map <T, T> {};

  GF<T> *gf;
  Term terms;

  explicit Poly(GF<T> *gf);
  void clear();
  void copy(Poly<T> *src);
  T degree();
  T lead();
  bool is_zero();
  T get(T exponent);
  void set(T exponent, T coef);
  void _neg(Poly<T> *result, Poly<T> *a);
  void _add(Poly<T> *result, Poly<T> *a, Poly<T> *b);
  void _sub(Poly<T> *result, Poly<T> *a, Poly<T> *b);
  void _mul(Poly<T> *result, Poly<T> *a, Poly<T> *b);
  void _div(Poly<T> *q, Poly<T> *r, Poly<T> *n, Poly<T> *d);
  void _derivative(Poly<T> *result, Poly<T> *a);
  void neg();
  void add(Poly<T> *b);
  void sub(Poly<T> *b);
  void mul(Poly<T> *b);
  void div(Poly<T> *d);
  void mod(Poly<T> *d);
  void derivative();
  T eval(T x);
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
  assert(gf->check(coef));

  typename Term::const_iterator it = terms.find(exponent);

  if (it == terms.end() && coef == 0)
    return;
  terms[exponent] = coef;
}

template <typename T>
void Poly<T>::_neg(Poly<T> *result, Poly<T> *a)
{
  result->clear();

  Poly<T> b(gf);
  sub(result, &b, a);
}

template <typename T>
void Poly<T>::_add(Poly<T> *result, Poly<T> *a, Poly<T> *b)
{
  result->clear();

  T max = gf->max(a->degree(), b->degree());

  for (int i = max; i >= 0; i--)
    result->set(i,
                gf->add(a->get(i), b->get(i)));
}

template <typename T>
void Poly<T>::_sub(Poly<T> *result, Poly<T> *a, Poly<T> *b)
{
  result->clear();

  T max = gf->max(a->degree(), b->degree());

  for (int i = max; i >= 0; i--)
    result->set(i,
                gf->sub(a->get(i), b->get(i)));
}

template <typename T>
void Poly<T>::_mul(Poly<T> *result, Poly<T> *a, Poly<T> *b)
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
void Poly<T>::_div(Poly<T> *q, Poly<T> *r, Poly<T> *n, Poly<T> *d)
{
  Poly<T> _q(gf), _r(gf);

  if (d->is_zero())
    throw NTL_EX_DIV_BY_ZERO;

  _q.clear();
  _r.copy(n);
  while (!_r.is_zero() && (_r.degree() >= d->degree())) {
    Poly<T> _t(gf);
    _t.set(_r.degree() - d->degree(),
           gf->div(_r.lead(), d->lead()));
    _q.add(&_t);
    _t.mul(d);
    _r.sub(&_t);
  }

  if (q != nullptr)
    q->copy(&_q);

  if (r != nullptr)
    r->copy(&_r);
}

template <typename T>
void Poly<T>::_derivative(Poly<T> *result, Poly<T> *a)
{
  result->clear();

  for (int i = a->degree(); i > 0; i--)
    result->set(i - 1, gf->mul(a->get(i), i % gf->p));
}

template <typename T>
void Poly<T>::neg()
{
  Poly<T> a(gf), b(gf);
  b.copy(this);
  _sub(this, &a, &b);
}

template <typename T>
void Poly<T>::add(Poly<T> *b)
{
  Poly<T> a(gf);
  a.copy(this);
  _add(this, &a, b);
}

template <typename T>
void Poly<T>::sub(Poly<T> *b)
{
  Poly<T> a(gf);
  a.copy(this);
  _sub(this, &a, b);
}

template <typename T>
void Poly<T>::mul(Poly<T> *b)
{
  Poly<T> a(gf);
  a.copy(this);
  _mul(this, &a, b);
}

template <typename T>
void Poly<T>::div(Poly<T> *b)
{
  Poly<T> a(gf);
  a.copy(this);
  _div(this, NULL, &a, b);
}

template <typename T>
void Poly<T>::mod(Poly<T> *b)
{
  Poly<T> a(gf);
  a.copy(this);
  _div(NULL, this, &a, b);
}

template <typename T>
void Poly<T>::derivative()
{
  Poly<T> a(gf);
  a.copy(this);
  _derivative(this, &a);
}

template <typename T>
T Poly<T>::eval(T x)
{
  int i = degree();
  T result = get(i);

  while (i >= 1) {
    result = gf->add(gf->mul(result, x), get(i - 1));
    i--;
  }
  return result;
}

template <typename T>
void Poly<T>::dump()
{
  typename Term::const_reverse_iterator it = terms.rbegin();

  for (; it != terms.rend(); it++)
    std::cout << " " << it->second << "x^" << it->first;
  std::cout << "\n";
}
