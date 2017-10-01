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
  void copy(Poly<T> *src, T offset);
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
  void compute_g0_g1(T deg0, T deg2, Poly<T> *g0, Poly<T> *g1);
  void taylor_expand(std::vector<Poly<T>> *result, T n, T t);
  void taylor_expand_t2(Poly<T> *g0, Poly<T> *g1, T n, T deg=0);
  T eval(T x);
  bool equal(Poly<T> *f);
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
void Poly<T>::copy(Poly<T> *src, T offset)
{
  clear();

  for (int i = src->degree(); i >= 0; i--)
    set(i + offset, src->get(i));
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

/**
 * Taylor expansion at (x^t - x)
 *  Algorithm 1 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 * @param result: vector of hi(x) polynomials
 * @param n
 * @param t
 */
template <typename T>
void Poly<T>::taylor_expand(std::vector<Poly<T>> *result, T n, T t)
{
  // it supports only GF2N
  assert(gf->p == 2);
  assert(n >= 1);
  assert(t > 1);

  T deg = degree();
  assert(deg < n);

  if (n <= t) {
    Poly<T> tmp(gf);
    tmp.copy(this);
    result->push_back(tmp);
    return;
  }

  T k = gf->arith->log2(n/t - 1);
  T deg2 = gf->arith->exp2(k);
  T deg0 = t*deg2;
  Poly<T> g0(gf), g1(gf);
  compute_g0_g1(deg0, deg2, &g0, &g1);

  std::vector<Poly<T>> V1, V2;
  g0.taylor_expand(&V1, deg0, t);
  g1.taylor_expand(&V2, n - deg0, t);

  result->reserve(V1.size() + V2.size() ); // preallocate memory
  result->insert(result->end(), V1.begin(), V1.end() );
  result->insert(result->end(), V2.begin(), V2.end());
}

template <typename T>
void Poly<T>::compute_g0_g1(T deg0, T deg2, Poly<T> *g0, Poly<T> *g1) {
  // find k s.t. t*2^k < n <= 2*t*2^k
  T deg1 = deg0 - deg2;

  // find f0, f1, f2 s.t.
  // f = f0 + x^(t2^k)( f1 + x^( (t-1)2^k) f2 )
  Poly<T> f0(gf), f1(gf), f2(gf);
  T i, j;

  for (i = 0; i < deg0; i++)
    f0.set(i, get(i));
  for (j = 0; j < deg1; j++, i++)
    f1.set(j, get(i));
  for (j = 0; j < deg2; j++, i++)
    f2.set(j, get(i));

  // std::cout << "f0:"; f0.dump();
  // std::cout << "f1:"; f1.dump();
  // std::cout << "f2:"; f2.dump();

  Poly<T> h(gf);
  h.copy(&f1);
  h.add(&f2);
  g0->copy(&h, deg2);
  g0->add(&f0);
  g1->copy(&f2, deg1);
  g1->add(&h);
}


/**
 * Taylor expansion at (x^2 - x)
 *  Algorithm 1 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 * @param G0: poly of gi0 of hi(x) polynomials = gi0 + gi1 * x
 * @param G1: poly of gi1 of hi(x) polynomials = gi0 + gi1 * x
 * @param n
 * @param s_deg: start degree of G0 and G1
 */
template <typename T>
void Poly<T>::taylor_expand_t2(Poly<T> *G0, Poly<T> *G1, T n, T s_deg)
{
  // it supports only GF2N
  assert(gf->p == 2);
  assert(n >= 1);

  T deg = degree();
  assert(deg < n);
  if (n <= 2) {
    G0->set(s_deg, this->get(0));
    G1->set(s_deg, this->get(1));
    return;
  }

  T k = gf->arith->log2(n/2 - 1);
  T deg2 = gf->arith->exp2(k);
  T deg0 = 2*deg2;

  Poly<T> g0(gf), g1(gf);
  compute_g0_g1(deg0, deg2, &g0, &g1);

  g0.taylor_expand_t2(G0, G1, deg0, s_deg);
  g1.taylor_expand_t2(G0, G1, n - deg0, s_deg + deg2);
}

template <typename T>
bool Poly<T>::equal(Poly<T> *f) {
  T deg = degree();
  if (deg != f->degree())
    return false;
  for (T i = 0; i <= deg; i++)
    if (get(i) != f->get(i))
      return false;
  return true;
}

template <typename T>
void Poly<T>::dump()
{
  typename Term::const_reverse_iterator it = terms.rbegin();

  for (; it != terms.rend(); it++)
    std::cout << " " << it->second << "x^" << it->first;
  std::cout << "\n";
}
