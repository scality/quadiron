/* -*- mode: c++ -*- */
#ifndef __NTL_POLY_H__
#define __NTL_POLY_H__

#include <map>

#include "core.h"
#include "gf.h"

template <typename T>
class Vec;

template <typename T>
class Poly {
  private:
    struct Term : std::map<int, T> {
    };
    GF<T>* gf;
    GF<T>* sub_field;
    Term terms;

  public:
    explicit Poly(GF<T>* gf);
    void clear();
    void copy(Poly<T>* src);
    void copy(Poly<T>* src, T offset);
    int degree();
    T lead();
    bool is_zero();
    T get(int exponent);
    void set(int exponent, T coef);
    void _neg(Poly<T>* result, Poly<T>* a);
    void _add(Poly<T>* result, Poly<T>* a, Poly<T>* b);
    void _sub(Poly<T>* result, Poly<T>* a, Poly<T>* b);
    void _mul(Poly<T>* result, Poly<T>* a, Poly<T>* b);
    void _div(Poly<T>* q, Poly<T>* r, Poly<T>* n, Poly<T>* d);
    void _gcd(Poly<T>* u, Poly<T>* v, Poly<T>* gcd);
    void _extended_gcd(
        Poly<T>* a,
        Poly<T>* b,
        Poly<T>* bezout_coef[2],
        Poly<T>* quotient_gcd[2],
        Poly<T>* gcd);
    void _derivative(Poly<T>* result, Poly<T>* a);
    void neg();
    void add(Poly<T>* b);
    void sub(Poly<T>* b);
    void mul(Poly<T>* b);
    void div(Poly<T>* d);
    void mod(Poly<T>* d);
    void derivative();
    T eval(T x);
    bool equal(Poly<T>* f);
    void to_vec(Vec<T>* vec);
    T to_num();
    void from_num(T x, int max_deg);
    void dump(std::ostream& dest);
    void dump();
};

template <typename T>
Poly<T>::Poly(GF<T>* gf)
{
    this->gf = gf;
    this->sub_field = gf->get_sub_field();
}

template <typename T>
void Poly<T>::clear()
{
    terms.clear();
}

template <typename T>
void Poly<T>::copy(Poly<T>* src)
{
    clear();

    for (int i = src->degree(); i >= 0; i--)
        set(i, src->get(i));
}

template <typename T>
void Poly<T>::copy(Poly<T>* src, T offset)
{
    clear();

    for (int i = src->degree(); i >= 0; i--)
        set(i + offset, src->get(i));
}

template <typename T>
int Poly<T>::degree()
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
T Poly<T>::get(int exponent)
{
    assert(exponent >= 0);

    typename Term::const_iterator it = terms.find(exponent);

    return (it == terms.end()) ? 0 : it->second;
}

template <typename T>
void Poly<T>::set(int exponent, T coef)
{
    assert(exponent >= 0);
    assert(gf->check(coef));

    typename Term::const_iterator it = terms.find(exponent);

    if (it == terms.end()) {
        if (coef == 0)
            return;
    } else {
        if (coef == 0) {
            terms.erase(it);
            return;
        }
    }
    terms[exponent] = coef;
}

template <typename T>
void Poly<T>::_neg(Poly<T>* result, Poly<T>* a)
{
    result->clear();

    Poly<T> b(gf);
    sub(result, &b, a);
}

template <typename T>
void Poly<T>::_add(Poly<T>* result, Poly<T>* a, Poly<T>* b)
{
    result->clear();

    int max = std::max(a->degree(), b->degree());

    for (int i = max; i >= 0; i--)
        result->set(i, gf->add(a->get(i), b->get(i)));
}

template <typename T>
void Poly<T>::_sub(Poly<T>* result, Poly<T>* a, Poly<T>* b)
{
    result->clear();

    int max = std::max(a->degree(), b->degree());

    for (int i = max; i >= 0; i--)
        result->set(i, gf->sub(a->get(i), b->get(i)));
}

template <typename T>
void Poly<T>::_mul(Poly<T>* result, Poly<T>* a, Poly<T>* b)
{
    result->clear();

    for (int i = a->degree(); i >= 0; i--)
        for (int j = b->degree(); j >= 0; j--)
            result->set(
                i + j,
                gf->add(result->get(i + j), gf->mul(a->get(i), b->get(j))));
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
void Poly<T>::_div(Poly<T>* q, Poly<T>* r, Poly<T>* n, Poly<T>* d)
{
    Poly<T> _q(gf), _r(gf);

    if (d->is_zero())
        throw NTL_EX_DIV_BY_ZERO;

    _q.clear();
    _r.copy(n);

    while (!_r.is_zero() && (_r.degree() >= d->degree())) {
        Poly<T> _t(gf);
        _t.set(_r.degree() - d->degree(), gf->div(_r.lead(), d->lead()));
        _q.add(&_t);
        _t.mul(d);
        _r.sub(&_t);
    }

    if (q != nullptr)
        q->copy(&_q);

    if (r != nullptr)
        r->copy(&_r);
}

/**
 * compute the GCD of 2 different polynomials (modern algorithm)
 *
 * @param u
 * @param v
 * @param gcd cannot be nullptr
 */
template <typename T>
void Poly<T>::_gcd(Poly<T>* u, Poly<T>* v, Poly<T>* gcd)
{
    Poly<T> r(sub_field);

    if (v->degree() == 0) {
        gcd->copy(u);
        return;
    }
    if (u->degree() == 0) {
        gcd->copy(v);
        return;
    }
    while (1) {
        r.copy(u);
        r.mod(v);
        if (r.degree() == 0) {
            gcd->copy(&r);
            return;
        }
        u->copy(v);
        v->copy(&r);
    }
}

/**
 * compute the extended euclidean algorithm
 *
 * @param a (e.g. the polynomial that define the field)
 * @param b (e.g. the element of the extension field)
 * @param bezout_coef can be nullptr
 * @param quotient_gcd can be nullptr
 * @param gcd can be nullptr
 */
template <typename T>
void Poly<T>::_extended_gcd(
    Poly<T>* a,
    Poly<T>* b,
    Poly<T>* bezout_coef[2],
    Poly<T>* quotient_gcd[2],
    Poly<T>* gcd)
{
    Poly<T> s(sub_field);
    Poly<T> old_s(sub_field);
    old_s.set(0, 1);
    Poly<T> t(sub_field);
    t.set(0, 1);
    Poly<T> old_t(sub_field);
    Poly<T> r(sub_field);
    r.copy(b);
    Poly<T> old_r(sub_field);
    old_r.copy(a);
    Poly<T> quotient(sub_field);
    Poly<T> tmp(sub_field);
    Poly<T> tmp2(sub_field);

    while (!r.is_zero()) {
        _div(&quotient, nullptr, &old_r, &r);

        tmp.copy(&r);
        tmp2.copy(&quotient);
        tmp2.mul(&r);
        r.copy(&old_r);
        r.sub(&tmp2);
        old_r.copy(&tmp);

        tmp.copy(&s);
        tmp2.copy(&quotient);
        tmp2.mul(&s);
        s.copy(&old_s);
        s.sub(&tmp2);
        old_s.copy(&tmp);

        tmp.copy(&t);
        tmp2.copy(&quotient);
        tmp2.mul(&t);
        t.copy(&old_t);
        t.sub(&tmp2);
        old_t.copy(&tmp);
    }

    if (old_r.degree() > 0) {
        // the primitive polynomial is not irreducible or
        // the element is multiple of the primitive polynomial
        throw(NTL_EX_INVAL);
    }

    if (bezout_coef) {
        bezout_coef[0]->copy(&old_s);
        bezout_coef[1]->copy(&old_t);
    }

    if (quotient_gcd) {
        quotient_gcd[0]->copy(&t);
        quotient_gcd[1]->copy(&s);
    }

    if (gcd) {
        gcd->copy(&old_r);
    }
}

template <typename T>
void Poly<T>::_derivative(Poly<T>* result, Poly<T>* a)
{
    T _card;

    if (sub_field)
        _card = sub_field->card();
    else
        _card = gf->card();

    result->clear();

    for (int i = a->degree(); i > 0; i--)
        result->set(i - 1, gf->mul(a->get(i), i % _card));
}

template <typename T>
void Poly<T>::neg()
{
    Poly<T> a(gf), b(gf);
    b.copy(this);
    _sub(this, &a, &b);
}

template <typename T>
void Poly<T>::add(Poly<T>* b)
{
    Poly<T> a(gf);
    a.copy(this);
    _add(this, &a, b);
}

template <typename T>
void Poly<T>::sub(Poly<T>* b)
{
    Poly<T> a(gf);
    a.copy(this);
    _sub(this, &a, b);
}

template <typename T>
void Poly<T>::mul(Poly<T>* b)
{
    Poly<T> a(gf);
    a.copy(this);
    _mul(this, &a, b);
}

template <typename T>
void Poly<T>::div(Poly<T>* b)
{
    Poly<T> a(gf);
    a.copy(this);
    _div(this, nullptr, &a, b);
}

template <typename T>
void Poly<T>::mod(Poly<T>* b)
{
    Poly<T> a(gf);
    a.copy(this);
    _div(nullptr, this, &a, b);
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
bool Poly<T>::equal(Poly<T>* f)
{
    T deg = degree();
    if (deg != f->degree())
        return false;
    for (T i = 0; i <= deg; i++)
        if (get(i) != f->get(i))
            return false;
    return true;
}

template <typename T>
void Poly<T>::to_vec(Vec<T>* vec)
{
    T deg = degree();
    assert(vec->get_n() > deg);
    int i;
    for (i = 0; i <= deg; i++)
        vec->set(i, get(i));
    for (; i < vec->get_n(); i++)
        vec->set(i, 0);
}

/**
 * convert a polynomial to numerical representation
 *
 *
 * @return
 */
template <typename T>
T Poly<T>::to_num()
{
    int i = degree();
    T result = 0;

    while (i >= 0) {
        result += get(i) * _exp<T>(gf->card(), i);
        i--;
    }
    return result;
}

/**
 * convert a numerical representation of a polynomial
 *
 * @param x
 */
template <typename T>
void Poly<T>::from_num(T x, int max_deg)
{
    clear(); // XXX
    for (int i = max_deg; i >= 0; i--) {
        T tmp = _exp<T>(gf->card(), i);
        T tmp2 = x / tmp;
        if (tmp2 != 0)
            set(i, tmp2);
        else
            set(i, 0);
        x -= tmp2 * tmp;
    }
}

template <typename T>
void Poly<T>::dump(std::ostream& dest)
{
    typename Term::const_reverse_iterator it = terms.rbegin();

    if (it == terms.rend()) {
        dest << "0";
        return;
    }
    bool first = true;
    for (; it != terms.rend(); it++) {
        if (first) {
            first = false;
        } else {
            dest << " + ";
        }
        dest << it->second;
        if (0 != it->first)
            dest << "x^" << it->first;
    }
}

template <typename T>
void Poly<T>::dump()
{
    dump(std::cout);
    std::cout << "\n";
}

#endif
