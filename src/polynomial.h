/* -*- mode: c++ -*- */
/*
 * Copyright 2017-2018 Scality
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef __NTTEC_POLYNOMIAL_H__
#define __NTTEC_POLYNOMIAL_H__

#include <iostream>
#include <map>

#include "arith.h"
#include "core.h"
#include "gf_base.h"

namespace nttec {

namespace vec {

template <typename T>
class Vector;

} // namespace vec

template <typename T>
class Polynomial {
  public:
    explicit Polynomial(gf::Field<T>* field);
    void clear();
    void copy(Polynomial<T>* src);
    void copy(Polynomial<T>* src, T offset);
    int degree();
    T lead();
    bool is_zero();
    T get(int exponent);
    void set(int exponent, T coef);
    void _neg(Polynomial<T>* result, Polynomial<T>* a);
    void _add(Polynomial<T>* result, Polynomial<T>* a, Polynomial<T>* b);
    void _sub(Polynomial<T>* result, Polynomial<T>* a, Polynomial<T>* b);
    void _mul(Polynomial<T>* result, Polynomial<T>* a, Polynomial<T>* b);
    void _div(
        Polynomial<T>* q,
        Polynomial<T>* r,
        Polynomial<T>* n,
        Polynomial<T>* d);
    void _gcd(Polynomial<T>* u, Polynomial<T>* v, Polynomial<T>* gcd);
    void _extended_gcd(
        Polynomial<T>* a,
        Polynomial<T>* b,
        Polynomial<T>* bezout_coef[2],
        Polynomial<T>* quotient_gcd[2],
        Polynomial<T>* gcd);
    void _derivative(Polynomial<T>* result, Polynomial<T>* a);
    void neg();
    void add(Polynomial<T>* b);
    void sub(Polynomial<T>* b);
    void mul(Polynomial<T>* b);
    void div(Polynomial<T>* d);
    void mod(Polynomial<T>* d);
    void derivative();
    T eval(T x);
    bool equal(Polynomial<T>* f);
    void to_vec(vec::Vector<T>* vec);
    T to_num();
    void from_num(T x, int max_deg);
    void dump(std::ostream& dest);
    void dump();

  private:
    struct Term : std::map<int, T> {
    };
    gf::Field<T>* field;
    gf::Field<T>* sub_field;
    Term terms;
};

template <typename T>
Polynomial<T>::Polynomial(gf::Field<T>* field)
{
    this->field = field;
    this->sub_field = field->get_sub_field();
}

template <typename T>
void Polynomial<T>::clear()
{
    terms.clear();
}

template <typename T>
void Polynomial<T>::copy(Polynomial<T>* src)
{
    clear();

    for (int i = src->degree(); i >= 0; i--)
        set(i, src->get(i));
}

template <typename T>
void Polynomial<T>::copy(Polynomial<T>* src, T offset)
{
    clear();

    for (int i = src->degree(); i >= 0; i--)
        set(i + offset, src->get(i));
}

template <typename T>
int Polynomial<T>::degree()
{
    return (terms.rbegin() == terms.rend()) ? 0 : terms.rbegin()->first;
}

template <typename T>
T Polynomial<T>::lead()
{
    return get(degree());
}

template <typename T>
bool Polynomial<T>::is_zero()
{
    return degree() == 0 && get(0) == 0;
}

template <typename T>
T Polynomial<T>::get(int exponent)
{
    assert(exponent >= 0);

    typename Term::const_iterator it = terms.find(exponent);

    return (it == terms.end()) ? 0 : it->second;
}

template <typename T>
void Polynomial<T>::set(int exponent, T coef)
{
    assert(exponent >= 0);
    assert(field->check(coef));

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
void Polynomial<T>::_neg(Polynomial<T>* result, Polynomial<T>* a)
{
    result->clear();

    Polynomial<T> b(field);
    sub(result, &b, a);
}

template <typename T>
void Polynomial<T>::_add(
    Polynomial<T>* result,
    Polynomial<T>* a,
    Polynomial<T>* b)
{
    result->clear();

    int max = std::max(a->degree(), b->degree());

    for (int i = max; i >= 0; i--)
        result->set(i, field->add(a->get(i), b->get(i)));
}

template <typename T>
void Polynomial<T>::_sub(
    Polynomial<T>* result,
    Polynomial<T>* a,
    Polynomial<T>* b)
{
    result->clear();

    int max = std::max(a->degree(), b->degree());

    for (int i = max; i >= 0; i--)
        result->set(i, field->sub(a->get(i), b->get(i)));
}

template <typename T>
void Polynomial<T>::_mul(
    Polynomial<T>* result,
    Polynomial<T>* a,
    Polynomial<T>* b)
{
    result->clear();

    for (int i = a->degree(); i >= 0; i--)
        for (int j = b->degree(); j >= 0; j--)
            result->set(
                i + j,
                field->add(
                    result->get(i + j), field->mul(a->get(i), b->get(j))));
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
void Polynomial<T>::_div(
    Polynomial<T>* q,
    Polynomial<T>* r,
    Polynomial<T>* n,
    Polynomial<T>* d)
{
    Polynomial<T> _q(field), _r(field);

    if (d->is_zero()) {
        throw DomainError("divisor is zero");
    }
    _q.clear();
    _r.copy(n);

    while (!_r.is_zero() && (_r.degree() >= d->degree())) {
        Polynomial<T> _t(field);
        _t.set(_r.degree() - d->degree(), field->div(_r.lead(), d->lead()));
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
void Polynomial<T>::_gcd(Polynomial<T>* u, Polynomial<T>* v, Polynomial<T>* gcd)
{
    Polynomial<T> r(sub_field);

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
void Polynomial<T>::_extended_gcd(
    Polynomial<T>* a,
    Polynomial<T>* b,
    Polynomial<T>* bezout_coef[2],
    Polynomial<T>* quotient_gcd[2],
    Polynomial<T>* gcd)
{
    Polynomial<T> s(sub_field);
    Polynomial<T> old_s(sub_field);
    old_s.set(0, 1);
    Polynomial<T> t(sub_field);
    t.set(0, 1);
    Polynomial<T> old_t(sub_field);
    Polynomial<T> r(sub_field);
    r.copy(b);
    Polynomial<T> old_r(sub_field);
    old_r.copy(a);
    Polynomial<T> quotient(sub_field);
    Polynomial<T> tmp(sub_field);
    Polynomial<T> tmp2(sub_field);

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
        throw DomainError(
            "irreducible primitive polynomial or "
            "multiple of the primitive polynomial");
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
void Polynomial<T>::_derivative(Polynomial<T>* result, Polynomial<T>* a)
{
    T _card;

    if (sub_field)
        _card = sub_field->card();
    else
        _card = field->card();

    result->clear();

    for (int i = a->degree(); i > 0; i--)
        result->set(i - 1, field->mul(a->get(i), i % _card));
}

template <typename T>
void Polynomial<T>::neg()
{
    Polynomial<T> a(field), b(field);
    b.copy(this);
    _sub(this, &a, &b);
}

template <typename T>
void Polynomial<T>::add(Polynomial<T>* b)
{
    Polynomial<T> a(field);
    a.copy(this);
    _add(this, &a, b);
}

template <typename T>
void Polynomial<T>::sub(Polynomial<T>* b)
{
    Polynomial<T> a(field);
    a.copy(this);
    _sub(this, &a, b);
}

template <typename T>
void Polynomial<T>::mul(Polynomial<T>* b)
{
    Polynomial<T> a(field);
    a.copy(this);
    _mul(this, &a, b);
}

template <typename T>
void Polynomial<T>::div(Polynomial<T>* b)
{
    Polynomial<T> a(field);
    a.copy(this);
    _div(this, nullptr, &a, b);
}

template <typename T>
void Polynomial<T>::mod(Polynomial<T>* b)
{
    Polynomial<T> a(field);
    a.copy(this);
    _div(nullptr, this, &a, b);
}

template <typename T>
void Polynomial<T>::derivative()
{
    Polynomial<T> a(field);
    a.copy(this);
    _derivative(this, &a);
}

template <typename T>
T Polynomial<T>::eval(T x)
{
    int i = degree();
    T result = get(i);

    while (i >= 1) {
        result = field->add(field->mul(result, x), get(i - 1));
        i--;
    }
    return result;
}

template <typename T>
bool Polynomial<T>::equal(Polynomial<T>* f)
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
void Polynomial<T>::to_vec(vec::Vector<T>* vec)
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
T Polynomial<T>::to_num()
{
    int i = degree();
    T result = 0;

    while (i >= 0) {
        result += get(i) * arith::exp<T>(field->card(), i);
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
void Polynomial<T>::from_num(T x, int max_deg)
{
    clear(); // XXX
    for (int i = max_deg; i >= 0; i--) {
        T tmp = arith::exp<T>(field->card(), i);
        T tmp2 = x / tmp;
        if (tmp2 != 0)
            set(i, tmp2);
        else
            set(i, 0);
        x -= tmp2 * tmp;
    }
}

template <typename T>
void Polynomial<T>::dump(std::ostream& dest)
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
void Polynomial<T>::dump()
{
    dump(std::cout);
    std::cout << "\n";
}

} // namespace nttec

#endif
