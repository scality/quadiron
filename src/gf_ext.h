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
#ifndef __NTTEC_GF_EXT_H__
#define __NTTEC_GF_EXT_H__

#include "arith.h"
#include "gf_base.h"

namespace nttec {

template <typename T>
class Polynomial;

namespace gf {

/** An extension Galois Field GF(q<sup>n</sup>) where `q` is a prime number. */
template <typename T>
class Extension : public gf::Field<T> {
  public:
    explicit Extension(T p, int n);
    ~Extension();
    T card(void);
    T card_minus_one(void);
    T neg(T a);
    T add(T a, T b);
    T sub(T a, T b);
    T mul(T a, T b);
    T div(T a, T b);
    T inv(T a);
    T exp(T a, T b);
    T log(T a, T b);
    T weak_rand(void);

  private:
    gf::Field<T>* sub_field;
    int max_deg;
    Polynomial<T>* primitive_poly;
};

template <typename T>
Extension<T>::Extension(T p, int n) : gf::Field<T>(p, n)
{
    // XXX shall check that p is prime
    sub_field = new Prime<T>(p);
    max_deg = 31; // XXX

    // XXX temp
    if (p == 2) {
        if (n == 8) {
            primitive_poly = new Polynomial<T>(sub_field);
            primitive_poly->set(8, 1);
            primitive_poly->set(4, 1);
            primitive_poly->set(3, 1);
            primitive_poly->set(1, 1);
            primitive_poly->set(0, 1);
        }
    } else if (p == 3) {
        if (n == 2) {
            // take irreducible polynomial from GF(3) of degree 2 P(X)=X^2+1
            primitive_poly = new Polynomial<T>(sub_field);
            primitive_poly->set(2, 1);
            primitive_poly->set(0, 1);
        } else if (n == 3) {
            // take irreducible polynomial from GF(3) of degree 3
            // P(X)=X^3+2X^2+1  another irreducible polynomial from GF(3) of
            // degree 3 is P(X)=X^3+2X+1
            primitive_poly = new Polynomial<T>(sub_field);
            primitive_poly->set(3, 1);
            primitive_poly->set(2, 2);
            // primitive_poly->set(1, 2);
            primitive_poly->set(0, 1);
        }
    } else {
        // XXX generate irreducible polynomials
        assert(false);
    }
}

template <typename T>
Extension<T>::~Extension()
{
    delete primitive_poly;
    delete sub_field;
}

template <typename T>
T Extension<T>::card(void)
{
    return arith::exp<T>(this->p, this->n);
}

template <typename T>
T Extension<T>::card_minus_one(void)
{
    return card() - 1;
}

template <typename T>
T Extension<T>::neg(T a)
{
    assert(this->check(a));

    return sub(0, a);
}

template <typename T>
T Extension<T>::add(T a, T b)
{
    assert(this->check(a));
    assert(this->check(b));

    Polynomial<T> _a(sub_field);
    _a.from_num(a, max_deg);
    Polynomial<T> _b(sub_field);
    _b.from_num(b, max_deg);

    _a.add(&_b);

    return _a.to_num();
}

template <typename T>
T Extension<T>::sub(T a, T b)
{
    assert(this->check(a));
    assert(this->check(b));

    Polynomial<T> _a(sub_field);
    _a.from_num(a, max_deg);
    Polynomial<T> _b(sub_field);
    _b.from_num(b, max_deg);

    _a.sub(&_b);

    return _a.to_num();
}

template <typename T>
T Extension<T>::mul(T a, T b)
{
    assert(this->check(a));
    assert(this->check(b));

    Polynomial<T> _a(sub_field);
    _a.from_num(a, max_deg);
    Polynomial<T> _b(sub_field);
    _b.from_num(b, max_deg);

    _a.mul(&_b);
    _a.mod(primitive_poly);

    return _a.to_num();
}

template <typename T>
T Extension<T>::div(T a, T b)
{
    assert(this->check(a));
    assert(this->check(b));

    Polynomial<T> _a(sub_field);
    _a.from_num(a, max_deg);
    Polynomial<T> _b(sub_field);
    _b.from_num(b, max_deg);

    _a.div(&_b);
    _a.mod(primitive_poly);

    return _a.to_num();
}

template <typename T>
T Extension<T>::inv(T a)
{
    assert(this->check(a));

    Polynomial<T> _u(sub_field);
    _u.copy(primitive_poly);
    Polynomial<T> _v(sub_field);
    _v.from_num(a, max_deg);

    Polynomial<T> _b1(sub_field);
    Polynomial<T> _b2(sub_field);
    Polynomial<T>* bezout[2];
    bezout[0] = &_b1;
    bezout[1] = &_b2;
    Polynomial<T> _gcd(sub_field);

    _gcd._extended_gcd(&_u, &_v, bezout, nullptr, &_gcd);

    assert(_gcd.degree() == 0);
    _gcd.set(0, sub_field->inv(_gcd.get(0)));

    _b2.mod(primitive_poly);
    _b2.mul(&_gcd);

    return _b2.to_num();
}

template <typename T>
T Extension<T>::exp(T a, T b)
{
    assert(this->check(a));
    assert(this->check(b));

    return gf::Field<T>::exp_quick(a, b);
}

template <typename T>
T Extension<T>::log(T a, T b)
{
    assert(this->check(a));

    return gf::Field<T>::log_naive(a, b);
}

template <typename T>
T Extension<T>::weak_rand(void)
{
    Polynomial<T> _a(sub_field);
    T num;

retry:
    for (int i = 0; i < this->n; i++)
        _a.set(i, arith::weak_rand0<T>(this->p));

    num = _a.to_num();
    if (0 == num)
        goto retry;

    return num;
}

} // namespace gf
} // namespace nttec

#endif
