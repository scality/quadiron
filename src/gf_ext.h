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
#ifndef __QUAD_GF_EXT_H__
#define __QUAD_GF_EXT_H__

#include "arith.h"
#include "gf_base.h"

namespace quadiron {

template <typename T>
class Polynomial;

namespace gf {

/** An extension Galois Field GF(q<sup>n</sup>) where `q` is a prime number. */
template <typename T>
class Extension : public gf::Field<T> {
  public:
    explicit Extension(T p, int n);
    ~Extension() = default;
    void init();
    void find_primitive_root();
    T card(void) const override;
    T card_minus_one(void) const override;
    T neg(T a) const override;
    T add(T a, T b) const override;
    T sub(T a, T b) const override;
    T mul(T a, T b) const override;
    T div(T a, T b) const override;
    T inv(T a) const override;
    T exp(T a, T b) const override;
    T log(T a, T b) const override;
    T weak_rand(void) const override;

  private:
    int max_deg;
    std::unique_ptr<Polynomial<T>> primitive_poly;
};

/**
 * @note set the 3rd argument as false to deactivate the computing of primitive
 * root on creating the base class
 */
template <typename T>
Extension<T>::Extension(T p, int n) : gf::Field<T>(p, n, false)
{
    // XXX shall check that p is prime
    max_deg = 31; // XXX

    // XXX temp
    if (p == 2) {
        if (n == 8) {
            primitive_poly = std::unique_ptr<Polynomial<T>>(
                new Polynomial<T>(*(this->sub_field)));
            primitive_poly->set(8, 1);
            primitive_poly->set(4, 1);
            primitive_poly->set(3, 1);
            primitive_poly->set(1, 1);
            primitive_poly->set(0, 1);
        }
    } else if (p == 3) {
        if (n == 2) {
            // take irreducible polynomial from GF(3) of degree 2 P(X)=X^2+1
            primitive_poly = std::unique_ptr<Polynomial<T>>(
                new Polynomial<T>(*(this->sub_field)));
            primitive_poly->set(2, 1);
            primitive_poly->set(0, 1);
        } else if (n == 3) {
            // take irreducible polynomial from GF(3) of degree 3
            // P(X)=X^3+2X^2+1  another irreducible polynomial from GF(3) of
            // degree 3 is P(X)=X^3+2X+1
            primitive_poly = std::unique_ptr<Polynomial<T>>(
                new Polynomial<T>(*(this->sub_field)));
            primitive_poly->set(3, 1);
            primitive_poly->set(2, 2);
            // primitive_poly->set(1, 2);
            primitive_poly->set(0, 1);
        }
    } else {
        // XXX generate irreducible polynomials
        assert(false);
    }

    // computing primitive root
    this->init();
}

template <typename T>
void Extension<T>::init()
{
    // compute factors of order
    RingModN<T>::compute_factors_of_order();

    // compute root of unity
    this->find_primitive_root();
}

template <typename T>
void Extension<T>::find_primitive_root()
{
    if (this->root) {
        return;
    }

    T h = Extension<T>::card_minus_one();
    if (h == 1) {
        this->root = 1;
        return;
    }

    T nb = 2;
    bool ok;
    typename std::vector<T>::size_type i;

    while (nb <= h) {
        ok = true;
        // check nb^divisor == 1
        // std::cout << "checking.." << nb << std::endl;
        for (i = 0; i != this->proper_divisors->size(); ++i) {
            // std::cout << nb << "^" << this->proper_divisors->at(i) << "\n";
            // std::cout << this->exp(nb, this->proper_divisors->at(i)) <<
            // std::endl;
            if (Extension<T>::exp(nb, this->proper_divisors->at(i)) == 1) {
                ok = false;
                break;
            }
        }
        if (ok) {
            this->root = nb;
            break;
        }
        nb++;
    }

    if (!this->root) {
        assert(false); // Root not found.
    }
}

template <typename T>
T Extension<T>::card(void) const
{
    return arith::exp<T>(this->p, this->n);
}

template <typename T>
T Extension<T>::card_minus_one(void) const
{
    return Extension<T>::card() - 1;
}

template <typename T>
T Extension<T>::neg(T a) const
{
    assert(this->check(a));

    return sub(0, a);
}

template <typename T>
T Extension<T>::add(T a, T b) const
{
    assert(this->check(a));
    assert(this->check(b));

    Polynomial<T> _a(*(this->sub_field));
    _a.from_num(a, max_deg);
    Polynomial<T> _b(*(this->sub_field));
    _b.from_num(b, max_deg);

    _a.add(&_b);

    return _a.to_num();
}

template <typename T>
T Extension<T>::sub(T a, T b) const
{
    assert(this->check(a));
    assert(this->check(b));

    Polynomial<T> _a(*(this->sub_field));
    _a.from_num(a, max_deg);
    Polynomial<T> _b(*(this->sub_field));
    _b.from_num(b, max_deg);

    _a.sub(&_b);

    return _a.to_num();
}

template <typename T>
T Extension<T>::mul(T a, T b) const
{
    assert(this->check(a));
    assert(this->check(b));

    Polynomial<T> _a(*(this->sub_field));
    _a.from_num(a, max_deg);
    Polynomial<T> _b(*(this->sub_field));
    _b.from_num(b, max_deg);

    _a.mul(&_b);
    _a.mod(primitive_poly.get());

    return _a.to_num();
}

template <typename T>
T Extension<T>::div(T a, T b) const
{
    assert(this->check(a));
    assert(this->check(b));

    Polynomial<T> _a(*(this->sub_field));
    _a.from_num(a, max_deg);
    Polynomial<T> _b(*(this->sub_field));
    _b.from_num(b, max_deg);

    _a.div(&_b);
    _a.mod(primitive_poly.get());

    return _a.to_num();
}

template <typename T>
T Extension<T>::inv(T a) const
{
    assert(this->check(a));

    Polynomial<T> _u(*(this->sub_field));
    _u.copy(primitive_poly.get());
    Polynomial<T> _v(*(this->sub_field));
    _v.from_num(a, max_deg);

    Polynomial<T> _b1(*(this->sub_field));
    Polynomial<T> _b2(*(this->sub_field));
    Polynomial<T>* bezout[2];
    bezout[0] = &_b1;
    bezout[1] = &_b2;
    Polynomial<T> _gcd(*(this->sub_field));

    _gcd._extended_gcd(&_u, &_v, bezout, nullptr, &_gcd);

    assert(_gcd.degree() == 0);
    _gcd.set(0, this->sub_field->inv(_gcd.get(0)));

    _b2.mod(primitive_poly.get());
    _b2.mul(&_gcd);

    return _b2.to_num();
}

template <typename T>
T Extension<T>::exp(T a, T b) const
{
    assert(Extension<T>::check(a));
    assert(Extension<T>::check(b));

    return gf::Field<T>::exp_quick(a, b);
}

template <typename T>
T Extension<T>::log(T a, T b) const
{
    assert(this->check(a));

    return gf::Field<T>::log_naive(a, b);
}

template <typename T>
T Extension<T>::weak_rand(void) const
{
    std::uniform_int_distribution<uint32_t> dis(0, this->p - 1);
    Polynomial<T> _a(*(this->sub_field));
    T num;

retry:
    for (int i = 0; i < this->n; i++)
        _a.set(i, dis(prng()));

    num = _a.to_num();
    if (0 == num)
        goto retry;

    return num;
}

} // namespace gf
} // namespace quadiron

#endif
