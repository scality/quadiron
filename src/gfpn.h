/* -*- mode: c++ -*- */
#ifndef __NTTEC_GFPN_H__
#define __NTTEC_GFPN_H__

#include "gf.h"

namespace nttec {

template <typename T>
class Poly;

namespace gf {

template <typename T>
class GFPN : public GF<T> {
  private:
    GF<T>* sub_field;
    int max_deg;
    Poly<T>* prim_poly;

  public:
    explicit GFPN(T p, int n);
    ~GFPN();
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
};

template <typename T>
GFPN<T>::GFPN(T p, int n) : GF<T>(p, n)
{
    // XXX shall check that p is prime
    sub_field = new GFP<T>(p);
    max_deg = 31; // XXX

    // XXX temp
    if (p == 2) {
        if (n == 8) {
            prim_poly = new Poly<T>(sub_field);
            prim_poly->set(8, 1);
            prim_poly->set(4, 1);
            prim_poly->set(3, 1);
            prim_poly->set(1, 1);
            prim_poly->set(0, 1);
        }
    } else if (p == 3) {
        if (n == 2) {
            // take irreducible polynomial from GF(3) of degree 2 P(X)=X^2+1
            prim_poly = new Poly<T>(sub_field);
            prim_poly->set(2, 1);
            prim_poly->set(0, 1);
        } else if (n == 3) {
            // take irreducible polynomial from GF(3) of degree 3
            // P(X)=X^3+2X^2+1  another irreducible polynomial from GF(3) of
            // degree 3 is P(X)=X^3+2X+1
            prim_poly = new Poly<T>(sub_field);
            prim_poly->set(3, 1);
            prim_poly->set(2, 2);
            // prim_poly->set(1, 2);
            prim_poly->set(0, 1);
        }
    } else {
        // XXX generate irreducible polynomials
        assert(false);
    }
}

template <typename T>
GFPN<T>::~GFPN()
{
    delete prim_poly;
    delete sub_field;
}

template <typename T>
T GFPN<T>::card(void)
{
    return arith::_exp<T>(this->p, this->n);
}

template <typename T>
T GFPN<T>::card_minus_one(void)
{
    return card() - 1;
}

template <typename T>
T GFPN<T>::neg(T a)
{
    assert(this->check(a));

    return sub(0, a);
}

template <typename T>
T GFPN<T>::add(T a, T b)
{
    assert(this->check(a));
    assert(this->check(b));

    Poly<T> _a(sub_field);
    _a.from_num(a, max_deg);
    Poly<T> _b(sub_field);
    _b.from_num(b, max_deg);

    _a.add(&_b);

    return _a.to_num();
}

template <typename T>
T GFPN<T>::sub(T a, T b)
{
    assert(this->check(a));
    assert(this->check(b));

    Poly<T> _a(sub_field);
    _a.from_num(a, max_deg);
    Poly<T> _b(sub_field);
    _b.from_num(b, max_deg);

    _a.sub(&_b);

    return _a.to_num();
}

template <typename T>
T GFPN<T>::mul(T a, T b)
{
    assert(this->check(a));
    assert(this->check(b));

    Poly<T> _a(sub_field);
    _a.from_num(a, max_deg);
    Poly<T> _b(sub_field);
    _b.from_num(b, max_deg);

    _a.mul(&_b);
    _a.mod(prim_poly);

    return _a.to_num();
}

template <typename T>
T GFPN<T>::div(T a, T b)
{
    assert(this->check(a));
    assert(this->check(b));

    Poly<T> _a(sub_field);
    _a.from_num(a, max_deg);
    Poly<T> _b(sub_field);
    _b.from_num(b, max_deg);

    _a.div(&_b);
    _a.mod(prim_poly);

    return _a.to_num();
}

template <typename T>
T GFPN<T>::inv(T a)
{
    assert(this->check(a));

    Poly<T> _u(sub_field);
    _u.copy(prim_poly);
    Poly<T> _v(sub_field);
    _v.from_num(a, max_deg);

    Poly<T> _b1(sub_field);
    Poly<T> _b2(sub_field);
    Poly<T>* bezout[2];
    bezout[0] = &_b1;
    bezout[1] = &_b2;
    Poly<T> _gcd(sub_field);

    _gcd._extended_gcd(&_u, &_v, bezout, nullptr, &_gcd);

    assert(_gcd.degree() == 0);
    _gcd.set(0, sub_field->inv(_gcd.get(0)));

    _b2.mod(prim_poly);
    _b2.mul(&_gcd);

    return _b2.to_num();
}

template <typename T>
T GFPN<T>::exp(T a, T b)
{
    assert(this->check(a));
    assert(this->check(b));

    return GF<T>::exp_quick(a, b);
}

template <typename T>
T GFPN<T>::log(T a, T b)
{
    assert(this->check(a));

    return GF<T>::log_naive(a, b);
}

template <typename T>
T GFPN<T>::weak_rand(void)
{
    Poly<T> _a(sub_field);
    T num;

retry:
    for (int i = 0; i < this->n; i++)
        _a.set(i, arith::_weak_rand0<T>(this->p));

    num = _a.to_num();
    if (0 == num)
        goto retry;

    return num;
}

} // namespace gf
} // namespace nttec

#endif
