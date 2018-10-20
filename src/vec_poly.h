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
#ifndef __QUAD_VEC_POLY_H__
#define __QUAD_VEC_POLY_H__

#include "gf_base.h"
#include "gf_nf4.h"
#include "vec_vector.h"

namespace quadiron {
namespace vec {

/** A vector of size \f$n\f$ represents a polynomial \f$P(X)\f$ of degree
 * \f$(n-1)\f$
 *
 * Vector: \f$v[0], v[1], ..., v[n-1]\f$ represents the polynomial:
 *  \f$P(X):= v[0] + v[1] X + ... + v[n-1] X^n\f$
 */
template <typename T>
class Poly : public Vector<T> {
  public:
    Poly(const gf::Field<T>& field, int n);
    Poly(const Poly<T>& a);
    int get_deg() const;
    void set_deg(int exponent);
    const T& get(int exponent) const override;
    void set(int exponent, T coef) override;
    void derivative();
    void derivative_nf4();
    T eval(T x);
    void mul(Poly<T>* b, int deg_out);
    void mul_to_x_plus_coef(T coef);
    void neg() override;
    void zero_fill() override;
    void dump() const override;

  private:
    const gf::Field<T>* field;
    T field_characteristic;
    T* buf;
    int degree;
};

template <typename T>
Poly<T>::Poly(const gf::Field<T>& field, int n) : Vector<T>(field, n)
{
    this->field = &field;
    this->field_characteristic = field.get_p();
    this->buf = this->get_mem();
    this->degree = -1;
}

template <typename T>
Poly<T>::Poly(const Poly<T>& a) : Vector<T>(a.get_gf(), a.get_n())
{
    this->field = (gf::Field<T>*)this->rn;
    this->buf = this->get_mem();
    this->field_characteristic = field->get_p();

    std::copy_n(a.get_mem(), this->n, buf);
    this->degree = a.get_deg();
}

template <typename T>
int Poly<T>::get_deg() const
{
    return degree;
}

template <typename T>
void Poly<T>::set_deg(int deg)
{
    assert(deg >= 0 && deg < this->n);
    this->degree = deg;
}

template <typename T>
inline const T& Poly<T>::get(int exponent) const
{
    assert(exponent >= 0 && exponent < this->n);

    return buf[exponent];
}

template <typename T>
void Poly<T>::set(int exponent, T coef)
{
    assert(exponent >= 0 && exponent < this->n);
    assert(field->check(coef));

    buf[exponent] = coef;
    if (exponent > degree && coef > 0)
        degree = exponent;
}

/** Compute the derivative of a polynomial
 *
 * \f{eqnarray*}{
 * P(X) &= a_0 + a_1 X + a_2 X^2 + ... + a_n X^n \\
 * P'(X)&=       a_1   + a_2 X   + ... + a_n X^{n-1}
 * \f}
 *
 * Note, derivate of \f$X^n\f$ is defined as sum of \f$n\f$ polynomials
 * \f$X^{n-1}\f$, not a multiplication of \f$n\f$ and \f$X^{n-1}\f$. It can be
 * expressed as below:
 * \f{eqnarray*}{
 *  derivative(X^n) &= (1 + 1 + ... + 1) X^{n-1} \\
 *                  &= (n \% p) * X^{n-1}
 * \f}
 * where \f$p\f$ is the characteristic of the field.
 */
template <typename T>
void Poly<T>::derivative()
{
    for (int deg = 1; deg <= degree; ++deg) {
        buf[deg - 1] = field->mul(deg % field_characteristic, buf[deg]);
    }
    buf[degree] = 0;
    degree--;
}

/** Derivative particularly for NF4 */
template <typename T>
void Poly<T>::derivative_nf4()
{
    for (int deg = 1; deg <= degree; ++deg) {
        const T _deg = field->replicate(deg % field_characteristic);
        buf[deg - 1] = field->mul(_deg, buf[deg]);
    }
    buf[degree] = 0;
    degree--;
}

template <typename T>
T Poly<T>::eval(T x)
{
    int i = degree;
    T result = buf[i];

    while (i >= 1) {
        result = field->add(field->mul(result, x), buf[i - 1]);
        i--;
    }
    return result;
}

/** Multiply to a polynomial \f$b\f$ with a limited degree of result
 *
 * Degree of result polynomial is limited by a given value, i.e. elements of
 * degree greater than this value are ignored. It reduces operations in cases
 * the degree limit is smaller than degree of input polynomials
 *
 * @param b - polynomial to multiply
 * @param deg_out - degree of result
 */
template <typename T>
void Poly<T>::mul(Poly<T>* b, int deg_out)
{
    assert(deg_out >= 0 && deg_out < this->n);

    for (int deg = deg_out; deg >= 0; deg--) {
        T val = 0;
        for (int deg_a = 0; deg_a <= deg && deg_a <= degree; deg_a++) {
            const T val_a = buf[deg_a];
            if (val_a > 0) {
                const int deg_b = deg - deg_a;
                const T val_b = b->get(deg_b);
                if (val_b > 0)
                    val = field->add(val, field->mul(val_a, val_b));
            }
        }
        buf[deg] = val;
    }
}

/** Multiply to polynomial \f$(X + coef)\f$
 *
 * \f{eqnarray*}{
 * P(X) = &a_0 + a_1 X + ... + a_n X^n \\
 * Q(X) = &coef + X \\
 * R(X) = &P(X) \times Q(X) \\
 *      = &a_0 * coef + a_1 * coef X + ... + a_n X^n + \\
 *        &a_0 X + a_1 X^2 + ... + a_{n-1} X^n + a_n X^{n+1} \\
 *      = &a_0 * coef + sum_{i=1}^{n}(a_{i-1} + a_i * coef) X^i + a_n X^{n+1}
 * \f}
 *
 * @param coef - coef of the multiplied polynomial (X + coef)
 */
template <typename T>
void Poly<T>::mul_to_x_plus_coef(T coef)
{
    assert(field->check(coef));

    const T val = buf[degree];
    for (int deg = degree; deg > 0; deg--) {
        buf[deg] = field->add(buf[deg - 1], field->mul(buf[deg], coef));
    }
    buf[0] = field->mul(buf[0], coef);
    set(degree + 1, val);
}

template <typename T>
void Poly<T>::neg()
{
    for (int i = 0; i <= degree; i++)
        buf[i] = field->neg(buf[i]);
}

template <typename T>
void Poly<T>::zero_fill()
{
    std::memset(buf, 0, this->n * sizeof(*buf));
}

template <typename T>
void Poly<T>::dump() const
{
    for (int i = 0; i <= degree; i++) {
        if (i) {
            std::cout << " + ";
        }
        std::cout << buf[i] << " x^" << i;
    }
    std::cout << "\n";
}

} // namespace vec
} // namespace quadiron

#endif
