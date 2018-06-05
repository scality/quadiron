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
#ifndef __NTTEC_GF_BIN_EXT_H__
#define __NTTEC_GF_BIN_EXT_H__

#include <limits>

#include "exceptions.h"
#include "gf_base.h"

namespace nttec {
namespace gf {

/** An extension Galois Field extended from GF(2)
 * @note set the 3rd argument as false to deactivate the computing of primitive
 * root on creating the base class
 */
template <typename T>
class BinExtension : public gf::Field<T> {
  public:
    explicit BinExtension(T n);
    ~BinExtension();
    void init();
    void find_primitive_root();
    T card(void) const override;
    T card_minus_one(void) const override;
    bool check(T a) const override;
    T neg(T a) const override;
    T add(T a, T b) const override;
    T sub(T a, T b) const override;
    T mul(T a, T b) const override;
    T div(T a, T b) const override;
    T inv(T a) const override;
    T exp(T a, T b) const override;
    T log(T a, T b) const override;
    void hadamard_mul(int n, T* x, T* y) const override;
    void hadamard_mul_doubled(int n, T* x, T* y) const override;
    void add_doubled(int n, T* x, T* y) const override;

  private:
    T n;
    T my_card;
    T sgroup_nb;
    T primitive_poly;
    T first_bit;
    T tab_nb;
    T* gflog = nullptr;
    T* gfilog = nullptr;
    T*** gfsplit = nullptr; // (n/4-1)*256*256 elements
    T* mask = nullptr;
    bool restricted = false;
    T _mul_log(T a, T b) const;
    T _mul_split(T a, T b) const;
    T _mul_by_two(T x) const;
    T _shift_left(T x, T shift) const;
    T _deg_of(T x, T max_deg) const;
    T _div_log(T a, T b) const;
    T _div_by_inv(T a, T b) const;
    T _inv_by_div(T a) const;
    T _inv_ext_gcd(T a) const;
    int mul_type;
    int div_type;
    int inv_type;
    void init_mask(void);
    void setup_tables(void);
    void setup_split_tables(void);
};

enum MulType { MUL_LOG_TAB, SPLIT_8_8 };
enum DivType { DIV_LOG_TAB, DIV_BY_INV };
enum InvType { INV_BY_DIV, INV_EXT_GCD };

template <typename T>
BinExtension<T>::BinExtension(T n) : gf::Field<T>(2, n, false)
{
    this->n = n;
    if (n / 8 > sizeof(T)) {
        assert(false);
    }
    if (n == 1)
        this->primitive_poly = 0x1;
    else if (n == 2)
        this->primitive_poly = 0x7;
    else if (n == 3)
        this->primitive_poly = 0xb;
    else if (n == 4)
        this->primitive_poly = 0x13;
    else if (n == 8)
        // alternative 0x11b, original one: 0x11d
        this->primitive_poly = 0x11d;
    else if (n == 16)
        // an alternative: 0x1002b
        this->primitive_poly = T(0x1100b);
    // Following irreducible polynomials are from Gadiel Seroussi's paper:
    //  "Table of Low-Weight Binary Irreducible Polynomials"
    else if (n == 32)
        // pentanomial x^32 + x^7 + x^3 + x^2 + 1
        this->primitive_poly = 0x8d;
    else if (n == 64)
        // pentanomial x^64 + x^4 + x^3 + x + 1
        this->primitive_poly = 0x1b;
    else if (n == 128)
        // pentanomial x^128 + x^7 + x^2 + x + 1
        this->primitive_poly = 0x87;
    else
        assert(false); // not supported

    if (n / 8 == sizeof(T)) {
        // restricted case
        this->restricted = true;
    }

    this->mask = new T[n];
    init_mask();

    this->tab_nb = n / 4 - 1;
    this->first_bit = this->mask[n - 1];
    if (this->restricted) {
        if (n <= 64)
            this->my_card = std::numeric_limits<T>::max();
        else {
            this->my_card = 0;
            for (T i = 0; i < n; i++)
                this->my_card += this->mask[i];
        }
    } else
        this->my_card = 2 * this->mask[n - 1];

    if (n <= 16) {
        setup_tables();
        this->mul_type = MUL_LOG_TAB;
        this->div_type = DIV_LOG_TAB;
        this->inv_type = INV_BY_DIV;
    } else {
        // currently only for (8, 8) split
        this->mul_type = SPLIT_8_8;
        this->div_type = DIV_BY_INV;
        this->inv_type = INV_EXT_GCD;
        this->sgroup_nb = n / 8;
        setup_split_tables();
    }

    // computing primitive root
    this->init();
}

template <typename T>
BinExtension<T>::~BinExtension()
{
    if (gflog != nullptr)
        delete[] gflog;
    if (gfilog != nullptr)
        delete[] gfilog;
    if (gfsplit != nullptr) {
        for (T t = 0; t < tab_nb; t++) {
            for (T i = 0; i < 256; i++) {
                delete[] gfsplit[t][i];
            }
            delete[] gfsplit[t];
        }
        delete[] gfsplit;
    }
    if (mask)
        delete[] mask;
}

template <typename T>
void BinExtension<T>::init()
{
    // compute factors of order
    RingModN<T>::compute_factors_of_order();

    // compute root of unity
    this->find_primitive_root();
}

template <typename T>
void BinExtension<T>::find_primitive_root()
{
    if (this->root) {
        return;
    }

    T h = BinExtension<T>::card_minus_one();
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
            if (BinExtension<T>::exp(nb, this->proper_divisors->at(i)) == 1) {
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
void BinExtension<T>::init_mask(void)
{
    T i, v;

    v = 1;
    for (i = 0; i < n; i++, v *= 2) {
        mask[i] = v;
    }
}

template <typename T>
void BinExtension<T>::setup_tables(void)
{
    T b, log;

    gflog = new T[my_card];
    gfilog = new T[my_card];
    b = 1;
    for (log = 0; log < my_card - 1; log++) {
        gflog[b] = log;
        gfilog[log] = b;
        b = b << 1;
        if (b & my_card)
            b = b ^ primitive_poly;
    }
}

template <typename T>
inline T BinExtension<T>::_mul_by_two(T x) const
{
    if ((x & first_bit) > 0) {
        return ((x ^ first_bit) * 2) ^ primitive_poly;
    }
    return x * 2;
}

template <typename T>
inline T BinExtension<T>::_shift_left(T x, T shift) const
{
    if (shift == 0) {
        return x;
    } else if (shift == 1) {
        return _mul_by_two(x);
    } else if (shift > 1) {
        return _shift_left(_mul_by_two(x), shift - 1);
    } else {
        throw InvalidArgument("shift-left of negative nb");
    }
}

/**
 * Setup tables for split multiplications
 *  There are (n/4 - 1) tables each of 256 x 256 elements of GF(2^n), hence each
 *    of size (n/8 * 2^16) bytes. Total memory for tables is 2n(n-4) KB.
 *  Table t = 0, .., (n/4-2) contains multiplication of a * b * x^(8t) for
 *    a, b = 0, ..., 256
 */
template <typename T>
void BinExtension<T>::setup_split_tables(void)
{
    T i, j, t;
    T x;
    T base;

    // alloc
    gfsplit = new T**[tab_nb];
    for (t = 0; t < tab_nb; t++) {
        gfsplit[t] = new T*[256];
        for (i = 0; i < 256; i++) {
            gfsplit[t][i] = new T[256];
        }
    }

    base = 1; // x^0
    for (t = 0; t < tab_nb; t++) {
        // setup for a = 0 and b = 0
        for (j = 0; j < 256; j++) {
            gfsplit[t][0][j] = 0;
            gfsplit[t][j][0] = 0;
        }
        // setup gfsplit[t][i][1]
        gfsplit[t][1][1] = base;
        for (i = 2; i < 256; i++) {
            if (i & 1) {
                x = gfsplit[t][i ^ 1][1];
                gfsplit[t][i][1] = x ^ base;
            } else {
                x = gfsplit[t][i >> 1][1];
                gfsplit[t][i][1] = _mul_by_two(x);
            }
        }
        // setup gfsplit[t][i][j]
        for (i = 1; i < 256; i++) {
            x = gfsplit[t][i][1];
            for (j = 2; j < 256; j++) {
                if (j & 1) {
                    gfsplit[t][i][j] = gfsplit[t][i][j ^ 1] ^ x;
                } else {
                    gfsplit[t][i][j] = _mul_by_two(gfsplit[t][i][j >> 1]);
                }
            }
        }
        // update base for next table: base * x^8
        for (i = 0; i < 8; i++)
            base = _mul_by_two(base);
    }
}

template <typename T>
T BinExtension<T>::card(void) const
{
    return my_card;
}

template <typename T>
T BinExtension<T>::card_minus_one(void) const
{
    if (this->restricted)
        return this->my_card;
    return (this->my_card - 1);
}

template <typename T>
bool BinExtension<T>::check(T a) const
{
    if (this->restricted)
        return (a >= 0 && a <= my_card);
    return (a >= 0 && a < my_card);
}

template <typename T>
T BinExtension<T>::neg(T a) const
{
    assert(check(a));

    return sub(0, a);
}

template <typename T>
T BinExtension<T>::add(T a, T b) const
{
    assert(check(a));
    assert(check(b));

    return a ^ b;
}

template <typename T>
T BinExtension<T>::sub(T a, T b) const
{
    assert(check(a));
    assert(check(b));

    return a ^ b;
}

template <typename T>
T BinExtension<T>::mul(T a, T b) const
{
    assert(check(a));
    assert(check(b));

    switch (mul_type) {
    case SPLIT_8_8:
        return _mul_split(a, b);
    case MUL_LOG_TAB:
    default:
        return _mul_log(a, b);
    }
}

template <typename T>
T BinExtension<T>::_mul_log(T a, T b) const
{
    assert(check(a));
    assert(check(b));

    T sum_log;

    if (a == 0 || b == 0)
        return 0;
    sum_log = gflog[a] + gflog[b];
    if (sum_log >= my_card - 1)
        sum_log -= my_card - 1;

    return gfilog[sum_log];
}

template <typename T>
T BinExtension<T>::_mul_split(T a, T b) const
{
    assert(check(a));
    assert(check(b));

    T tb, product;
    T mask;
    T i, j;

    product = 0;
    mask = 0xff;
    for (i = 0; i < sgroup_nb; i++) {
        tb = b;
        for (j = 0; j < sgroup_nb; j++) {
            product ^= gfsplit[i + j][a & mask][tb & mask];
            tb >>= 8;
        }
        a >>= 8;
    }
    return product;
}

template <typename T>
T BinExtension<T>::div(T a, T b) const
{
    assert(check(a));
    assert(check(b));

    switch (div_type) {
    case DIV_BY_INV:
        return _div_by_inv(a, b);
    case DIV_LOG_TAB:
    default:
        return _div_log(a, b);
    }
}

template <typename T>
T BinExtension<T>::_div_by_inv(T a, T b) const
{
    assert(check(a));
    assert(check(b));

    T inv = _inv_ext_gcd(b);
    return mul(a, inv);
}

template <typename T>
T BinExtension<T>::_div_log(T a, T b) const
{
    assert(check(a));
    assert(check(b));

    int diff_log;

    if (a == 0)
        return 0;
    if (b == 0)
        return -1;
    diff_log = gflog[a] - gflog[b];
    if (diff_log < 0)
        diff_log += my_card - 1;

    return gfilog[diff_log];
}

template <typename T>
T BinExtension<T>::inv(T a) const
{
    assert(check(a));

    switch (inv_type) {
    case INV_EXT_GCD:
        return _inv_ext_gcd(a);
    case INV_BY_DIV:
    default:
        return _inv_by_div(a);
    }
}

template <typename T>
T BinExtension<T>::_inv_by_div(T a) const
{
    assert(check(a));

    return div(1, a);
}

template <typename T>
T BinExtension<T>::exp(T a, T b) const
{
    assert(BinExtension<T>::check(a));
    assert(BinExtension<T>::check(b));

    return gf::Field<T>::exp_quick(a, b);
}

template <typename T>
T BinExtension<T>::log(T a, T b) const
{
    assert(check(a));

    T result;
    T tmp = a;
    if (b == 1)
        return 0;
    for (result = 1; result < my_card; result++) {
        if (tmp == b)
            return result;
        tmp = mul(tmp, a);
    }

    throw NoSolution("solution not found");
}

template <typename T>
inline T BinExtension<T>::_deg_of(T a, T max_deg) const
{
    T deg = max_deg;
    while ((mask[deg] & a) == 0)
        deg--;
    return deg;
}

/**
 * Implement Algo 2.48 in "Guide to Elliptic Curve Cryptography",
 *  Darrel Hankerson, Scott Vanstone, Alfred Menezes
 */
template <typename T>
T BinExtension<T>::_inv_ext_gcd(T x) const
{
    T uv[2];
    T g[2];
    int deg[2];
    int a, b;
    int j;

    a = 0;
    b = 1;
    uv[a] = x;
    uv[b] = 0;
    g[a] = 1;
    g[b] = 0;

    deg[a] = _deg_of(uv[a], n - 1);
    deg[b] = n;
    while (uv[a] != 1) {
        j = deg[a] - deg[b];
        if (j < 0) {
            a ^= 1;
            b ^= 1;
            j = abs(j);
        }
        uv[a] ^= _shift_left(uv[b], j);
        g[a] ^= _shift_left(g[b], j);
        deg[a] = _deg_of(uv[a], n - 1);
        deg[b] = _deg_of(uv[b], n - 1);
    }
    return g[a];
}

template <typename T>
void BinExtension<T>::hadamard_mul(int n, T* x, T* y) const
{
    // this->hadamard_mul(n, x, y);
    for (int i = 0; i < n; i++) {
        x[i] = mul(x[i], y[i]);
    }
}

template <typename T>
void BinExtension<T>::hadamard_mul_doubled(int n, T* x, T* y) const
{
    const int half = n / 2;
    T* x_next = x + half;

    // multiply y to the first half of `x`
    for (int i = 0; i < half; i++) {
        x[i] = mul(x[i], y[i]);
    }

    // multiply y to the second half of `x`
    for (int i = 0; i < half; i++) {
        x_next[i] = mul(x_next[i], y[i]);
    }
}

template <typename T>
void BinExtension<T>::add_doubled(int n, T* x, T* y) const
{
    const int half = n / 2;
    T* x_next = x + half;

    // add y to the first half of `x`
    for (int i = 0; i < half; i++) {
        x[i] = add(x[i], y[i]);
    }

    // add y to the second half of `x`
    for (int i = 0; i < half; i++) {
        x_next[i] = add(x_next[i], y[i]);
    }
}
} // namespace gf
} // namespace nttec

#endif
