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
#ifndef __NTTEC_GF_NF4_H__
#define __NTTEC_GF_NF4_H__

#include <iostream>

#include "gf_base.h"
#include "gf_prime.h"
#include "vec_vector.h"

#define MASK16 0xFFFF
#define MASK32 0xFFFFFFFF

namespace nttec {
namespace gf {

/** A group of `n` elements of GF(F<sub>4</sub>=65537).
 *
 * It allows performing `n` operations on GF(F<sub>4</sub>) in parallel.
 */
template <typename T>
class NF4 : public gf::Field<T> {
  public:
    explicit NF4(unsigned n);
    ~NF4();
    T card(void);
    T card_minus_one(void);
    bool check(T a);
    T neg(T a);
    T add(T a, T b);
    T sub(T a, T b);
    T mul(T a, T b);
    T div(T a, T b);
    T inv(T a);
    T exp(T a, T b);
    T log(T a, T b);
    T weak_rand_tuple(void);
    T weak_rand(void);
    T get_unit(void);
    T replicate(T a);
    T pack(T a);
    T pack(T a, uint32_t flag);
    GroupedValues<T> unpack(T a);
    T get_nth_root(T n);
    void compute_omegas(vec::Vector<T>* W, int n, T w);
    gf::Field<uint32_t>* get_sub_field();
    void hadamard_mul(int n, T* x, T* y);
    void add(int n, T* x, T* y);

  private:
    T unit;
    T q;
    T h;
    gf::Field<uint32_t>* sub_field;
    bool check_n(unsigned n);
    void init(void);

    T expand16(uint16_t* arr);
    T expand32(uint32_t* arr);
    // to debug
    void show_arr(uint32_t* arr);
};

template <typename T>
NF4<T>::NF4(unsigned n) : gf::Field<T>(T(65537), n)
{
    sub_field = new gf::Prime<uint32_t>(T(65537));

    if (!check_n(n)) {
        // not supported yet
        assert(false && "Input n is not supported for NF4");
    }

    unit = 0;
    q = 0;
    h = 0;
    init();
}

template <typename T>
NF4<T>::~NF4()
{
    delete sub_field;
}

/*
 * Each element of sub-field GFP(65537) is stored in 32bits
 * Hence, n <= sizeof(T) / 4
 */
template <typename T>
bool NF4<T>::check_n(unsigned n)
{
    return (n <= sizeof(T) / 4);
}

template <typename T>
T NF4<T>::expand16(uint16_t* arr)
{
    T c = arr[this->n - 1];
    for (int i = this->n - 2; i >= 0; i--) {
        c = (c << 16) | arr[i];
    }
    return c;
}

template <typename T>
T NF4<T>::expand32(uint32_t* arr)
{
    T c = arr[this->n - 1];
    for (int i = this->n - 2; i >= 0; i--) {
        c = ((c << 16) << 16) | arr[i];
    }
    return c;
}

template <typename T>
void NF4<T>::show_arr(uint32_t* arr)
{
    std::cout << "\t arr: ";
    for (int i = 0; i < this->n; i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

template <typename T>
T NF4<T>::replicate(T a)
{
    T b = a;
    for (int i = 1; i < this->n; i++) {
        b = ((b << 16) << 16) | a;
    }
    return b;
}

template <typename T>
void NF4<T>::init(void)
{
    unit = replicate(1);
    q = replicate(T(65537));
    h = replicate(T(65536));
}

template <typename T>
T NF4<T>::card(void)
{
    return q;
}

template <typename T>
T NF4<T>::card_minus_one(void)
{
    return h;
}

template <typename T>
T NF4<T>::get_unit(void)
{
    return unit;
}

template <typename T>
T NF4<T>::neg(T a)
{
    return sub(0, a);
}

template <typename T>
T NF4<T>::add(T a, T b)
{
    uint32_t arr[this->n];

    arr[0] = ((uint32_t)(a & MASK32) + (uint32_t)(b & MASK32)) % 65537;
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        arr[i] = ((uint32_t)(a & MASK32) + (uint32_t)(b & MASK32)) % 65537;
    }

    T c = expand32(arr);

    return c;
}

template <typename T>
T NF4<T>::sub(T a, T b)
{
    uint32_t arr[this->n];
    uint32_t ae, be;

    ae = (a & MASK32);
    be = (b & MASK32);
    if (ae >= be)
        arr[0] = ae - be;
    else
        arr[0] = 65537 + ae - be;
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        ae = (uint32_t)(a & MASK32);
        be = (uint32_t)(b & MASK32);
        if (ae >= be)
            arr[i] = ae - be;
        else
            arr[i] = 65537 + ae - be;
    }

    T c = expand32(arr);

    return c;
}

template <typename T>
T NF4<T>::mul(T a, T b)
{
    uint32_t arr[this->n];
    uint64_t ae;
    uint32_t be;

    ae = (uint64_t)(a & MASK32);
    be = (uint32_t)(b & MASK32);
    if (ae == 65536 && be == 65536)
        arr[0] = 1;
    else
        arr[0] = (ae * be) % 65537;
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        ae = (uint64_t)(a & MASK32);
        be = (uint32_t)(b & MASK32);
        if (ae == 65536 && be == 65536)
            arr[i] = 1;
        else
            arr[i] = (ae * be) % 65537;
    }

    T c = expand32(arr);
    return c;
}

template <typename T>
T NF4<T>::div(T a, T b)
{
    uint32_t arr[this->n];

    arr[0] = sub_field->div((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        arr[i] = sub_field->div((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
    }

    T c = expand32(arr);
    return c;
}

template <typename T>
T NF4<T>::inv(T a)
{
    uint32_t arr[this->n];

    arr[0] = sub_field->inv((uint32_t)(a & MASK32));
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        arr[i] = sub_field->inv((uint32_t)(a & MASK32));
    }

    T c = expand32(arr);
    return c;
}

template <typename T>
T NF4<T>::exp(T a, T b)
{
    uint32_t arr[this->n];

    arr[0] = sub_field->exp((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        arr[i] = sub_field->exp((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
    }

    T c = expand32(arr);
    return c;
}

template <typename T>
T NF4<T>::log(T a, T b)
{
    uint32_t arr[this->n];

    arr[0] = sub_field->log((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        arr[i] = sub_field->log((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
    }

    T c = expand32(arr);
    return c;
}

template <typename T>
T NF4<T>::weak_rand_tuple()
{
    T c = sub_field->weak_rand();
    for (int i = 1; i < this->n; i++) {
        c = ((c << 16) << 16) | sub_field->weak_rand();
    }
    return c;
}

template <typename T>
T NF4<T>::weak_rand(void)
{
    T c = weak_rand_tuple();
    return unpack(c).values;
}

/**
 * Pack of n numbers each of 16 bits into n numbers each of 32 bits
 */
template <typename T>
T NF4<T>::pack(T a)
{
    uint32_t arr[this->n];
    arr[0] = (uint32_t)(a & MASK16);
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16);
        arr[i] = (uint32_t)(a & MASK16);
    }

    T c = expand32(arr);
    return c;
}

/**
 * Pack of n numbers each of 16 bits into n numbers each of 32 bits
 *  If flag contains 2^i, ith number == 65537
 */
template <typename T>
T NF4<T>::pack(T a, uint32_t flag)
{
    uint32_t arr[this->n];
    if (flag & 1)
        arr[0] = 65536;
    else
        arr[0] = (uint32_t)(a & MASK16);
    for (int i = 1; i < this->n; i++) {
        flag >>= 1;
        a = (a >> 16);
        if (flag & 1)
            arr[i] = 65536;
        else
            arr[i] = (uint32_t)(a & MASK16);
    }

    T c = expand32(arr);
    return c;
}

/**
 * Unpack of n numbers each of 32 bits into n numbers each of 16 bits
 *  If ith number == 65537, add 2^i to flag to mark it, and consider this number
 *  as zero
 */
template <typename T>
GroupedValues<T> NF4<T>::unpack(T a)
{
    GroupedValues<T> b = GroupedValues<T>();
    uint32_t flag = 0;
    uint32_t ae;
    uint16_t arr[this->n];

    ae = (uint32_t)(a & MASK32);
    if (ae == 65536) {
        flag |= 1;
        arr[0] = 0;
    } else
        arr[0] = (uint16_t)ae;
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        ae = (uint32_t)(a & MASK32);
        if (ae == 65536) {
            flag |= (1 << i);
            arr[i] = 0;
        } else
            arr[i] = ae;
    }

    b.flag = flag;
    b.values = expand16(arr);
    return b;
}

// Use for fft
template <typename T>
T NF4<T>::get_nth_root(T n)
{
    T sub_w = sub_field->get_nth_root(n);
    T w = replicate(sub_w);
    return w;
}

template <typename T>
bool NF4<T>::check(T a)
{
    return (a >= 0 && a < q);
}

/**
 * Compute the different powers of the root of unity into a vector
 *
 * @param W output vector (must be of length n+1)
 * @param n
 * @param w n-th root of unity
 */
template <typename T>
void NF4<T>::compute_omegas(vec::Vector<T>* W, int n, T w)
{
    for (int i = 0; i < n; i++) {
        W->set(i, this->exp(w, replicate(i)));
    }
}

template <typename T>
gf::Field<uint32_t>* NF4<T>::get_sub_field()
{
    return sub_field;
}

template <typename T>
void NF4<T>::hadamard_mul(int n, T* x, T* y)
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
void NF4<T>::add(int n, T* x, T* y)
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

#ifdef NTTEC_USE_SIMD
/* Operations are vectorized by SIMD */

template <>
__uint128_t NF4<__uint128_t>::expand16(uint16_t* arr);

template <>
__uint128_t NF4<__uint128_t>::expand32(uint32_t* arr);

template <>
__uint128_t NF4<__uint128_t>::add(__uint128_t a, __uint128_t b);

template <>
__uint128_t NF4<__uint128_t>::sub(__uint128_t a, __uint128_t b);

template <>
__uint128_t NF4<__uint128_t>::mul(__uint128_t a, __uint128_t b);

template <>
__uint128_t NF4<__uint128_t>::pack(__uint128_t a);

template <>
__uint128_t NF4<__uint128_t>::pack(__uint128_t a, uint32_t flag);

template <>
GroupedValues<__uint128_t> NF4<__uint128_t>::unpack(__uint128_t a);

template <>
void NF4<__uint128_t>::hadamard_mul(int n, __uint128_t* x, __uint128_t* y);

template <>
void NF4<__uint128_t>::add(int n, __uint128_t* x, __uint128_t* y);

#endif // #ifdef NTTEC_USE_SIMD

} // namespace gf
} // namespace nttec

#endif
