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
#ifndef __QUAD_GF_NF4_H__
#define __QUAD_GF_NF4_H__

#include <iostream>
#include <vector>

#include "gf_base.h"
#include "gf_prime.h"
#include "vec_vector.h"

#define MASK16 0xFFFF
#define MASK32 0xFFFFFFFF

namespace quadiron {
namespace gf {

/** A group of `n` elements of GF(F<sub>4</sub>=65537).
 *
 * It allows performing `n` operations on GF(F<sub>4</sub>) in parallel.
 */
template <typename T>
class NF4 : public gf::Field<T> {
  public:
    using gf::Field<T>::neg;

    NF4(NF4&& other)
    noexcept
        : Field<T>(std::move(other)), unit(other.unit), q(other.q), h(other.h),
          sub_field(std::move(other.sub_field))
    {
    }

    ~NF4() = default;
    T card(void) const override;
    T card_minus_one(void) const override;
    T get_inv_n_mod_p(int n) const override;
    bool check(T a) const override;
    T neg(T a) const override;
    T add(T a, T b) const override;
    T sub(T a, T b) const override;
    T mul(T a, T b) const override;
    T div(T a, T b) const override;
    T inv(T a) const override;
    T exp(T a, T b) const override;
    T log(T a, T b) const override;
    T unpacked_rand(void) const;
    T rand(void) const override;
    T get_unit(void) const override;
    T replicate(T a) const override;
    T pack(T a) const;
    T pack(T a, uint32_t flag) const;
    GroupedValues<T> unpack(T a) const;
    void unpack(T a, GroupedValues<T>& b) const;
    T get_nth_root(T n) const override;
    void compute_omegas(vec::Vector<T>& W, int n, T w) const override;
    const gf::Field<uint32_t>& get_sub_field() const;
    void hadamard_mul(int n, T* x, T* y) const override;

  private:
    T unit;
    T q;
    T h;
    std::unique_ptr<gf::Field<uint32_t>> sub_field;

    // Scratch space for arithmetic operations.
    // Note: those are not part of the externally visible state of the class.
    mutable std::vector<uint16_t> scratch16;
    mutable std::vector<uint32_t> scratch32;

    bool check_n(unsigned n);
    explicit NF4(unsigned n);

    template <typename Class, typename... Args>
    friend Class create(Args... args);

    template <typename Base, typename Class, typename... Args>
    friend std::unique_ptr<Base> alloc(Args... args);

    T expand16(uint16_t* arr) const;
    T expand32(uint32_t* arr) const;
    // to debug
    void show_arr(uint32_t* arr);
};

template <typename T>
NF4<T>::NF4(unsigned n) : gf::Field<T>(T(65537), n)
{
    this->isNF4 = true;
    sub_field = gf::alloc<gf::Field<uint32_t>, gf::Prime<uint32_t>>(T(65537));

    if (!check_n(n)) {
        // not supported yet
        assert(false && "Input n is not supported for NF4");
    }

    unit = NF4<T>::replicate(1);
    q = NF4<T>::replicate(T(65537));
    h = NF4<T>::replicate(T(65536));

    scratch16.reserve(this->n);
    scratch32.reserve(this->n);
}

template <typename T>
inline T NF4<T>::get_inv_n_mod_p(int n) const
{
    return replicate(gf::Field<T>::get_inv_n_mod_p(n));
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
inline T NF4<T>::expand16(uint16_t* arr) const
{
    T c = arr[this->n - 1];
    for (int i = this->n - 2; i >= 0; i--) {
        c = (c << 16) | arr[i];
    }
    return c;
}

template <typename T>
inline T NF4<T>::expand32(uint32_t* arr) const
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
inline T NF4<T>::replicate(T a) const
{
    T b = a;
    for (int i = 1; i < this->n; i++) {
        b = ((b << 16) << 16) | a;
    }
    return b;
}

template <typename T>
inline T NF4<T>::card(void) const
{
    return q;
}

template <typename T>
inline T NF4<T>::card_minus_one(void) const
{
    return h;
}

template <typename T>
inline T NF4<T>::get_unit(void) const
{
    return unit;
}

template <typename T>
inline T NF4<T>::neg(T a) const
{
    return sub(0, a);
}

template <typename T>
inline T NF4<T>::add(T a, T b) const
{
    scratch32[0] =
        (narrow_cast<uint32_t>(a) + narrow_cast<uint32_t>(b)) % 65537;
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        scratch32[i] =
            (narrow_cast<uint32_t>(a) + narrow_cast<uint32_t>(b)) % 65537;
    }

    T c = expand32(scratch32.data());

    return c;
}

template <typename T>
inline T NF4<T>::sub(T a, T b) const
{
    uint32_t ae, be;

    ae = narrow_cast<uint32_t>(a);
    be = narrow_cast<uint32_t>(b);
    scratch32[0] = ae >= be ? ae - be : 65537 + ae - be;
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        ae = narrow_cast<uint32_t>(a);
        be = narrow_cast<uint32_t>(b);
        scratch32[i] = ae >= be ? ae - be : 65537 + ae - be;
    }

    T c = expand32(scratch32.data());

    return c;
}

template <typename T>
inline T NF4<T>::mul(T a, T b) const
{
    uint64_t ae;
    uint32_t be;

    ae = static_cast<uint64_t>(a & MASK32);
    be = narrow_cast<uint32_t>(b);
    scratch32[0] = (ae == 65536 && be == 65536) ? 1 : (ae * be) % 65537;
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        ae = static_cast<uint64_t>(a & MASK32);
        be = narrow_cast<uint32_t>(b);
        scratch32[i] = (ae == 65536 && be == 65536) ? 1 : (ae * be) % 65537;
    }

    T c = expand32(scratch32.data());
    return c;
}

template <typename T>
inline T NF4<T>::div(T a, T b) const
{
    scratch32[0] =
        sub_field->div(narrow_cast<uint32_t>(a), narrow_cast<uint32_t>(b));
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        scratch32[i] =
            sub_field->div(narrow_cast<uint32_t>(a), narrow_cast<uint32_t>(b));
    }

    T c = expand32(scratch32.data());
    return c;
}

template <typename T>
inline T NF4<T>::inv(T a) const
{
    scratch32[0] = sub_field->inv(narrow_cast<uint32_t>(a));
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        scratch32[i] = sub_field->inv(narrow_cast<uint32_t>(a));
    }

    T c = expand32(scratch32.data());
    return c;
}

template <typename T>
inline T NF4<T>::exp(T a, T b) const
{
    scratch32[0] =
        sub_field->exp(narrow_cast<uint32_t>(a), narrow_cast<uint32_t>(b));
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        scratch32[i] =
            sub_field->exp(narrow_cast<uint32_t>(a), narrow_cast<uint32_t>(b));
    }

    T c = expand32(scratch32.data());
    return c;
}

template <typename T>
inline T NF4<T>::log(T a, T b) const
{
    scratch32[0] =
        sub_field->log(narrow_cast<uint32_t>(a), narrow_cast<uint32_t>(b));
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        b = (b >> 16) >> 16;
        scratch32[i] =
            sub_field->log(narrow_cast<uint32_t>(a), narrow_cast<uint32_t>(b));
    }

    T c = expand32(scratch32.data());
    return c;
}

template <typename T>
T NF4<T>::rand() const
{
    T c = sub_field->rand();
    for (int i = 1; i < this->n; i++) {
        c = ((c << 16) << 16) | sub_field->rand();
    }
    return c;
}

template <typename T>
T NF4<T>::unpacked_rand(void) const
{
    T c = rand();
    return unpack(c).values;
}

/**
 * Pack of n numbers each of 16 bits into n numbers each of 32 bits
 */
template <typename T>
inline T NF4<T>::pack(T a) const
{
    scratch32[0] = static_cast<uint32_t>(a & MASK16);
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16);
        scratch32[i] = static_cast<uint32_t>(a & MASK16);
    }

    T c = expand32(scratch32.data());
    return c;
}

/**
 * Pack of n numbers each of 16 bits into n numbers each of 32 bits
 *  If flag contains 2^i, ith number == 65537
 */
template <typename T>
inline T NF4<T>::pack(T a, uint32_t flag) const
{
    scratch32[0] = (flag & 1) ? 65536 : static_cast<uint32_t>(a & MASK16);
    for (int i = 1; i < this->n; i++) {
        flag >>= 1;
        a = (a >> 16);
        scratch32[i] = (flag & 1) ? 65536 : static_cast<uint32_t>(a & MASK16);
    }

    T c = expand32(scratch32.data());
    return c;
}

/**
 * Unpack of n numbers each of 32 bits into n numbers each of 16 bits
 *  If ith number == 65537, add 2^i to flag to mark it, and consider this number
 *  as zero
 */
template <typename T>
inline GroupedValues<T> NF4<T>::unpack(T a) const
{
    GroupedValues<T> b = GroupedValues<T>();
    uint32_t flag = 0;
    uint32_t ae;

    ae = narrow_cast<uint32_t>(a);
    if (ae == 65536) {
        flag |= 1;
        scratch16[0] = 0;
    } else {
        scratch16[0] = narrow_cast<uint16_t>(ae);
    }
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        ae = narrow_cast<uint32_t>(a);
        if (ae == 65536) {
            flag |= (1 << i);
            scratch16[i] = 0;
        } else {
            scratch16[i] = ae;
        }
    }

    b.flag = flag;
    b.values = expand16(scratch16.data());
    return b;
}

template <typename T>
inline void NF4<T>::unpack(T a, GroupedValues<T>& b) const
{
    uint32_t flag = 0;
    uint32_t ae;

    ae = narrow_cast<uint32_t>(a);
    if (ae == 65536) {
        flag |= 1;
        scratch16[0] = 0;
    } else {
        scratch16[0] = narrow_cast<uint16_t>(ae);
    }
    for (int i = 1; i < this->n; i++) {
        a = (a >> 16) >> 16;
        ae = narrow_cast<uint32_t>(a);
        if (ae == 65536) {
            flag |= (1 << i);
            scratch16[i] = 0;
        } else {
            scratch16[i] = ae;
        }
    }

    b.flag = flag;
    b.values = expand16(scratch16.data());
}

// Use for fft
template <typename T>
inline T NF4<T>::get_nth_root(T n) const
{
    T sub_w = sub_field->get_nth_root(n);
    T w = replicate(sub_w);
    return w;
}

template <typename T>
bool NF4<T>::check(T a) const
{
    return (a >= 0 && a < q);
}

/** Compute the different powers of the root of unity into a vector.
 *
 * @param W output vector (must be of length n+1)
 * @param n the length of the output vector
 * @param w n-th root of unity
 */
template <typename T>
inline void NF4<T>::compute_omegas(vec::Vector<T>& W, int n, T w) const
{
    for (int i = 0; i < n; i++) {
        W.set(i, this->exp(w, replicate(i)));
    }
}

template <typename T>
const gf::Field<uint32_t>& NF4<T>::get_sub_field() const
{
    return *(sub_field.get());
}

template <typename T>
inline void NF4<T>::hadamard_mul(int n, T* x, T* y) const
{
    for (int i = 0; i < n; i++) {
        x[i] = mul(x[i], y[i]);
    }
}

#ifdef QUADIRON_USE_SIMD
/* Operations are vectorized by SIMD */

template <>
__uint128_t NF4<__uint128_t>::expand16(uint16_t* arr) const;

template <>
__uint128_t NF4<__uint128_t>::expand32(uint32_t* arr) const;

template <>
__uint128_t NF4<__uint128_t>::add(__uint128_t a, __uint128_t b) const;

template <>
__uint128_t NF4<__uint128_t>::sub(__uint128_t a, __uint128_t b) const;

template <>
__uint128_t NF4<__uint128_t>::mul(__uint128_t a, __uint128_t b) const;

template <>
__uint128_t NF4<__uint128_t>::pack(__uint128_t a) const;

template <>
__uint128_t NF4<__uint128_t>::pack(__uint128_t a, uint32_t flag) const;

template <>
GroupedValues<__uint128_t> NF4<__uint128_t>::unpack(__uint128_t a) const;

template <>
void NF4<__uint128_t>::unpack(__uint128_t a, GroupedValues<__uint128_t>& b)
    const;

template <>
void NF4<__uint128_t>::hadamard_mul(int n, __uint128_t* x, __uint128_t* y)
    const;

#endif // #ifdef QUADIRON_USE_SIMD

} // namespace gf
} // namespace quadiron

#endif
