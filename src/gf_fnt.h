/* -*- mode: c++ -*- */
/*
 * Copyright 2017-2019 Scality
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
#ifndef __QUAD_GF_FNT_H__
#define __QUAD_GF_FNT_H__

#include "gf_prime.h"

namespace quadiron {
namespace gf {

template <typename T>
T Fermat_num()
{
    return 65537;
}

template <>
inline uint16_t Fermat_num<uint16_t>()
{
    return 257;
}
template <>
inline uint32_t Fermat_num<uint32_t>()
{
    return 65537;
}

/** A Galois Field whose order is a prime Fermat number. */
template <typename T>
class Fnt : public gf::Prime<T> {
  public:
    Fnt(Fnt&&) = default;
    T add(T a, T b) const override;
    T sub(T a, T b) const override;
    T mul(T a, T b) const override;
    T mul_bounded(T a, T b) const override;

  private:
    explicit Fnt();
    unsigned half_word_bits_nb;

    template <typename Class, typename... Args>
    friend Class create(Args... args);

    template <typename Base, typename Class, typename... Args>
    friend std::unique_ptr<Base> alloc(Args... args);
};

template <typename T>
Fnt<T>::Fnt() : gf::Prime<T>(Fermat_num<T>())
{
    unsigned sizeofT = sizeof(T);
    assert(sizeofT == 2 || sizeofT == 4);

    half_word_bits_nb = 4 * sizeofT;
}

template <typename T>
inline T Fnt<T>::add(T a, T b) const
{
    assert(this->check(a));
    assert(this->check(b));

    T c = a + b;
    return (c < Fermat_num<T>()) ? c : c - Fermat_num<T>();
}

template <typename T>
inline T Fnt<T>::sub(T a, T b) const
{
    assert(this->check(a));
    assert(this->check(b));

    return (b > a) ? Fermat_num<T>() - b + a : a - b;
}

template <typename T>
inline T Fnt<T>::mul(T a, T b) const
{
    assert(this->check(a));
    assert(this->check(b));

    if (a == Fermat_num<T>() && b == Fermat_num<T>()) {
        return 1;
    }
    return mul_bounded(a, b);
}

/**
 * Modular multiplication for packed unsigned 32-bit integers
 *
 * @note We assume that at least `a` or `b` is less than `q-1` so it's
 * not necessary to verify overflow on multiplying elements
 *
 * @param a input register
 * @param b input register
 * @return (a * b) mod q
 */
template <typename T>
inline T Fnt<T>::mul_bounded(T a, T b) const
{
    assert(this->check(a));
    assert(this->check(b));

    T c = a * b;
    T l = static_cast<HalfSizeVal<T>>(c);
    c = c >> half_word_bits_nb;
    T h = static_cast<HalfSizeVal<T>>(c);
    T res = sub(l, h);

    return res;
}

} // namespace gf
} // namespace quadiron

#endif
