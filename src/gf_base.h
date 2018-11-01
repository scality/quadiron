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
#ifndef __QUAD_GF_BASE_H__
#define __QUAD_GF_BASE_H__

#include "arith.h"
#include "core.h"
#include "gf_ring.h"

namespace quadiron {

/** Galois Fields handling. */
namespace gf {

template <typename T>
class Prime;

/** Generic class for Galois Fields of q=p<sup>n</sup>. */
template <typename T>
class Field : public RingModN<T> {
  public:
    Field(Field&& other) noexcept
        : RingModN<T>(std::move(other)), p(other.p), n(other.n),
          sub_field(std::move(other.sub_field))
    {
    }
    virtual ~Field() = default;

    const Field<T>& get_sub_field() const;
    T get_p() const;
    int get_n() const;
    virtual T get_inv_n_mod_p(int n) const;
    virtual T get_unit() const;

  protected:
    T p;
    int n;
    std::unique_ptr<Prime<T>> sub_field;

    Field(T p, int n);

    template <typename Class, typename... Args>
    friend Class create(Args... args);

    template <typename Base, typename Class, typename... Args>
    friend std::unique_ptr<Base> alloc(Args... args);
};

template <typename T>
Field<T>::Field(T p, int n) : RingModN<T>(arith::exp<T>(p, n))
{
    // XXX: Shall check that `p` is prime.
    this->p = p;
    this->n = n;
    if (n == 1) {
        this->sub_field = std::unique_ptr<Prime<T>>(nullptr);
    } else {
        this->sub_field = gf::alloc<gf::Prime<T>, gf::Prime<T>>(p);
    }
}

/**
 * return the field in which is based the extension field (or the field
 * itself if n == 1)
 */
template <typename T>
const Field<T>& Field<T>::get_sub_field() const
{
    if (this->sub_field)
        return *sub_field;
    else
        return *this;
}

template <typename T>
T Field<T>::get_p() const
{
    return p;
}

template <typename T>
int Field<T>::get_n() const
{
    return n;
}

template <typename T>
T Field<T>::get_inv_n_mod_p(int n) const
{
    return this->inv(n) % p;
}

template <typename T>
T Field<T>::get_unit() const
{
    return 1;
}

} // namespace gf
} // namespace quadiron

#endif
