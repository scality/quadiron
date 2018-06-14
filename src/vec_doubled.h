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
#ifndef __NTTEC_VEC_DOUBLED_H__
#define __NTTEC_VEC_DOUBLED_H__

#include "vec_vector.h"

namespace nttec {
namespace vec {

/** A vector of size `n` virtually extented to `2n`.
 *
 * Example:
 *
 * In memory, the vec::Doubled looks like:
 *
 * v[0]  | v[1]   | … | v[n-1]
 * ----- | ------ | - | --------
 * v1    |  v2    | … | vN
 *
 * But it behaves likes:
 *
 * v[0]  | v[1]   | … | v[n-1] | v[n]  | v[n+1] | … | v[2n-1]
 * ----- | ------ | - | ------ | ----- | ------ | - | -----------
 * v1    |  v2    | … | vN     | v1    |  v2    | … | vN
 */
template <typename T>
class Doubled : public Vector<T> {
  public:
    explicit Doubled(const Vector<T>& vec);
    const int get_n(void) const override;
    const T& get(int i) const override;

  private:
    const Vector<T>* vec;
};

template <typename T>
Doubled<T>::Doubled(const Vector<T>& vec)
    : Vector<T>(vec.get_gf(), vec.get_n(), vec.get_mem(), vec.get_mem_len())
{
    this->vec = &vec;
}

template <typename T>
const int Doubled<T>::get_n(void) const
{
    return vec->get_n() * 2;
}

template <typename T>
inline const T& Doubled<T>::get(int i) const
{
    int vec_n = vec->get_n();

    assert(i >= 0 && i < 2 * vec_n);

    if (i < vec_n)
        return vec->get(i);
    else
        return vec->get(i - vec_n);
}

} // namespace vec
} // namespace nttec

#endif
