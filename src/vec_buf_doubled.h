/* -*- mode: c++ -*- */
/*
 * Copyright 2017-2018 the NTTEC authors
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
#ifndef __NTTEC_VEC_BUF_DOUBLED_H__
#define __NTTEC_VEC_BUF_DOUBLED_H__

#include "vec_buffers.h"

namespace nttec {
namespace vec {

/** A vector of `n` buffers virtually extented to `2n`.
 *
 * Physically, BuffersDoubled contains `n` buffers but behaves as if it
 * contained `2n` buffers.
 *
 * Example:
 *
 * We have a set of `N` independent buffers:
 *
 * buf1    | buf2    | … | bufN
 * ------- | ------- | - | -------
 * buf1[0] | buf2[0] | … | bufN[0]
 * buf1[1] | buf2[1] | … | bufN[1]
 * …       | …       | … | …
 *
 * In memory, the BuffersDoubled looks like:
 *
 * v[0]  | v[1]   | … | v[n-1]
 * ----- | ------ | - | --------
 * &buf1 |  &buf2 | … | &bufN
 *
 * But it behaves likes:
 *
 * v[0]  | v[1]   | … | v[n-1] | v[n]  | v[n+1] | … | v[2n-1]
 * ----- | ------ | - | ------ | ----- | ------ | - | -----------
 * &buf1 |  &buf2 | … | &bufN  | &buf1 |  &buf2 | … | &bufN
 */
template <typename T>
class BuffersDoubled : public Buffers<T> {
  public:
    explicit BuffersDoubled(Buffers<T>* vec);
    T* get(int i);

  private:
    Buffers<T>* vec;
};

template <typename T>
BuffersDoubled<T>::BuffersDoubled(Buffers<T>* vec)
    : Buffers<T>(2 * vec->get_n(), vec->get_size(), vec->get_mem())
{
    this->vec = vec;
}

template <typename T>
T* BuffersDoubled<T>::get(int i)
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
