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
#ifndef __NTTEC_VEC_ZERO_EXT_H__
#define __NTTEC_VEC_ZERO_EXT_H__

#include "vec_vector.h"

namespace nttec {
namespace vec {

/** A vector of size `n` virtually extented with zeros.
 *
 * In memory, the vec::ZeroExtended looks like:
 *
 * v[0]  | v[1]   | … | v[n-1]
 * ----- | ------ | - | --------
 * v1    |  v2    | … | vN
 *
 * But it behaves likes:
 *
 * v[0]  | v[1]   | … | v[n-1] | v[n]  | v[n+1] | …
 * ----- | ------ | - | ------ | ----- | ------ | -----
 * v1    |  v2    | … | vN     | 0     |  0     | 0
 */
template <typename T>
class ZeroExtended : public Vector<T> {
  public:
    ZeroExtended(Vector<T>* vec, int n);
    int get_n(void);
    T get(int i);

  private:
    Vector<T>* vec;
    int n;
};

template <typename T>
ZeroExtended<T>::ZeroExtended(Vector<T>* vec, int n)
    : Vector<T>(vec->rn, n, vec->get_mem(), vec->get_mem_len())
{
    this->vec = vec;
    this->n = n;
}

template <typename T>
int ZeroExtended<T>::get_n(void)
{
    return n;
}

template <typename T>
T ZeroExtended<T>::get(int i)
{
    assert(i >= 0 && i < n);

    return (i < vec->get_n()) ? vec->get(i) : 0;
}

} // namespace vec
} // namespace nttec

#endif
