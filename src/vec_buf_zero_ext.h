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
#ifndef __NTTEC_VEC_BUF_ZERO_EXT_H__
#define __NTTEC_VEC_BUF_ZERO_EXT_H__

#include <cstring>
#include <vector>

#include "vec_buffers.h"

namespace nttec {
namespace vec {

/** A vector of `n` buffers virtually extented with a `zero_chunk`.
 *
 * A `zero_chunk` is a buffers of zeros.
 *
 * Example:
 *
 * We have a set of `N` independent buffers (+ a buffer of zeros):
 *
 * buf1    | buf2    | … | bufN    | zero_chunk
 * ------- | ------- | - | ------- | ------------
 * buf1[0] | buf2[0] | … | bufN[0] | 0
 * buf1[1] | buf2[1] | … | bufN[1] | 0
 * …       | …       | … | …       | 0
 *
 * In memory, the BuffersZeroExtended looks like:
 *
 * v[0]  | v[1]   | … | v[n-1]
 * ----- | ------ | - | --------
 * &buf1 |  &buf2 | … | &bufN
 *
 * But it behaves likes:
 *
 * v[0]  | v[1]   | … | v[n-1] | v[n]        | v[n+1]       | …
 * ----- | ------ | - | ------ | ----------- | ------------ | -----
 * &buf1 |  &buf2 | … | &bufN  | &zero_chunk |  &zero_chunk | &zero_chunk
 */
template <typename T>
class BuffersZeroExtended : public Buffers<T> {
  public:
    BuffersZeroExtended(const Buffers<T>& vec, int n);
    ~BuffersZeroExtended();

  private:
    const Buffers<T>* vec;
    int vec_n;
    T* zero_chunk = nullptr;
};

template <typename T>
BuffersZeroExtended<T>::BuffersZeroExtended(const Buffers<T>& vec, int n)
    : Buffers<T>(n, vec.get_size(), vec.get_mem())
{
    this->vec = &vec;
    vec_n = vec.get_n();
    assert(n >= vec_n);
    // overwrite by a new vector
    this->mem =
        new std::vector<T*>(vec.get_mem()->begin(), vec.get_mem()->end());
    if (n > vec_n) {
        zero_chunk = aligned_allocate<T>(this->size);
        std::memset(zero_chunk, 0, this->size * sizeof(T));
        this->mem->resize(n, zero_chunk);
    }
}

template <typename T>
BuffersZeroExtended<T>::~BuffersZeroExtended()
{
    if (zero_chunk != nullptr) {
        aligned_deallocate<T>(zero_chunk);
    }
    this->mem->shrink_to_fit();
    delete this->mem;
}

} // namespace vec
} // namespace nttec

#endif
