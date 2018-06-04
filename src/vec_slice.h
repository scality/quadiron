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
#ifndef __NTTEC_VEC_SLICE_H__
#define __NTTEC_VEC_SLICE_H__

#include "vec_vector.h"

namespace nttec {
namespace vec {

/** A slice of an existing vector.
 *
 * Example:
 *
 * If you have the following Vector `vec`:
 *
 * v[0]  | v[1]   | v[2]   | v[3]   | v[4]   | v[5]
 * ----- | ------ | ------ | ------ | ------ | -----
 * v1    |  v2    |  v3    |  v4    |  v5    |  v6
 *
 * Then, Slice(vec, 3, 2) is:
 *
 * s[0]   | s[1]   | s[2]
 * ------ | ------ | -----
 *  v3    |  v4    |  v5
 */
template <typename T>
class Slice : public Vector<T> {
  public:
    explicit Slice(Vector<T>* vec, int n = 0, int offset = 0);
    T get(int i) const override;
    void set(int i, T val) override;
    void set_map(int offset);
    int get_offset(void);

  private:
    T* mem_offset;
    int base_len;
    int offset;
};

template <typename T>
Slice<T>::Slice(Vector<T>* vec, int n, int offset)
    : Vector<T>(
          vec->rn,
          n,
          vec->get_mem() + offset,
          vec->get_mem_len() > offset ? vec->get_mem_len() - offset : 0)
{
    this->mem_offset = this->get_mem();
    this->base_len = vec->get_n();
    this->offset = offset;
    this->n =
        (n + this->offset < this->base_len)
            ? n
            : (this->base_len > this->offset ? this->base_len - this->offset
                                             : 0);
}

template <typename T>
inline T Slice<T>::get(int i) const
{
    assert(i >= 0 && i < this->n);

    return mem_offset[i];
}

template <typename T>
void Slice<T>::set_map(int offset)
{
    this->mem_offset = this->mem_offset - this->offset + offset;
    this->offset = offset;
    int n = this->n;
    this->n =
        (n + this->offset < this->base_len)
            ? n
            : (this->base_len > this->offset ? this->base_len - this->offset
                                             : 0);
    this->set_mem(this->mem_offset, this->n);
}

template <typename T>
inline void Slice<T>::set(int i, T val)
{
    assert(i >= 0 && i < this->n);

    mem_offset[i] = val;
}

template <typename T>
int Slice<T>::get_offset()
{
    return this->offset;
}

} // namespace vec
} // namespace nttec

#endif
