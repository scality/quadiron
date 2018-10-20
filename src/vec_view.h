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
#ifndef __QUAD_VEC_VIEW_H__
#define __QUAD_VEC_VIEW_H__

#include "vec_vector.h"

namespace quadiron {
namespace vec {

/** A view of an existing vector.
 *
 * Example:
 *
 * If you have the following Vector `vec`:
 *
 * v[0]  | v[1]   | v[2]   | v[3]   | v[4]   | v[5]  | v[6]  | v[7]  | v[8]
 * ----- | ------ | ------ | ------ | ------ | ----- | ----- | ----- | -----
 * v1    |  v2    |  v3    |  v4    |  v5    |  v6   |  v7   |  v8   |  v9
 *
 * Then, View(vec, 6, 1, 2) is:
 *
 * v[0]   | v[1]   | v[2]
 * ------ | ------ | -----
 *  v2    |  v4    |  v6
 */
template <typename T>
class View : public Vector<T> {
  public:
    explicit View(Vector<T>* vec, int n = 0, int offset = 0, int step = 1);
    int get_n(void) const override;
    const T& get(int i) const override;
    void set(int i, T val) override;
    void set_map(int offset, int step);
    void set_len(int n);
    void set_vec(Vector<T>* vec);

  private:
    Vector<T>* vec;
    int offset;
    int step;
    int N;
};

template <typename T>
View<T>::View(Vector<T>* vec, int n, int offset, int step)
    : Vector<T>(vec->get_gf(), n, vec->get_mem(), vec->get_mem_len())
{
    this->vec = vec;
    this->offset = offset;
    this->step = step;
    this->n = n;
    this->N = vec->get_n();
    assert(this->n <= this->N);
}

template <typename T>
int View<T>::get_n(void) const
{
    return this->n;
}

template <typename T>
inline const T& View<T>::get(int i) const
{
    assert(i >= 0 && i < this->n);

    T loc = (offset + step * i) % this->N;
    return vec->get(loc);
}

template <typename T>
void View<T>::set_map(int offset, int step)
{
    this->offset = offset;
    this->step = step;
}

template <typename T>
void View<T>::set_len(int n)
{
    assert(this->n <= this->N);
    this->n = n;
}

template <typename T>
void View<T>::set_vec(Vector<T>* vec)
{
    assert(this->n <= vec->get_n());
    this->vec = vec;
    this->rn = vec->rn;
    this->N = vec->get_n();
}

template <typename T>
void View<T>::set(int i, T val)
{
    assert(i >= 0 && i < this->n);

    T loc = (offset + step * i) % this->N;
    vec->set(loc, val);
}

} // namespace vec
} // namespace quadiron

#endif
