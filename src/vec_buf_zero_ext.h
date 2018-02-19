/* -*- mode: c++ -*- */
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
    BuffersZeroExtended(Buffers<T>* vec, int n);
    ~BuffersZeroExtended();

  private:
    Buffers<T>* vec;
    int vec_n;
    T* zero_chunk = nullptr;
};

template <typename T>
BuffersZeroExtended<T>::BuffersZeroExtended(Buffers<T>* vec, int n)
    : Buffers<T>(n, vec->get_size(), vec->get_mem())
{
    this->vec = vec;
    vec_n = vec->get_n();
    assert(n >= vec_n);
    // overwrite by a new vector
    this->mem =
        new std::vector<T*>(vec->get_mem()->begin(), vec->get_mem()->end());
    if (n > vec_n) {
        zero_chunk = new T[this->size];
        std::memset(zero_chunk, 0, this->size * sizeof(T));
        this->mem->resize(n, zero_chunk);
    }
}

template <typename T>
BuffersZeroExtended<T>::~BuffersZeroExtended()
{
    if (zero_chunk != nullptr) {
        delete[] zero_chunk;
    }
    this->mem->shrink_to_fit();
    delete this->mem;
}

} // namespace vec
} // namespace nttec

#endif
