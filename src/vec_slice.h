/* -*- mode: c++ -*- */
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
    T get(int i);
    void set(int i, T val);
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
inline T Slice<T>::get(int i)
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
