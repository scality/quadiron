/* -*- mode: c++ -*- */
#ifndef __NTL_VMVEC_H__
#define __NTL_VMVEC_H__

#include "vec.h"

/**
 * Virtual vector on a base vector from an offset
 *  vmvec[i] = base[i + offset]
 *
 */
template <typename T>
class VmVec : public Vec<T> {
  private:
    T* mem_offset;
    int base_len;
    int offset;

  public:
    explicit VmVec(Vec<T>* vec, int n = 0, int offset = 0);
    T get(int i);
    void set(int i, T val);
    void set_map(int offset);
    int get_offset(void);
};

template <typename T>
VmVec<T>::VmVec(Vec<T>* vec, int n, int offset)
    : Vec<T>(
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
inline T VmVec<T>::get(int i)
{
    assert(i >= 0 && i < this->n);

    return mem_offset[i];
}

template <typename T>
void VmVec<T>::set_map(int offset)
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
inline void VmVec<T>::set(int i, T val)
{
    assert(i >= 0 && i < this->n);

    mem_offset[i] = val;
}

template <typename T>
int VmVec<T>::get_offset()
{
    return this->offset;
}

#endif
