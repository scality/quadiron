/* -*- mode: c++ -*- */
#ifndef __NTTEC_VVECP_H__
#define __NTTEC_VVECP_H__

#include "vec.h"

/**
 * Virtual vector that padds zero_chunk beyond managed vector length
 */
template <typename T>
class VVecp : public Vecp<T> {
  private:
    Vecp<T>* vec;
    int vec_n;
    T* zero_chunk = nullptr;

  public:
    VVecp(Vecp<T>* vec, int n);
    ~VVecp();
};

template <typename T>
VVecp<T>::VVecp(Vecp<T>* vec, int n)
    : Vecp<T>(n, vec->get_size(), vec->get_mem())
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
VVecp<T>::~VVecp()
{
    if (zero_chunk != nullptr) {
        delete[] zero_chunk;
    }
    this->mem->shrink_to_fit();
    delete this->mem;
}

#endif
