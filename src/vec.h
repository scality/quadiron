/* -*- mode: c++ -*- */
#ifndef __NTTEC_VEC_H__
#define __NTTEC_VEC_H__

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <iostream>

#include "v2vec.h"

namespace nttec {

template <typename T>
class Poly;

template <typename T>
class RN;

/** Low-level wrappers/helpers around memory buffer. */
namespace vec {

template <typename T>
class V2Vec;

template <typename Ts, typename Td, typename Tw>
inline void
pack_next(std::vector<Ts*>* src, std::vector<Td*>* dest, int n, size_t size)
{
    int i;
    std::vector<Tw*> tmp(n, nullptr);
    for (i = 0; i < n; i++) {
        tmp[i] = static_cast<Tw*>(static_cast<void*>(src->at(i)));
        std::copy_n(tmp[i], size, dest->at(i));
    }
    tmp.shrink_to_fit();
}

/*
 * Cast buffers of source to a type corresponding to 'word_size'
 * Copy casted elements to buffers of destination
 *
 *  sizeof(Ts) <= sizeof(Td)
 *
 * @param src: source vector
 * @param dest: destination vector
 * @param n: number of buffers per vector
 * @param size: number of elements per destination buffer
 * @param word_size: number of bytes used to store data in each element of
 *  destination buffers
 * @return
 */
template <typename Ts, typename Td>
inline void pack(
    std::vector<Ts*>* src,
    std::vector<Td*>* dest,
    int n,
    size_t size,
    size_t word_size)
{
    assert(sizeof(Td) >= word_size);
    assert(word_size % sizeof(Ts) == 0);
    // get only word_size bytes from each element
    if (word_size == 1) {
        pack_next<Ts, Td, uint8_t>(src, dest, n, size);
    } else if (word_size == 2) {
        pack_next<Ts, Td, uint16_t>(src, dest, n, size);
    } else if (word_size == 4) {
        pack_next<Ts, Td, uint32_t>(src, dest, n, size);
    } else if (word_size == 8) {
        pack_next<Ts, Td, uint64_t>(src, dest, n, size);
    } else if (word_size == 16) {
        pack_next<Ts, Td, __uint128_t>(src, dest, n, size);
    }
}

template <typename Ts, typename Td, typename Tw>
inline void
unpack_next(std::vector<Ts*>* src, std::vector<Td*>* dest, int n, size_t size)
{
    int i;
    std::vector<Tw*> tmp(n, nullptr);
    for (i = 0; i < n; i++) {
        tmp[i] = static_cast<Tw*>(static_cast<void*>(dest->at(i)));
        std::copy_n(src->at(i), size, tmp[i]);
    }
    tmp.shrink_to_fit();
}

/*
 * Cast buffers of destination to a type corresponding to 'word_size'
 * Copy elements from source buffers to casted buffers
 *
 *  sizeof(Ts) >= sizeof(Td)
 *
 * @param src: source vector
 * @param dest: destination vector
 * @param n: number of buffers per vector
 * @param size: number of elements per source buffer
 * @param word_size: number of bytes used to store data in each element of
 *  source buffers
 * @return
 */
template <typename Ts, typename Td>
inline void unpack(
    std::vector<Ts*>* src,
    std::vector<Td*>* dest,
    int n,
    size_t size,
    size_t word_size)
{
    assert(sizeof(Ts) >= word_size);
    assert(word_size % sizeof(Td) == 0);
    // get only word_size bytes from each element
    if (word_size == 1) {
        unpack_next<Ts, Td, uint8_t>(src, dest, n, size);
    } else if (word_size == 2) {
        unpack_next<Ts, Td, uint16_t>(src, dest, n, size);
    } else if (word_size == 4) {
        unpack_next<Ts, Td, uint32_t>(src, dest, n, size);
    } else if (word_size == 8) {
        unpack_next<Ts, Td, uint64_t>(src, dest, n, size);
    } else if (word_size == 16) {
        unpack_next<Ts, Td, __uint128_t>(src, dest, n, size);
    }
}

template <typename T>
class Vec {
  private:
    T* mem;
    int mem_len;
    bool new_mem;

  protected:
    int n;

  public:
    RN<T>* rn;
    Vec(RN<T>* rn, int n, T* mem = nullptr, int mem_len = 0);
    virtual ~Vec();
    virtual int get_n(void);
    int get_mem_len(void);
    void zero_fill(void);
    void fill(T val);
    virtual void set(int i, T val);
    virtual T get(int i);
    T* get_mem();
    void set_mem(T* mem, int mem_len);
    virtual bool is_v2vec();
    void mul_scalar(T scalar);
    void mul_beta(T beta);
    void hadamard_mul(Vec<T>* v);
    void hadamard_mul(V2Vec<T>* v);
    void add(Vec<T>* v);
    void add(V2Vec<T>* v);
    void add(Vec<T>* v, int offset);
    void add_mutual(Vec<T>* v);
    void add_mutual(Vec<T>* v, int offset);
    void add_mutual(Vec<T>* v, int offset, int len);
    void copy(Vec<T>* v);
    void copy(Vec<T>* v, int n);
    void copy(Vec<T>* v, int n, int offset);
    void copy(Vec<T>* v, int n, int dest_offset, int src_offset);
    bool eq(Vec<T>* v);
    void to_poly(Poly<T>* poly);
    virtual void dump(void);
};

template <typename T>
Vec<T>::Vec(RN<T>* rn, int n, T* mem, int mem_len)
{
    this->rn = rn;
    this->n = n;
    if (mem == nullptr) {
        this->mem = new T[n];
        this->mem_len = n;
        this->new_mem = true;
    } else {
        this->mem = mem;
        this->mem_len = mem_len;
        this->new_mem = false;
    }
}

template <typename T>
Vec<T>::~Vec()
{
    if (new_mem)
        delete[] this->mem;
}

template <typename T>
inline int Vec<T>::get_n(void)
{
    return this->n;
}

template <typename T>
inline int Vec<T>::get_mem_len(void)
{
    return this->mem_len;
}

template <typename T>
void Vec<T>::zero_fill(void)
{
    int i;

    for (i = 0; i < n; i++)
        set(i, 0);
}

template <typename T>
void Vec<T>::fill(T val)
{
    int i;

    for (i = 0; i < n; i++)
        mem[i] = val;
}

template <typename T>
inline void Vec<T>::set(int i, T val)
{
    assert(i >= 0 && i < n);

    mem[i] = val;
}

template <typename T>
inline T Vec<T>::get(int i)
{
    assert(i >= 0 && i < n);

    return mem[i];
}

template <typename T>
inline T* Vec<T>::get_mem()
{
    return mem;
}

template <typename T>
inline void Vec<T>::set_mem(T* mem, int mem_len)
{
    if (new_mem)
        delete[] this->mem;
    new_mem = false;
    this->mem = mem;
    this->mem_len = mem_len;
}

template <typename T>
bool Vec<T>::is_v2vec()
{
    return false;
}

/**
 * Multiplication of a vector by a scalar
 *
 * @param scalar
 */
template <typename T>
void Vec<T>::mul_scalar(T scalar)
{
    for (int i = 0; i < n; i++)
        set(i, rn->mul(get(i), scalar));
}

/**
 * Multiplication of ith element of a vector by a scalar = beta^i
 *
 * @param scalar
 */
template <typename T>
void Vec<T>::mul_beta(T beta)
{
    T coef = beta;
    for (int i = 1; i < n; i++) {
        mem[i] = rn->mul(mem[i], coef);
        coef = rn->mul(coef, beta);
    }
}

/**
 * entrywise product
 *
 * @param v
 */
template <typename T>
void Vec<T>::hadamard_mul(Vec<T>* v)
{
    assert(n == v->get_n());
    T* src = v->get_mem();
    for (int i = 0; i < n; i++)
        mem[i] = rn->mul(mem[i], src[i]);
}

/**
 * entrywise product
 *
 * @param v
 */
template <typename T>
void Vec<T>::hadamard_mul(V2Vec<T>* v)
{
    assert(n == v->get_n());
    T* src = v->get_mem();
    int i;
    int j;
    for (i = 0; i < n / 2; i++)
        mem[i] = rn->mul(mem[i], src[i]);
    for (j = 0; i < n; i++, j++)
        mem[i] = rn->mul(mem[i], src[j]);
}

template <>
void Vec<uint32_t>::hadamard_mul(V2Vec<uint32_t>* v);
template <>
void Vec<uint64_t>::hadamard_mul(V2Vec<uint64_t>* v);

template <typename T>
void Vec<T>::add(Vec<T>* v)
{
    assert(n == v->get_n());

    T* src = v->get_mem();
    for (int i = 0; i < n; i++)
        mem[i] = rn->add(mem[i], src[i]);
}

template <typename T>
void Vec<T>::add(V2Vec<T>* v)
{
    assert(n == v->get_n());
    T* src = v->get_mem();
    int i;
    int j;
    for (i = 0; i < n / 2; i++)
        mem[i] = rn->add(mem[i], src[i]);
    for (j = 0; i < n; i++, j++)
        mem[i] = rn->add(mem[i], src[j]);
}

template <typename T>
void Vec<T>::add(Vec<T>* v, int offset)
{
    assert(n >= v->get_n() + offset);

    T* src = v->get_mem();
    T* dest = mem + offset;

    for (int i = 0; i < v->get_n(); i++)
        dest[i] = rn->add(dest[i], src[i]);
}

template <typename T>
void Vec<T>::add_mutual(Vec<T>* v)
{
    int len = v->get_n();
    assert(n >= len);
    T* src = v->get_mem();
    T* dest = this->mem;
    for (int i = 0; i < len; i++)
        dest[i] = rn->add(dest[i], src[i]);
}

template <typename T>
void Vec<T>::add_mutual(Vec<T>* v, int offset)
{
    int len = v->get_n();
    assert(len == 0 || n - offset >= len);
    T* src = v->get_mem();
    T* dest = this->mem + offset;
    for (int i = 0; i < len; i++)
        dest[i] = rn->add(dest[i], src[i]);
}

template <typename T>
void Vec<T>::add_mutual(Vec<T>* v, int offset, int len)
{
    assert(len == 0 || n - offset >= len);
    assert(v->get_n() >= len);
    T* src = v->get_mem();
    T* dest = this->mem + offset;
    for (int i = 0; i < len; i++)
        dest[i] = rn->add(dest[i], src[i]);
}

template <>
void Vec<uint32_t>::add(V2Vec<uint32_t>* v);
template <>
void Vec<uint64_t>::add(V2Vec<uint64_t>* v);

template <typename T>
void Vec<T>::copy(Vec<T>* v)
{
    assert(v->get_mem_len() <= this->mem_len);
    std::copy_n(v->get_mem(), v->get_mem_len(), this->mem);
}

template <typename T>
void Vec<T>::copy(Vec<T>* v, int n)
{
    assert(n <= this->mem_len);
    int v_mem_len = v->get_mem_len();
    if (v_mem_len >= n)
        std::copy_n(v->get_mem(), n, this->mem);
    else {
        std::copy_n(v->get_mem(), v_mem_len, this->mem);
        std::memset(this->mem + v_mem_len, 0, sizeof(T) * (n - v_mem_len));
    }
}

template <typename T>
void Vec<T>::copy(Vec<T>* v, int n, int offset)
{
    assert(n + offset <= this->mem_len);
    int v_mem_len = v->get_mem_len();
    T* dest = this->mem + offset;
    if (v_mem_len >= n)
        std::copy_n(v->get_mem(), n, dest);
    else {
        std::copy_n(v->get_mem(), v_mem_len, dest);
        std::memset(dest + v_mem_len, 0, sizeof(T) * (n - v_mem_len));
    }
}

template <typename T>
void Vec<T>::copy(Vec<T>* v, int n, int dest_offset, int src_offset)
{
    assert(n + dest_offset <= this->mem_len);
    int src_len = v->get_mem_len() - src_offset;
    T* dest = this->mem + dest_offset;
    T* src = v->get_mem() + src_offset;
    if (src_len >= n)
        std::copy_n(src, n, dest);
    else {
        std::copy_n(src, src_len, dest);
        std::memset(dest + src_len, 0, sizeof(T) * (n - src_len));
    }
}

template <typename T>
bool Vec<T>::eq(Vec<T>* v)
{
    if (v->get_n() != this->n)
        return false;

    for (int i = 0; i < v->get_n(); i++) {
        if (get(i) != v->get(i))
            return false;
    }

    return true;
}

template <typename T>
void Vec<T>::to_poly(Poly<T>* poly)
{
    poly->clear();
    for (int i = 0; i < this->n; i++) {
        poly->set(i, get(i));
    }
}

template <typename T>
void Vec<T>::dump(void)
{
    std::cout << "( ";
    for (int i = 0; i < n; i++)
        std::cout << get(i) << " ";
    std::cout << ")\n";
}

} // namespace vec
} // namespace nttec

#endif
