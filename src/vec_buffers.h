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
#ifndef __NTTEC_VEC_BUFFERS_H__
#define __NTTEC_VEC_BUFFERS_H__

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

#include "core.h"

namespace nttec {
namespace vec {

/* Available cases of allocating memory
 * NO_ALLOC:    do not allocate any memory
 * SLICE_ALLOC: allocate only a vector of pointers each points to an allocated
 *              memory
 * ZERO_EXTEND: allocate a zero buffers and push back to the vector memory
 * COMBINED:    allocate only a vector of pointers that are combined from
 *              two vectors of two input buffers vector
 * FULL_ALLOC:  fully allocate memory including a vector of pointers each for a
 *              memory of `size` elements
 */
enum MEM_ALLOC_CASES {
    NO_ALLOC = 0,
    SLICE_ALLOC,
    ZERO_EXTEND,
    COMBINED,
    FULL_ALLOC,
};

/** A vector of `n` buffers (array of T).
 *
 * A Buffers contains pointers to `n` buffers such as
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
 * In memory, the Buffers looks like:
 *
 * v[0]  | v[1]   | … | v[n-1]
 * ----- | ------ | - | --------
 * &buf1 |  &buf2 | … | &bufN
 *
 * A Buffers can:
 * - owns all the memory (the vector of buffers and the buffers themselve).
 * - own the vector and use existing buffers (only the vector is allocated).
 * - nothing (a shallow copy of another Buffers).
 */
template <typename T>
class Buffers {
  public:
    Buffers(int n, size_t size, std::vector<T*>* mem = nullptr);
    Buffers(const Buffers<T>& vec, int n = 0);
    Buffers(const Buffers<T>& vec, int begin, int end);
    Buffers(const Buffers<T>& vec1, const Buffers<T>& vec2);
    virtual ~Buffers();
    virtual int get_n(void) const;
    virtual size_t get_size(void) const;
    int get_mem_len(void);
    void zero_fill(void);
    void fill(int i, T value);
    virtual void set(int i, T* buf);
    virtual T* get(int i) const;
    std::vector<T*>* get_mem() const;
    void set_mem(std::vector<T*>* mem);
    void copy(const Buffers<T>& v);
    void copy(int i, T* buf);
    void separate_even_odd();
    void separate_even_odd(const Buffers<T>& even, const Buffers<T>& odd);
    bool eq(const Buffers<T>& v);
    virtual void dump(void);

  protected:
    std::vector<T*>* mem = nullptr;
    int mem_len;
    size_t size;
    int n;

  private:
    int mem_alloc_case = FULL_ALLOC;
    T* zeros = nullptr;
};

/**
 * Constructor of Buffers
 *
 * @param n - size of vector of pointers, stored in `mem`
 * @param size - number of elements of each memory pointed by a pointer of `mem`
 * @param mem - null pointer or a pointer to a vector of buffers
 */
template <typename T>
Buffers<T>::Buffers(int n, size_t size, std::vector<T*>* mem)
{
    this->n = n;
    this->size = size;
    this->mem_len = n * size;
    if (mem == nullptr) {
        this->mem_alloc_case = FULL_ALLOC;
        this->mem = new std::vector<T*>(n, nullptr);
        size_t total_size = n * size;
        this->mem->at(0) = aligned_allocate<T>(total_size);
        for (int i = 1; i < n; i++) {
            this->mem->at(i) = this->mem->at(i - 1) + size;
        }
    } else {
        this->mem_alloc_case = NO_ALLOC;
        this->mem = mem;
    }
}

/**
 * Constructor of Buffers by cloning partly from a given vector
 *
 * @param vec - a given Buffers instance
 * @param n - number of buffers pointed by the constructed vector
 */
template <typename T>
Buffers<T>::Buffers(const Buffers<T>& vec, int n)
{
    assert(n >= 0);
    int i;
    int vec_n = vec.get_n();

    this->size = vec.get_size();
    this->n = (n == 0) ? vec_n : n;
    this->mem = new std::vector<T*>(this->n, nullptr);
    this->mem_len = n * size;

    this->mem_alloc_case = FULL_ALLOC;

    size_t total_size = n * size;
    this->mem->at(0) = aligned_allocate<T>(total_size);
    for (int i = 1; i < n; i++) {
        this->mem->at(i) = this->mem->at(i - 1) + size;
    }

    int copy_len = (this->n <= vec_n) ? this->n : vec_n;
    for (i = 0; i < copy_len; i++) {
        std::copy_n(vec.get(i), this->size, this->mem->at(i));
    }

    if (this->n > vec_n) { // padding zeros
        for (i = vec_n; i < this->n; i++) {
            std::memset(this->mem->at(i), 0, this->size * sizeof(T));
        }
    }
}

/**
 * Constructor of Buffers by slicing from a given vector. Buffers are sliced
 * from `begin` (inclusive) to `end` (exclusive). If `end` is out of the input
 * vector, zero buffers are padded.
 *
 * @param vec - a given Buffers instance
 * @param begin - index of buffer pointed by vec, that is the 1st buffer
 *                pointed by the output Buffers
 * @param end - index of buffer pointed by vec, that limits buffers
 *                pointed by the output Buffers
 */
template <typename T>
Buffers<T>::Buffers(const Buffers<T>& vec, int begin, int end)
{
    assert(begin >= 0 && begin < end);

    this->n = end - begin;
    this->size = vec.get_size();
    this->mem_len = this->n * this->size;
    std::vector<T*>* vec_mem = vec.get_mem();

    // slice from input buffers
    if (end <= vec.get_n()) {
        this->mem_alloc_case = SLICE_ALLOC;

        this->mem = new std::vector<T*>(
            vec_mem->begin() + begin, vec_mem->begin() + end);
    } else { // slice and padding zeros
        this->mem_alloc_case = ZERO_EXTEND;

        this->zeros = aligned_allocate<T>(this->size);
        std::memset(this->zeros, 0, this->size * sizeof(T));

        this->mem =
            new std::vector<T*>(vec_mem->begin() + begin, vec_mem->end());
        this->mem->insert(this->mem->end(), end - vec.get_n(), this->zeros);
    }
}

template <typename T>
Buffers<T>::Buffers(const Buffers<T>& vec1, const Buffers<T>& vec2)
{
    assert(vec1.get_size() == vec2.get_size());

    int n1 = vec1.get_n();
    int n2 = vec2.get_n();

    this->n = n1 + n2;
    this->size = vec1.get_size();
    this->mem_len = this->n * this->size;

    this->mem_alloc_case = COMBINED;

    this->mem = new std::vector<T*>();
    this->mem->insert(
        this->mem->end(), vec1.get_mem()->begin(), vec1.get_mem()->end());
    this->mem->insert(
        this->mem->end(), vec2.get_mem()->begin(), vec2.get_mem()->end());
}

template <typename T>
Buffers<T>::~Buffers()
{
    if (this->mem_alloc_case != NO_ALLOC && this->mem != nullptr) {
        if (this->mem_alloc_case == FULL_ALLOC) {
            aligned_deallocate<T>(this->mem->at(0));
        } else if (this->mem_alloc_case == ZERO_EXTEND) {
            aligned_deallocate<T>(this->zeros);
        }
        this->mem->shrink_to_fit();
        delete this->mem;
    }
}

template <typename T>
inline int Buffers<T>::get_n(void) const
{
    return n;
}

template <typename T>
inline size_t Buffers<T>::get_size(void) const
{
    return size;
}

template <typename T>
inline int Buffers<T>::get_mem_len(void)
{
    return mem_len;
}

template <typename T>
void Buffers<T>::zero_fill(void)
{
    for (int i = 0; i < n; i++)
        std::memset(mem->at(i), 0, size * sizeof(T));
}

template <typename T>
void Buffers<T>::fill(int i, T value)
{
    std::memset(mem->at(i), value, size * sizeof(T));
}

template <typename T>
inline void Buffers<T>::set(int i, T* buf)
{
    assert(i >= 0 && i < n);

    if ((mem_alloc_case == 0) && (this->mem->at(i) != nullptr))
        aligned_deallocate<T>(this->mem->at(i));

    this->mem->at(i) = buf;
}

template <typename T>
inline T* Buffers<T>::get(int i) const
{
    assert(i >= 0 && i < n);
    return this->mem->at(i);
}

template <typename T>
inline std::vector<T*>* Buffers<T>::get_mem() const
{
    return mem;
}

template <typename T>
inline void Buffers<T>::set_mem(std::vector<T*>* mem)
{
    this->mem = mem;
}

template <typename T>
void Buffers<T>::copy(const Buffers<T>& v)
{
    assert(v.get_n() == n);
    assert(v.get_size() <= size);
    size_t v_size = v.get_size();
    for (int i = 0; i < n; i++)
        std::copy_n(v.get(i), v_size, this->mem->at(i));
}

template <typename T>
void Buffers<T>::copy(int i, T* buf)
{
    std::copy_n(buf, this->size, this->mem->at(i));
}

template <typename T>
void Buffers<T>::separate_even_odd()
{
    std::vector<T*> _mem(n, nullptr);
    int half = n / 2;
    int j = 0;
    int i;
    for (i = 0; i < n; i += 2) {
        _mem[j] = get(i);            // even
        _mem[j + half] = get(i + 1); // odd
        j++;
    }
    for (i = 0; i < n; i++) {
        mem->at(i) = _mem[i];
    }
    _mem.shrink_to_fit();
}

template <typename T>
void Buffers<T>::separate_even_odd(
    const Buffers<T>& even,
    const Buffers<T>& odd)
{
    int j = 0;
    int i;
    std::vector<T*>* even_mem = even.get_mem();
    std::vector<T*>* odd_mem = odd.get_mem();
    for (i = 0; i < n; i += 2) {
        even_mem->at(j) = get(i);    // even
        odd_mem->at(j) = get(i + 1); // odd
        j++;
    }
}

template <typename T>
bool Buffers<T>::eq(const Buffers<T>& v)
{
    if ((v.get_n() != n) || (v.get_size() != size))
        return false;

    for (int i = 0; i < n; i++) {
        T* a = get(i);
        T* b = v.get(i);
        for (size_t j = 0; j < size; j++) {
            if (a[j] != b[j])
                return false;
        }
    }

    return true;
}

template <typename T>
void Buffers<T>::dump(void)
{
    for (int i = 0; i < n; i++) {
        std::cout << "\n\t" << i << ": ";
        for (size_t j = 0; j < size - 1; j++) {
            std::cout << (get(i))[j] << "-";
        }
        if (size > 0) {
            std::cout << (get(i))[size - 1];
        }
    }
    std::cout << "\n)\n";
}

} // namespace vec
} // namespace nttec

#endif
