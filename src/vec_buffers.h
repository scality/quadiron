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
#ifndef __QUAD_VEC_BUFFERS_H__
#define __QUAD_VEC_BUFFERS_H__

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

#include "core.h"
#include "simd/simd.h"

namespace quadiron {
namespace vec {

template <typename T>
class Vector;

// Forward declarations.
template <typename>
class Buffers;
template <typename T>
bool operator==(const Buffers<T>& lhs, const Buffers<T>& rhs);

/// Available cases of allocating memory
enum class BufMemAlloc {
    /// Do not allocate any memory
    NONE = 0,

    /// Allocate only a vector of pointers each points to an allocated memory
    SLICE,

    /// Allocate a zero buffers and push back to the vector memory
    ZERO_EXTEND,

    /// Allocate only a vector of pointers that are from two vectors of two
    /// input buffers vector
    COMBINED,

    /// Fully allocate memory including a vector of pointers each for a
    /// memory of `size` elements
    FULL,
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
class Buffers final {
  public:
    Buffers(int n, size_t size);
    Buffers(int n, size_t size, const std::vector<T*>& mem);
    Buffers(const Buffers<T>& vec, int n = 0);
    Buffers(const Buffers<T>& vec, int begin, int end);
    Buffers(const Buffers<T>& vec1, const Buffers<T>& vec2);
    Buffers(const Buffers<T>& vec, const Vector<T>& map, unsigned n);
    ~Buffers();
    int get_n(void) const;
    size_t get_size(void) const;
    size_t get_mem_len(void);
    void zero_fill(void);
    void fill(int i, T value);
    T* get(int i);
    const T* get(int i) const;
    const std::vector<T*>& get_mem() const;
    void set_mem(std::vector<T*>* mem);
    void copy(const Buffers<T>& v);
    void copy(int i, T* buf);
    friend bool operator==<T>(const Buffers<T>& lhs, const Buffers<T>& rhs);
    void dump(void);
    void swap(unsigned i, unsigned j);

  protected:
    std::vector<T*> mem;
    size_t mem_len;
    size_t size;
    int n;

  private:
    simd::AlignedAllocator<T> allocator;
    BufMemAlloc mem_alloc_case = BufMemAlloc::FULL;
    T* zeros = nullptr;
};

/**
 * Constructor of Buffers
 *
 * @param n - size of vector of pointers, stored in `mem`
 * @param size - number of elements of each memory pointed by a pointer of `mem`
 */
template <typename T>
Buffers<T>::Buffers(int n, size_t size)
{
    this->n = n;
    this->size = size;
    this->mem_len = n * size;

    this->mem_alloc_case = BufMemAlloc::FULL;
    mem.reserve(n);
    for (int i = 0; i < n; i++) {
        mem.push_back(this->allocator.allocate(size));
    }
}

/**
 * Constructor of Buffers
 *
 * @param n - size of vector of pointers, stored in `mem`
 * @param size - number of elements of each memory pointed by a pointer of `mem`
 * @param mem - a vector of buffers
 */
template <typename T>
Buffers<T>::Buffers(int n, size_t size, const std::vector<T*>& mem)
{
    this->n = n;
    this->size = size;
    this->mem_len = n * size;
    this->mem_alloc_case = BufMemAlloc::NONE;
    this->mem = mem;
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
    this->mem_len = n * size;

    this->mem_alloc_case = BufMemAlloc::FULL;

    mem.reserve(n);
    for (i = 0; i < this->n; i++) {
        mem.push_back(this->allocator.allocate(size));
    }

    int copy_len = (this->n <= vec_n) ? this->n : vec_n;
    for (i = 0; i < copy_len; i++) {
        std::copy_n(vec.get(i), this->size, mem[i]);
    }

    if (this->n > vec_n) { // padding zeros
        for (i = vec_n; i < this->n; i++) {
            std::memset(mem[i], 0, this->size * sizeof(T));
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
    const std::vector<T*> vec_mem = vec.get_mem();

    mem.reserve(this->n);
    // slice from input buffers
    if (end <= vec.get_n()) {
        this->mem_alloc_case = BufMemAlloc::SLICE;

        mem.insert(mem.end(), vec_mem.begin() + begin, vec_mem.begin() + end);
    } else { // slice and padding zeros
        this->mem_alloc_case = BufMemAlloc::ZERO_EXTEND;

        this->zeros = this->allocator.allocate(size);
        std::memset(this->zeros, 0, this->size * sizeof(T));

        mem.insert(mem.end(), vec_mem.begin() + begin, vec_mem.end());
        mem.insert(mem.end(), end - vec.get_n(), zeros);
    }
}

/**
 * Constructor of Buffers by concatinating two buffers.
 *
 * @param vec1 - a given Buffers instance
 * @param vec2 - a given Buffers instance
 */
template <typename T>
Buffers<T>::Buffers(const Buffers<T>& vec1, const Buffers<T>& vec2)
{
    assert(vec1.get_size() == vec2.get_size());

    int n1 = vec1.get_n();
    int n2 = vec2.get_n();

    this->n = n1 + n2;
    this->size = vec1.get_size();
    this->mem_len = this->n * this->size;

    this->mem_alloc_case = BufMemAlloc::COMBINED;

    mem.reserve(this->n);
    mem.insert(mem.end(), vec1.get_mem().begin(), vec1.get_mem().end());
    mem.insert(mem.end(), vec2.get_mem().begin(), vec2.get_mem().end());
}

/**
 * Constructor of Buffers whose elements are shuffled from a given Buffers.
 *
 * An extending is necessary if the output length is greater than that given
 * Buffers' length.
 *
 * The map will map every element of input vector to output:
 *  map[i] = j => output[j] = input[i]
 *
 * Ex1: input [p0, p1, p2, p3, p4] and map [3, 2, 0, 1]
 * Output: [p2, p3, p1, p0]
 *
 * Ex2: input [p0, p1, p2, p3, p4] and map [3, 6, 1, 5, 4]
 * Output: [0, p2, 0, p0, p4, p3, p1] where `0` is an all-zero buffer.
 *
 * @param vec - a given Buffers instance of `m` elements
 * @param map - a vector of `m` elements
 * @param n - output vector length
 */
template <typename T>
Buffers<T>::Buffers(
    const Buffers<T>& vec,
    const vec::Vector<T>& map,
    unsigned n)
{
    const unsigned map_len = map.get_n();
    const unsigned vec_n = vec.get_n();
    assert(map_len <= n);

    this->n = n;
    this->size = vec.get_size();
    this->mem_len = this->n * this->size;

    const std::vector<T*> vec_mem = vec.get_mem();
    // output is sliced & shuffled from `vec`
    mem.reserve(this->n);
    if (vec_n >= n) {
        this->mem_alloc_case = BufMemAlloc::SLICE;
    } else { // output is zero-extended & shuffled from `vec`
        this->mem_alloc_case = BufMemAlloc::ZERO_EXTEND;

        this->zeros = this->allocator.allocate(size);
        std::memset(this->zeros, 0, this->size * sizeof(T));

        for (unsigned i = 0; i < n; ++i) {
            mem.push_back(zeros);
        }
    }
    for (unsigned i = 0; i < map_len; ++i) {
        mem[map.get(i)] = vec_mem[i];
    }
}

template <typename T>
Buffers<T>::~Buffers()
{
    if (this->mem_alloc_case != BufMemAlloc::NONE && mem.size() > 0) {
        if (this->mem_alloc_case == BufMemAlloc::FULL) {
            for (int i = 0; i < n; i++) {
                this->allocator.deallocate(mem[i], size);
            }
        } else if (this->mem_alloc_case == BufMemAlloc::ZERO_EXTEND) {
            this->allocator.deallocate(this->zeros, size);
        }
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
inline size_t Buffers<T>::get_mem_len(void)
{
    return mem_len;
}

template <typename T>
void Buffers<T>::zero_fill(void)
{
    for (int i = 0; i < n; i++)
        std::memset(mem[i], 0, size * sizeof(T));
}

template <typename T>
void Buffers<T>::fill(int i, T value)
{
    std::fill_n(mem[i], size, value);
}

template <typename T>
inline T* Buffers<T>::get(int i)
{
    assert(i >= 0 && i < n);
    return mem[i];
}

template <typename T>
inline const T* Buffers<T>::get(int i) const
{
    assert(i >= 0 && i < n);
    return mem[i];
}

template <typename T>
inline const std::vector<T*>& Buffers<T>::get_mem() const
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
        std::copy_n(v.get(i), v_size, mem[i]);
}

template <typename T>
void Buffers<T>::copy(int i, T* buf)
{
    std::copy_n(buf, this->size, mem[i]);
}

template <typename T>
bool operator==(const Buffers<T>& lhs, const Buffers<T>& rhs)
{
    if (lhs.n != rhs.n || lhs.size != rhs.size) {
        return false;
    }
    for (int i = 0; i < lhs.n; i++) {
        const T* lhs_vec = lhs.get(i);
        const T* rhs_vec = rhs.get(i);

        for (size_t j = 0; j < lhs.size; j++) {
            if (lhs_vec[j] != rhs_vec[j]) {
                return false;
            }
        }
    }
    return true;
}

template <typename T>
void Buffers<T>::swap(unsigned i, unsigned j)
{
    using std::swap;
    swap(mem[i], mem[j]);
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
} // namespace quadiron

#endif
