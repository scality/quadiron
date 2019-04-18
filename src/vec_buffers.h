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

/** A vector of `n` DATA buffers (array of T) and an OPTIONAL vector of `n` META
 * buffers (array of `uint8_t`). A meta buffer corresponds to a data buffer.
 *Every `s`-byte data element of a data buffer corresponds to a `s`-bit element
 *of the meta buffer, where `s = sizeof(T)`. A pair of data and meta elements
 *can represent an integer of `8*s + s`-bits
 *
 * Each data element points to a buffer of size `m * sizeof(T)` bytes.
 * Hence, each meta element points to a buffer of size
 * `bsize = m * sizeof(T) / 8` bytes that should be an integer, i.e.
 * m * sizeof(T) % 8 = 0. This condition is not difficult to achieved as a
 * convenient `size` can be always chosen with a negligible impact on its
 * application.
 *
 * Example:
 *
 * We have two sets each of `N` independent buffers:
 *
 * data1    | data2    | … | dataN
 * -------- | -------- | - | -------
 * data1[0] | data2[0] | … | dataN[0]
 * data1[1] | data2[1] | … | dataN[1]
 * …        | …        | … | …
 *
 *
 * meta1    | meta2    | … | metaN
 * -------- | -------- | - | -------
 * meta1[0] | meta2[0] | … | metaN[0]
 * meta1[1] | meta2[1] | … | metaN[1]
 * …        | …        | … | …
 *
 * In memory, the Buffers looks like:
 *
 * u[0]   | u[1]    | … | u[n-1]
 * ------ | ------- | - | --------
 * &data1 |  &data2 | … | &dataN
 *
 * v[0]   | v[1]    | … | v[n-1]
 * ------ | ------- | - | --------
 * &meta1 |  &meta2 | … | &metaN
 *
 * A Buffers can:
 * - owns all the memory (the vector of buffers and the buffers themselve).
 * - own the vector and use existing buffers (only the vector is allocated).
 * - nothing (a shallow copy of another Buffers).
 */
template <typename T>
class Buffers final {
  public:
    Buffers(int n, size_t size, bool has_meta = false);
    Buffers(
        int n,
        size_t size,
        const std::vector<T*>& mem,
        const std::vector<uint8_t*>* meta = nullptr);
    Buffers(const Buffers<T>& vec, int n = 0);
    Buffers(const Buffers<T>& vec, int begin, int end);
    Buffers(const Buffers<T>& vec1, const Buffers<T>& vec2);
    Buffers(const Buffers<T>& vec, const Vector<T>& map, unsigned n);
    ~Buffers();
    bool has_meta() const;
    int get_n(void) const;
    size_t get_size(void) const;
    size_t get_meta_size(void) const;
    size_t get_mem_len(void);
    void zero_fill(void);
    void fill(int i, T value);
    T* get(int i);
    const T* get(int i) const;
    uint8_t* get_meta(int i);
    const uint8_t* get_meta(int i) const;
    const std::vector<T*>& get_mem() const;
    const std::vector<uint8_t*>& get_meta() const;
    void set_mem(std::vector<T*>* mem);
    void copy(const Buffers<T>& v);
    void copy(const Buffers<T>& v, int src_id, int dest_id);
    void radix2_fft_prepare(
        const Buffers<T>& input,
        const std::vector<T>& scrambler,
        unsigned data_len,
        unsigned group_len);
    void radix2_fft_inv_prepare(const Buffers<T>& input);
    friend bool operator==<T>(const Buffers<T>& lhs, const Buffers<T>& rhs);
    void dump(void);
    void swap(unsigned i, unsigned j);

  protected:
    std::vector<T*> mem;
    std::vector<uint8_t*> meta;
    size_t mem_len;
    size_t size;
    size_t meta_size = 0;
    int n;

  private:
    simd::AlignedAllocator<T> allocator;
    simd::AlignedAllocator<uint8_t> allocator_meta;
    BufMemAlloc mem_alloc_case = BufMemAlloc::FULL;
    T* zeros = nullptr;
    uint8_t* zeros_meta = nullptr;
    bool m_meta = false;

    void compute_meta_size();
    void allocate_meta(bool init_zero = false);
};

/**
 * Constructor of Buffers
 *
 * @param n - size of vector of pointers, stored in `mem`
 * @param size - number of elements of each memory pointed by a pointer of `mem`
 */
template <typename T>
Buffers<T>::Buffers(int n, size_t size, bool has_meta)
{
    this->n = n;
    this->size = size;
    this->mem_len = n * size;

    this->mem_alloc_case = BufMemAlloc::FULL;
    mem.reserve(n);
    for (int i = 0; i < n; i++) {
        mem.push_back(this->allocator.allocate(size));
    }

    if (has_meta) {
        this->m_meta = has_meta;
        this->compute_meta_size();
        this->allocate_meta(true);
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
Buffers<T>::Buffers(
    int n,
    size_t size,
    const std::vector<T*>& mem,
    const std::vector<uint8_t*>* meta)
{
    this->n = n;
    this->size = size;
    this->mem_len = n * size;
    this->mem_alloc_case = BufMemAlloc::NONE;
    this->mem = mem;

    if (meta) {
        this->m_meta = true;
        this->compute_meta_size();
        this->meta = *meta;
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

    this->m_meta = vec.has_meta();
    if (this->m_meta) {
        this->compute_meta_size();
        this->allocate_meta();

        for (i = 0; i < copy_len; i++) {
            std::copy_n(vec.get_meta(i), meta_size, meta[i]);
        }

        if (this->n > vec_n) { // padding zeros
            for (i = vec_n; i < this->n; i++) {
                std::fill_n(meta[i], meta_size, 0);
            }
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

    this->m_meta = vec.has_meta();
    if (this->m_meta) {
        const std::vector<uint8_t*> vec_meta = vec.get_meta();
        this->compute_meta_size();
        meta.reserve(this->n);
        // slice from input buffers
        if (end <= vec.get_n()) {
            meta.insert(
                meta.end(), vec_meta.begin() + begin, vec_meta.begin() + end);
        } else { // slice and padding zeros
            this->zeros_meta = this->allocator_meta.allocate(meta_size);
            std::fill_n(this->zeros_meta, meta_size, 0);

            meta.insert(meta.end(), vec_meta.begin() + begin, vec_meta.end());
            meta.insert(meta.end(), end - vec.get_n(), zeros_meta);
        }
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
    assert(vec1.has_meta() == vec2.has_meta());

    this->n = vec1.get_n() + vec2.get_n();
    this->size = vec1.get_size();
    this->mem_len = this->n * this->size;

    this->mem_alloc_case = BufMemAlloc::COMBINED;

    mem.reserve(this->n);
    mem.insert(mem.end(), vec1.get_mem().begin(), vec1.get_mem().end());
    mem.insert(mem.end(), vec2.get_mem().begin(), vec2.get_mem().end());

    this->m_meta = vec1.has_meta();
    if (this->m_meta) {
        this->compute_meta_size();
        meta.reserve(this->n);
        meta.insert(meta.end(), vec1.get_meta().begin(), vec1.get_meta().end());
        meta.insert(meta.end(), vec2.get_meta().begin(), vec2.get_meta().end());
    }
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

    this->m_meta = vec.has_meta();
    if (this->m_meta) {
        this->compute_meta_size();

        const std::vector<uint8_t*> vec_meta = vec.get_meta();
        // output is sliced & shuffled from `vec`
        meta.reserve(this->n);
        if (vec_n < n) { // output is zero-extended & shuffled from `vec`
            this->zeros_meta = this->allocator_meta.allocate(meta_size);
            std::fill_n(this->zeros_meta, meta_size, 0);

            for (unsigned i = 0; i < n; ++i) {
                meta.push_back(zeros_meta);
            }
        }
        for (unsigned i = 0; i < map_len; ++i) {
            meta[map.get(i)] = vec_meta[i];
        }
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
            if (this->m_meta) {
                for (int i = 0; i < n; i++) {
                    this->allocator_meta.deallocate(meta[i], meta_size);
                }
            }
        } else if (this->mem_alloc_case == BufMemAlloc::ZERO_EXTEND) {
            this->allocator.deallocate(this->zeros, size);
            if (this->m_meta) {
                this->allocator_meta.deallocate(this->zeros_meta, meta_size);
            }
        }
    }
}

template <typename T>
inline void Buffers<T>::compute_meta_size()
{
    const size_t bytes = size * sizeof(T);
    meta_size = bytes / CHAR_BIT;
    while (meta_size * CHAR_BIT < bytes) {
        meta_size++;
    }
}

template <typename T>
inline void Buffers<T>::allocate_meta(bool init_zero)
{
    meta.reserve(n);
    for (int i = 0; i < n; i++) {
        meta.push_back(this->allocator_meta.allocate(meta_size));
    }

    if (init_zero) {
        for (int i = 0; i < n; i++) {
            std::fill_n(meta[i], meta_size, 0);
        }
    }
}

template <typename T>
inline bool Buffers<T>::has_meta(void) const
{
    return m_meta;
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
inline size_t Buffers<T>::get_meta_size(void) const
{
    return meta_size;
}

template <typename T>
inline size_t Buffers<T>::get_mem_len(void)
{
    return mem_len;
}

template <typename T>
void Buffers<T>::zero_fill(void)
{
    for (int i = 0; i < n; i++) {
        std::fill_n(mem[i], size, 0);
    }
    reset_meta();
}

template <typename T>
void Buffers<T>::fill(int i, T value)
{
    std::fill_n(mem[i], size, value);
    if (m_meta) {
        std::fill_n(meta[i], meta_size, 0);
    }
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
inline uint8_t* Buffers<T>::get_meta(int i)
{
    assert(i >= 0 && i < n);
    return meta[i];
}

template <typename T>
inline const uint8_t* Buffers<T>::get_meta(int i) const
{
    assert(i >= 0 && i < n);
    return meta[i];
}

template <typename T>
inline const std::vector<T*>& Buffers<T>::get_mem() const
{
    return mem;
}

template <typename T>
inline const std::vector<uint8_t*>& Buffers<T>::get_meta() const
{
    return meta;
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
    for (int i = 0; i < n; ++i) {
        std::copy_n(v.get(i), v_size, mem[i]);
    }
    if (v.has_meta()) {
        meta_size = v.get_meta_size();
        for (int i = 0; i < n; ++i) {
            std::copy_n(v.get_meta(i), meta_size, meta[i]);
        }
    }
}

template <typename T>
void Buffers<T>::copy(const Buffers<T>& v, int src_id, int dest_id)
{
    assert(m_meta == v.has_meta());

    std::copy_n(v.get(src_id), size, mem[dest_id]);

    if (m_meta) {
        std::copy_n(v.get_meta(src_id), meta_size, meta[dest_id]);
    }
}

template <typename T>
bool operator==(const Buffers<T>& lhs, const Buffers<T>& rhs)
{
    if (lhs.n != rhs.n || lhs.get_size() != rhs.get_size()
        || lhs.get_meta_size() != rhs.get_meta_size()
        || lhs.has_meta() != rhs.has_meta()) {
        return false;
    }
    for (int i = 0; i < lhs.n; i++) {
        const T* lhs_vec = lhs.get(i);
        const T* rhs_vec = rhs.get(i);

        for (size_t j = 0; j < lhs.get_size(); j++) {
            if (lhs_vec[j] != rhs_vec[j]) {
                return false;
            }
        }
    }

    if (lhs.has_meta()) {
        for (int i = 0; i < lhs.n; i++) {
            const uint8_t* lhs_meta = lhs.get_meta(i);
            const uint8_t* rhs_meta = rhs.get_meta(i);

            for (size_t j = 0; j < lhs.get_meta_size(); j++) {
                if (lhs_meta[j] != rhs_meta[j]) {
                    return false;
                }
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
    if (m_meta) {
        swap(meta[i], meta[j]);
    }
}

/** Perform a preparatiion in radix-2 FFT algorithm
 *
 * @param input - input buffers
 * @param scrambler - a vector stores bit-reversed values
 * @param data_len - a conven length
 * @param group_len - number of elements in which a same operations will be
 * performed
 */
template <typename T>
void Buffers<T>::radix2_fft_prepare(
    const Buffers<T>& input,
    const std::vector<T>& scrambler,
    unsigned data_len,
    unsigned group_len)
{
    const unsigned input_len = input.get_n();
    const std::vector<T*>& i_mem = input.get_mem();
    const std::vector<uint8_t*>& i_meta = input.get_meta();

    for (unsigned idx = 0; idx < input_len; ++idx) {
        // set output  = scramble(input), i.e. bit reversal ordering
        for (unsigned i = scrambler[idx]; i < scrambler[idx] + group_len; ++i) {
            std::copy_n(i_mem[idx], size, mem[i]);
            m_meta ? std::copy_n(i_meta[idx], meta_size, meta[i]) : nullptr;
        }
    }
    for (unsigned idx = input_len; idx < data_len; ++idx) {
        // set output  = scramble(input), i.e. bit reversal ordering
        for (unsigned i = scrambler[idx]; i < scrambler[idx] + group_len; ++i) {
            std::fill_n(mem[i], size, 0);
            m_meta ? std::fill_n(meta[i], meta_size, 0) : nullptr;
        }
    }
}

/** Perform a preparatiion in radix-2 inverse FFT algorithm
 *
 * @param input - input buffers
 * performed
 */
template <typename T>
void Buffers<T>::radix2_fft_inv_prepare(const Buffers<T>& input)
{
    const unsigned len = this->n;
    const unsigned input_len = input.get_n();
    const std::vector<T*>& i_mem = input.get_mem();
    const std::vector<uint8_t*>& i_meta = input.get_meta();

    unsigned i;
    for (i = 0; i < input_len; ++i) {
        std::copy_n(i_mem[i], size, mem[i]);
        m_meta ? std::copy_n(i_meta[i], meta_size, meta[i]) : nullptr;
    }

    if (input_len < len) {
        const unsigned input_len_power_2 = arith::ceil2<unsigned>(input_len);
        for (; i < input_len_power_2; ++i) {
            std::fill_n(mem[i], size, 0);
            m_meta ? std::fill_n(meta[i], meta_size, 0) : nullptr;
        }
    }
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

    if (m_meta) {
        std::cout << "\nMeta:\n";

        for (int i = 0; i < n; i++) {
            std::cout << "\n\t" << i << ": ";
            for (size_t j = 0; j < meta_size - 1; j++) {
                std::cout << static_cast<unsigned>((get_meta(i))[j]) << "-";
            }
            if (size > 0) {
                std::cout << static_cast<unsigned>(
                    (get_meta(i))[meta_size - 1]);
            }
        }
    }

    std::cout << "\n\n";
}

} // namespace vec
} // namespace quadiron

#endif
