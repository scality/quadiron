/* -*- mode: c++ -*- */
#ifndef __NTTEC_VECP_H__
#define __NTTEC_VECP_H__

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

namespace nttec {
namespace vec {

/*
 * Wrapper of a vector of `n` buffers each of `size` elements of type T
 *  std::vector<T*> vec(n)
 *    vec[i] = T[size];
 */

template <typename T>
class Vecp {
  private:
    /* Three cases of allocating memory
     * 0: fully allocate memory including a vector of pointers each for a memory
     *    of `size` elements
     * 1: allocate only a vector of pointers each points to an allocated memory
     * 2: do not allocate any memory
     */
    int mem_alloc_case = 0;

  protected:
    std::vector<T*>* mem = nullptr;
    int mem_len;
    size_t size;
    int n;

  public:
    Vecp(int n, size_t size, std::vector<T*>* mem = nullptr);
    Vecp(Vecp<T>* vec, int n = 0);
    Vecp(Vecp<T>* vec, int begin, int end);
    virtual ~Vecp();
    virtual int get_n(void);
    virtual size_t get_size(void);
    int get_mem_len(void);
    void zero_fill(void);
    virtual void set(int i, T* buf);
    virtual T* get(int i);
    std::vector<T*>* get_mem();
    void set_mem(std::vector<T*>* mem);
    void copy(Vecp<T>* v);
    void separate_even_odd();
    void separate_even_odd(Vecp<T>* even, Vecp<T>* odd);
    bool eq(Vecp<T>* v);
    virtual void dump(void);
};

/**
 * Constructor of Vecp
 *
 * @param n - size of vector of pointers, stored in `mem`
 * @param size - number of elements of each memory pointed by a pointer of `mem`
 * @param mem - null pointer or a pointer to a vector of buffers
 */
template <typename T>
Vecp<T>::Vecp(int n, size_t size, std::vector<T*>* mem)
{
    this->n = n;
    this->size = size;
    this->mem_len = n * size;
    if (mem == nullptr) {
        this->mem = new std::vector<T*>(n, nullptr);
        for (int i = 0; i < n; i++) {
            this->mem->at(i) = new T[size];
        }
    } else {
        this->mem_alloc_case = 2;
        this->mem = mem;
    }
}

/**
 * Constructor of Vecp by cloning partly from a given vector
 *
 * @param vec - a given Vecp instance
 * @param n - number of buffers pointed by the constructed vector
 */
template <typename T>
Vecp<T>::Vecp(Vecp<T>* vec, int n)
{
    assert(n >= 0);
    int i;
    int vec_n = vec->get_n();

    this->size = vec->get_size();
    this->n = (n == 0) ? vec_n : n;
    this->mem = new std::vector<T*>(this->n, nullptr);
    this->mem_len = n * size;

    for (i = 0; i < this->n; i++) {
        this->mem->at(i) = new T[this->size];
    }

    int copy_len = (this->n <= vec_n) ? this->n : vec_n;
    for (i = 0; i < copy_len; i++) {
        std::copy_n(vec->get(i), this->size, this->mem->at(i));
    }

    if (this->n > vec_n) { // padding zeros
        for (i = vec_n; i < this->n; i++) {
            std::memset(this->mem->at(i), 0, this->size * sizeof(T));
        }
    }
}

/**
 * Constructor of Vecp by slicing from a given vector. Buffers are sliced from
 *  `begin` (inclusive) to `end` (exclusive)
 *
 * @param vec - a given Vecp instance
 * @param begin - index of buffer pointed by vec, that is the 1st buffer
 *                pointed by the output Vecp
 * @param end - index of buffer pointed by vec, that limits buffers
 *                pointed by the output Vecp
 */
template <typename T>
Vecp<T>::Vecp(Vecp<T>* vec, int begin, int end)
{
    assert(begin >= 0 && begin < end);
    assert(end <= vec->get_n());

    this->n = end - begin;
    this->size = vec->get_size();
    this->mem_len = this->n * this->size;
    this->mem_alloc_case = 1;
    this->mem = new std::vector<T*>(
        vec->get_mem()->begin() + begin, vec->get_mem()->begin() + end);
}

template <typename T>
Vecp<T>::~Vecp()
{
    if (this->mem_alloc_case < 2 && this->mem != nullptr) {
        if (this->mem_alloc_case == 0) {
            for (int i = 0; i < n; i++) {
                delete[] this->mem->at(i);
            }
        }
        this->mem->shrink_to_fit();
        delete this->mem;
    }
}

template <typename T>
inline int Vecp<T>::get_n(void)
{
    return n;
}

template <typename T>
inline size_t Vecp<T>::get_size(void)
{
    return size;
}

template <typename T>
inline int Vecp<T>::get_mem_len(void)
{
    return mem_len;
}

template <typename T>
void Vecp<T>::zero_fill(void)
{
    for (int i = 0; i < n; i++)
        std::memset(mem->at(i), 0, size * sizeof(T));
}

template <typename T>
inline void Vecp<T>::set(int i, T* buf)
{
    assert(i >= 0 && i < n);

    if ((mem_alloc_case == 0) && (this->mem->at(i) != nullptr))
        delete[] this->mem->at(i);

    this->mem->at(i) = buf;
}

template <typename T>
inline T* Vecp<T>::get(int i)
{
    assert(i >= 0 && i < n);
    return this->mem->at(i);
}

template <typename T>
inline std::vector<T*>* Vecp<T>::get_mem()
{
    return mem;
}

template <typename T>
inline void Vecp<T>::set_mem(std::vector<T*>* mem)
{
    this->mem = mem;
}

template <typename T>
void Vecp<T>::copy(Vecp<T>* v)
{
    assert(v->get_n() == n);
    assert(v->get_size() <= size);
    size_t v_size = v->get_size();
    for (int i = 0; i < n; i++)
        std::copy_n(v->get(i), v_size, this->mem->at(i));
}

template <typename T>
void Vecp<T>::separate_even_odd()
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
void Vecp<T>::separate_even_odd(Vecp<T>* even, Vecp<T>* odd)
{
    int j = 0;
    int i;
    std::vector<T*>* even_mem = even->get_mem();
    std::vector<T*>* odd_mem = odd->get_mem();
    for (i = 0; i < n; i += 2) {
        even_mem->at(j) = get(i);    // even
        odd_mem->at(j) = get(i + 1); // odd
        j++;
    }
}

template <typename T>
bool Vecp<T>::eq(Vecp<T>* v)
{
    if ((v->get_n() != n) || (v->get_size() != size))
        return false;

    for (int i = 0; i < n; i++) {
        T* a = get(i);
        T* b = v->get(i);
        for (size_t j = 0; j < size; j++) {
            if (a[j] != b[j])
                return false;
        }
    }

    return true;
}

template <typename T>
void Vecp<T>::dump(void)
{
    for (int i = 0; i < n; i++) {
        std::cout << "\n\t" << i << ": ";
        for (size_t j = 0; j < size; j++) {
            std::cout << unsigned((get(i))[j]) << "-";
        }
    }
    std::cout << "\n)\n";
}

} // namespace vec
} // namespace nttec

#endif
