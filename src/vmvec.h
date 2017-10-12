/* -*- mode: c++ -*- */
#pragma once

/**
 * Virtual vector on a base vector from an offset and a step
 *
 */
template<typename T>
class VmVec : public Vec<T>
{
 private:
  Vec<T> *vec;
  int offset;
  int step;
  int N;
  int loc(int i);
 public:
  explicit VmVec(Vec<T> *vec, int n = 0, int offset = 0, int step = 1);
  int get_n(void);
  T get(int i);
  void set(int i, T val);
  void set_map(int offset, int step);
  void set_len(int n);
  void set_vec(Vec<T> *vec);
};

template <typename T>
VmVec<T>::VmVec(Vec<T> *vec, int n, int offset, int step) :
  Vec<T>(vec->rn, n, vec->get_mem() + offset, vec->get_mem_len())
{
  this->vec = vec;
  this->offset = offset;
  this->step = step;
  this->N = vec->get_n();
  if (offset >= this->N)
    this->n = 0;
  else
    this->n = n;
  assert(this->n <= this->N);
}

template <typename T>
inline int VmVec<T>::loc(int i)
{
  assert(i >= 0 && i < this->n);
  return offset + step * i;
}

template <typename T>
int VmVec<T>::get_n(void)
{
  return this->n;
}

template <typename T>
inline T VmVec<T>::get(int i)
{
  int j = loc(i);
  if (j >= N)
    return 0;
  return vec->get(j);
}

template <typename T>
void VmVec<T>::set_map(int offset, int step)
{
  this->offset = offset;
  this->step = step;
}

template <typename T>
void VmVec<T>::set_len(int n)
{
  assert(this->n <= this->N);
  this->n = n;
}

template <typename T>
void VmVec<T>::set_vec(Vec<T> *vec)
{
  assert(this->n <= vec->get_n());
  this->vec = vec;
  this->rn = vec->rn;
  this->N = vec->get_n();
}

template <typename T>
inline void VmVec<T>::set(int i, T val)
{
  int j = loc(i);
  if (j < N)
    vec->set(j, val);
}
