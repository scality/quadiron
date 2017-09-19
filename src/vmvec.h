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
VmVec<T>::VmVec(Vec<T> *vec, int n, int offset, int step) : Vec<T>(vec->gf, n)
{
  this->vec = vec;
  this->offset = offset;
  this->step = step;
  this->n = n;
  this->N = vec->n;
  assert(this->n <= this->N);
}

template <typename T>
int VmVec<T>::get_n(void)
{
  return this->n;
}

template <typename T>
T VmVec<T>::get(int i)
{
  assert(i >= 0 && i < this->n);

  T loc = (offset + step * i) % this->N;
  return vec->get(loc);
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
  assert(this->n <= vec->n);
  this->vec = vec;
  this->gf = vec->gf;
  this->N = vec->n;
}

template <typename T>
void VmVec<T>::set(int i, T val)
{
  assert(i >= 0 && i < this->n);

  T loc = (offset + step * i) % this->N;
  vec->set(loc, val);
}
