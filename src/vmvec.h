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
  explicit VmVec(Vec<T> *vec, int n, int offset, int step);
  int get_n(void);
  T get(int i);
  void set(int i, T val);
  Vec<T>* toVec();
  void import(Vec<T>* v);
};

template <typename T>
VmVec<T>::VmVec(Vec<T> *vec, int n, int offset, int step) : Vec<T>(NULL, 0)
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
void VmVec<T>::set(int i, T val)
{
  assert(i >= 0 && i < this->n);

  T loc = (offset + step * i) % this->N;
  vec->set(loc, val);
}

template <typename T>
Vec<T>* VmVec<T>::toVec()
{
  Vec<T>* v = new Vec<T>(this->gf, this->n);
  for (int i = 0; i < this->n; i++)
    v->set(i, get(i));
  return v;
}

template <typename T>
void VmVec<T>::import(Vec<T> *v)
{
  assert(v->n == this->n);

  for (int i = 0; i < this->n; i++)
    this->set(i, v->get(i));
}
