/* -*- mode: c++ -*- */
#pragma once

template<typename T>
class Vec
{
 private:
  T *mem;
 public:
  GF<T> *gf;
  int n;
  Vec(GF<T> *gf, int n);
  ~Vec();
  void zero_fill(void);
  void set(int i, T val);
  virtual T get(int i);
  void mul_scalar(T scalar);
  void dump(void);
};

template <typename T>
Vec<T>::Vec(GF<T> *gf, int n)
{
  this->gf = gf;
  this->n = n;
  this->mem = new T[n];
}

template <typename T>
Vec<T>::~Vec()
{
  delete[] this->mem;
}

template <typename T>
void Vec<T>::zero_fill(void)
{
  int i;

  for (i = 0;i < n;i++)
    set(i, 0);
}

template <typename T>
void Vec<T>::set(int i, T val)
{
  assert(i >= 0 && i < n);
  
  mem[i] = val;
}

template <typename T>
T Vec<T>::get(int i)
{
  assert(i >= 0 && i < n);
  
  return mem[i];
}

/** 
 * Multiplication of a vector by a scalar
 * 
 * @param scalar 
 */
template <typename T>
void Vec<T>::mul_scalar(T scalar)
{
  for (int i = 0; i < n;i++)
    set(i, gf->mul(get(i), scalar));
}

template <typename T>
void Vec<T>::dump(void)
{
  std::cout << "( ";
  for (int i = 0;i < n;i++)
    std::cout << get(i) << " ";
  std::cout << ")\n";
}
