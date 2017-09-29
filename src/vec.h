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
  virtual int get_n(void);
  void zero_fill(void);
  virtual void set(int i, T val);
  virtual T get(int i);
  virtual T *get_mem();
  virtual bool is_v2vec();
  void mul_scalar(T scalar);
  void hadamard_mul(Vec<T> *v);
  void add(Vec<T> *v);
  void copy(Vec<T> *v);
  void copy(Vec<T> *v, int n);
  bool eq(Vec<T> *v);
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
int Vec<T>::get_n(void)
{
  return n;
}

template <typename T>
void Vec<T>::zero_fill(void)
{
  int i;

  for (i = 0; i < n; i++)
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

template <typename T>
T *Vec<T>::get_mem()
{
  return mem;
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
    set(i, gf->mul(get(i), scalar));
}

/**
 * entrywise product
 *
 * @param v
 */
template <typename T>
void Vec<T>::hadamard_mul(Vec<T> *v)
{
  assert(n == v->get_n());

  for (int i = 0; i < n; i++)
    set(i, gf->mul(get(i), v->get(i)));
}

template <>
void Vec<uint32_t>::hadamard_mul(Vec<uint32_t> *v);
template <>
void Vec<uint64_t>::hadamard_mul(Vec<uint64_t> *v);

template <typename T>
void Vec<T>::add(Vec<T> *v)
{
  assert(n == v->get_n());

  for (int i = 0; i < n; i++)
    set(i, gf->add(get(i), v->get(i)));
}

template <>
void Vec<uint32_t>::add(Vec<uint32_t> *v);
template <>
void Vec<uint64_t>::add(Vec<uint64_t> *v);

template <typename T>
void Vec<T>::copy(Vec<T> *v)
{
  assert(v->n == this->n);

  std::copy_n(v->mem, v->n, this->mem);
}

template <typename T>
void Vec<T>::copy(Vec<T> *v, int n)
{
  assert(n <= this->n);

  std::copy_n(v->mem, n, this->mem);
}

template <typename T>
bool Vec<T>::eq(Vec<T> *v)
{
  if (v->n != this->n)
    return false;

  for (int i = 0; i < v->n; i++) {
    if (get(i) != v->get(i))
      return false;
  }

  return true;
}

template <typename T>
void Vec<T>::dump(void)
{
  std::cout << "( ";
  for (int i = 0; i < n; i++)
    std::cout << get(i) << " ";
  std::cout << ")\n";
}
