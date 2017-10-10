/* -*- mode: c++ -*- */
#pragma once

template<typename T>
class Poly;

template<typename T>
class Vec
{
 private:
  T *mem;
  int mem_len;
  bool new_mem;
 protected:
  int n;
 public:
  MG<T> *mg;
  Vec(MG<T> *mg, int n, T* mem=NULL, int mem_len=0);
  ~Vec();
  virtual int get_n(void);
  int get_mem_len(void);
  void zero_fill(void);
  virtual void set(int i, T val);
  virtual T get(int i);
  virtual T *get_mem();
  virtual bool is_v2vec();
  void mul_scalar(T scalar);
  void hadamard_mul(Vec<T> *v);
  void add(Vec<T> *v);
  void add(Vec<T> *v, int offset);
  void add_mutual(Vec<T> *v);
  void add_mutual(Vec<T> *v, int offset);
  void add_mutual(Vec<T> *v, int offset, int len);
  void copy(Vec<T> *v);
  void copy(Vec<T> *v, int n);
  bool eq(Vec<T> *v);
  void to_poly(Poly<T> *poly);
  virtual void dump(void);
};

template <typename T>
Vec<T>::Vec(MG<T> *mg, int n, T* mem, int mem_len)
{
  this->mg = mg;
  this->n = n;
  if (mem == NULL) {
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
  if (new_mem) delete[] this->mem;
}

template <typename T>
int Vec<T>::get_n(void)
{
  return n;
}

template <typename T>
int Vec<T>::get_mem_len(void)
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
    set(i, mg->mul(get(i), scalar));
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
    set(i, mg->mul(get(i), v->get(i)));
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
    set(i, mg->add(get(i), v->get(i)));
}

template <typename T>
void Vec<T>::add(Vec<T> *v, int offset)
{
  assert(n >= v->get_n() + offset);

  int j;
  for (int i = 0; i < v->get_n(); i++) {
    j = i + offset;
    set(j, mg->add(get(j), v->get(i)));
  }
}

template <typename T>
void Vec<T>::add_mutual(Vec<T> *v)
{
  assert(n >= v->get_n());
  for (int i = 0; i < v->get_n(); i++)
    set(i, mg->add(get(i), v->get(i)));
}

template <typename T>
void Vec<T>::add_mutual(Vec<T> *v, int offset)
{
  assert(v->get_n() == 0 || n - offset >= v->get_n());
  int j;
  for (int i = 0; i < v->get_n(); i++) {
    j = i + offset;
    set(j, mg->add(get(j), v->get(i)));
  }
}

template <typename T>
void Vec<T>::add_mutual(Vec<T> *v, int offset, int len)
{
  assert(len == 0 || n - offset >= len);
  assert(v->get_n() >= len);
  int j;
  for (int i = 0; i < len; i++) {
    j = i + offset;
    set(j, mg->add(get(j), v->get(i)));
  }
}

template <>
void Vec<uint32_t>::add(Vec<uint32_t> *v);
template <>
void Vec<uint64_t>::add(Vec<uint64_t> *v);

template <typename T>
void Vec<T>::copy(Vec<T> *v)
{
  assert(v->get_mem_len() <= this->mem_len);
  std::copy_n(v->get_mem(), v->get_mem_len(), this->mem);
}

template <typename T>
void Vec<T>::copy(Vec<T> *v, int n)
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
bool Vec<T>::eq(Vec<T> *v)
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
void Vec<T>::to_poly(Poly<T> *poly)
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
