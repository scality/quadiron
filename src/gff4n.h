/* -*- mode: c++ -*- */
#pragma once
#include "ntl.h"

const uint32_t F4 = 65537;

template<typename T>
class Poly;

template<typename T>
class GFF4N : public GF<T>
{
private:
  T unit;
  T q;
  T h;
  GF<T> *sub_field;
  bool check_n(int n);
  void init(void);

public:
  explicit GFF4N(int n);
  ~GFF4N();
  T card(void);
  T card_minus_one(void);
  T neg(T a);
  T add(T a, T b);
  T sub(T a, T b);
  T mul(T a, T b);
  T div(T a, T b);
  T inv(T a);
  T exp(T a, T b);
  T log(T a, T b);
  T weak_rand(void);
  T get_unit(void);
};

template <typename T>
GFF4N<T>::GFF4N(int n) : GF<T>(65537, n)
{
  sub_field = new GFP<T>(65537);

  if (!check_n(n)) {
    // not supported yet
    assert("Input n is not supported for GFF4N");
  }

  unit = 0;
  q = 0;
  h = 0;
  init();
}

template <typename T>
GFF4N<T>::~GFF4N()
{
  delete sub_field;
}

/*
 * Each element of sub-field GFP(65537) is stored in 32bits
 * Hence, n <= sizeof(T) / 4
 */
template <typename T>
bool GFF4N<T>::check_n(int n)
{
  return (n <= sizeof(T) / 4);
}

template <typename T>
void GFF4N<T>::init(void)
{
  unit = 1;
  q = 65537;
  h = 65536;
  for (int i = 1; i < this->n; i++) {
    unit += (unit << 32);
    q += (q << 32);
    h += (h << 32);
  }
}

// FIXME: n >= sizeof(T) / sizeof(F4)
template <typename T>
T GFF4N<T>::card(void)
{
  return q;
}

template <typename T>
T GFF4N<T>::card_minus_one(void)
{
  return h;
}

template <typename T>
T GFF4N<T>::get_unit(void)
{
  return unit;
}

template <typename T>
T GFF4N<T>::neg(T a)
{
  return sub(0, a);
}

template <typename T>
T GFF4N<T>::add(T a, T b)
{
  T c = 0;
  uint32_t mask = 0xffffffff;
  int shift = 0;
  for (int i = 0; i < this->n; i++) {
    T tmp = sub_field->add(a & mask, b & mask);
    tmp <<= shift;
    c += tmp;
    shift += 32;
    a >>= 32;
    b >>= 32;
  }
  return c;
}

template <typename T>
T GFF4N<T>::sub(T a, T b)
{
  T c = 0;
  uint32_t mask = 0xffffffff;
  int shift = 0;
  for (int i = 0; i < this->n; i++) {
    T tmp = sub_field->sub(a & mask, b & mask);
    tmp <<= shift;
    c += tmp;
    shift += 32;
    a >>= 32;
    b >>= 32;
  }
  return c;
}

template <typename T>
T GFF4N<T>::mul(T a, T b)
{
  T c = 0;
  uint32_t mask = 0xffffffff;
  int shift = 0;
  for (int i = 0; i < this->n; i++) {
    T tmp = sub_field->mul(a & mask, b & mask);
    tmp <<= shift;
    c += tmp;
    shift += 32;
    a >>= 32;
    b >>= 32;
  }
  return c;
}

template <typename T>
T GFF4N<T>::div(T a, T b)
{
  T c = 0;
  uint32_t mask = 0xffffffff;
  int shift = 0;
  for (int i = 0; i < this->n; i++) {
    T tmp = sub_field->div(a & mask, b & mask);
    tmp <<= shift;
    c += tmp;
    shift += 32;
    a >>= 32;
    b >>= 32;
  }
  return c;
}

template <typename T>
T GFF4N<T>::inv(T a)
{
  T c = 0;
  uint32_t mask = 0xffffffff;
  int shift = 0;
  for (int i = 0; i < this->n; i++) {
    T tmp = sub_field->inv(a & mask);
    tmp <<= shift;
    c += tmp;
    shift += 32;
    a >>= 32;
  }
  return c;
}

template <typename T>
T GFF4N<T>::exp(T a, T b)
{
  T c = 0;
  uint32_t mask = 0xffffffff;
  int shift = 0;
  for (int i = 0; i < this->n; i++) {
    T tmp = sub_field->exp(a & mask, b & mask);
    tmp <<= shift;
    c += tmp;
    shift += 32;
    a >>= 32;
    b >>= 32;
  }
  return c;
}

template <typename T>
T GFF4N<T>::log(T a, T b)
{
  T c = 0;
  uint32_t mask = 0xffffffff;
  int shift = 0;
  for (int i = 0; i < this->n; i++) {
    T tmp = sub_field->log(a & mask, b & mask);
    tmp <<= shift;
    c += tmp;
    shift += 32;
    a >>= 32;
    b >>= 32;
  }
  return c;
}

template <typename T>
T GFF4N<T>::weak_rand(void)
{
  T c = 0;
  uint32_t mask = 0xffffffff;
  int shift = 0;
  for (int i = 0; i < this->n; i++) {
    T tmp = sub_field->weak_rand();
    tmp <<= shift;
    c += tmp;
    shift += 32;
  }
  return c;
}
