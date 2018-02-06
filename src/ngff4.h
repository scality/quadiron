/* -*- mode: c++ -*- */
#ifndef __NTL_NGFF4_H__
#define __NTL_NGFF4_H__

#include "gf.h"
#include "gfp.h"

#define MASK16  0xFFFF
#define MASK32  0xFFFFFFFF

template<typename T>
class Poly;

template<typename T>
class NGFF4 : public GF<T>
{
private:
  T unit;
  T q;
  T h;
  GF<uint32_t> *sub_field;
  bool check_n(unsigned n);
  void init(void);

  T expand16(uint16_t *arr);
  T expand32(uint32_t *arr);
  // to debug
  void show_arr(uint32_t *arr);

public:
  explicit NGFF4(unsigned n);
  ~NGFF4();
  T card(void);
  T card_minus_one(void);
  bool check(T a);
  T neg(T a);
  T add(T a, T b);
  T sub(T a, T b);
  T mul(T a, T b);
  T div(T a, T b);
  T inv(T a);
  T exp(T a, T b);
  T log(T a, T b);
  T weak_rand_tuple(void);
  T weak_rand(void);
  T get_unit(void);
  T replicate(T a);
  T pack(T a);
  T pack(T a, uint32_t flag);
  compT<T> unpack(T a);
  T get_nth_root(T n);
  void compute_omegas(Vec<T> *W, int n, T w);
  GF<uint32_t> *get_sub_field();
};

template <typename T>
NGFF4<T>::NGFF4(unsigned n) : GF<T>(65537, n)
{
  sub_field = new GFP<uint32_t>(65537);

  if (!check_n(n)) {
    // not supported yet
    assert("Input n is not supported for NGFF4");
  }

  unit = 0;
  q = 0;
  h = 0;
  init();
}

template <typename T>
NGFF4<T>::~NGFF4()
{
  delete sub_field;
}

/*
 * Each element of sub-field GFP(65537) is stored in 32bits
 * Hence, n <= sizeof(T) / 4
 */
template <typename T>
bool NGFF4<T>::check_n(unsigned n)
{
  return (n <= sizeof(T) / 4);
}

template <typename T>
inline T NGFF4<T>::expand16(uint16_t *arr)
{
  T c = arr[this->n - 1];
  for (int i = this->n - 2; i >= 0; i--) {
    c = (c << 16) | arr[i];
  }
  return c;
}

template <typename T>
inline T NGFF4<T>::expand32(uint32_t *arr)
{
  T c = arr[this->n - 1];
  for (int i = this->n - 2; i >= 0; i--) {
    c = ((c << 16) << 16) | arr[i];
  }
  return c;
}

template <typename T>
void NGFF4<T>::show_arr(uint32_t *arr)
{
  std::cout << "\t arr: ";
  for (int i = 0; i < this->n; i++) {
    std::cout << arr[i] << " ";
  }
  std::cout << std::endl;
}

template <typename T>
T NGFF4<T>::replicate(T a)
{
  T b = a;
  for (int i = 1; i < this->n; i++) {
    b = ((b << 16) << 16) | a;
  }
  return b;
}

template <typename T>
void NGFF4<T>::init(void)
{
  unit = replicate(1);
  q = replicate(65537);
  h = replicate(65536);
}

template <typename T>
T NGFF4<T>::card(void)
{
  return q;
}

template <typename T>
T NGFF4<T>::card_minus_one(void)
{
  return h;
}

template <typename T>
T NGFF4<T>::get_unit(void)
{
  return unit;
}

template <typename T>
T NGFF4<T>::neg(T a)
{
  return sub(0, a);
}

template <typename T>
T NGFF4<T>::add(T a, T b)
{
  uint32_t arr[this->n];

  arr[0] = ((uint32_t)(a & MASK32) + (uint32_t)(b & MASK32)) % 65537;
  for (int i = 1; i < this->n; i++) {
    a = (a >> 16) >> 16;
    b = (b >> 16) >> 16;
    arr[i] = ((uint32_t)(a & MASK32) + (uint32_t)(b & MASK32)) % 65537;
  }

  T c = expand32(arr);

  return c;
}

template <typename T>
T NGFF4<T>::sub(T a, T b)
{
  uint32_t arr[this->n];
  uint32_t ae, be;

  ae = (a & MASK32);
  be = (b & MASK32);
  if (ae >= be)
    arr[0] = ae - be;
  else
    arr[0] = 65537 + ae - be;
  for (int i = 1; i < this->n; i++) {
    a = (a >> 16) >> 16;
    b = (b >> 16) >> 16;
    ae = (uint32_t)(a & MASK32);
    be = (uint32_t)(b & MASK32);
    if (ae >= be)
      arr[i] = ae - be;
    else
      arr[i] = 65537 + ae - be;
  }

  T c = expand32(arr);

  return c;
}

template <typename T>
T NGFF4<T>::mul(T a, T b)
{
  uint32_t arr[this->n];
  uint64_t ae;
  uint32_t be;

  ae = (uint64_t)(a & MASK32);
  be = (uint32_t)(b & MASK32);
  if (ae == 65536 && be == 65536)
    arr[0] = 1;
  else
    arr[0] = (ae * be) % 65537;
  for (int i = 1; i < this->n; i++) {
    a = (a >> 16) >> 16;
    b = (b >> 16) >> 16;
    ae = (uint64_t)(a & MASK32);
    be = (uint32_t)(b & MASK32);
    if (ae == 65536 && be == 65536)
      arr[i] = 1;
    else
      arr[i] = (ae * be) % 65537;
  }

  T c = expand32(arr);
  return c;
}

template <typename T>
T NGFF4<T>::div(T a, T b)
{
  uint32_t arr[this->n];

  arr[0] = sub_field->div((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
  for (int i = 1; i < this->n; i++) {
    a = (a >> 16) >> 16;
    b = (b >> 16) >> 16;
    arr[i] = sub_field->div((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
  }

  T c = expand32(arr);
  return c;
}

template <typename T>
T NGFF4<T>::inv(T a)
{
  uint32_t arr[this->n];

  arr[0] = sub_field->inv((uint32_t)(a & MASK32));
  for (int i = 1; i < this->n; i++) {
    a = (a >> 16) >> 16;
    arr[i] = sub_field->inv((uint32_t)(a & MASK32));
  }

  T c = expand32(arr);
  return c;
}

template <typename T>
T NGFF4<T>::exp(T a, T b)
{
  uint32_t arr[this->n];

  arr[0] = sub_field->exp((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
  for (int i = 1; i < this->n; i++) {
    a = (a >> 16) >> 16;
    b = (b >> 16) >> 16;
    arr[i] = sub_field->exp((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
  }

  T c = expand32(arr);
  return c;
}

template <typename T>
T NGFF4<T>::log(T a, T b)
{
  uint32_t arr[this->n];

  arr[0] = sub_field->log((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
  for (int i = 1; i < this->n; i++) {
    a = (a >> 16) >> 16;
    b = (b >> 16) >> 16;
    arr[i] = sub_field->log((uint32_t)(a & MASK32), (uint32_t)(b & MASK32));
  }

  T c = expand32(arr);
  return c;
}

template <typename T>
T NGFF4<T>::weak_rand_tuple()
{
  T c = sub_field->weak_rand();
  for (int i = 1; i < this->n; i++) {
    c = ((c << 16) << 16) | sub_field->weak_rand();
  }
  return c;
}

template <typename T>
T NGFF4<T>::weak_rand(void)
{
  T c = weak_rand_tuple();
  return unpack(c).val;
}

/**
 * Pack of n numbers each of 16 bits into n numbers each of 32 bits
 */
template <typename T>
T NGFF4<T>::pack(T a)
{
  uint32_t arr[this->n];
  arr[0] = (uint32_t)(a & MASK16);
  for (int i = 1; i < this->n; i++) {
    a = (a >> 16);
    arr[i] = (uint32_t)(a & MASK16);
  }

  T c = expand32(arr);
  return c;
}

/**
 * Pack of n numbers each of 16 bits into n numbers each of 32 bits
 *  If flag contains 2^i, ith number == 65537
 */
template <typename T>
T NGFF4<T>::pack(T a, uint32_t flag)
{
  uint32_t arr[this->n];
  if (flag & 1)
    arr[0] = 65536;
  else
    arr[0] = (uint32_t)(a & MASK16);
  for (int i = 1; i < this->n; i++) {
    flag >>= 1;
    a = (a >> 16);
    if (flag & 1)
      arr[i] = 65536;
    else
      arr[i] = (uint32_t)(a & MASK16);
  }

  T c = expand32(arr);
  return c;
}

/**
 * Unpack of n numbers each of 32 bits into n numbers each of 16 bits
 *  If ith number == 65537, add 2^i to flag to mark it, and consider this number
 *  as zero
 */
template <typename T>
compT<T> NGFF4<T>::unpack(T a)
{
  compT<T> b = compT<T>();
  uint32_t flag = 0;
  uint32_t ae;
  uint16_t arr[this->n];

  ae = (uint32_t)(a & MASK32);
  if (ae == 65536) {
    flag |= 1;
    arr[0] = 0;
  } else
    arr[0] = (uint16_t)ae;
  for (int i = 1; i < this->n; i++) {
    a = (a >> 16) >> 16;
    ae = (uint32_t)(a & MASK32);
    if (ae == 65536) {
      flag |= (1 << i);
      arr[i] = 0;
    } else
      arr[i] = ae;
  }

  b.flag = flag;
  b.val = expand16(arr);
  return b;
}

// Use for fft
template <typename T>
T NGFF4<T>::get_nth_root(T n)
{
  T sub_w = sub_field->get_nth_root(n);
  T w = replicate(sub_w);
  return w;
}

template <typename T>
bool NGFF4<T>::check(T a)
{
  return (a >= 0 && a < q);
}

/**
 * Compute the different powers of the root of unity into a vector
 *
 * @param W output vector (must be of length n+1)
 * @param n
 * @param w n-th root of unity
 */
template <typename T>
void NGFF4<T>::compute_omegas(Vec<T> *W, int n, T w)
{
  for (int i = 0; i < n; i++) {
    W->set(i, this->exp(w, replicate(i)));
  }
}

template <typename T>
GF<uint32_t> *NGFF4<T>::get_sub_field()
{
  return sub_field;
}

#endif
