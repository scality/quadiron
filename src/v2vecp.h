/* -*- mode: c++ -*- */
#pragma once

/**
 * Virtual (v->get_n()*2) x 1 vertical vector for the need of
 * Cooley-Tukey algorithm
 *
 * | v_0 |
 * | v_1 |
 * | ...
 * | v_b |
 * | v_0 |
 * | v_1 |
 * | ...
 * | v_b |
 */
template<typename T>
class V2Vecp : public Vecp<T>
{
 private:
  Vecp<T> *vec;
 public:
  explicit V2Vecp(Vecp<T> *vec);
  T* get(int i);
};

template <typename T>
V2Vecp<T>::V2Vecp(Vecp<T> *vec) :
  Vecp<T>(2*vec->get_n(), vec->get_size(), vec->get_mem())
{
  this->vec = vec;
}

template <typename T>
T* V2Vecp<T>::get(int i)
{
  int vec_n = vec->get_n();
  assert(i >= 0 && i < 2*vec_n);

  if (i < vec_n)
    return vec->get(i);
  else
    return vec->get(i - vec_n);
}
