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
class V2Vec : public Vec<T>
{
 private:
  Vec<T> *vec;
 public:
  explicit V2Vec(Vec<T> *vec);
  int get_n(void);
  T get(int i);
};

template <typename T>
V2Vec<T>::V2Vec(Vec<T> *vec) : Vec<T>(NULL, 0)
{
  this->vec = vec;
}

template <typename T>
int V2Vec<T>::get_n(void)
{
  return vec->n * 2;
}

template <typename T>
T V2Vec<T>::get(int i)
{
  assert(i >= 0 && i < 2*vec->n);

  if (i < vec->n)
    return vec->get(i);
  else
    return vec->get(i - vec->n);
}

