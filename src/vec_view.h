/* -*- mode: c++ -*- */
#ifndef __NTTEC_VEC_VIEW_H__
#define __NTTEC_VEC_VIEW_H__

#include "vec_vector.h"

namespace nttec {
namespace vec {

/** A view of an existing vector.
 *
 * Example:
 *
 * If you have the following Vector `vec`:
 *
 * v[0]  | v[1]   | v[2]   | v[3]   | v[4]   | v[5]  | v[6]  | v[7]  | v[8]
 * ----- | ------ | ------ | ------ | ------ | ----- | ----- | ----- | -----
 * v1    |  v2    |  v3    |  v4    |  v5    |  v6   |  v7   |  v8   |  v9
 *
 * Then, View(vec, 6, 1, 2) is:
 *
 * v[0]   | v[1]   | v[2]
 * ------ | ------ | -----
 *  v2    |  v4    |  v6
 */
template <typename T>
class View : public Vector<T> {
  public:
    explicit View(Vector<T>* vec, int n = 0, int offset = 0, int step = 1);
    int get_n(void);
    T get(int i);
    void set(int i, T val);
    void set_map(int offset, int step);
    void set_len(int n);
    void set_vec(Vector<T>* vec);

  private:
    Vector<T>* vec;
    int offset;
    int step;
    int N;
};

template <typename T>
View<T>::View(Vector<T>* vec, int n, int offset, int step)
    : Vector<T>(vec->rn, n, vec->get_mem(), vec->get_mem_len())
{
    this->vec = vec;
    this->offset = offset;
    this->step = step;
    this->n = n;
    this->N = vec->get_n();
    assert(this->n <= this->N);
}

template <typename T>
int View<T>::get_n(void)
{
    return this->n;
}

template <typename T>
T View<T>::get(int i)
{
    assert(i >= 0 && i < this->n);

    T loc = (offset + step * i) % this->N;
    return vec->get(loc);
}

template <typename T>
void View<T>::set_map(int offset, int step)
{
    this->offset = offset;
    this->step = step;
}

template <typename T>
void View<T>::set_len(int n)
{
    assert(this->n <= this->N);
    this->n = n;
}

template <typename T>
void View<T>::set_vec(Vector<T>* vec)
{
    assert(this->n <= vec->get_n());
    this->vec = vec;
    this->rn = vec->rn;
    this->N = vec->get_n();
}

template <typename T>
void View<T>::set(int i, T val)
{
    assert(i >= 0 && i < this->n);

    T loc = (offset + step * i) % this->N;
    vec->set(loc, val);
}

} // namespace vec
} // namespace nttec

#endif
