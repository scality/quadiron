/* -*- mode: c++ -*- */
#ifndef __NTTEC_VEC_DOUBLED_H__
#define __NTTEC_VEC_DOUBLED_H__

#include "vec_vector.h"

namespace nttec {
namespace vec {

/** A vector of size `n` virtually extented to `2n`.
 *
 * Example:
 *
 * In memory, the vec::Doubled looks like:
 *
 * v[0]  | v[1]   | … | v[n-1]
 * ----- | ------ | - | --------
 * v1    |  v2    | … | vN
 *
 * But it behaves likes:
 *
 * v[0]  | v[1]   | … | v[n-1] | v[n]  | v[n+1] | … | v[2n-1]
 * ----- | ------ | - | ------ | ----- | ------ | - | -----------
 * v1    |  v2    | … | vN     | v1    |  v2    | … | vN
 */
template <typename T>
class Doubled : public Vector<T> {
  public:
    explicit Doubled(Vector<T>* vec);
    int get_n(void);
    T get(int i);
    bool is_v2vec();

  private:
    Vector<T>* vec;
};

template <typename T>
bool Doubled<T>::is_v2vec()
{
    return true;
}

template <typename T>
Doubled<T>::Doubled(Vector<T>* vec)
    : Vector<T>(vec->rn, vec->get_n(), vec->get_mem(), vec->get_mem_len())
{
    this->vec = vec;
}

template <typename T>
int Doubled<T>::get_n(void)
{
    return vec->get_n() * 2;
}

template <typename T>
T Doubled<T>::get(int i)
{
    int vec_n = vec->get_n();

    assert(i >= 0 && i < 2 * vec_n);

    if (i < vec_n)
        return vec->get(i);
    else
        return vec->get(i - vec_n);
}

} // namespace vec
} // namespace nttec

#endif
