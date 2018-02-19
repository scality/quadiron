/* -*- mode: c++ -*- */
#ifndef __NTTEC_VEC_ZERO_EXT_H__
#define __NTTEC_VEC_ZERO_EXT_H__

#include "vec_vector.h"

namespace nttec {
namespace vec {

/** A vector of size `n` virtually extented with zeros.
 *
 * In memory, the vec::ZeroExtended looks like:
 *
 * v[0]  | v[1]   | … | v[n-1]
 * ----- | ------ | - | --------
 * v1    |  v2    | … | vN
 *
 * But it behaves likes:
 *
 * v[0]  | v[1]   | … | v[n-1] | v[n]  | v[n+1] | …
 * ----- | ------ | - | ------ | ----- | ------ | -----
 * v1    |  v2    | … | vN     | 0     |  0     | 0
 */
template <typename T>
class ZeroExtended : public Vector<T> {
  public:
    ZeroExtended(Vector<T>* vec, int n);
    int get_n(void);
    T get(int i);

  private:
    Vector<T>* vec;
    int n;
};

template <typename T>
ZeroExtended<T>::ZeroExtended(Vector<T>* vec, int n)
    : Vector<T>(vec->rn, n, vec->get_mem(), vec->get_mem_len())
{
    this->vec = vec;
    this->n = n;
}

template <typename T>
int ZeroExtended<T>::get_n(void)
{
    return n;
}

template <typename T>
T ZeroExtended<T>::get(int i)
{
    assert(i >= 0 && i < n);

    return (i < vec->get_n()) ? vec->get(i) : 0;
}

} // namespace vec
} // namespace nttec

#endif
