/* -*- mode: c++ -*- */
#ifndef __NTTEC_VEC_BUF_DOUBLED_H__
#define __NTTEC_VEC_BUF_DOUBLED_H__

#include "vec_buffers.h"

namespace nttec {
namespace vec {

/** A vector of `n` buffers virtually extented to `2n`.
 *
 * Physically, BuffersDoubled contains `n` buffers but behaves as if it
 * contained `2n` buffers.
 *
 * Example:
 *
 * We have a set of `N` independent buffers:
 *
 * buf1    | buf2    | … | bufN
 * ------- | ------- | - | -------
 * buf1[0] | buf2[0] | … | bufN[0]
 * buf1[1] | buf2[1] | … | bufN[1]
 * …       | …       | … | …
 *
 * In memory, the BuffersDoubled looks like:
 *
 * v[0]  | v[1]   | … | v[n-1]
 * ----- | ------ | - | --------
 * &buf1 |  &buf2 | … | &bufN
 *
 * But it behaves likes:
 *
 * v[0]  | v[1]   | … | v[n-1] | v[n]  | v[n+1] | … | v[2n-1]
 * ----- | ------ | - | ------ | ----- | ------ | - | -----------
 * &buf1 |  &buf2 | … | &bufN  | &buf1 |  &buf2 | … | &bufN
 */
template <typename T>
class BuffersDoubled : public Buffers<T> {
  public:
    explicit BuffersDoubled(Buffers<T>* vec);
    T* get(int i);

  private:
    Buffers<T>* vec;
};

template <typename T>
BuffersDoubled<T>::BuffersDoubled(Buffers<T>* vec)
    : Buffers<T>(2 * vec->get_n(), vec->get_size(), vec->get_mem())
{
    this->vec = vec;
}

template <typename T>
T* BuffersDoubled<T>::get(int i)
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
