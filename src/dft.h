/* -*- mode: c++ -*- */
#ifndef __NTTEC_DFT_H__
#define __NTTEC_DFT_H__

#include "gf.h"
#include "vec.h"

namespace nttec {

/** Various Fast Fourier Transform implementations. */
namespace fft {

/*
 * DFT over the field gf and over vectors of size n with w as n-th
 * root of unity
 */
template <typename T>
class DFT {
  protected:
    gf::GF<T>* gf;
    int n;
    T inv_n_mod_p;
    vec::Vec<T>* vec_inv_n = nullptr;
    DFT(gf::GF<T>* gf, int n);

  public:
    virtual ~DFT();
    int get_n();
    gf::GF<T>* get_gf();
    virtual void fft(vec::Vec<T>* output, vec::Vec<T>* input) = 0;
    virtual void ifft(vec::Vec<T>* output, vec::Vec<T>* input) = 0;
    virtual void fft_inv(vec::Vec<T>* output, vec::Vec<T>* input) = 0;
};

template <typename T>
DFT<T>::DFT(gf::GF<T>* gf, int n)
{
    this->gf = gf;
    this->n = n;
    this->inv_n_mod_p = gf->inv(n) % gf->get_sub_field()->card();

    this->vec_inv_n = new vec::Vec<T>(this->gf, n);
    for (int i = 0; i < n; i++) {
        this->vec_inv_n->set(i, this->inv_n_mod_p);
    }
}

template <typename T>
DFT<T>::~DFT()
{
    if (vec_inv_n != nullptr)
        delete vec_inv_n;
}

template <typename T>
int DFT<T>::get_n()
{
    return n;
}

template <typename T>
gf::GF<T>* DFT<T>::get_gf()
{
    return gf;
}

} // namespace fft
} // namespace nttec

#endif
