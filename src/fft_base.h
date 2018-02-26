/* -*- mode: c++ -*- */
#ifndef __NTTEC_FFT_BASE_H__
#define __NTTEC_FFT_BASE_H__

#include "gf_base.h"
#include "vec_buffers.h"
#include "vec_vector.h"

namespace nttec {

/** Various Fast Fourier Transform (FFT) implementations. */
namespace fft {

/** Base class for Fourier Transform on Galois Fields.
 *
 *  The Fourier Transform is applied on vector of size `n` with ω as
 *  n<sup>th</sup> root of unity.
 */
template <typename T>
class FourierTransform {
  public:
    virtual ~FourierTransform();
    int get_n();
    gf::Field<T>* get_gf();
    /** Compute the Fourier Transform. */
    virtual void fft(vec::Vector<T>* output, vec::Vector<T>* input) = 0;
    virtual void fft(vec::Buffers<T>* output, vec::Buffers<T>* input){};
    /** Compute the Inverse Fourier Transform. */
    virtual void ifft(vec::Vector<T>* output, vec::Vector<T>* input) = 0;
    virtual void ifft(vec::Buffers<T>* output, vec::Buffers<T>* input){};
    /** Compute the summation for the inverse FFT formula. */
    virtual void fft_inv(vec::Vector<T>* output, vec::Vector<T>* input) = 0;
    virtual void fft_inv(vec::Buffers<T>* output, vec::Buffers<T>* input){};

  protected:
    gf::Field<T>* gf;
    int n;
    T inv_n_mod_p;
    vec::Vector<T>* vec_inv_n = nullptr;
    FourierTransform(gf::Field<T>* gf, int n);
};

template <typename T>
FourierTransform<T>::FourierTransform(gf::Field<T>* gf, int n)
{
    this->gf = gf;
    this->n = n;
    this->inv_n_mod_p = gf->inv(n) % gf->get_sub_field()->card();

    this->vec_inv_n = new vec::Vector<T>(this->gf, n);
    for (int i = 0; i < n; i++) {
        this->vec_inv_n->set(i, this->inv_n_mod_p);
    }
}

template <typename T>
FourierTransform<T>::~FourierTransform()
{
    if (vec_inv_n != nullptr)
        delete vec_inv_n;
}

template <typename T>
int FourierTransform<T>::get_n()
{
    return n;
}

template <typename T>
gf::Field<T>* FourierTransform<T>::get_gf()
{
    return gf;
}

} // namespace fft
} // namespace nttec

#endif
