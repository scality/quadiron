/* -*- mode: c++ -*- */
#ifndef __NTTEC_FFT_2_H__
#define __NTTEC_FFT_2_H__

#include "fft_base.h"
#include "gf_base.h"
#include "vec_buffers.h"
#include "vec_vector.h"

namespace nttec {
namespace fft {

/** FFT implementation specialized for FFT size=2.
 *
 *  Since w<sup>2</sup> = 1 we have either w = -1 or w = 1 (we exclude)
 *
 *    And w<sup>-1</sup> = w = -1
 */
template <typename T>
class Size2 : public FourierTransform<T> {
  public:
    explicit Size2(gf::Field<T>* gf);
    ~Size2();
    void fft(vec::Vector<T>* output, vec::Vector<T>* input);
    void ifft(vec::Vector<T>* output, vec::Vector<T>* input);
    void fft_inv(vec::Vector<T>* output, vec::Vector<T>* input);
    void fft(vec::Buffers<T>* output, vec::Buffers<T>* input);
    void ifft(vec::Buffers<T>* output, vec::Buffers<T>* input);
    void fft_inv(vec::Buffers<T>* output, vec::Buffers<T>* input);
};

template <typename T>
Size2<T>::Size2(gf::Field<T>* gf) : FourierTransform<T>(gf, 2)
{
}

template <typename T>
Size2<T>::~Size2()
{
}

template <typename T>
void Size2<T>::fft(vec::Vector<T>* output, vec::Vector<T>* input)
{
    T a = input->get(0);
    T b = input->get(1);
    output->set(0, this->gf->add(a, b));
    output->set(1, this->gf->sub(a, b));
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void Size2<T>::fft_inv(vec::Vector<T>* output, vec::Vector<T>* input)
{
    fft(output, input);
}

template <typename T>
void Size2<T>::ifft(vec::Vector<T>* output, vec::Vector<T>* input)
{
    fft_inv(output, input);
    if (this->inv_n_mod_p > 1)
        output->mul_scalar(this->inv_n_mod_p);
}

template <typename T>
void Size2<T>::fft(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    output->copy(input);
    size_t buf_len = input->get_size();
    this->gf->add_two_bufs(input->get(1), output->get(0), buf_len);
    this->gf->sub_two_bufs(
        input->get(0), output->get(1), output->get(1), buf_len);
}

/*
 * This function performs an inverse DFT formular without a multiplication to
 * the coefficient (n^(-1) mod p)
 *
 */
template <typename T>
void Size2<T>::fft_inv(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    fft(output, input);
}

template <typename T>
void Size2<T>::ifft(vec::Buffers<T>* output, vec::Buffers<T>* input)
{
    fft_inv(output, input);
    if (this->inv_n_mod_p > 1)
        this->gf->mul_vec_to_vecp(this->vec_inv_n, output, output);
}

} // namespace fft
} // namespace nttec

#endif
