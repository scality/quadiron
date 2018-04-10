/* -*- mode: c++ -*- */
/*
 * Copyright 2017-2018 Scality
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef __NTTEC_FEC_CONTEXT_H__
#define __NTTEC_FEC_CONTEXT_H__

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <vector>
#include <sys/time.h>

#include "fft_base.h"
#include "gf_base.h"
#include "vec_poly.h"
#include "vec_zero_ext.h"

namespace nttec {

/** Forward Error Correction code implementations. */
namespace fec {

/** A class for context of decoding
 */
template <typename T>
class DecodeContext {
  public:
    DecodeContext(
        gf::Field<T>* gf,
        fft::FourierTransform<T>& fft,
        fft::FourierTransform<T>& fft_2k,
        const int k,
        const int n,
        const size_t size = 0)
    {
        this->k = k;
        this->n = n;
        this->size = size;
        this->gf = gf;
        this->fft = &fft;
        this->fft_2k = &fft_2k;
        this->vx_zero = -1;

        this->len_2k = this->gf->get_code_len_high_compo(2 * this->k);

        A = new vec::Poly<T>(this->gf, this->n);
        A_fft_2k = new vec::Poly<T>(this->gf, this->len_2k);
        inv_A_i = new vec::Poly<T>(this->gf, this->n);
        S = new vec::Poly<T>(this->gf, k);

        poly1_n = new vec::Poly<T>(this->gf, n);
        poly2_n = new vec::Poly<T>(this->gf, n);

        poly1_2k = new vec::Poly<T>(this->gf, len_2k);
        poly2_2k = new vec::Poly<T>(this->gf, len_2k);

        // zero-out all polynomials as they are used in full-length for FFT
        A->zero();
        inv_A_i->zero();
        S->zero();

        poly1_n->zero();
        poly2_n->zero();
    }

    ~DecodeContext()
    {
        if (A)
            delete A;
        if (A_fft_2k)
            delete A_fft_2k;
        if (inv_A_i)
            delete inv_A_i;
        if (S)
            delete S;
        if (poly1_n)
            delete poly1_n;
        if (poly2_n)
            delete poly2_n;
    }

    unsigned get_len_2k()
    {
        return len_2k;
    }

    vec::Poly<T>* get_A()
    {
        return A;
    }

    vec::Poly<T>* get_A_fft_2k()
    {
        return A_fft_2k;
    }

    vec::Poly<T>* get_inv_A_i()
    {
        return inv_A_i;
    }

    vec::Poly<T>* get_poly1_n()
    {
        return poly1_n;
    }

    vec::Poly<T>* get_poly2_n()
    {
        return poly2_n;
    }

    vec::Poly<T>* get_poly1_2k()
    {
        return poly1_2k;
    }

    vec::Poly<T>* get_poly2_2k()
    {
        return poly2_2k;
    }

    vec::Poly<T>* get_S()
    {
        return S;
    }

    vec::Vector<T>* get_frag_ids()
    {
        return fragments_ids;
    }

    void set_frag_ids(vec::Vector<T>* fragments_ids)
    {
        this->fragments_ids = fragments_ids;
    }

    int get_vx_zero()
    {
        return vx_zero;
    }

    void set_vx_zero(int value)
    {
        this->vx_zero = value;
    }

    void init(vec::Vector<T>* fragments_ids, const vec::Vector<T>& vx)
    {
        this->fragments_ids = fragments_ids;

        // compute A(x) = prod_j(x-x_j)
        A->set(0, 1);
        for (int i = 0; i < k; ++i) {
            A->mul_to_x_plus_coef(this->gf->sub(0, vx.get(i)));
        }

        // compute A'(x) since A_i(x_i) = A'_i(x_i)
        vec::Poly<T> _A(*A);
        _A.derivative();

        // compute A_i(x_i)
        this->fft->fft(inv_A_i, &_A);

        // compute 1/(x_i * A_i(x_i))
        // we care only about elements corresponding to fragments_ids
        for (unsigned i = 0; i < k; ++i) {
            unsigned j = fragments_ids->get(i);
            if (i != vx_zero) {
                inv_A_i->set(
                    j,
                    this->gf->inv(this->gf->mul(inv_A_i->get(j), vx.get(i))));
            } else {
                inv_A_i->set(j, this->gf->inv(inv_A_i->get(j)));
            }
        }

        // compute FFT(A) of length 2k
        if (this->fft_2k) {
            vec::ZeroExtended<T> A_2k(A, len_2k);
            this->fft_2k->fft(A_fft_2k, &A_2k);
        }
    }

    void find_vx_zero(vec::Vector<T>* betas)
    {
        vx_zero = -1;
        // vector x=(x_0, x_1, ..., x_k-1)
        for (int i = 0; i < k; ++i) {
            if (betas->get(fragments_ids->get(i)) == 0) {
                vx_zero = i;
                break;
            }
        }
    }

    void dump()
    {
        std::cout << "Dump context:\n";
        std::cout << "A:";
        A->dump();
        std::cout << "inv_A_i init:";
        inv_A_i->dump();
        std::cout << "inv_A_i:";
        inv_A_i->dump();
    }

  public:
    vec::Vector<T>* fragments_ids;

  private:
    unsigned k;
    unsigned n;
    unsigned len_2k;
    size_t size;
    int vx_zero;
    gf::Field<T>* gf;
    fft::FourierTransform<T>* fft;
    fft::FourierTransform<T>* fft_2k;

    vec::Poly<T>* A = nullptr;
    vec::Poly<T>* A_fft_2k = nullptr;
    vec::Poly<T>* inv_A_i = nullptr;
    vec::Poly<T>* S = nullptr;

    vec::Poly<T>* poly1_n = nullptr;
    vec::Poly<T>* poly2_n = nullptr;
    vec::Poly<T>* poly1_2k = nullptr;
    vec::Poly<T>* poly2_2k = nullptr;
};

} // namespace fec
} // namespace nttec

#endif
