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
#include "gf_nf4.h"
#include "vec_poly.h"
#include "vec_zero_ext.h"

namespace nttec {

namespace fec {

/* List of vectors handled by context */
enum class CtxVec { A_FFT_2K = 0, INV_A_I, N1, N2, V2K1, V2K2 };

/* List of polynomials handled by context */
enum class CtxPoly { A = 0, S };

/** A class for context of decoding
 */
template <typename T>
class DecodeContext {
  public:
    DecodeContext(
        const gf::Field<T>& gf,
        fft::FourierTransform<T>& fft,
        fft::FourierTransform<T>& fft_2k,
        const vec::Vector<T>& fragments_ids,
        const vec::Vector<T>& vx,
        const int k,
        const int n,
        int vx_zero = -1,
        const size_t size = 0)
    {
        this->k = k;
        this->n = n;
        this->size = size;
        this->gf = &gf;
        this->fft = &fft;
        this->fft_2k = &fft_2k;
        this->vx_zero = vx_zero;

        this->len_2k = this->gf->get_code_len_high_compo(2 * this->k);

        this->fragments_ids = &fragments_ids;

        A = std::unique_ptr<vec::Poly<T>>(new vec::Poly<T>(gf, n));
        A_fft_2k =
            std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, len_2k));
        inv_A_i = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, n));
        S = std::unique_ptr<vec::Poly<T>>(new vec::Poly<T>(gf, k));

        vec1_n = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, n));
        vec2_n = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, n));

        vec1_2k =
            std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, len_2k));
        vec2_2k =
            std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, len_2k));

        // zero-out all polynomials as they are used in full-length for FFT
        A->zero_fill();
        inv_A_i->zero_fill();
        S->zero_fill();

        vec1_n->zero_fill();
        vec2_n->zero_fill();

        init(vx);
    }

    ~DecodeContext() {}

    unsigned get_len_2k() const
    {
        return len_2k;
    }

    const vec::Vector<T>& get_fragments_id() const
    {
        return *fragments_ids;
    }

    vec::Vector<T>& get_vector(CtxVec type) const
    {
        switch (type) {
        case CtxVec::A_FFT_2K:
            return *A_fft_2k;
        case CtxVec::INV_A_I:
            return *inv_A_i;
        case CtxVec::N1:
            return *vec1_n;
        case CtxVec::N2:
            return *vec2_n;
        case CtxVec::V2K1:
            return *vec1_2k;
        case CtxVec::V2K2:
            return *vec2_2k;
        }
    }

    vec::Poly<T>& get_poly(CtxPoly type) const
    {
        switch (type) {
        case CtxPoly::A:
            return *A;
        case CtxPoly::S:
            return *S;
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

    void dump() const
    {
        std::cout << "Dump context:\n";
        std::cout << "A:";
        A->dump();
        std::cout << "inv_A_i init:";
        inv_A_i->dump();
        std::cout << "inv_A_i:";
        inv_A_i->dump();
    }

  private:
    void init(const vec::Vector<T>& vx)
    {
        // compute A(x) = prod_j(x-x_j)
        if (gf->isNF4) {
            A->set(0, this->gf->replicate(1));
        } else {
            A->set(0, 1);
        }

        for (int i = 0; i < k; ++i) {
            A->mul_to_x_plus_coef(this->gf->sub(0, vx.get(i)));
        }

        // compute A'(x) since A_i(x_i) = A'_i(x_i)
        vec::Poly<T> _A(*A);
        if (gf->isNF4) {
            _A.derivative_nf4();
        } else {
            _A.derivative();
        }

        // compute A_i(x_i)
        this->fft->fft(inv_A_i.get(), &_A);

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
            vec::ZeroExtended<T> A_2k(A.get(), len_2k);
            this->fft_2k->fft(A_fft_2k.get(), &A_2k);
        }
    }

  public:
    int vx_zero;

  private:
    unsigned k;
    unsigned n;
    unsigned len_2k;
    size_t size;
    const gf::Field<T>* gf;
    fft::FourierTransform<T>* fft;
    fft::FourierTransform<T>* fft_2k;

    const vec::Vector<T>* fragments_ids;

    std::unique_ptr<vec::Poly<T>> A = nullptr;
    std::unique_ptr<vec::Vector<T>> A_fft_2k = nullptr;
    std::unique_ptr<vec::Vector<T>> inv_A_i = nullptr;
    std::unique_ptr<vec::Poly<T>> S = nullptr;

    std::unique_ptr<vec::Vector<T>> vec1_n = nullptr;
    std::unique_ptr<vec::Vector<T>> vec2_n = nullptr;
    std::unique_ptr<vec::Vector<T>> vec1_2k = nullptr;
    std::unique_ptr<vec::Vector<T>> vec2_2k = nullptr;
};

} // namespace fec
} // namespace nttec

#endif
