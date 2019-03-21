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
#ifndef __QUAD_FEC_CONTEXT_H__
#define __QUAD_FEC_CONTEXT_H__

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <vector>
#include <sys/time.h>

#include "fft_base.h"
#include "gf_base.h"
#include "gf_nf4.h"
#include "property.h"
#include "vec_poly.h"
#include "vec_zero_ext.h"

namespace quadiron {

namespace fec {

/* List of vectors handled by context */
enum class CtxVec { A_FFT_2K = 0, INV_A_I, N1, N2, V2K1, V2K2 };

/* List of polynomials handled by context */
enum class CtxPoly { A = 0, S };

/* List of buffers handled by context */
enum class CtxBuf { N1 = 0, N2, K1, B2K1, B2K2 };

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
        std::vector<Properties>& input_props,
        const vec::Vector<T>& vx,
        const int k,
        const int n,
        int vx_zero = -1,
        const size_t size = 0,
        vec::Buffers<T>* output = nullptr)
        : props_indices(input_props.size(), 0)
    {
        this->k = k;
        this->n = n;
        this->size = size;
        this->gf = &gf;
        this->fft = &fft;
        this->fft_2k = &fft_2k;
        this->vx_zero = vx_zero;

        this->len_2k = this->gf->get_code_len_high_compo(2 * this->k);
        this->max_n_2k = (this->n > this->len_2k) ? this->n : this->len_2k;

        this->fragments_ids = &fragments_ids;

        for (auto& props : input_props) {
            // Sort properties on the basis of location of pairs in ascending
            // order.
            props.sort();
        }

        A = std::make_unique<vec::Poly<T>>(gf, n);
        A_fft_2k = std::make_unique<vec::Vector<T>>(gf, len_2k);
        inv_A_i = std::make_unique<vec::Vector<T>>(gf, k);
        S = std::make_unique<vec::Poly<T>>(gf, k);

        // zero-out all polynomials as they are used in full-length for FFT
        A->zero_fill();
        S->zero_fill();

        if (this->size == 0) {
            vec1_n = std::make_unique<vec::Vector<T>>(gf, n);
            vec2_n = std::make_unique<vec::Vector<T>>(gf, n);

            vec1_2k = std::make_unique<vec::Vector<T>>(gf, len_2k);
            vec2_2k = std::make_unique<vec::Vector<T>>(gf, len_2k);

            vec1_n->zero_fill();
            vec2_n->zero_fill();
        } else {
            assert(output != nullptr);

            // Buffers each of which is fully allocated
            // Buffer of length `len_2k`
            buf1_2k = std::make_unique<vec::Buffers<T>>(len_2k, size);
            // Buffer of length `max_n_2k - k`
            bNmK = std::make_unique<vec::Buffers<T>>(max_n_2k - k, size);

            // Buffers that are derived from the two above ones
            // Buffer sliced from `k` first elements of `buf1_2k`
            buf1_k = std::make_unique<vec::Buffers<T>>(*buf1_2k, 0, k);
            // An `n`-length buffer that is zero-extended and shuffled from
            // `buf1_k`
            buf1_n =
                std::make_unique<vec::Buffers<T>>(*buf1_k, fragments_ids, n);
            // A `max_n_2k`-length buffer combined from `output` and `bNmK`
            buf_max_n_2k = std::make_unique<vec::Buffers<T>>(*output, *bNmK);
            // An `n`-length buffer sliced from `buf_max_n_2k`
            buf2_n = std::make_unique<vec::Buffers<T>>(*buf_max_n_2k, 0, n);
            // An `len_2k`-length buffer sliced from `buf_max_n_2k`
            buf2_2k =
                std::make_unique<vec::Buffers<T>>(*buf_max_n_2k, 0, len_2k);
        }
        init(vx);
    }

    ~DecodeContext() = default;

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
        // To quell an overzealous `Wreturn-type` from GCC.
        default:
            throw InvalidArgument("invalid Vec type");
        }
    }

    vec::Poly<T>& get_poly(CtxPoly type) const
    {
        switch (type) {
        case CtxPoly::A:
            return *A;
        case CtxPoly::S:
            return *S;
        // To quell an overzealous `Wreturn-type` from GCC.
        default:
            throw InvalidArgument("invalid Poly type");
        }
    }
    vec::Buffers<T>& get_buffer(CtxBuf type) const
    {
        switch (type) {
        case CtxBuf::N1:
            return *buf1_n;
        case CtxBuf::N2:
            return *buf2_n;
        case CtxBuf::K1:
            return *buf1_k;
        case CtxBuf::B2K1:
            return *buf1_2k;
        case CtxBuf::B2K2:
            return *buf2_2k;
        // To quell an overzealous `Wreturn-type` from GCC.
        default:
            throw InvalidArgument("invalid Buf type");
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

        for (unsigned i = 0; i < k; ++i) {
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
        vec::Vector<T> _A_fft(*gf, n);
        this->fft->fft(_A_fft, _A);

        // compute 1/(x_i * A_i(x_i))
        // we care only about elements corresponding to fragments_ids
        for (int i = 0; i < static_cast<int>(k); ++i) {
            unsigned j = fragments_ids->get(i);
            if (i != vx_zero) {
                inv_A_i->set(
                    i, this->gf->inv(this->gf->mul(_A_fft.get(j), vx.get(i))));
            } else {
                inv_A_i->set(i, this->gf->inv(_A_fft.get(j)));
            }
        }

        // compute FFT(A) of length 2k
        if (this->fft_2k) {
            vec::ZeroExtended<T> A_2k(*A, len_2k);
            this->fft_2k->fft(*A_fft_2k, A_2k);
        }
    }

  public:
    int vx_zero;
    std::vector<size_t> props_indices;

  private:
    unsigned k;
    unsigned n;
    unsigned len_2k;
    unsigned max_n_2k;
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

    // Buffers each of which is fully allocated
    // Buffer of length `len_2k`
    std::unique_ptr<vec::Buffers<T>> buf1_2k = nullptr;
    // Buffer of length `max_n_2k - k`
    std::unique_ptr<vec::Buffers<T>> bNmK = nullptr;

    // Buffers that are derived from the two above ones
    // Buffer sliced from `k` first elements of `buf1_2k`
    std::unique_ptr<vec::Buffers<T>> buf1_k = nullptr;
    // An `n`-length buffer shuffled from a zero-extended of `k`-length buffer
    std::unique_ptr<vec::Buffers<T>> buf1_n = nullptr;
    // A `max_n_2k`-length buffer combined from `output` and `bNmK`
    std::unique_ptr<vec::Buffers<T>> buf_max_n_2k = nullptr;
    // An `n`-length buffer sliced from `buf_max_n_2k`
    std::unique_ptr<vec::Buffers<T>> buf2_n = nullptr;
    // An `len_2k`-length buffer sliced from `buf_max_n_2k`
    std::unique_ptr<vec::Buffers<T>> buf2_2k = nullptr;
};

} // namespace fec
} // namespace quadiron

#endif
