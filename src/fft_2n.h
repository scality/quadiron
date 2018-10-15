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
#ifndef __QUAD_FFT_2N_H__
#define __QUAD_FFT_2N_H__

#include "arith.h"
#include "fft_2.h"
#include "fft_base.h"
#include "fft_single.h"
#include "gf_base.h"
#include "vec_vector.h"
#include "vec_zero_ext.h"

/** Compute bit-reversed number for a given number
 *
 * @param x - input number
 * @param len - upper bound of the number's range. It is a power of 2
 * @param bit_len - log2 of len
 * @return - bit-reversed number of `x`
 */
inline unsigned reverse_bitwise(unsigned x, unsigned len, unsigned bit_len)
{
    unsigned mask = len;
    unsigned res = 0;
    unsigned bit_num = bit_len;
    while (mask > 0) {
        unsigned char bit = (x & mask) >> bit_num;
        res |= (bit << (bit_len - 1 - bit_num)); // NOLINT
        --bit_num;
        mask >>= 1;
    }
    return res;
}

namespace quadiron {
namespace fft {

/** Implementation of the radix-2 FFT
 *
 * It uses bit-reversal permutation algorithm that is originally described in
 * Algorithm 9.5.5 in \cite primenumbers
 */
template <typename T>
class Radix2 : public FourierTransform<T> {
  public:
    Radix2(
        const gf::Field<T>& gf,
        int n,
        int data_len = 0,
        size_t pkt_size = 0);
    ~Radix2() = default;
    void fft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void ifft(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void fft_inv(vec::Vector<T>& output, vec::Vector<T>& input) override;
    void fft(vec::Buffers<T>& output, vec::Buffers<T>& input) override;
    void ifft(vec::Buffers<T>& output, vec::Buffers<T>& input) override;
    void fft_inv(vec::Buffers<T>& output, vec::Buffers<T>& input) override;

  private:
    void init_bitrev();
    void bit_rev_permute(vec::Vector<T>& vec);
    void bit_rev_permute(vec::Buffers<T>& vec);
    void butterfly_ct_step(
        vec::Buffers<T>& buf,
        T r,
        unsigned start,
        unsigned m,
        unsigned step);
    void butterfly_ct_two_layers_step(
        vec::Buffers<T>& buf,
        unsigned start,
        unsigned m);
    void butterfly_gs_step(
        vec::Buffers<T>& buf,
        T r,
        unsigned start,
        unsigned m,
        unsigned step);
    void butterfly_gs_step_simple(
        vec::Buffers<T>& buf,
        T r,
        unsigned start,
        unsigned m,
        unsigned step);

    // Only used for non-vectorized elements
    void butterfly_ct_step_slow(
        vec::Buffers<T>& buf,
        T coef,
        unsigned start,
        unsigned m,
        unsigned step,
        size_t offset = 0);
    void butterfly_gs_step_slow(
        vec::Buffers<T>& buf,
        T coef,
        unsigned start,
        unsigned m,
        unsigned step,
        size_t offset = 0);
    void butterfly_gs_step_simple_slow(
        vec::Buffers<T>& buf,
        T coef,
        unsigned start,
        unsigned m,
        unsigned step,
        size_t offset = 0);

    unsigned data_len; // number of real input elements
    T card;
    T card_minus_one;
    T w;
    T inv_w;
    size_t pkt_size;
    size_t buf_size;

    std::unique_ptr<T[]> rev = nullptr;
    std::unique_ptr<vec::Vector<T>> W = nullptr;
    std::unique_ptr<vec::Vector<T>> inv_W = nullptr;
    T* vec_W;
};

/**
 * n-th root will be constructed with primitive root
 *
 * @param gf
 * @param n FFT length, for now must be a power of 2
 * @param m length of input vector without zero padding. It allows shorterning
 *  operation cycles
 * @param pkt_size size of packet, i.e. number of symbols per chunk will be
 *  received and processed at a time
 *
 * @return
 */
template <typename T>
Radix2<T>::Radix2(const gf::Field<T>& gf, int n, int data_len, size_t pkt_size)
    : FourierTransform<T>(gf, n)
{
    w = gf.get_nth_root(n);
    inv_w = gf.inv(w);
    this->pkt_size = pkt_size;
    this->data_len = data_len > 0 ? data_len : n;
    buf_size = pkt_size * sizeof(T);

    W = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, n));
    inv_W = std::unique_ptr<vec::Vector<T>>(new vec::Vector<T>(gf, n));
    gf.compute_omegas(*W, n, w);
    gf.compute_omegas(*inv_W, n, inv_w);

    vec_W = W->get_mem();

    card = this->gf->card();
    card_minus_one = this->gf->card_minus_one();

    rev = std::unique_ptr<T[]>(new T[n]);
    init_bitrev();
}

template <typename T>
void Radix2<T>::init_bitrev()
{
    unsigned len = this->n;
    unsigned log_n = arith::log2<T>(len);

    // init bit-reversed indices
    for (unsigned i = 0; i < len; ++i) {
        rev[i] = reverse_bitwise(i, len, log_n);
    }
}

template <typename T>
void Radix2<T>::bit_rev_permute(vec::Vector<T>& vec)
{
    for (unsigned i = 0; i < (unsigned)this->n; ++i) {
        if (rev[i] < i) {
            vec.swap(i, rev[i]);
        }
    }
}

template <typename T>
void Radix2<T>::bit_rev_permute(vec::Buffers<T>& vec)
{
    for (unsigned i = 0; i < (unsigned)this->n; ++i) {
        if (rev[i] < i) {
            vec.swap(i, rev[i]);
        }
    }
}

/** Perform decimation-in-time FFT
 *
 * Process:
 * - Output is initialized as the input in the bit-reversal ordering
 * - Cooley-Tukey butterfly is performed on output
 * It leads to a O(N*logN) complexity.
 * @note: we consider particular cases where input is padded by zeros that
 * always happens in erasure encoding process. Normally, FFT process starts
 * performing butterfly operations on groups of size 2, then 4, 8 etc.
 * Thanks to zero padding, we can start from a group size where there is at
 * most a non-zero element of input in the group, all other elements of the
 * group are zero. Indeed, after performing butterfly operation on the
 * group, all elements equal to the non-zero one.
 *
 * Given FFT length N, number of data length K, output is initialized in
 * the following clever way:
 * - Each group contains N/K elements
 * - Perform on groups containing non-zero element of input: copy input's
 * element at the bit-reversed index to all elements of correspondent group
 * of output
 * - For other groups, initialize them normally
 * It leads to a O(K*logN) complexity
 *
 * @param output - output vector
 * @param input - input vector
 */
template <typename T>
void Radix2<T>::fft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    const unsigned len = this->n;
    const unsigned input_len = input.get_n();
    // to support FFT on input vectors of length greater than from `data_len`
    const unsigned group_len =
        (input_len > data_len) ? len / input_len : len / data_len;

    for (unsigned idx = 0; idx < input_len; ++idx) {
        // set output  = scramble(input), i.e. bit reversal ordering
        const T a = input.get(idx);
        const unsigned end = rev[idx] + group_len;
        for (unsigned i = rev[idx]; i < end; ++i) {
            output.set(i, a);
        }
    }
    for (unsigned idx = input_len; idx < data_len; ++idx) {
        // set output  = scramble(input), i.e. bit reversal ordering
        const unsigned end = rev[idx] + group_len;
        for (unsigned i = rev[idx]; i < end; ++i) {
            output.set(i, 0);
        }
    }
    // perform butterfly operations
    for (unsigned m = group_len; m < len; m *= 2) {
        const unsigned doubled_m = 2 * m;
        const unsigned ratio = len / doubled_m;
        for (unsigned j = 0; j < m; ++j) {
            const T r = W->get(j * ratio);
            for (unsigned i = j; i < len; i += doubled_m) {
                const T a = output.get(i);
                const T b = this->gf->mul(r, output.get(i + m));
                output.set(i, this->gf->add(a, b));
                output.set(i + m, this->gf->sub(a, b));
            }
        }
    }
}

/** Perform decimation-in-frequency FFT or inverse FFT
 *
 * @note: input vector is in reversed-bit order
 *
 * Process:
 * - Genteleman-Sande butterfly is performed on input
 * - Output is copied from input at bit-reversed indices
 * It leads to a O(N*logN) complexity.
 *
 * @param output - output vector
 * @param input - input vector
 */
template <typename T>
void Radix2<T>::fft_inv(vec::Vector<T>& output, vec::Vector<T>& input)
{
    const unsigned len = this->n;

    output.copy(&input);

    for (unsigned m = len / 2; m >= 1; m /= 2) {
        unsigned doubled_m = 2 * m;
        for (unsigned j = 0; j < m; ++j) {
            T r = inv_W->get(j * len / doubled_m);
            for (unsigned i = j; i < len; i += doubled_m) {
                const T a = output.get(i);
                const T b = output.get(i + m);
                output.set(i, this->gf->add(a, b));
                output.set(i + m, this->gf->mul(r, this->gf->sub(a, b)));
            }
        }
    }

    // reversion of elements of output to return values on the natural order
    bit_rev_permute(output);
}

template <typename T>
void Radix2<T>::ifft(vec::Vector<T>& output, vec::Vector<T>& input)
{
    fft_inv(output, input);

    /*
     * We need to divide output to `N` for the inverse formular
     */
    output.mul_scalar(this->inv_n_mod_p);
}

/** Perform decimation-in-time FFT
 *
 * @param output - output buffers
 * @param input - input buffers
 */
template <typename T>
void Radix2<T>::fft(vec::Buffers<T>& output, vec::Buffers<T>& input)
{
    const unsigned len = this->n;
    const unsigned input_len = input.get_n();
    // to support FFT on input vectors of length greater than from `data_len`
    const unsigned group_len =
        (input_len > data_len) ? len / input_len : len / data_len;

    const std::vector<T*>& i_mem = input.get_mem();
    const std::vector<T*>& o_mem = output.get_mem();

    for (unsigned idx = 0; idx < input_len; ++idx) {
        // set output  = scramble(input), i.e. bit reversal ordering
        for (unsigned i = rev[idx]; i < rev[idx] + group_len; ++i) {
            memcpy(o_mem[i], i_mem[idx], buf_size);
        }
    }
    for (unsigned idx = input_len; idx < data_len; ++idx) {
        // set output  = scramble(input), i.e. bit reversal ordering
        for (unsigned i = rev[idx]; i < rev[idx] + group_len; ++i) {
            memset(o_mem[i], 0, buf_size);
        }
    }

    // ----------------------
    // Two layers at a time
    // ----------------------
    unsigned m = group_len;
    const unsigned end = len / 2;
    for (; m < end; m <<= 2) {
        for (unsigned j = 0; j < m; ++j) {
            butterfly_ct_two_layers_step(output, j, m);
        }
    }
    if (m < len) {
        assert(m == end);
        // perform the last butterfly operations
        for (unsigned j = 0; j < m; ++j) {
            const T r = W->get(j);
            butterfly_ct_step(output, r, j, m, len);
        }
    }
}

// for each pair (P, Q) = (buf[i], buf[i + m]):
// P = P + c * Q
// Q = P - c * Q
template <typename T>
void Radix2<T>::butterfly_ct_step(
    vec::Buffers<T>& buf,
    T r,
    unsigned start,
    unsigned m,
    unsigned step)
{
    butterfly_ct_step_slow(buf, r, start, m, step);
}

/**
 * Butterly CT on two-layers at a time
 *
 * For each quadruple
 * (P, Q, R, S) = (buf[i], buf[i + m], buf[i + 2 * m], buf[i + 3 * m])
 * First layer: butterfly on (P, Q) and (R, S) for step = 2 * m
 *      coef r1 = W[start * n / (2 * m)]
 *      P = P + r1 * Q
 *      Q = P - r1 * Q
 *      R = R + r1 * S
 *      S = R - r1 * S
 * Second layer: butterfly on (P, R) and (Q, S) for step = 4 * m
 *      coef r2 = W[start * n / (4 * m)]
 *      coef r3 = W[(start + m) * n / (4 * m)]
 *      P = P + r2 * R
 *      R = P - r2 * R
 *      Q = Q + r3 * S
 *      S = Q - r3 * S
 *
 * @param buf - working buffers
 * @param start - index of buffer among `m` ones
 * @param m - current group size
 */
template <typename T>
void Radix2<T>::butterfly_ct_two_layers_step(
    vec::Buffers<T>& buf,
    unsigned start,
    unsigned m)
{
    const unsigned step = m << 2;
    //  ---------
    // First layer
    //  ---------
    const T r1 = W->get(start * this->n / m / 2);
    // first pair
    butterfly_ct_step(buf, r1, start, m, step);
    // second pair
    butterfly_ct_step(buf, r1, start + 2 * m, m, step);
    //  ---------
    // Second layer
    //  ---------
    // first pair
    const T r2 = W->get(start * this->n / m / 4);
    butterfly_ct_step(buf, r2, start, 2 * m, step);
    // second pair
    const T r3 = W->get((start + m) * this->n / m / 4);
    butterfly_ct_step(buf, r3, start + m, 2 * m, step);
}

template <typename T>
void Radix2<T>::butterfly_ct_step_slow(
    vec::Buffers<T>& buf,
    T coef,
    unsigned start,
    unsigned m,
    unsigned step,
    size_t offset)
{
    for (int i = start; i < this->n; i += step) {
        T* a = buf.get(i);
        T* b = buf.get(i + m);
        // perform butterfly operation for Cooley-Tukey FFT algorithm
        for (size_t j = offset; j < this->pkt_size; ++j) {
            T x = this->gf->mul(coef, b[j]);
            b[j] = this->gf->sub(a[j], x);
            a[j] = this->gf->add(a[j], x);
        }
    }
}

/** Perform decimation-in-frequency FFT or inverse FFT
 *
 * Input buffer is in reversed-bit order. Hence butterfly operations can
 * be performed directly on this order. But the output must be re-order in the
 * natural order by using swap operations. A swap operation is to swap only two
 * buffers' pointer, hence it's not costly. However, it will modify the output
 * buffers. To avoid that, we reverse the output twice, before and after the
 * butterfly operation.
 *
 * @param output - output buffers
 * @param input - input buffers
 */
template <typename T>
void Radix2<T>::fft_inv(vec::Buffers<T>& output, vec::Buffers<T>& input)
{
    const unsigned len = this->n;
    const unsigned input_len = input.get_n();

    // 1st reversion of elements of output
    bit_rev_permute(output);

    // copy input to output
    const std::vector<T*>& i_mem = input.get_mem();
    const std::vector<T*>& o_mem = output.get_mem();
    unsigned i;
    for (i = 0; i < input_len; ++i) {
        memcpy(o_mem[i], i_mem[i], buf_size);
    }

    unsigned m = len / 2;

    if (input_len < len) {
        const unsigned input_len_power_2 =
            arith::get_smallest_power_of_2<unsigned>(input_len);
        for (; i < input_len_power_2; ++i) {
            memset(o_mem[i], 0, buf_size);
        }

        // For Q are zeros only => Q = c * P
        for (; m >= input_len; m /= 2) {
            unsigned doubled_m = 2 * m;
            for (unsigned j = 0; j < m; ++j) {
                const T r = inv_W->get(j * len / doubled_m);
                butterfly_gs_step_simple(output, r, j, m, doubled_m);
            }
        }
    }

    // Next, normal butterlfy GS is performed
    for (; m >= 1; m /= 2) {
        unsigned doubled_m = 2 * m;
        for (unsigned j = 0; j < m; ++j) {
            const T r = inv_W->get(j * len / doubled_m);
            butterfly_gs_step(output, r, j, m, doubled_m);
        }
    }

    // 2nd reversion of elements of output to return its natural order
    bit_rev_permute(output);
}

// for each pair (P, Q) = (buf[i], buf[i + m]):
// Q = c * P
template <typename T>
void Radix2<T>::butterfly_gs_step_simple(
    vec::Buffers<T>& buf,
    T coef,
    unsigned start,
    unsigned m,
    unsigned step)
{
    butterfly_gs_step_simple_slow(buf, coef, start, m, step);
}

// for each pair (P, Q) = (buf[i], buf[i + m]):
// P = P + Q
// Q = c * (P - Q)
template <typename T>
void Radix2<T>::butterfly_gs_step(
    vec::Buffers<T>& buf,
    T coef,
    unsigned start,
    unsigned m,
    unsigned step)
{
    butterfly_gs_step_slow(buf, coef, start, m, step);
}

template <typename T>
void Radix2<T>::butterfly_gs_step_slow(
    vec::Buffers<T>& buf,
    T coef,
    unsigned start,
    unsigned m,
    unsigned step,
    size_t offset)
{
    for (int i = start; i < this->n; i += step) {
        T* a = buf.get(i);
        T* b = buf.get(i + m);
        // perform butterfly operation for Cooley-Tukey FFT algorithm
        for (size_t j = offset; j < this->pkt_size; ++j) {
            T x = this->gf->sub(a[j], b[j]);
            a[j] = this->gf->add(a[j], b[j]);
            b[j] = this->gf->mul(coef, x);
        }
    }
}

template <typename T>
void Radix2<T>::butterfly_gs_step_simple_slow(
    vec::Buffers<T>& buf,
    T coef,
    unsigned start,
    unsigned m,
    unsigned step,
    size_t offset)
{
    for (int i = start; i < this->n; i += step) {
        T* a = buf.get(i);
        T* b = buf.get(i + m);
        // perform butterfly operation for Cooley-Tukey FFT algorithm
        for (size_t j = offset; j < this->pkt_size; ++j) {
            b[j] = this->gf->mul(coef, a[j]);
        }
    }
}

template <typename T>
void Radix2<T>::ifft(vec::Buffers<T>& output, vec::Buffers<T>& input)
{
    fft_inv(output, input);

    /*
     * We need to divide output to `N` for the inverse formular
     */
    this->gf->mul_vec_to_vecp(*(this->vec_inv_n), output, output);
}

} // namespace fft
} // namespace quadiron

#endif
