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
#ifndef __NTTEC_FEC_BASE_H__
#define __NTTEC_FEC_BASE_H__

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <vector>
#include <sys/time.h>

#include "fft_base.h"
#include "gf_base.h"
#include "misc.h"
#include "property.h"
#include "vec_buffers.h"
#include "vec_poly.h"
#include "vec_vector.h"

namespace nttec {

/** Forward Error Correction code implementations. */
namespace fec {

static inline timeval tick()
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return tv;
}

static inline uint64_t hrtime_usec(timeval begin)
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return 1000000 * (tv.tv_sec - begin.tv_sec) + tv.tv_usec - begin.tv_usec;
}

/*
 * Get and cast mem of Buffers<Ts> to a vector of Td*
 */
template <typename Ts, typename Td>
std::vector<Td*>* cast_mem_of_vecp(vec::Buffers<Ts>* s)
{
    int i;
    int n = s->get_n();

    // std::cout << "\ninput: "; s->dump();

    std::vector<Ts*>* mem_s = s->get_mem();
    std::vector<Td*>* mem_d = new std::vector<Td*>(n, nullptr);
    for (i = 0; i < n; i++) {
        mem_d->at(i) = static_cast<Td*>(static_cast<void*>(mem_s->at(i)));
    }

    return mem_d;
}

enum class FecType {
    /** Systematic code
     *
     * Take n_data input, generate n_parities outputs.
     */
    SYSTEMATIC,

    /** Non-systematic code
     *
     * Take n_data input, generate n_data+n_parities outputs.
     */
    NON_SYSTEMATIC
};

/** Base class for Forward Error Correction (FEC) codes. */
template <typename T>
class FecCode {
  public:
    FecType type;
    unsigned word_size;
    unsigned n_data;
    unsigned n_parities;
    unsigned code_len;
    unsigned n_outputs;
    size_t pkt_size; // packet size, i.e. number of words per packet
    size_t buf_size; // packet size in bytes

    // Length of operating codeword. It's calculated by derived Class
    // FIXME: move n to protected
    T n;

    uint64_t total_encode_cycles = 0;
    uint64_t n_encode_ops = 0;
    uint64_t total_decode_cycles = 0;
    uint64_t n_decode_ops = 0;

    uint64_t total_enc_usec = 0;
    uint64_t total_dec_usec = 0;

    FecCode(
        FecType type,
        unsigned word_size,
        unsigned n_data,
        unsigned n_parities,
        size_t pkt_size = 8);
    virtual ~FecCode();

    /**
     * Return the number actual parities for SYSTEMATIC it is exactly
     * n_parities, for NON_SYSTEMATIC it maybe at least n_data+n_parities (but
     * sometimes more).
     * @return
     */
    virtual int get_n_outputs() = 0;

    virtual void encode(
        vec::Vector<T>* output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* words) = 0;
    virtual void encode(
        vec::Buffers<T>* output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>* words){};
    virtual void decode_add_data(int fragment_index, int row) = 0;
    virtual void decode_add_parities(int fragment_index, int row) = 0;
    virtual void decode_build(void) = 0;

    /**
     * Decode a vector of words
     *
     * @param props properties bound to parity fragments
     * @param offset offset in the data fragments
     * @param output original data (must be of n_data length)
     * @param fragments_ids identifiers of available fragments
     * @param words input words if SYSTEMATIC must be n_data,
     *  if NON_SYSTEMATIC get_n_outputs()
     */
    virtual void decode(
        vec::Vector<T>* output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words);

    bool readw(T* ptr, std::istream* stream);
    bool writew(T val, std::ostream* stream);

    bool read_pkt(char* pkt, std::istream* stream);
    bool write_pkt(char* pkt, std::ostream* stream);

    void encode_bufs(
        std::vector<std::istream*> input_data_bufs,
        std::vector<std::ostream*> output_parities_bufs,
        std::vector<Properties>& output_parities_props);

    void encode_packet(
        std::vector<std::istream*> input_data_bufs,
        std::vector<std::ostream*> output_parities_bufs,
        std::vector<Properties>& output_parities_props);

    bool decode_bufs(
        std::vector<std::istream*> input_data_bufs,
        std::vector<std::istream*> input_parities_bufs,
        const std::vector<Properties>& input_parities_props,
        std::vector<std::ostream*> output_data_bufs);

    gf::Field<T>* get_gf()
    {
        return this->gf;
    }

    void reset_stats_enc()
    {
        total_encode_cycles = 0;
        n_encode_ops = 0;
        total_enc_usec = 0;
    }

    void reset_stats_dec()
    {
        total_decode_cycles = 0;
        n_decode_ops = 0;
        total_dec_usec = 0;
    }

  protected:
    // primitive nth root of unity
    T r;
    gf::Field<T>* gf = nullptr;
    std::unique_ptr<fft::FourierTransform<T>> fft = nullptr;
    std::unique_ptr<fft::FourierTransform<T>> fft_full = nullptr;
    // This vector MUST be initialized by derived Class using multiplicative FFT
    std::unique_ptr<vec::Vector<T>> inv_r_powers = nullptr;

    // pure abstract methods that will be defined in derived class
    virtual void check_params() = 0;
    virtual void init_gf() = 0;
    virtual void init_fft() = 0;
    virtual void init_others() = 0;

    // This function will called in constructor of every derived class
    void fec_init()
    {
        // check compatible parameters
        check_params();
        // create the field
        init_gf();
        // create FFT
        init_fft();
        // init other parameters dedicated for derived class
        init_others();
    }

    virtual void decode_prepare(
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words,
        vec::Vector<T>* vx,
        int* vx_zero = 0);

    virtual void decode_lagrange(
        vec::Vector<T>* output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words,
        vec::Vector<T>* vx,
        int vx_zero = 0);

    virtual void decode_vec_lagrange(
        vec::Vector<T>* output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words,
        vec::Vector<T>* vx,
        int vx_zero = 0);
};

/**
 * Create an encoder
 *
 * @param word_size in bytes
 * @param n_data
 * @param n_parities
 */
template <typename T>
FecCode<T>::FecCode(
    FecType type,
    unsigned word_size,
    unsigned n_data,
    unsigned n_parities,
    size_t pkt_size)
{
    assert(type == FecType::SYSTEMATIC || type == FecType::NON_SYSTEMATIC);

    this->type = type;
    this->word_size = word_size;
    this->n_data = n_data;
    this->n_parities = n_parities;
    this->code_len = n_data + n_parities;
    this->n_outputs =
        (type == FecType::SYSTEMATIC) ? this->n_parities : this->code_len;
    this->pkt_size = pkt_size;
    this->buf_size = pkt_size * word_size;
}

template <typename T>
FecCode<T>::~FecCode()
{
}

template <typename T>
inline bool FecCode<T>::readw(T* ptr, std::istream* stream)
{
    if (word_size == 1) {
        uint8_t c;
        if (stream->read((char*)&c, 1)) {
            *ptr = c;
            return true;
        }
    } else if (word_size == 2) {
        uint16_t s;
        if (stream->read((char*)&s, 2)) {
            *ptr = s;
            return true;
        }
    } else if (word_size == 4) {
        uint32_t s;
        if (stream->read((char*)&s, 4)) {
            *ptr = s;
            return true;
        }
    } else if (word_size == 8) {
        uint64_t s;
        if (stream->read((char*)&s, 8)) {
            *ptr = s;
            return true;
        }
    } else if (word_size == 16) {
        __uint128_t s;
        if (stream->read((char*)&s, 16)) {
            *ptr = s;
            return true;
        }
    } else {
        assert(false && "no such size");
    }
    return false;
}

template <typename T>
inline bool FecCode<T>::writew(T val, std::ostream* stream)
{
    if (word_size == 1) {
        uint8_t c = val;
        if (stream->write((char*)&c, 1))
            return true;
    } else if (word_size == 2) {
        uint16_t s = val;
        if (stream->write((char*)&s, 2))
            return true;
    } else if (word_size == 4) {
        uint32_t s = val;
        if (stream->write((char*)&s, 4))
            return true;
    } else if (word_size == 8) {
        uint64_t s = val;
        if (stream->write((char*)&s, 8))
            return true;
    } else if (word_size == 16) {
        __uint128_t s = val;
        if (stream->write((char*)&s, 16))
            return true;
    } else {
        assert(false && "no such size");
    }
    return false;
}

template <typename T>
inline bool FecCode<T>::read_pkt(char* pkt, std::istream* stream)
{
    if (stream->read(pkt, buf_size)) {
        return true;
    }
    return false;
}

template <typename T>
inline bool FecCode<T>::write_pkt(char* pkt, std::ostream* stream)
{
    if (stream->write(pkt, buf_size))
        return true;
    return false;
}

/**
 * Encode buffers
 *
 * @param input_data_bufs must be exactly n_data
 * @param output_parities_bufs must be exactly get_n_outputs() (set nullptr when
 * not missing/wanted)
 * @param output_parities_props must be exactly get_n_outputs() specific
 * properties that the called is supposed to store along with parities
 *
 * @note all streams must be of equal size
 */
template <typename T>
void FecCode<T>::encode_bufs(
    std::vector<std::istream*> input_data_bufs,
    std::vector<std::ostream*> output_parities_bufs,
    std::vector<Properties>& output_parities_props)
{
    bool cont = true;
    off_t offset = 0;

    assert(input_data_bufs.size() == n_data);
    assert(output_parities_bufs.size() == n_outputs);
    assert(output_parities_props.size() == n_outputs);

    vec::Vector<T> words(gf, n_data);
    vec::Vector<T> output(gf, get_n_outputs());

    reset_stats_enc();

    while (true) {
        words.zero_fill();
        for (unsigned i = 0; i < n_data; i++) {
            T tmp;
            if (!readw(&tmp, input_data_bufs[i])) {
                cont = false;
                break;
            }
            words.set(i, tmp);
        }
        if (!cont)
            break;

        // std::cout << "words at " << offset << ": "; words.dump();

        timeval t1 = tick();
        uint64_t start = rdtsc();
        encode(&output, output_parities_props, offset, &words);
        uint64_t end = rdtsc();
        uint64_t t2 = hrtime_usec(t1);

        // std::cout << "output: "; output.dump();

        total_enc_usec += t2;
        total_encode_cycles += (end - start) / word_size;
        n_encode_ops++;

        for (unsigned i = 0; i < n_outputs; i++) {
            T tmp = output.get(i);
            writew(tmp, output_parities_bufs[i]);
        }
        offset += word_size;
    }
}

template <typename T>
void FecCode<T>::encode_packet(
    std::vector<std::istream*> input_data_bufs,
    std::vector<std::ostream*> output_parities_bufs,
    std::vector<Properties>& output_parities_props)
{
    assert(input_data_bufs.size() == n_data);
    assert(output_parities_bufs.size() == n_outputs);
    assert(output_parities_props.size() == n_outputs);

    bool cont = true;
    off_t offset = 0;
    bool full_word_size = (word_size == sizeof(T));

    vec::Buffers<uint8_t> words_char(n_data, buf_size);
    std::vector<uint8_t*>* words_mem_char = words_char.get_mem();
    std::vector<T*>* words_mem_T = nullptr;
    if (full_word_size)
        words_mem_T = cast_mem_of_vecp<uint8_t, T>(&words_char);
    vec::Buffers<T> words(n_data, pkt_size, words_mem_T);
    words_mem_T = words.get_mem();

    int output_len = get_n_outputs();

    vec::Buffers<T> output(output_len, pkt_size);
    std::vector<T*>* output_mem_T = output.get_mem();
    std::vector<uint8_t*>* output_mem_char = nullptr;
    if (full_word_size)
        output_mem_char = cast_mem_of_vecp<T, uint8_t>(&output);
    vec::Buffers<uint8_t> output_char(output_len, buf_size, output_mem_char);
    output_mem_char = output_char.get_mem();

    reset_stats_enc();

    while (true) {
        // TODO: get number of read bytes -> true buf size
        // words.zero_fill();
        for (unsigned i = 0; i < n_data; i++) {
            if (!read_pkt((char*)(words_mem_char->at(i)), input_data_bufs[i])) {
                cont = false;
                break;
            }
        }
        if (!cont)
            break;

        if (!full_word_size)
            vec::pack<uint8_t, T>(
                words_mem_char, words_mem_T, n_data, pkt_size, word_size);

        timeval t1 = tick();
        uint64_t start = rdtsc();
        encode(&output, output_parities_props, offset, &words);
        uint64_t end = rdtsc();
        uint64_t t2 = hrtime_usec(t1);

        total_enc_usec += t2;
        total_encode_cycles += (end - start) / buf_size;
        n_encode_ops++;

        if (!full_word_size)
            vec::unpack<T, uint8_t>(
                output_mem_T, output_mem_char, output_len, pkt_size, word_size);

        for (unsigned i = 0; i < n_outputs; i++) {
            write_pkt((char*)(output_mem_char->at(i)), output_parities_bufs[i]);
        }
        offset += buf_size;
    }
}

/**
 * Decode buffers
 *
 * @param input_data_bufs if SYSTEMATIC must be exactly n_data otherwise it is
 * unused (use nullptr when missing)
 * @param input_parities_bufs if SYSTEMATIC must be exactly n_parities otherwise
 * get_n_outputs() (use nullptr when missing)
 * @param input_parities_props if SYSTEMATIC must be exactly n_parities
 * otherwise get_n_outputs() caller is supposed to provide specific information
 * bound to parities
 * @param output_data_bufs must be exactly n_data (use nullptr when not
 * missing/wanted)
 *
 * @note All streams must be of equal size
 *
 * @return true if decode succeded, else false
 */
template <typename T>
bool FecCode<T>::decode_bufs(
    std::vector<std::istream*> input_data_bufs,
    std::vector<std::istream*> input_parities_bufs,
    const std::vector<Properties>& input_parities_props,
    std::vector<std::ostream*> output_data_bufs)
{
    off_t offset = 0;
    bool cont = true;

    unsigned fragment_index = 0;

    if (type == FecType::SYSTEMATIC) {
        assert(input_data_bufs.size() == n_data);
    }
    assert(input_parities_bufs.size() == n_outputs);
    assert(input_parities_props.size() == n_outputs);
    assert(output_data_bufs.size() == n_data);

    reset_stats_dec();

    if (type == FecType::SYSTEMATIC) {
        for (unsigned i = 0; i < n_data; i++) {
            if (input_data_bufs[i] != nullptr) {
                decode_add_data(fragment_index, i);
                fragment_index++;
            }
        }
        // data is in clear so nothing to do
        if (fragment_index == n_data)
            return true;
    }

    if (fragment_index < n_data) {
        // finish with parities available
        for (unsigned i = 0; i < code_len; i++) {
            if (input_parities_bufs[i] != nullptr) {
                decode_add_parities(fragment_index, i);
                fragment_index++;
                // stop when we have enough parities
                if (fragment_index == n_data)
                    break;
            }
        }
        // unable to decode
        if (fragment_index < n_data)
            return false;
    }

    decode_build();

    int n_words = code_len;
    if (type == FecType::SYSTEMATIC) {
        n_words = n_data;
    }

    vec::Vector<T> words(gf, n_words);
    vec::Vector<T> fragments_ids(gf, n_words);
    vec::Vector<T> output(gf, n_data);

    while (true) {
        words.zero_fill();
        fragment_index = 0;
        if (type == FecType::SYSTEMATIC) {
            for (unsigned i = 0; i < n_data; i++) {
                if (input_data_bufs[i] != nullptr) {
                    T tmp;
                    if (!readw(&tmp, input_data_bufs[i])) {
                        cont = false;
                        break;
                    }
                    fragments_ids.set(fragment_index, i);
                    words.set(fragment_index, tmp);
                    fragment_index++;
                }
            }
            if (!cont)
                break;
            // stop when we have enough parities
            if (fragment_index == n_data)
                break;
        }
        for (unsigned i = 0; i < n_outputs; i++) {
            if (input_parities_bufs[i] != nullptr) {
                T tmp;
                if (!readw(&tmp, input_parities_bufs[i])) {
                    cont = false;
                    break;
                }
                fragments_ids.set(fragment_index, i);
                words.set(fragment_index, tmp);
                fragment_index++;
                // stop when we have enough parities
                if (fragment_index == n_data)
                    break;
            }
        }
        if (!cont)
            break;

        timeval t1 = tick();
        uint64_t start = rdtsc();
        decode(&output, input_parities_props, offset, &fragments_ids, &words);
        uint64_t end = rdtsc();
        uint64_t t2 = hrtime_usec(t1);

        total_dec_usec += t2;
        total_decode_cycles += (end - start) / word_size;
        n_decode_ops++;

        for (unsigned i = 0; i < n_data; i++) {
            if (output_data_bufs[i] != nullptr) {
                T tmp = output.get(i);
                writew(tmp, output_data_bufs[i]);
            }
        }

        offset += word_size;
    }

    return true;
}

/**
 * Perform a Lagrange interpolation to find the coefficients of the
 * polynomial
 *
 * @note If all fragments are available ifft(words) is enough
 *
 * @param output must be exactly n_data
 * @param props special values dictionary must be exactly n_data
 * @param offset used to locate special values
 * @param fragments_ids unused
 * @param words v=(v_0, v_1, ..., v_k-1) k must be exactly n_data
 */
template <typename T>
void FecCode<T>::decode(
    vec::Vector<T>* output,
    const std::vector<Properties>& props,
    off_t offset,
    vec::Vector<T>* fragments_ids,
    vec::Vector<T>* words)
{
    int vx_zero;
    // vector x=(x_0, x_1, ..., x_k-1)
    vec::Vector<T> vx(this->gf, this->n_data);

    // prepare for decoding
    decode_prepare(props, offset, fragments_ids, words, &vx, &vx_zero);

    // Lagrange interpolation
    decode_vec_lagrange(
        output, props, offset, fragments_ids, words, &vx, vx_zero);
}

/* Prepare for decoding
 * It supports for FEC using multiplicative FFT over FNT
 */
template <typename T>
void FecCode<T>::decode_prepare(
    const std::vector<Properties>& props,
    off_t offset,
    vec::Vector<T>* fragments_ids,
    vec::Vector<T>* words,
    vec::Vector<T>* vx,
    int* vx_zero)
{
    int k = this->n_data; // number of fragments received
    // vector x=(x_0, x_1, ..., x_k-1)
    for (int i = 0; i < k; i++) {
        vx->set(i, this->gf->exp(r, fragments_ids->get(i)));
    }

    for (int i = 0; i < k; i++) {
        const int j = fragments_ids->get(i);
        auto data = props[j].get(ValueLocation(offset, j));

        // Check if the symbol is a special case whick is marked by "@".
        // Note: this check is necessary when word_size is not large enough to
        // cover all symbols of the field.
        // Following check is used for FFT over FNT where the single special
        // case symbol equals card - 1
        if (data && *data == "@") {
            words->set(i, this->gf->card() - 1);
        }
    }
}

/* Lagrange interpolation
 * It supports for FEC using multiplicative FFT over GF differing NF4
 */
template <typename T>
void FecCode<T>::decode_lagrange(
    vec::Vector<T>* output,
    const std::vector<Properties>& props,
    off_t offset,
    vec::Vector<T>* fragments_ids,
    vec::Vector<T>* words,
    vec::Vector<T>* vx,
    int vx_zero)
{
    if (inv_r_powers == nullptr) {
        throw LogicError("The pointer 'inv_r_powers' must be initialized");
    }

    Polynomial<T> A(this->gf);
    Polynomial<T> _A(this->gf);
    Polynomial<T> N_p(this->gf);
    Polynomial<T> S(this->gf);

    int k = this->n_data; // number of fragments received

    // compute A(x) = prod_j(x-x_j)
    A.set(0, 1);
    for (int i = 0; i < k; i++) {
        A.mul_to_x_plus_coef(this->gf->sub(0, vx->get(i)));
    }

    // compute A'(x) since A_i(x_i) = A'_i(x_i)
    _A.copy(&A);
    _A.derivative();

    // evaluate n_i=v_i/A'_i(x_i)
    vec::Vector<T> _n(this->gf, k);
    for (int i = 0; i < k; i++) {
        _n.set(i, this->gf->div(words->get(i), _A.eval(vx->get(i))));
    }

    // compute N'(x) = sum_i{n_i * x^z_i}
    for (int i = 0; i <= k - 1; i++) {
        N_p.set(fragments_ids->get(i), _n.get(i));
    }

    // We have to find the numerator of the following expression:
    // P(x)/A(x) = sum_i=0_k-1(n_i/(x-x_i)) mod x^n
    // using Taylor series we rewrite the expression into
    // P(x)/A(x) = -sum_i=0_k-1(sum_j=0_n-1(n_i*x_i^(-j-1)*x^j))
    for (T i = 0; i <= k - 1; i++) {
        S.set(i, N_p.eval(inv_r_powers->get(i + 1)));
    }
    S.neg();
    S.mul(&A, k - 1);
    // No need to mod x^n since only last n_data coefs are obtained
    // output is n_data length
    for (unsigned i = 0; i < this->n_data; i++)
        output->set(i, S.get(i));
}

template <typename T>
void FecCode<T>::decode_vec_lagrange(
    vec::Vector<T>* output,
    const std::vector<Properties>& props,
    off_t offset,
    vec::Vector<T>* fragments_ids,
    vec::Vector<T>* words,
    vec::Vector<T>* vx,
    int vx_zero)
{
    if (this->fft == nullptr) {
        throw LogicError("FEC base: FFT must be initialized");
    }
    if (this->fft_full == nullptr) {
        throw LogicError("FEC base: FFT must be initialized");
    }
    int k = this->n_data; // number of fragments received

    vec::Poly<T> A(this->gf, n);
    vec::Poly<T> _A_fft(this->gf, n);
    vec::Poly<T> N_p(this->gf, n);
    vec::Poly<T> N_p_ifft(this->gf, n);
    vec::Poly<T> S(this->gf, k);

    // compute A(x) = prod_j(x-x_j)
    A.zero();
    A.set(0, 1);
    for (int i = 0; i < k; i++) {
        A.mul_to_x_plus_coef(this->gf->sub(0, vx->get(i)));
    }

    // compute A'(x) since A_i(x_i) = A'_i(x_i)
    vec::Poly<T> _A(A);
    _A.derivative();
    this->fft->fft(&_A_fft, &_A);

    // evaluate n_i=v_i/A'_i(x_i)
    vec::Vector<T> _n(this->gf, k);
    for (int i = 0; i < k; i++) {
        _n.set(
            i, this->gf->div(words->get(i), _A_fft.get(fragments_ids->get(i))));
    }

    // compute N'(x) = sum_i{n_i * x^z_i}
    N_p.zero();
    for (int i = 0; i <= k - 1; i++) {
        N_p.set(fragments_ids->get(i), _n.get(i));
    }
    this->fft_full->fft_inv(&N_p_ifft, &N_p);

    // We have to find the numerator of the following expression:
    // P(x)/A(x) = sum_i=0_k-1(n_i/(x-x_i)) mod x^n
    // using Taylor series we rewrite the expression into
    // P(x)/A(x) = -sum_i=0_k-1(sum_j=0_n-1(n_i*x_i^(-j-1)*x^j))
    for (int i = 0; i <= k - 1; i++) {
        S.set(i, N_p_ifft.get(i + 1));
    }
    S.neg();
    S.mul(&A, k - 1);
    // No need to mod x^n since only last n_data coefs are obtained
    // output is n_data length
    for (unsigned i = 0; i < this->n_data; i++)
        output->set(i, S.get(i));
}

} // namespace fec
} // namespace nttec

#endif
