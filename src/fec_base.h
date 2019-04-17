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
#ifndef __QUAD_FEC_BASE_H__
#define __QUAD_FEC_BASE_H__

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <vector>
#include <sys/time.h>

#include "fec_context.h"
#include "fft_base.h"
#include "gf_base.h"
#include "misc.h"
#include "property.h"
#include "vec_buffers.h"
#include "vec_cast.h"
#include "vec_poly.h"
#include "vec_slice.h"
#include "vec_vector.h"

#ifdef QUADIRON_USE_SIMD

#include "simd.h"

#endif // #ifdef QUADIRON_USE_SIMD

namespace quadiron {

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
    virtual ~FecCode() = default;

    /** Return the number of output parts.
     *
     * @return `n_parities` for SYSTEMATIC and at least `n_data + n_parities`
     * for NON_SYSTEMATIC
     */
    virtual int get_n_outputs() = 0;

    virtual void encode(
        vec::Vector<T>& output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>& words) = 0;
    virtual void
    encode_post_process(vec::Vector<T>&, std::vector<Properties>&, off_t){};
    virtual void encode(
        vec::Buffers<T>&,
        std::vector<Properties>&,
        off_t,
        vec::Buffers<T>&){};
    virtual void
    encode_post_process(vec::Buffers<T>&, std::vector<Properties>&, off_t){};
    virtual void decode_add_data(int /* fragment_index */, int /* row */){};
    virtual void decode_add_parities(int /* fragment_index */, int /* row */){};
    virtual void decode_build(void){};

    /**
     * Decode a vector of words
     *
     * @param context decoding context
     * @param output original data (must be of n_data length)
     * @param props properties bound to parity fragments
     * @param offset offset in the data fragments
     * @param words input words if SYSTEMATIC must be n_data,
     *  if NON_SYSTEMATIC get_n_outputs()
     */
    virtual void decode(
        DecodeContext<T>& context,
        vec::Vector<T>& output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>& words);

    virtual void decode(
        DecodeContext<T>& context,
        vec::Buffers<T>& output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>& words);

    bool readw(T* ptr, std::istream* stream);
    bool writew(T val, std::ostream* stream);

    bool read_pkt(char* pkt, std::istream& stream);
    bool write_pkt(char* pkt, std::ostream& stream, size_t bytes);

    void encode_streams_horizontal(
        std::vector<std::istream*> input_data_bufs,
        std::vector<std::ostream*> output_parities_bufs,
        std::vector<Properties>& output_parities_props);

    void encode_streams_vertical(
        std::vector<std::istream*> input_data_bufs,
        std::vector<std::ostream*> output_parities_bufs,
        std::vector<Properties>& output_parities_props);

    bool decode_streams_horizontal(
        std::vector<std::istream*> input_data_bufs,
        std::vector<std::istream*> input_parities_bufs,
        std::vector<Properties>& input_parities_props,
        std::vector<std::ostream*> output_data_bufs);

    virtual std::unique_ptr<DecodeContext<T>> init_context_dec(
        vec::Vector<T>& fragments_ids,
        std::vector<Properties>& input_props,
        size_t size = 0,
        vec::Buffers<T>* output = nullptr);

    bool decode_streams_vertical(
        std::vector<std::istream*> input_data_bufs,
        std::vector<std::istream*> input_parities_bufs,
        std::vector<Properties>& input_parities_props,
        std::vector<std::ostream*> output_data_bufs);

    void encode_blocks_vertical(
        std::vector<uint8_t*> data_bufs,
        std::vector<uint8_t*> parities_bufs,
        std::vector<Properties>& parities_props,
        std::vector<bool> wanted_idxs,
        size_t block_size_bytes);

    bool decode_blocks_vertical(
        std::vector<uint8_t*> data_bufs,
        std::vector<uint8_t*> parities_bufs,
        std::vector<Properties>& parities_props,
        std::vector<int> missing_idxs,
        std::vector<bool> wanted_idxs,
        size_t block_size_bytes);

    const gf::Field<T>& get_gf()
    {
        return *gf;
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
    std::unique_ptr<gf::Field<T>> gf = nullptr;
    std::unique_ptr<fft::FourierTransform<T>> fft = nullptr;
    std::unique_ptr<fft::FourierTransform<T>> fft_2k = nullptr;
    // This vector MUST be initialized by derived Class using multiplicative FFT
    std::unique_ptr<vec::Vector<T>> inv_r_powers = nullptr;
    // This vector MUST be initialized by derived Class using multiplicative FFT
    std::unique_ptr<vec::Vector<T>> r_powers = nullptr;
    // buffers for intermediate symbols used for systematic FNT
    std::unique_ptr<vec::Buffers<T>> dec_inter_codeword;

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
        DecodeContext<T>& context,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>& words);

    virtual void decode_apply(
        DecodeContext<T>& context,
        vec::Vector<T>& output,
        vec::Vector<T>& words);

    virtual void decode_prepare(
        DecodeContext<T>& context,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>& words);

    virtual void decode_apply(
        DecodeContext<T>& context,
        vec::Buffers<T>& output,
        vec::Buffers<T>& words);
};

/// Create an encoder.
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
inline bool FecCode<T>::readw(T* ptr, std::istream* stream)
{
    if (word_size == 1) {
        uint8_t c;
        if (stream->read(reinterpret_cast<char*>(&c), sizeof(c))) {
            *ptr = c;
            return true;
        }
    } else if (word_size == 2) {
        uint16_t s;
        if (stream->read(reinterpret_cast<char*>(&s), sizeof(s))) {
            *ptr = s;
            return true;
        }
    } else if (word_size == 4) {
        uint32_t s;
        if (stream->read(reinterpret_cast<char*>(&s), sizeof(s))) {
            *ptr = s;
            return true;
        }
    } else if (word_size == 8) {
        uint64_t s;
        if (stream->read(reinterpret_cast<char*>(&s), sizeof(s))) {
            *ptr = s;
            return true;
        }
    } else if (word_size == 16) {
        __uint128_t s;
        if (stream->read(reinterpret_cast<char*>(&s), sizeof(s))) {
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
        if (stream->write(reinterpret_cast<char*>(&c), sizeof(c)))
            return true;
    } else if (word_size == 2) {
        uint16_t s = val;
        if (stream->write(reinterpret_cast<char*>(&s), sizeof(s)))
            return true;
    } else if (word_size == 4) {
        uint32_t s = val;
        if (stream->write(reinterpret_cast<char*>(&s), sizeof(s)))
            return true;
    } else if (word_size == 8) {
        uint64_t s = val;
        if (stream->write(reinterpret_cast<char*>(&s), sizeof(s)))
            return true;
    } else if (word_size == 16) {
        __uint128_t s = val;
        if (stream->write(reinterpret_cast<char*>(&s), sizeof(s)))
            return true;
    } else {
        assert(false && "no such size");
    }
    return false;
}

template <typename T>
inline bool FecCode<T>::read_pkt(char* pkt, std::istream& stream)
{
    return static_cast<bool>(stream.read(pkt, buf_size));
}

template <typename T>
inline bool FecCode<T>::write_pkt(char* pkt, std::ostream& stream, size_t bytes)
{
    return static_cast<bool>(stream.write(pkt, bytes));
}

/** Encode streams
 *
 * @param input_data_bufs vector size must be exactly n_data
 * @param output_parities_bufs vector size must be exactly n_parities
 * (set entries no nullptr when not wanted)
 * @param output_parities_props vector size must be exactly n_parities
 * (set entries no nullptr when not wanted)
 *
 * @note all streams must be of equal size
 */
template <typename T>
void FecCode<T>::encode_streams_horizontal(
    std::vector<std::istream*> input_data_bufs,
    std::vector<std::ostream*> output_parities_bufs,
    std::vector<Properties>& output_parities_props)
{
    bool cont = true;
    off_t offset = 0;

    assert(input_data_bufs.size() == n_data);
    assert(output_parities_bufs.size() == n_outputs);
    assert(output_parities_props.size() == n_outputs);

    vec::Vector<T> words(*(this->gf), n_data);
    vec::Vector<T> output(*(this->gf), get_n_outputs());

    // clear property vectors
    for (auto& props : output_parities_props) {
        props.clear();
    }

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

        timeval t1 = tick();
        uint64_t start = hw_timer();
        encode(output, output_parities_props, offset, words);
        uint64_t end = hw_timer();
        uint64_t t2 = hrtime_usec(t1);

        total_enc_usec += t2;
        total_encode_cycles += (end - start) / word_size;
        n_encode_ops++;

        for (unsigned i = 0; i < n_outputs; i++) {
            T tmp = output.get(i);
            writew(tmp, output_parities_bufs[i]);
        }
        offset++;
    }
}

template <typename T>
void FecCode<T>::encode_streams_vertical(
    std::vector<std::istream*> input_data_bufs,
    std::vector<std::ostream*> output_parities_bufs,
    std::vector<Properties>& output_parities_props)
{
    assert(input_data_bufs.size() == n_data);
    assert(output_parities_bufs.size() == n_outputs);
    assert(output_parities_props.size() == n_outputs);

    // clear property vectors
    for (auto& props : output_parities_props) {
        props.clear();
    }

    bool cont = true;
    off_t offset = 0;

    // vector of buffers storing data read from chunk
    vec::Buffers<char> words_char(n_data, buf_size);
    const std::vector<char*> words_mem_char = words_char.get_mem();
    // vector of buffers storing data that are performed in encoding, i.e. FFT
    vec::Buffers<T> words(n_data, pkt_size);
    const std::vector<T*> words_mem_T = words.get_mem();

    int output_len = get_n_outputs();

    // vector of buffers storing data that are performed in encoding, i.e. FFT
    vec::Buffers<T> output(output_len, pkt_size);
    const std::vector<T*> output_mem_T = output.get_mem();
    // vector of buffers storing data in output chunk
    vec::Buffers<char> output_char(output_len, buf_size);
    const std::vector<char*> output_mem_char = output_char.get_mem();

    reset_stats_enc();

    // Number of bytes would be read from each input stream
    // We suppose that these stream returns the same quantity of data.
    size_t read_bytes = buf_size;

    while (cont) {
        for (unsigned i = 0; i < n_data; i++) {
            if (!read_pkt(words_mem_char.at(i), *(input_data_bufs[i]))) {
                read_bytes = input_data_bufs[i]->gcount();
                // Zero-out trailing part
                std::fill_n(
                    words_mem_char.at(i) + read_bytes,
                    buf_size - read_bytes,
                    0);

                cont = false;
            }
        }

        if (read_bytes == 0) {
            break;
        }

        vec::pack<char, T>(
            words_mem_char, words_mem_T, n_data, pkt_size, word_size);

        timeval t1 = tick();
        uint64_t start = hw_timer();
        encode(output, output_parities_props, offset, words);
        uint64_t end = hw_timer();
        uint64_t t2 = hrtime_usec(t1);

        total_enc_usec += t2;
        total_encode_cycles += (end - start) / buf_size;
        n_encode_ops++;

        vec::unpack<T, char>(
            output_mem_T, output_mem_char, output_len, pkt_size, word_size);

        for (unsigned i = 0; i < n_outputs; i++) {
            write_pkt(
                output_mem_char.at(i), *(output_parities_bufs[i]), read_bytes);
        }
        offset += pkt_size;
    }
}

/** Decode streams
 *
 * @param input_data_bufs vector size must be exactly n_data
 * (set entries to nullptr when missing)
 * @param input_parities_bufs vector size must be exactly n_parities
 * (set entries to nullptr when missing)
 * @param input_parities_props vector size must be exactly n_parities
 * (set entries to nullptr when missing)
 * @param output_data_bufs vector size must be exactly n_data
 * (set entries to nullptr when not wanted)
 *
 * @note All streams must be of equal size
 *
 * @return true if decode succeded, else false
 */
template <typename T>
bool FecCode<T>::decode_streams_horizontal(
    std::vector<std::istream*> input_data_bufs,
    std::vector<std::istream*> input_parities_bufs,
    std::vector<Properties>& input_parities_props,
    std::vector<std::ostream*> output_data_bufs)
{
    off_t offset = 0;
    bool cont = true;

    unsigned fragment_index = 0;
    unsigned parity_index = 0;
    unsigned avail_data_nb = 0;

    if (type == FecType::SYSTEMATIC) {
        assert(input_data_bufs.size() == n_data);
    }
    assert(input_parities_bufs.size() == n_outputs);
    assert(input_parities_props.size() == n_outputs);
    assert(output_data_bufs.size() == n_data);

    reset_stats_dec();
    // ids of received fragments, from 0 to codelen-1
    vec::Vector<T> fragments_ids(*(this->gf), n_data);

    if (type == FecType::SYSTEMATIC) {
        for (unsigned i = 0; i < n_data; i++) {
            if (input_data_bufs[i] != nullptr) {
                decode_add_data(fragment_index, i);
                fragments_ids.set(fragment_index, i);
                fragment_index++;
            }
        }
        avail_data_nb = fragment_index;
        // data is in clear so nothing to do
        if (fragment_index == n_data)
            return true;
    }

    vec::Vector<T> avail_parity_ids(*(this->gf), n_data - avail_data_nb);

    if (fragment_index < n_data) {
        // finish with parities available
        for (unsigned i = 0; i < n_outputs; i++) {
            if (input_parities_bufs[i] != nullptr) {
                decode_add_parities(fragment_index, i);
                unsigned j = (type == FecType::SYSTEMATIC) ? n_data + i : i;
                fragments_ids.set(fragment_index, j);
                avail_parity_ids.set(parity_index, i);
                fragment_index++;
                parity_index++;
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

    vec::Vector<T> words(*(this->gf), n_words);
    vec::Vector<T> output(*(this->gf), n_data);

    std::unique_ptr<DecodeContext<T>> context =
        init_context_dec(fragments_ids, input_parities_props);
    while (true) {
        words.zero_fill();
        if (type == FecType::SYSTEMATIC) {
            for (unsigned i = 0; i < avail_data_nb; ++i) {
                unsigned data_idx = fragments_ids.get(i);
                T tmp;
                if (!readw(&tmp, input_data_bufs[data_idx])) {
                    cont = false;
                    break;
                }
                words.set(i, tmp);
            }
            if (!cont)
                break;
        }
        for (unsigned i = 0; i < n_data - avail_data_nb; ++i) {
            unsigned parity_idx = avail_parity_ids.get(i);
            T tmp;
            if (!readw(&tmp, input_parities_bufs[parity_idx])) {
                cont = false;
                break;
            }
            words.set(avail_data_nb + i, tmp);
        }
        if (!cont)
            break;

        timeval t1 = tick();
        uint64_t start = hw_timer();
        decode(*context, output, input_parities_props, offset, words);
        uint64_t end = hw_timer();
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

        offset++;
    }

    return true;
}

/**
 * Perform a Lagrange interpolation to find the coefficients of the
 * polynomial. The algorithm is mainly extracted from the paper \cite fnt-rs.
 *
 * @note If all fragments are available ifft(words) is enough
 *
 * We have to find the first \f$k\f$ coefficients the numerator of the following
 * expression:
 * \f[
 *  \frac{P(x)}{A(x)} = (\sum_{i=0}^{k-1}\frac{n_i}{x-x_i}) \% x^k
 * \f]
 * where \f$x_i\f$, \f$v_i\f$ are the index and value of the \f$i\f$th
 * received fragments, and
 * \f{eqnarray*}{
 *  &n_i = \frac{v_i}{{A'_i}(x_i)} \\
 *  &A(x) = \prod_{i=0}^{k-1} (x - x_i)
 * \f}
 *
 * Using Taylor series
 * \f[ \frac{1}{x_i - x} = \sum \frac{x^j}{{x_i}^{j+1}} \f]
 * we rewrite the expression into
 * \f{eqnarray*}{
 * \frac{P(x)}{A(x)}
 *  &= -\sum_{i=0}^{k-1}
 *      \sum_{j=0}^{k-1}(
 *          \frac{n_i}{x_i} \times {x_i}^{-j} \times x^j) \\
 *  &= -\sum_{j=0}^{k-1}
 *      x^j \times
 *      \sum_{i=0}^{k-1}(
 *          \frac{n_i}{x_i} \times {x_i}^{-j}) \\
 * \f}
 *
 * As \f$x_i\f$ is a power of the root of unity order \f$n\f$, i.e.
 * \f$x_i = r^{z_i}\f$, we define the following polynomial
 * \f[
 *  N(x) = \sum_{i=0}^{k-1}{\frac{n_i}{x_i} \times x^{z_i}}
 * \f]
 * Hence,
 * \f{eqnarray*}{
 * \frac{P(x)}{A(x)}
 *  &= -\sum_{j=0}^{k-1} N(r^{-j}) x^j
 * \f}
 * Note, these \f$N'(r^{-j})\f$ are the first \f$k\f$ elements of
 * \f$iFFT_n(N(x))\f$.
 *
 * The last step is to perform multiply \f$A(x)\f$ to
 * \f$Q(x) := \sum_{j=0}^{k-1} N(r^{-j}) x^j\f$. It can be done by using the
 * convolution theorem:
 * \f[ P(x) = iFFT_{2k}( FFT_{2k}(A(x)) \cdot FFT_{2k}(Q(x)) ) \f]
 * where \f$FFT_{2k}, iFFT_{2k}\f$ are FFT and inverse FFT of length \f$2k\f$.
 *
 * @param context calculated for given indices of received fragments
 * @param output must be exactly n_data
 * @param props special values dictionary must be exactly n_data
 * @param offset used to locate special values
 * @param words \f$v=(v_0, v_1, ..., v_{k-1})\f$ \f$k\f$ must be exactly n_data
 */
template <typename T>
void FecCode<T>::decode(
    DecodeContext<T>& context,
    vec::Vector<T>& output,
    const std::vector<Properties>& props,
    off_t offset,
    vec::Vector<T>& words)
{
    // prepare for decoding
    decode_prepare(context, props, offset, words);

    // Lagrange interpolation
    decode_apply(context, output, words);
}

/* Initialize context for decoding
 * It supports for FEC using multiplicative FFT over FNT
 */
template <typename T>
std::unique_ptr<DecodeContext<T>> FecCode<T>::init_context_dec(
    vec::Vector<T>& fragments_ids,
    std::vector<Properties>& input_props,
    size_t size,
    vec::Buffers<T>* output)
{
    if (this->inv_r_powers == nullptr) {
        throw LogicError("FEC base: vector (inv_r)^i must be initialized");
    }
    if (this->r_powers == nullptr) {
        throw LogicError("FEC base: vector r^i must be initialized");
    }
    if (this->fft == nullptr) {
        throw LogicError("FEC base: FFT must be initialized");
    }

    int k = this->n_data; // number of fragments received
    // vector x=(x_0, x_1, ..., x_k-1)
    vec::Vector<T> vx(*(this->gf), k);
    for (int i = 0; i < k; ++i) {
        vx.set(i, r_powers->get(fragments_ids.get(i)));
    }

    return std::make_unique<DecodeContext<T>>(
        *gf,
        *fft,
        *fft_2k,
        fragments_ids,
        input_props,
        vx,
        n_data,
        n,
        -1,
        size,
        output);
}

/* Prepare for decoding
 * It supports for FEC using multiplicative FFT over FNT
 */
template <typename T>
void FecCode<T>::decode_prepare(
    DecodeContext<T>& context,
    const std::vector<Properties>& props,
    off_t offset,
    vec::Vector<T>& words)
{
    const vec::Vector<T>& fragments_ids = context.get_fragments_id();
    for (unsigned i = 0; i < this->n_data; ++i) {
        const int j = fragments_ids.get(i);
        if (props[j].is_marked(context.props_indices[j], offset)) {
            // Check if the symbol is a special case whick is marked by
            // `OOR_MARK`, i.e. true. Note: this check is necessary when
            // word_size is not large enough to cover all symbols of the field.
            // Following check is used for FFT over FNT where the single special
            // case symbol equals card - 1
            if (props[j].marker(context.props_indices[j]) == OOR_MARK) {
                words.set(i, this->gf->card() - 1);
            }
            context.props_indices.at(j)++;
        }
    }
}

/**
 * Perform decoding on received chunks by using pre-calculated context.
 *
 * @param context calculated for given indices of received fragments
 * @param output must be exactly n_data
 * @param words vector \f$v=(v_0, v_1, ..., v_{k-1})\f$, \f$k\f$ must be exactly
 * n_data
 */
template <typename T>
void FecCode<T>::decode_apply(
    DecodeContext<T>& context,
    vec::Vector<T>& output,
    vec::Vector<T>& words)
{
    unsigned len_2k = context.get_len_2k();

    const vec::Vector<T>& fragments_ids = context.get_fragments_id();
    vec::Vector<T>& inv_A_i = context.get_vector(CtxVec::INV_A_I);
    vec::Vector<T>& A_fft_2k = context.get_vector(CtxVec::A_FFT_2K);
    vec::Vector<T>& vec1_n = context.get_vector(CtxVec::N1);
    vec::Vector<T>& vec2_n = context.get_vector(CtxVec::N2);
    vec::Vector<T>& vec1_2k = context.get_vector(CtxVec::V2K1);
    vec::Vector<T>& vec2_2k = context.get_vector(CtxVec::V2K2);

    unsigned k = this->n_data; // number of fragments received

    // compute N(x) and stored in vec1_n
    vec1_n.zero_fill();
    for (unsigned i = 0; i <= k - 1; ++i) {
        vec1_n.set(
            fragments_ids.get(i), this->gf->mul(words.get(i), inv_A_i.get(i)));
    }

    // compute vec2_n = FFT(vec1_n)
    this->fft->fft_inv(vec2_n, vec1_n);

    // vec_tmp_2k: first k elements from vec2_n
    //             last (len_2k - k) elements are padded
    vec::Slice<T> vec_tmp_k(&vec2_n, k);
    vec::ZeroExtended<T> vec_tmp_2k(vec_tmp_k, len_2k);

    // compute FFT_2k(Q(x))
    this->fft_2k->fft(vec1_2k, vec_tmp_2k);

    // multiply FFT_2k(A(x)) to FFT_2k(Q(x)), results are stored in vec1_2k
    vec1_2k.hadamard_mul(&A_fft_2k);

    // compute iFFT_{2k}( FFT_{2k}(A(x)) \cdot FFT_{2k}(Q(x)) )
    this->fft_2k->ifft(vec2_2k, vec1_2k);

    // perform negative
    vec2_2k.neg();

    // get decoded symbols are the first k elements of vec2_2k
    for (unsigned i = 0; i < k; ++i)
        output.set(i, vec2_2k.get(i));
}

/********** Decoding over vec::PolyBuf **********/

/** Decode streams
 *
 * @param input_data_bufs vector size must be exactly n_data
 * (set entries to nullptr when missing)
 * @param input_parities_bufs vector size must be exactly n_parities
 * (set entries to nullptr when missing)
 * @param input_parities_props vector size must be exactly n_parities
 * (set entries to nullptr when missing)
 * @param output_data_bufs vector size must be exactly n_data
 * (set entries to nullptr when not wanted)
 *
 * @pre All streams must be of equal size
 *
 * @return true if decode succeeded, else false
 */
template <typename T>
bool FecCode<T>::decode_streams_vertical(
    std::vector<std::istream*> input_data_bufs,
    std::vector<std::istream*> input_parities_bufs,
    std::vector<Properties>& input_parities_props,
    std::vector<std::ostream*> output_data_bufs)
{
    bool cont = true;
    off_t offset = 0;

    unsigned fragment_index = 0;
    unsigned parity_index = 0;
    unsigned avail_data_nb = 0;

    if (type == FecType::SYSTEMATIC) {
        assert(input_data_bufs.size() == n_data);
    }
    assert(input_parities_bufs.size() == n_outputs);
    assert(input_parities_props.size() == n_outputs);
    assert(output_data_bufs.size() == n_data);

    // ids of received fragments, from 0 to codelen-1
    vec::Vector<T> fragments_ids(*(this->gf), n_data);

    if (type == FecType::SYSTEMATIC) {
        for (unsigned i = 0; i < n_data; i++) {
            if (input_data_bufs[i] != nullptr) {
                decode_add_data(fragment_index, i);
                fragments_ids.set(fragment_index, i);
                fragment_index++;
            }
        }
        avail_data_nb = fragment_index;
        // data is in clear so nothing to do
        if (fragment_index == n_data)
            return true;
    }

    vec::Vector<T> avail_parity_ids(*(this->gf), n_data - avail_data_nb);

    if (fragment_index < n_data) {
        // finish with parities available
        for (unsigned i = 0; i < n_outputs; i++) {
            if (input_parities_bufs[i] != nullptr) {
                decode_add_parities(fragment_index, i);
                unsigned j = (type == FecType::SYSTEMATIC) ? n_data + i : i;
                fragments_ids.set(fragment_index, j);
                avail_parity_ids.set(parity_index, i);
                fragment_index++;
                parity_index++;
                // stop when we have enough parities
                if (fragment_index == n_data)
                    break;
            }
        }
        // unable to decode
        if (fragment_index < n_data)
            return false;
    }
    fragments_ids.sort();

    decode_build();

    // vector of buffers storing data read from chunk
    vec::Buffers<char> words_char(n_data, buf_size);
    const std::vector<char*> words_mem_char = words_char.get_mem();
    // vector of buffers storing data that are performed in encoding, i.e. FFT
    vec::Buffers<T> words(n_data, pkt_size);
    const std::vector<T*> words_mem_T = words.get_mem();

    int output_len = n_data;

    // vector of buffers storing data that are performed in decoding, i.e. FFT
    vec::Buffers<T> output(output_len, pkt_size);
    const std::vector<T*> output_mem_T = output.get_mem();
    // vector of buffers storing data in output chunk
    vec::Buffers<char> output_char(output_len, buf_size);
    const std::vector<char*> output_mem_char = output_char.get_mem();

    std::unique_ptr<DecodeContext<T>> context = init_context_dec(
        fragments_ids, input_parities_props, pkt_size, &output);

    reset_stats_dec();

    // Number of bytes would be read from each input stream
    // We suppose that these stream returns the same quantity of data.
    size_t read_bytes = buf_size;

    while (cont) {
        if (type == FecType::SYSTEMATIC) {
            for (unsigned i = 0; i < avail_data_nb; i++) {
                unsigned data_idx = fragments_ids.get(i);
                if (!read_pkt(
                        words_mem_char.at(i), *(input_data_bufs[data_idx]))) {
                    read_bytes = input_data_bufs[data_idx]->gcount();
                    // Zero-out trailing part
                    std::fill_n(
                        words_mem_char.at(i) + read_bytes,
                        buf_size - read_bytes,
                        0);

                    cont = false;
                }
            }
        }
        for (unsigned i = 0; i < n_data - avail_data_nb; ++i) {
            unsigned parity_idx = avail_parity_ids.get(i);
            if (!read_pkt(
                    words_mem_char.at(avail_data_nb + i),
                    *(input_parities_bufs[parity_idx]))) {
                read_bytes = input_parities_bufs[parity_idx]->gcount();
                // Zero-out trailing part
                std::fill_n(
                    words_mem_char.at(avail_data_nb + i) + read_bytes,
                    buf_size - read_bytes,
                    0);

                cont = false;
            }
        }

        if (read_bytes == 0) {
            break;
        }

        vec::pack<char, T>(
            words_mem_char, words_mem_T, n_data, pkt_size, word_size);

        timeval t1 = tick();
        uint64_t start = hw_timer();
        decode(*context, output, input_parities_props, offset, words);
        uint64_t end = hw_timer();
        uint64_t t2 = hrtime_usec(t1);

        total_dec_usec += t2;
        total_decode_cycles += (end - start) / word_size;
        n_decode_ops++;

        vec::unpack<T, char>(
            output_mem_T, output_mem_char, output_len, pkt_size, word_size);

        for (unsigned i = 0; i < n_data; i++) {
            if (output_data_bufs[i] != nullptr) {
                write_pkt(
                    output_mem_char.at(i), *(output_data_bufs[i]), read_bytes);
            }
        }
        offset += pkt_size;
    }

    return true;
}

/** Encode blocks
 *
 * @param data_bufs vector size must be exactly n_data
 * (set entries to nullptr when missing)
 * @param parities_bufs vector size must be exactly n_outputs
 * (set entries to nullptr when not wanted)
 * @param parities_props vector size must be exactly n_outputs
 * (set entries to nullptr when not wanted)
 * @param wanted_idxs bool array of missing_idxs of len n_outputs indicating
 * wanted (value 1) or not wanted fragments (value 0) - wanted blocks MUST BE
 * allocated by caller
 * @param block_size_bytes the block size in bytes
 *
 * @pre All blocks must be of equal size
 */
template <typename T>
void FecCode<T>::encode_blocks_vertical(
    std::vector<uint8_t*> data_bufs,
    std::vector<uint8_t*> parities_bufs,
    std::vector<Properties>& parities_props,
    std::vector<bool> wanted_idxs,
    size_t block_size_bytes)
{
    assert(data_bufs.size() == n_data);
    assert(parities_bufs.size() == n_outputs);
    assert(parities_props.size() == n_outputs);

    // clear property vectors
    for (auto& props : parities_props) {
        props.clear();
    }

    size_t offset = 0;
    size_t block_size = block_size_bytes / word_size;

    // vector of buffers storing data read from chunk
    vec::Buffers<uint8_t> words_char(n_data, buf_size);
    const std::vector<uint8_t*> words_mem_char = words_char.get_mem();
    // vector of buffers storing data that are performed in encoding, i.e. FFT
    vec::Buffers<T> words(n_data, pkt_size);
    const std::vector<T*> words_mem_T = words.get_mem();

    int output_len = get_n_outputs();

    // vector of buffers storing data that are performed in encoding, i.e. FFT
    vec::Buffers<T> output(output_len, pkt_size);
    const std::vector<T*> output_mem_T = output.get_mem();
    // vector of buffers storing data in output chunk
    vec::Buffers<uint8_t> output_char(output_len, buf_size);
    const std::vector<uint8_t*> output_mem_char = output_char.get_mem();

    reset_stats_enc();

    while (offset < block_size) {
        size_t remain_size = block_size - offset;
        size_t copy_size = std::min(pkt_size, remain_size);
        for (unsigned i = 0; i < n_data; i++) {
            memcpy(
                reinterpret_cast<char*>(words_mem_char.at(i)),
                data_bufs[i] + offset * word_size,
                copy_size * word_size);
        }

        // Zero-out trailing part of data
        if (copy_size < pkt_size) {
            const size_t copy_bytes = copy_size * word_size;
            const size_t trailing_bytes = buf_size - copy_bytes;
            for (unsigned i = 0; i < n_data; i++) {
                memset(
                    reinterpret_cast<char*>(words_mem_char.at(i)) + copy_bytes,
                    0,
                    trailing_bytes);
            }
        }

        vec::pack<uint8_t, T>(
            words_mem_char, words_mem_T, n_data, pkt_size, word_size);

        timeval t1 = tick();
        uint64_t start = hw_timer();
        encode(output, parities_props, offset, words);
        uint64_t end = hw_timer();
        uint64_t t2 = hrtime_usec(t1);

        total_enc_usec += t2;
        total_encode_cycles += (end - start) / (copy_size * word_size);
        n_encode_ops++;

        vec::unpack<T, uint8_t>(
            output_mem_T, output_mem_char, output_len, pkt_size, word_size);

        for (unsigned i = 0; i < n_outputs; i++) {
            if (wanted_idxs[i]) {
                memcpy(
                    parities_bufs[i] + offset * word_size,
                    reinterpret_cast<char*>(output_mem_char.at(i)),
                    copy_size * word_size);
            }
        }
        offset += pkt_size;
    }
}

/** Decode blocks
 *
 * @param data_bufs vector size must be exactly n_data
 * (set entries to nullptr when missing)
 * @param parities_bufs vector size must be exactly n_parities
 * (set entries to nullptr when missing)
 * @param parities_props vector size must be exactly n_parities
 * (set entries to nullptr when missing)
 * @param missing_idxs array of missing indexes of vector size code_len
 * indicating presence (value 1) or absence of fragments (value 0). It applies
 * for both data and parities
 * - non missing blocks MUST BE allocated by caller
 * @param wanted_idxs bool array of wanted indexes of vector size n_data
 * indicating wanted (value 1) or not wanted fragments (value 0). It applies
 * only for data
 * - wanted blocks MUST BE allocated
 * by caller
 * @param block_size_bytes the block size in bytes
 *
 * @pre All blocks must be of equal size
 *
 * @return true if decode succeeded, else false
 */
template <typename T>
bool FecCode<T>::decode_blocks_vertical(
    std::vector<uint8_t*> data_bufs,
    std::vector<uint8_t*> parities_bufs,
    std::vector<Properties>& parities_props,
    std::vector<int> missing_idxs,
    std::vector<bool> wanted_idxs,
    size_t block_size_bytes)
{
    size_t offset = 0;
    size_t block_size = block_size_bytes / word_size;

    unsigned fragment_index = 0;
    unsigned parity_index = 0;
    unsigned avail_data_nb = 0;

    if (type == FecType::SYSTEMATIC) {
        assert(data_bufs.size() == n_data);
    }
    assert(parities_bufs.size() == n_outputs);
    assert(parities_props.size() == n_outputs);

    // ids of received fragments, from 0 to codelen-1
    vec::Vector<T> fragments_ids(*(this->gf), n_data);

    if (type == FecType::SYSTEMATIC) {
        for (unsigned i = 0; i < n_data; i++) {
            if (!missing_idxs[i]) {
                decode_add_data(fragment_index, i);
                fragments_ids.set(fragment_index, i);
                fragment_index++;
            }
            avail_data_nb = fragment_index;
            // data is in clear so nothing to do
            if (fragment_index == n_data)
                return true;
        }
    }

    vec::Vector<T> avail_parity_ids(*(this->gf), n_data - avail_data_nb);

    if (fragment_index < n_data) {
        // finish with parities available
        for (unsigned i = 0; i < n_outputs; i++) {
            unsigned j = (type == FecType::SYSTEMATIC) ? n_data + i : i;
            if (!missing_idxs[j]) {
                decode_add_parities(fragment_index, i);
                fragments_ids.set(fragment_index, j);
                avail_parity_ids.set(parity_index, i);
                fragment_index++;
                parity_index++;
                // stop when we have enough parities
                if (fragment_index == n_data)
                    break;
            }
        }
        // unable to decode
        if (fragment_index < n_data)
            return false;
    }
    fragments_ids.sort();

    decode_build();

    // vector of buffers storing data read from chunk
    vec::Buffers<uint8_t> words_char(n_data, buf_size);
    const std::vector<uint8_t*> words_mem_char = words_char.get_mem();
    // vector of buffers storing data that are performed in encoding, i.e. FFT
    vec::Buffers<T> words(n_data, pkt_size);
    const std::vector<T*> words_mem_T = words.get_mem();

    int output_len = n_data;

    // vector of buffers storing data that are performed in decoding, i.e. FFT
    vec::Buffers<T> output(output_len, pkt_size);
    const std::vector<T*> output_mem_T = output.get_mem();
    // vector of buffers storing data in output chunk
    vec::Buffers<uint8_t> output_char(output_len, buf_size);
    const std::vector<uint8_t*> output_mem_char = output_char.get_mem();

    std::unique_ptr<DecodeContext<T>> context =
        init_context_dec(fragments_ids, parities_props, pkt_size, &output);

    reset_stats_dec();

    while (offset < block_size) {
        size_t remain_size = block_size - offset;
        size_t copy_size = std::min(pkt_size, remain_size);
        if (type == FecType::SYSTEMATIC) {
            for (unsigned i = 0; i < avail_data_nb; i++) {
                unsigned data_idx = fragments_ids.get(i);
                memcpy(
                    reinterpret_cast<char*>(words_mem_char.at(i)),
                    data_bufs[data_idx] + offset * word_size,
                    copy_size * word_size);
            }
        }
        for (unsigned i = 0; i < n_data - avail_data_nb; ++i) {
            unsigned parity_idx = avail_parity_ids.get(i);
            memcpy(
                reinterpret_cast<char*>(words_mem_char.at(avail_data_nb + i)),
                parities_bufs[parity_idx] + offset * word_size,
                copy_size * word_size);
        }

        // Zero-out trailing part of data
        if (copy_size < pkt_size) {
            const size_t copy_bytes = copy_size * word_size;
            const size_t trailing_bytes = buf_size - copy_bytes;
            for (unsigned i = 0; i < n_data; i++) {
                memset(
                    reinterpret_cast<char*>(words_mem_char.at(i)) + copy_bytes,
                    0,
                    trailing_bytes);
            }
        }

        vec::pack<uint8_t, T>(
            words_mem_char, words_mem_T, n_data, pkt_size, word_size);

        timeval t1 = tick();
        uint64_t start = hw_timer();
        decode(*context, output, parities_props, offset, words);
        uint64_t end = hw_timer();
        uint64_t t2 = hrtime_usec(t1);

        total_dec_usec += t2;
        total_decode_cycles += (end - start) / word_size;
        n_decode_ops++;

        vec::unpack<T, uint8_t>(
            output_mem_T, output_mem_char, output_len, pkt_size, word_size);

        for (unsigned i = 0; i < n_data; i++) {
            if (wanted_idxs[i]) {
                memcpy(
                    data_bufs[i] + offset * word_size,
                    reinterpret_cast<char*>(output_mem_char.at(i)),
                    copy_size * word_size);
            }
        }
        offset += pkt_size;
    }

    return true;
}

/**
 * Perform a Lagrange interpolation to find the coefficients of the
 * polynomial
 *
 * @note If all fragments are available ifft(words) is enough
 *
 * @param context decoding context
 * @param output must be exactly n_data
 * @param offset used to locate special values
 * @param words vector \f$v=(v_0, v_1, ..., v_{k-1})\f$, \f$k\f$ must be exactly
 * n_data
 */
template <typename T>
void FecCode<T>::decode(
    DecodeContext<T>& context,
    vec::Buffers<T>& output,
    const std::vector<Properties>& props,
    off_t offset,
    vec::Buffers<T>& words)
{
    // prepare for decoding
    decode_prepare(context, props, offset, words);

    // Lagrange interpolation
    decode_apply(context, output, words);

    if (type == FecType::SYSTEMATIC) {
        this->fft->fft(*dec_inter_codeword, output);
        for (unsigned i = 0; i < this->n_data; i++) {
            output.copy(*dec_inter_codeword, i, i);
        }
    }
}

/* Prepare for decoding
 * It supports for FEC using multiplicative FFT over FNT
 */
template <typename T>
void FecCode<T>::decode_prepare(
    DecodeContext<T>& context,
    const std::vector<Properties>& props,
    off_t offset,
    vec::Buffers<T>& words)
{
    // FIXME: could we integrate this preparation into vec::pack?
    // It will reduce a loop on all data
    const vec::Vector<T>& fragments_ids = context.get_fragments_id();
    off_t offset_max = offset + pkt_size;

    T thres = (this->gf->card() - 1);
    for (unsigned i = 0; i < this->n_data; i++) {
        unsigned frag_id = fragments_ids.get(i);
        if (type == FecType::SYSTEMATIC && frag_id < this->n_data) {
            continue;
        }
        T* chunk = words.get(i);
        if (type == FecType::SYSTEMATIC) {
            frag_id -= this->n_data;
        }

        // loop over marked symbols
        while (props[frag_id].in_range(
            context.props_indices.at(frag_id), offset, offset_max)) {
            const size_t loc_offset =
                props[frag_id].location(context.props_indices.at(frag_id));
            // As loc.offset := offset + j
            const size_t j = (loc_offset - offset);

            // Check if the symbol is a special case whick is marked by
            // `OOR_MARK`.
            // Note: this check is necessary when word_size is not large
            // enough to cover all symbols of the field. Following check is
            // used for FFT over FNT where the single special case symbol
            // equals card - 1
            if (props[frag_id].marker(context.props_indices.at(frag_id))
                == OOR_MARK) {
                chunk[j] = thres;
            }
            context.props_indices.at(frag_id)++;
        }
    }
}

/**
 * Perform a Lagrange interpolation to find the coefficients of the
 * polynomial
 *
 * @note If all fragments are available ifft(words) is enough
 *
 * @param context decoding context
 * @param output must be exactly n_data
 * @param words vector \f$v=(v_0, v_1, ..., v_{k-1})\f$, \f$k\f$ must be exactly
 * n_data
 */
template <typename T>
void FecCode<T>::decode_apply(
    DecodeContext<T>& context,
    vec::Buffers<T>& output,
    vec::Buffers<T>& words)
{
    vec::Vector<T>& inv_A_i = context.get_vector(CtxVec::INV_A_I);
    vec::Vector<T>& A_fft_2k = context.get_vector(CtxVec::A_FFT_2K);

    vec::Buffers<T>& buf1_k = context.get_buffer(CtxBuf::K1);
    vec::Buffers<T>& buf1_n = context.get_buffer(CtxBuf::N1);
    vec::Buffers<T>& buf2_n = context.get_buffer(CtxBuf::N2);
    vec::Buffers<T>& buf1_2k = context.get_buffer(CtxBuf::B2K1);
    vec::Buffers<T>& buf2_2k = context.get_buffer(CtxBuf::B2K2);

    // compute N'(x) = sum_i{n_i * x^z_i}
    // where n_i=v_i/A'_i(x_i)
    this->gf->mul_vec_to_vecp(inv_A_i, words, buf1_k);

    // compute buf2_n
    this->fft->fft_inv(buf2_n, buf1_n);

    this->fft_2k->fft(buf1_2k, output);

    // multiply FFT(A) and buf2_2k
    this->gf->mul_vec_to_vecp(A_fft_2k, buf1_2k, buf1_2k);

    this->fft_2k->ifft(buf2_2k, buf1_2k);

    // negatize output
    this->gf->neg(output);
}

} // namespace fec
} // namespace quadiron

#endif
