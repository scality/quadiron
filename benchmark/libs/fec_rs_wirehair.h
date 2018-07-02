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
#ifndef __NTTEC_FEC_RS_WIREHAIR_H__
#define __NTTEC_FEC_RS_WIREHAIR_H__

// use <wirehair.h> instead when linking against installed
#include "wirehair.h"

#include "fec_base.h"
#include "vec_vector.h"

namespace nttec {
namespace fec {

typedef unsigned char u8;
/** Reed-Solomon (RS) Erasure code over GF(2<sup>n</sup>) (Cauchy or
 *  Vandermonde) using the WireHair library
 */
template <typename T>
class RsWH : public FecCode<T> {
  public:
    RsMatrixType mat_type;

    RsWH(
        unsigned word_size,
        unsigned n_data,
        unsigned n_parities,
        RsMatrixType type,
        size_t pkt_size = 8 * 1024)
        : FecCode<T>(
              FecType::SYSTEMATIC,
              word_size,
              n_data,
              n_parities,
              pkt_size)
    {
        mat_type = type;
        this->fec_init();
    }

    ~RsWH() = default;

    inline void check_params() override
    {
        if (this->n_data > 64000) {
            assert(false); // not support yet
            exit(1);
        }
    }

    inline void init_gf() override {}

    inline void init_fft() override {}

    inline void init_others() override
    {
        const WirehairResult initResult = wirehair_init();

        if (initResult != Wirehair_Success) {
            std::cout << "!!! Wirehair initialization failed: " << initResult
                      << std::endl;
            return;
        }

        message_size = this->n_data * this->buf_size;

        size = this->n_data * this->pkt_size;
        message = std::unique_ptr<std::vector<T>>(new std::vector<T>(size));
    }

    int get_n_outputs() override
    {
        return this->n_parities;
    }

    void encode(
        vec::Vector<T>& output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>& words) override
    {
    }

    void encode(
        vec::Buffers<T>& output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>& words) override
    {
        size_t block_size = this->buf_size;
        std::vector<T*>* vec_output = output.get_mem();
        std::vector<T*>* vec_data = words.get_mem();

        WirehairCodec encoder = wirehair_encoder_create(
            nullptr, vec_data->at(0), message_size, block_size);

        if (!encoder) {
            std::cout << "!!! Failed to create encoder" << std::endl;
        }

        uint32_t writeLen = 0;

        for (unsigned idx = 0; idx < this->n_parities; ++idx) {
            const unsigned blockId = idx + this->n_data;

            WirehairResult encodeResult = wirehair_encode(
                encoder, blockId, vec_output->at(idx), block_size, &writeLen);

            if (encodeResult != Wirehair_Success) {
                std::cout << "wirehair_encode failed" << std::endl;
            }
        }

        wirehair_free(encoder);
    }

    void decode_add_data(int fragment_index, int row) override {}

    void decode_add_parities(int fragment_index, int row) override {}

    void decode_build() override {}

    std::unique_ptr<DecodeContext<T>>
    init_context_dec(vec::Vector<T>& fragments_ids, size_t size) override
    {
        this->fragments_ids = &fragments_ids;
        std::unique_ptr<DecodeContext<T>> context;
        return context;
    }

    void decode(
        const DecodeContext<T>& context,
        vec::Vector<T>& output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>& words) override
    {
    }

    void decode(
        const DecodeContext<T>& context,
        vec::Buffers<T>& output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Buffers<T>& words) override
    {
        size_t block_size = this->buf_size;
        std::vector<T*>* vec_data = words.get_mem();
        std::vector<T*>* vec_output = output.get_mem();

        WirehairCodec decoder =
            wirehair_decoder_create(nullptr, message_size, block_size);

        if (!decoder) {
            std::cout << "!!! Failed to create decoder" << std::endl;
            return;
        }

        unsigned receivedFragsNb = fragments_ids->get_n();
        bool decoding_success = false;

        for (unsigned idx = 0; idx < receivedFragsNb; ++idx) {
            unsigned blockId = fragments_ids->get(idx);

            WirehairResult decodeResult = wirehair_decode(
                decoder, blockId, vec_data->at(idx), block_size);

            if (decodeResult == Wirehair_Success) {
                decoding_success = true;
                break;
            }
        }

        if (decoding_success) {
            WirehairResult recoverResult =
                wirehair_recover(decoder, vec_output->at(0), message_size);
                // wirehair_recover(decoder, &message->at(0), message_size);

            if (recoverResult != Wirehair_Success) {
                std::cout << "wirehair_recover failed" << std::endl;
            }
        } else {
            std::cout << "wirehair_decode failed" << std::endl;
        }

        wirehair_free(decoder);
    }

  private:
    size_t message_size;
    size_t size;
    vec::Vector<T>* fragments_ids;
    std::unique_ptr<std::vector<T>> message = nullptr;
};

} // namespace fec
} // namespace nttec

#endif
