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
#ifndef __NTTEC_FEC_RS_CCAT_H__
#define __NTTEC_FEC_RS_CCAT_H__

// use <wirehair.h> instead when linking against installed
#include "ccat.h"
#include "CCatCpp.h"

#include "fec_base.h"
#include "vec_vector.h"

/// C++ convenience wrapper around the C API
/// Override OnRecoveredData() to receive data
class CCatWrapper
{
public:
    // Initialize and pick window size in data length
    bool Initialize(unsigned WindowPackets = 16)
    {
        Destroy();

        CCatSettings settings;
        settings.AppContextPtr = this;
        settings.OnRecoveredData = [](CCatOriginal original, void* context)
        {
            CCatWrapper* thiz = (CCatWrapper*)context;
            thiz->OnRecoveredData(original);
        };
        settings.WindowMsec = 100;
        settings.WindowPackets = WindowPackets;

        CCatResult result = ccat_create(&settings, &Codec);
        if (result != CCat_Success)
        {
            Error = true;
            return false;
        }

        Error = false;
        return true;
    }

    bool IsError() const
    {
        return Error;
    }

    void Destroy()
    {
        if (Codec)
        {
            ccat_destroy(Codec);
            Codec = nullptr;
        }
    }

    virtual ~CCatWrapper()
    {
        Destroy();
    }

    // Received data handlers:

    void OnOriginal(const CCatOriginal& original)
    {
        CCatResult result = ccat_decode_original(Codec, &original);
        if (result != CCat_Success)
            Error = true;
    }

    void OnRecovery(const CCatRecovery& recovery)
    {
        CCatResult result = ccat_decode_recovery(Codec, &recovery);
        if (result != CCat_Success)
            Error = true;
    }

    // Outgoing data:

    void SendOriginal(const CCatOriginal& original)
    {
        CCatResult result = ccat_encode_original(Codec, &original);
        if (result != CCat_Success)
            Error = true;
    }

    bool SendRecovery(CCatRecovery& recovery)
    {
        CCatResult result = ccat_encode_recovery(Codec, &recovery);
        if (result != CCat_Success)
            Error = true;
        return result == CCat_Success;
    }

protected:
    virtual void OnRecoveredData(const CCatOriginal& original)
    {
        // Default does nothing
        (void)original;
    }

    bool Error = false;
    CCatCodec Codec = nullptr;
};

class CCatWrapperDecoder : public CCatWrapper
{
public:
    uint64_t RecoveredPackets = 0;
    uint64_t OriginalPackets = 0;

    void OnRecoveredData(const CCatOriginal& original) override
    {
        std::cout << "CCatWrapperDecoder OnRecoveredData\n";
        std::cout << "\tSequenceNumber " << original.SequenceNumber << "\n";
        std::cout << "\tBytes " << original.Bytes << "\n";
        ++RecoveredPackets;
        ++OriginalPackets;
    }
};

namespace nttec {
namespace fec {

typedef unsigned char u8;
/** Reed-Solomon (RS) Erasure code over GF(2<sup>n</sup>) (Cauchy or
 *  Vandermonde) using the WireHair library
 */
template <typename T>
class RsCCat : public FecCode<T> {
  public:
    RsMatrixType mat_type;

    RsCCat(
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

    ~RsCCat() = default;

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
        // std::cout << "init encoder..";
        encoder = std::unique_ptr<CCatWrapper>(new CCatWrapper());
        if (!encoder->Initialize(this->n_data))
        {
            std::cout << "Failed to initialize encoder\n";
        }
        std::cout << "done\n";

        // std::cout << "init decoder..";
        decoder = std::unique_ptr<CCatWrapperDecoder>(new CCatWrapperDecoder());
        if (!decoder->Initialize(this->n_data))
        {
            std::cout << "Failed to initialize decoder\n";
        }
        // std::cout << "done\n";
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
        std::vector<T*>* vec_output = output.get_mem();
        std::vector<T*>* vec_data = words.get_mem();

        size_t block_id = this->n_data * offset;
        // std::cout << "encode from "  << block_id << " IsError " << encoder->IsError() << std::endl;

        for (unsigned idx = 0; idx < this->n_data; ++idx, ++block_id) {
            CCatOriginal block;
            block.Data = reinterpret_cast<uint8_t*>(vec_data->at(idx));
            block.Bytes = this->buf_size;
            block.SequenceNumber = block_id;
            encoder->SendOriginal(block);
            // std::cout << "send original "  << block_id << " result " << encoder->IsError() << std::endl;
        }

        for (unsigned idx = 0; idx < this->n_parities; ++idx) {
            CCatRecovery block;
            encoder->SendRecovery(block);
            // std::cout << "Parity " << idx << "\n";
            // std::cout << "\tCount " << block.Count << "\n";
            // std::cout << "\tSequenceStart " << block.SequenceStart << "\n";
            // std::cout << "\tBytes " << block.Bytes << "\n";
            // std::cout << "\tRecoveryRow " << block.RecoveryRow << "\n";
            // std::cout << "\tResult " << encoder->IsError() << std::endl;
            std::memcpy(vec_output->at(idx), block.Data, this->buf_size);
        }

        // std::cout << "WORDS"; words.dump();
        // std::cout << "ENCODED"; output.dump();
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
        // std::vector<T*>* vec_output = output.get_mem();
        std::vector<T*>* vec_data = words.get_mem();

        // std::cout << "Received: "; words.dump();

        const unsigned receivedNb = fragments_ids->get_n();

        for (unsigned i = 0; i < receivedNb; ++i) {
            const unsigned block_id = fragments_ids->get(i);
            // std::cout << "frag "  << block_id << std::endl;
            if (block_id < this->n_data) {
                CCatOriginal block;
                block.Data = reinterpret_cast<uint8_t*>(vec_data->at(i));
                block.Bytes = this->buf_size;
                block.SequenceNumber = block_id;

                decoder->OnOriginal(block);
            } else {
                CCatRecovery block;
                block.Data = reinterpret_cast<uint8_t*>(vec_data->at(i));
                block.Bytes = this->buf_size;
                block.SequenceStart = 0;
                // block.SequenceNumber = block_id;

                decoder->OnRecovery(block);
            }
            // std::cout << "Result decoder " << decoder->IsError() << std::endl;
        }
    }

  private:
    vec::Vector<T>* fragments_ids;
    std::unique_ptr<CCatWrapper> encoder = nullptr;
    std::unique_ptr<CCatWrapperDecoder> decoder = nullptr;
};

} // namespace fec
} // namespace nttec

#endif
