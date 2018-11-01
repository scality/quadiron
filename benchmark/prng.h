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
#ifndef __QUAD_BENCH_PRNG_H__
#define __QUAD_BENCH_PRNG_H__

#include <cstdint>
#include <cstdlib>

// CRC-32C (iSCSI) polynomial in reversed bit order
#define POLY 0x82f63b78

/*
 *  From http://www.pcg-random.org/
 *
 * *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
 * Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
 */
class PRNG {
  private:
    uint64_t state;
    uint64_t inc;

  public:
    PRNG(uint64_t initstate = 0, uint64_t initseq = 0)
    {
        this->state = initstate;
        this->inc = initseq;
    }

    void srandom(uint64_t initstate, uint64_t initseq)
    {
        state = 0U;
        inc = (initseq << 1u) | 1u;
        _rand();
        state += initstate;
        _rand();
    }

    uint32_t _rand()
    {
        uint64_t oldstate = state;
        state = oldstate * 6364136223846793005ULL + inc;
        uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
        uint32_t rot = oldstate >> 59u;
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

    void gen_chunk(void* chunk, size_t size)
    {
        uint8_t* buffer = static_cast<uint8_t*>(chunk);
        if (size < 4) {
            buffer[0] = static_cast<uint8_t>(_rand()); // narrow_cast
            for (size_t i = 1; i < size; i++) {
                buffer[i] = buffer[0];
            }
        } else {
            uint32_t crc = 0;
            crc = ~crc;
            for (size_t i = 4; i < size; i++) {
                uint8_t val = static_cast<uint8_t>(_rand()); // narrow_cast
                buffer[i] = val;
                crc ^= val;
                for (int k = 0; k < 8; k++)
                    crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
            }
            *reinterpret_cast<uint32_t*>(buffer) = ~crc;
        }
    }

    uint32_t get_crc(void* chunk, size_t size)
    {
        if (size < 4)
            return 0;
        return *static_cast<uint32_t*>(chunk);
    }

    bool check_chunk(void* chunk, size_t size)
    {
        uint8_t* buffer = static_cast<uint8_t*>(chunk);
        if (size < 4) {
            uint8_t val = buffer[0];
            for (size_t i = 1; i < size; i++) {
                if (buffer[i] != val)
                    return false;
            }
        } else {
            uint32_t crc = 0;
            crc = ~crc;
            for (size_t i = 4; i < size; i++) {
                crc ^= buffer[i];
                for (int k = 0; k < 8; k++)
                    crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
            }
            uint32_t _crc = *reinterpret_cast<uint32_t*>(buffer);
            if (_crc != ~crc)
                return false;
        }
        return true;
    }
};

#endif
