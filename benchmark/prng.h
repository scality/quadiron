/* -*- mode: c++ -*- */
#ifndef __NTTEC_BENCH_PRNG_H__
#define __NTTEC_BENCH_PRNG_H__

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
        rand();
        state += initstate;
        rand();
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
        uint8_t* buffer = (uint8_t*)chunk;
        if (size < 4) {
            buffer[0] = (uint8_t)_rand();
            for (size_t i = 1; i < size; i++) {
                buffer[i] = buffer[0];
            }
        } else {
            uint32_t crc = 0;
            crc = ~crc;
            for (size_t i = 4; i < size; i++) {
                uint8_t val = (uint8_t)_rand();
                buffer[i] = val;
                crc ^= val;
                for (int k = 0; k < 8; k++)
                    crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
            }
            *(uint32_t*)buffer = ~crc;
        }
    }

    uint32_t get_crc(void* chunk, size_t size)
    {
        if (size < 4)
            return 0;
        return *(uint32_t*)chunk;
    }

    bool check_chunk(void* chunk, size_t size)
    {
        uint8_t* buffer = (uint8_t*)chunk;
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
            uint32_t _crc = *(uint32_t*)buffer;
            if (_crc != ~crc)
                return false;
        }
        return true;
    }
};

#endif
