/*
 * Copyright 2017-2019 Scality
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
#ifndef __QUAD_BENCH_COMP_PERF_BASE_H__
#define __QUAD_BENCH_COMP_PERF_BASE_H__

#include <benchmark/benchmark.h>

#include "vec_buffers.h"

namespace vec = quadiron::vec;

template <typename T>
class PerfBase : public benchmark::Fixture {
  public:
    T q;
    unsigned word_size;

    void randomize_data(const vec::Buffers<T>& vec, int max = 0)
    {
        const T upper = (max == 0) ? q : max + 1;
        std::uniform_int_distribution<T> dis(0, upper - 1);

        const size_t vec_n = vec.get_n();
        const size_t vec_size = vec.get_size();

        const std::vector<T*>& mem = vec.get_mem();
        for (size_t i = 0; i < vec_n; ++i) {
            for (size_t j = 0; j < vec_size; j++) {
                mem[i][j] = dis(quadiron::prng());
            }
        }
    }

    void randomize_byte_buffer(uint8_t* buf, size_t size)
    {
        std::uniform_int_distribution<> dis(0, 256);
        for (size_t i = 0; i < size; ++i) {
            buf[i] = dis(quadiron::prng());
        }
    }
};

#endif
