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
#ifndef __QUAD_BENCH_COMP_PERF_FEC_FNT_H__
#define __QUAD_BENCH_COMP_PERF_FEC_FNT_H__

#include "fec_rs_fnt.h"
#include "perf_fec_base.h"

template <typename T>
class PerfFecFnt : public PerfFec<T> {
  public:
    PerfFecFnt()
    {
        if (sizeof(T) == 2) {
            this->q = 257;
            this->word_size = 1;
        } else if (sizeof(T) == 4) {
            this->q = static_cast<T>(65537);
            this->word_size = 2;
        } else {
            throw "Wrong TypeParam for PerfFecFnt tests";
        }

        this->prepare_codec_streams();
    }

    fec::FecType get_fec_type(unsigned type)
    {
        return (type == FEC_NON_SYSTEMATIC) ? fec::FecType::NON_SYSTEMATIC
                                            : fec::FecType::SYSTEMATIC;
    }

    inline bool
    skip_with_error(benchmark::State& st, size_t n_len, size_t n_data)
    {
        if (n_len <= n_data) {
            st.SkipWithError("k should be less than n");
            return true;
        }

        if (2 * n_data >= this->q) {
            st.SkipWithError("RS Fnt has not supported k >= q/2 yet");
            return true;
        }

        return false;
    }

    void enc_dec(benchmark::State& st)
    {
        const size_t pkt_size = simd::countof<T>() * st.range(0);
        const size_t n_len = st.range(1);
        const unsigned inv_rate = 1 + arith::log2<T>(st.range(2));
        const fec::FecType type = get_fec_type(st.range(3));

        // k = inv_rate * n / MAX_INV_RATE
        const size_t n_data = n_len * inv_rate / MAX_INV_RATE;
        const size_t n_parities = n_len - n_data;

        bool err = skip_with_error(st, n_len, n_data);
        if (err) {
            return;
        }

        fec::RsFnt<T> fec(type, this->word_size, n_data, n_parities, pkt_size);

        this->base_enc_dec(st, fec, n_len, n_data, pkt_size);
    }

    void enc_dec_blocks(benchmark::State& st)
    {
        const size_t pkt_size = simd::countof<T>() * st.range(0);
        const size_t block_size = st.range(1);
        const size_t n_len = st.range(2);
        const unsigned inv_rate = 1 + arith::log2<T>(st.range(3));
        const fec::FecType type = get_fec_type(st.range(4));

        // k = inv_rate * n / MAX_INV_RATE
        const size_t n_data = n_len * inv_rate / MAX_INV_RATE;
        const size_t n_parities = n_len - n_data;

        bool err = skip_with_error(st, n_len, n_data);
        if (err) {
            return;
        }

        fec::RsFnt<T> fec(type, this->word_size, n_data, n_parities, pkt_size);

        this->base_enc_dec_blocks(st, fec, n_len, n_data, block_size, pkt_size);
    }

    void enc_dec_streams(benchmark::State& st)
    {
        const size_t pkt_size = simd::countof<T>() * st.range(0);
        const size_t n_len = st.range(1);
        const unsigned inv_rate = 1 + arith::log2<T>(st.range(2));
        const fec::FecType type = get_fec_type(st.range(3));

        // k = inv_rate * n / MAX_INV_RATE
        const size_t n_data = n_len * inv_rate / MAX_INV_RATE;
        const size_t n_parities = n_len - n_data;

        bool err = skip_with_error(st, n_len, n_data);
        if (err) {
            return;
        }

        fec::RsFnt<T> fec(type, this->word_size, n_data, n_parities, pkt_size);

        this->base_enc_dec_streams(st, fec, n_len, n_data, pkt_size);
    }
};

BENCHMARK_TEMPLATE_DEFINE_F(PerfFecFnt, EncDec16, uint16_t)
(benchmark::State& st)
{
    this->enc_dec(st);
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfFecFnt, EncDec16)
    ->RangeMultiplier(2)
    ->Ranges({{MIN_VEC_LEN, MAX_VEC_LEN} /* Pkt_size in #registers */,
              {MIN_CODE_LEN, 256} /* Code length */,
              {1, 1 << (MAX_INV_RATE - 2)} /* Inversed rate */,
              {FEC_NON_SYSTEMATIC, FEC_SYSTEMATIC}});

BENCHMARK_TEMPLATE_DEFINE_F(PerfFecFnt, EncDecBuf16, uint16_t)
(benchmark::State& st)
{
    this->enc_dec_blocks(st);
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfFecFnt, EncDecBuf16)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(2)
    ->Ranges({{MIN_VEC_LEN, MAX_VEC_LEN} /* Pkt_size in #registers */,
              {MIN_BLOCK_SIZE, MAX_BLOCK_SIZE} /* Block size */,
              {MIN_CODE_LEN, 256} /* Code length */,
              {1, 1 << (MAX_INV_RATE - 2)} /* Inversed rate */,
              {FEC_NON_SYSTEMATIC, FEC_SYSTEMATIC}});

BENCHMARK_TEMPLATE_DEFINE_F(PerfFecFnt, EncDecStream16, uint16_t)
(benchmark::State& st)
{
    this->enc_dec_streams(st);
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfFecFnt, EncDecStream16)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(2)
    ->Ranges({{MIN_VEC_LEN, MAX_VEC_LEN} /* Pkt_size in #registers */,
              {MIN_CODE_LEN, 256} /* Code length */,
              {1, 1 << (MAX_INV_RATE - 2)} /* Inversed rate */,
              {FEC_NON_SYSTEMATIC, FEC_SYSTEMATIC}});

BENCHMARK_TEMPLATE_DEFINE_F(PerfFecFnt, EncDec32, uint32_t)
(benchmark::State& st)
{
    this->enc_dec(st);
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfFecFnt, EncDec32)
    ->RangeMultiplier(2)
    ->Ranges({{MIN_VEC_LEN, MAX_VEC_LEN} /* Pkt_size in #registers */,
              {MIN_CODE_LEN, MAX_CODE_LEN} /* Code length */,
              {1, 1 << (MAX_INV_RATE - 2)} /* Inversed rate */,
              {FEC_NON_SYSTEMATIC, FEC_SYSTEMATIC}});

BENCHMARK_TEMPLATE_DEFINE_F(PerfFecFnt, EncDecBuf32, uint32_t)
(benchmark::State& st)
{
    this->enc_dec_blocks(st);
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfFecFnt, EncDecBuf32)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(2)
    ->Ranges({{MIN_VEC_LEN, MAX_VEC_LEN} /* Pkt_size in #registers */,
              {MIN_BLOCK_SIZE, MAX_BLOCK_SIZE} /* Block size */,
              {MIN_CODE_LEN, MAX_CODE_LEN} /* Code length */,
              {1, 1 << (MAX_INV_RATE - 2)} /* Inversed rate */,
              {FEC_NON_SYSTEMATIC, FEC_SYSTEMATIC}});

BENCHMARK_TEMPLATE_DEFINE_F(PerfFecFnt, EncDecStream32, uint32_t)
(benchmark::State& st)
{
    this->enc_dec_streams(st);
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfFecFnt, EncDecStream32)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(2)
    ->Ranges({{MIN_VEC_LEN, MAX_VEC_LEN} /* Pkt_size in #registers */,
              {MIN_CODE_LEN, MAX_CODE_LEN} /* Code length */,
              {1, 1 << (MAX_INV_RATE - 2)} /* Inversed rate */,
              {FEC_NON_SYSTEMATIC, FEC_SYSTEMATIC}});

#endif
