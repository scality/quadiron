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
#ifndef __QUAD_BENCH_COMP_PERF_FFT_H__
#define __QUAD_BENCH_COMP_PERF_FFT_H__

#include <chrono>
#include <functional>
#include <vector>

#include "fft_2n.h"
#include "gf_prime.h"
#include "perf_base.h"
#include "simd/simd.h"
#include "vec_buffers.h"

namespace gf = quadiron::gf;
namespace fft = quadiron::fft;
namespace simd = quadiron::simd;

// Code length
constexpr size_t MIN_FFT_LEN = 32;
constexpr size_t MAX_FFT_LEN = 1024;
// Maximal ratio of output length on input length
// out_len / in_len = {1, .., MAX_RATIO_OUT_IN}
constexpr unsigned MAX_RATIO_OUT_IN = 4;

template <typename T>
class PerfRadix2Fft : public PerfBase<T> {
  public:
    std::unique_ptr<gf::Field<T>> gf = nullptr;

    PerfRadix2Fft()
    {
        if (sizeof(T) == 2) {
            this->q = 257;
            this->word_size = 1;
        } else if (sizeof(T) == 4) {
            this->q = static_cast<T>(65537);
            this->word_size = 2;
        } else {
            throw "Wrong TypeParam for PerfRadix2Fft tests";
        }

        gf = gf::alloc<gf::Field<T>, gf::Prime<T>>(this->q);
    }

    void fft(benchmark::State& st, bool ifft = false)
    {
        const size_t pkt_size = simd::countof<T>() * st.range(0);
        const size_t ratio = st.range(1);
        const size_t fft_len = st.range(2);

        if (pkt_size <= 0) {
            st.SkipWithError("First argument for packet size must be positive");
        }
        if (ratio <= 0) {
            st.SkipWithError(
                "Second argument for ratio output/input length must be "
                "positive");
        }
        if (fft_len <= 0) {
            st.SkipWithError("Thrid argument for FFT length must be positive");
        }
        if (fft_len < ratio) {
            st.SkipWithError(
                "FFT length is too small vs ratio output/input length");
        }

        const size_t input_len = fft_len / ratio;
        fft::Radix2<T> fft_2n(*gf, fft_len, input_len, pkt_size);

        vec::Buffers<T> input(input_len, pkt_size);
        vec::Buffers<T> output(fft_len, pkt_size);
        this->randomize_data(input);

        fft::OpCounter counter;
        std::chrono::duration<double, std::micro> elapsed_us{};

        if (ifft) {
            counter = fft_2n.ifft_op_counter(input_len);
            auto start = std::chrono::high_resolution_clock::now();
            for (auto _ : st) {
                fft_2n.ifft(output, input);
            }
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_us += std::chrono::duration_cast<
                std::chrono::duration<double, std::micro>>(end - start);
        } else {
            counter = fft_2n.fft_op_counter(input_len);
            auto start = std::chrono::high_resolution_clock::now();
            for (auto _ : st) {
                fft_2n.fft(output, input);
            }
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_us += std::chrono::duration_cast<
                std::chrono::duration<double, std::micro>>(end - start);
        }

        const size_t bytes = fft_len * pkt_size * this->word_size;
        st.counters.insert(
            {{"(1) FFT length", fft_len},
             {"(2) Input length", input_len},
             {"(3) Packet size",
              benchmark::Counter(
                  pkt_size,
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)},
             {"(4) Butterfly", counter.butterfly},
             {"(4-) Add", counter.add},
             {"(4-) Sub", counter.sub},
             {"(4-) Neg", counter.neg},
             {"(4-) Mul", counter.mul},
             {"(5) Elapsed time (us)", elapsed_us.count() / st.iterations()},
             {"(6) Throughput",
              benchmark::Counter(
                  st.iterations() * bytes,
                  benchmark::Counter::kIsRate,
                  benchmark::Counter::OneK::kIs1024)}});
    }
};

BENCHMARK_TEMPLATE_DEFINE_F(PerfRadix2Fft, Fft16, uint16_t)
(benchmark::State& st)
{
    this->fft(st);
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfRadix2Fft, Fft16)
    ->RangeMultiplier(2)
    ->Ranges({{MIN_VEC_LEN, MAX_VEC_LEN} /* Number of vectors per buffer */,
              {1, MAX_RATIO_OUT_IN} /* Ratio out over in length */,
              {MIN_FFT_LEN, 256} /* FFT length */});

BENCHMARK_TEMPLATE_DEFINE_F(PerfRadix2Fft, iFft16, uint16_t)
(benchmark::State& st)
{
    this->fft(st, true /* inverse FFT */);
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfRadix2Fft, iFft16)
    ->RangeMultiplier(2)
    ->Ranges({{MIN_VEC_LEN, MAX_VEC_LEN} /* Number of vectors per buffer */,
              {1, MAX_RATIO_OUT_IN} /* Ratio out over in length */,
              {MIN_FFT_LEN, 256} /* FFT length */});

BENCHMARK_TEMPLATE_DEFINE_F(PerfRadix2Fft, Fft32, uint32_t)
(benchmark::State& st)
{
    this->fft(st);
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfRadix2Fft, Fft32)
    ->RangeMultiplier(2)
    ->Ranges({{MIN_VEC_LEN, MAX_VEC_LEN} /* Number of vectors per buffer */,
              {1, MAX_RATIO_OUT_IN} /* Ratio out over in length */,
              {MIN_FFT_LEN, MAX_FFT_LEN} /* FFT length */});

BENCHMARK_TEMPLATE_DEFINE_F(PerfRadix2Fft, iFft32, uint32_t)
(benchmark::State& st)
{
    this->fft(st, true /* inverse FFT */);
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfRadix2Fft, iFft32)
    ->RangeMultiplier(2)
    ->Ranges({{MIN_VEC_LEN, MAX_VEC_LEN} /* Number of vectors per buffer */,
              {1, MAX_RATIO_OUT_IN} /* Ratio out over in length */,
              {MIN_FFT_LEN, MAX_FFT_LEN} /* FFT length */});

#endif
