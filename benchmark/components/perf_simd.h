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

#ifndef __QUAD_BENCH_COMP_PERF_SIMD_H__
#define __QUAD_BENCH_COMP_PERF_SIMD_H__

#include <chrono>
#include <functional>
#include <vector>

#include "perf_base.h"
#include "vec_buffers.h"

#ifdef QUADIRON_USE_SIMD

#include "simd.h"
#include "simd/simd.h"
#include "simd_fnt.h"

namespace simd = quadiron::simd;

// Buffer size range: number of vectors (registers)
constexpr size_t MIN_BUF_LEN = 1;
constexpr size_t MAX_BUF_LEN = 64;
// Iterations nb
constexpr size_t ITERATIONS_NB = 10'000;

template <typename T>
class PerfSimd : public PerfBase<T> {
  public:
    size_t n;
    size_t size;
    std::unique_ptr<vec::Buffers<T>> workload = nullptr;

    PerfSimd()
    {
        if (sizeof(T) == 2) {
            this->q = 257;
            this->word_size = 1;
        } else if (sizeof(T) == 4) {
            this->q = static_cast<T>(65537);
            this->word_size = 2;
        } else {
            throw "Wrong TypeParam for PerfSimd tests";
        }

        n = 2;
        size = simd::countof<T>() * MAX_BUF_LEN;
        workload = std::make_unique<vec::Buffers<T>>(n, size);
        this->randomize_data(*workload);
    }

    template <typename TFunc>
    void field_arith(benchmark::State& st, const TFunc& f)
    {
        const size_t vec_len = st.range(0);

        const std::vector<T*>& mem = workload->get_mem();
        simd::VecType* vec_x = reinterpret_cast<simd::VecType*>(mem[0]);
        simd::VecType* vec_y = reinterpret_cast<simd::VecType*>(mem[1]);

        std::chrono::duration<double, std::nano> elapsed_ns{};
        for (auto _ : st) {
            auto start = std::chrono::high_resolution_clock::now();
            for (unsigned i = 0; i < vec_len; ++i) {
                simd::VecType x = simd::load_to_reg(&vec_x[i]);
                simd::VecType y = simd::load_to_reg(&vec_y[i]);

                f(x, y);

                simd::store_to_mem(&vec_x[i], x);
            }
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_ns += std::chrono::duration_cast<
                std::chrono::duration<double, std::nano>>(end - start);
        }

        const double elapsed = elapsed_ns.count() / st.iterations();
        st.counters.insert({{"(1) Vector len", vec_len},
                            {"(2) Packet size",
                             benchmark::Counter(
                                 simd::countof<T>() * vec_len,
                                 benchmark::Counter::kDefaults,
                                 benchmark::Counter::OneK::kIs1024)},
                            {"(3) Elapsed time (ns)", elapsed},
                            {"(4) Average elapsed time (ns)",
                             elapsed / static_cast<double>(vec_len)}});
    }

    template <typename TFunc>
    void butterfly_op(benchmark::State& st, const TFunc& f)
    {
        const size_t vec_len = st.range(0);
        const std::vector<T*>& mem = workload->get_mem();

        simd::VecType* vec_x = reinterpret_cast<simd::VecType*>(mem[0]);
        simd::VecType* vec_y = reinterpret_cast<simd::VecType*>(mem[1]);

        // to randomize coefficient
        std::uniform_int_distribution<T> dis(1, this->q - 2);

        std::chrono::duration<double, std::nano> elapsed_ns{};
        for (auto _ : st) {
            T coef = dis(quadiron::prng());

            const simd::CtGsCase ct_case = simd::get_case<T>(coef, this->q);
            const simd::VecType c = simd::set_one(coef);

            auto start = std::chrono::high_resolution_clock::now();
            for (unsigned i = 0; i < vec_len; ++i) {
                simd::VecType x = simd::load_to_reg(&vec_x[i]);
                simd::VecType y = simd::load_to_reg(&vec_y[i]);

                f(ct_case, c, x, y);

                simd::store_to_mem(&vec_x[i], x);
                simd::store_to_mem(&vec_y[i], y);
            }
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_ns += std::chrono::duration_cast<
                std::chrono::duration<double, std::nano>>(end - start);
        }

        const size_t pkt_size = vec_len * simd::countof<T>();
        const size_t bytes = pkt_size * this->word_size;
        const double elapsed = elapsed_ns.count() / st.iterations();
        st.counters.insert(
            {{"(1) Vector len", vec_len},
             {"(2) Packet size",
              benchmark::Counter(
                  pkt_size,
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)},
             {"(3) Elapsed time (ns)", elapsed},
             {"(4) Average elapsed time (ns)",
              elapsed / static_cast<double>(vec_len)},
             {"(5) Throughput (B/s)",
              benchmark::Counter(
                  1'000'000'000 * static_cast<double>(bytes) / elapsed,
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)}});
    }
};

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, Add16, uint16_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::add<uint16_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, Add16)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, Sub16, uint16_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::sub<uint16_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, Sub16)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, Mul16, uint16_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mul<uint16_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, Mul16)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, ModAdd16, uint16_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mod_add<uint16_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, ModAdd16)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, ModSub16, uint16_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mod_sub<uint16_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, ModSub16)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, ModMul16, uint16_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mod_mul<uint16_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, ModMul16)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, ButterflyCt16, uint16_t)
(benchmark::State& st)
{
    this->butterfly_op(
        st,
        [](simd::CtGsCase ct_case,
           const simd::VecType& c,
           simd::VecType& x,
           simd::VecType& y) {
            simd::butterfly_ct<uint16_t>(ct_case, c, x, y);
        });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, ButterflyCt16)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, Add32, uint32_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::add<uint32_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, Add32)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, Sub32, uint32_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::sub<uint32_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, Sub32)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, Mul32, uint32_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mul<uint32_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, Mul32)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, ModAdd32, uint32_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mod_add<uint32_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, ModAdd32)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, ModSub32, uint32_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mod_sub<uint32_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, ModSub32)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, ModMul32, uint32_t)(benchmark::State& st)
{
    this->field_arith(st, [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mod_mul<uint32_t>(x, y);
    });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, ModMul32)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

BENCHMARK_TEMPLATE_DEFINE_F(PerfSimd, ButterflyCt32, uint32_t)
(benchmark::State& st)
{
    this->butterfly_op(
        st,
        [](simd::CtGsCase ct_case,
           const simd::VecType& c,
           simd::VecType& x,
           simd::VecType& y) {
            simd::butterfly_ct<uint32_t>(ct_case, c, x, y);
        });
}
// NOLINTNEXTLINE(cert-err58-cpp)
BENCHMARK_REGISTER_F(PerfSimd, ButterflyCt32)
    ->RangeMultiplier(2)
    ->Range(MIN_BUF_LEN, MAX_BUF_LEN);

#endif

#endif
