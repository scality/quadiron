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
#include <vector>

#include <functional>
#include <gtest/gtest.h>

#include "arith.h"
#include "core.h"
#include "fec_rs_fnt.h"
#include "fft_2n.h"
#include "gf_prime.h"
#include "misc.h"
#include "vec_buffers.h"

namespace vec = quadiron::vec;
namespace gf = quadiron::gf;
namespace fft = quadiron::fft;
namespace fec = quadiron::fec;

#ifdef QUADIRON_USE_SIMD

#include "simd.h"
#include "simd/simd.h"
#include "simd_fnt.h"

namespace simd = quadiron::simd;

template <typename T>
void show(simd::VecType val)
{
    const size_t n = simd::countof<T>();
    T buffer[n];
    simd::store_to_mem(reinterpret_cast<simd::VecType*>(buffer), val);
    for (unsigned i = 0; i < n; i++) {
        std::cout << unsigned(buffer[i]) << " ";
    }
    std::cout << "\n";
}

template <typename T>
void dump(T* buf, size_t bytes)
{
    const size_t nb = bytes / sizeof(T);
    for (size_t i = 0; i < nb; ++i) {
        std::cout << unsigned(buf[i]) << " ";
    }
    std::cout << "\n";
}

template <typename T>
class PerfFnt : public ::testing::Test {
  public:
    PerfFnt()
    {
        if (sizeof(T) == 2) {
            this->q = 257;
            this->word_size = 1;
        } else if (sizeof(T) == 4) {
            this->q = static_cast<T>(65537);
            this->word_size = 2;
        } else {
            throw "Wrong TypeParam for PerfFnt tests";
        }

        this->distribution =
            std::make_unique<std::uniform_int_distribution<uint32_t>>(0, q - 1);
    }

    void
    buf_rand_data(vec::Buffers<T>& vec, bool has_meta = false, int _max = 0)
    {
        const T max = (_max == 0) ? std::numeric_limits<T>::max() : _max;
        std::uniform_int_distribution<T> dis(0, max - 1);

        const std::vector<T*>& mem = vec.get_mem();
        const size_t size = vec.get_size();
        const size_t n = vec.get_n();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < size; j++) {
                mem[i][j] = dis(quadiron::prng());
            }
        }

        if (has_meta) {
            const std::vector<uint8_t*>& meta = vec.get_meta();
            const size_t meta_size = vec.get_meta_size();
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < meta_size; ++j) {
                    meta[i][j] = static_cast<uint8_t>(dis(quadiron::prng()));
                }
            }
        }
    }

    simd::VecType rand_vec(T lower = 0, T upper_bound = 0)
    {
        const size_t n = simd::countof<T>();
        T buf[n];
        simd::VecType* vec = reinterpret_cast<simd::VecType*>(buf);

        T bound = upper_bound ? upper_bound : q;
        bound -= lower;

        for (unsigned i = 0; i < n; i++) {
            buf[i] = lower + distribution->operator()(quadiron::prng()) % bound;
        }

        return vec[0];
    }

    template <typename Tx>
    void gen_rand_data(Tx* vec, size_t len)
    {
        for (size_t i = 0; i < len; i++) {
            vec[i] = distribution->operator()(quadiron::prng());
        }
    }

    simd::VecType copy(simd::VecType x)
    {
        const size_t n = simd::countof<T>();
        T buf[n];
        T val[n];
        simd::VecType* vec = reinterpret_cast<simd::VecType*>(buf);

        simd::store_to_mem(reinterpret_cast<simd::VecType*>(val), x);
        std::copy_n(val, n, buf);

        return vec[0];
    }

    template <typename TFunc>
    void core_op_perf_single(const std::string& text, const TFunc& f)
    {
        const size_t len = simd::countof<T>();
        // Init a Buffers to obtain aligned memory
        quadiron::vec::Buffers<T> data_buf(2, len);
        T* x = data_buf.get(0);
        T* y = data_buf.get(1);

        for (unsigned i = 0; i < len; ++i) {
            x[i] =
                1
                + (distribution->operator()(quadiron::prng()) % (this->q - 1));
            y[i] =
                1
                + (distribution->operator()(quadiron::prng()) % (this->q - 1));
        }

        simd::VecType* vec_x = reinterpret_cast<simd::VecType*>(x);
        simd::VecType* vec_y = reinterpret_cast<simd::VecType*>(y);

        uint64_t start = quadiron::hw_timer();
        for (unsigned i = 0; i < iters_nb; ++i) {
            simd::VecType _x = simd::load_to_reg(vec_x);
            simd::VecType _y = simd::load_to_reg(vec_y);

            f(_x, _y);

            simd::store_to_mem(vec_x, _x);
        }
        uint64_t end = quadiron::hw_timer();
        double avg_cycles_nb =
            static_cast<double>(end - start) / static_cast<double>(iters_nb);
        std::cout << "#CPU cycles of " << text << ": " << avg_cycles_nb << "\n";
    }

    template <typename TFunc>
    void core_op_perf(const std::string& text, const TFunc& f)
    {
        std::cout << text << "\n";
        std::cout << "\tVectors nb\t#CPU cycles\n";
        for (auto vec_len : arr_vec_len) {
            const size_t len = vec_len * simd::countof<T>();

            // Init a Buffers to obtain aligned memory
            quadiron::vec::Buffers<T> data_buf(2, len);
            T* buf_x = data_buf.get(0);
            T* buf_y = data_buf.get(1);
            gen_rand_data(buf_x, len);
            gen_rand_data(buf_y, len);

            simd::VecType* data_x = reinterpret_cast<simd::VecType*>(buf_x);
            simd::VecType* data_y = reinterpret_cast<simd::VecType*>(buf_y);

            uint64_t start = quadiron::hw_timer();
            for (unsigned i = 0; i < iters_nb; ++i) {
                for (size_t j = 0; j < vec_len; ++j) {
                    simd::VecType x = simd::load_to_reg(&data_x[j]);
                    simd::VecType y = simd::load_to_reg(&data_y[j]);

                    f(x, y);

                    simd::store_to_mem(&data_x[j], x);
                }
            }
            uint64_t end = quadiron::hw_timer();
            double avg_cycles_nb = static_cast<double>(end - start)
                                   / static_cast<double>(iters_nb)
                                   / static_cast<double>(vec_len);

            std::cout << "\t" << vec_len << "\t\t" << avg_cycles_nb << "\n";
        }
        std::cout << "\n";
    }

    template <typename TFunc>
    void butterfly_perf(const std::string& text, const TFunc& f)
    {
        std::cout << text << "\n";
        std::cout << "\tVectors nb\t#CPU cycles\n";
        for (auto vec_len : arr_vec_len) {
            const size_t len = vec_len * simd::countof<T>();

            // Init a Buffers to obtain aligned memory
            quadiron::vec::Buffers<T> data_buf(2, len);
            T* buf_x = data_buf.get(0);
            T* buf_y = data_buf.get(1);
            gen_rand_data(buf_x, len);
            gen_rand_data(buf_y, len);

            simd::VecType* data_x = reinterpret_cast<simd::VecType*>(buf_x);
            simd::VecType* data_y = reinterpret_cast<simd::VecType*>(buf_y);

            T coef = 1
                     + this->distribution->operator()(quadiron::prng())
                           % (this->q - 2);
            const simd::CtGsCase ct_case = simd::get_case<T>(coef, this->q);
            const simd::VecType c = simd::set_one(coef);

            uint64_t start = quadiron::hw_timer();
            for (unsigned i = 0; i < iters_nb; ++i) {
                for (size_t j = 0; j < vec_len; ++j) {
                    simd::VecType x = simd::load_to_reg(&data_x[j]);
                    simd::VecType y = simd::load_to_reg(&data_y[j]);

                    f(ct_case, c, x, y);

                    simd::store_to_mem(&data_x[j], x);
                    simd::store_to_mem(&data_y[j], y);
                }
            }
            uint64_t end = quadiron::hw_timer();

            double avg_cycles_nb = static_cast<double>(end - start)
                                   / static_cast<double>(iters_nb)
                                   / static_cast<double>(vec_len);

            std::cout << "\t" << vec_len << "\t\t" << avg_cycles_nb << "\n";
        }
        std::cout << "\n";
    }

    template <typename TFunc>
    void fft_perf(const std::string& text, size_t fft_len, const TFunc& f)
    {
        std::cout << text << " of length " << fft_len << "\n";
        std::cout << "\tVectors nb\t#CPU cycles\n";

        for (auto vec_len : arr_vec_len) {
            const size_t len = vec_len * simd::countof<T>();

            vec::Buffers<T> input(fft_len, len, true);
            vec::Buffers<T> output(fft_len, len, true);

            buf_rand_data(input);
            buf_rand_data(output);

            uint64_t start = quadiron::hw_timer();
            for (unsigned i = 0; i < iters_nb; ++i) {
                (i % 2) ? f(output, input) : f(input, output);
            }
            uint64_t end = quadiron::hw_timer();

            double avg_cycles_nb = static_cast<double>(end - start)
                                   / static_cast<double>(iters_nb)
                                   / static_cast<double>(vec_len);

            std::cout << "\t" << vec_len << "\t\t" << avg_cycles_nb << "\n";
        }
    }

    void fnt_perf(
        const std::string& text,
        fec::FecCode<T>& fec,
        size_t fft_len,
        int n_data,
        size_t vec_len)
    {
        const size_t len = vec_len * simd::countof<T>();
        const bool has_meta = true;

        vec::Buffers<T> data_frags(n_data, len, has_meta);
        vec::Buffers<T> encoded_frags(fft_len, len, has_meta);

        // It's necessary to set `data_frags` all zeros for `RsNf4` as
        // `data_frags` has not meta
        data_frags.zero_fill();

        std::vector<quadiron::Properties> props(fec.get_n_outputs());

        buf_rand_data(data_frags);

        uint64_t start = quadiron::hw_timer();
        for (unsigned i = 0; i < iters_nb; ++i) {
            for (int i = 0; i < fec.get_n_outputs(); i++) {
                props[i] = quadiron::Properties();
            }
            data_frags.reset_meta();

            fec.encode(encoded_frags, props, 0, data_frags);
        }
        uint64_t end = quadiron::hw_timer();

        double avg_cycles_nb = static_cast<double>(end - start)
                               / static_cast<double>(iters_nb)
                               / static_cast<double>(vec_len);

        std::cout << "\t" << text << "\t\t" << n_data << "\t\t" << fft_len - n_data << "\t\t"
                  << len << "\t\t" << avg_cycles_nb << "\n";
    }

    T q;
    unsigned word_size;
    std::unique_ptr<std::uniform_int_distribution<uint32_t>> distribution;
    std::vector<size_t> arr_vec_len =
        {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384};
    size_t iters_nb = 1e1;
    std::vector<size_t> arr_fft_len = {16, 32, 64, 128, 256};
    std::vector<size_t> arr_k = {8, 16, 32, 64, 128};
};

using AllTypes = ::testing::Types<uint32_t>;
TYPED_TEST_CASE(PerfFnt, AllTypes);

TYPED_TEST(PerfFnt, PerfSimdSingle) // NOLINT
{
    this->core_op_perf_single(
        "Add", [](simd::VecType& x, const simd::VecType& y) {
            x = simd::add<TypeParam>(x, y);
        });

    this->core_op_perf_single(
        "Sub", [](simd::VecType& x, const simd::VecType& y) {
            x = simd::sub<TypeParam>(x, y);
        });

    this->core_op_perf_single(
        "Mul", [](simd::VecType& x, const simd::VecType& y) {
            x = simd::mul<TypeParam>(x, y);
        });

    this->core_op_perf_single(
        "Min", [](simd::VecType& x, const simd::VecType& y) {
            x = simd::min<TypeParam>(x, y);
        });
}

TYPED_TEST(PerfFnt, PerfSimdBuf) // NOLINT
{
    this->core_op_perf("Add", [](simd::VecType& x, const simd::VecType& y) {
        x = simd::add<TypeParam>(x, y);
    });

    this->core_op_perf("Sub", [](simd::VecType& x, const simd::VecType& y) {
        x = simd::sub<TypeParam>(x, y);
    });

    this->core_op_perf("Mul", [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mul<TypeParam>(x, y);
    });

    this->core_op_perf("Min", [](simd::VecType& x, const simd::VecType& y) {
        x = simd::min<TypeParam>(x, y);
    });
}

TYPED_TEST(PerfFnt, PerfModBuf) // NOLINT
{
    this->core_op_perf("ModAdd", [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mod_add<TypeParam>(x, y);
    });

    this->core_op_perf("ModSub", [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mod_sub<TypeParam>(x, y);
    });

    this->core_op_perf("ModMul", [](simd::VecType& x, const simd::VecType& y) {
        x = simd::mod_mul<TypeParam>(x, y);
    });
}

TYPED_TEST(PerfFnt, PerfPackUnpack) // NOLINT
{
    std::cout << "Pack & Unpack"
              << "\n";
    std::cout << "\tVectors nb\t#CPU cycles\n";
    for (const auto vec_len : this->arr_vec_len) {
        const size_t len = vec_len * simd::countof<TypeParam>();

        // Init a Buffers to obtain aligned memory
        quadiron::vec::Buffers<TypeParam> data_buf(1, len);
        TypeParam* buf_data = data_buf.get(0);
        this->gen_rand_data(buf_data, len);

        quadiron::vec::Buffers<simd::MetaType> meta_buf(1, vec_len);
        simd::MetaType* buf_meta = meta_buf.get(0);
        this->gen_rand_data(buf_meta, vec_len);

        simd::VecType* data = reinterpret_cast<simd::VecType*>(buf_data);
        simd::MetaType* meta = buf_meta;

        uint64_t start = quadiron::hw_timer();
        for (unsigned i = 0; i < this->iters_nb; ++i) {
            for (size_t j = 0; j < vec_len; ++j) {
                simd::VecType lo, hi;

                simd::VecType x = simd::load_to_reg(&data[j]);

                simd::unpack<TypeParam>(meta[j], x, hi, lo);
                simd::pack<TypeParam>(lo, hi, x, meta[j]);

                simd::store_to_mem(&data[j], x);
            }
        }
        uint64_t end = quadiron::hw_timer();
        double avg_cycles_nb = static_cast<double>(end - start)
                               / static_cast<double>(this->iters_nb)
                               / static_cast<double>(vec_len);

        std::cout << "\t" << vec_len << "\t\t" << avg_cycles_nb << "\n";
    }
    std::cout << "\n";
}

TYPED_TEST(PerfFnt, PerfButterfly) // NOLINT
{
    this->butterfly_perf(
        "Butterfly_CT",
        [](simd::CtGsCase ct_case,
           const simd::VecType& c,
           simd::VecType& x,
           simd::VecType& y) {
            simd::butterfly_ct<TypeParam>(ct_case, c, x, y);
        });

    this->butterfly_perf(
        "Butterfly_GS",
        [](simd::CtGsCase ct_case,
           const simd::VecType& c,
           simd::VecType& x,
           simd::VecType& y) {
            simd::butterfly_gs<TypeParam>(ct_case, c, x, y);
        });
}

TYPED_TEST(PerfFnt, PerfFftRadix2) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));

    for (auto fft_len : this->arr_fft_len) {
        fft::Radix2<TypeParam> fft_2n(gf, fft_len);

        this->fft_perf(
            "FFT",
            fft_len,
            [&fft_2n](
                vec::Buffers<TypeParam>& output,
                vec::Buffers<TypeParam>& input) { fft_2n.fft(output, input); });
    }
}

TYPED_TEST(PerfFnt, PerfFntEnc) // NOLINT
{
    std::cout << "FNT performance\n";
    std::cout << "\tType\t\tk\t\tm\t\tpkt_size\t\t#CPU cycles\n";

    for (auto fft_len : this->arr_fft_len) {
        for (auto n_data : this->arr_k) {
            if (n_data >= fft_len) {
                continue;
            }
            const int n_parities = fft_len - n_data;

            for (auto vec_len : this->arr_vec_len) {
                fec::RsFnt<TypeParam> fec(
                    fec::FecType::NON_SYSTEMATIC,
                    this->word_size,
                    n_data,
                    n_parities,
                    vec_len);

                this->fnt_perf("Enc", fec, fft_len, n_data, vec_len);
            }
        }
    }
}
#endif
