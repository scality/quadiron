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
#include "misc.h"
#include "vec_buffers.h"

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
class SimdTestFnt : public ::testing::Test {
  public:
    SimdTestFnt()
    {
        if (sizeof(T) == 2) {
            this->q = 257;
        } else if (sizeof(T) == 4) {
            this->q = static_cast<T>(65537);
        } else {
            throw "Wrong TypeParam for SimdTestFnt tests";
        }

        this->distribution =
            std::make_unique<std::uniform_int_distribution<uint32_t>>(0, q - 1);
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

    bool is_equal(simd::VecType x, simd::VecType y)
    {
        return simd::is_zero(simd::bit_xor(x, y));
    }

    simd::VecType from_msb8_mask(simd::MetaType meta)
    {
        const size_t n = simd::countof<uint8_t>();
        uint8_t buf[n];
        simd::VecType* vec = reinterpret_cast<simd::VecType*>(buf);

        for (unsigned i = 0; i < n; ++i) {
            buf[i] = meta & 1;
            meta >>= 1;
        }

        return vec[0];
    }

    simd::VecType mod_mul(simd::VecType x, simd::VecType y)
    {
        const size_t n = simd::countof<T>();
        T _x[n];
        T _y[n];
        T _z[n];
        simd::store_to_mem(reinterpret_cast<simd::VecType*>(_x), x);
        simd::store_to_mem(reinterpret_cast<simd::VecType*>(_y), y);
        for (unsigned i = 0; i < n; i++) {
            _z[i] = (quadiron::DoubleSizeVal<T>(_x[i]) * _y[i]) % q;
        }

        simd::VecType* vec = reinterpret_cast<simd::VecType*>(_z);

        return vec[0];
    }

    /* Butterfly Cooley-Tukey operation
     * x <- x + c * y
     * y <- x - c * y
     */
    void butterfly_ct(simd::VecType c, simd::VecType& x, simd::VecType& y)
    {
        const size_t n = simd::countof<T>();
        T c_buf[n];
        T x_buf[n];
        T y_buf[n];

        simd::store_to_mem(reinterpret_cast<simd::VecType*>(c_buf), c);
        simd::store_to_mem(reinterpret_cast<simd::VecType*>(x_buf), x);
        simd::store_to_mem(reinterpret_cast<simd::VecType*>(y_buf), y);

        for (unsigned i = 0; i < n; ++i) {
            T mul = (quadiron::DoubleSizeVal<T>(c_buf[i]) * y_buf[i]) % q;
            T u = (x_buf[i] + mul) % q;
            T v = x_buf[i] >= mul ? x_buf[i] - mul : q + x_buf[i] - mul;

            x_buf[i] = u;
            y_buf[i] = v;
        }

        x = simd::load_to_reg(reinterpret_cast<simd::VecType*>(x_buf));
        y = simd::load_to_reg(reinterpret_cast<simd::VecType*>(y_buf));
    }

    /* Butterfly Genteleman-Sande operation
     * x <- x + y
     * y <- c * (x - y)
     */
    void butterfly_gs(simd::VecType c, simd::VecType& x, simd::VecType& y)
    {
        const size_t n = simd::countof<T>();
        T c_buf[n];
        T x_buf[n];
        T y_buf[n];

        simd::store_to_mem(reinterpret_cast<simd::VecType*>(c_buf), c);
        simd::store_to_mem(reinterpret_cast<simd::VecType*>(x_buf), x);
        simd::store_to_mem(reinterpret_cast<simd::VecType*>(y_buf), y);

        for (unsigned i = 0; i < n; ++i) {
            T sub = x_buf[i] >= y_buf[i] ? x_buf[i] - y_buf[i]
                                         : q + x_buf[i] - y_buf[i];
            T u = (x_buf[i] + y_buf[i]) % q;
            T v = (quadiron::DoubleSizeVal<T>(c_buf[i]) * sub) % q;
            x_buf[i] = u;
            y_buf[i] = v;
        }

        x = simd::load_to_reg(reinterpret_cast<simd::VecType*>(x_buf));
        y = simd::load_to_reg(reinterpret_cast<simd::VecType*>(y_buf));
    }

    /* Butterfly Genteleman-Sande simple operation where y = 0
     * y <- c * x
     */
    void butterfly_simple_gs(simd::VecType c, simd::VecType& x)
    {
        const size_t n = simd::countof<T>();
        T c_buf[n];
        T x_buf[n];

        simd::store_to_mem(reinterpret_cast<simd::VecType*>(c_buf), c);
        simd::store_to_mem(reinterpret_cast<simd::VecType*>(x_buf), x);

        for (unsigned i = 0; i < n; ++i) {
            x_buf[i] = (quadiron::DoubleSizeVal<T>(c_buf[i]) * x_buf[i]) % q;
        }

        x = simd::load_to_reg(reinterpret_cast<simd::VecType*>(x_buf));
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
        std::cout << "Average nb of CPU cycles of " << text << ": "
                  << avg_cycles_nb << "\n";
    }

    template <typename TFunc>
    void core_op_perf(const std::string& text, const TFunc& f)
    {
        std::cout << text << "\n";
        std::cout << "\tVectors nb\t\tAverage nb of CPU cycles\n";
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
            ;

            std::cout << "\t" << vec_len << "\t\t" << avg_cycles_nb << "\n";
        }
        std::cout << "\n";
    }

    template <typename TFunc>
    void butterfly_perf(const std::string& text, const TFunc& f)
    {
        std::cout << text << "\n";
        std::cout << "\tVectors nb\t\tAverage nb of CPU cycles\n";
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
            ;

            std::cout << "\t" << vec_len << "\t\t" << avg_cycles_nb << "\n";
        }
        std::cout << "\n";
    }

    T q;
    std::unique_ptr<std::uniform_int_distribution<uint32_t>> distribution;
    std::vector<size_t> arr_vec_len =
        {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384};
    size_t iters_nb = 1e3;
};

using AllTypes = ::testing::Types<uint16_t, uint32_t>;
TYPED_TEST_CASE(SimdTestFnt, AllTypes);

TYPED_TEST(SimdTestFnt, TestMetaToMask) // NOLINT
{
    simd::MetaType meta = static_cast<simd::MetaType>(0x89abcdefUL);
    simd::VecType got = simd::from_msb8_mask(meta);
    simd::VecType expected = this->from_msb8_mask(meta);

    ASSERT_TRUE(this->is_equal(expected, got));
}

TYPED_TEST(SimdTestFnt, TestModAddSub) // NOLINT
{
    for (unsigned i = 0; i < 100; ++i) {
        simd::VecType x = this->rand_vec(this->q / 2);
        simd::VecType y = this->rand_vec(this->q / 2);

        simd::VecType u = simd::mod_add<TypeParam>(x, y);
        simd::VecType v = simd::mod_sub<TypeParam>(u, x);
        simd::VecType z = simd::mod_add<TypeParam>(v, x);

        ASSERT_TRUE(this->is_equal(y, v));
        ASSERT_TRUE(this->is_equal(u, z));
    }
}

TYPED_TEST(SimdTestFnt, TestModNeg) // NOLINT
{
    for (unsigned i = 0; i < 100; ++i) {
        simd::VecType x = this->rand_vec();

        simd::VecType y = simd::mod_neg<TypeParam>(x);
        simd::VecType u = simd::mod_sub<TypeParam>(simd::zero(), x);

        ASSERT_TRUE(this->is_equal(u, y));
    }
}

TYPED_TEST(SimdTestFnt, TestModMul) // NOLINT
{
    for (unsigned i = 0; i < 100; ++i) {
        simd::VecType x = this->rand_vec(0, this->q - 1);
        simd::VecType y = this->rand_vec();

        // check mod_mul
        simd::VecType u = simd::mod_mul<TypeParam>(x, y);
        simd::VecType v = this->mod_mul(x, y);
        ASSERT_TRUE(this->is_equal(v, u));

        // check mod_mul_safe
        simd::VecType a = simd::card_minus_one<TypeParam>();
        simd::VecType b = simd::card_minus_one<TypeParam>();
        simd::VecType c = this->mod_mul(a, b);
        simd::VecType d = simd::mod_mul<TypeParam>(a, b);
        simd::VecType e = simd::mod_mul_safe<TypeParam>(a, b);

        ASSERT_FALSE(this->is_equal(c, d));
        ASSERT_TRUE(this->is_equal(c, e));
    }
}

TYPED_TEST(SimdTestFnt, TestButterflyCt) // NOLINT
{
    for (unsigned i = 0; i < 100; ++i) {
        std::vector<TypeParam> r_values = {
            1,
            static_cast<TypeParam>(
                this->distribution->operator()(quadiron::prng())),
            static_cast<TypeParam>(this->q - 1)};

        for (const TypeParam r : r_values) {
            const simd::CtGsCase ct_case =
                simd::get_case<TypeParam>(r, this->q);

            simd::VecType c = simd::set_one(r);

            simd::VecType x = this->rand_vec();
            simd::VecType y = this->rand_vec();
            simd::VecType x_expected = this->copy(x);
            simd::VecType y_expected = this->copy(y);

            this->butterfly_ct(c, x_expected, y_expected);
            simd::butterfly_ct<TypeParam>(ct_case, c, x, y);

            ASSERT_TRUE(this->is_equal(x_expected, x));
            ASSERT_TRUE(this->is_equal(y_expected, y));
        }
    }
}

TYPED_TEST(SimdTestFnt, TestButterflyGs) // NOLINT
{
    for (unsigned i = 0; i < 100; ++i) {
        std::vector<TypeParam> r_values = {
            1,
            static_cast<TypeParam>(
                this->distribution->operator()(quadiron::prng())),
            static_cast<TypeParam>(this->q - 1)};

        for (const TypeParam r : r_values) {
            const simd::CtGsCase ct_case =
                simd::get_case<TypeParam>(r, this->q);

            simd::VecType c = simd::set_one(r);

            simd::VecType x = this->rand_vec();
            simd::VecType y = this->rand_vec();
            simd::VecType x_expected = this->copy(x);
            simd::VecType y_expected = this->copy(y);

            this->butterfly_gs(c, x_expected, y_expected);
            simd::butterfly_gs<TypeParam>(ct_case, c, x, y);

            ASSERT_TRUE(this->is_equal(x_expected, x));
            ASSERT_TRUE(this->is_equal(y_expected, y));

            this->butterfly_simple_gs(c, x_expected);
            simd::butterfly_simple_gs<TypeParam>(ct_case, c, x);

            ASSERT_TRUE(this->is_equal(x_expected, x));
        }
    }
}

TYPED_TEST(SimdTestFnt, PerfSimdSingle) // NOLINT
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

TYPED_TEST(SimdTestFnt, PerfSimdBuf) // NOLINT
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

TYPED_TEST(SimdTestFnt, PerfModBuf) // NOLINT
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

TYPED_TEST(SimdTestFnt, PerfPackUnpack) // NOLINT
{
    std::cout << "Pack & Unpack"
              << "\n";
    std::cout << "\tVectors nb\t\tAverage nb of CPU cycles\n";
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

TYPED_TEST(SimdTestFnt, PerfButterfly) // NOLINT
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

#endif
