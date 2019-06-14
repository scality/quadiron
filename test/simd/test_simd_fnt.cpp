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
#include <vector>
#include <functional>

#include <gtest/gtest.h>

#include "arith.h"
#include "core.h"
#include "misc.h"

#ifdef QUADIRON_USE_SIMD

#include "simd.h"
#include "simd/simd.h"
#include "simd_fnt.h"

namespace simd = quadiron::simd;

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

        T max = upper_bound ? upper_bound : q - 1;
        std::uniform_int_distribution<uint32_t> rand(lower, max);

        for (unsigned i = 0; i < n; i++) {
            buf[i] = rand(quadiron::prng());
        }

        return vec[0];
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

    T q;
    std::unique_ptr<std::uniform_int_distribution<uint32_t>> distribution;
};

using AllTypes = ::testing::Types<uint16_t, uint32_t>;
TYPED_TEST_CASE(SimdTestFnt, AllTypes);

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

#endif
