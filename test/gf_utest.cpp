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
#include <gtest/gtest.h>

#include "gf_bin_ext.h"
#include "gf_nf4.h"

namespace gf = quadiron::gf;

template <typename T>
class GfTestCommon : public ::testing::Test {
  public:
    void test_negation(const gf::Field<T>& gf)
    {
        for (int i = 0; i < 100; i++) {
            const T x = gf.weak_rand();
            const T y = gf.neg(x);
            ASSERT_EQ(gf.add(x, y), 0);
        }
    }

    void test_reciprocal(const gf::Field<T>& gf)
    {
        int n_found = 0;

        for (int i = 0; i < 100; i++) {
            const T x = gf.weak_rand();
            try {
                const T y = gf.inv(x);
                ASSERT_EQ(gf.mul(x, y), gf.get_unit());
                n_found++;
            } catch (const quadiron::Exception&) {
                continue;
            }
        }
        ASSERT_GT(n_found, 0);
    }

    void test_log(const gf::Field<T>& gf)
    {
        int n_found = 0;

        for (int i = 0; i < 1000; i++) {
            const T x = gf.weak_rand();
            const T y = gf.weak_rand();
            try {
                const T z = gf.log(x, y);
                const T t = gf.exp(x, z);

                ASSERT_EQ(t, y);
                n_found++;
            } catch (...) {
                continue;
            }
        }
        ASSERT_GT(n_found, 0);
    }

    void test_find_primitive_root(gf::Field<T>* gf)
    {
        gf->find_primitive_root();
        ASSERT_TRUE(gf->check_primitive_root(gf->get_root()));
    }

    void test_get_order(const gf::Field<T>& gf)
    {
        const T h = gf.card_minus_one();

        for (int i = 0; i < 1000; i++) {
            const T x = gf.weak_rand();
            const T order = gf.get_order(x);

            ASSERT_EQ(gf.exp(x, order), 1);
            ASSERT_EQ(h % order, 0);
        }
    }

    void test_get_nth_root(const gf::Field<T>& gf)
    {
        for (int i = 0; i < 1000; i++) {
            const T x = gf.weak_rand();
            const T nth_root = gf.get_nth_root(x);

            ASSERT_EQ(gf.exp(nth_root, x), 1);
        }
    }
};

using AllTypes = ::testing::Types<uint32_t, uint64_t, __uint128_t>;
TYPED_TEST_CASE(GfTestCommon, AllTypes);

TYPED_TEST(GfTestCommon, TestGfNf4) // NOLINT
{
    for (unsigned n = 1; n <= sizeof(TypeParam) / 4; n++) {
        quadiron::prng().seed(time(0));

        auto gf(gf::create<gf::NF4<TypeParam>>(n));
        this->test_negation(gf);
        this->test_negation(gf);

        // Test pack/unpack.
        for (int i = 0; i < 100; i++) {
            const TypeParam x = gf.weak_rand_tuple();
            const quadiron::GroupedValues<TypeParam> z = gf.unpack(x);
            const TypeParam y = gf.pack(z.values, z.flag);
            ASSERT_EQ(x, y);
        }
    }
}

TYPED_TEST(GfTestCommon, TestGf2n) // NOLINT
{
    for (TypeParam n = 8; n <= 8 * sizeof(TypeParam); n *= 2) {
        quadiron::prng().seed(time(0));

        auto gf(gf::create<gf::BinExtension<TypeParam>>(n));
        this->test_negation(gf);
        this->test_reciprocal(gf);
        // this->test_log(gf);
        if (n < 128) { // It works slowly in GF2N(128).
            this->test_get_order(gf);
            this->test_get_nth_root(gf);
            this->test_find_primitive_root(&gf);
        }
    }
}

template <typename T>
class GfTestNo128 : public GfTestCommon<T> {
};

using No128 = ::testing::Types<uint32_t, uint64_t>;
TYPED_TEST_CASE(GfTestNo128, No128);

TYPED_TEST(GfTestNo128, TestGf5) // NOLINT
{
    quadiron::prng().seed(time(0));

    auto gf(gf::create<gf::Prime<TypeParam>>(5));
    this->test_negation(gf);
    this->test_reciprocal(gf);
    this->test_log(gf);
    this->test_get_order(gf);
    this->test_get_nth_root(gf);
    this->test_find_primitive_root(&gf);
}

TYPED_TEST(GfTestNo128, TestGf256) // NOLINT
{
    quadiron::prng().seed(time(0));

    auto gf(gf::create<gf::BinExtension<TypeParam>>(8));
    this->test_negation(gf);
    this->test_reciprocal(gf);
    this->test_log(gf);
    this->test_get_order(gf);
    this->test_get_nth_root(gf);
    this->test_find_primitive_root(&gf);
}
