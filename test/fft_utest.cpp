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

#include "fft_2n.h"
#include "fft_add.h"
#include "fft_ct.h"
#include "fft_gt.h"
#include "fft_naive.h"
#include "fft_single.h"
#include "gf_bin_ext.h"
#include "gf_prime.h"

namespace fft = quadiron::fft;
namespace gf = quadiron::gf;

template <typename T>
class FftTest : public ::testing::Test {
  public:
    const unsigned q = 65537;
    const unsigned code_len = 32;

    quadiron::vec::Vector<T>
    random_vec(const gf::Field<T>& gf, int size, unsigned to_init)
    {
        quadiron::vec::Vector<T> vec(gf, size);
        vec.zero_fill();

        for (unsigned i = 0; i < to_init; i++) {
            vec.set(i, gf.weak_rand());
        }

        return vec;
    }

    // Taylor expansion on (x^t - x).
    void test_taylor_expand(const gf::Field<T>& gf, fft::Additive<T>* fft)
    {
        for (int i = 0; i < 1000; i++) {
            const int t = 2 + gf.weak_rand() % (fft->get_n() - 1);
            const int n = t + gf.weak_rand() % (fft->get_n() - t + 1);

            quadiron::vec::Vector<T> v1 = this->random_vec(gf, n, n);

            int m = n / t;
            while (m * t < n) {
                m++;
            }
            quadiron::vec::Vector<T> v2(gf, t * m);

            fft->taylor_expand(v2, v1, n, t);
            quadiron::vec::Vector<T> _v1(gf, n);
            fft->inv_taylor_expand(_v1, v2, t);
            ASSERT_EQ(_v1, v1);
        }
    }

    // Taylor expansion on (x^2 - x).
    void test_taylor_expand_t2(const gf::Field<T>& gf, fft::Additive<T>* fft)
    {
        for (int i = 0; i < 1000; i++) {
            const int n = fft->get_n();
            quadiron::vec::Vector<T> v1 = this->random_vec(gf, n, n);

            fft->taylor_expand_t2(v1, n, true);
            quadiron::vec::Vector<T> _v1(gf, n);
            fft->inv_taylor_expand_t2(_v1);
            ASSERT_EQ(_v1, v1);
        }
    }

    void test_fft_codec(
        const gf::Field<T>& gf,
        fft::FourierTransform<T>* fft,
        int n_data)
    {
        quadiron::vec::Vector<T> _v(gf, fft->get_n());
        quadiron::vec::Vector<T> v2(gf, fft->get_n());
        for (int j = 0; j < 1000; j++) {
            quadiron::vec::Vector<T> v =
                this->random_vec(gf, fft->get_n(), n_data);

            fft->fft(_v, v);
            fft->ifft(v2, _v);
            ASSERT_EQ(v, v2);
        }
    }

    void test_fft_1vs1(
        const gf::Field<T>& gf,
        fft::FourierTransform<T>* fft1,
        fft::FourierTransform<T>* fft2,
        int n_data)
    {
        ASSERT_EQ(fft1->get_n(), fft2->get_n());

        quadiron::vec::Vector<T> fft_1(gf, fft1->get_n());
        quadiron::vec::Vector<T> fft_2(gf, fft1->get_n());
        quadiron::vec::Vector<T> ifft_1(gf, fft1->get_n());
        quadiron::vec::Vector<T> ifft_2(gf, fft1->get_n());
        for (int j = 0; j < 100; j++) {
            quadiron::vec::Vector<T> v =
                this->random_vec(gf, fft1->get_n(), n_data);

            fft1->fft(fft_1, v);
            fft2->fft(fft_2, v);

            ASSERT_EQ(fft_1, fft_2);
            fft1->ifft(ifft_1, fft_1);
            fft2->ifft(ifft_2, fft_2);

            ASSERT_EQ(ifft_1, ifft_2);
            ASSERT_EQ(ifft_1, v);
        }
    }
};

using TestedTypes = ::testing::Types<uint32_t, uint64_t>;
TYPED_TEST_CASE(FftTest, TestedTypes);

TYPED_TEST(FftTest, TestGcd) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(97));
    quadiron::SignedDoubleSizeVal<TypeParam> bezout[2];

    ASSERT_EQ(gf.inv(20), 34);
    for (int i = 0; i < 100; i++) {
        const TypeParam x = gf.weak_rand();

        ASSERT_EQ(
            quadiron::arith::extended_gcd<TypeParam>(97, x, bezout, nullptr),
            1);
        TypeParam y = gf.inv(x);
        if (bezout[1] < 0) {
            bezout[1] = bezout[1] + gf.card();
        }
        ASSERT_EQ(bezout[1], y);
    }
}

TYPED_TEST(FftTest, TestQuadraticResidues) // NOLINT
{
    auto gf32(gf::create<gf::Prime<TypeParam>>(32));
    for (int i = 0; i < 32; i++) {
        ASSERT_TRUE(gf32.is_quadratic_residue(gf32.exp(i, 2)));
    }

    auto gf7(gf::create<gf::Prime<TypeParam>>(7));
    ASSERT_TRUE(gf7.is_quadratic_residue(2));
    ASSERT_FALSE(gf7.is_quadratic_residue(5));

    auto gf8(gf::create<gf::Prime<TypeParam>>(8));
    ASSERT_TRUE(gf8.is_quadratic_residue(1));
    ASSERT_FALSE(gf8.is_quadratic_residue(3));
}

TYPED_TEST(FftTest, TestFftNaive) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();

    ASSERT_EQ(quadiron::arith::jacobi<TypeParam>(R, this->q), -1);

    // With this encoder we cannot exactly satisfy users request,
    // we need to pad n = minimal divisor of (q-1)
    // that is at least (n_parities + n_data).
    const unsigned n = gf.get_code_len(this->code_len);

    // Compute root of order n-1 such as r^(n-1) mod q == 1.
    const unsigned r = gf.get_nth_root(n);

    fft::Naive<TypeParam> fft(gf, n, r);
    this->test_fft_codec(gf, &fft, this->code_len);
}

TYPED_TEST(FftTest, TestNaiveVsFft2kVec) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();

    ASSERT_EQ(quadiron::arith::jacobi<TypeParam>(R, this->q), -1);

    // With this encoder we cannot exactly satisfy users request,
    // we need to pad n = minimal divisor of (q-1)
    // that is at least (n_parities + n_data).
    const unsigned n = gf.get_code_len(this->code_len);

    // Compute root of order n-1 such as r^(n-1) mod q == 1.
    const unsigned r = gf.get_nth_root(n);

    fft::Naive<TypeParam> fft_naive(gf, n, r);
    fft::Radix2<TypeParam> fft_2n(gf, n, n);

    this->test_fft_1vs1(gf, &fft_naive, &fft_2n, this->code_len);
}
TYPED_TEST(FftTest, TestFft2kVec) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();

    ASSERT_EQ(quadiron::arith::jacobi<TypeParam>(R, this->q), -1);

    // With this encoder we cannot exactly satisfy users request,
    // we need to pad n = minimal divisor of (q-1)
    // that is at least (n_parities + n_data).
    const unsigned n = gf.get_code_len(this->code_len);

    for (unsigned data_len = 2; data_len <= n; data_len *= 2) {
        fft::Radix2<TypeParam> fft(gf, n, data_len);
        quadiron::vec::Vector<TypeParam> _v(gf, fft.get_n());
        quadiron::vec::Vector<TypeParam> v2(gf, fft.get_n());
        for (unsigned len = 2; len < n; len *= 2) {
            for (int j = 0; j < 100; j++) {
                quadiron::vec::Vector<TypeParam> v =
                    this->random_vec(gf, len, len);

                fft.fft(_v, v);
                fft.ifft(v2, _v);
                quadiron::vec::Slice<TypeParam> _v2(&v2, len);
                ASSERT_EQ(v, _v2);
            }
        }
    }
}

TYPED_TEST(FftTest, TestFft2kVecp) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();
    const size_t size = 4;

    ASSERT_EQ(quadiron::arith::jacobi<TypeParam>(R, this->q), -1);

    // With this encoder we cannot exactly satisfy users request,
    // we need to pad n = minimal divisor of (q-1)
    // that is at least (n_parities + n_data).
    const unsigned n = gf.get_code_len(this->code_len);

    for (unsigned data_len = 2; data_len <= n; data_len *= 2) {
        fft::Radix2<TypeParam> fft(gf, n, data_len, size);

        const int vec_n = fft.get_n();
        quadiron::vec::Buffers<TypeParam> v2(vec_n, size);
        quadiron::vec::Buffers<TypeParam> _v2(vec_n, size);
        for (unsigned len = 2; len <= n; len *= 2) {
            quadiron::vec::Buffers<TypeParam> v(len, size);
            quadiron::vec::Buffers<TypeParam> _v(_v2, 0, len);
            for (int j = 0; j < 100; j++) {
                for (unsigned i = 0; i < len; i++) {
                    TypeParam* mem = v.get(i);
                    for (size_t u = 0; u < size; u++) {
                        mem[u] = gf.weak_rand();
                    }
                }
                fft.fft(v2, v);
                fft.ifft(_v2, v2);

                ASSERT_EQ(v, _v);
            }
        }
    }
}

TYPED_TEST(FftTest, TestNaiveVsFft2kVecp) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();
    const size_t size = 2;

    ASSERT_EQ(quadiron::arith::jacobi<TypeParam>(R, this->q), -1);

    // With this encoder we cannot exactly satisfy users request,
    // we need to pad n = minimal divisor of (q-1)
    // that is at least (n_parities + n_data).
    const unsigned n = gf.get_code_len(this->code_len);

    // Compute root of order n-1 such as r^(n-1) mod q == 1.
    const unsigned r = gf.get_nth_root(n);

    fft::Naive<TypeParam> fft_naive(gf, n, r, size);
    fft::Radix2<TypeParam> fft_2n(gf, n, n, size);

    ASSERT_EQ(fft_naive.get_n(), fft_2n.get_n());

    quadiron::vec::Buffers<TypeParam> v(n, size);
    quadiron::vec::Buffers<TypeParam> fft1(n, size);
    quadiron::vec::Buffers<TypeParam> fft2(n, size);
    quadiron::vec::Buffers<TypeParam> ifft1(n, size);
    quadiron::vec::Buffers<TypeParam> ifft2(n, size);
    for (int j = 0; j < 100; j++) {
        for (unsigned i = 0; i < n; i++) {
            TypeParam* mem = v.get(i);
            for (size_t u = 0; u < size; u++) {
                mem[u] = gf.weak_rand();
            }
        }

        fft_naive.fft(fft1, v);
        fft_2n.fft(fft2, v);

        ASSERT_EQ(fft1, fft2);

        fft_naive.ifft(ifft1, fft1);
        fft_2n.ifft(ifft2, fft2);

        ASSERT_EQ(ifft1, ifft2);
        ASSERT_EQ(ifft1, v);
    }
}

TYPED_TEST(FftTest, TestFftGt) // NOLINT
{
    auto gf(gf::create<gf::BinExtension<TypeParam>>(16));

    // With this encoder we cannot exactly satisfy users request,
    // we need to pad n = minimal divisor of (q-1)
    // that is at least (this->n_parities + n_data).
    const TypeParam n = gf.get_code_len(this->code_len);

    fft::GoodThomas<TypeParam> fft(gf, n);
    this->test_fft_codec(gf, &fft, this->code_len);
}

TYPED_TEST(FftTest, TestFftCtGfp) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));

    // With this encoder we cannot exactly satisfy users request,
    // We need to pad n = minimal divisor of (q-1)
    // that is at least (n_parities + n_data).
    const TypeParam n = gf.get_code_len(this->code_len);

    fft::CooleyTukey<TypeParam> fft(gf, n);
    this->test_fft_codec(gf, &fft, this->code_len);
}

TYPED_TEST(FftTest, TestFftCtGf2n) // NOLINT
{
    const size_t max_n = 8 * sizeof(TypeParam);
    for (size_t gf_n = 4; gf_n <= 128 && gf_n <= max_n; gf_n *= 2) {
        auto gf(gf::create<gf::BinExtension<TypeParam>>(gf_n));

        unsigned len = this->code_len;
        if (gf.card_minus_one() <= this->code_len) {
            len = gf.weak_rand();
        }

        // With this encoder we cannot exactly satisfy users request,
        // we need to pad n = minimal divisor of (q-1)
        // that is at least (n_parities + n_data).
        const TypeParam n = gf.get_code_len(len);

        fft::CooleyTukey<TypeParam> fft(gf, n);
        this->test_fft_codec(gf, &fft, len);
    }
}

TYPED_TEST(FftTest, TestFftAdd) // NOLINT
{
    for (size_t gf_n = 4; gf_n <= 128 && gf_n <= 8 * sizeof(TypeParam);
         gf_n *= 2) {
        auto gf(gf::create<gf::BinExtension<TypeParam>>(gf_n));

        unsigned len = this->code_len;
        if (gf.card_minus_one() <= this->code_len) {
            len = gf_n;
        }

        // n is power of 2 and at least n_data + n_parities.
        const int n = quadiron::arith::ceil2<TypeParam>(len);
        const int m = quadiron::arith::log2<TypeParam>(n);
        fft::Additive<TypeParam> fft(gf, m);

        this->test_taylor_expand(gf, &fft);
        this->test_taylor_expand_t2(gf, &fft);
        this->test_fft_codec(gf, &fft, len);
    }
}

TYPED_TEST(FftTest, TestFftNaive2) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();

    ASSERT_EQ(quadiron::arith::jacobi<TypeParam>(R, this->q), -1);

    // The test is spefically for length-8 FFT
    const unsigned n = 8;

    // compute root of order n-1 such as r^(n-1) mod q == 1
    const unsigned r = gf.get_nth_root(n);

    fft::Naive<TypeParam> fft(gf, n, r);

    quadiron::vec::Vector<TypeParam> v(gf, fft.get_n());
    quadiron::vec::Vector<TypeParam> _v(gf, fft.get_n());
    quadiron::vec::Vector<TypeParam> v2(gf, fft.get_n());
    v.zero_fill();
    v.set(0, 27746);
    v.set(1, 871);
    v.set(2, 49520);

    fft.fft(_v, v);
    ASSERT_EQ(_v.get(0), 12600);
    ASSERT_EQ(_v.get(1), 27885);
    ASSERT_EQ(_v.get(2), 17398);
    ASSERT_EQ(_v.get(3), 4624);
    ASSERT_EQ(_v.get(4), 10858);
    ASSERT_EQ(_v.get(5), 36186);
    ASSERT_EQ(_v.get(6), 4591);
    ASSERT_EQ(_v.get(7), 42289);

    fft.ifft(v2, _v);
    ASSERT_EQ(v, v2);
}

TYPED_TEST(FftTest, TestFft2Gfp) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));

    fft::Size2<TypeParam> fft(gf);
    this->test_fft_codec(gf, &fft, 2);
}

TYPED_TEST(FftTest, TestFftSingleGfp) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const int n = quadiron::arith::ceil2<TypeParam>(this->code_len);
    fft::Single<TypeParam> fft(gf, n);

    this->test_fft_codec(gf, &fft, 1);
}

TYPED_TEST(FftTest, TestFftGf2n) // NOLINT
{
    const size_t max_n = 8 * sizeof(TypeParam);
    for (TypeParam gf_n = 4; gf_n <= 128 && gf_n <= max_n; gf_n *= 2) {
        auto gf(gf::create<gf::BinExtension<TypeParam>>(gf_n));
        const TypeParam R = gf.get_primitive_root();

        ASSERT_EQ(gf.exp(R, gf.card_minus_one()), 1);

        unsigned len = this->code_len;
        if (gf.card_minus_one() <= this->code_len) {
            len = gf.weak_rand();
        }

        // With this encoder we cannot exactly satisfy users request,
        // we need to pad.
        const unsigned n = gf.get_code_len(len);

        const TypeParam r = gf.get_nth_root(n);
        ASSERT_EQ(gf.exp(r, n), 1);

        fft::Naive<TypeParam> fft(gf, n, r);
        this->test_fft_codec(gf, &fft, len);
    }
}
