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
#include "fft_large.h"
#include "fft_naive.h"
#include "fft_single.h"
#include "gf_bin_ext.h"
#include "gf_prime.h"
#include "misc.h"
#include "vec_poly.h"

namespace fft = quadiron::fft;
namespace gf = quadiron::gf;
namespace arith = quadiron::arith;
namespace vec = quadiron::vec;

template <typename T>
class FftTest : public ::testing::Test {
  public:
    const unsigned q = 65537;
    const std::vector<unsigned> code_lengths = {2, 4, 32, 64};
    const std::vector<int> vec_t = {2, 3, 4, 8, 10};
    const std::vector<int> vec_n_minus_t = {-1, 0, 1, 2, 10};
    const std::vector<int> vec_n_power_2 = {1, 2, 4, 8, 16, 32};

    vec::Vector<T>
    random_vec(const gf::Field<T>& gf, int size, unsigned to_init)
    {
        vec::Vector<T> vec(gf, size);
        vec.zero_fill();

        for (unsigned i = 0; i < to_init; i++) {
            vec.set(i, gf.rand());
        }

        return vec;
    }

    // Taylor expansion on (x^t - x).
    void test_taylor_expand(const gf::Field<T>& gf, fft::Additive<T>* fft)
    {
        for (int i = 0; i < 100; i++) {
            for (auto const& t : this->vec_t) {
                for (auto const& n_minus_t : this->vec_n_minus_t) {
                    int n = t + n_minus_t;
                    if (n > fft->get_n()) {
                        n = fft->get_n();
                    }

                    vec::Vector<T> v1(this->random_vec(gf, n, n));

                    int m = n / t;
                    while (m * t < n) {
                        m++;
                    }
                    vec::Vector<T> v2(gf, t * m);

                    fft->taylor_expand(v2, v1, n, t);
                    vec::Vector<T> _v1(gf, n);
                    fft->inv_taylor_expand(_v1, v2, t);
                    ASSERT_EQ(_v1, v1);
                }
            }
        }
    }

    // Taylor expansion on (x^2 - x).
    void test_taylor_expand_t2(const gf::Field<T>& gf, fft::Additive<T>* fft)
    {
        for (int i = 0; i < 100; i++) {
            for (auto n : this->vec_n_power_2) {
                if (n > fft->get_n()) {
                    n = fft->get_n();
                }

                vec::Vector<T> v1(this->random_vec(gf, n, n));

                fft->taylor_expand_t2(v1, n, true);
                vec::Vector<T> _v1(gf, n);
                fft->inv_taylor_expand_t2(_v1);
                ASSERT_EQ(_v1, v1);
            }
        }
    }

    void test_fft_codec(
        const gf::Field<T>& gf,
        fft::FourierTransform<T>* fft,
        int n_data)
    {
        vec::Vector<T> _v(gf, fft->get_n());
        vec::Vector<T> v2(gf, fft->get_n());
        for (int j = 0; j < 1000; j++) {
            vec::Vector<T> v(this->random_vec(gf, fft->get_n(), n_data));

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

        vec::Vector<T> fft_1(gf, fft1->get_n());
        vec::Vector<T> fft_2(gf, fft1->get_n());
        vec::Vector<T> ifft_1(gf, fft1->get_n());
        vec::Vector<T> ifft_2(gf, fft1->get_n());
        for (int j = 0; j < 100; j++) {
            vec::Vector<T> v(this->random_vec(gf, fft1->get_n(), n_data));

            fft1->fft(fft_1, v);
            fft2->fft(fft_2, v);

            ASSERT_EQ(fft_1, fft_2);
            fft1->ifft(ifft_1, fft_1);
            fft2->ifft(ifft_2, fft_2);

            ASSERT_EQ(ifft_1, ifft_2);
            ASSERT_EQ(ifft_1, v);
        }
    }

    /** Convert a number into a vector of digits padded with zeros
     *
     * @param vec destination vector
     * @param n vector size
     * @param num number in string representation
     */
    void big_num_to_vec(vec::Vector<T>& vec, size_t n, std::string num)
    {
        size_t i;
        const size_t len = num.length();
        for (i = 0; i < len; i++) {
            vec.set(i, num[len - i - 1] - '0');
        }
        for (; i < n; i++) {
            vec.set(i, 0);
        }
    }

    /** Compute sum of 2 large numbers in string representation
     *
     * @param str1 first string
     * @param str2 second string
     *
     * @return the sum string
     */
    std::string big_num_add(std::string str1, std::string str2)
    {
        if (str1.length() > str2.length()) {
            std::swap(str1, str2);
        }

        std::string str;

        size_t n1 = str1.length();
        size_t n2 = str2.length();

        std::reverse(str1.begin(), str1.end());
        std::reverse(str2.begin(), str2.end());

        int carry = 0;
        for (size_t i = 0; i < n1; i++) {
            int sum = ((str1[i] - '0') + (str2[i] - '0') + carry);
            str.push_back(sum % 10 + '0');
            carry = sum / 10;
        }

        for (size_t i = n1; i < n2; i++) {
            int sum = ((str2[i] - '0') + carry);
            str.push_back(sum % 10 + '0');
            carry = sum / 10;
        }

        if (carry) {
            str.push_back(carry + '0');
        }

        std::reverse(str.begin(), str.end());

        return str;
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
        const TypeParam x = gf.rand();

        ASSERT_EQ(arith::extended_gcd<TypeParam>(97, x, bezout, nullptr), 1);
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

    ASSERT_EQ(arith::jacobi<TypeParam>(R, this->q), -1);

    for (auto const& code_len : this->code_lengths) {
        // With this encoder we cannot exactly satisfy users request,
        // we need to pad n = minimal divisor of (q-1)
        // that is at least (n_parities + n_data).
        const unsigned n = gf.get_code_len(code_len);

        // Compute root of order n-1 such as r^(n-1) mod q == 1.
        const unsigned r = gf.get_nth_root(n);

        fft::Naive<TypeParam> fft(gf, n, r);
        this->test_fft_codec(gf, &fft, code_len);
    }
}

TYPED_TEST(FftTest, TestNaiveVsFft2kVec) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();

    ASSERT_EQ(arith::jacobi<TypeParam>(R, this->q), -1);

    for (auto const& code_len : this->code_lengths) {
        // With this encoder we cannot exactly satisfy users request,
        // we need to pad n = minimal divisor of (q-1)
        // that is at least (n_parities + n_data).
        const unsigned n = gf.get_code_len(code_len);

        // Compute root of order n-1 such as r^(n-1) mod q == 1.
        const unsigned r = gf.get_nth_root(n);

        fft::Naive<TypeParam> fft_naive(gf, n, r);
        fft::Radix2<TypeParam> fft_2n(gf, n, n);

        this->test_fft_1vs1(gf, &fft_naive, &fft_2n, code_len);
    }
}
TYPED_TEST(FftTest, TestFft2kVec) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();

    ASSERT_EQ(arith::jacobi<TypeParam>(R, this->q), -1);

    for (auto const& code_len : this->code_lengths) {
        // With this encoder we cannot exactly satisfy users request,
        // we need to pad n = minimal divisor of (q-1)
        // that is at least (n_parities + n_data).
        const unsigned n = gf.get_code_len(code_len);

        for (unsigned data_len = 2; data_len <= n; data_len *= 2) {
            fft::Radix2<TypeParam> fft(gf, n, data_len);
            vec::Vector<TypeParam> _v(gf, fft.get_n());
            vec::Vector<TypeParam> v2(gf, fft.get_n());
            for (unsigned len = 2; len < n; len *= 2) {
                for (int j = 0; j < 100; j++) {
                    vec::Vector<TypeParam> v(this->random_vec(gf, len, len));

                    fft.fft(_v, v);
                    fft.ifft(v2, _v);
                    vec::Slice<TypeParam> _v2(&v2, len);
                    ASSERT_EQ(v, _v2);
                }
            }
        }
    }
}

TYPED_TEST(FftTest, TestFft2kVecp) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();
    const size_t size = 4;

    ASSERT_EQ(arith::jacobi<TypeParam>(R, this->q), -1);

    for (auto const& code_len : this->code_lengths) {
        // With this encoder we cannot exactly satisfy users request,
        // we need to pad n = minimal divisor of (q-1)
        // that is at least (n_parities + n_data).
        const unsigned n = gf.get_code_len(code_len);

        for (unsigned data_len = 2; data_len <= n; data_len *= 2) {
            fft::Radix2<TypeParam> fft(gf, n, data_len, size);

            const int vec_n = fft.get_n();
            vec::Buffers<TypeParam> v2(vec_n, size);
            vec::Buffers<TypeParam> _v2(vec_n, size);
            for (unsigned len = 2; len <= n; len *= 2) {
                vec::Buffers<TypeParam> v(len, size);
                vec::Buffers<TypeParam> _v(_v2, 0, len);
                for (int j = 0; j < 100; j++) {
                    for (unsigned i = 0; i < len; i++) {
                        TypeParam* mem = v.get(i);
                        for (size_t u = 0; u < size; u++) {
                            mem[u] = gf.rand();
                        }
                    }
                    fft.fft(v2, v);
                    fft.ifft(_v2, v2);

                    ASSERT_EQ(v, _v);
                }
            }
        }
    }
}

TYPED_TEST(FftTest, TestNaiveVsFft2kVecp) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();
    const size_t size = 2;

    ASSERT_EQ(arith::jacobi<TypeParam>(R, this->q), -1);

    for (auto const& code_len : this->code_lengths) {
        // With this encoder we cannot exactly satisfy users request,
        // we need to pad n = minimal divisor of (q-1)
        // that is at least (n_parities + n_data).
        const unsigned n = gf.get_code_len(code_len);

        // Compute root of order n-1 such as r^(n-1) mod q == 1.
        const unsigned r = gf.get_nth_root(n);

        fft::Naive<TypeParam> fft_naive(gf, n, r, size);
        fft::Radix2<TypeParam> fft_2n(gf, n, n, size);

        ASSERT_EQ(fft_naive.get_n(), fft_2n.get_n());

        vec::Buffers<TypeParam> v(n, size);
        vec::Buffers<TypeParam> fft1(n, size);
        vec::Buffers<TypeParam> fft2(n, size);
        vec::Buffers<TypeParam> ifft1(n, size);
        vec::Buffers<TypeParam> ifft2(n, size);
        for (int j = 0; j < 100; j++) {
            for (unsigned i = 0; i < n; i++) {
                TypeParam* mem = v.get(i);
                for (size_t u = 0; u < size; u++) {
                    mem[u] = gf.rand();
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
}

TYPED_TEST(FftTest, TestFftGt) // NOLINT
{
    auto gf(gf::create<gf::BinExtension<TypeParam>>(16));

    for (auto const& code_len : this->code_lengths) {
        // With this encoder we cannot exactly satisfy users request,
        // we need to pad n = minimal divisor of (q-1)
        // that is at least (this->n_parities + n_data).
        const TypeParam n = gf.get_code_len(code_len);

        fft::GoodThomas<TypeParam> fft(gf, n);
        this->test_fft_codec(gf, &fft, code_len);
    }
}

TYPED_TEST(FftTest, TestFftCtGfp) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));

    for (auto const& code_len : this->code_lengths) {
        // With this encoder we cannot exactly satisfy users request,
        // We need to pad n = minimal divisor of (q-1)
        // that is at least (n_parities + n_data).
        const TypeParam n = gf.get_code_len(code_len);

        fft::CooleyTukey<TypeParam> fft(gf, n);
        this->test_fft_codec(gf, &fft, code_len);
    }
}

TYPED_TEST(FftTest, TestFftCtGf2n) // NOLINT
{
    const size_t max_n = 8 * sizeof(TypeParam);
    for (size_t gf_n = 4; gf_n <= 128 && gf_n <= max_n; gf_n *= 2) {
        auto gf(gf::create<gf::BinExtension<TypeParam>>(gf_n));

        for (auto const& code_len : this->code_lengths) {
            unsigned len = code_len;
            if (gf.card_minus_one() <= code_len) {
                len = gf.rand();
            }

            // With this encoder we cannot exactly satisfy users request,
            // we need to pad n = minimal divisor of (q-1)
            // that is at least (n_parities + n_data).
            const TypeParam n = gf.get_code_len(len);

            fft::CooleyTukey<TypeParam> fft(gf, n);
            this->test_fft_codec(gf, &fft, len);
        }
    }
}

TYPED_TEST(FftTest, TestFftAdd) // NOLINT
{
    for (size_t gf_n = 4; gf_n <= 128 && gf_n <= 8 * sizeof(TypeParam);
         gf_n *= 2) {
        auto gf(gf::create<gf::BinExtension<TypeParam>>(gf_n));

        for (auto const& code_len : this->code_lengths) {
            unsigned len = code_len;
            if (gf.card_minus_one() <= code_len) {
                len = gf_n;
            }

            // n is power of 2 and at least n_data + n_parities.
            const int n = arith::ceil2<TypeParam>(len);
            const int m = arith::log2<TypeParam>(n);
            fft::Additive<TypeParam> fft(gf, m);

            this->test_taylor_expand(gf, &fft);
            this->test_taylor_expand_t2(gf, &fft);
            this->test_fft_codec(gf, &fft, len);
        }
    }
}

TYPED_TEST(FftTest, TestFftNaive2) // NOLINT
{
    auto gf(gf::create<gf::Prime<TypeParam>>(this->q));
    const unsigned R = gf.get_primitive_root();

    ASSERT_EQ(arith::jacobi<TypeParam>(R, this->q), -1);

    // The test is spefically for length-8 FFT
    const unsigned n = 8;

    // compute root of order n-1 such as r^(n-1) mod q == 1
    const unsigned r = gf.get_nth_root(n);

    fft::Naive<TypeParam> fft(gf, n, r);

    vec::Vector<TypeParam> v(gf, fft.get_n());
    vec::Vector<TypeParam> _v(gf, fft.get_n());
    vec::Vector<TypeParam> v2(gf, fft.get_n());
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

    for (auto const& code_len : this->code_lengths) {
        const int n = arith::ceil2<TypeParam>(code_len);
        fft::Single<TypeParam> fft(gf, n);

        this->test_fft_codec(gf, &fft, 1);
    }
}

TYPED_TEST(FftTest, TestFftGf2n) // NOLINT
{
    const size_t max_n = 8 * sizeof(TypeParam);
    for (TypeParam gf_n = 4; gf_n <= 128 && gf_n <= max_n; gf_n *= 2) {
        auto gf(gf::create<gf::BinExtension<TypeParam>>(gf_n));
        const TypeParam R = gf.get_primitive_root();

        ASSERT_EQ(gf.exp(R, gf.card_minus_one()), 1);

        for (auto const& code_len : this->code_lengths) {
            unsigned len = code_len;
            if (gf.card_minus_one() <= code_len) {
                len = gf.rand();
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
}

TYPED_TEST(FftTest, TestFftLarge) // NOLINT
{
    const unsigned n = 256;
    const unsigned q = 7681;

    std::unique_ptr<gf::Field<TypeParam>> gf = nullptr;
    std::unique_ptr<vec::Vector<TypeParam>> a = nullptr;
    std::unique_ptr<vec::Vector<TypeParam>> _a1 = nullptr;
    std::unique_ptr<vec::Vector<TypeParam>> _a2 = nullptr;
    std::unique_ptr<fft::FourierTransform<TypeParam>> fft = nullptr;
    std::unique_ptr<fft::FourierTransform<TypeParam>> fft2 = nullptr;

    gf = gf::alloc<gf::Field<TypeParam>, gf::Prime<TypeParam>>(q);

    TypeParam r = gf->get_nth_root(n);
    int l = arith::log2<TypeParam>(n);
    fft = std::make_unique<fft::Naive<TypeParam>>(*gf, n, r);
    fft2 = std::make_unique<fft::Large<TypeParam>>(*gf, l, r);

    a = std::make_unique<vec::Vector<TypeParam>>(*gf, n);
    _a1 = std::make_unique<vec::Vector<TypeParam>>(*gf, n);
    _a2 = std::make_unique<vec::Vector<TypeParam>>(*gf, n);
    a->rand();

    fft->fft(*_a1, *a);
    fft2->fft(*_a2, *a);

    ASSERT_EQ(*_a1 == *_a2, true);
}

TYPED_TEST(FftTest, TestFftEquivalence) // NOLINT
{
    std::vector<uint64_t> qs({7681, 12289, 65537});
    std::vector<uint64_t> ns({256, 512, 1024});

    for (auto q : qs) {
        for (auto n : ns) {
            ASSERT_EQ(arith::is_prime(q), true);

            std::unique_ptr<gf::Field<TypeParam>> gf = nullptr;
            std::unique_ptr<vec::Vector<TypeParam>> a = nullptr;
            std::unique_ptr<vec::Vector<TypeParam>> _a1 = nullptr;
            std::unique_ptr<vec::Vector<TypeParam>> _a2 = nullptr;
            std::unique_ptr<vec::Vector<TypeParam>> _a3 = nullptr;
            std::unique_ptr<fft::FourierTransform<TypeParam>> fft = nullptr;
            std::unique_ptr<fft::FourierTransform<TypeParam>> fft2 = nullptr;
            std::unique_ptr<fft::FourierTransform<TypeParam>> fft3 = nullptr;

            gf = gf::alloc<gf::Field<TypeParam>, gf::Prime<TypeParam>>(q);

            TypeParam r = gf->get_nth_root(n);
            int l = arith::log2<TypeParam>(n);
            fft = std::make_unique<fft::Naive<TypeParam>>(*gf, n, r);
            try {
                fft2 = std::make_unique<fft::Large<TypeParam>>(*gf, l, r);
            } catch (quadiron::Exception&) {
                fft2 = nullptr;
            }
            try {
                fft3 =
                    std::make_unique<fft::Radix2<TypeParam>>(*gf, n, n, 1024);
            } catch (quadiron::Exception&) {
                fft3 = nullptr;
            }

            a = std::make_unique<vec::Vector<TypeParam>>(*gf, n);
            _a1 = std::make_unique<vec::Vector<TypeParam>>(*gf, n);
            _a2 = std::make_unique<vec::Vector<TypeParam>>(*gf, n);
            _a3 = std::make_unique<vec::Vector<TypeParam>>(*gf, n);
            a->rand();

            uint64_t t1 = 0;
            uint64_t t2 = 0;
            uint64_t t3 = 0;

            t1 = quadiron::hw_timer();
            fft->fft(*_a1, *a);
            t1 = quadiron::hw_timer() - t1;
            if (fft2) {
                t2 = quadiron::hw_timer();
                fft2->fft(*_a2, *a);
                t2 = quadiron::hw_timer() - t2;
            }
            if (fft3) {
                t3 = quadiron::hw_timer();
                fft3->fft(*_a3, *a);
                t3 = quadiron::hw_timer() - t3;
            }
            if (fft2) {
                ASSERT_EQ(*_a1 == *_a2, true);
            }
            if (fft3) {
                ASSERT_EQ(*_a1 == *_a3, true);
            }
        }
    }
}

TYPED_TEST(FftTest, TestFftLarge2) // NOLINT
{
    // Schönhage-Strassen algorithm
    // Example taken from Pierre Meunier's book
    if (sizeof(TypeParam) < 8)
        return;

    const int base = 10;
    const int p = 14; // we could multiply integers of 2^p digits

    uint64_t l = p + 1;

    // choose 2 prime numbers of the form p=a.2^n+1
    // because if x is not a quadratic residue then w=x^a is
    // a 2^n-th principal root of unity in GF_p
    uint64_t a1 = 2;
    uint64_t a2 = 5;
    uint64_t p1 = a1 * arith::exp<TypeParam>(2, 15) + 1;
    uint64_t p2 = a2 * arith::exp<TypeParam>(2, 15) + 1;
    ASSERT_EQ(arith::is_prime<TypeParam>(p1), true);
    ASSERT_EQ(arith::is_prime<TypeParam>(p2), true);

    // ensure their product is bounded (b-1)^2*2^(n-1) < m
    uint64_t m = p1 * p2;
    // check overflow
    assert(m / p1 == p2);
    assert(
        arith::exp<TypeParam>((base - 1), 2) * arith::exp<TypeParam>(p, 2) < m);

    // find x so it is not a quadratic residue in GF_p1 and GF_p2
    assert(arith::jacobi<TypeParam>(3, p1) == arith::jacobi<TypeParam>(p1, 3));
    assert(arith::jacobi<TypeParam>(p1, 3) == arith::jacobi<TypeParam>(2, 3));
    assert(arith::jacobi<TypeParam>(3, p2) == arith::jacobi<TypeParam>(p2, 3));
    assert(arith::jacobi<TypeParam>(p2, 3) == arith::jacobi<TypeParam>(2, 3));
    assert(arith::jacobi<TypeParam>(2, 3) == -1);
    // which means x=3 is not a quadratic residue in GF_p1 and GF_p2

    // therefore we can compute 2^n-th roots of unity in GF_p1 and GF_p2
    uint64_t w1 = arith::exp<TypeParam>(3, a1);
    uint64_t w2 = arith::exp<TypeParam>(3, a2);
    assert(w1 == 9);
    assert(w2 == 243);

    // find root of unity in GF_p1p2
    uint64_t _a[2];
    uint64_t _n[2];
    _a[0] = w1;
    _n[0] = p1;
    _a[1] = w2;
    _n[1] = p2;
    uint64_t w = arith::chinese_remainder<uint64_t>(2, _a, _n);
    ASSERT_EQ(w, 25559439);

    std::unique_ptr<gf::Field<TypeParam>> gf_m =
        gf::alloc<gf::Field<TypeParam>, gf::Prime<TypeParam>>(m);
    fft::Large<TypeParam> fft(*gf_m, l, 25559439);

    // parse the big numbers
    std::string X = "1236548787985654354598651354984132468";
    std::string Y = "745211515185321545554545854598651354984132468";

    std::unique_ptr<vec::Vector<TypeParam>> _X = nullptr;
    _X = std::make_unique<vec::Vector<TypeParam>>(*gf_m, fft.get_n());
    this->big_num_to_vec(*_X, fft.get_n(), X);
    std::unique_ptr<vec::Vector<TypeParam>> _Y = nullptr;
    _Y = std::make_unique<vec::Vector<TypeParam>>(*gf_m, fft.get_n());
    this->big_num_to_vec(*_Y, fft.get_n(), Y);

    std::unique_ptr<vec::Vector<TypeParam>> sfX = nullptr;
    sfX = std::make_unique<vec::Vector<TypeParam>>(*gf_m, fft.get_n());
    std::unique_ptr<vec::Vector<TypeParam>> sfY = nullptr;
    sfY = std::make_unique<vec::Vector<TypeParam>>(*gf_m, fft.get_n());
    std::unique_ptr<vec::Vector<TypeParam>> _XY = nullptr;
    _XY = std::make_unique<vec::Vector<TypeParam>>(*gf_m, fft.get_n());
    std::unique_ptr<vec::Vector<TypeParam>> sfXY = nullptr;
    sfXY = std::make_unique<vec::Vector<TypeParam>>(*gf_m, fft.get_n());

    fft.fft(*sfX, *_X);
    fft.fft(*sfY, *_Y);

    for (int i = 0; i <= fft.get_n() - 1; i++) {
        quadiron::DoubleSizeVal<TypeParam> val =
            quadiron::DoubleSizeVal<TypeParam>(sfX->get(i)) * sfY->get(i);
        _XY->set(i, val % m);
    }

    fft.ifft(*sfXY, *_XY);

    // carry propagation
    std::string z("0");
    for (int i = 0; i <= fft.get_n() - 1; i++) {
        std::string b(std::to_string(sfXY->get(i)));
        if (b != std::string("0")) {
            b.append(i, '0');
            z = this->big_num_add(z, b);
        }
    }

    ASSERT_EQ(
        z,
        std::string("921490395895362412399910100421159322")
            + "712298564831565484737491129935640058571771024");
}
