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
#include <random>

#include <gtest/gtest.h>

#include "quadiron.h"

namespace arith = quadiron::arith;

template <typename T>
class ArithTestCommon : public ::testing::Test {
  public:
    T max;
    std::uniform_int_distribution<uint32_t> uniform_dist_max;

    ArithTestCommon()
    {
        if (sizeof(T) <= sizeof(uint64_t))
            max = std::numeric_limits<T>::max();
        else {
            max = 0;
            for (size_t i = 0; i < sizeof(T) * 8; i++) {
                max += static_cast<T>(1) << i;
            }
        }
        uniform_dist_max = std::uniform_int_distribution<uint32_t>(
            1, static_cast<uint32_t>(max - 1));

        quadiron::prng().seed(time(0));
    }
    ~ArithTestCommon() = default;

    void check_all_primes(const std::vector<T>& primes, bool distinct)
    {
        ASSERT_GT(primes.size(), 0);
        ASSERT_TRUE(arith::is_prime<T>(primes.at(0)));

        for (size_t j = 1; j != primes.size(); ++j) {
            ASSERT_TRUE(arith::is_prime<T>(primes.at(j)));
            if (distinct) {
                ASSERT_NE(primes.at(j - 1), primes.at(j));
            }
        }
    }

    void check_divisors(T nb, const std::vector<T>& divisors, bool proper)
    {
        if (proper && arith::is_prime<T>(nb)) {
            ASSERT_EQ(divisors.size(), 0);
        } else {
            ASSERT_GT(divisors.size(), 0);
        }

        for (auto v : divisors) {
            ASSERT_EQ(nb % v, 0);
            if (proper) {
                ASSERT_TRUE(arith::is_prime<T>(nb / v));
            }
        }
    }

    void
    check_prime_divisors(T nb, const std::vector<T>& divisors, bool coprime)
    {
        ASSERT_GT(divisors.size(), 0);

        for (size_t i = 0; i != divisors.size(); ++i) {
            ASSERT_EQ(nb % divisors.at(i), 0);

            if (!coprime) {
                ASSERT_TRUE(arith::is_prime<T>(divisors.at(i)));
            } else {
                for (size_t j = i + 1; j != divisors.size(); ++j) {
                    ASSERT_NE(divisors.at(i), divisors.at(j));
                }
            }
        }
    }
};

using AllTypes = ::testing::Types<uint32_t, uint64_t, __uint128_t>;
TYPED_TEST_CASE(ArithTestCommon, AllTypes);

TYPED_TEST(ArithTestCommon, TestBasicOperations) // NOLINT
{
    ASSERT_EQ(arith::sqrt<TypeParam>(2025), 45);

    ASSERT_TRUE(arith::is_prime<TypeParam>(2));
    ASSERT_TRUE(arith::is_prime<TypeParam>(3));
    ASSERT_TRUE(arith::is_prime<TypeParam>(13));
    ASSERT_TRUE(!arith::is_prime<TypeParam>(4));
    ASSERT_TRUE(!arith::is_prime<TypeParam>(15));
}

TYPED_TEST(ArithTestCommon, TestReciprocal) // NOLINT
{
    const TypeParam sub_max = arith::sqrt<TypeParam>(this->max);
    std::uniform_int_distribution<uint32_t> dis(1, sub_max - 1);

    for (int i = 0; i < 1000; i++) {
        const TypeParam x = dis(quadiron::prng());
        const TypeParam y = arith::exp<TypeParam>(x, 2);
        const TypeParam z = arith::sqrt<TypeParam>(y);

        ASSERT_EQ(z, x);
    }
}

TYPED_TEST(ArithTestCommon, TestJacobi) // NOLINT
{
    ASSERT_EQ(arith::jacobi<TypeParam>(1001, 9907), -1);
    ASSERT_EQ(arith::jacobi<TypeParam>(19, 45), 1);
    ASSERT_EQ(arith::jacobi<TypeParam>(8, 21), -1);
    ASSERT_EQ(arith::jacobi<TypeParam>(5, 21), 1);
    ASSERT_EQ(arith::jacobi<TypeParam>(47, 221), -1);
    ASSERT_EQ(arith::jacobi<TypeParam>(2, 221), -1);
}

TYPED_TEST(ArithTestCommon, TestFactorDistinctPrime) // NOLINT
{
    for (int i = 0; i < 1000; i++) {
        const TypeParam x = this->uniform_dist_max(quadiron::prng());
        const std::vector<TypeParam> primes =
            arith::factor_distinct_prime<TypeParam>(x);

        this->check_all_primes(primes, true);
    }
}

TYPED_TEST(ArithTestCommon, TestFactorPrime) // NOLINT
{
    for (int i = 0; i < 1000; i++) {
        std::vector<TypeParam> primes;
        std::vector<int> exponent;
        const TypeParam x = this->uniform_dist_max(quadiron::prng());

        arith::factor_prime<TypeParam>(x, &primes, &exponent);
        this->check_all_primes(primes, true);
        // Check primes exponent.
        ASSERT_GT(primes.size(), 0);
        ASSERT_GT(exponent.size(), 0);

        TypeParam y = 1;
        for (size_t j = 0; j != primes.size(); ++j) {
            y *= arith::exp<TypeParam>(primes.at(j), exponent.at(j));
        }
        ASSERT_EQ(y, x);
    }
}

TYPED_TEST(ArithTestCommon, TestGetProperDivisor) // NOLINT
{
    for (int i = 0; i < 1000; i++) {
        const TypeParam x = this->uniform_dist_max(quadiron::prng());
        const std::vector<TypeParam> divisors =
            arith::get_proper_divisors<TypeParam>(x);

        this->check_divisors(x, divisors, true);
    }
}

TYPED_TEST(ArithTestCommon, TestGetProperDivisor2) // NOLINT
{
    for (int i = 0; i < 1000; i++) {
        const TypeParam x = this->uniform_dist_max(quadiron::prng());
        const std::vector<TypeParam> factors =
            arith::factor_distinct_prime<TypeParam>(x);
        const std::vector<TypeParam> divisors =
            arith::get_proper_divisors<TypeParam>(x, factors);

        this->check_divisors(x, divisors, true);
    }
}

TYPED_TEST(ArithTestCommon, TestGetAllDivisor) // NOLINT
{
    for (int i = 0; i < 1000; i++) {
        const TypeParam x = this->uniform_dist_max(quadiron::prng());
        const std::vector<TypeParam> divisors =
            arith::get_all_divisors<TypeParam>(x);

        this->check_divisors(x, divisors, false);
    }
}

TYPED_TEST(ArithTestCommon, TestGetCodeLength) // NOLINT
{
    for (int i = 0; i < 1000; i++) {
        const TypeParam order = this->uniform_dist_max(quadiron::prng());
        std::uniform_int_distribution<uint32_t> dis(1, order - 1);
        const TypeParam n = dis(quadiron::prng());
        const TypeParam len = arith::get_code_len<TypeParam>(order, n);

        ASSERT_EQ(order % len, 0);
        ASSERT_GE(len, n);
    }
}

TYPED_TEST(ArithTestCommon, TestGetCodeLengthHighCompo) // NOLINT
{
    for (int i = 0; i < 1000; i++) {
        const TypeParam order = this->uniform_dist_max(quadiron::prng());
        std::uniform_int_distribution<uint32_t> dis(1, order - 1);
        const TypeParam n = dis(quadiron::prng());
        const TypeParam len =
            arith::get_code_len_high_compo<TypeParam>(order, n);

        ASSERT_EQ(order % len, 0);
        ASSERT_GE(len, n);
    }
}

TYPED_TEST(ArithTestCommon, TestGetCodeLengthHighCompo2) // NOLINT
{
    for (int i = 0; i < 1000; i++) {
        const TypeParam order = this->uniform_dist_max(quadiron::prng());
        std::uniform_int_distribution<uint32_t> dis(1, order - 1);
        const TypeParam n = dis(quadiron::prng());
        const std::vector<TypeParam> factors =
            arith::get_prime_factors<TypeParam>(order);
        const TypeParam len =
            arith::get_code_len_high_compo<TypeParam>(factors, n);

        ASSERT_EQ(order % len, 0);
        ASSERT_GE(len, n);
    }
}

TYPED_TEST(ArithTestCommon, TestGetCoprimeFactors) // NOLINT
{
    for (int i = 0; i < 1000; i++) {
        const TypeParam n = this->uniform_dist_max(quadiron::prng());
        const std::vector<TypeParam> divisors =
            arith::get_coprime_factors<TypeParam>(n);

        this->check_prime_divisors(n, divisors, true);
    }
}

TYPED_TEST(ArithTestCommon, TestGetPrimeFactors) // NOLINT
{
    for (int i = 0; i < 1000; i++) {
        const TypeParam n = this->uniform_dist_max(quadiron::prng());
        const std::vector<TypeParam> divisors =
            arith::get_prime_factors<TypeParam>(n);

        this->check_prime_divisors(n, divisors, false);
    }
}

template <typename T>
class ArithTestNo128 : public ::testing::Test {
};

using No128 = ::testing::Types<uint32_t, uint64_t>;
TYPED_TEST_CASE(ArithTestNo128, No128);

/* References:
 * http://www.math.unm.edu/~loring/links/discrete_f05/remainder.pdf
 * http://gauss.math.luc.edu/greicius/Math201/Fall2012/Lectures/ChineseRemainderThm.article.pdf
 */
TYPED_TEST(ArithTestNo128, TestChineseRemainder) // NOLINT
{
    TypeParam a[4];
    TypeParam n[4];
    TypeParam omega;

    a[0] = 4;
    n[0] = 107;
    a[1] = 2;
    n[1] = 74;
    omega = arith::chinese_remainder<TypeParam>(2, a, n);
    ASSERT_EQ(omega, 5996);

    a[0] = 6;
    n[0] = 7;
    a[1] = 4;
    n[1] = 8;
    omega = arith::chinese_remainder<TypeParam>(2, a, n);
    ASSERT_EQ(omega, 20);

    a[0] = 3;
    n[0] = 4;
    a[1] = 0;
    n[1] = 6;
    omega = arith::chinese_remainder<TypeParam>(2, a, n);
    // no solution XXX detect it
}

TYPED_TEST(ArithTestNo128, TestExtendedGcd) // NOLINT
{
    quadiron::SignedDoubleSizeVal<TypeParam> bezout[2];

    // Not explicitely related to GF(97).
    ASSERT_EQ(arith::extended_gcd<TypeParam>(240, 46, nullptr, nullptr), 2);
    ASSERT_EQ(arith::extended_gcd<TypeParam>(54, 24, nullptr, nullptr), 6);
    ASSERT_EQ(arith::extended_gcd<TypeParam>(210, 45, nullptr, nullptr), 15);
    ASSERT_EQ(arith::extended_gcd<TypeParam>(97, 20, bezout, nullptr), 1);
    ASSERT_TRUE(bezout[0] == -7 && bezout[1] == 34);
}

// Schönhage-Strassen algorithm (example taken from Pierre Meunier's book).
TEST(ArithTest, TestBignumMultiplication) // NOLINT
{
    const int b = 10; // Base.
    const int p = 14; // We could multiply integers of 2^p digits.

    // Choose two prime numbers of the form p=a.2^n+1
    // Because if `x` is not a quadratic residue then w=x^a is
    // a 2^n-th principal root of unity in GF_p
    const uint64_t a1 = 2;
    const uint64_t a2 = 5;
    const uint64_t p1 = a1 * arith::exp<uint64_t>(2, 15) + 1;
    const uint64_t p2 = a2 * arith::exp<uint64_t>(2, 15) + 1;
    ASSERT_TRUE(arith::is_prime<uint64_t>(p1));
    ASSERT_TRUE(arith::is_prime<uint64_t>(p2));

    // Ensure their product is bounded (b-1)^2*2^(n-1) < m.
    const uint64_t m = p1 * p2;
    // Check overflow.
    ASSERT_EQ(m / p1, p2);
    ASSERT_LT(arith::exp<uint64_t>((b - 1), 2) * arith::exp<uint64_t>(p, 2), m);

    // Find `x` so it is not a quadratic residue in GF_p1 and GF_p2.²
    ASSERT_EQ(arith::jacobi<uint64_t>(3, p1), arith::jacobi<uint64_t>(p1, 3));
    ASSERT_EQ(arith::jacobi<uint64_t>(p1, 3), arith::jacobi<uint64_t>(2, 3));
    ASSERT_EQ(arith::jacobi<uint64_t>(3, p2), arith::jacobi<uint64_t>(p2, 3));
    ASSERT_EQ(arith::jacobi<uint64_t>(p2, 3), arith::jacobi<uint64_t>(2, 3));
    ASSERT_EQ(arith::jacobi<uint64_t>(2, 3), -1);
    // Which means x=3 is not a quadratic residue in GF_p1 and GF_p2.

    // Therefore we can compute 2^n-th roots of unity in GF_p1 and GF_p2
    const uint64_t w1 = arith::exp<uint64_t>(3, a1);
    const uint64_t w2 = arith::exp<uint64_t>(3, a2);
    ASSERT_EQ(w1, 9);
    ASSERT_EQ(w2, 243);

    // find root of unity in GF_p1p2
    uint64_t _a[2] = {w1, w2};
    uint64_t _n[2] = {p1, p2};
    const uint64_t w = arith::chinese_remainder<uint64_t>(2, _a, _n);
    ASSERT_EQ(w, 25559439);
}
