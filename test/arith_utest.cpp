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

#include "nttec.h"

template <typename T>
class ArithUtest {
  public:
    T max;
    std::uniform_int_distribution<uint32_t> uniform_dist_max;

    ArithUtest()
    {
        if (sizeof(T) < 16)
            max = std::numeric_limits<T>::max();
        else {
            max = 0;
            for (size_t i = 0; i < sizeof(T) * 8; i++)
                max += (T)1 << i;
        }
        uniform_dist_max = std::uniform_int_distribution<uint32_t>(
            1, static_cast<uint32_t>(max - 1));
    }

    ~ArithUtest() {}

    void test_basic_ops()
    {
        std::cout << "test_basic_ops\n";
        assert(nttec::arith::sqrt<T>(2025) == 45);
        assert(nttec::arith::is_prime<T>(2));
        assert(nttec::arith::is_prime<T>(3));
        assert(nttec::arith::is_prime<T>(13));
        assert(!nttec::arith::is_prime<T>(4));
        assert(!nttec::arith::is_prime<T>(15));
    }

    void test_reciprocal()
    {
        std::cout << "test_reciprocal\n";

        int i;
        T sub_max = nttec::arith::sqrt<T>(max);
        std::uniform_int_distribution<uint32_t> dis(1, sub_max - 1);

        for (i = 0; i < 1000; i++) {
            T x, y, z;

            x = dis(nttec::prng());
            // std::cout << "x=" << x << "\n";
            y = nttec::arith::exp<T>(x, 2);
            // std::cout << "exp(x)=" << y << "\n";
            z = nttec::arith::sqrt<T>(y);
            // std::cout << "z=" << z << "\n";
            assert(z == x);
        }
    }

    /**
     * http://www.math.unm.edu/~loring/links/discrete_f05/remainder.pdf
     * http://gauss.math.luc.edu/greicius/Math201/Fall2012/Lectures/ChineseRemainderThm.article.pdf
     *
     * @param gf
     */
    void test_chinese_remainder()
    {
        std::cout << "test_chinese_remainder\n";

        T a[4];
        T n[4];
        T omega;

        a[0] = 4;
        n[0] = 107;
        a[1] = 2;
        n[1] = 74;
        omega = nttec::arith::chinese_remainder<T>(2, a, n);
        assert(omega == 5996);

        a[0] = 6;
        n[0] = 7;
        a[1] = 4;
        n[1] = 8;
        omega = nttec::arith::chinese_remainder<T>(2, a, n);
        assert(omega == 20);

        a[0] = 3;
        n[0] = 4;
        a[1] = 0;
        n[1] = 6;
        omega = nttec::arith::chinese_remainder<T>(2, a, n);
        // no solution XXX detect it
    }

    void test_jacobi()
    {
        std::cout << "test_jacobi\n";
        assert(nttec::arith::jacobi<T>(1001, 9907) == -1);
        assert(nttec::arith::jacobi<T>(19, 45) == 1);
        assert(nttec::arith::jacobi<T>(8, 21) == -1);
        assert(nttec::arith::jacobi<T>(5, 21) == 1);
        assert(nttec::arith::jacobi<T>(47, 221) == -1);
        assert(nttec::arith::jacobi<T>(2, 221) == -1);
    }

    /**
     * Schönhage-Strassen algorithm
     * Example taken from Pierre Meunier's book
     *
     * @param gf
     */
    void test_mul_bignum()
    {
        // current test is only for uint64_t
        if (sizeof(T) != 8)
            return;
        std::cout << "test_mul_bignum\n";

        int b = 10; // base
        int p = 14; // we could multiply integers of 2^p digits
        // T max_digits = nttec::arith::exp<T>(2, p);
        // std::cerr << "p=" << p << " max_digits=" << max_digits << "\n";

        // T l = p + 1;
        // std::cerr << "l=" << l << "\n";

        // choose 2 prime numbers of the form p=a.2^n+1
        // because if x is not a quadratic residue then w=x^a is
        // a 2^n-th principal root of unity in GF_p
        T a1 = 2;
        T a2 = 5;
        T p1 = a1 * nttec::arith::exp<T>(2, 15) + 1;
        T p2 = a2 * nttec::arith::exp<T>(2, 15) + 1;
        // std::cerr << "p1=" << p1 << " p2=" << p2 << "\n";
        assert(nttec::arith::is_prime<T>(p1));
        assert(nttec::arith::is_prime<T>(p2));

        // ensure their product is bounded (b-1)^2*2^(n-1) < m
        T m = p1 * p2;
        // check overflow
        assert(m / p1 == p2);
        // std::cerr << " m=" << m << "\n";
        assert(
            nttec::arith::exp<T>((b - 1), 2) * nttec::arith::exp<T>(p, 2) < m);

        // find x so it is not a quadratic residue in GF_p1 and GF_p2
        assert(
            nttec::arith::jacobi<T>(3, p1) == nttec::arith::jacobi<T>(p1, 3));
        assert(nttec::arith::jacobi<T>(p1, 3) == nttec::arith::jacobi<T>(2, 3));
        assert(
            nttec::arith::jacobi<T>(3, p2) == nttec::arith::jacobi<T>(p2, 3));
        assert(nttec::arith::jacobi<T>(p2, 3) == nttec::arith::jacobi<T>(2, 3));
        assert(nttec::arith::jacobi<T>(2, 3) == -1);
        // which means x=3 is not a quadratic residue in GF_p1 and GF_p2

        // therefore we can compute 2^n-th roots of unity in GF_p1 and GF_p2
        T w1 = nttec::arith::exp<T>(3, a1);
        T w2 = nttec::arith::exp<T>(3, a2);
        // std::cerr << "w1=" << w1 << " w2=" << w2 << "\n";
        assert(w1 == 9);
        assert(w2 == 243);

        // find root of unity in GF_p1p2
        T _a[2];
        T _n[2];
        _a[0] = w1;
        _n[0] = p1;
        _a[1] = w2;
        _n[1] = p2;
        T w = nttec::arith::chinese_remainder<T>(2, _a, _n);
        // std::cerr << " w=" << w << "\n";
        assert(w == 25559439);
    }

    void test_ext_gcd()
    {
        std::cout << "test_ext_gcd\n";
        nttec::SignedDoubleSizeVal<T> bezout[2];

        // not explicitely related to GF(97)
        assert(2 == nttec::arith::extended_gcd<T>(240, 46, nullptr, nullptr));
        assert(6 == nttec::arith::extended_gcd<T>(54, 24, nullptr, nullptr));
        assert(15 == nttec::arith::extended_gcd<T>(210, 45, nullptr, nullptr));
        //
        assert(1 == nttec::arith::extended_gcd<T>(97, 20, bezout, nullptr));
        assert(bezout[0] == -7 && bezout[1] == 34);
    }

    void check_all_primes(std::vector<T>* primes, bool distinct)
    {
        assert(primes->size() > 0);

        typename std::vector<T>::size_type j;

        assert(nttec::arith::is_prime<T>(primes->at(0)));

        for (j = 1; j != primes->size(); ++j) {
            // std::cout << j << ": " << primes->at(j) << "\n";
            assert(nttec::arith::is_prime<T>(primes->at(j)));
            if (distinct)
                assert(primes->at(j - 1) != primes->at(j));
        }
    }

    void check_primes_exponent(
        T nb,
        std::vector<T>* primes,
        std::vector<int>* exponent)
    {
        assert(primes->size() > 0);
        assert(exponent->size() > 0);

        typename std::vector<int>::size_type j;
        T y = 1;
        for (j = 0; j != primes->size(); ++j) {
            y *= nttec::arith::exp<T>(primes->at(j), exponent->at(j));
        }
        assert(y == nb);
    }

    void test_factor_distinct_prime()
    {
        std::cout << "test_factor_distinct_prime\n";

        int i;

        for (i = 0; i < 1000; i++) {
            const T x = uniform_dist_max(nttec::prng());
            std::vector<T> primes;

            nttec::arith::factor_distinct_prime<T>(x, &primes);
            check_all_primes(&primes, true);
        }
    }

    void test_factor_prime()
    {
        std::cout << "test_factor_prime\n";

        int i;
        std::vector<T> primes;
        std::vector<int> exponent;

        for (i = 0; i < 1000; i++) {
            const T x = uniform_dist_max(nttec::prng());

            nttec::arith::factor_prime<T>(x, &primes, &exponent);
            check_all_primes(&primes, true);
            check_primes_exponent(x, &primes, &exponent);
            primes.clear();
            exponent.clear();
        }
    }

    void check_divisors(T nb, std::vector<T>* divisors, bool proper)
    {
        if (proper && nttec::arith::is_prime<T>(nb))
            assert(divisors->size() == 0);
        else
            assert(divisors->size() > 0);

        typename std::vector<T>::size_type i;
        for (i = 0; i != divisors->size(); ++i) {
            // std::cout << i << ": " << divisors->at(i) << "\n";
            assert(nb % divisors->at(i) == 0);
            if (proper)
                assert(nttec::arith::is_prime<T>(nb / divisors->at(i)));
        }
    }

    void test_get_proper_divisors()
    {
        std::cout << "test_get_proper_divisors\n";

        int i;
        std::vector<T> divisors;
        for (i = 0; i < 1000; i++) {
            const T x = uniform_dist_max(nttec::prng());

            nttec::arith::get_proper_divisors<T>(x, &divisors);
            check_divisors(x, &divisors, true);
            divisors.clear();
        }
    }

    void test_get_proper_divisors_2()
    {
        std::cout << "test_get_proper_divisors_2\n";

        int i;
        std::vector<T> divisors;
        for (i = 0; i < 1000; i++) {
            const T x = uniform_dist_max(nttec::prng());

            std::vector<T> factors;
            nttec::arith::factor_distinct_prime<T>(x, &factors);
            nttec::arith::get_proper_divisors<T>(x, &factors, &divisors);
            check_divisors(x, &divisors, true);
            divisors.clear();
        }
    }

    void test_get_all_divisors()
    {
        std::cout << "test_get_all_divisors\n";

        int i;
        std::vector<T> divisors;
        for (i = 0; i < 1000; i++) {
            const T x = uniform_dist_max(nttec::prng());

            nttec::arith::get_all_divisors<T>(x, &divisors);
            check_divisors(x, &divisors, false);
            divisors.clear();
        }
    }

    void test_get_code_len()
    {
        std::cout << "test_get_code_len\n";

        int i;
        for (i = 0; i < 1000; i++) {
            const T order = uniform_dist_max(nttec::prng());
            std::uniform_int_distribution<uint32_t> dis(1, order - 1);
            const T n = dis(nttec::prng());
            const T len = nttec::arith::get_code_len<T>(order, n);

            assert(order % len == 0);
            assert(len >= n);
        }
    }

    void test_get_code_len_high_compo()
    {
        std::cout << "test_get_code_len_high_compo\n";

        int i;
        for (i = 0; i < 1000; i++) {
            const T order = uniform_dist_max(nttec::prng());
            std::uniform_int_distribution<uint32_t> dis(1, order - 1);
            const T n = dis(nttec::prng());
            const T len = nttec::arith::get_code_len_high_compo<T>(order, n);

            assert(order % len == 0);
            assert(len >= n);
        }
    }

    void test_get_code_len_high_compo_2()
    {
        std::cout << "test_get_code_len_high_compo_2\n";

        int i;
        for (i = 0; i < 1000; i++) {
            const T order = uniform_dist_max(nttec::prng());
            std::uniform_int_distribution<uint32_t> dis(1, order - 1);
            const T n = dis(nttec::prng());

            std::vector<T> factors;
            nttec::arith::get_prime_factors<T>(order, &factors);
            const T len = nttec::arith::get_code_len_high_compo<T>(&factors, n);
            // std::cout << "len=" << len << "\n";
            assert(order % len == 0);
            assert(len >= n);
        }
    }

    void check_prime_divisors(T nb, std::vector<T>* divisors, bool coprime)
    {
        assert(divisors->size() > 0);

        typename std::vector<T>::size_type i, j;
        for (i = 0; i != divisors->size(); ++i) {
            // std::cout << i << ": " << divisors->at(i) << "\n";
            assert(nb % divisors->at(i) == 0);
            if (!coprime)
                assert(nttec::arith::is_prime<T>(divisors->at(i)));
            else {
                for (j = i + 1; j != divisors->size(); ++j) {
                    assert(divisors->at(i) != divisors->at(j));
                }
            }
        }
    }

    void test_get_coprime_factors()
    {
        std::cout << "test_get_coprime_factors\n";

        std::vector<T> divisors;
        int i;
        for (i = 0; i < 1000; i++) {
            const T n = uniform_dist_max(nttec::prng());

            nttec::arith::get_coprime_factors<T>(n, &divisors);
            check_prime_divisors(n, &divisors, true);
            divisors.clear();
        }
    }

    void test_get_prime_factors()
    {
        std::cout << "test_get_prime_factors\n";

        std::vector<T> divisors;
        int i;
        for (i = 0; i < 1000; i++) {
            const T n = uniform_dist_max(nttec::prng());

            nttec::arith::get_prime_factors<T>(n, &divisors);
            check_prime_divisors(n, &divisors, false);
            divisors.clear();
        }
    }

    void arith_utest()
    {
        std::cout << "arith_utest with sizeof(T)=" << sizeof(T) << "\n";

        nttec::prng().seed(time(0));

        test_basic_ops();
        test_reciprocal();
        test_chinese_remainder();
        test_jacobi();
        test_mul_bignum();
        test_ext_gcd();
        test_factor_distinct_prime();
        test_factor_prime();
        test_get_proper_divisors();
        test_get_all_divisors();
        test_get_code_len();
        test_get_code_len_high_compo();
        test_get_coprime_factors();
        test_get_prime_factors();
    }

    void arith_utest_no256()
    {
        std::cout << "arith_utest_no256\n";

        nttec::prng().seed(time(0));

        test_basic_ops();
        test_reciprocal();
        test_jacobi();
        test_factor_distinct_prime();
        test_factor_prime();
        test_get_proper_divisors();
        test_get_all_divisors();
        test_get_code_len();
        test_get_code_len_high_compo();
        test_get_coprime_factors();
        test_get_prime_factors();
    }
};

void arith_utest()
{
    ArithUtest<uint32_t> ArithUtest_uint32;
    ArithUtest_uint32.arith_utest();
    ArithUtest<uint64_t> ArithUtest_uint64;
    ArithUtest_uint64.arith_utest();
    ArithUtest<__uint128_t> ArithUtest_uint128;
    ArithUtest_uint128.arith_utest_no256();
}
