/* -*- mode: c++ -*- */
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
#ifndef __QUAD_ARITH_H__
#define __QUAD_ARITH_H__

#include <cassert>
#include <cstdlib>
#include <vector>

#include "big_int.h"
#include "core.h"
#include "exceptions.h"

namespace quadiron {

/** Base/core arithmetical functions of QuadIron. */
namespace arith {

template <typename T>
T sqrt(T n);
template <typename T>
T exp(T base, T exponent);
template <typename T>
T exp_mod(T base, T exponent, T modulus);
template <typename T>
bool is_power_of_2(int x);
template <typename T>
T ceil2(int x);
template <typename T>
int log2(int x);
template <typename T>
int exp2(int x);
template <typename T>
SignedDoubleSizeVal<T> extended_gcd(
    SignedDoubleSizeVal<T> a,
    SignedDoubleSizeVal<T> b,
    SignedDoubleSizeVal<T> bezout_coef[2],
    SignedDoubleSizeVal<T> quotient_gcd[2]);
template <typename T>
T chinese_remainder(int n_mod, T a[], T n[]);
template <typename T>
int jacobi(SignedDoubleSizeVal<T> n, SignedDoubleSizeVal<T> m);
template <typename T>
bool solovay_strassen1(T a, T n);
template <typename T>
bool solovay_strassen(T n, T (*rand_func)(T n));
template <typename T>
bool is_prime(T n);
template <typename T>
T gcd(T u, T v);
template <typename T>
std::vector<T> factor_distinct_prime(T n);
template <typename T>
void factor_prime(T nb, std::vector<T>* primes, std::vector<int>* exponent);
template <typename T>
std::vector<T> get_proper_divisors(T n);
template <typename T>
std::vector<T> get_proper_divisors(T n, const std::vector<T>& primes);
template <typename T>
std::vector<T> get_all_divisors(T n);
template <typename T>
T get_code_len(T order, T n);
template <typename T>
T get_code_len_high_compo(T order, T n);
template <typename T>
T get_code_len_high_compo(const std::vector<T>& factors, T n);
template <typename T>
std::vector<T> get_coprime_factors(T nb);
template <typename T>
std::vector<T> get_prime_factors(T nb);
template <typename T>
std::vector<T> get_prime_factors(
    const std::vector<T>& primes,
    const std::vector<int>& exponent);

/** Compute the integer square root.
 *
 * @param[in] n a value
 * @return the integer square root of `n`
 *
 * @pre `n` must be strictly positive.
 */
// Taken from Wikipedia.
template <typename T>
T sqrt(T n)
{
    T place = static_cast<T>(1) << (sizeof(T) * CHAR_BIT - 2);

    while (place > n) {
        place /= 4;
    }

    T root = 0;
    while (place) {
        if (n >= root + place) {
            n -= root + place;
            root += place * 2;
        }
        root /= 2;
        place /= 4;
    }

    return root;
}

/** Compute the regular exponentation.
 *
 * @param[in] base the base `b`
 * @param[in] exponent the exponent `e`
 * @return the value of \f$b^e\f$
 */
template <typename T>
T exp(T base, T exponent)
{
    if (exponent == 0) {
        return 1;
    }
    T result = 1;
    while (true) {
        if (exponent % 2 == 1) {
            result *= base;
        }
        exponent >>= 1;
        if (exponent == 0) {
            break;
        }
        base *= base;
    }

    return result;
}

/** Compute the modular exponentation.
 *
 * @param[in] base the base `b`
 * @param[in] exponent the exponent `e`
 * @param[in] modulus the modulus `m`
 * @return the value of \f$b^e \mod m\f$
 */
// Taken from Applied Cryptography by Bruce Schneier.
template <typename T>
T exp_mod(T base, T exponent, T modulus)
{
    if (modulus == 1) {
        return 0;
    }

    T result = 1;
    base = base % modulus;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = (result * base) % modulus;
        }
        exponent = exponent >> 1;
        base = (base * base) % modulus;
    }

    return result;
}

/** Test if `n` is a power of two.
 *
 * @param[in] n a number
 * @return true if `n` is a power of two.
 */
template <typename T>
bool is_power_of_2(int n)
{
    return n > 0 && !(n & (n - 1));
}

/** Round up to the next power of two.
 *
 * @param[in] n a number
 * @return a power of two greater than or equal to `n`
 */
template <typename T>
T ceil2(int n)
{
    if (is_power_of_2<T>(n)) {
        return n;
    }
    return exp2<T>(log2<T>(n) + 1);
}

/** Compute the binary logarithm.
 *
 * @param[in] n a number
 * @return the value of \f$log_2(n)\f$
 *
 * @pre `n` must be strictly positive.
 */
template <typename T>
int log2(int n)
{
    if (n == 0) {
        throw DomainError("pole error: n is zero");
    }
    if (n == 1) {
        return 0;
    }

    int result = 0;
    while (n > 1) {
        n >>= 1;
        result++;
    }

    return result;
}

/** Compute \f$2^n\f$.
 *
 * @param[in] n a number
 * @return the value of \f$2^n\f$
 */
template <typename T>
int exp2(int n)
{
    return static_cast<T>(1) << n;
}

/** Compute the Extended Euclidean algorithm.
 *
 * Perform ax+by = gcd(a, b)
 *   If a and b are coprime then:
 *         x = inv(a) mod b
 *         y = inv(b) mod a
 *
 * @param[in] a a number
 * @param[in] b a number
 * @param[out] bezout_coef computed bezout coefficients, might be nullptr
 * @param[out] quotient_gcd computed quotient by the GCD, might be nullptr
 *
 * @todo take care of the signs of input
 *
 * @return the GCD
 */
// Taken from Wikipedia.
template <typename T>
SignedDoubleSizeVal<T> extended_gcd(
    SignedDoubleSizeVal<T> a,
    SignedDoubleSizeVal<T> b,
    SignedDoubleSizeVal<T> bezout_coef[2],
    SignedDoubleSizeVal<T> quotient_gcd[2])
{
    SignedDoubleSizeVal<T> s = 0;
    SignedDoubleSizeVal<T> old_s = 1;
    SignedDoubleSizeVal<T> t = 1;
    SignedDoubleSizeVal<T> old_t = 0;
    SignedDoubleSizeVal<T> r = b;
    SignedDoubleSizeVal<T> old_r = a;
    SignedDoubleSizeVal<T> quotient;
    SignedDoubleSizeVal<T> tmp;

    while (0 != r) {
        quotient = old_r / r;

        tmp = r;
        r = old_r - quotient * r;
        old_r = tmp;

        tmp = s;
        s = old_s - quotient * s;
        old_s = tmp;

        tmp = t;
        t = old_t - quotient * t;
        old_t = tmp;
    }

    if (bezout_coef) {
        bezout_coef[0] = old_s;
        bezout_coef[1] = old_t;
    }

    if (quotient_gcd) {
        quotient_gcd[0] = t;
        quotient_gcd[1] = s;
    }

    return old_r;
}

/** Apply the Chinese remainder theorem.
 *
 * @param[in] n_mod number of moduli
 * @param[in] a the a's (integers)
 * @param[in] n the n's (moduli)
 *
 * @todo Check if there is a solution.
 *
 * @return the solution.
 */
// Taken from Wikipedia.
template <typename T>
T chinese_remainder(int n_mod, T a[], T n[])
{
    SignedDoubleSizeVal<T> acc;
    SignedDoubleSizeVal<T> x;
    SignedDoubleSizeVal<T>* N = new SignedDoubleSizeVal<T>[n_mod];
    SignedDoubleSizeVal<T>* M = new SignedDoubleSizeVal<T>[n_mod];

    acc = 1;
    for (int i = 0; i < n_mod; i++) {
        acc *= n[i];
    }

    for (int i = 0; i < n_mod; i++) {
        SignedDoubleSizeVal<T> bezout[2];

        N[i] = acc / n[i];
        extended_gcd<T>(N[i], n[i], bezout, nullptr);
        M[i] = bezout[0]; // XXX % n[i];
    }

    x = 0;
    for (int i = 0; i < n_mod; i++) {
        x += N[i] * M[i] * a[i];
    }

    delete[] N;
    delete[] M;

    x = x % acc;

    if (x < 0) {
        x = acc + x;
    }

    return x;
}

/** Compute the Jacobi symbol of the numbers `n` and `m`.
 *
 * The Jacobi symbol is defined as follow:
 *    n    |  0 if n % m == 0
 *  ( _) = |  1 if n % m != 0 and for some x, n % m == x^2 (quadratic residue)
 *    m    | -1 if n % m != 0 and there is no such x
 *
 * @param[in] n a value
 * @param[in] m a modulus
 *
 * @return the Jacobi symbol
 * @throw DomainError if b is even or b is negative
 */
// Taken from https://groups.google.com/forum/#!topic/sci.crypt/v9_cNF06XjU
template <typename T>
int jacobi(SignedDoubleSizeVal<T> n, SignedDoubleSizeVal<T> m)
{
    SignedDoubleSizeVal<T> t;
    int jac;

    if ((m % 2) == 0) {
        throw DomainError("m must be odd");
    }
    if (m < 0) {
        throw DomainError("m must be positive");
    }

    jac = 1;
    while (m != 1) {
        // if the gcd of $n$ and $m$ is $>1$ Jacobi returns $0$
        if ((n == 0) || ((n % 2 == 0) && (m % 2 == 0))) {
            jac = 0;
            m = 1;
        }
        // $J(n,2*m) = J(n,m) * J(n,2) = J(n,m) * (-1)^{(n^2-1)/8}$
        else if (m % 2 == 0) {
            if ((n % 8 == 3) || (n % 8 == 5))
                jac = -jac;
            m = m / 2;
        }
        // $J(2*n,m) = J(n,m) * J(2,m) = J(n,m) * (-1)^{(m^2-1)/8}$
        else if (n % 2 == 0) {
            if ((m % 8 == 3) || (m % 8 == 5))
                jac = -jac;
            n = n / 2;
        }
        // $J(-n,m) = J(n,m) * J(-1,m) = J(n,m) * (-1)^{(m-1)/2}$
        else if (n < 0) {
            if (m % 4 == 3)
                jac = -jac;
            n = -n;
        }
        // $J(n,m) = J(m,n) * (-1)^{(n-1)*(m-1)/4}$ (quadratic reciprocity)
        else {
            if ((n % 4 == 3) && (m % 4 == 3))
                jac = -jac;
            t = n;
            n = m % n;
            m = t;
        }
    }

    return jac;
}

template <typename T>
bool solovay_strassen1(T a, T n)
{
    const T exponent = (n - 1) / 2;
    const T res = exp_mod(a, exponent, n);
    int j = jacobi<T>(a, n);
    T _res;
    if (j < 0) {
        _res = res - n;
    } else {
        _res = res;
    }
    return j == _res;
}

/** Perform the Solovay–Strassen primality test.
 *
 * @param n check if `n` is prime
 * @param rand_func return a random number between 0 and n-1
 * @return true if `n` is (probably) prime else false.
 */
template <typename T>
bool solovay_strassen(T n, T (*rand_func)(T n))
{
    for (int i = 0; i < 100; i++) {
        T a = rand_func(n);
        if (!solovay_strassen1(a, n))
            return false;
    }
    return true;
}

/** Brute force prime checking until sqrt(n)
 *
 * @param[in] n a number
 * @return true if `n` is a prime number, otherwise false
 */
template <typename T>
bool is_prime(T n)
{
    if (n == 2)
        return true;

    T root = sqrt<T>(n);
    for (T i = 2; i <= root; i++) {
        if (n % i == 0)
            return false;
    }
    return true;
}

/** Factor a given number into distinct primes.
 *
 * @param nb the number to factorize.
 * @return the list of prime factor
 *
 * @note Only primes are stored, their exponent is ignored.
 */
template <typename T>
std::vector<T> factor_distinct_prime(T nb)
{
    std::vector<T> primes;
    T last_found = 1;

    while (nb % 2 == 0) {
        if (last_found != 2) {
            primes.push_back(2);
            last_found = 2;
        }
        nb = nb / 2;
    }
    // `nb` must be odd at this point => we can skip one element.
    assert((nb & 1) == 1);
    for (T i = 3; i <= sqrt<T>(nb); i = i + 2) {
        // While `i` divides `nb`, get `i` and divide `nb`.
        while (nb % i == 0) {
            if (last_found != i) {
                primes.push_back(i);
                last_found = i;
            }
            nb = nb / i;
        }
    }
    // Handle the case when `nb` is a prime number greater than 2.
    if (nb > 2) {
        primes.push_back(nb);
    }

    return primes;
}

/// Get the proper divisors of a number.
template <typename T>
std::vector<T> get_proper_divisors(T nb)
{
    const std::vector<T> primes = factor_distinct_prime<T>(nb);
    return get_proper_divisors(nb, primes);
}

/// Get proper divisors of a number from its factored distinct primes.
template <typename T>
std::vector<T> get_proper_divisors(T nb, const std::vector<T>& primes)
{
    std::vector<T> divisors;

    for (auto prime : primes) {
        if (prime < nb) {
            divisors.push_back(nb / prime);
        }
    }

    return divisors;
}

/// Factor a given number into primes.
template <typename T>
void factor_prime(T nb, std::vector<T>* primes, std::vector<int>* exponent)
{
    int occurence = 0;

    assert(primes != nullptr);
    assert(exponent != nullptr);

    while (nb % 2 == 0) {
        occurence++;
        if (occurence == 1) {
            primes->push_back(2);
        }
        nb = nb / 2;
    }

    if (occurence > 0) {
        exponent->push_back(occurence);
        occurence = 0;
    }
    // `nb` must be odd at this point => we can skip one element.
    assert((nb & 1) == 1);
    for (T i = 3; i <= sqrt<T>(nb); i = i + 2) {
        // While `i` divides `nb`, get `i` and divide `nb`.
        while (nb % i == 0) {
            occurence++;
            if (occurence == 1) {
                primes->push_back(i);
            }
            nb = nb / i;
        }
        if (occurence > 0) {
            exponent->push_back(occurence);
            occurence = 0;
        }
    }
    // Handle the case when `nb` is a prime number greater than 2.
    if (nb > 2) {
        primes->push_back(nb);
        exponent->push_back(1);
    }
}

/** Modern Euclidean algorithm implementation
 *
 * This is an implementation of the Algorithm A, p.322 of "The Art of Computer
 * Programming" by Donald Knuth
 *
 * @param[in] u a number
 * @param[in] v a number
 * @return the GCD(u, v)
 */
template <typename T>
T gcd(T u, T v)
{
    T r;
    if (v == 0)
        return u;
    if (u == 0)
        return v;
    while (1) {
        r = u % v;
        if (r == 0)
            return v;
        u = v;
        v = r;
    }
}

/// Compute all divisors of a number.
template <typename T>
std::vector<T> get_all_divisors(T n)
{
    std::vector<T> divisors;
    std::vector<T> n_div_i;

    for (T i = 1; i <= sqrt<T>(n); i++) {
        // While `i` divides `nb`, get `i` and divide `nb/i`.
        if (n % i == 0) {
            divisors.push_back(i);
            n_div_i.push_back(n / i);
        }
    }
    divisors.insert(divisors.end(), n_div_i.rbegin(), n_div_i.rend());

    return divisors;
}

/** Find smallest number that is ≥ `n` and divisible by `order`.
 *
 * @param[in] order a prime number
 * @param[in] n the lower bound for the result
 * @return the code length
 */
template <typename T>
T get_code_len(T order, T n)
{
    if (order % n == 0) {
        return n;
    }
    if (order < n) {
        assert(false);
    }
    T nb_sqrt = sqrt<T>(order);
    T i;
    if (n > nb_sqrt) {
        for (i = order / n; i >= 1; i--) {
            // if i divides n, return n/i
            if (order % i == 0) {
                return order / i;
            }
        }
        assert(false);
    }
    for (i = n; i <= nb_sqrt; i++) {
        // if i divides n, return i
        if (order % i == 0) {
            return i;
        }
    }
    // order is prime
    return order;
}

/** Find smallest number that is highly composite, ≥ `n` and divisible by
 * `order`.
 *
 * @param[in] order a prime number
 * @param[in] n the lower bound for the result
 * @return the code length
 */
template <typename T>
T get_code_len_high_compo(T order, T n)
{
    assert(order >= n);
    return get_code_len_high_compo(get_prime_factors<T>(order), n);
}

/** Find smallest number that is highly composite and ≥ `n`
 *
 * @param[in] factors vector of all prime factors of a given order
 * @param[in] n the lower bound for the result
 * @return the code length
 */
template <typename T>
T get_code_len_high_compo(const std::vector<T>& factors, T n)
{
    T x = 1;
    typename std::vector<T>::size_type i, j;
    // forward to get a divisor of (q-1) >= n and of highly composited
    for (i = 0; i != factors.size(); ++i) {
        x *= factors[i];
        if (x >= n) {
            // backward to get smaller number
            for (j = 0; j != i + 1 && j != factors.size(); j++) {
                x /= factors[j];
                if (x < n) {
                    return x * factors[j];
                }
            }
        }
    }
    assert(false);
    return 0;
}

/** Get all the coprime factors of a number.
 *
 * @param nb the number to be factored
 * @return the co-prime divisors of `nb`
 */
template <typename T>
std::vector<T> get_coprime_factors(T nb)
{
    std::vector<T> coprimes;
    std::vector<T> primes;
    std::vector<int> exponents;

    factor_prime<T>(nb, &primes, &exponents);
    assert(primes.size() == exponents.size());

    coprimes.reserve(primes.size());
    for (size_t i = 0; i != primes.size(); ++i) {
        coprimes.push_back(exp<T>(primes[i], exponents[i]));
    }

    return coprimes;
}

/** Get all the prime factors of a number.
 *
 * @param nb the number to be factored
 * @return the prime factors of `nb`.
 */
template <typename T>
std::vector<T> get_prime_factors(T nb)
{
    std::vector<T> primes;
    std::vector<int> exponent;

    factor_prime<T>(nb, &primes, &exponent);
    return get_prime_factors<T>(primes, exponent);
}

/** Get all prime factors of a number from pre-computed primes and exponents.
 *
 * @param primes primes
 * @param exponents exponents for each prime
 * @return the prime factors of `nb`.
 *
 * @pre `primes` and `exponents` must have the same size.
 */
template <typename T>
std::vector<T> get_prime_factors(
    const std::vector<T>& primes,
    const std::vector<int>& exponents)
{
    std::vector<T> factors;
    assert(primes.size() == exponents.size());

    for (size_t i = 0; i != primes.size(); ++i) {
        for (int j = 0; j != exponents[i]; ++j) {
            factors.push_back(primes[i]);
        }
    }
    return factors;
}

} // namespace arith
} // namespace quadiron

#endif
