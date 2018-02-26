/* -*- mode: c++ -*- */
/*
 * Copyright 2017-2018 the NTTEC authors
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
#ifndef __NTTEC_ARITH_H__
#define __NTTEC_ARITH_H__

#include <cassert>
#include <cstdlib>
#include <vector>

#include "big_int.h"
#include "core.h"
#include "exceptions.h"

namespace nttec {

template <typename T>
using DoubleSizeVal = typename DoubleSize<T>::T;

template <typename T>
using SignedDoubleSizeVal = typename SignedDoubleSize<T>::T;

/** Base/core arithmetical functions of NTTEC. */
namespace arith {

template <class T>
T sqrt(T n);
template <class T>
T exp(T base, T exponent);
template <class T>
T exp_mod(T base, T exponent, T modulus);
template <class T>
bool is_power_of_2(int x);
template <class T>
T get_smallest_power_of_2(int x);
template <class T>
int log2(int x);
template <class T>
int exp2(int x);
template <class T>
SignedDoubleSizeVal<T> extended_gcd(
    SignedDoubleSizeVal<T> a,
    SignedDoubleSizeVal<T> b,
    SignedDoubleSizeVal<T> bezout_coef[2],
    SignedDoubleSizeVal<T> quotient_gcd[2]);
template <class T>
T chinese_remainder(int n_mod, T a[], T n[]);
template <class T>
int jacobi(SignedDoubleSizeVal<T> n, SignedDoubleSizeVal<T> m);
template <class T>
bool solovay_strassen1(T a, T n);
template <class T>
bool solovay_strassen(T n);
template <class T>
bool is_prime(T n);
template <class T>
T weak_rand(T max);
template <class T>
T weak_rand0(T max);
template <class T>
T gcd(T u, T v);
template <class T>
void factor_distinct_prime(T n, std::vector<T>* output);
template <class T>
void factor_prime(T nb, std::vector<T>* primes, std::vector<int>* exponent);
template <class T>
void get_proper_divisors(T n, std::vector<T>* output);
template <class T>
void get_proper_divisors(T n, std::vector<T>* primes, std::vector<T>* output);
template <class T>
void get_all_divisors(T n, std::vector<T>* output);
template <class T>
T get_code_len(T order, T n);
template <class T>
T get_code_len_high_compo(T order, T n);
template <class T>
T get_code_len_high_compo(std::vector<T>* factors, T n);
template <class T>
void get_coprime_factors(T nb, std::vector<T>* output);
template <class T>
void get_prime_factors(T nb, std::vector<T>* output);
template <class T>
void get_prime_factors_final(
    std::vector<T>* primes,
    std::vector<int>* exponent,
    std::vector<T>* output);

/**
 * integer square root (from Wikipedia)
 *
 * @param n
 *
 * @return
 */
template <class T>
T sqrt(T remainder)
{
    // calculated by precompiler = same runtime as: place = 0x40000000
    T place = (T)1 << (sizeof(T) * 8 - 2);
    while (place > remainder)
        place /= 4; // optimized by complier as place >>= 2

    T root = 0;
    while (place) {
        if (remainder >= root + place) {
            remainder -= root + place;
            root += place * 2;
        }
        root /= 2;
        place /= 4;
    }

    return root;
}

/**
 * Regular exponentation
 *
 * @param gf
 * @param base
 * @param exponent
 * @param modulus
 *
 * @return
 */
template <class T>
T exp(T base, T exponent)
{
    if (exponent == 0)
        return 1;
    T result = 1;
    while (true) {
        if (exponent % 2 == 1)
            result *= base;
        exponent >>= 1;
        if (exponent == 0)
            break;
        base *= base;
    }

    return result;
}

/**
 * Modular exponentation taken from Applied Cryptography by Bruce Schneier.
 *
 * @param gf
 * @param base
 * @param exponent
 * @param modulus
 *
 * @return
 */
template <class T>
T exp_mod(T base, T exponent, T modulus)
{
    if (1 == modulus)
        return 0;

    T result = 1;
    base = base % modulus;
    while (exponent > 0) {
        if (exponent % 2 == 1)
            result = (result * base) % modulus;
        exponent = exponent >> 1;
        base = (base * base) % modulus;
    }

    return result;
}

/**
 * check if x is a power of 2
 *
 * @param x
 *
 * @return
 */
template <class T>
bool is_power_of_2(int x)
{
    return x > 0 && !(x & (x - 1));
}

/**
 * Get a smallest power of 2 and >= x
 *
 * @param x
 *
 * @return
 */
template <class T>
T get_smallest_power_of_2(int x)
{
    if (is_power_of_2<T>(x))
        return x;
    return exp2<T>(log2<T>(x) + 1);
}

/**
 * Compute log2(x)
 *
 * @param exponent
 *
 * @return
 */
template <class T>
int log2(int x)
{
    int result = 0;

    if (x == 0) {
        throw DomainError("pole error: x is zero");
    }

    if (x == 1)
        return 0;

    while (x > 1) {
        x >>= 1;
        result++;
    }

    return result;
}

/**
 * Compute 2^x
 *
 * @param x
 *
 * @return
 */
template <class T>
int exp2(int x)
{
    return 1 << x;
}

/**
 * Extended Euclidean algorithm (from Wikipedia)
 * Perform ax+by = gcd(a, b)
 *   If a and b are coprime then:
 *         x = inv(a) mod b
 *         y = inv(b) mod a
 *
 * @param a
 * @param b
 * @param bezout_coef[output] computed bezout coefficients, might be nullptr
 * @param quotient_gcd[output] computed quotient by the GCD, might be nullptr
 *
 * XXX take care of the signs of input
 *
 * @return the GCD
 */
template <class T>
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
 * Implementation from Wikipedia.
 *
 * @param n_mod number of moduli
 * @param a the a's (integers)
 * @param n the n's (moduli)
 *
 * @todo Check if there is a solution.
 *
 * @return the solution.
 */
template <class T>
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

/**
 * Compute the jacobi symbol of the 2 numbers
 * https://groups.google.com/forum/#!topic/sci.crypt/v9_cNF06XjU
 *
 * The jacobi symbol is defined as follow:
 *    n    |  0 if n % m == 0
 *  ( _) = |  1 if n % m != 0 and for some x, n % m == x^2 (quadratic residue)
 *    m    | -1 if n % m != 0 and there is no such x
 *
 * @param n
 * @param m
 *
 * @return the jacobi symbol
 * @throw DomainError if b is even or b is negative
 */
template <class T>
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

template <class T>
bool solovay_strassen1(T a, T n)
{
    const T exponent = (n - 1) / 2;
    const T res = exp_mod(a, exponent, n);
    int j = jacobi(a, n);
    if (j < 0) {
        res = res - n;
    }
    return j == res;
}

/** Perform the Solovay–Strassen primality test.
 *
 * @param n check if `n` is prime
 * @return true if `n` is (probably) prime else false.
 */
template <class T>
bool solovay_strassen(T n)
{
    for (int i = 0; i < 100; i++) {
        T a = weak_rand<T>(n);
        if (!solovay_strassen1(a, n))
            return false;
    }
    return true;
}

/**
 * Brute force prime checking until sqrt(n)
 *
 * @param n
 *
 * @return
 */
template <class T>
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

/**
 * Returns a number n such as 0 < n < max
 *
 * @param max
 *
 * @return
 */
template <class T>
T weak_rand(T max)
{
    T r;
retry:
    r = rand() % max;
    if (0 == r)
        goto retry;
    return r;
}

/**
 * Returns a number n such as 0 <= n < max
 *
 * @param max
 *
 * @return
 */
template <class T>
T weak_rand0(T max)
{
    T r;
    r = rand() % max;
    return r;
}

/*
 * A given number `n` is factored into primes-> Only primes are stored, their
 *  exponent is ignored.
 */
template <class T>
void factor_distinct_prime(T nb, std::vector<T>* output)
{
    T last_found = 1;
    while (nb % 2 == 0) {
        if (last_found != 2) {
            output->push_back(2);
            last_found = 2;
        }
        nb = nb / 2;
    }
    // n must be odd at this point.  So we can skip one element
    for (T i = 3; i <= sqrt<T>(nb); i = i + 2) {
        // While i divides n, get i and divide n
        while (nb % i == 0) {
            if (last_found != i) {
                output->push_back(i);
                last_found = i;
            }
            nb = nb / i;
        }
    }
    // This condition is to handle the case when n
    // is a prime number greater than 2
    if (nb > 2) {
        output->push_back(nb);
    }
}

/*
 * Get proper divisors of a given number from its factored distince primes
 */
template <class T>
void get_proper_divisors(T nb, std::vector<T>* output)
{
    std::vector<T> input;
    typename std::vector<T>::iterator it;
    factor_distinct_prime<T>(nb, &input);
    // std::cout << "nb: " << nb << std::endl;
    for (it = input.begin(); it != input.end(); ++it) {
        if (*it < nb) {
            output->push_back(nb / (*it));
            // std::cout << *it << std::endl;
        }
    }
}

/*
 * Get proper divisors of a given number from its factored distince primes
 */
template <class T>
void get_proper_divisors(T nb, std::vector<T>* primes, std::vector<T>* output)
{
    assert(primes != nullptr);
    typename std::vector<T>::iterator it;
    // std::cout << "nb: " << nb << std::endl;
    for (it = primes->begin(); it != primes->end(); ++it) {
        if (*it < nb) {
            output->push_back(nb / (*it));
            // std::cout << *it << std::endl;
        }
    }
}

/*
 * A given number `n` is factored into primes
 */
template <class T>
void factor_prime(T nb, std::vector<T>* primes, std::vector<int>* exponent)
{
    // std::cout << nb << ": ";
    int occurence = 0;
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
    // n must be odd at this point.  So we can skip one element
    for (T i = 3; i <= sqrt<T>(nb); i = i + 2) {
        // While i divides n, get i and divide n
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
    // This condition is to handle the case when n
    // is a prime number greater than 2
    if (nb > 2) {
        primes->push_back(nb);
        exponent->push_back(1);
    }
    // typename std::vector<int>::size_type i;
    // for (i = 0; i != primes->size(); i++) {
    //   std::cout << primes->at(i) << "^" << exponent->at(i) << " ";
    // }
    // std::cout << std::endl;
}

/**
 * Modern Euclidean algorithm
 * Implementation of Algorithm A in "The Art of Computing", p.322, Donald Knuth
 * @param u
 * @param v
 * @return the GCD(u, v)
 */
template <class T>
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

/**
 * compute all divisors of a number
 *
 */
template <class T>
void get_all_divisors(T n, std::vector<T>* output)
{
    typename std::vector<T> tmp;

    T nb_sqrt = sqrt<T>(n);
    for (T i = 1; i <= nb_sqrt; i++) {
        // While i divides n, get i and n/i
        if (n % i == 0) {
            // std::cout << n << ":" << i << " " << n/i << std::endl;
            output->push_back(i);
            tmp.push_back(n / i);
        }
    }
    output->insert(output->end(), tmp.rbegin(), tmp.rend());
}

/**
 * find smallest number is
 *  - at least n
 *  - divisible by (order)
 *
 * @param order - a prime number
 * @param n
 *
 * @return root
 */
template <class T>
T get_code_len(T order, T n)
{
    if (order % n == 0)
        return n;
    if (order < n)
        assert(false);
    T nb_sqrt = sqrt<T>(order);
    T i;
    if (n > nb_sqrt) {
        for (i = order / n; i >= 1; i--) {
            // if i divides n, return n/i
            if (order % i == 0)
                return order / i;
        }
        assert(false);
    }
    for (i = n; i <= nb_sqrt; i++) {
        // if i divides n, return i
        if (order % i == 0)
            return i;
    }
    // order is prime
    return order;
}

/**
 * find smallest number is
 *  - highly composited
 *  - at least n
 *  - divisible by (order)
 *
 * @param order - a prime number
 * @param n
 *
 * @return root
 */
template <class T>
T get_code_len_high_compo(T order, T n)
{
    if (order < n)
        assert(false);

    std::vector<T> factors;
    get_prime_factors<T>(order, &factors);
    T x = 1;
    typename std::vector<T>::size_type i, j;
    // forward to get a divisor of (q-1) >= n and of highly composited
    for (i = 0; i != factors.size(); ++i) {
        x *= factors.at(i);
        if (x >= n) {
            // backward to get smaller number
            for (j = 0; j != i + 1 && j != factors.size(); j++) {
                x /= factors.at(j);
                if (x < n)
                    return x * factors.at(j);
            }
        }
    }
    assert(false);
    return 0;
}

/**
 * find smallest number is
 *  - highly composited
 *  - at least n
 *  - divisible by (order)
 *
 * @param factors - vector of all prime factors of order. Replication of each
 * factor equals to its exponent
 * @param n
 *
 * @return code length
 */
template <class T>
T get_code_len_high_compo(std::vector<T>* factors, T n)
{
    assert(factors != nullptr);

    T x = 1;
    typename std::vector<T>::size_type i, j;
    // forward to get a divisor of (q-1) >= n and of highly composited
    for (i = 0; i != factors->size(); ++i) {
        x *= factors->at(i);
        if (x >= n) {
            // backward to get smaller number
            for (j = 0; j != i + 1 && j != factors->size(); j++) {
                x /= factors->at(j);
                if (x < n)
                    return x * factors->at(j);
            }
        }
    }
    assert(false);
    return 0;
}

/**
 * get all coprime factors of a number
 *
 * @param nb - a number to be refactored
 * @param output - vector of co-prime divisors of nb
 * @return
 */
template <class T>
void get_coprime_factors(T nb, std::vector<T>* output)
{
    std::vector<T> primes;
    std::vector<int> exponent;
    factor_prime<T>(nb, &primes, &exponent);

    typename std::vector<T>::size_type i;
    for (i = 0; i != primes.size(); ++i) {
        output->push_back(exp<T>(primes[i], exponent[i]));
    }
}

/**
 * get all prime factors of a number
 *
 * @param nb - a number to be refactored
 * @param output - vector of co-prime divisors of nb
 * @return
 */
template <class T>
void get_prime_factors(T nb, std::vector<T>* output)
{
    std::vector<T> primes;
    std::vector<int> exponent;
    factor_prime<T>(nb, &primes, &exponent);
    get_prime_factors_final<T>(&primes, &exponent, output);
}

/**
 * get all prime factors of a number
 *
 * @param primes - vector of primes
 * @param primes - vector of exponents each for prime
 * @param output - vector of co-prime divisors of nb
 */
template <class T>
void get_prime_factors_final(
    std::vector<T>* primes,
    std::vector<int>* exponent,
    std::vector<T>* output)
{
    typename std::vector<T>::size_type i;
    for (i = 0; i != primes->size(); ++i)
        for (int j = 0; j != exponent->at(i); ++j) {
            output->push_back(primes->at(i));
        }
}

} // namespace arith
} // namespace nttec

#endif
