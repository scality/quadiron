/* -*- mode: c++ -*- */
#pragma once

template<typename T>
using DoubleT = typename Double<T>::T;

template<typename T>
using SignedDoubleT = typename SignedDouble<T>::T;

/**
 * Class for arithmetics
 *
 */
template<typename T>
class Arith
{
 public:
  Arith() {};
  virtual ~Arith() {};

 public:
  T sqrt(T n);
  T exp(T base, T exponent);
  T exp_mod(T base, T exponent, T modulus);
  bool is_power_of_2(int x);
  T get_smallest_power_of_2(int x);
  int log2(int x);
  int exp2(int x);
  SignedDoubleT<T> extended_gcd(SignedDoubleT<T> a, SignedDoubleT<T> b,
    SignedDoubleT<T> bezout_coef[2], SignedDoubleT<T> quotient_gcd[2]);
  T chinese_remainder(int n_mod, T a[], T n[]);
  int jacobi(SignedDoubleT<T> n, SignedDoubleT<T> m);
  bool solovay_strassen1(T a, T n);
  bool solovay_strassen(T n);
  bool is_prime(T n);
  T weak_rand(T max);
  T gcd(T u, T v);
  void factor_distinct_prime(T n, std::vector<T> *output);
  void factor_prime(T nb, std::vector<T> *primes, std::vector<T> *exponent);
  void get_proper_divisors(T n, std::vector<T> *output);
  void get_proper_divisors(T n, std::vector<T> *primes, std::vector<T> *output);
  void compute_all_divisors(T nb, std::vector<T> *output);
  T get_code_len(T order, T n);
  T get_code_len_high_compo(T order, T n);
  T get_code_len_high_compo(std::vector<T> *factors, T n);
  void get_coprime_factors(T nb, std::vector<T> *output);
  void get_prime_factors(T nb, std::vector<T> *output);
  void get_prime_factors_final(std::vector<T> *primes,
    std::vector<T> *exponent, std::vector<T> *output);
};

/**
 * integer square root (from Wikipedia)
 *
 * @param n
 *
 * @return
 */
template <typename T>
T Arith<T>::sqrt(T remainder)
{
  // calculated by precompiler = same runtime as: place = 0x40000000
  T place = (T)1 << (sizeof (T) * 8 - 2);
  while (place > remainder)
    place /= 4;  // optimized by complier as place >>= 2

  T root = 0;
  while (place)
    {
      if (remainder >= root+place)
        {
          remainder -= root+place;
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
template <typename T>
T Arith<T>::exp(T base, T exponent)
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
template <typename T>
T Arith<T>::exp_mod(T base, T exponent, T modulus)
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
template <typename T>
bool Arith<T>::is_power_of_2(int x)
{
  return x > 0 && !(x & (x-1));
}

/**
 * Get a smallest power of 2 and >= x
 *
 * @param x
 *
 * @return
 */
template <typename T>
T Arith<T>::get_smallest_power_of_2(int x)
{
  if (is_power_of_2(x))
    return x;
  return exp2(log2(x) + 1);
}

/**
 * Compute log2(x)
 *
 * @param exponent
 *
 * @return
 */
template <typename T>
int Arith<T>::log2(int x)
{
  int result = 0;

  if (x == 0)
    throw NTL_EX_INVAL;

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
template <typename T>
int Arith<T>::exp2(int x)
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
 * @param bezout_coef[output] computed bezout coefficients, might be NULL
 * @param quotient_gcd[output] computed quotient by the GCD, might be NULL
 *
 * XXX take care of the signs of input
 *
 * @return the GCD
 */
template <typename T>
SignedDoubleT<T> Arith<T>::extended_gcd(SignedDoubleT<T> a, SignedDoubleT<T> b,
  SignedDoubleT<T> bezout_coef[2], SignedDoubleT<T> quotient_gcd[2])
{
  SignedDoubleT<T> s = 0;
  SignedDoubleT<T> old_s = 1;
  SignedDoubleT<T> t = 1;
  SignedDoubleT<T> old_t = 0;
  SignedDoubleT<T> r = b;
  SignedDoubleT<T> old_r = a;
  SignedDoubleT<T> quotient;
  SignedDoubleT<T> tmp;

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

/**
 * Chinese remainder theorem (from Wikipedia)
 *
 * @param n_mod number of moduli
 * @param a the a's (integers)
 * @param n the n's (moduli)
 *
 * XXX check if there is a solution
 *
 * @return the solution
 * @throw NTL_EX_NO_SOLUTION
 */
template <typename T>
T Arith<T>::chinese_remainder(int n_mod, T a[], T n[])
{
  int i;
  SignedDoubleT<T> _N, x;
  SignedDoubleT<T> *N = new SignedDoubleT<T>[n_mod];
  SignedDoubleT<T> *M = new SignedDoubleT<T>[n_mod];

  _N = 1;
  for (i = 0; i < n_mod; i++)
    _N *= n[i];

  for (i = 0; i < n_mod; i++) {
    SignedDoubleT<T> bezout[2];

    N[i] = _N / n[i];
    extended_gcd(N[i], n[i], bezout, NULL);
    M[i] = bezout[0];  // XXX % n[i];
  }

  x = 0;
  for (i = 0; i < n_mod; i++)
    x += N[i] * M[i] * a[i];

  delete[] N;
  delete[] M;

  x = x % _N;

  if (x < 0)
    x = _N + x;

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
 * @throw NT_EX_INVAL if b is even or b is negative
 */
template <typename T>
int Arith<T>::jacobi(SignedDoubleT<T> n, SignedDoubleT<T> m)
{
  SignedDoubleT<T> t;
  int jac;

  // m must be odd
  if ((m % 2) == 0)
    throw NTL_EX_INVAL;

  // m must be positive
  if (m < 0)
    throw NTL_EX_INVAL;

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
bool Arith<T>::solovay_strassen1(T a, T n)
{
  T _n = (n - 1) / 2;
  T _j = exp_mod(a, _n, n);
  int j = jacobi(a, n);
  if (j < 0)
    _j = _j - n;
  return j == _j;
}

/**
 * Perform the Solvay-Strassen primality test
 *
 * @param n check if n is prime
 *
 * @return true if n is (probably) prime else false
 */
template <typename T>
bool Arith<T>::solovay_strassen(T n)
{
  int ok = 0;
  for (int i = 0; i < 100; i++) {
    T a = weak_rand(n);
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
template <typename T>
bool Arith<T>::is_prime(T n)
{
  if (n == 2)
    return true;

  T root = this->sqrt(n);
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
template <typename T>
T Arith<T>::weak_rand(T max)
{
  T r;
 retry:
  r = rand() % max;
  if (0 == r)
    goto retry;
  return r;
}

/*
 * A given number `n` is factored into primes-> Only primes are stored, their
 *  exponent is ignored.
 */
template<typename T>
void Arith<T>::factor_distinct_prime(T nb, std::vector<T> *output)
{
  T last_found = 1;
  while (nb % 2 == 0) {
    if (last_found != 2) {
      output->push_back(2);
      last_found = 2;
    }
    nb = nb/2;
  }
  // n must be odd at this point.  So we can skip one element
  for (T i = 3; i <= this->sqrt(nb); i = i + 2) {
    // While i divides n, get i and divide n
    while (nb % i == 0) {
      if (last_found != i) {
        output->push_back(i);
        last_found = i;
      }
      nb = nb/i;
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
template<typename T>
void Arith<T>::get_proper_divisors(T nb, std::vector<T> *output)
{
  std::vector<T> input;
  typename std::vector<T>::iterator it;
  factor_distinct_prime(nb, &input);
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
template<typename T>
void Arith<T>::get_proper_divisors(T nb, std::vector<T> *primes,
  std::vector<T> *output)
{
  assert(primes != NULL);
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
template<typename T>
void Arith<T>::factor_prime(T nb, std::vector<T> *primes,
  std::vector<T> *exponent)
{
  // std::cout << nb << ": ";
  T occurence = 0;
  while (nb % 2 == 0) {
    occurence++;
    if (occurence == 1) {
      primes->push_back(2);
    }
    nb = nb/2;
  }
  if (occurence > 0) {
    exponent->push_back(occurence);
    occurence = 0;
  }
  // n must be odd at this point.  So we can skip one element
  for (T i = 3; i <= this->sqrt(nb); i = i + 2) {
    // While i divides n, get i and divide n
    while (nb % i == 0) {
      occurence++;
      if (occurence == 1) {
        primes->push_back(i);
      }
      nb = nb/i;
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
template <typename T>
T Arith<T>::gcd(T u, T v)
{
  T r;
  if (v == 0) return u;
  if (u == 0) return v;
  while (1) {
    r = u % v;
    if (r == 0) return v;
    u = v;
    v = r;
  }
}

/**
 * compute all divisors of a number
 *
 */
template <typename T>
void Arith<T>::compute_all_divisors(T nb, std::vector<T>* output)
{
  typename std::vector<T> tmp;

  T nb_sqrt = this->sqrt(nb);
  for (T i = 1; i <= nb_sqrt; i++) {
    // While i divides n, get i and n/i
    if (nb % i == 0) {
      // std::cout << nb << ":" << i << " " << nb/i << std::endl;
      output->push_back(i);
      tmp.push_back(nb / i);
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
template <typename T>
T Arith<T>::get_code_len(T order, T n)
{
  if (order % n == 0) return n;
  if (order < n) assert(false);
  T nb_sqrt = this->sqrt(order);
  T i;
  if (n > nb_sqrt) {
    for (i = order/n; i >= 1; i--) {
      // if i divides n, return n/i
      if (order % i == 0)
        return order/i;
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
template <typename T>
T Arith<T>::get_code_len_high_compo(T order, T n)
{
  if (order < n) assert(false);

  std::vector<T> factors;
  get_prime_factors(order, &factors);
  T x = 1;
  typename std::vector<T>::size_type i, j;
  // forward to get a divisor of (q-1) >= n and of highly composited
  for (i = 0; i != factors.size(); ++i) {
    x *= factors.at(i);
    if (x >= n) {
      // backward to get smaller number
      for (j = 0; j != i+1 && j != factors.size(); j++) {
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
template <typename T>
T Arith<T>::get_code_len_high_compo(std::vector<T> *factors, T n)
{
  assert(factors != NULL);

  T x = 1;
  typename std::vector<T>::size_type i, j;
  // forward to get a divisor of (q-1) >= n and of highly composited
  for (i = 0; i != factors->size(); ++i) {
    x *= factors->at(i);
    if (x >= n) {
      // backward to get smaller number
      for (j = 0; j != i+1 && j != factors->size(); j++) {
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
template <typename T>
void Arith<T>::get_coprime_factors(T nb, std::vector<T> *output)
{
  std::vector<T> primes;
  std::vector<T> exponent;
  factor_prime(nb, &primes, &exponent);

  typename std::vector<T>::size_type i;
  for (i = 0; i != primes.size(); ++i) {
    output->push_back(exp(primes[i], exponent[i]));
  }
}

/**
 * get all prime factors of a number
 *
 * @param nb - a number to be refactored
 * @param output - vector of co-prime divisors of nb
 * @return
 */
template <typename T>
void Arith<T>::get_prime_factors(T nb, std::vector<T> *output)
{
  std::vector<T> primes;
  std::vector<T> exponent;
  factor_prime(nb, &primes, &exponent);
  get_prime_factors_final(&primes, &exponent, output);
}

/**
 * get all prime factors of a number
 *
 * @param primes - vector of primes
 * @param primes - vector of exponents each for prime
 * @param output - vector of co-prime divisors of nb
 * @return
 */
template <typename T>
void Arith<T>::get_prime_factors_final(std::vector<T> *primes,
  std::vector<T> *exponent, std::vector<T> *output)
{
  typename std::vector<T>::size_type i, j;
  for (i = 0; i != primes->size(); ++i)
    for (j = 0; j != exponent->at(i); ++j) {
      output->push_back(primes->at(i));
    }
}
