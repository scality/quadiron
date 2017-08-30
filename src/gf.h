/* -*- mode: c++ -*- */
#pragma once

template<typename T>
using DoubleT = typename Double<T>::T;

template<typename T>
using SignedDoubleT = typename SignedDouble<T>::T;

template<typename T>
class Vec;

/**
 * Generic class for Galois Fields of q=p^n
 *
 * @note Functions starting with _ are convenience function that does
 *       not call GF specific operations
 *
 * @param p primer number
 * @param n exponent
 */
template<typename T>
class GF
{
 public:
  T p;
  T n;
  T root;

 protected:
  GF(T p, T n);

 public:
  virtual T card(void) = 0;
  virtual T card_minus_one(void) = 0;
  virtual bool check(T a) = 0;
  virtual T max(T a, T b) = 0;
  virtual T min(T a, T b) = 0;
  virtual T neg(T a) = 0;
  virtual T add(T a, T b) = 0;
  virtual T sub(T a, T b) = 0;
  virtual T mul(T a, T b) = 0;
  virtual T div(T a, T b) = 0;
  virtual T inv(T a) = 0;
  virtual T exp(T a, T b) = 0;
  virtual T log(T a, T b) = 0;
  T _card();
  T exp_naive(T base, T exponent);
  T exp_quick(T base, T exponent);
  T _sqrt(T n);
  T _exp(T base, T exponent);
  T _exp_mod(T base, T exponent, T modulus);
  T log_naive(T base, T exponent);
  bool _is_power_of_2(int x);
  int _log2(int x);
  int _exp2(int x);
  SignedDoubleT<T> _extended_gcd(SignedDoubleT<T> a, SignedDoubleT<T> b,
    SignedDoubleT<T> bezout_coef[2], SignedDoubleT<T> quotient_gcd[2]);
  T _chinese_remainder(int n_mod, T a[], T n[]);
  bool is_quadratic_residue(T q);
  int _jacobi(SignedDoubleT<T> n, SignedDoubleT<T> m);
  bool _solovay_strassen1(T a, T n);
  bool _solovay_strassen(T n);
  bool _is_prime(T n);
  T get_nth_root(int n, int R);
  void compute_omegas(Vec<T> *W, int n, T w);
  T _weak_rand(T max);
  T weak_rand(void);
  void _factor_distinct_prime(T n, std::vector<T> *output);
  void _factor_prime(T nb, std::vector<T> *primes, std::vector<T> *exponent);
  void _get_proper_divisors(T n, std::vector<T> *output);
  bool is_prime_root(T nb);
  void find_prime_root();
  bool check_prime_root(T nb);
  bool check_prime_root_naive(T nb);
  T do_step_get_order(T x, T h, std::vector<T> *primes,
    std::vector<T> *exponent);
  T get_order(T x);
};

template <typename T>
GF<T>::GF(T p, T n)
{
  this->p = p;
  this->n = n;
  this->root = 0;
}

template <typename T>
T GF<T>::_card()
{
  return _exp(p, n);
}

/**
 * Naive exponentiation in the field
 *
 * @param gf
 * @param base
 * @param exponent
 *
 * @return
 */
template <typename T>
T GF<T>::exp_naive(T base, T exponent)
{
  T result;
  T i;

  if (0 == exponent)
    return 1;

  if (1 == exponent)
    return base;

  result = base;
  for (i = 1; i < exponent; i++)
    result = this->mul(result, base);

  return result;
}

/**
 * Quick exponentiation in the field
 *
 * @param gf
 * @param base
 * @param exponent
 *
 * @return
 */
template <typename T>
T GF<T>::exp_quick(T base, T exponent)
{
  T result;
  T i;

  if (0 == exponent)
    return 1;

  if (1 == exponent)
    return base;

  T tmp = this->exp_quick(base, exponent / 2);
  result = this->mul(tmp, tmp);
  if (exponent % 2 == 1)
    result = this->mul(result, base);
  return result;
}

/**
 * integer square root (from Wikipedia)
 *
 * @param n
 *
 * @return
 */
template <typename T>
T GF<T>::_sqrt(T remainder)
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
T GF<T>::_exp(T base, T exponent)
{
  T result = 1;
  while (exponent > 0) {
    if (exponent % 2 == 1)
      result *= base;
    exponent >>= 1;
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
T GF<T>::_exp_mod(T base, T exponent, T modulus)
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
 * Naive brute force algorithm
 *
 * @param gf
 * @param base
 * @param exponent
 *
 * @throw NTL_EX_NOT_FOUND if result is not found
 *
 * return
 */
template <typename T>
T GF<T>::log_naive(T base, T exponent)
{
  T result;

  for (result = 1; result < card(); result++) {
    if (exp(base, result) == exponent)
      return result;
  }

  // not found
  throw NTL_EX_NOT_FOUND;
}

/**
 * check if x is a power of 2
 *
 * @param x
 *
 * @return
 */
template <typename T>
bool GF<T>::_is_power_of_2(int x)
{
  return x > 0 && !(x & (x-1));
}

/**
 * Compute log2(x)
 *
 * @param exponent
 *
 * @return
 */
template <typename T>
int GF<T>::_log2(int x)
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
int GF<T>::_exp2(int x)
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
SignedDoubleT<T> GF<T>::_extended_gcd(SignedDoubleT<T> a, SignedDoubleT<T> b,
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
T GF<T>::_chinese_remainder(int n_mod, T a[], T n[])
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
    _extended_gcd(N[i], n[i], bezout, NULL);
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
 * check if q is a quadractic residue in GF_p^n
 *
 * q is a quadratic residue, if x^2 == q % n
 *
 * @param q
 *
 * @return boolean
 */
template <typename T>
bool GF<T>::is_quadratic_residue(T q)
{
  T i;

  for (i = 0; i < card(); i++) {
    if (exp(i, 2) == q)
      return true;
  }

  return false;
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
int GF<T>::_jacobi(SignedDoubleT<T> n, SignedDoubleT<T> m)
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
bool GF<T>::_solovay_strassen1(T a, T n)
{
  T _n = (n - 1) / 2;
  T _j = _exp_mod(a, _n, n);
  int j = _jacobi(a, n);
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
bool GF<T>::_solovay_strassen(T n)
{
  int ok = 0;
  for (int i = 0; i < 100; i++) {
    T a = _weak_rand(n);
    if (!_solovay_strassen1(a, n))
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
bool GF<T>::_is_prime(T n)
{
  if (n == 2)
    return true;

  T root = _sqrt(n);
  for (T i = 2; i <= root; i++) {
    if (n % i == 0)
      return false;
  }
  return true;
}

/**
 * compute root of order n-1 such as r^(n-1) mod q == 1
 * As suggested by the paper:
 * FNT-based Reed-Solomon Erasure Codes
 * by Alexandre Soro and Jerome Lacan
 *
 * @param R primitive root
 * @param n must be a power of 2
 *
 * @return root
 */
template <typename T>
T GF<T>::get_nth_root(int n, int R)
{
  assert(_is_power_of_2(n));
  u_int l = _log2(n);
  assert(l <= 16);
  mpz_class _R = R;
  std::stringstream my_card;
  my_card << card();
  mpz_class _q(my_card.str());
  mpz_class _r;
  mpz_powm_ui(_r.get_mpz_t(), _R.get_mpz_t(), _exp2(16-l), _q.get_mpz_t());
  // std::cout << "n=" << n << " l=" << l << " r=" << _r << " q=" << _q << "\n";
  return _r.get_ui();
}

/**
 * Compute the different powers of the root of unity into a vector
 *
 * @note cache the result in a file called W<w>.cache
 *
 * @note XXX not reentrant
 *
 * @param W output vector (must be of length n+1)
 * @param n
 * @param w n-th root of unity
 */
template <typename T>
void GF<T>::compute_omegas(Vec<T> *W, int n, T w)
{
  std::ostringstream filename;

  filename << "W" << w << ".cache";

  if (-1 == access(filename.str().c_str(), F_OK)) {
    std::ofstream file;
    file.open(filename.str().c_str(), std::ios::out);
    for (int i = 0; i <= n; i++) {
      W->set(i, exp(w, i));
      file << W->get(i) << "\n";
    }
  } else {
    std::ifstream file;
    int i = 0;
    file.open(filename.str().c_str(), std::ios::in);
    T tmp;
    while (file >> tmp) {
      W->set(i, tmp);
      i++;
    }
    assert(i == n + 1);
  }
}

/**
 * Returns a number n such as 0 < n < max
 *
 * @param max
 *
 * @return
 */
template <typename T>
T GF<T>::_weak_rand(T max)
{
  T r;
 retry:
  r = rand() % max;
  if (0 == r)
    goto retry;
  return r;
}

template <typename T>
T GF<T>::weak_rand(void)
{
  return _weak_rand(card());
}

/*
 * A given number `n` is factored into primes-> Only primes are stored, their
 *  exponent is ignored.
 */
template<typename T>
void GF<T>::_factor_distinct_prime(T nb, std::vector<T> *output)
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
  for (T i = 3; i <= sqrt(nb); i = i + 2) {
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
void GF<T>::_get_proper_divisors(T nb, std::vector<T> *output)
{
  std::vector<T> *input = new std::vector<T>();
  typename std::vector<T>::iterator it;
  _factor_distinct_prime(nb, input);
  // std::cout << "nb: " << nb << std::endl;
  for (it = input->begin(); it != input->end(); ++it) {
    if (*it < nb) {
      output->push_back(nb / (*it));
      // std::cout << *it << std::endl;
    }
  }
}

/*
 * Check if a number is primitive root or not, i.e. check its
 *  order y = q - 1, i.e. x^i != 1 for every i < q-1
 *
 * Note, order of a number must be a divisor of (q-1)
 *
 * Algorithm to checking whether x is primitive root or not
 *  Methodology: checking its order is not a proper divisor of (q-1)
 * 1. Prime factorization (q - 1) = p1^r1 * p2^r2 * ... pm^rm
 * 2. Get "necessary" proper divisors of (q - 1):
 *      D = {(q-1)/p1, (q-1)/p2, .., (q-1)/pm }
 * 3. x is a primitive root if for every divisor d of D:
 *                          x^d != 1
 * Proof:
 *  Suppose that x satisfies the algorithm but its order is a proper divisor y
 *  of (q-1). This order can be expressed as (q-1)/y where
 *              y = (p1^s1 * p2^s2 * ... * pm^sm)
 *  and, x^((q-1)/y) = 1
 *  Withouth loss of generality, we suppose s1 >= 1. We have
 *      x^((q-1)/y) = 1 => (x^((q-1)/y))^(y/p1) = 1 <=> x^((q-1)/p1) = 1.
 *  Contradict to Step 3 of the algorith.
 */
template <typename T>
bool GF<T>::is_prime_root(T nb)
{
  std::vector<T>* divisors = new std::vector<T>();
  typename std::vector<T>::iterator divisor;
  T h = card_minus_one();
  // get necessary proper divisors of h
  _get_proper_divisors(h, divisors);
  // check nb^divisor == 1
  // std::cout << "checking.." << nb << std::endl;
  for (divisor = divisors->begin(); divisor != divisors->end(); ++divisor) {
    // std::cout << nb << "^" << *divisor << "=" << exp(nb, *divisor) << std::endl;
    if (exp(nb, *divisor) == 1) {
      return false;
    }
  }
  return true;
}

/*
 * Find primitive root of finite field. A number x is a primitive root if its
 *  order y = q - 1, i.e. x^i != 1 for every i < q-1
 *
 * Use the algorithm of checking if a number is primitive root or not
 */
template <typename T>
void GF<T>::find_prime_root()
{
  std::vector<T>* divisors = new std::vector<T>();
  typename std::vector<T>::iterator divisor;
  T h = card_minus_one();
  // get all divisors of h
  _get_proper_divisors(h, divisors);
  T nb = 2;
  bool ok;
  while (nb <= h) {
    ok = true;
    // check nb^divisor == 1
    // std::cout << "checking.." << nb << std::endl;
    for (divisor = divisors->begin(); divisor != divisors->end(); ++divisor) {
      // std::cout << nb << "^" << *divisor << "=" << exp(nb, *divisor) << std::endl;
      if (exp(nb, *divisor) == 1) {
        ok = false;
        break;
      }
    }
    if (ok) {
      this->root = nb;
      return;
    }
    nb++;
  }
  if (!this->root)
    assert(false);  // not found root
}

/*
 * A given number `n` is factored into primes->
 */
template<typename T>
void GF<T>::_factor_prime(T nb, std::vector<T> *primes,
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
  for (T i = 3; i <= sqrt(nb); i = i + 2) {
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

template <typename T>
T GF<T>::do_step_get_order(T x, T h, std::vector<T> *primes,
  std::vector<T> *exponent)
{
  T y, p, r;
  while (!primes->empty()) {
    p = primes->back();
    r = exponent->back();
    primes->pop_back();
    exponent->pop_back();
    y = h/p;
    if (exp(x, y) != 1) {
      // remove (p, r)
      while (r > 1) {
        y /= p;
        r--;
      }
      continue;
    }
    // exp(x, y) == 1
    if (r > 1) {
      primes->push_back(p);
      exponent->push_back(r - 1);
    }
    // do next
    return do_step_get_order(x, y, primes, exponent);
  }
  return h;
}

/*
 * Get order of an element x of GF(q)
 *  The order d is the smallest divisor of (q-1) such that x^d = 1
 * Pseudocode:
 *  Prime factorisation of (q-1) = p1^r1 * p2^r2 * ... * pm^rm
 *  D = { (p1, r1), (p2, r2), .. , (pm, rm) }
 *  h = q - 1
 *   function do_step_get_order(x, h, D) {
 *     for (pi, ri) in D {
 *       y = h/pi
 *       newD = D \ (pi, ri)
 *       if (x^y != 1) {
 *           y /= pi^(ri - 1)
 *           continue;
 *       }
 *       // x^y == 1
 *       if (ri > 1)
 *         newD = newD + (pi, ri - 1)
 *       return do_step_get_order(x, y, newD)
 *     }
 *     return h;
 *   }
 */
template <typename T>
T GF<T>::get_order(T x)
{
  std::vector<T> *primes = new std::vector<T>();
  std::vector<T> *exponent = new std::vector<T>();
  T h = card_minus_one();
  _factor_prime(h, primes, exponent);
  T order = do_step_get_order(x, h, primes, exponent);
  if (order == 1) return h;
  return order;
}

/*
 * Check whether a number is a primitive root
 */
template <typename T>
bool GF<T>::check_prime_root(T nb)
{
  T h = card_minus_one();
  return (get_order(nb) == h);
}

/*
 * Check naively whether a number is a primitive root
 */
template <typename T>
bool GF<T>::check_prime_root_naive(T nb)
{
  T h = card_minus_one();
  T i = 1;
  T tmp = nb;
  while (i < h) {
    // std::cout << i << ":" << tmp << std::endl;
    if (tmp == 1) return false;
    tmp = mul(tmp, nb);
    i++;
  }
  // std::cout << "last one" << tmp << std::endl;
  if (tmp == 1) return true;
  return false;
}
