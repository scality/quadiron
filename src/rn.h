/* -*- mode: c++ -*- */
#pragma once

template<typename T>
class Vec;

/**
 * Generic class for the Ring of integers modulo N
 *
 * @param n cardinal
 */
template<typename T>
class RN
{
private:
  T _card;
  T root;

public:
  RN(T card);
  virtual ~RN();
  virtual T card(void);
  virtual T card_minus_one(void) ;
  virtual bool check(T a) ;
  virtual T neg(T a) ;
  virtual T add(T a, T b) ;
  virtual T sub(T a, T b) ;
  virtual T mul(T a, T b) ;
  virtual T div(T a, T b) ;
  T inv_bezout(T a);
  virtual T inv(T a) ;
  virtual T exp(T a, T b) ;
  virtual T log(T a, T b) ;
  T exp_naive(T base, T exponent);
  T exp_quick(T base, T exponent);
  T log_naive(T base, T exponent);
  void compute_factors_of_order();
  bool is_quadratic_residue(T q);
  virtual void compute_omegas(Vec<T> *W, int n, T w);
  void compute_omegas_cached(Vec<T> *W, int n, T w);
  virtual T weak_rand(void);
  bool is_prime_root(T nb);
  void find_prime_root();
  T get_root();
  T get_prime_root();
  bool check_prime_root(T nb);
  bool check_order_naive(T nb, T order);
  T do_step_get_order(T x, T h, std::vector<T> *primes,
    std::vector<int> *exponent);
  T get_order(T x);
  virtual T get_nth_root(T n);
  T get_code_len(T n);
  T get_code_len_high_compo(T n);
 private:
  bool compute_factors_of_order_done = false;
  std::vector<T>* primes = NULL;
  std::vector<int>* exponent = NULL;
  std::vector<T>* all_primes_factors = NULL;
  std::vector<T>* proper_divisors = NULL;
};

template <typename T>
RN<T>::RN(T card)
{
  this->_card = card;
  this->root = 0;

  this->primes = new std::vector<T>();
  this->exponent = new std::vector<int>();
  this->all_primes_factors = new std::vector<T>();
  this->proper_divisors = new std::vector<T>();
}

template <typename T>
RN<T>::~RN()
{
  if (primes) delete primes;
  if (exponent) delete exponent;
  if (all_primes_factors) delete all_primes_factors;
  if (proper_divisors) delete proper_divisors;
}

template <typename T>
T RN<T>::card(void)
{
  return this->_card;
}

template <typename T>
T RN<T>::card_minus_one(void)
{
  return this->_card - 1;
}

template <typename T>
bool RN<T>::check(T a)
{
  return (a >= 0 && a < this->_card);
}

template <typename T>
T RN<T>::neg(T a)
{
  assert(check(a));

  return sub(0, a);
}

template <typename T>
T RN<T>::add(T a, T b)
{
  assert(check(a));
  assert(check(b));

  T c = a + b;
  if (c >= this->_card)
    return c - this->_card;
  else
    return c;
}

template <typename T>
T RN<T>::sub(T a, T b)
{
  assert(check(a));
  assert(check(b));

  if (a >= b)
    return a - b;
  else
    return this->_card - (b - a);
}

template <typename T>
T RN<T>::mul(T a, T b)
{
  assert(check(a));
  assert(check(b));

  return T((DoubleT<T>(a) * b) % this->_card);
}

template <typename T>
T RN<T>::div(T a, T b)
{
  assert(check(a));
  assert(check(b));

  T inv_b = inv(b);

  return mul(a, inv_b);
}

/**
 * Inverse by Bezout
 *
 * @param a
 *
 * @return
 */
template <typename T>
T RN<T>::inv_bezout(T a)
{
  assert(check(a));

  SignedDoubleT<T> x = a;
  SignedDoubleT<T> n = this->_card;
  SignedDoubleT<T> bezout[2];

  _extended_gcd<T>(x, n, bezout, NULL);
  if (bezout[0] < 0)
    bezout[0] = this->_card + bezout[0];
  return bezout[0];
}

template <typename T>
T RN<T>::inv(T a)
{
  return inv_bezout(a);
}

template <typename T>
T RN<T>::exp(T a, T b)
{
  assert(check(a));
  assert(check(b));

  return RN<T>::exp_quick(a, b);
}

template <typename T>
T RN<T>::log(T a, T b)
{
  assert(check(a));

  return RN<T>::log_naive(a, b);
}

/**
 * Naive exponentiation in the group
 *
 * @param base
 * @param exponent
 *
 * @return
 */
template <typename T>
T RN<T>::exp_naive(T base, T exponent)
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
 * Quick exponentiation in the group
 *
 * @param base
 * @param exponent
 *
 * @return
 */
template <typename T>
T RN<T>::exp_quick(T base, T exponent)
{
  T result;

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
 * Naive brute force algorithm in the group
 *
 * @param base
 * @param exponent
 *
 * @throw NTL_EX_NOT_FOUND if result is not found
 *
 * return
 */
template <typename T>
T RN<T>::log_naive(T base, T exponent)
{
  T result;

  for (result = 1; result < card(); result++) {
    if (exp(base, result) == exponent)
      return result;
  }

  // not found
  throw NTL_EX_NOT_FOUND;
}

template <typename T>
void RN<T>::compute_factors_of_order()
{
  if (this->compute_factors_of_order_done) return;

  T h = this->card_minus_one();
  // prime factorisation of order, i.e. order = p_i^e_i where
  //  p_i, e_i are ith element of this->primes and this->exponent
  _factor_prime<T>(h, this->primes, this->exponent);
  // store all primes in a vector. A prime is replicated according to its
  //  exponent
  _get_prime_factors_final<T>(this->primes, this->exponent,
                              this->all_primes_factors);
  // calculate all proper divisor of order. A proper divisor = order/p_i for
  //  each prime divisor of order.
  _get_proper_divisors<T>(h, this->primes, this->proper_divisors);

  this->compute_factors_of_order_done = true;
}

/**
 * check if q is a quadractic residue
 *
 * q is a quadratic residue, if x^2 == q % n
 *
 * @param q
 *
 * @return boolean
 */
template <typename T>
bool RN<T>::is_quadratic_residue(T q)
{
  T i;

  for (i = 0; i < this->card(); i++) {
    if (this->exp(i, 2) == q)
      return true;
  }

  return false;
}

/**
 * Compute the different powers of the root of unity into a vector
 *
 * @param W output vector (must be of length n+1)
 * @param n
 * @param w n-th root of unity
 */
template <typename T>
void RN<T>::compute_omegas(Vec<T> *W, int n, T w)
{
  for (int i = 0; i <= n; i++) {
    W->set(i, this->exp(w, i));
  }
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
void RN<T>::compute_omegas_cached(Vec<T> *W, int n, T w)
{
  std::ostringstream filename;

  filename << "W" << w << ".cache";

  if (-1 == access(filename.str().c_str(), F_OK)) {
    std::ofstream file;
    file.open(filename.str().c_str(), std::ios::out);
    for (int i = 0; i <= n; i++) {
      W->set(i, this->exp(w, i));
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

template <typename T>
T RN<T>::weak_rand(void)
{
  return _weak_rand<T>(this->card());
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
bool RN<T>::is_prime_root(T nb)
{
  bool ok = true;
  typename std::vector<T>::size_type i;
  T h = this->card_minus_one();
  if (!this->compute_factors_of_order_done) {
    this->compute_factors_of_order();
  }
  // check nb^divisor == 1
  // std::cout << "checking.." << nb << std::endl;
  for (i = 0; i != this->proper_divisors->size(); ++i) {
    // std::cout << nb << "^" << this->proper_divisors->at(i) << "=";
    // std::cout << exp(nb, this->proper_divisors->at(i)) << std::endl;
    if (this->exp(nb, this->proper_divisors->at(i)) == 1) {
      ok = false;
      break;
    }
  }
  return ok;
}

/*
 * Find primitive root of finite field. A number x is a primitive root if its
 *  order y = q - 1, i.e. x^i != 1 for every i < q-1
 *
 * Use the algorithm of checking if a number is primitive root or not
 */
template <typename T>
void RN<T>::find_prime_root()
{
  if (this->root) return;
  T nb = 2;
  bool ok;
  typename std::vector<T>::size_type i;
  T h = this->card_minus_one();
  if (!this->compute_factors_of_order_done) {
    this->compute_factors_of_order();
  }
  while (nb <= h) {
    ok = true;
    // check nb^divisor == 1
    // std::cout << "checking.." << nb << std::endl;
    for (i = 0; i != this->proper_divisors->size(); ++i) {
      // std::cout << nb << "^" << proper_divisors->at(i) << "=";
      // std::cout << this->exp(nb, this->proper_divisors->at(i)) << std::endl;
      if (this->exp(nb, this->proper_divisors->at(i)) == 1) {
        ok = false;
        break;
      }
    }
    if (ok) {
      this->root = nb;
      break;
    }
    nb++;
  }

  if (!this->root)
    assert(false);  // not found root
}

template <typename T>
T RN<T>::get_root()
{
  return root;
}

template <typename T>
T RN<T>::do_step_get_order(T x, T h, std::vector<T> *primes,
  std::vector<int> *exponent)
{
  T y, p;
  int r;
  while (!primes->empty()) {
    p = primes->back();
    r = exponent->back();
    primes->pop_back();
    exponent->pop_back();
    y = h/p;
    if (this->exp(x, y) != 1) {
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
 * Get order of an element x of RN(q)
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
T RN<T>::get_order(T x)
{
  if (x == 0 || x == 1) return 1;
  T h = this->card_minus_one();
  if (!this->compute_factors_of_order_done)
    _factor_prime<T>(h, this->primes, this->exponent);
  std::vector<T> _primes(*primes);
  std::vector<int> _exponent(*exponent);
  T order = do_step_get_order(x, h, &_primes, &_exponent);

  if (order == 1) return h;
  return order;
}

/*
 * Check whether a number is a primitive root
 */
template <typename T>
bool RN<T>::check_prime_root(T nb)
{
  T h = this->card_minus_one();
  return (get_order(nb) == h);
}

/*
 * Check naively order of a number
 */
template <typename T>
bool RN<T>::check_order_naive(T nb, T order)
{
  if (this->exp(nb, order) != 1) return false;

  T i = 1;
  T tmp = nb;
  while (i < order - 1) {
    // std::cout << i << ":" << tmp << std::endl;
    if (tmp == 1) return false;
    tmp = this->mul(tmp, nb);
    i++;
  }
  return true;
}

/**
 * compute root of order n: g^((q-1)/d) where d = gcd(n, q-1)
 *
 * @param n
 *
 * @return root
 */
template <typename T>
T RN<T>::get_nth_root(T n)
{
  T q_minus_one = this->card_minus_one();
  T d = _gcd<T>(n, q_minus_one);
  if (!this->root)
    this->find_prime_root();
  T nth_root = this->exp(this->root, q_minus_one/d);
  return nth_root;
}

/**
 * return prime root

 * @return root
 */
template <typename T>
T RN<T>::get_prime_root()
{
  if (!this->root)
    this->find_prime_root();
  return this->root;
}

/**
 * find smallest number is
 *  - at least n
 *  - divisible by (q-1)
 *
 * @param n
 *
 * @return root
 */
template <typename T>
T RN<T>::get_code_len(T n)
{
  T nb = this->card_minus_one();
  if (nb < n) assert(false);

  return _get_code_len<T>(nb, n);
}

/**
 * find smallest number is
 *  - highly composited
 *  - at least n
 *  - divisible by (q-1)
 *
 * @param n
 *
 * @return root
 */
template <typename T>
T RN<T>::get_code_len_high_compo(T n)
{
  T nb = this->card_minus_one();
  if (nb < n) assert(false);

  if (!this->compute_factors_of_order_done) {
    this->compute_factors_of_order();
  }
  return _get_code_len_high_compo<T>(this->all_primes_factors, n);
}
