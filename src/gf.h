/* -*- mode: c++ -*- */
#pragma once

template<typename T>
class Arith;

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
  Arith<T> *arith;

  GF(T p, T n);
  virtual ~GF();

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
  T log_naive(T base, T exponent);
  bool is_quadratic_residue(T q);
  void compute_omegas(Vec<T> *W, int n, T w);
  void compute_omegas_cached(Vec<T> *W, int n, T w);
  T weak_rand(void);
  bool is_prime_root(T nb);
  void find_prime_root();
  T get_prime_root();
  bool check_prime_root(T nb);
  bool check_order_naive(T nb, T order);
  T do_step_get_order(T x, T h, std::vector<T> *primes,
    std::vector<T> *exponent);
  T get_order(T x);
  T get_nth_root(T n);
  T get_code_len(T n);
  T get_code_len_high_compo(T n);
  void compute_all_divisors_h();
 private:
  std::vector<T>* _all_divisors_h = NULL;
};

template <typename T>
GF<T>::GF(T p, T n)
{
  this->p = p;
  this->n = n;
  this->root = 0;
  arith = new Arith<T>();
}

template <typename T>
GF<T>::~GF()
{
  delete arith;
  delete _all_divisors_h;
}

template <typename T>
T GF<T>::_card()
{
  return arith->exp(p, n);
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
 * Compute the different powers of the root of unity into a vector
 *
 * @param W output vector (must be of length n+1)
 * @param n
 * @param w n-th root of unity
 */
template <typename T>
void GF<T>::compute_omegas(Vec<T> *W, int n, T w)
{
  for (int i = 0; i <= n; i++) {
    W->set(i, exp(w, i));
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
void GF<T>::compute_omegas_cached(Vec<T> *W, int n, T w)
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

template <typename T>
T GF<T>::weak_rand(void)
{
  return arith->weak_rand(card());
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
  bool ok = true;
  std::vector<T>* divisors = new std::vector<T>();
  typename std::vector<T>::iterator divisor;
  T h = card_minus_one();
  // get necessary proper divisors of h
  arith->get_proper_divisors(h, divisors);
  // check nb^divisor == 1
  // std::cout << "checking.." << nb << std::endl;
  for (divisor = divisors->begin(); divisor != divisors->end(); ++divisor) {
    // std::cout << nb << "^" << *divisor << "=" << exp(nb, *divisor) << std::endl;
    if (exp(nb, *divisor) == 1) {
      ok = false;
      break;
    }
  }

  delete divisors;
  return ok;
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
  if (this->root) return;
  std::vector<T>* divisors = new std::vector<T>();
  typename std::vector<T>::iterator divisor;
  T h = card_minus_one();
  // get all divisors of h
  arith->get_proper_divisors(h, divisors);
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
      break;
    }
    nb++;
  }

  delete divisors;

  if (!this->root)
    assert(false);  // not found root
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
  if (x == 0 || x == 1) return 1;
  std::vector<T> *primes = new std::vector<T>();
  std::vector<T> *exponent = new std::vector<T>();
  T h = card_minus_one();
  arith->factor_prime(h, primes, exponent);
  T order = do_step_get_order(x, h, primes, exponent);

  delete primes;
  delete exponent;

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
 * Check naively order of a number
 */
template <typename T>
bool GF<T>::check_order_naive(T nb, T order)
{
  if (exp(nb, order) != 1) return false;

  T i = 1;
  T tmp = nb;
  while (i < order - 1) {
    // std::cout << i << ":" << tmp << std::endl;
    if (tmp == 1) return false;
    tmp = mul(tmp, nb);
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
T GF<T>::get_nth_root(T n)
{
  T q_minus_one = this->card_minus_one();
  T d = arith->gcd(n, q_minus_one);
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
T GF<T>::get_prime_root()
{
  if (!this->root)
    this->find_prime_root();
  return this->root;
}

/**
 * compute all divisors of q-1
 *
 */
template <typename T>
void GF<T>::compute_all_divisors_h()
{
  if (this->_all_divisors_h != NULL) return;
  this->_all_divisors_h = new std::vector<T>();
  arith->compute_all_divisors(this->card_minus_one(), this->_all_divisors_h);
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
T GF<T>::get_code_len(T n)
{
  T nb = card_minus_one();
  if (nb < n) assert(false);

  return arith->get_code_len(nb, n);
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
T GF<T>::get_code_len_high_compo(T n)
{
  T nb = card_minus_one();
  if (nb < n) assert(false);

  return arith->get_code_len_high_compo(nb, n);
}
