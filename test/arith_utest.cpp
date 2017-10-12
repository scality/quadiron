
#include "ntl.h"

template<typename T>
class ArithUtest
{
 public:
  T max;
  Arith<T> *arith;
  ArithUtest() {
    if (sizeof(T) < 16)
      max = std::numeric_limits<T>::max();
    else {
      max = 0;
      for (int i = 0; i < sizeof(T)*8; i++)
        max += (T)1 << i;
    }
    arith = new Arith<T>();
  }

  ~ArithUtest() {
    delete arith;
  }

  void test_basic_ops()
  {
    std::cout << "test_basic_ops\n";
    assert(arith->sqrt(2025) == 45);
    assert(arith->is_prime(2));
    assert(arith->is_prime(3));
    assert(arith->is_prime(13));
    assert(!arith->is_prime(4));
    assert(!arith->is_prime(15));
  }

  void test_reciprocal()
  {
    std::cout << "test_reciprocal\n";
    int i;
    T sub_max = arith->sqrt(max);
    for (i = 0; i < 1000; i++) {
      T x, y, z;

      // std::cout << "i=" << i << "\n";

      x = arith->weak_rand(sub_max);
      // std::cout << "x=" << x << "\n";
      y = arith->exp(x, 2);
      // std::cout << "exp(x)=" << y << "\n";
      z = arith->sqrt(y);
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
    omega = arith->chinese_remainder(2, a, n);
    assert(omega == 5996);

    a[0] = 6;
    n[0] = 7;
    a[1] = 4;
    n[1] = 8;
    omega = arith->chinese_remainder(2, a, n);
    assert(omega == 20);

    a[0] = 3;
    n[0] = 4;
    a[1] = 0;
    n[1] = 6;
    omega = arith->chinese_remainder(2, a, n);
    // no solution XXX detect it
  }

  void test_jacobi()
  {
    std::cout << "test_jacobi\n";
    assert(arith->jacobi(1001, 9907) == -1);
    assert(arith->jacobi(19, 45) == 1);
    assert(arith->jacobi(8, 21) == -1);
    assert(arith->jacobi(5, 21) == 1);
    assert(arith->jacobi(47, 221) == -1);
    assert(arith->jacobi(2, 221) == -1);
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
    if (sizeof(T) != 8) return;
    std::cout << "test_mul_bignum\n";

    int b = 10;  // base
    int p = 14;  // we could multiply integers of 2^p digits
    T max_digits = arith->exp(2, p);
    // std::cerr << "p=" << p << " max_digits=" << max_digits << "\n";

    T l = p + 1;
    // std::cerr << "l=" << l << "\n";

    // choose 2 prime numbers of the form p=a.2^n+1
    // because if x is not a quadratic residue then w=x^a is
    // a 2^n-th principal root of unity in GF_p
    T a1 = 2;
    T a2 = 5;
    T p1 = a1 * arith->exp(2, 15) + 1;
    T p2 = a2 * arith->exp(2, 15) + 1;
    // std::cerr << "p1=" << p1 << " p2=" << p2 << "\n";
    assert(arith->is_prime(p1));
    assert(arith->is_prime(p2));

    // ensure their product is bounded (b-1)^2*2^(n-1) < m
    T m = p1 * p2;
    // check overflow
    assert(m/p1 == p2);
    // std::cerr << " m=" << m << "\n";
    assert(arith->exp((b - 1), 2) * arith->exp(p, 2) < m);

    // find x so it is not a quadratic residue in GF_p1 and GF_p2
    assert(arith->jacobi(3, p1) == arith->jacobi(p1, 3));
    assert(arith->jacobi(p1, 3) == arith->jacobi(2, 3));
    assert(arith->jacobi(3, p2) == arith->jacobi(p2, 3));
    assert(arith->jacobi(p2, 3) == arith->jacobi(2, 3));
    assert(arith->jacobi(2, 3) == -1);
    // which means x=3 is not a quadratic residue in GF_p1 and GF_p2

    // therefore we can compute 2^n-th roots of unity in GF_p1 and GF_p2
    T w1 = arith->exp(3, a1);
    T w2 = arith->exp(3, a2);
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
    T w = arith->chinese_remainder(2, _a, _n);
    // std::cerr << " w=" << w << "\n";
    assert(w == 25559439);
  }

  void test_ext_gcd()
  {
    std::cout << "test_ext_gcd\n";
    SignedDoubleT<T> bezout[2];

    // not explicitely related to GF(97)
    assert(2 == arith->extended_gcd(240, 46, NULL, NULL));
    assert(6 == arith->extended_gcd(54, 24, NULL, NULL));
    assert(15 == arith->extended_gcd(210, 45, NULL, NULL));
    //
    assert(1 == arith->extended_gcd(97, 20, bezout, NULL));
    assert(bezout[0] == -7 && bezout[1] == 34);
  }

  void check_all_primes(std::vector<T> *primes, bool distinct) {
    assert(primes->size() > 0);

    typename std::vector<T>::size_type j;

    assert(arith->is_prime(primes->at(0)));

    for (j = 1; j != primes->size(); ++j) {
      // std::cout << j << ": " << primes->at(j) << "\n";
      assert(arith->is_prime(primes->at(j)));
      if (distinct)
        assert(primes->at(j-1) != primes->at(j));
    }
  }

  void check_primes_exponent(T nb, std::vector<T> *primes,
    std::vector<T> *exponent) {
    assert(primes->size() > 0);
    assert(exponent->size() > 0);

    typename std::vector<T>::size_type j;
    T y = 1;
    for (j = 0; j != primes->size(); ++j) {
      y *= arith->exp(primes->at(j), exponent->at(j));
    }
    assert(y == nb);
  }

  void test_factor_distinct_prime()
  {
    std::cout << "test_factor_distinct_prime\n";

    int i;
    for (i = 0; i < 1000; i++) {
      T x;
      std::vector<T> primes;
      // std::cout << "i=" << i << "\n";
      x = arith->weak_rand(max);
      // std::cout << "x=" << x << "\n";
      arith->factor_distinct_prime(x, &primes);
      check_all_primes(&primes, true);
    }
  }

  void test_factor_prime()
  {
    std::cout << "test_factor_prime\n";

    int i;
    std::vector<T> primes;
    std::vector<T> exponent;
    for (i = 0; i < 1000; i++) {
      T x;
      // std::cout << "i=" << i << "\n";
      x = arith->weak_rand(max);
      // std::cout << "x=" << x << "\n";
      arith->factor_prime(x, &primes, &exponent);
      check_all_primes(&primes, true);
      check_primes_exponent(x, &primes, &exponent);
      primes.clear();
      exponent.clear();
    }
  }

  void check_divisors(T nb, std::vector<T> *divisors, bool proper) {
    if (proper && arith->is_prime(nb))
      assert(divisors->size() == 0);
    else
      assert(divisors->size() > 0);

    typename std::vector<T>::size_type i;
    for (i = 0; i != divisors->size(); ++i) {
      // std::cout << i << ": " << divisors->at(i) << "\n";
      assert(nb % divisors->at(i) == 0);
      if (proper)
        assert(arith->is_prime(nb / divisors->at(i)));
    }
  }

  void test_get_proper_divisors()
  {
    std::cout << "test_get_proper_divisors\n";

    int i;
    std::vector<T> divisors;
    for (i = 0; i < 1000; i++) {
      T x;
      // std::cout << "i=" << i << "\n";
      x = arith->weak_rand(max);
      // std::cout << "x=" << x << "\n";
      arith->get_proper_divisors(x, &divisors);
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
      T x;
      // std::cout << "i=" << i << "\n";
      x = arith->weak_rand(max);
      std::vector<T> factors;
      arith->factor_distinct_prime(x, &factors);
      // std::cout << "x=" << x << "\n";
      arith->get_proper_divisors(x, &factors, &divisors);
      check_divisors(x, &divisors, true);
      divisors.clear();
    }
  }

  void test_compute_all_divisors()
  {
    std::cout << "test_compute_all_divisors\n";

    int i;
    std::vector<T> divisors;
    for (i = 0; i < 1000; i++) {
      T x;
      // std::cout << "i=" << i << "\n";
      x = arith->weak_rand(max);
      // std::cout << "x=" << x << "\n";
      arith->compute_all_divisors(x, &divisors);
      check_divisors(x, &divisors, false);
      divisors.clear();
    }
  }

  void test_get_code_len()
  {
    std::cout << "test_get_code_len\n";

    int i;
    for (i = 0; i < 1000; i++) {
      T order, n;
      // std::cout << "i=" << i << "\n";
      order = arith->weak_rand(max);
      // std::cout << "order=" << order << "\n";
      n = arith->weak_rand(order);
      // std::cout << "n=" << n << "\n";
      T len = arith->get_code_len(order, n);
      // std::cout << "len=" << len << "\n";
      assert(order % len == 0);
      assert(len >= n);
    }
  }

  void test_get_code_len_high_compo()
  {
    std::cout << "test_get_code_len_high_compo\n";

    int i;
    for (i = 0; i < 1000; i++) {
      T order, n;
      // std::cout << "i=" << i << "\n";
      order = arith->weak_rand(max);
      // std::cout << "order=" << order << "\n";
      n = arith->weak_rand(order);
      // std::cout << "n=" << n << "\n";
      T len = arith->get_code_len_high_compo(order, n);
      // std::cout << "len=" << len << "\n";
      assert(order % len == 0);
      assert(len >= n);
    }
  }

  void test_get_code_len_high_compo_2()
  {
    std::cout << "test_get_code_len_high_compo_2\n";

    int i;
    for (i = 0; i < 1000; i++) {
      T order, n;
      // std::cout << "i=" << i << "\n";
      order = arith->weak_rand(max);
      // std::cout << "order=" << order << "\n";
      n = arith->weak_rand(order);
      // std::cout << "n=" << n << "\n";
      std::vector<T> factors;
      arith->get_prime_factors(order, &factors);
      T len = arith->get_code_len_high_compo(&factors, n);
      // std::cout << "len=" << len << "\n";
      assert(order % len == 0);
      assert(len >= n);
    }
  }

  void check_prime_divisors(T nb, std::vector<T> *divisors, bool coprime) {
    assert(divisors->size() > 0);

    typename std::vector<T>::size_type i, j;
    for (i = 0; i != divisors->size(); ++i) {
      // std::cout << i << ": " << divisors->at(i) << "\n";
      assert(nb % divisors->at(i) == 0);
      if (!coprime)
        assert(arith->is_prime(divisors->at(i)));
      else {
        for (j = i+1; j != divisors->size(); ++j) {
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
      T n;
      // std::cout << "i=" << i << "\n";
      n = arith->weak_rand(max);
      // std::cout << "n=" << n << "\n";
      arith->get_coprime_factors(n, &divisors);
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
      T n;
      // std::cout << "i=" << i << "\n";
      n = arith->weak_rand(max);
      // std::cout << "n=" << n << "\n";
      arith->get_prime_factors(n, &divisors);
      check_prime_divisors(n, &divisors, false);
      divisors.clear();
    }
  }

  void arith_utest()
  {
    std::cout << "arith_utest with sizeof(T)=" << sizeof(T) << "\n";

    srand(time(0));

    test_basic_ops();
    test_reciprocal();
    test_chinese_remainder();
    test_jacobi();
    test_mul_bignum();
    test_ext_gcd();
    test_factor_distinct_prime();
    test_factor_prime();
    test_get_proper_divisors();
    test_compute_all_divisors();
    test_get_code_len();
    test_get_code_len_high_compo();
    test_get_coprime_factors();
    test_get_prime_factors();
  }

  void arith_utest_no256()
  {
    std::cout << "arith_utest_no256\n";

    srand(time(0));

    test_basic_ops();
    test_reciprocal();
    test_jacobi();
    test_factor_distinct_prime();
    test_factor_prime();
    test_get_proper_divisors();
    test_compute_all_divisors();
    test_get_code_len();
    test_get_code_len_high_compo();
    test_get_coprime_factors();
    test_get_prime_factors();
  }

  void arith_utest_mpz()
  {
    std::cout << "arith_utest_mpz\n";

    srand(time(0));

    test_basic_ops();
    test_chinese_remainder();
    test_mul_bignum();
    test_jacobi();
    test_ext_gcd();
  }
};

template class Arith<uint32_t>;
template class Arith<uint64_t>;
template class Arith<mpz_class>;

void arith_utest()
{
  ArithUtest<uint32_t> ArithUtest_uint32;
  ArithUtest_uint32.arith_utest();
  ArithUtest<uint64_t> ArithUtest_uint64;
  ArithUtest_uint64.arith_utest();
  ArithUtest<__uint128_t> ArithUtest_uint128;
  ArithUtest_uint128.arith_utest_no256();
  ArithUtest<mpz_class> ArithUtest_mpz;
  ArithUtest_mpz.arith_utest_mpz();
}
