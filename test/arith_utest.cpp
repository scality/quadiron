
#include "ntl.h"

template<typename T>
class ArithUtest
{
 public:
  T max;
  Arith<T> *arith;
  ArithUtest() {
    max = std::numeric_limits<T>::max();
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
    std::vector<T> *primes = new std::vector<T>();
    for (i = 0; i < 1000; i++) {
      T x;
      // std::cout << "i=" << i << "\n";
      x = arith->weak_rand(max);
      // std::cout << "x=" << x << "\n";
      arith->factor_distinct_prime(x, primes);
      check_all_primes(primes, true);
      primes->clear();
    }
  }

  void test_factor_prime()
  {
    std::cout << "test_factor_prime\n";

    int i;
    std::vector<T> *primes = new std::vector<T>();
    std::vector<T> *exponent = new std::vector<T>();
    for (i = 0; i < 1000; i++) {
      T x;
      // std::cout << "i=" << i << "\n";
      x = arith->weak_rand(max);
      // std::cout << "x=" << x << "\n";
      arith->factor_prime(x, primes, exponent);
      check_all_primes(primes, true);
      check_primes_exponent(x, primes, exponent);
      primes->clear();
      exponent->clear();
    }

    delete primes;
    delete exponent;
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
    std::vector<T> *divisors = new std::vector<T>();
    for (i = 0; i < 1000; i++) {
      T x;
      // std::cout << "i=" << i << "\n";
      x = arith->weak_rand(max);
      // std::cout << "x=" << x << "\n";
      arith->get_proper_divisors(x, divisors);
      check_divisors(x, divisors, true);
      divisors->clear();
    }
    delete divisors;
  }

  void test_compute_all_divisors()
  {
    std::cout << "test_compute_all_divisors\n";

    int i;
    std::vector<T> *divisors = new std::vector<T>();
    for (i = 0; i < 1000; i++) {
      T x;
      // std::cout << "i=" << i << "\n";
      x = arith->weak_rand(max);
      // std::cout << "x=" << x << "\n";
      arith->compute_all_divisors(x, divisors);
      check_divisors(x, divisors, false);
      divisors->clear();
    }
    delete divisors;
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

    std::vector<T> *divisors = new std::vector<T>();
    int i;
    for (i = 0; i < 1000; i++) {
      T n;
      // std::cout << "i=" << i << "\n";
      n = arith->weak_rand(max);
      // std::cout << "n=" << n << "\n";
      arith->get_coprime_factors(n, divisors);
      check_prime_divisors(n, divisors, true);
      divisors->clear();
    }
    delete divisors;
  }

  void test_get_prime_factors()
  {
    std::cout << "test_get_prime_factors\n";

    std::vector<T> *divisors = new std::vector<T>();
    int i;
    for (i = 0; i < 1000; i++) {
      T n;
      // std::cout << "i=" << i << "\n";
      n = arith->weak_rand(max);
      // std::cout << "n=" << n << "\n";
      arith->get_prime_factors(n, divisors);
      check_prime_divisors(n, divisors, false);
      divisors->clear();
    }
    delete divisors;
  }

  void arith_utest()
  {
    std::cout << "arith_utest\n";

    srand(time(0));

    test_basic_ops();
    test_reciprocal();
    test_chinese_remainder();
    test_jacobi();
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
