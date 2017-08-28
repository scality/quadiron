
#include "ntl.h"

template<typename T>
class GFUtest
{
 public:
  void test_basic_ops()
  {
    assert(__gf64._sqrt(2025) == 45);
    assert(__gf64._is_prime(2));
    assert(__gf64._is_prime(3));
    assert(__gf64._is_prime(13));
    assert(!__gf64._is_prime(4));
    assert(!__gf64._is_prime(15));
  }

  void test_negation(GF<T> *gf)
  {
    int i;

    for (i = 0; i < 100; i++) {
      T x, y;

      // std::cout << "i=" << i << "\n";

      x = gf->weak_rand();
      // std::cout << "x=" << x << "\n";
      y = gf->neg(x);
      // std::cout << "inv(x)=" << y << "\n";
      assert(gf->add(x, y) == 0);
    }
  }

  void test_reciprocal(GF<T> *gf)
  {
    int i;

    for (i = 0; i < 100; i++) {
      T x, y;

      // std::cout << "i=" << i << "\n";

      x = gf->weak_rand();
      // std::cout << "x=" << x << "\n";
      y = gf->inv(x);
      // std::cout << "inv(x)=" << y << "\n";
      assert(gf->mul(x, y) == 1);
    }
  }

  void test_log(GF<T> *gf)
  {
    int i;

    for (i = 0; i < 1000; i++) {
      T x, y, z, t;

      // std::cout << "i=" << i << "\n";
      // std::cout << gf->card() << "\n";
      x = gf->weak_rand();
      y = gf->weak_rand();
      // std::cout << "x=" << x << "\n";
      // std::cout << "y=" << y << "\n";
      try {
        z = gf->log(x, y);
        // std::cout << "z=" << z << "\n";
      } catch (...) {
        // std::cout << "not found\n";
        continue;
      }
      t = gf->exp(x, z);
      // std::cout << "t=" << t << "\n";
      assert(t == y);
    }
  }

  void test_find_prime_root(GF<T> *gf)
  {
    gf->find_prime_root();
    // std::cout << "root " << gf->root << std::endl;
    assert(gf->check_prime_root(gf->root));
  }

  void test_negation_gf5()
  {
    std::cout << "test_negation_gf5\n";

    GFP<T> gf5(5);
    test_negation(&gf5);
  }

  void test_reciprocal_gf5()
  {
    std::cout << "test_reciprocal_gf5\n";

    GFP<T> gf5(5);
    test_reciprocal(&gf5);
  }

  void test_prime_root_gf5()
  {
    std::cout << "test_prime_root_gf5\n";

    GFP<T> gf5(5);
    test_find_prime_root(&gf5);
  }

  void test_log_gf()
  {
    std::cout << "test_log_gf\n";

    GFP<T> gf(32);
    test_log(&gf);
  }

  void test_prime_root_gf32()
  {
    std::cout << "test_prime_root_gf32\n";

    GFP<T> gf(32);
    test_find_prime_root(&gf);
  }

  void test_negation_gf256()
  {
    std::cout << "test_negation_gf256\n";

    GF2N<T> gf256(8);
    test_negation(&gf256);
  }

  void test_reciprocal_gf256()
  {
    std::cout << "test_reciprocal_gf256\n";

    GF2N<T> gf256(8);
    test_reciprocal(&gf256);
  }

  void test_log_gf256()
  {
    std::cout << "test_log_gf256\n";

    GF2N<T> gf256(8);
    test_log(&gf256);
  }

  void test_prime_root_gf256()
  {
    std::cout << "test_prime_root_gf256\n";

    GF2N<T> gf256(8);
    test_find_prime_root(&gf256);
  }

  void test_misc_gf29()
  {
    GFP<T> gf29(29);

    std::cout << "test_misc_gf29\n";

    assert(gf29.neg(28) == 1);
    assert(gf29.neg(18) == 11);
    assert(gf29.neg(4) == 25);

    assert(gf29.add(28, 1) == 0);
  }

  void test_negation_gf2_bign(T n)
  {
    std::cout << "test_negation_gf(2^" << n << ")\n";

    GF2N<T> gf2n(n);
    test_negation(&gf2n);
  }

  void test_reciprocal_gf2_bign(T n)
  {
    std::cout << "test_reciprocal_gf(2^" << n << ")\n";

    GF2N<T> gf2n(n);
    test_reciprocal(&gf2n);
  }

  void test_log_gf2_bign(T n)
  {
    std::cout << "test_log_gf(2^" << n << ")\n";

    GF2N<T> gf2n(n);
    test_log(&gf2n);
  }

  void test_prime_root_gf2(T n)
  {
    std::cout << "test_prime_root_gf(2^" << n << ")\n";

    GF2N<T> gf2n(n);
    test_find_prime_root(&gf2n);
  }

  void gf_utest()
  {
    std::cout << "gf_utest\n";

    srand(time(0));

    test_basic_ops();
    test_negation_gf5();
    test_reciprocal_gf5();
    test_prime_root_gf5();
    test_log_gf();
    // test_prime_root_gf32();
    test_negation_gf256();
    test_reciprocal_gf256();
    test_log_gf256();
    test_misc_gf29();
    test_negation_gf256();
    test_prime_root_gf256();
  }

  void gf_utest_2_bign(T n)
  {
    std::cout << "gf_utest 2^" << n << ")\n";

    srand(time(0));

    test_negation_gf2_bign(n);
    test_reciprocal_gf2_bign(n);
    // test_prime_root_gf2(n);
    // test_log_gf2_bign(n);
  }

  void gf_utest_nogf2n()
  {
    std::cout << "gf_utest\n";

    srand(time(0));

    test_basic_ops();
    test_negation_gf5();
    test_reciprocal_gf5();
    test_log_gf();
    test_misc_gf29();
  }
};

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;

template class GF<uint64_t>;
template class GFP<uint64_t>;
template class GF2N<uint64_t>;

template class GF<mpz_class>;
template class GFP<mpz_class>;

void gf_utest()
{
  GFUtest<uint32_t> gfutest_uint32;
  gfutest_uint32.gf_utest();
  gfutest_uint32.test_prime_root_gf2(8);
  gfutest_uint32.test_prime_root_gf2(16);
  gfutest_uint32.gf_utest_2_bign(32);
  GFUtest<uint64_t> gfutest_uint64;
  gfutest_uint64.gf_utest();
  gfutest_uint64.gf_utest_2_bign(32);
  gfutest_uint64.gf_utest_2_bign(64);
  gfutest_uint64.test_prime_root_gf2(8);
  gfutest_uint64.test_prime_root_gf2(16);
  gfutest_uint64.test_prime_root_gf2(32);
  GFUtest<__uint128_t> gfutest_uint128;
  // gfutest_uint128.gf_utest(); // gfp(n) does not work for uint128
  gfutest_uint128.gf_utest_2_bign(32);
  gfutest_uint128.gf_utest_2_bign(64);
  gfutest_uint128.gf_utest_2_bign(128);
  GFUtest<mpz_class> gfutest_mpz;
  gfutest_mpz.gf_utest_nogf2n();  // XXX gf2n broken for now
}
