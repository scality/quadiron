
#include "ntl.h"

template<typename T>
class GFUtest
{
 public:
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
    assert(gf->check_prime_root(gf->get_root()));
  }

  void test_get_order(GF<T> *gf)
  {
    int i;
    T x;
    T order;
    T h = gf->card_minus_one();
    for (i = 0; i < 1000; i++) {
      // std::cout << "i=" << i << "\n";
      // std::cout << gf->card() << "\n";
      x = gf->weak_rand();
      order = gf->get_order(x);
      // std::cout << "x=" << x << " " << order << std::endl;
      assert(gf->exp(x, order) == 1);
      assert(h % order == 0);
    }
  }

  void test_get_nth_root(GF<T> *gf)
  {
    int i;
    T x;
    T nth_root;
    T h = gf->card_minus_one();
    for (i = 0; i < 1000; i++) {
      // std::cout << "i=" << i << "\n";
      // std::cout << gf->card() << "\n";
      x = gf->weak_rand();
      nth_root = gf->get_nth_root(x);
      // std::cout << "x=" << x << " " << nth_root << std::endl;
      assert(gf->exp(nth_root, x) == 1);
    }
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

  void test_log_gf5()
  {
    std::cout << "test_log_gf5\n";

    GFP<T> gf5(5);
    test_log(&gf5);
  }

  void test_prime_root_gf5()
  {
    std::cout << "test_prime_root_gf5\n";

    GFP<T> gf5(5);
    test_find_prime_root(&gf5);
    test_get_order(&gf5);
    test_get_nth_root(&gf5);
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
    test_get_order(&gf256);
    test_get_nth_root(&gf256);
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
    test_get_order(&gf2n);
    test_get_nth_root(&gf2n);
  }

  void gf_utest()
  {
    std::cout << "gf_utest\n";

    srand(time(0));

    test_negation_gf5();
    test_reciprocal_gf5();
    test_log_gf5();
    test_prime_root_gf5();
    test_negation_gf256();
    test_reciprocal_gf256();
    test_log_gf256();
    test_negation_gf256();
    test_prime_root_gf256();
  }

  void gf_utest_2_bign(T n)
  {
    std::cout << "gf_utest 2^" << n << "\n";

    srand(time(0));

    test_negation_gf2_bign(n);
    test_reciprocal_gf2_bign(n);
    if (n < 128)  // it works slowly in GF2N(128)
      test_prime_root_gf2(n);
    // test_log_gf2_bign(n);
  }

  void gf_utest_nogf2n()
  {
    std::cout << "gf_utest_nogf2n\n";

    srand(time(0));

    test_negation_gf5();
    test_reciprocal_gf5();
    test_log_gf5();
    test_prime_root_gf5();
  }

  void gf_utest_2_n()
  {
    T max_n = 8 * sizeof(T);
    std::cout << "gf_utest_2_n for max_n=" << max_n << "\n";
    for (int i = 8; i <=max_n; i *= 2)
      gf_utest_2_bign(i);
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
  gfutest_uint32.gf_utest_2_n();
  GFUtest<uint64_t> gfutest_uint64;
  gfutest_uint64.gf_utest();
  gfutest_uint64.gf_utest_2_n();
  GFUtest<__uint128_t> gfutest_uint128;
  // gfutest_uint128.gf_utest(); // gfp(n) does not work for uint128
  gfutest_uint128.gf_utest_2_n();
  GFUtest<mpz_class> gfutest_mpz;
  gfutest_mpz.gf_utest_nogf2n();  // XXX gf2n broken for now
}
