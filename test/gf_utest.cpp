
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
    
    for (i = 0; i < 100;i++) {
      T x, y;
      
      //std::cout << "i=" << i << "\n";
      
      x = gf->weak_rand();
      //std::cout << "x=" << x << "\n";
      y = gf->neg(x);
      //std::cout << "inv(x)=" << y << "\n";
      assert(gf->add(x, y) == 0);
    }
  }

  void test_reciprocal(GF<T> *gf)
  {
    int i;
    
    for (i = 0; i < 100;i++) {
      T x, y;
      
      //std::cout << "i=" << i << "\n";
      
      x = gf->weak_rand();
      //std::cout << "x=" << x << "\n";
      y = gf->inv(x);
      //std::cout << "inv(x)=" << y << "\n";
      assert(gf->mul(x, y) == 1);
    }
  }

  void test_log(GF<T> *gf)
  {
    int i;
    
    for (i = 0; i < 1000;i++) {
      T x, y, z, t;
      
      //std::cout << "i=" << i << "\n";
      //std::cout << gf->card() << "\n";
      x = gf->weak_rand();
      y = gf->weak_rand();
      //std::cout << "x=" << x << "\n";
      //std::cout << "y=" << y << "\n";
      try {
        z = gf->log(x, y);
        //std::cout << "z=" << z << "\n";
      } catch (...) {
        //std::cout << "not found\n";
        continue ;
      }
      t = gf->exp(x, z);
      //std::cout << "t=" << t << "\n";
      assert(t == y);
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

  void test_log_gf()
  {
    std::cout << "test_log_gf\n";

    GFP<T> gf(32);
    test_log(&gf);
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

  void test_misc_gf29()
  {
    GFP<T> gf29(29);
    
    std::cout << "test_misc_gf29\n";
    
    assert(gf29.neg(28) == 1);
    assert(gf29.neg(18) == 11);
    assert(gf29.neg(4) == 25);
    
    assert(gf29.add(28, 1) == 0);
  }

  void test_negation_gf2_32()
  {
    std::cout << "test_negation_gf(2^32)\n";

    GF2N<T> gf2_32(32);
    test_negation(&gf2_32);
  }

  void test_reciprocal_gf2_32()
  {
    std::cout << "test_reciprocal_gf(2^32\n";

    GF2N<T> gf2_32(32);
    test_reciprocal(&gf2_32);
  }

  void test_log_gf2_32()
  {
    std::cout << "test_log_gf(2^32\n";

    GF2N<T> gf2_32(32);
    test_log(&gf2_32);
  }

  void gf_utest()
  {
    std::cout << "gf_utest\n";
    
    srand(time(0));

    test_basic_ops();
    test_negation_gf5();
    test_reciprocal_gf5();
    test_log_gf();
    test_negation_gf256();
    test_reciprocal_gf256();
    test_log_gf256();
    test_misc_gf29();
  }

  void gf_utest_2_32()
  {
    std::cout << "gf_utest\n";

    srand(time(0));

    test_negation_gf2_32();
    test_reciprocal_gf2_32();
    // test_log_gf2_32();
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
  GFUtest<uint64_t> gfutest_uint64;
  gfutest_uint64.gf_utest();
  gfutest_uint64.gf_utest_2_32();
  GFUtest<mpz_class> gfutest_mpz;
  gfutest_mpz.gf_utest_nogf2n(); //XXX gf2n broken for now
}
