
#include "ntl.h"
#include "gf.cpp"
#include "gfp.cpp"
#include "gf2n.cpp"

template<typename T>
class GFUtest
{
public:

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
      assert(gf->add(x, y) == gf->zero());
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
      assert(gf->mul(x, y) == gf->one());
    }
  }

  void test_log(GF<T> *gf)
  {
    int i;
    
    for (i = 0; i < 100;i++) {
      T x, y;
      
      //std::cout << "i=" << i << "\n";
      x = gf->weak_rand();
      //std::cout << "x=" << x << "\n";
      try {
        y = gf->log(x);
      } catch (...) {
        //std::cout << "not found\n";
        continue ;
      }
      //std::cout << "log(x)=" << y << "\n";
      assert(gf->pow(y) == x);
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
  
  void gf_utest()
  {
    std::cout << "gf_utest\n";
    
    srand(time(0));

    test_negation_gf5();
    test_reciprocal_gf5();
    test_log_gf5();
    test_negation_gf256();
    test_reciprocal_gf256();
    test_log_gf256();
    test_misc_gf29();
  }

  void gf_utest_nogf2n()
  {
    std::cout << "gf_utest\n";
    
    srand(time(0));

    test_negation_gf5();
    test_reciprocal_gf5();
    test_log_gf5();
    test_misc_gf29();
  }

};

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;

template class GF<mpz_class>;
template class GFP<mpz_class>;

void gf_utest()
{
  GFUtest<uint32_t> gfutest_uint32;
  gfutest_uint32.gf_utest();
  GFUtest<mpz_class> gfutest_mpz;
  gfutest_mpz.gf_utest_nogf2n(); //XXX gf2n broken for now
}
