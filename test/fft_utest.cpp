
#include "ntl.h"
#include "gf.cpp"
#include "gfp.cpp"
#include "gf2n.cpp"

template<typename T>
class FFTUtest
{
public:

  void test_gcd(GF<T> *gf)
  {
    T bezout[2];

    //not explicitely related to GF(97)
    assert(2 == gf->extended_gcd(240, 46, NULL, NULL));
    assert(6 == gf->extended_gcd(54, 24, NULL, NULL));
    assert(15 == gf->extended_gcd(210, 45, NULL, NULL));
    //
    assert(1 == gf->extended_gcd(97, 20, bezout, NULL));
    assert(bezout[0] == -7 && bezout[1] == 34);
    assert(gf->inv(20) == 34);
    //
    int i;
    for (i = 0;i < 100;i++) {
      T x = gf->weak_rand();
      assert(1 == gf->extended_gcd(97, x, bezout, NULL));
      std::cerr << bezout[0] << "*" << 97 << " " << bezout[1] << "*" << x << "=1\n";
      std::cerr << gf->inv(60) << "\n";
      //assert(bezout[1] == gf->inv(60));
    }
  }

  void test_gcd_gf97()
  {
    std::cout << "test_gcd_gf97\n";

    GFP<T> gf(97);
    test_gcd(&gf);
  }

  /** 
   * Example taken from Pierre Meunier's book
   * 
   * @param gf 
   */
  void test_chinese_remainder(GF<T> *gf)
  {
    T p1 = 2 * gf->exp(2, 15) + 1;
    T p2 = 5 * gf->exp(2, 15) + 1;
    T m = p1 * p2;
    T a[2];
    T n[2];
    a[0] = 9;
    a[1] = 243;
    n[0] = p1;
    n[1] = p2;
    T omega = gf->chinese_remainder(2, a, n);
    std::cerr << "p1=" << p1 << " p2=" << p2 << " m=" << m << " omega=" << omega << "\n";
    assert(omega == 25559439);
  }

  void test_chinese_remainder_gf5()
  {
    std::cout << "test_chinese_remainder_gf5\n";

    GFP<T> gf5(5);
    test_chinese_remainder(&gf5);
  }

  void fft_utest()
  {
    std::cout << "fft_utest\n";

    test_gcd_gf97();
    test_chinese_remainder_gf5();
  }
};

template class GF<uint64_t>;
template class GFP<uint64_t>;
template class GF2N<uint64_t>;

template class GF<mpz_class>;
template class GFP<mpz_class>;

void fft_utest()
{
  //FFTUtest<uint64_t> fftutest_uint64;
  //fftutest_uint64.fft_utest();
  //FFTUtest<mpz_class> fftutest_mpz;
  //fftutest_mpz.fft_utest();
}
