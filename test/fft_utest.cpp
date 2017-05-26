
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
    SignedDoubleT<T> bezout[2];

    //not explicitely related to GF(97)
    assert(2 == gf->_extended_gcd(240, 46, NULL, NULL));
    assert(6 == gf->_extended_gcd(54, 24, NULL, NULL));
    assert(15 == gf->_extended_gcd(210, 45, NULL, NULL));
    //
    assert(1 == gf->_extended_gcd(97, 20, bezout, NULL));
    assert(bezout[0] == -7 && bezout[1] == 34);
    assert(gf->inv(20) == 34);
    //
    int i;
    for (i = 0;i < 100;i++) {
      T x = gf->weak_rand();
      assert(1 == gf->_extended_gcd(97, x, bezout, NULL));
      //std::cerr << bezout[0] << "*" << 97 << " " << bezout[1] << "*" << x << "=1\n";
      T y = gf->inv(x);
      //std::cerr << "inv(x)=" << y << "\n";
      if (bezout[1] < 0)
        bezout[1] = bezout[1] + gf->card();
      assert(bezout[1] == y);
    }
  }

  void test_gcd_gf97()
  {
    std::cout << "test_gcd_gf97\n";

    GFP<T> gf(97);
    test_gcd(&gf);
  }

  /** 
   * http://www.math.unm.edu/~loring/links/discrete_f05/remainder.pdf
   * http://gauss.math.luc.edu/greicius/Math201/Fall2012/Lectures/ChineseRemainderThm.article.pdf
   * 
   * @param gf 
   */
  void test_chinese_remainder1(GF<T> *gf)
  {
    T a[4];
    T n[4];
    T omega;

    a[0] = 4;
    n[0] = 107;
    a[1] = 2;
    n[1] = 74;
    omega = gf->_chinese_remainder(2, a, n);
    assert(omega == 5996);

    a[0] = 6;
    n[0] = 7;
    a[1] = 4;
    n[1] = 8;
    omega = gf->_chinese_remainder(2, a, n);
    assert(omega == 20);

    a[0] = 3;
    n[0] = 4;
    a[1] = 0;
    n[1] = 6;
    omega = gf->_chinese_remainder(2, a, n);
    //no solution
  }

  /** 
   * Example taken from Pierre Meunier's book
   * 
   * @param gf 
   */
  void test_chinese_remainder2(GF<T> *gf)
  {
    T p1 = 2 * gf->_exp(2, 15) + 1;
    T p2 = 5 * gf->_exp(2, 15) + 1;
    DoubleT<T> m = p1 * p2;
    T a[2];
    T n[2];
    a[0] = 9;
    n[0] = p1;
    a[1] = 243;
    n[1] = p2;
    T omega = gf->_chinese_remainder(2, a, n);
    std::cerr << "p1=" << p1 << " p2=" << p2 << " m=" << m << " omega=" << omega << "\n";
    assert(omega == 25559439);
  }

  void test_chinese_remainder_gf5()
  {
    std::cout << "test_chinese_remainder_gf5\n";

    GFP<T> gf5(5);
    test_chinese_remainder1(&gf5);
    test_chinese_remainder2(&gf5);
  }

  void fft_utest()
  {
    std::cout << "fft_utest\n";

    test_gcd_gf97();
    test_chinese_remainder_gf5();
  }
};

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;

//template class GF<uint64_t>;
//template class GFP<uint64_t>;
//template class GF2N<uint64_t>;

//template class GF<mpz_class>;
//template class GFP<mpz_class>;

void fft_utest()
{
  FFTUtest<uint32_t> fftutest_uint32;
  fftutest_uint32.fft_utest();
  //FFTUtest<uint64_t> fftutest_uint64;
  //fftutest_uint64.fft_utest();
  //FFTUtest<mpz_class> fftutest_mpz;
  //fftutest_mpz.fft_utest();
}
