
#include "ntl.h"

template<typename T>
class FFTUtest
{
public:

  void test_gcd1(GF<T> *gf)
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

  void test_gcd()
  {
    std::cout << "test_gcd\n";

    GFP<T> gf(97);
    test_gcd1(&gf);
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

    GFP<T> gf5(5);

    T a[4];
    T n[4];
    T omega;

    a[0] = 4;
    n[0] = 107;
    a[1] = 2;
    n[1] = 74;
    omega = gf5._chinese_remainder(2, a, n);
    assert(omega == 5996);

    a[0] = 6;
    n[0] = 7;
    a[1] = 4;
    n[1] = 8;
    omega = gf5._chinese_remainder(2, a, n);
    assert(omega == 20);

    a[0] = 3;
    n[0] = 4;
    a[1] = 0;
    n[1] = 6;
    omega = gf5._chinese_remainder(2, a, n);
    //no solution XXX detect it
  }

  void test_quadratic_residues()
  {
    std::cout << "test_quadratic_residues\n";

    GFP<T> gf32(32);
    int i;
    for (i = 0;i < 32;i++) {
      assert(gf32.is_quadratic_residue(gf32.exp(i, 2)));
    }

    GFP<T> gf7(7);
    assert(gf7.is_quadratic_residue(2));
    assert(!gf7.is_quadratic_residue(5));

    GFP<T> gf8(8);
    assert(gf8.is_quadratic_residue(1));
    assert(!gf8.is_quadratic_residue(3));
  }

  void test_jacobi()
  {
    assert(__gf64._jacobi(1001, 9907) == -1);
    assert(__gf64._jacobi(19, 45) == 1);
    assert(__gf64._jacobi(8, 21) == -1);
    assert(__gf64._jacobi(5, 21) == 1);
    assert(__gf64._jacobi(47, 221) == -1);
    assert(__gf64._jacobi(2, 221) == -1);
  }
  
  /** 
   * convert a number into a vector of digits padded with zeros
   * 
   * @param num 
   * 
   * @return 
   */
  Vec<T> *_convert_string2vec(GF<T> *gf, int N, char num[])
  {
    int i;
    Vec<T> *vec = new Vec<T>(gf, N);
    int len = strlen(num);
    
    for (i = 0;i < len;i++) {
      vec->set(i, num[len - i - 1] - '0');
    }
    for (;i < N;i++) {  
      vec->set(i, 0);
    }
    
    return vec;
  }

  /** 
   * Schönhage-Strassen algorithm
   * Example taken from Pierre Meunier's book
   * 
   * @param gf 
   */
  void test_mul_bignum()
  {
    std::cout << "test_mul_bignum\n";

    int b = 10; //base
    int p = 14; //we could multiply integers of 2^p digits
    int max_digits = __gf64._exp(2, p);
    //std::cerr << "p=" << p << " max_digits=" << max_digits << "\n";

    uint64_t n = p + 1;
    uint64_t N = __gf64._exp(2, n);
    //std::cerr << "n=" << n << " N=" << N << "\n";

    //choose 2 prime numbers of the form p=a.2^n+1
    //because if x is not a quadratic residue then w=x^a is
    //a 2^n-th principal root of unity in GF_p
    uint64_t a1 = 2;
    uint64_t a2 = 5;
    uint64_t p1 = a1 * __gf64._exp(2, 15) + 1;
    uint64_t p2 = a2 * __gf64._exp(2, 15) + 1;
    //std::cerr << "p1=" << p1 << " p2=" << p2 << "\n";
    assert(__gf64._is_prime(p1));
    assert(__gf64._is_prime(p2));
    
    //ensure their product is bounded (b-1)^2*2^(n-1) < m
    uint64_t m = p1 * p2;
    //check overflow
    assert(m/p1 == p2);
    //std::cerr << " m=" << m << "\n";
    assert(__gf64._exp((b - 1), 2) * __gf64._exp(p, 2) < m);

    //find x so it is not a quadratic residue in GF_p1 and GF_p2
    assert(__gf64._jacobi(3, p1) == __gf64._jacobi(p1, 3));
    assert(__gf64._jacobi(p1, 3) == __gf64._jacobi(2, 3));
    assert(__gf64._jacobi(3, p2) == __gf64._jacobi(p2, 3));
    assert(__gf64._jacobi(p2, 3) == __gf64._jacobi(2, 3));
    assert(__gf64._jacobi(2, 3) == -1);
    //which means x=3 is not a quadratic residue in GF_p1 and GF_p2

    //therefore we can compute 2^n-th roots of unity in GF_p1 and GF_p2
    uint64_t w1 = __gf64._exp(3, a1);
    uint64_t w2 = __gf64._exp(3, a2);
    //std::cerr << "w1=" << w1 << " w2=" << w2 << "\n";
    assert(w1 == 9);
    assert(w2 == 243);

    //find root of unity in GF_p1p2
    uint64_t _a[2];
    uint64_t _n[2];
    _a[0] = w1;
    _n[0] = p1;
    _a[1] = w2;
    _n[1] = p2;
    uint64_t w = __gf64._chinese_remainder(2, _a, _n);
    //std::cerr << " w=" << w << "\n";
    assert(w == 25559439);

    GFP<T> gf_m(m);
    FFT<T> fft(&gf_m, n, 25559439);

    //parse the big numbers
    char X[] = "1236548787985654354598651354984132468";
    char Y[] = "745211515185321545554545854598651354984132468";

    Vec<T> *_X = _convert_string2vec(&gf_m, N, X);
    //_X->dump();
    Vec<T> *_Y = _convert_string2vec(&gf_m, N, Y);
    //_Y->dump();

    Vec<T> *sfX = new Vec<T>(&gf_m, N);
    Vec<T> *sfY = new Vec<T>(&gf_m, N);
    Vec<T> *_XY = new Vec<T>(&gf_m, N);
    Vec<T> *sfXY = new Vec<T>(&gf_m, N);

    fft.fft(sfX, _X);
    fft.fft(sfY, _Y);

    for (int i = 0;i <= N-1;i++) {
      DoubleT<T> val = DoubleT<T>(sfX->get(i)) * sfY->get(i);
      _XY->set(i, val % m);
    }

    fft.ifft(sfXY, _XY);

    //carry propagation
    mpz_class z = 0;
    for (int i = 0;i <= N-1;i++) {
      mpz_class t, b;
      b = 10;
      mpz_pow_ui(t.get_mpz_t(), b.get_mpz_t(), i);
      z += sfXY->get(i) * t;
    }

    //std::cout << z << "\n";
    assert(z.get_str() == "921490395895362412399910100421159322712298564831565484737491129935640058571771024");

    delete sfXY;
    delete _XY;
    delete sfX;
    delete sfY;
    delete _X;
    delete _Y;
  }

  void test_fft()
  {
    u_int n;
    u_int N;
    u_int r;
    u_int q = 65537;
    GFP<T> gf = GFP<T>(q);
    u_int R = 3; //primitive root
    u_int n_data = 3;
    u_int n_parities = 3;

    std::cout << "test_fft\n";

    assert(gf._jacobi(R, q) == -1);

    //with this encoder we cannot exactly satisfy users request, we need to pad
    n = __gf64._log2(n_data + n_parities) + 1;
    N = __gf64._exp(2, n);

    //compute root of order N-1 such as r^(N-1) mod q == 1
    mpz_class _r = __gfmpz._exp(R, __gfmpz._exp(2, 16-n)) % gf.p;
    r = _r.get_ui();

    //std::cerr << "n=" << n << "\n";
    //std::cerr << "N=" << N << "\n";
    //std::cerr << "r=" << r << "\n";

    FFT<T> fft = FFT<T>(&gf, n, r);
    
    for (int j = 0;j < 100000;j++) {
      Vec<T> v(&gf, N), _v(&gf, N), v2(&gf, N);
      v.zero_fill();
      for (int i = 0;i < n_data;i++)
        v.set(i, gf.weak_rand());
      //v.dump();
      fft.fft(&_v, &v);
      //_v.dump();
      fft.ifft(&v2, &_v);
      //v2.dump();
      assert(v.eq(&v2));
    }
  }

  void fft_utest_no_mul_bignum()
  {
    std::cout << "fft_utest\n";

    test_gcd();
    test_chinese_remainder();
    test_quadratic_residues();
    test_jacobi();
  }

  void fft_utest()
  {
    std::cout << "fft_utest\n";

    test_gcd();
    test_chinese_remainder();
    test_quadratic_residues();
    test_jacobi();
    test_fft();
    test_mul_bignum();
  }
};

template class Mat<uint32_t>;
template class Vec<uint32_t>;
template class FFT<uint32_t>;

template class Mat<uint64_t>;
template class Vec<uint64_t>;
template class FFT<uint64_t>;

template class Mat<mpz_class>;
template class Vec<mpz_class>;

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;

template class GF<uint64_t>;
template class GFP<uint64_t>;
template class GF2N<uint64_t>;

template class GF<mpz_class>;
template class GFP<mpz_class>;

void fft_utest()
{
  FFTUtest<uint32_t> fftutest_uint32;
  fftutest_uint32.fft_utest_no_mul_bignum();
  FFTUtest<uint64_t> fftutest_uint64;
  fftutest_uint64.fft_utest();
  FFTUtest<mpz_class> fftutest_mpz;
  fftutest_mpz.fft_utest_no_mul_bignum(); //too slow
}
