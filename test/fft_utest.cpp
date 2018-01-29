
#include "ntl.h"

template<typename T>
class FFTUtest
{
 public:
  void test_gcd1(GF<T> *gf)
  {
    SignedDoubleT<T> bezout[2];

    assert(gf->inv(20) == 34);
    //
    int i;
    for (i = 0; i < 100; i++) {
      T x = gf->weak_rand();
      assert(1 == _extended_gcd<T>(97, x, bezout, NULL));
      // std::cerr << bezout[0] << "*" << 97 << " " << bezout[1] << "*";
      // stdd:cerr << x << "=1\n";
      T y = gf->inv(x);
      // std::cerr << "inv(x)=" << y << "\n";
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

  void test_quadratic_residues()
  {
    std::cout << "test_quadratic_residues\n";

    GFP<T> gf32(32);
    int i;
    for (i = 0; i < 32; i++) {
      assert(gf32.is_quadratic_residue(gf32.exp(i, 2)));
    }

    GFP<T> gf7(7);
    assert(gf7.is_quadratic_residue(2));
    assert(!gf7.is_quadratic_residue(5));

    GFP<T> gf8(8);
    assert(gf8.is_quadratic_residue(1));
    assert(!gf8.is_quadratic_residue(3));
  }

  /**
   * convert a number into a vector of digits padded with zeros
   *
   * @param num
   *
   * @return
   */
  Vec<T> *_convert_string2vec(GF<T> *gf, int n, char num[])
  {
    int i;
    Vec<T> *vec = new Vec<T>(gf, n);
    int len = strlen(num);

    for (i = 0; i < len; i++) {
      vec->set(i, num[len - i - 1] - '0');
    }
    for (; i < n; i++) {
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
    if (sizeof(T) < 8)
      return;
    std::cout << "test_mul_bignum\n";

    int b = 10;  // base
    int p = 14;  // we could multiply integers of 2^p digits
    // int max_digits = _exp<T>(2, p);
    // std::cerr << "p=" << p << " max_digits=" << max_digits << "\n";

    uint64_t l = p + 1;
    // std::cerr << "l=" << l << "\n";

    // choose 2 prime numbers of the form p=a.2^n+1
    // because if x is not a quadratic residue then w=x^a is
    // a 2^n-th principal root of unity in GF_p
    uint64_t a1 = 2;
    uint64_t a2 = 5;
    uint64_t p1 = a1 * _exp<T>(2, 15) + 1;
    uint64_t p2 = a2 * _exp<T>(2, 15) + 1;
    // std::cerr << "p1=" << p1 << " p2=" << p2 << "\n";
    assert(_is_prime<T>(p1));
    assert(_is_prime<T>(p2));

    // ensure their product is bounded (b-1)^2*2^(n-1) < m
    uint64_t m = p1 * p2;
    // check overflow
    assert(m/p1 == p2);
    // std::cerr << " m=" << m << "\n";
    assert(_exp<T>((b - 1), 2) * _exp<T>(p, 2) < m);

    // find x so it is not a quadratic residue in GF_p1 and GF_p2
    assert(_jacobi<T>(3, p1) == _jacobi<T>(p1, 3));
    assert(_jacobi<T>(p1, 3) == _jacobi<T>(2, 3));
    assert(_jacobi<T>(3, p2) == _jacobi<T>(p2, 3));
    assert(_jacobi<T>(p2, 3) == _jacobi<T>(2, 3));
    assert(_jacobi<T>(2, 3) == -1);
    // which means x=3 is not a quadratic residue in GF_p1 and GF_p2

    // therefore we can compute 2^n-th roots of unity in GF_p1 and GF_p2
    uint64_t w1 = _exp<T>(3, a1);
    uint64_t w2 = _exp<T>(3, a2);
    // std::cerr << "w1=" << w1 << " w2=" << w2 << "\n";
    assert(w1 == 9);
    assert(w2 == 243);

    // find root of unity in GF_p1p2
    uint64_t _a[2];
    uint64_t _n[2];
    _a[0] = w1;
    _n[0] = p1;
    _a[1] = w2;
    _n[1] = p2;
    uint64_t w = _chinese_remainder<uint64_t>(2, _a, _n);
    // std::cerr << " w=" << w << "\n";
    assert(w == 25559439);

    GFP<T> gf_m(m);
    FFTLN<T> fft(&gf_m, l, 25559439);

    // parse the big numbers
    char X[] = "1236548787985654354598651354984132468";
    char Y[] = "745211515185321545554545854598651354984132468";

    Vec<T> *_X = _convert_string2vec(&gf_m, fft.get_n(), X);
    // _X->dump();
    Vec<T> *_Y = _convert_string2vec(&gf_m, fft.get_n(), Y);
    // _Y->dump();

    Vec<T> *sfX = new Vec<T>(&gf_m, fft.get_n());
    Vec<T> *sfY = new Vec<T>(&gf_m, fft.get_n());
    Vec<T> *_XY = new Vec<T>(&gf_m, fft.get_n());
    Vec<T> *sfXY = new Vec<T>(&gf_m, fft.get_n());

    fft.fft(sfX, _X);
    fft.fft(sfY, _Y);

    for (int i = 0; i <= fft.get_n()-1; i++) {
      DoubleT<T> val = DoubleT<T>(sfX->get(i)) * sfY->get(i);
      _XY->set(i, val % m);
    }

    fft.ifft(sfXY, _XY);

    // carry propagation
    mpz_class z = 0;
    mpz_class u;
    for (int i = 0; i <= fft.get_n()-1; i++) {
      mpz_class t, b;
      b = 10;
      mpz_pow_ui(t.get_mpz_t(), b.get_mpz_t(), i);
      u = mpz_class(std::to_string(sfXY->get(i)));
      z += u * t;
    }

    // std::cout << z << "\n";
    assert(z.get_str() == std::string("921490395895362412399910100421159322") +
      "712298564831565484737491129935640058571771024");

    delete sfXY;
    delete _XY;
    delete sfX;
    delete sfY;
    delete _X;
    delete _Y;
  }

  void test_fftn()
  {
    u_int n;
    u_int r;
    T q = 65537;
    GFP<T> gf = GFP<T>(q);
    u_int R = gf.get_prime_root();  // primitive root
    u_int n_data = 3;
    u_int n_parities = 3;

    std::cout << "test_fftn\n";

    assert(_jacobi<T>(R, q) == -1);

    // with this encoder we cannot exactly satisfy users request, we need to pad
    // n = minimal divisor of (q-1) that is at least (n_parities + n_data)
    n = gf.get_code_len(n_parities + n_data);

    // compute root of order n-1 such as r^(n-1) mod q == 1
    r = gf.get_nth_root(n);

    // std::cerr << "n=" << n << "\n";
    // std::cerr << "r=" << r << "\n";

    DFTN<T> fft = DFTN<T>(&gf, n, r);

    Vec<T> v(&gf, fft.get_n()), _v(&gf, fft.get_n()), v2(&gf, fft.get_n());
    for (int j = 0; j < 100000; j++) {
      v.zero_fill();
      for (int i = 0; i < n_data; i++)
        v.set(i, gf.weak_rand());
      // v.dump();
      fft.fft(&_v, &v);
      // _v.dump();
      fft.ifft(&v2, &_v);
      // v2.dump();
      assert(v.eq(&v2));
    }
  }

  void test_fft2k()
  {
    std::cout << "test_fft2k\n";
    test_fft2k_vec();
    test_fft2k_vecp();
  }

  void test_fft2k_vec()
  {
    u_int n;
    u_int q = 65537;
    GFP<T> gf = GFP<T>(q);
    u_int R = gf.get_prime_root();  // primitive root
    u_int n_data = 3;
    u_int n_parities = 3;

    std::cout << "test_fft2k_vec\n";

    assert(_jacobi<T>(R, q) == -1);

    // with this encoder we cannot exactly satisfy users request, we need to pad
    // n = minimal divisor of (q-1) that is at least (n_parities + n_data)
    n = gf.get_code_len(n_parities + n_data);

    // std::cerr << "n=" << n << "\n";

    FFT2K<T> fft = FFT2K<T>(&gf, n);

    Vec<T> v(&gf, fft.get_n()), _v(&gf, fft.get_n()), v2(&gf, fft.get_n());
    for (int j = 0; j < 100000; j++) {
      v.zero_fill();
      for (int i = 0; i < n_data; i++)
        v.set(i, gf.weak_rand());
      // v.dump();
      fft.fft(&_v, &v);
      // _v.dump();
      fft.ifft(&v2, &_v);
      // v2.dump();
      assert(v.eq(&v2));
    }
  }

  void test_fft2k_vecp()
  {
    u_int n;
    u_int q = 65537;
    GFP<T> gf = GFP<T>(q);
    u_int R = gf.get_prime_root();  // primitive root
    u_int n_data = 3;
    u_int n_parities = 3;
    size_t size = 4;

    std::cout << "test_fft2k_vecp\n";

    assert(_jacobi<T>(R, q) == -1);

    // with this encoder we cannot exactly satisfy users request, we need to pad
    // n = minimal divisor of (q-1) that is at least (n_parities + n_data)
    n = gf.get_code_len(n_parities + n_data);

    // std::cerr << "n=" << n << "\n";

    FFT2K<T> fft = FFT2K<T>(&gf, n);

    int vec_n = fft.get_n();

    Vecp<T> v(n_data, size), v2(vec_n, size), _v2(vec_n, size);
    VVecp<T> _v(&v, vec_n);
    for (int j = 0; j < 100000; j++) {
      for (int i = 0; i < n_data; i++) {
        T* mem = v.get(i);
        for (int u = 0; u < size; u++) {
          mem[u] = gf.weak_rand();
        }
      }
      // _v.dump();
      fft.fft(&v2, &_v);
      // v2.dump();
      fft.ifft(&_v2, &v2);
      // _v2.dump();
      assert(_v.eq(&_v2));
    }
  }

  void test_fftpf()
  {
    T n;
    GF2N<T> gf = GF2N<T>(4);
    T n_data = 3;
    T n_parities = 3;

    std::cout << "test_fftpf\n";

    // with this encoder we cannot exactly satisfy users request, we need to pad
    // n = minimal divisor of (q-1) that is at least (n_parities + n_data)
    n = gf.get_code_len(n_parities + n_data);

    // std::cerr << "n=" << n << "\n";

    FFTPF<T> fft = FFTPF<T>(&gf, n);

    Vec<T> v(&gf, fft.get_n()), _v(&gf, fft.get_n()), v2(&gf, fft.get_n());
    for (int j = 0; j < 10000; j++) {
      v.zero_fill();
      for (int i = 0; i < n_data; i++)
        v.set(i, gf.weak_rand());
        // v.dump();
      fft.fft(&_v, &v);
        // _v.dump();
      fft.ifft(&v2, &_v);
        // v2.dump();
      assert(v.eq(&v2));
    }
  }

  void test_fftct_gfp()
  {
    T n;
    T q = 65537;
    GFP<T> gf = GFP<T>(q);
    T n_data = 3;
    T n_parities = 3;

    std::cout << "test_fftct_gfp\n";

    // with this encoder we cannot exactly satisfy users request, we need to pad
    // n = minimal divisor of (q-1) that is at least (n_parities + n_data)
    n = gf.get_code_len(n_parities + n_data);

    // std::cerr << "n=" << n << "\n";

    FFTCT<T> fft = FFTCT<T>(&gf, n);

    Vec<T> v(&gf, fft.get_n()), _v(&gf, fft.get_n()), v2(&gf, fft.get_n());
    for (int j = 0; j < 10000; j++) {
      v.zero_fill();
      for (int i = 0; i < n_data; i++)
        v.set(i, gf.weak_rand());
        // v.dump();
      fft.fft(&_v, &v);
        // _v.dump();
      fft.ifft(&v2, &_v);
        // v2.dump();
      assert(v.eq(&v2));
    }
  }

  void test_fftct_gf2n()
  {
    T n;
    T n_data = 3;
    T n_parities = 3;
    for (int gf_n = 4; gf_n <= 128 && gf_n <= 8 * sizeof(T); gf_n *= 2) {
      GF2N<T> gf = GF2N<T>(gf_n);

      std::cout << "test_fftct_gf2n=" << gf_n << "\n";

      // with this encoder we cannot exactly satisfy users request, we need to pad
      // n = minimal divisor of (q-1) that is at least (n_parities + n_data)
      n = gf.get_code_len(n_parities + n_data);

      // std::cerr << "n=" << n << "\n";

      FFTCT<T> fft = FFTCT<T>(&gf, n);

      Vec<T> v(&gf, fft.get_n()), _v(&gf, fft.get_n()), v2(&gf, fft.get_n());
      for (int j = 0; j < 10000; j++) {
        v.zero_fill();
        for (int i = 0; i < n_data; i++)
          v.set(i, gf.weak_rand());
          // v.dump();
        fft.fft(&_v, &v);
          // _v.dump();
        fft.ifft(&v2, &_v);
          // v2.dump();
        assert(v.eq(&v2));
      }
    }
  }

  void run_taylor_expand(GF<T> *gf, FFTADD<T> *fft) {
    int t = 2 + gf->weak_rand() % (fft->get_n() - 2);
    int n = t + 1 + gf->weak_rand() % (fft->get_n() - t);
    Vec<T> v1(gf, n);
    int m = n/t;
    while (m*t < n) m++;
    Vec<T> v2(gf, t*m);
    for (int i = 0; i < n; i++)
      v1.set(i, gf->weak_rand());
    // v1.dump();
    fft->taylor_expand(&v2, &v1, n, t);
    Vec<T> _v1(gf, n);
    fft->inv_taylor_expand(&_v1, &v2, t);
    // _v1.dump();
    assert(_v1.eq(&v1));
  }

  // taylor expansion on (x^t - x)
  void test_taylor_expand(GF<T> *gf, FFTADD<T> *fft)
  {
    std::cout << "test_taylor_expand\n";
    for (int i = 0; i < 1000; i++)
      run_taylor_expand(gf, fft);
  }

  void run_taylor_expand_t2(GF<T> *gf, FFTADD<T> *fft) {
    int n = fft->get_n();
    Vec<T> v1(gf, n);
    for (int i = 0; i < n; i++)
      v1.set(i, gf->weak_rand());
    // v1.dump();
    fft->taylor_expand_t2(&v1, n, true);
    // v1.dump();
    Vec<T> _v1(gf, n);
    fft->inv_taylor_expand_t2(&_v1);
    // _v1.dump();
    assert(_v1.eq(&v1));
  }

  // taylor expansion on (x^2 - x)
  void test_taylor_expand_t2(GF<T> *gf, FFTADD<T> *fft)
  {
    std::cout << "test_taylor_expand_t2\n";
    for (int i = 0; i < 1000; i++)
      run_taylor_expand_t2(gf, fft);
  }

  void test_fftadd_codec(GF<T> *gf, FFTADD<T> *fft, int n_data)
  {
    std::cout << "test_fftadd_codec\n";
    Vec<T> v(gf, fft->get_n()), _v(gf, fft->get_n()), v2(gf, fft->get_n());
    for (int j = 0; j < 10000; j++) {
      v.zero_fill();
      for (int i = 0; i < n_data; i++)
        v.set(i, gf->weak_rand());
        // v.dump();
      fft->fft(&_v, &v);
        // _v.dump();
      fft->ifft(&v2, &_v);
        // v2.dump();
      assert(v.eq(&v2));
    }
  }

  void test_fftadd()
  {
    int n, m;
    int n_data = 3;
    int n_parities = 3;
    for (int gf_n = 4; gf_n <= 128 && gf_n <= 8 * sizeof(T); gf_n *= 2) {
      GF2N<T> gf = GF2N<T>(gf_n);
      std::cout << "test_fftadd_with_n=" << gf_n << "\n";
      // n is power of 2 and at least n_data + n_parities
      n = _get_smallest_power_of_2<T>(n_data + n_parities);
      m = _log2<T>(n);

      // std::cerr << "n=" << n << "\n";
      FFTADD<T> fft = FFTADD<T>(&gf, m);

      test_taylor_expand(&gf, &fft);
      test_taylor_expand_t2(&gf, &fft);
      test_fftadd_codec(&gf, &fft, n_data);
    }
  }

  void test_fft2_gfp()
  {
    GFP<T> gf = GFP<T>(3);
    T n_data = 1;

    std::cout << "test_fft2_gfp\n";

    // with this encoder we cannot exactly satisfy users request, we need to pad
    // n = minimal divisor of (q-1) that is at least (n_parities + n_data)
    // n = gf.get_code_len(n_parities + n_data);

    // std::cerr << "n=" << n << "\n";

    FFT2<T> fft = FFT2<T>(&gf);

    Vec<T> v(&gf, fft.get_n()), _v(&gf, fft.get_n()), v2(&gf, fft.get_n());
    for (int j = 0; j < 100000; j++) {
      v.zero_fill();
      for (int i = 0; i < n_data; i++)
        v.set(i, gf.weak_rand());
        // v.dump();
      fft.fft(&_v, &v);
        // _v.dump();
      fft.ifft(&v2, &_v);
        // v2.dump();
      assert(v.eq(&v2));
    }
  }

  void test_fft2()
  {
    u_int n;
    u_int r;
    u_int q = 65537;
    GFP<T> gf = GFP<T>(q);
    u_int R = gf.get_prime_root();  // primitive root
    u_int n_data = 3;
    u_int n_parities = 3;

    std::cout << "test_fft2\n";

    assert(_jacobi<T>(R, q) == -1);

    // with this encoder we cannot exactly satisfy users request, we need to pad
    // n = minimal divisor of (q-1) that is at least (n_parities + n_data)
    n = gf.get_code_len(n_parities + n_data);

    // compute root of order n-1 such as r^(n-1) mod q == 1
    r = gf.get_nth_root(n);

    // std::cerr << "r=" << r << "\n";

    DFTN<T> fft = DFTN<T>(&gf, n, r);

    Vec<T> v(&gf, fft.get_n()), _v(&gf, fft.get_n()), v2(&gf, fft.get_n());
    v.zero_fill();
    for (int i = 0; i < n_data; i++)
      v.set(i, gf.weak_rand());
    v.set(0, 27746);
    v.set(1, 871);
    v.set(2, 49520);
    // v.dump();
    fft.fft(&_v, &v);
    // _v.dump();
    assert(_v.get(0) == 12600);
    assert(_v.get(1) == 27885);
    assert(_v.get(2) == 17398);
    assert(_v.get(3) == 4624);
    assert(_v.get(4) == 10858);
    assert(_v.get(5) == 36186);
    assert(_v.get(6) == 4591);
    assert(_v.get(7) == 42289);
    fft.ifft(&v2, &_v);
    // v2.dump();
    assert(v.eq(&v2));
  }

  void test_fft_gf2n()
  {
    for (int n = 4; n <= 128 && n <= 8 * sizeof(T); n *= 2)
      test_fft_gf2n_with_n(n);
  }

  void test_fft_gf2n_with_n(int n)
  {
    T r;
    GF2N<T> gf = GF2N<T>(n);
    T R = gf.get_prime_root();
    T n_data = 3;
    T n_parities = 3;

    std::cout << "test_fft_gf2n_with_n=" << n << "\n";

    assert(gf.exp(R, gf.card_minus_one()) == 1);

    // with this encoder we cannot exactly satisfy users request, we need to pad
    n = gf.get_code_len(n_data + n_parities);

    r = gf.get_nth_root(n);
    assert(gf.exp(r, n) == 1);

    // std::cerr << "n=" << n << "\n";
    // std::cerr << "r=" << r << "\n";

    DFTN<T> fft = DFTN<T>(&gf, n, r);

    Vec<T> v(&gf, fft.get_n()), _v(&gf, fft.get_n()), v2(&gf, fft.get_n());
    for (T i = 0; i < 100000; i++) {
      v.zero_fill();
      for (int i = 0; i < n_data; i++)
        v.set(i, gf.weak_rand());
      // v.dump();
      fft.fft(&_v, &v);
      // _v.dump();
      fft.ifft(&v2, &_v);
      // v2.dump();
      assert(v.eq(&v2));
    }
  }

  void fft_utest_no_mul_bignum()
  {
    std::cout << "fft_utest_no_mul_bignum\n";

    test_gcd();
    test_quadratic_residues();
  }

  void fft_utest()
  {
    std::cout << "fft_utest\n";

    test_gcd();
    test_quadratic_residues();
    test_fftn();
    test_fft2k();
    test_fftpf();
    test_fftct_gfp();
    test_fftct_gf2n();
    test_fftadd();
    test_fft2();
    test_fft2_gfp();
    test_fft_gf2n();
    test_mul_bignum();
  }
};

template class Mat<uint32_t>;
template class Vec<uint32_t>;
template class DFT<uint32_t>;
template class FFTLN<uint32_t>;
template class FFTPF<uint32_t>;
template class FFTADD<uint32_t>;

template class Mat<uint64_t>;
template class Vec<uint64_t>;
template class DFT<uint64_t>;
template class FFTPF<uint64_t>;
template class FFTADD<uint64_t>;

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
  fftutest_uint32.fft_utest();
  FFTUtest<uint64_t> fftutest_uint64;
  fftutest_uint64.fft_utest();
  FFTUtest<mpz_class> fftutest_mpz;
  fftutest_mpz.fft_utest_no_mul_bignum();  // too slow
}
