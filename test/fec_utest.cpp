
#include "ntl.h"

template<typename T>
class FECUtest
{
 public:
  void test_fecfntrs()
  {
    u_int n_data = 3;
    u_int n_parities = 3;

    FECFNTRS<T> fec = FECFNTRS<T>(2, n_data, n_parities);
    GF<T> *gf = fec.get_gf();

    Vec<T> v(gf, n_data), _v(gf, fec.n), _v2(gf, n_data), f(gf, n_data),
      v2(gf, n_data);
    for (int j = 0; j < 10000; j++) {
      std::vector<KeyValue*>props(fec.n, nullptr);
      for (int i = 0; i < fec.n; i++)
        props[i] = new KeyValue();
      for (int i = 0; i < n_data; i++)
        v.set(i, gf->weak_rand());
      // v.dump();
      fec.encode(&_v, props, 0, &v);
      // _v.dump();
      _v2.copy(&_v, n_data);
      for (int i = 0; i < n_data; i++)
        f.set(i, i);
      fec.decode(&v2, props, 0, &f, &_v2);
      // v2.dump();
      assert(v.eq(&v2));
      for (int i = 0; i < fec.n; i++)
        delete props[i];
     }
  }

  void test_fecgf2nfftrs() {
    std::cout << "test_fecgf2nfftrs with sizeof(T)=" << sizeof(T) << "\n";
    for (int wordsize = 1; wordsize <= sizeof(T); wordsize *= 2)
      test_fecgf2nfftrs_with_wordsize(wordsize);
  }
  void test_fecgf2nfftrs_with_wordsize(int wordsize)
  {
    std::cout << "test_fecgf2nfftrs_with_wordsize=" << wordsize << "\n";

    u_int n_data = 3;
    u_int n_parities = 3;

    FECGF2NFFTRS<T> fec = FECGF2NFFTRS<T>(wordsize, n_data, n_parities);
    GF<T> *gf = fec.get_gf();

    std::vector<KeyValue*> props;
    Vec<T> v(gf, n_data), _v(gf, fec.n), _v2(gf, n_data), f(gf, n_data),
      v2(gf, n_data);
    for (int j = 0; j < 10000; j++) {
      for (int i = 0; i < n_data; i++)
        v.set(i, gf->weak_rand());
      // v.dump();
      fec.encode(&_v, props, 0, &v);
      // _v.dump();
      _v2.copy(&_v, n_data);
      for (int i = 0; i < n_data; i++)
        f.set(i, i);
      fec.decode(&v2, props, 0, &f, &_v2);
      // v2.dump();
      assert(v.eq(&v2));
     }
  }

  void fec_utest()
  {
    std::cout << "fec_utest\n";

    test_fecfntrs();
    test_fecgf2nfftrs();
  }
};

template class Mat<uint32_t>;
template class Vec<uint32_t>;
template class FEC<uint32_t>;
template class FECGF2NRS<uint32_t>;
template class FECFNTRS<uint32_t>;
template class FECGF2NFFTRS<uint32_t>;

template class Mat<uint64_t>;
template class Vec<uint64_t>;
template class FEC<uint64_t>;
template class FECGF2NRS<uint64_t>;
template class FECFNTRS<uint64_t>;
template class FECGF2NFFTRS<uint64_t>;

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

void fec_utest()
{
  FECUtest<uint32_t> fec_utest_uint32;
  fec_utest_uint32.fec_utest();
  FECUtest<uint64_t> fec_utest_uint64;
  fec_utest_uint64.fec_utest();
  FECUtest<__uint128_t> fec_utest_uint128;
  // fecfntrs does not work for uint128_t
  fec_utest_uint128.test_fecgf2nfftrs();  // it runs slowly over gf2n(128)
}
