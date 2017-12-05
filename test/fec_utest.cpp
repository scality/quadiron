
#include "ntl.h"

template<typename T>
class FECUtest
{
 public:
  void test_fecfntrs()
  {
    std::cout << "test_fecfntrs\n";
    u_int n_data = 3;
    u_int n_parities = 3;

    FECFNTRS<T> fec = FECFNTRS<T>(2, n_data, n_parities);
    run_test(&fec, fec.n, n_data, true);
  }

  void test_fecngff4rs()
  {
    std::cout << "test_fecngff4rs\n";
    u_int n_data = 3;
    u_int n_parities = 3;

    for (u_int word_size = 4; word_size <= sizeof(T)/2; word_size += 4) {
      std::cout << "test_fecngff4rs with word_size=" << word_size << "\n";
      FECNGFF4RS<T> fec = FECNGFF4RS<T>(word_size, n_data, n_parities);
      run_test(&fec, fec.n, n_data, true);
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
    run_test(&fec, fec.n, n_data);
  }

  void test_fecgf2nfftaddrs() {
    std::cout << "test_fecgf2nfftaddrs with sizeof(T)=" << sizeof(T) << "\n";
    for (int wordsize = 1; wordsize <= sizeof(T); wordsize *= 2)
      test_fecgf2nfftaddrs_with_wordsize(wordsize);
  }
  void test_fecgf2nfftaddrs_with_wordsize(int wordsize)
  {
    std::cout << "test_fecgf2nfftaddrs_with_wordsize=" << wordsize << "\n";

    u_int n_data = 3;
    u_int n_parities = 3;

    FECGF2NFFTADDRS<T> fec = FECGF2NFFTADDRS<T>(wordsize, n_data, n_parities);
    run_test(&fec, fec.n, n_data);
  }

  void test_fecgfpfftrs() {
    std::cout << "test_fecgfpfftrs with sizeof(T)=" << sizeof(T) << "\n";
    for (int word_size = 1; word_size <= 4 && word_size < sizeof(T);
      word_size *= 2)
      test_fecgfpfftrs_with_word_size(word_size);
  }
  void test_fecgfpfftrs_with_word_size(int word_size)
  {
    std::cout << "test_fecgfpfftrs_with_word_size=" << word_size << "\n";

    u_int n_data = 3;
    u_int n_parities = 3;

    FECGFPFFTRS<T> fec = FECGFPFFTRS<T>(word_size, n_data, n_parities);
    run_test(&fec, fec.n, n_data, true);
  }

  void run_test(FEC<T> *fec, int n, int n_data, bool propos_flag=false) {
    GF<T> *gf = fec->get_gf();

    Vec<T> v(gf, n_data), _v(gf, n), _v2(gf, n_data), f(gf, n_data),
      v2(gf, n_data);
    Vec<T> v_p(gf, n_data);
    std::vector<int> ids;
    for (int i = 0; i < n; i++)
      ids.push_back(i);
    std::vector<KeyValue*>props(n, nullptr);
    for (int j = 0; j < 10; j++) {
      if (propos_flag) {
        for (int i = 0; i < n; i++)
          props[i] = new KeyValue();
      }
      for (int i = 0; i < n_data; i++)
        v.set(i, gf->weak_rand());
      // FIXME: ngff4 will modify v after encode
      v_p.copy(&v);
      // std::cout << " v:"; v.dump();
      fec->encode(&_v, props, 0, &v);
      // std::cout << "_v:"; _v.dump();
      std::random_shuffle(ids.begin(), ids.end());
      for (int i = 0; i < n_data; i++) {
        f.set(i, ids.at(i));
        _v2.set(i, _v.get(ids.at(i)));
      }
      fec->decode(&v2, props, 0, &f, &_v2);
      // std::cout << "v2:"; v2.dump();
      assert(v_p.eq(&v2));
      if (propos_flag) {
        for (int i = 0; i < n; i++)
          delete props[i];
      }
    }
  }

  void fec_utest()
  {
    std::cout << "fec_utest\n";

    test_fecfntrs();
    test_fecngff4rs();
    test_fecgf2nfftrs();
    test_fecgf2nfftaddrs();
    test_fecgfpfftrs();
  }
};

template class Mat<uint32_t>;
template class Vec<uint32_t>;
template class FEC<uint32_t>;
template class FECGF2NRS<uint32_t>;
template class FECFNTRS<uint32_t>;
template class FECNGFF4RS<uint32_t>;
template class FECGF2NFFTRS<uint32_t>;

template class Mat<uint64_t>;
template class Vec<uint64_t>;
template class FEC<uint64_t>;
template class FECGF2NRS<uint64_t>;
template class FECFNTRS<uint64_t>;
template class FECNGFF4RS<uint64_t>;
template class FECGF2NFFTRS<uint64_t>;
template class FECGFPFFTRS<uint64_t>;

template class Mat<__uint128_t>;
template class Vec<__uint128_t>;
template class FEC<__uint128_t>;
template class FECGF2NRS<__uint128_t>;
// template class FECFNTRS<uint64_t>;
template class FECNGFF4RS<__uint128_t>;
template class FECGF2NFFTRS<__uint128_t>;
template class FECGFPFFTRS<__uint128_t>;

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
  fec_utest_uint128.test_fecngff4rs();
  // fecfntrs does not work for uint128_t
  fec_utest_uint128.test_fecgf2nfftrs();  // it runs slowly over gf2n(128)
}
