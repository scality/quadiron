#include "nttec.h"

template <typename T>
class FECUtest {
  public:
    void test_fecfntrs()
    {
        std::cout << "test_fecfntrs\n";
        unsigned n_data = 3;
        unsigned n_parities = 3;

        nttec::fec::FECFNTRS<T> fec(2, n_data, n_parities);
        run_test(&fec, fec.n, n_data, n_data + n_parities, true);
    }

    void test_fecngff4rs()
    {
        std::cout << "test_fecngff4rs\n";
        unsigned n_data = 3;
        unsigned n_parities = 3;

        for (int i = 1; i < nttec::arith::log2<T>(sizeof(T)); i++) {
            unsigned word_size = 1 << i;
            std::cout << "test_fecngff4rs with word_size=" << word_size << "\n";
            nttec::fec::FECNGFF4RS<T> fec(word_size, n_data, n_parities);
            run_test(&fec, fec.n, n_data, n_data + n_parities, true);
        }
    }

    void test_fecgf2nfftrs()
    {
        std::cout << "test_fecgf2nfftrs with sizeof(T)=" << sizeof(T) << "\n";
        for (size_t wordsize = 1; wordsize <= sizeof(T); wordsize *= 2)
            test_fecgf2nfftrs_with_wordsize(wordsize);
    }
    void test_fecgf2nfftrs_with_wordsize(int wordsize)
    {
        std::cout << "test_fecgf2nfftrs_with_wordsize=" << wordsize << "\n";

        unsigned n_data = 3;
        unsigned n_parities = 3;

        nttec::fec::FECGF2NFFTRS<T> fec(wordsize, n_data, n_parities);
        run_test(&fec, fec.n, n_data, n_data + n_parities);
    }

    void test_fecgf2nfftaddrs()
    {
        std::cout << "test_fecgf2nfftaddrs with sizeof(T)=" << sizeof(T)
                  << "\n";
        for (size_t wordsize = 1; wordsize <= sizeof(T); wordsize *= 2)
            test_fecgf2nfftaddrs_with_wordsize(wordsize);
    }
    void test_fecgf2nfftaddrs_with_wordsize(int wordsize)
    {
        std::cout << "test_fecgf2nfftaddrs_with_wordsize=" << wordsize << "\n";

        unsigned n_data = 3;
        unsigned n_parities = 3;

        nttec::fec::FECGF2NFFTADDRS<T> fec(wordsize, n_data, n_parities);
        run_test(&fec, fec.n, n_data, n_data + n_parities);
    }

    void test_fecgfpfftrs()
    {
        std::cout << "test_fecgfpfftrs with sizeof(T)=" << sizeof(T) << "\n";
        for (size_t word_size = 1; word_size <= 4 && word_size < sizeof(T);
             word_size *= 2)
            test_fecgfpfftrs_with_word_size(word_size);
    }
    void test_fecgfpfftrs_with_word_size(int word_size)
    {
        std::cout << "test_fecgfpfftrs_with_word_size=" << word_size << "\n";

        unsigned n_data = 3;
        unsigned n_parities = 3;

        nttec::fec::FECGFPFFTRS<T> fec(word_size, n_data, n_parities);
        run_test(&fec, fec.n, n_data, n_data + n_parities, true);
    }

    void run_test(
        nttec::fec::FEC<T>* fec,
        int n,
        int n_data,
        int code_len,
        bool propos_flag = false)
    {
        nttec::gf::GF<T>* gf = fec->get_gf();

        nttec::vec::Vec<T> v(gf, n_data);
        nttec::vec::Vec<T> _v(gf, n);
        nttec::vec::Vec<T> _v2(gf, n_data);
        nttec::vec::Vec<T> f(gf, n_data);
        nttec::vec::Vec<T> v2(gf, n_data);
        nttec::vec::Vec<T> v_p(gf, n_data);
        std::vector<int> ids;
        for (int i = 0; i < code_len; i++)
            ids.push_back(i);
        std::vector<nttec::Properties> props(code_len);
        for (int j = 0; j < 1000; j++) {
            if (propos_flag) {
                for (int i = 0; i < code_len; i++)
                    props[i] = nttec::Properties();
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

void fec_utest()
{
    FECUtest<uint32_t> fec_utest_uint32;
    fec_utest_uint32.fec_utest();
    FECUtest<uint64_t> fec_utest_uint64;
    fec_utest_uint64.fec_utest();
    FECUtest<__uint128_t> fec_utest_uint128;
    fec_utest_uint128.test_fecngff4rs();
    // fecfntrs does not work for uint128_t
    fec_utest_uint128.test_fecgf2nfftrs(); // it runs slowly over gf2n(128)
}
