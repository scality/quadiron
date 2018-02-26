#include "nttec.h"

template <typename T>
class FECUtest {
  public:
    void test_fec_rs_fnt()
    {
        std::cout << "test fec::RsFnt\n";
        unsigned n_data = 3;
        unsigned n_parities = 3;

        nttec::fec::RsFnt<T> fec(2, n_data, n_parities);
        run_test(&fec, fec.n, n_data, n_data + n_parities, true);
    }

    void test_fec_rs_nf4()
    {
        std::cout << "Test fec::RsNf4\n";
        unsigned n_data = 3;
        unsigned n_parities = 3;

        for (int i = 1; i < nttec::arith::log2<T>(sizeof(T)); i++) {
            unsigned word_size = 1 << i;
            std::cout << "Test fec::RsNgff4 with word_size=" << word_size
                      << '\n';
            nttec::fec::RsNf4<T> fec(word_size, n_data, n_parities);
            run_test(&fec, fec.n, n_data, n_data + n_parities, true);
        }
    }

    void test_fec_rs_gf2n_fft()
    {
        std::cout << "Test fec::RsGf2nFft with sizeof(T)=" << sizeof(T) << '\n';
        for (size_t wordsize = 1; wordsize <= sizeof(T); wordsize *= 2) {
            test_fec_rs_gf2n_fft_with_wordsize(wordsize);
        }
    }
    void test_fec_rs_gf2n_fft_with_wordsize(int wordsize)
    {
        std::cout << "Test fec::RsGf2nFft with wordsize=" << wordsize << '\n';

        unsigned n_data = 3;
        unsigned n_parities = 3;

        nttec::fec::RsGf2nFft<T> fec(wordsize, n_data, n_parities);
        run_test(&fec, fec.n, n_data, n_data + n_parities);
    }

    void test_fec_rs_gf2n_fft_add()
    {
        std::cout << "Test fec::RsGf2nFftAdd with sizeof(T)=" << sizeof(T)
                  << '\n';
        for (size_t wordsize = 1; wordsize <= sizeof(T); wordsize *= 2)
            test_fec_rs_gf2n_fft_add_with_wordsize(wordsize);
    }
    void test_fec_rs_gf2n_fft_add_with_wordsize(int wordsize)
    {
        std::cout << "Test fec::RsGf2nFftAdd with wordsize=" << wordsize
                  << '\n';

        unsigned n_data = 3;
        unsigned n_parities = 3;

        nttec::fec::RsGf2nFftAdd<T> fec(wordsize, n_data, n_parities);
        run_test(&fec, fec.n, n_data, n_data + n_parities);
    }

    void test_fec_rs_gfp_fft()
    {
        std::cout << "Test fec::RsGfpFft with sizeof(T)=" << sizeof(T) << '\n';
        for (size_t word_size = 1; word_size <= 4 && word_size < sizeof(T);
             word_size *= 2) {
            test_fec_rs_gfp_fft_with_word_size(word_size);
        }
    }
    void test_fec_rs_gfp_fft_with_word_size(int word_size)
    {
        std::cout << "Test fec::RsGfpFft with word size=" << word_size << '\n';

        unsigned n_data = 3;
        unsigned n_parities = 3;

        nttec::fec::RsGfpFft<T> fec(word_size, n_data, n_parities);
        run_test(&fec, fec.n, n_data, n_data + n_parities, true);
    }

    void run_test(
        nttec::fec::FecCode<T>* fec,
        int n,
        int n_data,
        int code_len,
        bool props_flag = false)
    {
        nttec::gf::Field<T>* gf = fec->get_gf();

        nttec::vec::Vector<T> v(gf, n_data);
        nttec::vec::Vector<T> _v(gf, n);
        nttec::vec::Vector<T> _v2(gf, n_data);
        nttec::vec::Vector<T> f(gf, n_data);
        nttec::vec::Vector<T> v2(gf, n_data);
        nttec::vec::Vector<T> v_p(gf, n_data);
        std::vector<int> ids;
        for (int i = 0; i < code_len; i++)
            ids.push_back(i);
        std::vector<nttec::Properties> props(code_len);
        for (int j = 0; j < 1000; j++) {
            if (props_flag) {
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
        std::cout << "Unit tests for Forward Error Correction codes\n";

        test_fec_rs_fnt();
        test_fec_rs_nf4();
        test_fec_rs_gf2n_fft();
        test_fec_rs_gf2n_fft_add();
        test_fec_rs_gfp_fft();
    }
};

void fec_utest()
{
    FECUtest<uint32_t> fec_utest_uint32;
    fec_utest_uint32.fec_utest();
    FECUtest<uint64_t> fec_utest_uint64;
    fec_utest_uint64.fec_utest();
    FECUtest<__uint128_t> fec_utest_uint128;
    fec_utest_uint128.test_fec_rs_nf4();
    // fec::RsFnt does not work for uint128_t
    fec_utest_uint128.test_fec_rs_gf2n_fft(); // it runs slowly over gf2n(128)
}
