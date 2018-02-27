/*
 * Copyright 2017-2018 the NTTEC authors
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#include "nttec.h"

template <typename T>
class FFTUtest {
  public:
    void test_gcd1(nttec::gf::Field<T>* gf)
    {
        nttec::SignedDoubleSizeVal<T> bezout[2];

        assert(gf->inv(20) == 34);
        //
        int i;
        for (i = 0; i < 100; i++) {
            T x = gf->weak_rand();
            assert(1 == nttec::arith::extended_gcd<T>(97, x, bezout, nullptr));
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

        nttec::gf::Prime<T> gf(97);
        test_gcd1(&gf);
    }

    void test_quadratic_residues()
    {
        std::cout << "test_quadratic_residues\n";

        nttec::gf::Prime<T> gf32(32);
        int i;
        for (i = 0; i < 32; i++) {
            assert(gf32.is_quadratic_residue(gf32.exp(i, 2)));
        }

        nttec::gf::Prime<T> gf7(7);
        assert(gf7.is_quadratic_residue(2));
        assert(!gf7.is_quadratic_residue(5));

        nttec::gf::Prime<T> gf8(8);
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
    nttec::vec::Vector<T>*
    _convert_string2vec(nttec::gf::Field<T>* gf, int n, char num[])
    {
        int i;
        nttec::vec::Vector<T>* vec = new nttec::vec::Vector<T>(gf, n);
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

        int b = 10; // base
        int p = 14; // we could multiply integers of 2^p digits
        // int max_digits = nttec::arith::exp<T>(2, p);
        // std::cerr << "p=" << p << " max_digits=" << max_digits << "\n";

        uint64_t l = p + 1;
        // std::cerr << "l=" << l << "\n";

        // choose 2 prime numbers of the form p=a.2^n+1
        // because if x is not a quadratic residue then w=x^a is
        // a 2^n-th principal root of unity in GF_p
        uint64_t a1 = 2;
        uint64_t a2 = 5;
        uint64_t p1 = a1 * nttec::arith::exp<T>(2, 15) + 1;
        uint64_t p2 = a2 * nttec::arith::exp<T>(2, 15) + 1;
        // std::cerr << "p1=" << p1 << " p2=" << p2 << "\n";
        assert(nttec::arith::is_prime<T>(p1));
        assert(nttec::arith::is_prime<T>(p2));

        // ensure their product is bounded (b-1)^2*2^(n-1) < m
        uint64_t m = p1 * p2;
        // check overflow
        assert(m / p1 == p2);
        // std::cerr << " m=" << m << "\n";
        assert(
            nttec::arith::exp<T>((b - 1), 2) * nttec::arith::exp<T>(p, 2) < m);

        // find x so it is not a quadratic residue in GF_p1 and GF_p2
        assert(
            nttec::arith::jacobi<T>(3, p1) == nttec::arith::jacobi<T>(p1, 3));
        assert(nttec::arith::jacobi<T>(p1, 3) == nttec::arith::jacobi<T>(2, 3));
        assert(
            nttec::arith::jacobi<T>(3, p2) == nttec::arith::jacobi<T>(p2, 3));
        assert(nttec::arith::jacobi<T>(p2, 3) == nttec::arith::jacobi<T>(2, 3));
        assert(nttec::arith::jacobi<T>(2, 3) == -1);
        // which means x=3 is not a quadratic residue in GF_p1 and GF_p2

        // therefore we can compute 2^n-th roots of unity in GF_p1 and GF_p2
        uint64_t w1 = nttec::arith::exp<T>(3, a1);
        uint64_t w2 = nttec::arith::exp<T>(3, a2);
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
        uint64_t w = nttec::arith::chinese_remainder<uint64_t>(2, _a, _n);
        // std::cerr << " w=" << w << "\n";
        assert(w == 25559439);

        nttec::gf::Prime<T> gf_m(m);
        nttec::fft::Large<T> fft(&gf_m, l, 25559439);

        // parse the big numbers
        char X[] = "1236548787985654354598651354984132468";
        char Y[] = "745211515185321545554545854598651354984132468";

        nttec::vec::Vector<T>* _X = _convert_string2vec(&gf_m, fft.get_n(), X);
        // _X->dump();
        nttec::vec::Vector<T>* _Y = _convert_string2vec(&gf_m, fft.get_n(), Y);
        // _Y->dump();

        nttec::vec::Vector<T>* sfX =
            new nttec::vec::Vector<T>(&gf_m, fft.get_n());
        nttec::vec::Vector<T>* sfY =
            new nttec::vec::Vector<T>(&gf_m, fft.get_n());
        nttec::vec::Vector<T>* _XY =
            new nttec::vec::Vector<T>(&gf_m, fft.get_n());
        nttec::vec::Vector<T>* sfXY =
            new nttec::vec::Vector<T>(&gf_m, fft.get_n());

        fft.fft(sfX, _X);
        fft.fft(sfY, _Y);

        for (int i = 0; i <= fft.get_n() - 1; i++) {
            nttec::DoubleSizeVal<T> val =
                nttec::DoubleSizeVal<T>(sfX->get(i)) * sfY->get(i);
            _XY->set(i, val % m);
        }

        fft.ifft(sfXY, _XY);

        // carry propagation
        mpz_class z = 0;
        mpz_class u;
        for (int i = 0; i <= fft.get_n() - 1; i++) {
            mpz_class t, b;
            b = 10;
            mpz_pow_ui(t.get_mpz_t(), b.get_mpz_t(), i);
            u = mpz_class(std::to_string(sfXY->get(i)));
            z += u * t;
        }

        // std::cout << z << "\n";
        assert(
            z.get_str()
            == std::string("921490395895362412399910100421159322")
                   + "712298564831565484737491129935640058571771024");

        delete sfXY;
        delete _XY;
        delete sfX;
        delete sfY;
        delete _X;
        delete _Y;
    }

    void test_fft_naive()
    {
        unsigned n;
        unsigned r;
        T q = 65537;
        nttec::gf::Prime<T> gf(q);
        unsigned R = gf.get_primitive_root();
        unsigned n_data = 3;
        unsigned n_parities = 3;

        std::cout << "Test fft::Naive\n";

        assert(nttec::arith::jacobi<T>(R, q) == -1);

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        n = gf.get_code_len(n_parities + n_data);

        // compute root of order n-1 such as r^(n-1) mod q == 1
        r = gf.get_nth_root(n);

        // std::cerr << "n=" << n << "\n";
        // std::cerr << "r=" << r << "\n";

        nttec::fft::Naive<T> fft(&gf, n, r);

        nttec::vec::Vector<T> v(&gf, fft.get_n());
        nttec::vec::Vector<T> _v(&gf, fft.get_n());
        nttec::vec::Vector<T> v2(&gf, fft.get_n());
        for (int j = 0; j < 100000; j++) {
            v.zero_fill();
            for (unsigned i = 0; i < n_data; i++)
                v.set(i, gf.weak_rand());
            // v.dump();
            fft.fft(&_v, &v);
            // _v.dump();
            fft.ifft(&v2, &_v);
            // v2.dump();
            assert(v.eq(&v2));
        }
    }

    void test_fft_2k()
    {
        std::cout << "Test fft::Radix2\n";
        test_fft_2k_vec();
        test_fft_2k_vecp();
    }

    void test_fft_2k_vec()
    {
        unsigned n;
        unsigned q = 65537;
        nttec::gf::Prime<T> gf(q);
        unsigned R = gf.get_primitive_root();
        unsigned n_data = 3;
        unsigned n_parities = 3;

        std::cout << "Test fft::Radix2 vec\n";

        assert(nttec::arith::jacobi<T>(R, q) == -1);

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        n = gf.get_code_len(n_parities + n_data);

        // std::cerr << "n=" << n << "\n";

        nttec::fft::Radix2<T> fft(&gf, n);

        nttec::vec::Vector<T> v(&gf, fft.get_n());
        nttec::vec::Vector<T> _v(&gf, fft.get_n());
        nttec::vec::Vector<T> v2(&gf, fft.get_n());
        for (int j = 0; j < 100000; j++) {
            v.zero_fill();
            for (unsigned i = 0; i < n_data; i++)
                v.set(i, gf.weak_rand());
            // v.dump();
            fft.fft(&_v, &v);
            // _v.dump();
            fft.ifft(&v2, &_v);
            // v2.dump();
            assert(v.eq(&v2));
        }
    }

    void test_fft_2k_vecp()
    {
        unsigned n;
        unsigned q = 65537;
        nttec::gf::Prime<T> gf(q);
        unsigned R = gf.get_primitive_root();
        unsigned n_data = 3;
        unsigned n_parities = 3;
        size_t size = 4;

        std::cout << "Test fft::Radix2 vecp\n";

        assert(nttec::arith::jacobi<T>(R, q) == -1);

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        n = gf.get_code_len(n_parities + n_data);

        // std::cerr << "n=" << n << "\n";

        nttec::fft::Radix2<T> fft(&gf, n);

        int vec_n = fft.get_n();

        nttec::vec::Buffers<T> v(n_data, size);
        nttec::vec::Buffers<T> v2(vec_n, size);
        nttec::vec::Buffers<T> _v2(vec_n, size);
        nttec::vec::BuffersZeroExtended<T> _v(&v, vec_n);
        for (int j = 0; j < 100000; j++) {
            for (unsigned i = 0; i < n_data; i++) {
                T* mem = v.get(i);
                for (size_t u = 0; u < size; u++) {
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

    void test_fft_gt()
    {
        T n;
        nttec::gf::BinExtension<T> gf(4);
        T n_data = 3;
        T n_parities = 3;

        std::cout << "Test fft::GoodThomas\n";

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        n = gf.get_code_len(n_parities + n_data);

        // std::cerr << "n=" << n << "\n";

        nttec::fft::GoodThomas<T> fft(&gf, n);

        nttec::vec::Vector<T> v(&gf, fft.get_n());
        nttec::vec::Vector<T> _v(&gf, fft.get_n());
        nttec::vec::Vector<T> v2(&gf, fft.get_n());
        for (int j = 0; j < 10000; j++) {
            v.zero_fill();
            for (T i = 0; i < n_data; i++)
                v.set(i, gf.weak_rand());
            // v.dump();
            fft.fft(&_v, &v);
            // _v.dump();
            fft.ifft(&v2, &_v);
            // v2.dump();
            assert(v.eq(&v2));
        }
    }

    void test_fft_ct_gfp()
    {
        T n;
        T q = 65537;
        nttec::gf::Prime<T> gf(q);
        T n_data = 3;
        T n_parities = 3;

        std::cout << "Test fft::CooleyTukey on GF(p)\n";

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        n = gf.get_code_len(n_parities + n_data);

        // std::cerr << "n=" << n << "\n";

        nttec::fft::CooleyTukey<T> fft(&gf, n);

        nttec::vec::Vector<T> v(&gf, fft.get_n());
        nttec::vec::Vector<T> _v(&gf, fft.get_n());
        nttec::vec::Vector<T> v2(&gf, fft.get_n());
        for (int j = 0; j < 10000; j++) {
            v.zero_fill();
            for (T i = 0; i < n_data; i++)
                v.set(i, gf.weak_rand());
            // v.dump();
            fft.fft(&_v, &v);
            // _v.dump();
            fft.ifft(&v2, &_v);
            // v2.dump();
            assert(v.eq(&v2));
        }
    }

    void test_fft_ct_gf2n()
    {
        T n;
        T n_data = 3;
        T n_parities = 3;
        for (size_t gf_n = 4; gf_n <= 128 && gf_n <= 8 * sizeof(T); gf_n *= 2) {
            nttec::gf::BinExtension<T> gf(gf_n);

            std::cout << "Test fft::CooleyTukey on GF(2n)=" << gf_n << "\n";

            // with this encoder we cannot exactly satisfy users request, we
            // need to pad n = minimal divisor of (q-1) that is at least
            // (n_parities + n_data)
            n = gf.get_code_len(n_parities + n_data);

            // std::cerr << "n=" << n << "\n";

            nttec::fft::CooleyTukey<T> fft(&gf, n);

            nttec::vec::Vector<T> v(&gf, fft.get_n()), _v(&gf, fft.get_n()),
                v2(&gf, fft.get_n());
            for (int j = 0; j < 10000; j++) {
                v.zero_fill();
                for (T i = 0; i < n_data; i++)
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

    void
    run_taylor_expand(nttec::gf::Field<T>* gf, nttec::fft::Additive<T>* fft)
    {
        int t = 2 + gf->weak_rand() % (fft->get_n() - 2);
        int n = t + 1 + gf->weak_rand() % (fft->get_n() - t);
        nttec::vec::Vector<T> v1(gf, n);
        int m = n / t;
        while (m * t < n)
            m++;
        nttec::vec::Vector<T> v2(gf, t * m);
        for (int i = 0; i < n; i++)
            v1.set(i, gf->weak_rand());
        // v1.dump();
        fft->taylor_expand(&v2, &v1, n, t);
        nttec::vec::Vector<T> _v1(gf, n);
        fft->inv_taylor_expand(&_v1, &v2, t);
        // _v1.dump();
        assert(_v1.eq(&v1));
    }

    // taylor expansion on (x^t - x)
    void
    test_taylor_expand(nttec::gf::Field<T>* gf, nttec::fft::Additive<T>* fft)
    {
        std::cout << "test_taylor_expand\n";
        for (int i = 0; i < 1000; i++)
            run_taylor_expand(gf, fft);
    }

    void
    run_taylor_expand_t2(nttec::gf::Field<T>* gf, nttec::fft::Additive<T>* fft)
    {
        int n = fft->get_n();
        nttec::vec::Vector<T> v1(gf, n);
        for (int i = 0; i < n; i++)
            v1.set(i, gf->weak_rand());
        // v1.dump();
        fft->taylor_expand_t2(&v1, n, true);
        // v1.dump();
        nttec::vec::Vector<T> _v1(gf, n);
        fft->inv_taylor_expand_t2(&_v1);
        // _v1.dump();
        assert(_v1.eq(&v1));
    }

    // taylor expansion on (x^2 - x)
    void
    test_taylor_expand_t2(nttec::gf::Field<T>* gf, nttec::fft::Additive<T>* fft)
    {
        std::cout << "test_taylor_expand_t2\n";
        for (int i = 0; i < 1000; i++)
            run_taylor_expand_t2(gf, fft);
    }

    void test_fftadd_codec(
        nttec::gf::Field<T>* gf,
        nttec::fft::Additive<T>* fft,
        int n_data)
    {
        std::cout << "test_fftadd_codec\n";
        nttec::vec::Vector<T> v(gf, fft->get_n());
        nttec::vec::Vector<T> _v(gf, fft->get_n());
        nttec::vec::Vector<T> v2(gf, fft->get_n());
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
        for (size_t gf_n = 4; gf_n <= 128 && gf_n <= 8 * sizeof(T); gf_n *= 2) {
            nttec::gf::BinExtension<T> gf(gf_n);
            std::cout << "test_fftadd_with_n=" << gf_n << "\n";
            // n is power of 2 and at least n_data + n_parities
            n = nttec::arith::get_smallest_power_of_2<T>(n_data + n_parities);
            m = nttec::arith::log2<T>(n);

            // std::cerr << "n=" << n << "\n";
            nttec::fft::Additive<T> fft(&gf, m);

            test_taylor_expand(&gf, &fft);
            test_taylor_expand_t2(&gf, &fft);
            test_fftadd_codec(&gf, &fft, n_data);
        }
    }

    void test_fft_2_gfp()
    {
        nttec::gf::Prime<T> gf(3);
        T n_data = 1;

        std::cout << "Test fft::Size2 GF(p)\n";

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data) n = gf.get_code_len(n_parities + n_data);

        // std::cerr << "n=" << n << "\n";

        nttec::fft::Size2<T> fft(&gf);

        nttec::vec::Vector<T> v(&gf, fft.get_n());
        nttec::vec::Vector<T> _v(&gf, fft.get_n());
        nttec::vec::Vector<T> v2(&gf, fft.get_n());
        for (int j = 0; j < 100000; j++) {
            v.zero_fill();
            for (T i = 0; i < n_data; i++)
                v.set(i, gf.weak_rand());
            // v.dump();
            fft.fft(&_v, &v);
            // _v.dump();
            fft.ifft(&v2, &_v);
            // v2.dump();
            assert(v.eq(&v2));
        }
    }

    void test_fft_single_gfp()
    {
        nttec::gf::Prime<T> gf(39);

        std::cout << "Test fft::Single on GF(p)\n";

        int n = 16;
        nttec::fft::Single<T> fft(&gf, n);

        nttec::vec::Vector<T> v(&gf, fft.get_n());
        nttec::vec::Vector<T> _v(&gf, fft.get_n());
        nttec::vec::Vector<T> v2(&gf, fft.get_n());
        for (int j = 0; j < 100000; j++) {
            v.zero_fill();
            v.set(0, gf.weak_rand());
            // v.dump();
            fft.fft(&_v, &v);
            // _v.dump();
            fft.ifft(&v2, &_v);
            // v2.dump();
            assert(v.eq(&v2));
        }
    }

    void test_fft_2()
    {
        unsigned n;
        unsigned r;
        unsigned q = 65537;
        nttec::gf::Prime<T> gf(q);
        unsigned R = gf.get_primitive_root();
        unsigned n_data = 3;
        unsigned n_parities = 3;

        std::cout << "Test fft::Size2\n";

        assert(nttec::arith::jacobi<T>(R, q) == -1);

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad n = minimal divisor of (q-1) that is at least (n_parities +
        // n_data)
        n = gf.get_code_len(n_parities + n_data);

        // compute root of order n-1 such as r^(n-1) mod q == 1
        r = gf.get_nth_root(n);

        // std::cerr << "r=" << r << "\n";

        nttec::fft::Naive<T> fft(&gf, n, r);

        nttec::vec::Vector<T> v(&gf, fft.get_n());
        nttec::vec::Vector<T> _v(&gf, fft.get_n());
        nttec::vec::Vector<T> v2(&gf, fft.get_n());
        v.zero_fill();
        for (T i = 0; i < n_data; i++)
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
        for (size_t n = 4; n <= 128 && n <= 8 * sizeof(T); n *= 2)
            test_fft_gf2n_with_n(n);
    }

    void test_fft_gf2n_with_n(int n)
    {
        T r;
        nttec::gf::BinExtension<T> gf(n);
        T R = gf.get_primitive_root();
        T n_data = 3;
        T n_parities = 3;

        std::cout << "Test fft::Naive on gf2n with n=" << n << "\n";

        assert(gf.exp(R, gf.card_minus_one()) == 1);

        // with this encoder we cannot exactly satisfy users request, we need to
        // pad
        n = gf.get_code_len(n_data + n_parities);

        r = gf.get_nth_root(n);
        assert(gf.exp(r, n) == 1);

        // std::cerr << "n=" << n << "\n";
        // std::cerr << "r=" << r << "\n";

        nttec::fft::Naive<T> fft(&gf, n, r);

        nttec::vec::Vector<T> v(&gf, fft.get_n());
        nttec::vec::Vector<T> _v(&gf, fft.get_n());
        nttec::vec::Vector<T> v2(&gf, fft.get_n());
        for (T i = 0; i < 100000; i++) {
            v.zero_fill();
            for (T i = 0; i < n_data; i++)
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
        test_fft_naive();
        test_fft_2k();
        test_fft_gt();
        test_fft_ct_gfp();
        test_fft_ct_gf2n();
        test_fftadd();
        test_fft_2();
        test_fft_2_gfp();
        test_fft_single_gfp();
        test_fft_gf2n();
        test_mul_bignum();
    }
};

void fft_utest()
{
    FFTUtest<uint32_t> fftutest_uint32;
    fftutest_uint32.fft_utest();
    FFTUtest<uint64_t> fftutest_uint64;
    fftutest_uint64.fft_utest();
    FFTUtest<mpz_class> fftutest_mpz;
    fftutest_mpz.fft_utest_no_mul_bignum(); // too slow
}
