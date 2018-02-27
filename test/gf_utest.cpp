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
class GFUtest {
  public:
    void simple_tests()
    {
        std::cout << "simple tests\n";

        nttec::gf::Extension<T> gf256(
            2, 8); // dumb way to declare a binary field but for sake of testing
        nttec::Polynomial<T> _x(gf256.get_sub_field());
        _x.set(6, 1);
        _x.set(4, 1);
        _x.set(1, 1);
        _x.set(0, 1);
        T x = gf256.inv(_x.to_num());
        _x.from_num(x, 31);
        assert(_x.degree() == 7);
        assert(_x.get(7) == 1);
        assert(_x.get(6) == 1);
        assert(_x.get(5) == 0);
        assert(_x.get(4) == 0);
        assert(_x.get(3) == 1);
        assert(_x.get(2) == 0);
        assert(_x.get(1) == 1);
        assert(_x.get(0) == 0);

        nttec::gf::Extension<T> gf27(3, 3);
        nttec::Polynomial<T> _y(gf27.get_sub_field());
        _y.set(2, 1);
        _y.set(0, 1);
        T y = gf27.inv(_y.to_num());
        _y.from_num(y, 31);
        assert(_y.degree() == 2);
        assert(_y.get(2) == 1);
        assert(_y.get(1) == 1);
        assert(_y.get(0) == 0);
    }

    void test_negation(nttec::gf::Field<T>* gf)
    {
        int i;

        for (i = 0; i < 100; i++) {
            T x, y;

            // std::cout << "i=" << i << "\n";

            x = gf->weak_rand();
            // std::cout << "x=" << x << "\n";
            y = gf->neg(x);
            // std::cout << "inv(x)=" << y << "\n";
            assert(gf->add(x, y) == 0);
        }
    }

    void test_negation_gf_nf4(nttec::gf::NF4<T>* gf)
    {
        int i;

        for (i = 0; i < 100; i++) {
            T x, y;

            // std::cout << "i=" << i << "\n";

            x = gf->weak_rand_tuple();
            // std::cout << "x=" << x << "\n";
            y = gf->neg(x);
            // std::cout << "inv(x)=" << y << "\n";
            assert(gf->add(x, y) == 0);
        }
    }

    void test_reciprocal(nttec::gf::Field<T>* gf)
    {
        int i;
        int n_found = 0;

        for (i = 0; i < 100; i++) {
            T x, y;

            // std::cout << "i=" << i << "\n";

            x = gf->weak_rand();
            // std::cout << "x=" << x << "\n";
            try {
                y = gf->inv(x);
            } catch (nttec::Exception e) {
                continue;
            }
            // std::cout << "inv(x)=" << y << "\n";
            assert(gf->mul(x, y) == 1);
            n_found++;
        }
        assert(n_found > 0);
    }

    void test_reciprocal_gf_nf4(nttec::gf::NF4<T>* gf)
    {
        int i;
        int n_found = 0;

        for (i = 0; i < 100; i++) {
            T x, y;

            // std::cout << "i=" << i << "\n";

            x = gf->weak_rand_tuple();
            // std::cout << "x=" << x << "\n";
            try {
                y = gf->inv(x);
            } catch (nttec::Exception e) {
                continue;
            }
            // std::cout << "inv(x)=" << y << "\n";
            assert(gf->mul(x, y) == gf->get_unit());
            n_found++;
        }
        assert(n_found > 0);
    }

    void test_log(nttec::gf::Field<T>* gf)
    {
        int i;
        int n_found = 0;

        for (i = 0; i < 1000; i++) {
            T x, y, z, t;

            // std::cout << "i=" << i << "\n";
            // std::cout << gf->card() << "\n";
            x = gf->weak_rand();
            y = gf->weak_rand();
            // std::cout << "x=" << x << "\n";
            // std::cout << "y=" << y << "\n";
            try {
                z = gf->log(x, y);
                // std::cout << "z=" << z << "\n";
            } catch (...) {
                // std::cout << "not found\n";
                continue;
            }
            t = gf->exp(x, z);
            // std::cout << "t=" << t << "\n";
            assert(t == y);
            n_found++;
        }
        assert(n_found > 0);
    }

    void test_pack_unpack(nttec::gf::NF4<T>* gf)
    {
        int i;

        for (i = 0; i < 100; i++) {
            T x, y;
            nttec::GroupedValues<T> z;

            // std::cout << "i=" << i << "\n";

            x = gf->weak_rand_tuple();
            // std::cout << "x=" << x << "\n";
            z = gf->unpack(x);
            // std::cout << "unpack(x)=" << z.values << "\n";
            y = gf->pack(z.values, z.flag);
            // std::cout << "pack(z)=" << y << "\n";
            assert(x == y);
        }
    }

    void test_find_prime_root(nttec::gf::Field<T>* gf)
    {
        gf->find_prime_root();
        // std::cout << "root " << gf->root << std::endl;
        assert(gf->check_prime_root(gf->get_root()));
    }

    void test_get_order(nttec::gf::Field<T>* gf)
    {
        int i;
        T x;
        T order;
        T h = gf->card_minus_one();
        for (i = 0; i < 1000; i++) {
            // std::cout << "i=" << i << "\n";
            // std::cout << gf->card() << "\n";
            x = gf->weak_rand();
            order = gf->get_order(x);
            assert(gf->exp(x, order) == 1);
            assert(h % order == 0);
        }
    }

    void test_get_nth_root(nttec::gf::Field<T>* gf)
    {
        int i;
        T x;
        T nth_root;
        for (i = 0; i < 1000; i++) {
            // std::cout << "i=" << i << "\n";
            // std::cout << gf->card() << "\n";
            x = gf->weak_rand();
            nth_root = gf->get_nth_root(x);
            // std::cout << "x=" << x << " " << nth_root << std::endl;
            assert(gf->exp(nth_root, x) == 1);
        }
    }

    void test_negation_gf5()
    {
        std::cout << "test_negation_gf5\n";

        nttec::gf::Prime<T> gf5(5);
        test_negation(&gf5);
    }

    void test_reciprocal_gf5()
    {
        std::cout << "test_reciprocal_gf5\n";

        nttec::gf::Prime<T> gf5(5);
        test_reciprocal(&gf5);
    }

    void test_log_gf5()
    {
        std::cout << "test_log_gf5\n";

        nttec::gf::Prime<T> gf5(5);
        test_log(&gf5);
    }

    void test_prime_root_gf5()
    {
        std::cout << "test_prime_root_gf5\n";

        nttec::gf::Prime<T> gf5(5);
        test_find_prime_root(&gf5);
        test_get_order(&gf5);
        test_get_nth_root(&gf5);
    }

    void test_negation_gf9()
    {
        std::cout << "test_negation_gf9\n";

        nttec::gf::Extension<T> gf9(3, 2);
        test_negation(&gf9);
    }

    void test_reciprocal_gf9()
    {
        std::cout << "test_reciprocal_gf9\n";

        nttec::gf::Extension<T> gf9(3, 2);
        test_reciprocal(&gf9);
    }

    void test_log_gf9()
    {
        std::cout << "test_log_gf9\n";

        nttec::gf::Extension<T> gf9(3, 2);
        test_log(&gf9);
    }

    void test_prime_root_gf9()
    {
        std::cout << "test_prime_root_gf9\n";

        nttec::gf::Extension<T> gf9(3, 2);
        test_find_prime_root(&gf9);
        test_get_order(&gf9);
        test_get_nth_root(&gf9);
    }

    void test_negation_gf256()
    {
        std::cout << "test_negation_gf256\n";

        nttec::gf::BinExtension<T> gf256(8);
        test_negation(&gf256);
    }

    void test_reciprocal_gf256()
    {
        std::cout << "test_reciprocal_gf256\n";

        nttec::gf::BinExtension<T> gf256(8);
        test_reciprocal(&gf256);
    }

    void test_log_gf256()
    {
        std::cout << "test_log_gf256\n";

        nttec::gf::BinExtension<T> gf256(8);
        test_log(&gf256);
    }

    void test_prime_root_gf256()
    {
        std::cout << "test_prime_root_gf256\n";

        nttec::gf::BinExtension<T> gf256(8);
        test_find_prime_root(&gf256);
        test_get_order(&gf256);
        test_get_nth_root(&gf256);
    }

    void test_negation_gf2_bign(T n)
    {
        std::cout << "test_negation_gf(2^" << n << ")\n";

        nttec::gf::BinExtension<T> gf2n(n);
        test_negation(&gf2n);
    }

    void test_reciprocal_gf2_bign(T n)
    {
        std::cout << "test_reciprocal_gf(2^" << n << ")\n";

        nttec::gf::BinExtension<T> gf2n(n);
        test_reciprocal(&gf2n);
    }

    void test_log_gf2_bign(T n)
    {
        std::cout << "test_log_gf(2^" << n << ")\n";

        nttec::gf::BinExtension<T> gf2n(n);
        test_log(&gf2n);
    }

    void test_prime_root_gf2(T n)
    {
        std::cout << "test_prime_root_gf(2^" << n << ")\n";

        nttec::gf::BinExtension<T> gf2n(n);
        test_find_prime_root(&gf2n);
        test_get_order(&gf2n);
        test_get_nth_root(&gf2n);
    }

    void gf_utest_gf_nf4_with_n(T n)
    {
        std::cout << "gf_utest_gf_nf4_with_n=" << n << "\n";

        srand(time(0));

        nttec::gf::NF4<T> gf(n);

        test_negation_gf_nf4(&gf);
        test_reciprocal_gf_nf4(&gf);
        test_pack_unpack(&gf);
    }

    void gf_utest_gf_nf4()
    {
        int max_n = sizeof(T) / 4;
        std::cout << "sizeofT=" << sizeof(T)
                  << "gf_utest_gf_nf4 for max_n=" << max_n << "\n";
        for (int n = 1; n <= max_n; n++) {
            gf_utest_gf_nf4_with_n(n);
        }
    }

    void gf_utest()
    {
        std::cout << "gf_utest\n";

        srand(time(0));

        simple_tests();
        test_negation_gf5();
        test_reciprocal_gf5();
        test_log_gf5();
        test_prime_root_gf5();
        test_negation_gf9();
        test_reciprocal_gf9();
        test_log_gf9();
        test_prime_root_gf9();
        test_negation_gf256();
        test_reciprocal_gf256();
        test_log_gf256();
        test_negation_gf256();
        test_prime_root_gf256();
        gf_utest_gf_nf4();
    }

    void gf_utest_gf_2_bign(T n)
    {
        std::cout << "gf_utest 2^" << n << "\n";

        srand(time(0));

        test_negation_gf2_bign(n);
        test_reciprocal_gf2_bign(n);
        if (n < 128) // it works slowly in GF2N(128)
            test_prime_root_gf2(n);
        // test_log_gf2_bign(n);
    }

    void gf_utest_nogf2n()
    {
        std::cout << "gf_utest_nogf2n\n";

        srand(time(0));

        test_negation_gf5();
        test_reciprocal_gf5();
        test_log_gf5();
        test_prime_root_gf5();
    }

    void gf_utest_gf_2_n()
    {
        T max_n = 8 * sizeof(T);
        std::cout << "gf_utest_gf_2_n for max_n=" << max_n << "\n";
        for (T i = 8; i <= max_n; i *= 2) {
            gf_utest_gf_2_bign(i);
        }
    }
};

void gf_utest()
{
    GFUtest<uint32_t> gfutest_uint32;
    gfutest_uint32.gf_utest();
    gfutest_uint32.gf_utest_gf_2_n();
    GFUtest<uint64_t> gfutest_uint64;
    gfutest_uint64.gf_utest();
    gfutest_uint64.gf_utest_gf_2_n();
    GFUtest<__uint128_t> gfutest_uint128;
    // gfutest_uint128.gf_utest(); // gfp(n) does not work for uint128
    gfutest_uint128.gf_utest_gf_nf4();
    gfutest_uint128.gf_utest_gf_2_n();
    GFUtest<mpz_class> gfutest_mpz;
    gfutest_mpz.gf_utest_nogf2n(); // XXX gf2n broken for now
}
