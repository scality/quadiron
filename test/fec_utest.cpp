/*
 * Copyright 2017-2018 Scality
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
#include <gtest/gtest.h>

#include "quadiron.h"

template <typename T>
class FecTestCommon : public ::testing::Test {
  public:
    const unsigned n_data = 3;
    const unsigned n_parities = 3;

    void run_test(
        quadiron::fec::FecCode<T>& fec,
        int n_data,
        int code_len,
        bool props_flag = false)
    {
        const quadiron::gf::Field<T>* gf = &(fec.get_gf());

        quadiron::vec::Vector<T> v(*gf, n_data);
        quadiron::vec::Vector<T> _v(*gf, fec.n);
        quadiron::vec::Vector<T> _v2(*gf, n_data);
        quadiron::vec::Vector<T> f(*gf, n_data);
        quadiron::vec::Vector<T> v2(*gf, n_data);
        quadiron::vec::Vector<T> v_p(*gf, n_data);
        std::vector<int> ids;

        for (int i = 0; i < code_len; i++) {
            ids.push_back(i);
        }

        std::vector<quadiron::Properties> props(code_len);
        for (int j = 0; j < 1000; j++) {
            if (props_flag) {
                for (int i = 0; i < code_len; i++) {
                    props[i] = quadiron::Properties();
                }
            }
            for (int i = 0; i < n_data; i++) {
                v.set(i, gf->weak_rand());
            }
            // FIXME: ngff4 will modify v after encode
            v_p.copy(&v);
            fec.encode(&_v, props, 0, &v);
            std::random_shuffle(ids.begin(), ids.end());
            for (int i = 0; i < n_data; i++) {
                f.set(i, ids.at(i));
                _v2.set(i, _v.get(ids.at(i)));
            }
            std::unique_ptr<quadiron::fec::DecodeContext<T>> context =
                fec.init_context_dec(f);
            fec.decode(*context, &v2, props, 0, &_v2);
            ASSERT_TRUE(v_p.eq(&v2));
        }
    }
};

using AllTypes = ::testing::Types<uint32_t, uint64_t, __uint128_t>;
TYPED_TEST_CASE(FecTestCommon, AllTypes);

TYPED_TEST(FecTestCommon, TestNf4) // NOLINT
{
    const int iter_count = quadiron::arith::log2<TypeParam>(sizeof(TypeParam));

    for (int i = 1; i < iter_count; i++) {
        const unsigned word_size = 1 << i;
        quadiron::fec::RsNf4<TypeParam> fec(
            word_size, this->n_data, this->n_parities);

        this->run_test(
            fec, this->n_data, this->n_data + this->n_parities, true);
    }
}

TYPED_TEST(FecTestCommon, TestGf2nFft) // NOLINT
{
    for (size_t wordsize = 1; wordsize <= sizeof(TypeParam); wordsize *= 2) {
        quadiron::fec::RsGf2nFft<TypeParam> fec(
            wordsize, this->n_data, this->n_parities);

        this->run_test(fec, this->n_data, this->n_data + this->n_parities);
    }
}

TYPED_TEST(FecTestCommon, TestGf2nFftAdd) // NOLINT
{
    for (size_t wordsize = 1; wordsize <= sizeof(TypeParam); wordsize *= 2) {
        quadiron::fec::RsGf2nFftAdd<TypeParam> fec(
            wordsize, this->n_data, this->n_parities);

        this->run_test(fec, this->n_data, this->n_data + this->n_parities);
    }
}

template <typename T>
class FecTestNo128 : public FecTestCommon<T> {
};

using No128 = ::testing::Types<uint32_t, uint64_t>;
TYPED_TEST_CASE(FecTestNo128, No128);

TYPED_TEST(FecTestNo128, TestFnt) // NOLINT
{
    quadiron::fec::RsFnt<TypeParam> fec(2, this->n_data, this->n_parities);
    this->run_test(fec, this->n_data, this->n_data + this->n_parities, true);
}

TYPED_TEST(FecTestNo128, TestGfpFft) // NOLINT
{
    for (size_t word_size = 1; word_size <= 4 && word_size < sizeof(TypeParam);
         word_size *= 2) {
        quadiron::fec::RsGfpFft<TypeParam> fec(
            word_size, this->n_data, this->n_parities);

        this->run_test(
            fec, this->n_data, this->n_data + this->n_parities, true);
    }
}
