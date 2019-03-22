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

namespace fec = quadiron::fec;
namespace arith = quadiron::arith;
namespace vec = quadiron::vec;

template <typename T>
class FecTestCommon : public ::testing::Test {
  public:
    const unsigned n_data = 3;
    const unsigned n_parities = 3;

    void
    run_test(fec::FecCode<T>& fec, bool props_flag = false, bool is_nf4 = false)
    {
        const int code_len = n_data + n_parities;

        const quadiron::gf::Field<T>& gf = fec.get_gf();
        const quadiron::gf::NF4<T>& nf4 =
            static_cast<const quadiron::gf::NF4<T>&>(gf);

        vec::Vector<T> data_frags(gf, n_data);
        vec::Vector<T> copied_data_frags(gf, n_data);
        vec::Vector<T> encoded_frags(gf, fec.n);
        vec::Vector<T> received_frags(gf, n_data);
        vec::Vector<T> decoded_frags(gf, n_data);
        std::vector<int> ids;
        vec::Vector<T> fragments_ids(gf, n_data);

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

            // gen_data(gf, data_frags, is_nf4);
            for (unsigned i = 0; i < n_data; i++) {
                const T v = is_nf4 ? nf4.unpacked_rand() : gf.rand();
                data_frags.set(i, v);
            }
            // FIXME: ngff4 will modify v after encode
            copied_data_frags.copy(&data_frags);

            fec.encode(encoded_frags, props, 0, data_frags);

            std::random_shuffle(ids.begin(), ids.end());
            for (unsigned i = 0; i < n_data; i++) {
                fragments_ids.set(i, ids.at(i));
                received_frags.set(i, encoded_frags.get(ids.at(i)));
            }
            std::unique_ptr<fec::DecodeContext<T>> context =
                fec.init_context_dec(fragments_ids, props);

            fec.decode(*context, decoded_frags, props, 0, received_frags);

            ASSERT_EQ(copied_data_frags, decoded_frags);
        }
    }
};

using AllTypes = ::testing::Types<uint32_t, uint64_t, __uint128_t>;
TYPED_TEST_CASE(FecTestCommon, AllTypes);

TYPED_TEST(FecTestCommon, TestNf4) // NOLINT
{
    const int iter_count = arith::log2<TypeParam>(sizeof(TypeParam));

    for (int i = 1; i < iter_count; i++) {
        const unsigned word_size = 1 << i;
        fec::RsNf4<TypeParam> fec(word_size, this->n_data, this->n_parities);

        this->run_test(fec, true, true);
    }
}

TYPED_TEST(FecTestCommon, TestGf2nFft) // NOLINT
{
    for (size_t wordsize = 1; wordsize <= sizeof(TypeParam); wordsize *= 2) {
        fec::RsGf2nFft<TypeParam> fec(wordsize, this->n_data, this->n_parities);

        this->run_test(fec);
    }
}

TYPED_TEST(FecTestCommon, TestGf2nFftAdd) // NOLINT
{
    for (size_t wordsize = 1; wordsize <= sizeof(TypeParam); wordsize *= 2) {
        fec::RsGf2nFftAdd<TypeParam> fec(
            wordsize, this->n_data, this->n_parities);

        this->run_test(fec);
    }
}

template <typename T>
class FecTestNo128 : public FecTestCommon<T> {
};

using No128 = ::testing::Types<uint32_t, uint64_t>;
TYPED_TEST_CASE(FecTestNo128, No128);

TYPED_TEST(FecTestNo128, TestFnt) // NOLINT
{
    for (unsigned word_size = 1; word_size <= 2; ++word_size) {
        fec::RsFnt<TypeParam> fec(
            fec::FecType::NON_SYSTEMATIC,
            word_size,
            this->n_data,
            this->n_parities);
        this->run_test(fec, true);
    }
}

TYPED_TEST(FecTestNo128, TestFntSys) // NOLINT
{
    for (unsigned word_size = 1; word_size <= 2; ++word_size) {
        fec::RsFnt<TypeParam> fec(
            fec::FecType::SYSTEMATIC,
            word_size,
            this->n_data,
            this->n_parities);
        this->run_test(fec, true);
    }
}

TYPED_TEST(FecTestNo128, TestGfpFft) // NOLINT
{
    for (size_t word_size = 1; word_size <= 4 && word_size < sizeof(TypeParam);
         word_size *= 2) {
        fec::RsGfpFft<TypeParam> fec(word_size, this->n_data, this->n_parities);

        this->run_test(fec, true);
    }
}
