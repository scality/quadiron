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
#include <random>

#include <gtest/gtest.h>

#include "quadiron.h"

namespace vec = quadiron::vec;
namespace gf = quadiron::gf;

template <typename T>
class BuffersTest : public ::testing::Test {
  public:
    std::unique_ptr<vec::Buffers<T>>
    gen_buffers_rand_data(int n, int size, int _max = 0)
    {
        std::mt19937 prng;
        T max_val = 65537;
        const int max = (_max == 0) ? max_val : _max;
        std::uniform_int_distribution<uint32_t> dis(0, max - 1);
        auto vec = std::make_unique<vec::Buffers<T>>(n, size);

        for (int i = 0; i < n; i++) {
            T* buf = quadiron::aligned_allocate<T>(size);
            for (int j = 0; j < size; j++) {
                buf[j] = dis(prng);
            }
            vec->set(i, buf);
        }

        return vec;
    }

    std::unique_ptr<vec::Vector<T>>
    gen_rand_vector(gf::RingModN<T>& gf, size_t n, T max_val)
    {
        std::vector<T> mem;
        mem.reserve(max_val);

        for (size_t i = 0; i < max_val; i++) {
            mem.push_back(i);
        }
        std::random_shuffle(mem.begin(), mem.end());

        auto vec = std::make_unique<vec::Vector<T>>(gf, n);
        for (size_t i = 0; i < n; i++) {
            vec->set(i, mem[i]);
        }
        return vec;
    }

    bool check_eq(const T* buf1, const T* buf2, size_t len)
    {
        return memcmp(buf1, buf2, len * sizeof(T)) == 0;
    }

    bool check_all_zeros(const T* buf, size_t len)
    {
        for (size_t i = 0; i < len; ++i) {
            if (buf[i] != 0) {
                return false;
            }
        }
        return true;
    }

    bool check_shuffled_bufs(
        const vec::Buffers<T>& input,
        const vec::Buffers<T>& output,
        const vec::Vector<T>& map)
    {
        const size_t input_len = input.get_n();
        const size_t output_len = output.get_n();
        const size_t map_len = map.get_n();

        const size_t size = input.get_size();

        std::vector<bool> check(output.get_n(), false);

        for (unsigned i = 0; i < map_len; ++i) {
            if (!check_eq(input.get(i), output.get(map.get(i)), size)) {
                return false;
            }
            check[map.get(i)] = true;
        }

        // Check zero-extended if it's necessary
        if (output_len > input_len) {
            for (size_t i = 0; i < output_len; ++i) {
                if (!check[i]) {
                    if (!check_all_zeros(output.get(i), size)) {
                        return false;
                    }
                }
            }
        }

        return true;
    }
};

using TestedTypes = ::testing::Types<uint32_t, uint64_t>;
TYPED_TEST_CASE(BuffersTest, TestedTypes);

TYPED_TEST(BuffersTest, TestConstructors) // NOLINT
{
    const int n = 16;
    const int begin = 5;
    const int end = 12;
    const int size = 4;

    auto vec1 = this->gen_buffers_rand_data(n, size);
    vec::Buffers<TypeParam> vec2(*vec1, begin, end);

    const std::vector<TypeParam*> mem1 = vec1->get_mem();
    const std::vector<TypeParam*> mem2 = vec2.get_mem();

    ASSERT_EQ(vec2.get_n(), end - begin);
    ASSERT_EQ(vec2.get_size(), vec1->get_size());

    for (int i = 0; i < end - begin; i++) {
        for (int j = 0; j < size; j++) {
            mem2.at(i)[j] = mem1.at(i + begin)[j];
        }
    }

    std::vector<TypeParam*> mem3(end - begin, nullptr);
    for (int i = 0; i < end - begin; i++) {
        mem3[i] = mem1.at(i + begin);
    }
    vec::Buffers<TypeParam> vec3(end - begin, size, mem3);

    ASSERT_EQ(vec2, vec3);

    auto gf(gf::create<gf::Prime<TypeParam>>(65537));

    // no-extension
    const int out_n_1 = n - 5;
    auto map_1 = this->gen_rand_vector(gf, out_n_1, out_n_1);
    vec::Buffers<TypeParam> vec4(*vec1, *map_1, out_n_1);
    ASSERT_TRUE(this->check_shuffled_bufs(*vec1, vec4, *map_1));

    // extension
    const int out_n_2 = n + 10;
    auto map_2 = this->gen_rand_vector(gf, n, out_n_2);
    vec::Buffers<TypeParam> vec5(*vec1, *map_2, out_n_2);
    ASSERT_TRUE(this->check_shuffled_bufs(*vec1, vec5, *map_2));
}

TYPED_TEST(BuffersTest, TestEvenOddSeparation) // NOLINT
{
    const int n = 8;
    const int size = 32;
    const int half = n / 2;
    auto vec1 = this->gen_buffers_rand_data(n, size);
    vec::Buffers<TypeParam> vec2(n, size);

    vec2.copy(*vec1);

    std::vector<TypeParam*>* even_mem =
        new std::vector<TypeParam*>(half, nullptr);
    std::vector<TypeParam*>* odd_mem =
        new std::vector<TypeParam*>(half, nullptr);
    vec::Buffers<TypeParam> i_even(half, size, *even_mem);
    vec::Buffers<TypeParam> i_odd(half, size, *odd_mem);
    vec1->separate_even_odd(i_even, i_odd);

    vec1->separate_even_odd();

    vec::Buffers<TypeParam> _i_even(*vec1, 0, half);
    vec::Buffers<TypeParam> _i_odd(*vec1, half, n);
    ASSERT_EQ(i_even, _i_even);
    ASSERT_EQ(i_odd, _i_odd);

    const std::vector<TypeParam*> mem1 = vec1->get_mem();
    const std::vector<TypeParam*> mem2 = vec2.get_mem();

    bool ok = true;
    for (int i = 0; i < n / 2; i += 2) {
        TypeParam* even1 = mem1.at(i);
        TypeParam* even2 = mem2.at(i * 2);
        TypeParam* odd1 = mem1.at(i + n / 2);
        TypeParam* odd2 = mem2.at(i * 2 + 1);
        for (int j = 0; j < size; j++) {
            if (even1[j] != even2[j] || odd1[j] != odd2[j]) {
                ok = false;
                i = n;
                break;
            }
        }
    }
    ASSERT_TRUE(ok);
}

TYPED_TEST(BuffersTest, TestZeroExtented) // NOLINT
{
    const int n = 8;
    const int size = 32;
    const int n1 = 4;
    const int n2 = 10;

    auto vec = this->gen_buffers_rand_data(n, size);
    vec::Buffers<TypeParam> vec1(*vec, n1);
    vec::Buffers<TypeParam> vec2(*vec, n2);

    vec::Buffers<TypeParam> _vec1(*vec, 0, n1);
    vec::Buffers<TypeParam> _vec2(*vec, n2);

    ASSERT_EQ(vec1, _vec1);
    ASSERT_EQ(vec2, _vec2);

    vec::Buffers<TypeParam> vec3(vec2, n1);
    ASSERT_EQ(vec3, vec1);
}

TYPED_TEST(BuffersTest, TestPackUnpack) // NOLINT
{
    const int iter_count = quadiron::arith::log2<TypeParam>(sizeof(TypeParam));

    for (int i = 0; i <= iter_count; i++) {
        const int word_size = quadiron::arith::exp2<TypeParam>(i);
        const int n = 8;
        const int size = 32;
        const int bytes_size = size * word_size;
        const TypeParam max = (static_cast<TypeParam>(1) << word_size) + 1;
        auto words = this->gen_buffers_rand_data(n, size, max);
        const std::vector<TypeParam*> mem_T = words->get_mem();

        // Pack manually from TypeParam to uint8_t.
        vec::Buffers<uint8_t> vec_char(n, bytes_size);
        const std::vector<uint8_t*> mem_char = vec_char.get_mem();
        for (int j = 0; j < n; j++) {
            int t = 0;
            TypeParam* buf_T = mem_T.at(j);
            uint8_t* buf_char = mem_char.at(j);

            for (int k = 0; k < size; k++) {
                TypeParam symb = buf_T[k];
                buf_char[t] = static_cast<uint8_t>(symb & 0xFF);

                t++;
                for (int u = 1; u < word_size; u++) {
                    symb >>= 8;
                    buf_char[t] = static_cast<uint8_t>(symb & 0xFF);
                    t++;
                }
            }
        }

        // Pack bufs of type uint8_t to bufs of type TypeParam.
        vec::Buffers<TypeParam> vec_T_tmp(n, size);
        const std::vector<TypeParam*> mem_T_tmp = vec_T_tmp.get_mem();
        vec::pack<uint8_t, TypeParam>(mem_char, mem_T_tmp, n, size, word_size);
        ASSERT_EQ(vec_T_tmp, *words);

        // Unpack bufs of type TypeParam to bufs of type uint8_t.
        vec::Buffers<uint8_t> vec_char_tmp(n, bytes_size);
        const std::vector<uint8_t*> mem_char_tmp = vec_char_tmp.get_mem();
        vec::unpack<TypeParam, uint8_t>(
            mem_T_tmp, mem_char_tmp, n, size, word_size);
        ASSERT_EQ(vec_char_tmp, vec_char);
    }
}
