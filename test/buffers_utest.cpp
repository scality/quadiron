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
#include "simd/simd.h"

namespace vec = quadiron::vec;
namespace gf = quadiron::gf;

template <typename T>
bool is_all_zeros(const T* buf, size_t len)
{
    for (size_t i = 0; i < len; ++i) {
        if (buf[i] != 0) {
            return false;
        }
    }
    return true;
}

template <typename T>
class BuffersTest : public ::testing::Test {
  public:
    quadiron::simd::AlignedAllocator<T> allocator;
    quadiron::simd::AlignedAllocator<uint8_t> allocator_meta;

    BuffersTest()
    {
        quadiron::prng().seed(time(0));
    }
    ~BuffersTest() = default;

    std::unique_ptr<vec::Buffers<T>>
    gen_buffers_rand_data(int n, int size, bool has_meta = false, int _max = 0)
    {
        const T max = (_max == 0) ? std::numeric_limits<T>::max() : _max;
        std::uniform_int_distribution<T> dis(0, max - 1);
        auto vec = std::make_unique<vec::Buffers<T>>(n, size, has_meta);

        const std::vector<T*> mem = vec->get_mem();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < size; j++) {
                mem[i][j] = dis(quadiron::prng());
            }
        }

        if (has_meta) {
            const std::vector<uint8_t*> meta = vec->get_meta();
            const size_t meta_size = vec->get_meta_size();
            for (int i = 0; i < n; ++i) {
                for (size_t j = 0; j < meta_size; ++j) {
                    meta[i][j] = static_cast<uint8_t>(dis(quadiron::prng()));
                }
            }
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
            if (!std::equal(
                    input.get(i),
                    input.get(i) + size,
                    output.get(map.get(i)))) {
                return false;
            }
            if (input.has_meta()
                && !std::equal(
                       input.get_meta(i),
                       input.get_meta(i) + input.get_meta_size(),
                       output.get_meta(map.get(i)))) {
                return false;
            }
            check[map.get(i)] = true;
        }

        // Check zero-extended if it's necessary
        if (output_len > input_len) {
            for (size_t i = 0; i < output_len; ++i) {
                if (!check[i]) {
                    if (!is_all_zeros(output.get(i), size)) {
                        return false;
                    }
                    if (input.has_meta()
                        && !is_all_zeros(
                               output.get_meta(i), output.get_meta_size())) {
                        return false;
                    }
                }
            }
        }

        return true;
    }
};

using TestedTypes = ::testing::Types<uint8_t, uint16_t, uint32_t, uint64_t>;
TYPED_TEST_CASE(BuffersTest, TestedTypes);

TYPED_TEST(BuffersTest, TestConstructors) // NOLINT
{
    const std::vector<bool> tests = {true, false};

    for (bool const& has_meta : tests) {
        const int n = 16;
        const int begin = 5;
        const int end = 12;
        const int size = 4;

        auto vec1 = this->gen_buffers_rand_data(n, size, has_meta);
        ASSERT_EQ(has_meta, vec1->has_meta());

        const size_t meta_size = vec1->get_meta_size();

        vec::Buffers<TypeParam> vec2(*vec1, begin, end);

        const std::vector<TypeParam*> mem1 = vec1->get_mem();
        const std::vector<TypeParam*> mem2 = vec2.get_mem();

        // Check Slice constructor
        ASSERT_EQ(vec2.get_n(), end - begin);
        ASSERT_EQ(vec2.get_size(), vec1->get_size());
        for (int i = begin, j = 0; i < end; ++i, ++j) {
            ASSERT_TRUE(std::equal(mem1[i], mem1[i] + size, mem2[j]));
        }

        if (has_meta) {
            ASSERT_EQ(vec2.has_meta(), vec1->has_meta());
            ASSERT_EQ(vec2.get_meta_size(), meta_size);

            const std::vector<uint8_t*> meta1 = vec1->get_meta();
            const std::vector<uint8_t*> meta2 = vec2.get_meta();

            for (int i = begin, j = 0; i < end; ++i, ++j) {
                ASSERT_TRUE(
                    std::equal(meta1[i], meta1[i] + meta_size, meta2[j]));
            }
        }

        // Check constructor from given mem and meta
        std::vector<TypeParam*> mem3(end - begin);
        for (int i = 0; i < end - begin; i++) {
            mem3[i] = this->allocator.allocate(size);
            std::copy_n(mem1[i + begin], size, mem3[i]);
        }

        if (has_meta) {
            const std::vector<uint8_t*> meta1 = vec1->get_meta();
            std::vector<uint8_t*> meta3(end - begin);
            for (int i = 0; i < end - begin; i++) {
                meta3[i] = this->allocator_meta.allocate(meta_size);
                std::copy_n(meta1[i + begin], meta_size, meta3[i]);
            }

            vec::Buffers<TypeParam> vec3(end - begin, size, mem3, &meta3);

            ASSERT_EQ(vec2, vec3);

        } else {
            vec::Buffers<TypeParam> vec3(end - begin, size, mem3);

            ASSERT_EQ(vec2, vec3);
        }

        auto gf(gf::create<gf::Prime<TypeParam>>(static_cast<TypeParam>(31)));

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
}

TYPED_TEST(BuffersTest, TestZeroExtented) // NOLINT
{
    const int n = 8;
    const int size = 32;
    const int n1 = 4;
    const int n2 = 10;

    const std::vector<bool> tests = {true, false};

    for (bool const& has_meta : tests) {
        auto vec = this->gen_buffers_rand_data(n, size, has_meta);
        // cloned constructor: no zero-padding as `n1 < n`
        vec::Buffers<TypeParam> vec1(*vec, n1);
        // cloned constructor: there are zero-padding as `n2 > n`
        vec::Buffers<TypeParam> vec2(*vec, n2);

        // slice from `vec` as `n1 < n`
        vec::Buffers<TypeParam> _vec1(*vec, 0, n1);
        // slice and zero-extended from `vec` as `n2 > n`
        vec::Buffers<TypeParam> _vec2(*vec, 0, n2);

        ASSERT_EQ(vec1, _vec1);
        ASSERT_EQ(vec2, _vec2);

        vec::Buffers<TypeParam> vec3(vec2, n);
        ASSERT_EQ(vec3, *vec);
    }
}

TYPED_TEST(BuffersTest, TestPackUnpack) // NOLINT
{
    const int iter_count = quadiron::arith::log2<TypeParam>(sizeof(TypeParam));
    const std::vector<bool> tests = {true, false};

    for (int i = 0; i <= iter_count; i++) {
        const int word_size = quadiron::arith::exp2<TypeParam>(i);
        const int n = 8;
        const int size = 32;
        const int bytes_size = size * word_size;
        const TypeParam max = (static_cast<TypeParam>(1) << word_size) + 1;

        for (bool const& has_meta : tests) {
            auto words = this->gen_buffers_rand_data(n, size, has_meta, max);

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
                        symb = static_cast<TypeParam>(symb) >> 8;
                        buf_char[t] = static_cast<uint8_t>(symb & 0xFF);
                        t++;
                    }
                }
            }

            // Pack bufs of type uint8_t to bufs of type TypeParam.
            vec::Buffers<TypeParam> vec_T_tmp(n, size);
            const std::vector<TypeParam*> mem_T_tmp = vec_T_tmp.get_mem();

            vec::pack<uint8_t, TypeParam>(
                mem_char, mem_T_tmp, n, size, word_size);

            // Unpack bufs of type TypeParam to bufs of type uint8_t.
            vec::Buffers<uint8_t> vec_char_tmp(n, bytes_size);
            const std::vector<uint8_t*> mem_char_tmp = vec_char_tmp.get_mem();
            vec::unpack<TypeParam, uint8_t>(
                mem_T_tmp, mem_char_tmp, n, size, word_size);

            for (int i = 0; i < n; ++i) {
                ASSERT_TRUE(
                    std::equal(mem_T_tmp[i], mem_T_tmp[i] + size, mem_T[i]));
                ASSERT_TRUE(std::equal(
                    mem_char_tmp[i],
                    mem_char_tmp[i] + bytes_size,
                    mem_char[i]));
            }
        }
    }
}

TYPED_TEST(BuffersTest, TestCalculateSize) // NOLINT
{
    // Convential size depends on type of word, i.e. TypeParam, hence we use an
    // array for expected values.
    typedef struct tuple {
        size_t size;                // in TypeParam words
        size_t size_alignment;      // in bytes
        size_t meta_size_alignment; // in bytes
        size_t expected[4]; // 4 values for 4 TypeParams [uint8_t, uint16_t,
                            // uint32_t, uint64_t]
    } tuple;

    std::vector<tuple> cases = {
        // sse
        {12, 16, 2, {16, 16, 12, 12}},
        {25, 16, 2, {32, 32, 28, 26}},
        // avx
        {12, 32, 4, {32, 16, 16, 12}},
        {25, 32, 4, {32, 32, 32, 28}},
        // whatever
        {11, 15, 6, {240, 120, 60, 30}},
    };

    for (auto const& t : cases) {
        const size_t new_size = vec::Buffers<TypeParam>::get_conv_size(
            t.size, t.size_alignment, t.meta_size_alignment);
        const size_t id = quadiron::arith::log2<TypeParam>(sizeof(TypeParam));

        ASSERT_EQ(t.expected[id], new_size);
    }
}
