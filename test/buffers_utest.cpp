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
#include <math.h>
#include <limits>

#include <gtest/gtest.h>

#include "quadiron.h"
#include "simd/simd.h"

namespace vec = quadiron::vec;
namespace gf = quadiron::gf;

template <typename T>
bool is_equal(const T* buf1, const T* buf2, size_t len)
{
    return memcmp(buf1, buf2, len * sizeof(T)) == 0;
}

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
std::unique_ptr<vec::Buffers<T>> help_gen_buffers_rand_data(
    quadiron::simd::AlignedAllocator<T>& allocator,
    quadiron::simd::AlignedAllocator<char>& char_allocator,
    int n,
    int size,
    int _max = 0)
{
    std::mt19937 prng;
    T max_val = std::numeric_limits<T>::max();
    const int max = (_max == 0) ? max_val : _max;
    std::uniform_int_distribution<uint32_t> dis(0, max - 1);
    auto vec = std::make_unique<vec::Buffers<T>>(n, size);

    // by default, each element refers to two words
    const int words_per_element = 2;
    const size_t bmap_size = ceil(size * words_per_element / 8);
    std::uniform_int_distribution<char> bmap_dis(-128, 127);

    for (int i = 0; i < n; i++) {
        T* buf = allocator.allocate(size);
        for (int j = 0; j < size; j++) {
            buf[j] = dis(prng);
        }
        vec->set(i, buf);
        char* bit_map = char_allocator.allocate(bmap_size);
        for (int j = 0; j < bmap_size; j++) {
            bit_map[j] = bmap_dis(prng);
        }
        vec->set_bit(i, bit_map);
    }

    return vec;
}

template <typename T>
class BuffersTest : public ::testing::Test {
  public:
    quadiron::simd::AlignedAllocator<T> allocator;
    quadiron::simd::AlignedAllocator<char> bmap_allocator;

    std::unique_ptr<vec::Buffers<T>>
    gen_buffers_rand_data(int n, int size, int _max = 0)
    {
        return help_gen_buffers_rand_data(
            allocator, bmap_allocator, n, size, _max);
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
        const size_t bmap_size = input.get_bmap_size();

        std::vector<bool> check(output.get_n(), false);

        for (unsigned i = 0; i < map_len; ++i) {
            if (!is_equal(input.get(i), output.get(map.get(i)), size)
                || !is_equal(
                       input.get_bit(i),
                       output.get_bit(map.get(i)),
                       bmap_size)) {
                return false;
            }
            check[map.get(i)] = true;
        }

        // Check zero-extended if it's necessary
        if (output_len > input_len) {
            for (size_t i = 0; i < output_len; ++i) {
                if (!check[i]) {
                    if (!is_all_zeros(output.get(i), size)
                        || !is_all_zeros(output.get_bit(i), bmap_size)) {
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

TYPED_TEST(BuffersTest, TestSliceConstructor) // NOLINT
{
    const int n = 16;
    const int begin = 5;
    const int end = 20;
    const int size = 4;

    // a Buffers with random data
    auto vec1 = this->gen_buffers_rand_data(n, size);
    const size_t bmap_size = vec1->get_bmap_size();

    // sliced Buffers
    vec::Buffers<TypeParam> vec2(*vec1, begin, end);

    const std::vector<TypeParam*> mem1 = vec1->get_mem();
    const std::vector<TypeParam*> mem2 = vec2.get_mem();

    const std::vector<char*> bmap1 = vec1->get_bmap();
    const std::vector<char*> bmap2 = vec2.get_bmap();


    ASSERT_EQ(vec2.get_n(), end - begin);
    ASSERT_EQ(vec2.get_size(), vec1->get_size());
    ASSERT_EQ(vec2.get_bmap_size(), vec1->get_bmap_size());

    // check common parts
    int i;
    for (i = 0; i < n - begin; ++i) {
        ASSERT_TRUE(
            is_equal(mem1[i + begin], mem2[i], size)
            && is_equal(bmap1[i + begin], bmap2[i], bmap_size));
    }
    // check padding parts
    for (; i < end - begin; ++i) {
        ASSERT_TRUE(
            is_all_zeros(mem2[i], size) && is_all_zeros(bmap2[i], bmap_size));
    }
}

TYPED_TEST(BuffersTest, TestCombinedConstructor) // NOLINT
{
    const int n1 = 16;
    const int n2 = 8;
    const int size = 4;

    // first Buffers with random data of length n1
    auto vec1 = this->gen_buffers_rand_data(n1, size);

    // second Buffers with random data of length n2
    auto vec2 = this->gen_buffers_rand_data(n2, size);
    // a Buffers combined from vec1 & vec2
    vec::Buffers<TypeParam> vec3(*vec1, *vec2);

    // a Buffers sliced of vec3 corresponding to vec1
    vec::Buffers<TypeParam> _vec1(vec3, 0, n1);
    ASSERT_EQ(*vec1, _vec1);

    // a Buffers sliced of vec3 corresponding to vec2
    vec::Buffers<TypeParam> _vec2(vec3, n1, n1 + n2);
    ASSERT_EQ(*vec2, _vec2);
}

TYPED_TEST(BuffersTest, TestMapConstructor) // NOLINT
{
    const int n = 16;
    const int size = 4;

    // a Buffers with random data of length n
    auto vec = this->gen_buffers_rand_data(n, size);

    auto gf(gf::create<gf::Prime<TypeParam>>(65537));

    // no-extension
    const int out_n_1 = n - 5;
    auto map_1 = this->gen_rand_vector(gf, out_n_1, out_n_1);
    vec::Buffers<TypeParam> vec1(*vec, *map_1, out_n_1);
    ASSERT_TRUE(this->check_shuffled_bufs(*vec, vec1, *map_1));

    // extension
    const int out_n_2 = n + 10;
    auto map_2 = this->gen_rand_vector(gf, n, out_n_2);
    vec::Buffers<TypeParam> vec2(*vec, *map_2, out_n_2);
    ASSERT_TRUE(this->check_shuffled_bufs(*vec, vec2, *map_2));
}

TYPED_TEST(BuffersTest, TestCastConstructor) // NOLINT
{
    const int n = 8;
    const int size = 32;
    const int size_T = size / sizeof(TypeParam);

    auto vec1 =
        help_gen_buffers_rand_data(this->bmap_allocator, this->bmap_allocator, n, size);
    const size_t bmap_size = vec1->get_bmap_size();

    vec::Buffers<TypeParam> vec2(*vec1);

    const std::vector<char*> mem1 = vec1->get_mem();
    const std::vector<TypeParam*> mem2 = vec2.get_mem();

    const std::vector<char*> bmap1 = vec1->get_bmap();
    const std::vector<char*> bmap2 = vec2.get_bmap();

    for (int i = 0; i < n; ++i) {
        TypeParam* _mem1_at_i = reinterpret_cast<TypeParam*>(mem1[i]);

        ASSERT_TRUE(
            is_equal(_mem1_at_i, mem2[i], size_T)
            && is_equal(bmap1[i], bmap2[i], bmap_size));
    }
}
