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
#include <array>

#include <gtest/gtest.h>

#include "quadiron.h"

template <typename T>
class VectorTest : public ::testing::Test {
};

using TestedTypes = ::testing::Types<uint32_t, uint64_t>;
TYPED_TEST_CASE(VectorTest, TestedTypes);

TYPED_TEST(VectorTest, TestHadamardMultiplication) // NOLINT
{
    const quadiron::gf::Prime<TypeParam> gfp(65537);
    const std::array<TypeParam, 16> vec1_data = {1,
                                                 64,
                                                 4096,
                                                 65533,
                                                 65281,
                                                 49153,
                                                 16,
                                                 1024,
                                                 65536,
                                                 65473,
                                                 61441,
                                                 4,
                                                 256,
                                                 16384,
                                                 65521,
                                                 64513};
    const std::array<TypeParam, 8> vec2_data = {
        6505, 23324, 7736, 43678, 27296, 50697, 59274, 48649};
    const std::array<TypeParam, 16> expected = {
        6505,
        50922,
        32285,
        21899,
        24683,
        61827,
        30866,
        8456,
        59032,
        14615,
        33252,
        43638,
        40854,
        3710,
        34671,
        57081,
    };
    quadiron::vec::Vector<TypeParam> vec1(gfp, 16);
    quadiron::vec::Vector<TypeParam> vec2(gfp, 8);
    quadiron::vec::Doubled<TypeParam> v2vec2(&vec2);

    for (unsigned i = 0; i < vec1_data.size(); ++i) {
        vec1.set(i, vec1_data[i]);
    }
    for (unsigned i = 0; i < vec2_data.size(); ++i) {
        vec2.set(i, vec2_data[i]);
    }

    vec1.hadamard_mul(&v2vec2);

    for (unsigned i = 0; i < expected.size(); ++i) {
        ASSERT_EQ(vec1.get(i), expected[i]);
    }
}

TYPED_TEST(VectorTest, TestAddition) // NOLINT
{
    const quadiron::gf::Prime<TypeParam> gfp(65537);
    const std::array<TypeParam, 8> vec1_data = {
        5459, 11947, 44310, 21807, 60078, 53590, 21227, 43730};
    const std::array<TypeParam, 4> vec2_data = {39466, 40329, 16012, 15149};
    const std::array<TypeParam, 8> expected = {
        44925,
        52276,
        60322,
        36956,
        34007,
        28382,
        37239,
        58879,
    };
    quadiron::vec::Vector<TypeParam> vec1(gfp, 8);
    quadiron::vec::Vector<TypeParam> vec2(gfp, 4);
    quadiron::vec::Doubled<TypeParam> v2vec2(&vec2);

    for (unsigned i = 0; i < vec1_data.size(); ++i) {
        vec1.set(i, vec1_data[i]);
    }
    for (unsigned i = 0; i < vec2_data.size(); ++i) {
        vec2.set(i, vec2_data[i]);
    }

    vec1.add(&v2vec2);

    for (unsigned i = 0; i < expected.size(); ++i) {
        ASSERT_EQ(vec1.get(i), expected[i]);
    }
}

TYPED_TEST(VectorTest, TestSlices) // NOLINT
{
    const quadiron::gf::Prime<TypeParam> gfp(65537);
    const int len = 20;
    const int len1 = 7;
    const int len2 = 5;
    const int offset1 = 2;
    const int offset2 = 3;
    const int len3 = (offset2 + len2 > len1) ? len1 - offset2 : len2;
    quadiron::vec::Vector<TypeParam> base_vec(gfp, len);

    for (int i = 0; i < len; i++) {
        base_vec.set(i, gfp.weak_rand());
    }
    // vmvec1 = base_vec[offset1, .., offset1 + len1 - 1]
    quadiron::vec::Slice<TypeParam> vmvec1(&base_vec, len1, offset1);
    // vmvec2 = vmvec1[offset2, .., min(offset2 + len2 - 1, len1 - 1)]
    quadiron::vec::Slice<TypeParam> vmvec2(&vmvec1, len2, offset2);
    // vmvec3 = base_vec[offset1 + offset2, .., offset1 + len1 - 1]
    quadiron::vec::Slice<TypeParam> vmvec3(&base_vec, len3, offset1 + offset2);

    ASSERT_TRUE(vmvec3.eq(&vmvec2));
}
