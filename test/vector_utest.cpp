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

namespace gf = quadiron::gf;
namespace vec = quadiron::vec;

template <typename T>
class VectorTest : public ::testing::Test {
};

using TestedTypes = ::testing::Types<uint32_t, uint64_t>;
TYPED_TEST_CASE(VectorTest, TestedTypes);

TYPED_TEST(VectorTest, TestSlices) // NOLINT
{
    const auto gfp(gf::create<gf::Prime<TypeParam>>(65537));
    const int len = 20;
    const int len1 = 7;
    const int len2 = 5;
    const int offset1 = 2;
    const int offset2 = 3;
    const int len3 = (offset2 + len2 > len1) ? len1 - offset2 : len2;
    vec::Vector<TypeParam> base_vec(gfp, len);

    for (int i = 0; i < len; i++) {
        base_vec.set(i, gfp.rand());
    }
    // vmvec1 = base_vec[offset1, .., offset1 + len1 - 1]
    vec::Slice<TypeParam> vmvec1(&base_vec, len1, offset1);
    // vmvec2 = vmvec1[offset2, .., min(offset2 + len2 - 1, len1 - 1)]
    vec::Slice<TypeParam> vmvec2(&vmvec1, len2, offset2);
    // vmvec3 = base_vec[offset1 + offset2, .., offset1 + len1 - 1]
    vec::Slice<TypeParam> vmvec3(&base_vec, len3, offset1 + offset2);

    ASSERT_EQ(vmvec3, vmvec2);
}
