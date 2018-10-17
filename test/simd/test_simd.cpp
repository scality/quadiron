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
#include <vector>

#include <gtest/gtest.h>

#include "simd/simd.h"

namespace simd = quadiron::simd;

TEST(SimdTest, TestCountOf) // NOLINT
{
    std::vector<std::size_t> expected;

    switch (simd::INSTRUCTION_SET) {
    case simd::InstructionSet::NONE:
        expected = {1, 1, 1, 1};
        break;
    case simd::InstructionSet::SSE:
        expected = {16, 8, 4, 2};
        break;
    case simd::InstructionSet::AVX:
        expected = {32, 16, 8, 4};
        break;
    }

    ASSERT_EQ(simd::countof<uint8_t>(), expected[0]);
    ASSERT_EQ(simd::countof<uint16_t>(), expected[1]);
    ASSERT_EQ(simd::countof<uint32_t>(), expected[2]);
    ASSERT_EQ(simd::countof<uint64_t>(), expected[3]);
}
