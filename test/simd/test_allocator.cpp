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
#include <cstdint>
#include <limits>

#include <gtest/gtest.h>

#include "simd/allocator.h"

namespace simd = quadiron::simd;

class SimdAllocatorTest : public ::testing::Test {
  public:
    simd::AlignedAllocator<int> allocator;
};

// Check that we gracefully handle zero-sized allocation.
TEST_F(SimdAllocatorTest, TestEmptyAlloc) // NOLINT
{
    int* ptr = this->allocator.allocate(0);

    ASSERT_TRUE(simd::addr_is_aligned(ptr));
    this->allocator.deallocate(ptr, 0);
}

// Check that we always return aligned memory.
TEST_F(SimdAllocatorTest, TestAlignment) // NOLINT
{
    for (int i = 1; i < 1000; ++i) {
        int* ptr = this->allocator.allocate(i);

        ASSERT_TRUE(simd::addr_is_aligned(ptr));
        this->allocator.deallocate(ptr, i);
    }
}

// Check that we gracefully handle zero-sized allocation.
TEST_F(SimdAllocatorTest, TestFreeNullPtr) // NOLINT
{
    // Like the default allocator, we should't crash on nullptr.
    this->allocator.deallocate(nullptr, 0);
}

// Check that we properly handle overflow.
TEST_F(SimdAllocatorTest, TestAllocateTooMuch) // NOLINT
{
    ASSERT_THROW(
        this->allocator.allocate(std::numeric_limits<std::size_t>::max()),
        std::bad_alloc);
}
