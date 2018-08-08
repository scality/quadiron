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

TEST(MatrixTest, TestInvertGf11) // NOLINT
{
    const quadiron::gf::Prime<uint32_t> gf11(11);
    quadiron::vec::Matrix<uint32_t> mat(gf11, 3, 3);

    mat.set(0, 0, 2);
    mat.set(0, 1, 1);
    mat.set(0, 2, 2);
    mat.set(1, 0, 1);
    mat.set(1, 1, 2);
    mat.set(1, 2, 9);
    mat.set(2, 0, 1);
    mat.set(2, 1, 2);
    mat.set(2, 2, 7);

    mat.inv();

    ASSERT_EQ(mat.get(0, 0), 8);
    ASSERT_EQ(mat.get(0, 1), 6);
    ASSERT_EQ(mat.get(0, 2), 1);
    ASSERT_EQ(mat.get(1, 0), 7);
    ASSERT_EQ(mat.get(1, 1), 9);
    ASSERT_EQ(mat.get(1, 2), 10);
    ASSERT_EQ(mat.get(2, 0), 0);
    ASSERT_EQ(mat.get(2, 1), 6);
    ASSERT_EQ(mat.get(2, 2), 5);
}

TEST(MatrixTest, TestInvertGf29) // NOLINT
{
    const quadiron::gf::Prime<uint32_t> gf29(29);
    quadiron::vec::Matrix<uint32_t> mat(gf29, 3, 3);

    mat.set(0, 0, 22);
    mat.set(0, 1, 27);
    mat.set(0, 2, 18);
    mat.set(1, 0, 18);
    mat.set(1, 1, 28);
    mat.set(1, 2, 5);
    mat.set(2, 0, 4);
    mat.set(2, 1, 17);
    mat.set(2, 2, 1);

    mat.inv();

    ASSERT_EQ(mat.get(0, 0), 1);
    ASSERT_EQ(mat.get(0, 1), 18);
    ASSERT_EQ(mat.get(0, 2), 8);
    ASSERT_EQ(mat.get(1, 0), 2);
    ASSERT_EQ(mat.get(1, 1), 8);
    ASSERT_EQ(mat.get(1, 2), 11);
    ASSERT_EQ(mat.get(2, 0), 20);
    ASSERT_EQ(mat.get(2, 1), 24);
    ASSERT_EQ(mat.get(2, 2), 14);
}
