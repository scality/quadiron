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
#include "quadiron_c.h"

template <typename T>
class QuadironCTest : public ::testing::Test {
  public:
    QuadironCTest()
    {
        quadiron::prng().seed(time(0));
    }
    ~QuadironCTest() = default;

    void randomize_buffer(uint8_t* buf, size_t size)
    {
        std::uniform_int_distribution<> dis(0, 256);
        for (u_int i = 0; i < size; i++)
            buf[i] = dis(quadiron::prng());
    }

    /** Generate all combinations of n given k
     *
     * @param n numbers of items
     * @param k number of elements in combinations
     *
     * @return the vector of combinations
     */
    std::vector<std::vector<int>> generate_combinations(int n, int k)
    {
        std::vector<std::vector<int>> combinations;
        std::vector<int> selected;
        std::vector<int> selector(n);
        std::fill(selector.begin(), selector.begin() + k, 1);
        do {
            for (int i = 0; i < n; i++) {
                if (selector[i]) {
                    selected.push_back(i);
                }
            }
            combinations.push_back(selected);
            selected.clear();
        } while (std::prev_permutation(selector.begin(), selector.end()));
        return combinations;
    }

    /** Convert list of indexes to list of boolean values
     *
     * @param src input idx_list
     * @param n_src number of entries in src
     * @param dst output idx list as boolean values (1/0)
     * @param n_dst number of entries in dst
     */
    void convert_idx_list(
        std::vector<int> src,
        int n_src,
        std::vector<int>& dst,
        int n_dst)
    {
        std::fill_n(dst.begin(), n_dst, 0);
        for (int i = 0; i < n_src; i++) {
            ASSERT_LT(src[i], n_dst);
            dst[src[i]] = 1;
        }
    }

    /** Test encode/decode/reconstruct
     *
     * This test will create n_data fragments and generate n_parities parities,
     * then decode with a set of missing fragments. Finaly reconstruct them.
     *
     * @param n_data number of data
     * @param n_parities number of parities
     * @param block_size size of block in bytes
     * @param systematic 1 if systematic else 0
     * @param missing_idxs vector of boolean vales indicating missing fragments
     * for decode and reconstruct
     *        - must be of length n_parities
     */
    void test_encode_decode_reconstruct(
        int n_data,
        int n_parities,
        size_t block_size,
        int systematic,
        std::vector<int> missing_idxs)
    {
        struct QuadironFnt32* inst =
            quadiron_fnt32_new(2, n_data, n_parities, systematic);
        size_t metadata_size =
            quadiron_fnt32_get_metadata_size(inst, block_size);
        std::vector<std::vector<uint8_t>> data(n_data);
        std::vector<uint8_t*> _data(n_data); // for C API
        std::vector<std::vector<uint8_t>> ref_data(n_data);
        std::vector<uint8_t*> _ref_data(n_data); // for C API
        int n_outputs;
        if (systematic) {
            n_outputs = n_parities;
        } else {
            n_outputs = n_data + n_parities;
        }
        std::vector<std::vector<uint8_t>> parity(n_outputs);
        std::vector<uint8_t*> _parity(n_outputs); // for C API
        std::vector<int> wanted_idxs(n_outputs);

        for (int i = 0; i < n_data; i++) {
            data.at(i).resize(block_size + metadata_size);
            _data[i] = data.at(i).data();
            ref_data.at(i).resize(block_size);
            _ref_data[i] = ref_data.at(i).data();
            randomize_buffer(_data[i] + metadata_size, block_size);
            std::copy_n(_data[i] + metadata_size, block_size, _ref_data[i]);
        }

        for (int i = 0; i < n_outputs; i++) {
            parity.at(i).resize(block_size + metadata_size);
            _parity[i] = parity.at(i).data();
        }

        // we want all parities
        std::fill_n(wanted_idxs.begin(), n_outputs, 1);

        ASSERT_EQ(
            quadiron_fnt32_encode(
                inst,
                _data.data(),
                _parity.data(),
                wanted_idxs.data(),
                block_size),
            0);

        for (int i = 0; i < n_data; i++) {
            if (missing_idxs[i]) {
                std::fill_n(_data[i], block_size + metadata_size, 0);
            }
        }

        for (int i = 0; i < n_parities; i++) {
            if (missing_idxs[n_data + i]) {
                std::fill_n(_parity[i], block_size + metadata_size, 0);
            }
        }

        for (int i = 0; i < n_data; i++) {
            if (missing_idxs[i]) {
                ASSERT_EQ(
                    quadiron_fnt32_reconstruct(
                        inst,
                        _data.data(),
                        _parity.data(),
                        missing_idxs.data(),
                        i,
                        block_size),
                    0);
            }
        }

        for (int i = 0; i < n_parities; i++) {
            if (missing_idxs[n_data + i]) {
                ASSERT_EQ(
                    quadiron_fnt32_reconstruct(
                        inst,
                        _data.data(),
                        _parity.data(),
                        missing_idxs.data(),
                        n_data + i,
                        block_size),
                    0);
            }
        }

        for (int i = 0; i < n_data; i++) {
            if (missing_idxs[i]) {
                std::fill_n(_data[i], block_size + metadata_size, 0);
            }
        }

        for (int i = 0; i < n_parities; i++) {
            if (missing_idxs[n_data + i]) {
                std::fill_n(_parity[i], block_size + metadata_size, 0);
            }
        }

        // FIXME: for non-systematic FNT, `quadiron_fnt32_decode` will
        // overwrite `_data` (that stores actually encoded fragments)
        // by decoded data. Hence, if we test `reconstruct` after, we could use
        // wrong fragments. A quick fix is to move the test of
        // `quadiron_fnt32_decode` at the end.
        ASSERT_EQ(
            quadiron_fnt32_decode(
                inst,
                _data.data(),
                _parity.data(),
                missing_idxs.data(),
                block_size),
            0);

        for (int i = 0; i < n_data; i++) {
            ASSERT_TRUE(std::equal(
                _ref_data[i],
                _ref_data[i] + block_size,
                _data[i] + metadata_size));
        }

        quadiron_fnt32_delete(inst);
    }

    void test_all_decodable_scenarios(int k, int m, int systematic)
    {
        for (int i = 0; i <= m; i++) {
            const auto combinations = generate_combinations(k + m, i);
            for (auto it = combinations.begin(); it != combinations.end();
                 it++) {
                std::vector<int> missing_idxs(k + m);
                convert_idx_list(*it, i, missing_idxs, k + m);
                test_encode_decode_reconstruct(
                    k, m, 10000, systematic, missing_idxs);
            }
        }
    }
};

using AllTypes = ::testing::Types<uint32_t>;
TYPED_TEST_CASE(QuadironCTest, AllTypes);

TYPED_TEST(QuadironCTest, TestBasicSys) // NOLINT
{
    this->test_all_decodable_scenarios(3, 3, 1);
}

TYPED_TEST(QuadironCTest, TestBasicNSys) // NOLINT
{
    this->test_all_decodable_scenarios(3, 3, 0);
}
