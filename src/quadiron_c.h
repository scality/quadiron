/* -*- mode: c++ -*- */
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
#ifndef __QUAD_QUADIRON_C_H__
#define __QUAD_QUADIRON_C_H__

#ifdef __cplusplus
extern "C" {
#endif

/** Create FNT FEC - This FEC is relatively complex because it requires storing
 * a metadata header.
 *
 * @param[in] word_size FNT only supports 1 or 2
 * @param[in] n_data number of data fragments
 * @param[in] n_parities number of parity fragments
 * @param[in] systematic if 1 then the code is systematic otherwise
 * non-systematic
 *
 * @return the FEC instance pointer
 */
struct QuadironFnt32*
quadiron_fnt32_new(int word_size, int n_data, int n_parities, int systematic);

/** Delete FEC
 *
 * @param[in,out] fecp the FEC instance pointer
 */
void quadiron_fnt32_delete(struct QuadironFnt32* fecp);

/** Return metadata size
 *
 * FNT requires to store specific information in headers and therefore caller
 * needs to preallocate metadata_size in blocks.
 *
 * @param[in] fecp the FEC instance
 * @param[in] block_size the metadata_size is computed acc/to the block_size
 *
 * @return the metadata size
 */
int quadiron_fnt32_get_metadata_size(
    struct QuadironFnt32* fecp,
    size_t block_size);

/** Encode blocks
 *
 * @param[in] fecp the FEC instance
 * @param[in] data must be exactly n_data
 * buffers must allocate block_size + metadata_size
 * @param[out] parity must be exactly n_outputs
 * - set entries to NULL when not wanted
 * - callers must allocate block_size + metadata_size when wanted
 * - n_outputs is n_parities if systematic, and n_data + n_parities if
 * non-systematic
 * @param[in] wanted_idxs array of length n_outputs indicating
 * the wish (value 1) or not (value 0) of parities
 * @param[in] block_size the block size in bytes
 *
 * @return 0 if encode succeeded, else -1
 */
int quadiron_fnt32_encode(
    struct QuadironFnt32* fecp,
    uint8_t** data,
    uint8_t** parity,
    int* wanted_idxs,
    size_t block_size);

/** Decode blocks
 *
 * @note For non-systematic codes parities must be provided as data and parities
 *
 * @param[in] fecp the FEC instance
 * @param[in,out] data must be exactly n_data
 * - callers must allocate block_size + metadata_size for both provided and
 * missing data
 * - it is not possible to choose not to decode a data (see reconstruct api for
 * that)
 * @param[in] parity must be exactly n_outputs
 * - set entries to NULL when missing
 * - n_outputs is n_parities if systematic, and n_data + n_parities if
 * non-systematic
 * - callers must allocate block_size + metadata_size for provided parities
 * @param[in] missing_idxs array of length code_len indicating
 * presence (value 1) or absence (value 0) of fragments (data and parities)
 * - code_len is n_data + n_parities
 * @param[in] block_size the block size in bytes
 *
 * @return 0 if decode succeeded, else -1
 */
int quadiron_fnt32_decode(
    struct QuadironFnt32* fecp,
    uint8_t** data,
    uint8_t** parity,
    int* missing_idxs,
    size_t block_size);

/** Reconstruct block
 *
 * @note For non-systematic codes parities must be provided as data and parities
 *
 * @param[in] fecp the FEC instance
 * @param[in,out] data must be exactly n_data
 * - set entries to NULL when missing
 * - callers must allocate block_size + metadata_size for provided and
 * destination_idx data
 * @param[in,out] parity must be exactly n_parities
 * - set entries to NULL when missing
 * - callers must allocate block_size + metadata_size
 * @param[in] missing_idxs array of missing_idxs of len code_len indicating
 * presence (value 1) or absence (value 0) of fragments (data and parities)
 * @param[in] destination_idx index of fragment to reconstruct (data or parity)
 * @param[in] block_size the block size in bytes
 *
 * @return 0 if reconstruct succeeded, else -1
 */
int quadiron_fnt32_reconstruct(
    struct QuadironFnt32* fecp,
    uint8_t** data,
    uint8_t** parity,
    int* missing_idxs,
    unsigned int destination_idx,
    size_t block_size);

/** Dump a buffer on stderr (debug function)
 *
 * @param[in] buf the buffer
 * @param[in] size the buffer size
 */
void quadiron_hex_dump(uint8_t* buf, size_t size);

#ifdef __cplusplus
}
#endif

#endif
