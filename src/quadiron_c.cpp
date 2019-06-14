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
#include "property.h"
#include "quadiron.h"
#include "quadiron_c.h"

extern "C" {

struct QuadironFnt32*
quadiron_fnt32_new(int word_size, int n_data, int n_parities, int systematic)
{
    const size_t pkt_size = 1024;

    if (word_size == 1 || word_size == 2) {
        return reinterpret_cast<struct QuadironFnt32*>(
            new quadiron::fec::RsFnt<uint32_t>(
                systematic ? quadiron::fec::FecType::SYSTEMATIC
                           : quadiron::fec::FecType::NON_SYSTEMATIC,
                word_size,
                n_data,
                n_parities,
                pkt_size));
    }

    return nullptr;
}

void quadiron_fnt32_delete(struct QuadironFnt32* fecp)
{
    delete reinterpret_cast<quadiron::fec::RsFnt<uint32_t>*>(fecp);
}

int quadiron_fnt32_get_metadata_size(
    struct QuadironFnt32* /* fecp */,
    size_t block_size)
{
    /*
     * We assume that a special value of 65536 may occur uniformly.
     * We count 4 bytes per special value.
     * We see large and roundup by 16 items.
     */
    return ((block_size / 65536) + 16) * 4;
}

int quadiron_fnt32_encode(
    struct QuadironFnt32* fecp,
    uint8_t** data,
    uint8_t** parity,
    int* wanted_idxs,
    size_t block_size)
{
    quadiron::fec::RsFnt<uint32_t>* fec =
        reinterpret_cast<quadiron::fec::RsFnt<uint32_t>*>(fecp);
    std::vector<uint8_t*> data_vec(fec->n_data, nullptr);
    std::vector<uint8_t*> parities_vec(fec->n_outputs, nullptr);
    std::vector<quadiron::Properties> parities_props(fec->n_outputs);
    std::vector<bool> wanted_idxs_vec(fec->n_outputs);
    int metadata_size = quadiron_fnt32_get_metadata_size(fecp, block_size);

    for (unsigned i = 0; i < fec->n_outputs; i++) {
        wanted_idxs_vec[i] = wanted_idxs[i] ? true : false;
    }

    if (fec->type == quadiron::fec::FecType::SYSTEMATIC) {
        for (unsigned i = 0; i < fec->n_data; i++) {
            data_vec[i] = data[i] + metadata_size;
        }
        for (unsigned i = 0; i < fec->n_parities; i++) {
            parities_vec[i] = parity[i] + metadata_size;
        }
    } else {
        for (unsigned i = 0; i < fec->n_data; i++) {
            data_vec[i] = data[i] + metadata_size;
            parities_vec[i] = data[i] + metadata_size;
        }
        for (unsigned i = 0; i < fec->n_parities; i++) {
            parities_vec[fec->n_data + i] = parity[i] + metadata_size;
        }
    }

    fec->encode_blocks_vertical(
        data_vec, parities_vec, parities_props, wanted_idxs_vec, block_size);

    if (fec->type == quadiron::fec::FecType::SYSTEMATIC) {
        quadiron::Properties null_prop;

        for (unsigned i = 0; i < fec->n_data; i++) {
            uint32_t* metadata = reinterpret_cast<uint32_t*>(data[i]);
            int ret = null_prop.fnt_serialize(metadata, metadata_size / 4);
            if (ret == -1) {
                return -1;
            }
        }
        for (unsigned i = 0; i < fec->n_parities; i++) {
            uint32_t* metadata = reinterpret_cast<uint32_t*>(parity[i]);
            int ret =
                parities_props[i].fnt_serialize(metadata, metadata_size / 4);
            if (ret == -1) {
                return -1;
            }
        }
    } else {
        for (unsigned i = 0; i < fec->n_data; i++) {
            uint32_t* metadata = reinterpret_cast<uint32_t*>(data[i]);
            int ret =
                parities_props[i].fnt_serialize(metadata, metadata_size / 4);
            if (ret == -1) {
                return -1;
            }
        }
        for (unsigned i = 0; i < fec->n_parities; i++) {
            uint32_t* metadata = reinterpret_cast<uint32_t*>(parity[i]);
            int ret = parities_props[fec->n_data + i].fnt_serialize(
                metadata, metadata_size / 4);
            if (ret == -1) {
                return -1;
            }
        }
    }

    return 0;
}

int quadiron_fnt32_decode(
    struct QuadironFnt32* fecp,
    uint8_t** data,
    uint8_t** parity,
    int* missing_idxs,
    size_t block_size)
{
    quadiron::fec::RsFnt<uint32_t>* fec =
        reinterpret_cast<quadiron::fec::RsFnt<uint32_t>*>(fecp);
    std::vector<uint8_t*> data_vec(fec->n_data, nullptr);
    std::vector<uint8_t*> parities_vec(fec->n_outputs, nullptr);
    std::vector<quadiron::Properties> parities_props(fec->n_outputs);
    std::vector<int> missing_idxs_vec(
        missing_idxs, missing_idxs + fec->code_len);
    std::vector<bool> wanted_idxs_vec(fec->n_data, true);
    int metadata_size = quadiron_fnt32_get_metadata_size(fecp, block_size);
    bool res;

    if (fec->type == quadiron::fec::FecType::SYSTEMATIC) {
        for (unsigned i = 0; i < fec->n_data; i++) {
            data_vec[i] = data[i] + metadata_size;
        }
        for (unsigned i = 0; i < fec->n_parities; i++) {
            if (!missing_idxs[fec->n_data + i]) {
                parities_vec[i] = parity[i] + metadata_size;
                uint32_t* metadata = reinterpret_cast<uint32_t*>(parity[i]);
                int ret = parities_props[i].fnt_deserialize(
                    metadata, metadata_size / 4);
                if (ret == -1)
                    return -1;
            }
        }
    } else {
        for (unsigned i = 0; i < fec->n_data; i++) {
            if (!missing_idxs[i]) {
                parities_vec[i] = data[i] + metadata_size;
                uint32_t* metadata = reinterpret_cast<uint32_t*>(data[i]);
                int ret = parities_props[i].fnt_deserialize(
                    metadata, metadata_size / 4);
                if (ret == -1)
                    return -1;
            }
            data_vec[i] = data[i] + metadata_size;
        }
        for (unsigned i = 0; i < fec->n_parities; i++) {
            if (!missing_idxs[fec->n_data + i]) {
                parities_vec[fec->n_data + i] = parity[i] + metadata_size;
                uint32_t* metadata = reinterpret_cast<uint32_t*>(parity[i]);
                int ret = parities_props[fec->n_data + i].fnt_deserialize(
                    metadata, metadata_size / 4);
                if (ret == -1)
                    return -1;
            }
        }
    }

    res = fec->decode_blocks_vertical(
        data_vec,
        parities_vec,
        parities_props,
        missing_idxs_vec,
        wanted_idxs_vec,
        block_size);
    if (!res)
        return -1;

    // reset metadata of data
    for (unsigned i = 0; i < fec->n_data; i++) {
        parities_props[i].clear();
        uint32_t* metadata = reinterpret_cast<uint32_t*>(data[i]);
        int ret = parities_props[i].fnt_serialize(metadata, metadata_size / 4);
        if (ret == -1) {
            return -1;
        }
    }

    return 0;
}

int quadiron_fnt32_reconstruct(
    struct QuadironFnt32* fecp,
    uint8_t** data,
    uint8_t** parity,
    int* missing_idxs,
    unsigned int destination_idx,
    size_t block_size)
{
    quadiron::fec::RsFnt<uint32_t>* fec =
        reinterpret_cast<quadiron::fec::RsFnt<uint32_t>*>(fecp);
    std::vector<uint8_t*> data_vec(fec->n_data);
    std::vector<uint8_t*> parities_vec(fec->n_outputs);
    std::vector<quadiron::Properties> parities_props(fec->n_outputs);
    std::vector<int> missing_idxs_vec(
        missing_idxs, missing_idxs + fec->code_len);
    std::vector<bool> wanted_data_vec(fec->n_data, false);
    std::vector<bool> wanted_idxs_vec(fec->n_outputs, false);
    int metadata_size = quadiron_fnt32_get_metadata_size(fecp, block_size);
    bool res;

    if (fec->type == quadiron::fec::FecType::SYSTEMATIC) {
        for (unsigned i = 0; i < fec->n_data; i++) {
            data_vec[i] = data[i] + metadata_size;
        }
        for (unsigned i = 0; i < fec->n_parities; i++) {
            if (!missing_idxs[fec->n_data + i]) {
                uint32_t* metadata = reinterpret_cast<uint32_t*>(parity[i]);
                int ret = parities_props[i].fnt_deserialize(
                    metadata, metadata_size / 4);
                if (ret == -1)
                    return -1;
            }
            parities_vec[i] = parity[i] + metadata_size;
        }
    } else {
        for (unsigned i = 0; i < fec->n_data; i++) {
            if (!missing_idxs[i]) {
                uint32_t* metadata = reinterpret_cast<uint32_t*>(data[i]);
                int ret = parities_props[i].fnt_deserialize(
                    metadata, metadata_size / 4);
                if (ret == -1)
                    return -1;
            }
            parities_vec[i] = data[i] + metadata_size;
        }
        for (unsigned i = 0; i < fec->n_parities; i++) {
            if (!missing_idxs[fec->n_data + i]) {
                uint32_t* metadata = reinterpret_cast<uint32_t*>(parity[i]);
                int ret = parities_props[fec->n_data + i].fnt_deserialize(
                    metadata, metadata_size / 4);
                if (ret == -1)
                    return -1;
            }
            parities_vec[fec->n_data + i] = parity[i] + metadata_size;
        }
    }

    int need_decode = 0;

    if (fec->type == quadiron::fec::FecType::SYSTEMATIC) {
        /*
         * Easy case where the target is a data then simply decode
         */
        if (destination_idx < fec->n_data) {
            wanted_idxs_vec[destination_idx] = true;
            res = fec->decode_blocks_vertical(
                data_vec,
                parities_vec,
                parities_props,
                missing_idxs_vec,
                wanted_idxs_vec,
                block_size);
            if (!res) {
                return -1;
            }
            for (unsigned i = 0; i < fec->n_data; i++) {
                quadiron::Properties null_prop;

                if (i == destination_idx) {
                    uint32_t* metadata = reinterpret_cast<uint32_t*>(data[i]);
                    int ret =
                        null_prop.fnt_serialize(metadata, metadata_size / 4);
                    if (ret == -1) {
                        return -1;
                    }
                }
            }
            return 0;
        }
    }

    /*
     * At this point we want a parity to be reconstructed
     * If systematic we may need to decode if a data is missing.
     * If non-systematic we always need to decode
     */
    std::vector<std::vector<uint8_t>> blocks(fec->n_data);
    if (fec->type == quadiron::fec::FecType::SYSTEMATIC) {
        for (unsigned i = 0; i < fec->n_data; i++) {
            if (missing_idxs[i]) {
                need_decode = 1;
                wanted_data_vec[i] = true;
                blocks.at(i).resize(block_size);
                data_vec[i] = blocks.at(i).data();
            }
        }
    } else {
        need_decode = 1;
        for (unsigned i = 0; i < fec->n_data; i++) {
            wanted_data_vec[i] = true;
            blocks.at(i).resize(block_size);
            data_vec[i] = blocks.at(i).data();
        }
    }

    if (need_decode) {
        res = fec->decode_blocks_vertical(
            data_vec,
            parities_vec,
            parities_props,
            missing_idxs_vec,
            wanted_data_vec,
            block_size);
        if (!res) {
            return -1;
        }
    }

    /*
     * At this point we have all data blocks.
     */
    if (fec->type == quadiron::fec::FecType::SYSTEMATIC) {
        wanted_idxs_vec[destination_idx - fec->n_data] = true;
    } else {
        wanted_idxs_vec[destination_idx] = true;
    }

    fec->encode_blocks_vertical(
        data_vec, parities_vec, parities_props, wanted_idxs_vec, block_size);

    if (fec->type == quadiron::fec::FecType::SYSTEMATIC) {
        for (unsigned i = 0; i < fec->n_parities; i++) {
            if (i == destination_idx - fec->n_data) {
                uint32_t* metadata = reinterpret_cast<uint32_t*>(parity[i]);
                int ret = parities_props[i].fnt_serialize(
                    metadata, metadata_size / 4);
                if (ret == -1) {
                    return -1;
                }
            }
        }
    } else {
        for (unsigned i = 0; i < fec->n_data; i++) {
            if (i == destination_idx) {
                uint32_t* metadata = reinterpret_cast<uint32_t*>(data[i]);
                int ret = parities_props[i].fnt_serialize(
                    metadata, metadata_size / 4);
                if (ret == -1) {
                    return -1;
                }
            }
        }
        for (unsigned i = 0; i < fec->n_parities; i++) {
            if (fec->n_data + i == destination_idx) {
                uint32_t* metadata = reinterpret_cast<uint32_t*>(parity[i]);
                int ret = parities_props[fec->n_data + i].fnt_serialize(
                    metadata, metadata_size / 4);
                if (ret == -1) {
                    return -1;
                }
            }
        }
    }

    return 0;
}

void quadiron_hex_dump(uint8_t* buf, size_t size)
{
    quadiron::hex_dump(std::cerr, buf, size, true);
}
}
