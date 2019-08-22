/*
 * Copyright 2017-2019 Scality
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
#ifndef __QUAD_BENCH_COMP_PERF_FEC_BASE_H__
#define __QUAD_BENCH_COMP_PERF_FEC_BASE_H__

#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include "arith.h"
#include "fec_base.h"
#include "gf_base.h"
#include "perf_base.h"
#include "simd/simd.h"

namespace arith = quadiron::arith;
namespace gf = quadiron::gf;
namespace fec = quadiron::fec;
namespace simd = quadiron::simd;

constexpr size_t KB = 1024;

// Pkt_size range: number of vectors (registers)
constexpr size_t MIN_VEC_LEN = 32;
constexpr size_t MAX_VEC_LEN = 64;
// Block size
constexpr size_t MIN_BLOCK_SIZE = 16 * KB;
constexpr size_t MAX_BLOCK_SIZE = 16 * KB;
// Fragment size
constexpr size_t FRAGMENT_SIZE_BYTES = 50 * KB;
// Code length
constexpr size_t MIN_CODE_LEN = 32;
constexpr size_t MAX_CODE_LEN = 1024;
// "Inversed" code rate
// k = {1, 2, .., MAX_INV_RATE-1} * n / MAX_INV_RATE
constexpr unsigned MAX_INV_RATE = 4;
// Fec type
constexpr unsigned FEC_NON_SYSTEMATIC = 1;
constexpr unsigned FEC_SYSTEMATIC = 2;

template <typename T>
class PerfFec : public PerfBase<T> {
  public:
    // to create files once
    static bool do_create_files;

    std::vector<std::string> d_files;
    std::vector<std::string> c_files;
    std::vector<std::string> r_files;

    size_t fragment_size = FRAGMENT_SIZE_BYTES;

    void create_files(
        std::vector<std::string>& files,
        size_t bytes,
        bool rand_fill = true)
    {
        std::vector<uint8_t> buf(bytes);
        const size_t files_nb = files.size();

        if (rand_fill) {
            for (size_t i = 0; i < files_nb; ++i) {
                this->randomize_byte_buffer(buf.data(), bytes);

                const std::string file_name = files[i];
                std::ofstream outfile(file_name, std::ofstream::binary);
                outfile.write(reinterpret_cast<char*>(buf.data()), bytes);
                outfile.close();
            }
        } else {
            for (size_t i = 0; i < files_nb; ++i) {
                const std::string file_name = files[i];
                std::ofstream outfile(file_name, std::ofstream::binary);
                outfile.seekp(bytes);
                outfile << 'a';
                outfile.close();
            }
        }
    }

    void create_files_name(
        const std::string& prefix,
        const std::string& suffix,
        size_t files_nb,
        std::vector<std::string>& files)
    {
        files.reserve(files_nb);
        for (size_t i = 0; i < files_nb; ++i) {
            std::string file_name =
                prefix + std::string("_") + std::to_string(i) + suffix;
            files.push_back(file_name);
        }
    }

    void prepare_codec_streams()
    {
        const size_t coded_files_nb = MAX_CODE_LEN;
        const size_t data_files_nb =
            coded_files_nb * (MAX_INV_RATE - 1) / MAX_INV_RATE;
        const size_t bytes = fragment_size;

        create_files_name("data", ".dat", data_files_nb, d_files);
        create_files_name("coded", ".dat", coded_files_nb, c_files);
        create_files_name("recovered", ".dat", data_files_nb, r_files);

        if (PerfFec<T>::do_create_files) {
            create_files(d_files, bytes);
            create_files(c_files, bytes, false /* not fill random value */);
            create_files(r_files, bytes, false /* not fill random value */);
            PerfFec<T>::do_create_files = false;
        }
    }

    void open_streams(
        const std::vector<std::string>& files,
        std::vector<std::istream*>& streams,
        bool in_and_out_flag = false)
    {
        std::ios_base::openmode mode = std::ios_base::in;
        if (in_and_out_flag) {
            mode |= std::ios_base::out;
        }
        for (size_t i = 0; i < streams.size(); ++i) {
            streams[i] = new std::ifstream(files[i], mode);
        }
    }

    void open_streams(
        const std::vector<std::string>& files,
        std::vector<std::ostream*>& streams,
        bool in_and_out_flag = false)
    {
        std::ios_base::openmode mode = std::ios_base::out;
        if (in_and_out_flag) {
            mode |= std::ios_base::in;
        }
        for (size_t i = 0; i < streams.size(); ++i) {
            streams[i] = new std::ofstream(files[i], mode);
        }
    }

    void close_streams(std::vector<std::istream*>& streams)
    {
        for (size_t i = 0; i < streams.size(); ++i) {
            (static_cast<std::ifstream*>(streams[i]))->close();
            delete streams[i];
        }
    }

    void close_streams(std::vector<std::ostream*>& streams)
    {
        for (size_t i = 0; i < streams.size(); ++i) {
            (static_cast<std::ofstream*>(streams[i]))->close();
            delete streams[i];
        }
    }

    void base_enc_dec(
        benchmark::State& st,
        fec::FecCode<T>& fec,
        size_t n_len,
        size_t n_data,
        size_t pkt_size)
    {
        const quadiron::gf::Field<T>& gf = fec.get_gf();
        const size_t n_outputs = fec.get_n_outputs();

        vec::Buffers<T> data_frags(n_data, pkt_size);
        vec::Buffers<T> encoded_frags(n_outputs, pkt_size);
        vec::Buffers<T> decoded_frags(n_data, pkt_size);

        this->randomize_data(data_frags);

        vec::Vector<T> fragments_ids(gf, n_data);
        std::vector<size_t> ids(n_len);
        for (size_t i = 0; i < n_len; ++i) {
            ids[i] = i;
        }
        std::vector<quadiron::Properties> props(n_outputs);
        for (size_t i = 0; i < n_outputs; ++i) {
            props[i] = quadiron::Properties();
        }

        const std::vector<T*>& data_mem = data_frags.get_mem();

        std::chrono::duration<double, std::micro> enc_elapsed_us{};
        std::chrono::duration<double, std::micro> dec_elapsed_us{};
        for (auto _ : st) {
            st.PauseTiming();
            for (size_t i = 0; i < n_outputs; ++i) {
                props[i].clear();
            }
            this->randomize_data(data_frags);
            st.ResumeTiming();

            auto enc_start = std::chrono::high_resolution_clock::now();
            fec.encode(encoded_frags, props, 0, data_frags);
            auto enc_end = std::chrono::high_resolution_clock::now();
            enc_elapsed_us += std::chrono::duration_cast<
                std::chrono::duration<double, std::micro>>(enc_end - enc_start);

            st.PauseTiming();
            std::random_shuffle(ids.begin(), ids.end());
            for (size_t i = 0; i < n_data; ++i) {
                fragments_ids.set(i, ids.at(i));
            }

            const std::vector<T*>& enc_mem = encoded_frags.get_mem();
            std::vector<T*> rec_mem(n_data);
            if (fec.type == fec::FecType::SYSTEMATIC) {
                for (size_t i = 0; i < n_data; ++i) {
                    const T id = fragments_ids.get(i);
                    if (id < n_data) {
                        rec_mem[i] = data_mem[id];
                    } else {
                        rec_mem[i] = enc_mem[id - n_data];
                    }
                }
            } else {
                for (size_t i = 0; i < n_data; ++i) {
                    rec_mem[i] = enc_mem[fragments_ids.get(i)];
                }
            }
            vec::Buffers<T> received_frags(n_data, pkt_size, rec_mem);

            std::unique_ptr<fec::DecodeContext<T>> context =
                fec.init_context_dec(
                    fragments_ids, props, pkt_size, &decoded_frags);

            st.ResumeTiming();

            auto dec_start = std::chrono::high_resolution_clock::now();
            fec.decode(*context, decoded_frags, props, 0, received_frags);
            auto dec_end = std::chrono::high_resolution_clock::now();
            dec_elapsed_us += std::chrono::duration_cast<
                std::chrono::duration<double, std::micro>>(dec_end - dec_start);
        }

        const size_t enc_bytes =
            st.iterations() * n_len * pkt_size * this->word_size;
        const size_t dec_bytes =
            st.iterations() * n_data * pkt_size * this->word_size;
        st.counters.insert(
            {{"|(1) Systematic", fec.type == fec::FecType::SYSTEMATIC},
             {"|(2) Code rate",
              static_cast<double>(n_data) / static_cast<double>(n_len)},
             {"|(3) Code len", n_len},
             {"|(4) Data len", n_data},
             {"|(5) Packet size",
              benchmark::Counter(
                  pkt_size,
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)},
             {"|(6) Encode time (us)",
              enc_elapsed_us.count() / st.iterations()},
             {"|(7) Encode speed (B/s)",
              benchmark::Counter(
                  1'000'000 * static_cast<double>(enc_bytes)
                      / enc_elapsed_us.count(),
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)},
             {"|(8) Decode time (us)",
              dec_elapsed_us.count() / st.iterations()},
             {"|(9) Decode speed (B/s)",
              benchmark::Counter(
                  1'000'000 * static_cast<double>(dec_bytes)
                      / dec_elapsed_us.count(),
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)}});
    }

    void base_enc_dec_blocks(
        benchmark::State& st,
        fec::FecCode<T>& fec,
        size_t n_len,
        size_t n_data,
        size_t block_size,
        size_t pkt_size)
    {
        const size_t n_outputs = fec.get_n_outputs();

        std::vector<uint8_t> k_data(n_data * block_size);
        std::vector<uint8_t> m_data(n_outputs * block_size);

        std::vector<uint8_t*> data_bufs(n_data);
        for (size_t i = 0; i < n_data; ++i) {
            data_bufs[i] = k_data.data() + i * block_size;
        }
        std::vector<uint8_t*> parities_bufs(n_outputs);
        for (size_t i = 0; i < n_outputs; ++i) {
            parities_bufs[i] = m_data.data() + i * block_size;
        }

        std::vector<quadiron::Properties> props(n_outputs);
        for (size_t i = 0; i < n_outputs; ++i) {
            props[i] = quadiron::Properties();
        }

        std::vector<int> missing_idx(n_len);

        std::vector<size_t> ids(n_len);
        for (size_t i = 0; i < n_len; ++i) {
            ids[i] = i;
        }

        std::vector<bool> wanted_idxs(n_len);
        std::fill_n(wanted_idxs.begin(), n_len, true);

        std::chrono::duration<double, std::micro> enc_elapsed_us{};
        std::chrono::duration<double, std::micro> dec_elapsed_us{};
        for (auto _ : st) {
            st.PauseTiming();
            for (size_t i = 0; i < n_outputs; ++i) {
                props[i].clear();
            }
            this->randomize_byte_buffer(k_data.data(), n_data * block_size);
            st.ResumeTiming();

            auto enc_start = std::chrono::high_resolution_clock::now();
            fec.encode_blocks_vertical(
                data_bufs, parities_bufs, props, wanted_idxs, block_size);
            auto enc_end = std::chrono::high_resolution_clock::now();
            enc_elapsed_us += std::chrono::duration_cast<
                std::chrono::duration<double, std::micro>>(enc_end - enc_start);

            st.PauseTiming();

            std::vector<uint8_t*> parities_for_dec(n_outputs);
            std::random_shuffle(ids.begin(), ids.end());

            if (fec.type == fec::FecType::SYSTEMATIC) {
                for (size_t i = 0; i < n_data; ++i) {
                    if (ids[i] >= n_data) {
                        const T id = ids[i] - n_data;
                        parities_for_dec[id] = parities_bufs[id];
                    }
                }
            } else {
                for (size_t i = 0; i < n_data; ++i) {
                    parities_for_dec[ids[i]] = parities_bufs[ids[i]];
                }
            }
            std::fill_n(missing_idx.begin(), n_len, 0);
            for (size_t i = n_data; i < n_len; ++i) {
                missing_idx[ids[i]] = 1;
            }

            st.ResumeTiming();

            auto dec_start = std::chrono::high_resolution_clock::now();
            fec.decode_blocks_vertical(
                data_bufs,
                parities_for_dec,
                props,
                missing_idx,
                wanted_idxs,
                block_size);
            auto dec_end = std::chrono::high_resolution_clock::now();
            dec_elapsed_us += std::chrono::duration_cast<
                std::chrono::duration<double, std::micro>>(dec_end - dec_start);
        }

        const size_t enc_bytes = st.iterations() * n_len * block_size;
        const size_t dec_bytes = st.iterations() * n_data * block_size;
        st.counters.insert(
            {{"|(1) Systematic", fec.type == fec::FecType::SYSTEMATIC},
             {"|(2) Code rate",
              static_cast<double>(n_data) / static_cast<double>(n_len)},
             {"|(3) Code len", n_len},
             {"|(4) Data len", n_data},
             {"|(5) Block size",
              benchmark::Counter(
                  block_size,
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)},
             {"|(6) Packet size",
              benchmark::Counter(
                  pkt_size,
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)},
             {"|(7) Encode time (us)",
              enc_elapsed_us.count() / st.iterations()},
             {"|(8) Encode speed (B/s)",
              benchmark::Counter(
                  1'000'000 * static_cast<double>(enc_bytes)
                      / enc_elapsed_us.count(),
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)},
             {"|(9) Decode time (us)",
              dec_elapsed_us.count() / st.iterations()},
             {"|(10) Decode speed (B/s)",
              benchmark::Counter(
                  1'000'000 * static_cast<double>(dec_bytes)
                      / dec_elapsed_us.count(),
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)}});
    }

    void base_enc_dec_streams(
        benchmark::State& st,
        fec::FecCode<T>& fec,
        size_t n_len,
        size_t n_data,
        size_t pkt_size)
    {
        const size_t n_outputs = fec.get_n_outputs();

        std::vector<std::istream*> data_streams(n_data);
        std::vector<std::ostream*> enc_streams(n_outputs);
        std::vector<std::ostream*> dec_streams(n_data);
        std::vector<std::istream*> all_rec_streams(n_outputs);

        std::vector<quadiron::Properties> props(n_outputs);
        for (size_t i = 0; i < n_outputs; ++i) {
            props[i] = quadiron::Properties();
        }

        std::vector<int> missing_idx(n_len);

        std::vector<size_t> ids(n_len);
        for (size_t i = 0; i < n_len; ++i) {
            ids[i] = i;
        }

        std::chrono::duration<double, std::micro> enc_elapsed_us{};
        std::chrono::duration<double, std::micro> dec_elapsed_us{};
        for (auto _ : st) {
            st.PauseTiming();

            for (size_t i = 0; i < n_outputs; ++i) {
                props[i].clear();
            }
            open_streams(d_files, data_streams);
            open_streams(c_files, enc_streams);

            st.ResumeTiming();

            auto enc_start = std::chrono::high_resolution_clock::now();
            fec.encode_streams_vertical(data_streams, enc_streams, props);
            auto enc_end = std::chrono::high_resolution_clock::now();
            enc_elapsed_us += std::chrono::duration_cast<
                std::chrono::duration<double, std::micro>>(enc_end - enc_start);

            st.PauseTiming();

            close_streams(data_streams);
            close_streams(enc_streams);

            open_streams(r_files, dec_streams);
            open_streams(c_files, all_rec_streams);

            std::random_shuffle(ids.begin(), ids.end());
            std::vector<std::istream*> rec_streams(n_len, nullptr);

            if (fec.type == fec::FecType::SYSTEMATIC) {
                open_streams(d_files, data_streams);
                for (size_t i = 0; i < n_data; ++i) {
                    const T id = ids[i];
                    if (id < n_data) {
                        rec_streams[id] = data_streams[id];
                    } else {
                        rec_streams[id] = all_rec_streams[id - n_data];
                    }
                }
            } else {
                for (size_t i = 0; i < n_data; ++i) {
                    rec_streams[ids[i]] = all_rec_streams[ids[i]];
                }
            }

            st.ResumeTiming();

            auto dec_start = std::chrono::high_resolution_clock::now();
            fec.decode_streams_vertical(
                data_streams, rec_streams, props, dec_streams);
            auto dec_end = std::chrono::high_resolution_clock::now();
            dec_elapsed_us += std::chrono::duration_cast<
                std::chrono::duration<double, std::micro>>(dec_end - dec_start);

            if (fec.type == fec::FecType::SYSTEMATIC) {
                close_streams(data_streams);
            }
            close_streams(dec_streams);
            close_streams(all_rec_streams);
        }

        const size_t enc_bytes = st.iterations() * n_len * fragment_size;
        const size_t dec_bytes = st.iterations() * n_data * fragment_size;
        st.counters.insert(
            {{"|(1) Systematic", fec.type == fec::FecType::SYSTEMATIC},
             {"|(2) Code rate",
              static_cast<double>(n_data) / static_cast<double>(n_len)},
             {"|(3) Code len", n_len},
             {"|(4) Data len", n_data},
             {"|(5) Fragment size",
              benchmark::Counter(
                  fragment_size,
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)},
             {"|(6) Packet size",
              benchmark::Counter(
                  pkt_size,
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)},
             {"|(7) Encode time (us)",
              enc_elapsed_us.count() / st.iterations()},
             {"|(8) Encode speed (B/s)",
              benchmark::Counter(
                  1'000'000 * static_cast<double>(enc_bytes)
                      / enc_elapsed_us.count(),
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)},
             {"|(9) Decode time (us)",
              dec_elapsed_us.count() / st.iterations()},
             {"|(10) Decode speed (B/s)",
              benchmark::Counter(
                  1'000'000 * static_cast<double>(dec_bytes)
                      / dec_elapsed_us.count(),
                  benchmark::Counter::kDefaults,
                  benchmark::Counter::OneK::kIs1024)}});
    }
};

// Initialize static member of class PerfFec
template <typename T>
bool PerfFec<T>::do_create_files = true;

#endif
