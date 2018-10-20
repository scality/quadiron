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
#include <string>

#include "benchmark.h"

template <typename T>
Benchmark<T>::Benchmark(Params_t* params)
{
    this->params = params;
    this->fec_type = params->fec_type;
    this->word_size = params->word_size;
    this->k = params->k;
    this->m = params->m;
    this->n = params->k + params->m;
    this->n_c = this->n; // if systematic, it will be updated in init()
    this->operation_on_packet = params->operation_on_packet;
    this->pkt_size = params->pkt_size;
    this->chunk_size = params->chunk_size;
    this->samples_nb = params->samples_nb;
    if (params->extra_param > -1) {
        this->extra_param = params->extra_param;
    }

    int error;
    if ((error = this->check_params()) < 0) {
        std::cerr << errors_desc.at(error) << std::endl;
        // NOLINTNEXTLINE(cert-err60-cpp)
        throw std::invalid_argument(errors_desc.at(error));
    }

    if ((error = this->init()) < 0) {
        std::cerr << errors_desc.at(error) << std::endl;
        // NOLINTNEXTLINE(cert-err60-cpp)
        throw std::invalid_argument(errors_desc.at(error));
    }
}

template <typename T>
Benchmark<T>::~Benchmark()
{
    if (fec != nullptr)
        delete fec;
    if (prng != nullptr)
        delete prng;
    if (c_chunks_id != nullptr)
        delete c_chunks_id;

    if (d_streams != nullptr) {
        for (int i = 0; i < k; i++) {
            delete d_streams->at(i);
        }
        delete d_streams;
    }
    if (c_streams != nullptr) {
        for (int i = 0; i < n_c; i++) {
            delete c_streams->at(i);
        }
        delete c_streams;
    }
    if (a_streams != nullptr) {
        for (int i = 0; i < n; i++) {
            delete a_streams->at(i);
        }
        delete a_streams;
    }
    if (r_streams != nullptr) {
        for (int i = 0; i < k; i++) {
            delete r_streams->at(i);
        }
        delete r_streams;
    }

    if (d_istreambufs != nullptr) {
        for (int i = 0; i < k; i++) {
            delete d_istreambufs->at(i);
        }
        delete d_istreambufs;
    }
    if (c_istreambufs != nullptr) {
        for (int i = 0; i < n_c; i++) {
            delete c_istreambufs->at(i);
        }
        delete c_istreambufs;
    }
    if (c_ostreambufs != nullptr) {
        for (int i = 0; i < n_c; i++) {
            delete c_ostreambufs->at(i);
        }
        delete c_ostreambufs;
    }
    if (r_ostreambufs != nullptr) {
        for (int i = 0; i < k; i++) {
            delete r_ostreambufs->at(i);
        }
        delete r_ostreambufs;
    }

    if (d_chunks != nullptr) {
        for (int i = 0; i < k; i++) {
            this->allocator.deallocate(d_chunks->at(i), chunk_size);
        }
        delete d_chunks;
    }
    if (c_chunks != nullptr) {
        for (int i = 0; i < n_c; i++) {
            this->allocator.deallocate(c_chunks->at(i), chunk_size);
        }
        delete c_chunks;
    }
    if (r_chunks != nullptr) {
        for (int i = 0; i < k; i++) {
            this->allocator.deallocate(r_chunks->at(i), chunk_size);
        }
        delete r_chunks;
    }

    if (enc_stats != nullptr)
        delete enc_stats;
    if (dec_stats != nullptr)
        delete dec_stats;
}

template <typename T>
int Benchmark<T>::init()
{
    switch (fec_type) {
    case EC_TYPE_RS_GF2N_V:
        fec = new quadiron::fec::RsGf2n<T>(
            word_size, k, m, quadiron::fec::RsMatrixType::VANDERMONDE);
        break;
    case EC_TYPE_RS_GF2N_C:
        fec = new quadiron::fec::RsGf2n<T>(
            word_size, k, m, quadiron::fec::RsMatrixType::CAUCHY);
        break;
    case EC_TYPE_RS_GF2N_FFT:
        fec = new quadiron::fec::RsGf2nFft<T>(word_size, k, m);
        break;
    case EC_TYPE_RS_GF2N_FFT_ADD:
        fec = new quadiron::fec::RsGf2nFftAdd<T>(word_size, k, m);
        break;
    case EC_TYPE_RS_GFP_FFT:
        fec = new quadiron::fec::RsGfpFft<T>(word_size, k, m);
        break;
    case EC_TYPE_RS_NF4:
        fec = new quadiron::fec::RsNf4<T>(word_size, k, m, pkt_size);
        break;
    case EC_TYPE_RS_FNT:
        fec = new quadiron::fec::RsFnt<T>(
            quadiron::fec::FecType::NON_SYSTEMATIC, word_size, k, m, pkt_size);
        break;
    case EC_TYPE_RS_FNT_SYS:
        fec = new quadiron::fec::RsFnt<T>(
            quadiron::fec::FecType::SYSTEMATIC, word_size, k, m, pkt_size);
        break;
    default:
        return ERR_FEC_TYPE_NOT_SUPPORTED;
    }

    this->systematic_ec = (fec->type == quadiron::fec::FecType::SYSTEMATIC);
    if (this->systematic_ec) {
        this->n_c = this->m;
    }

    this->prng = new PRNG(time(nullptr));

    // Allocate memory for data
    int i;
    d_chunks = new std::vector<uint8_t*>(k);
    c_chunks = new std::vector<uint8_t*>(n_c);
    r_chunks = new std::vector<uint8_t*>(k);

    for (i = 0; i < k; i++) {
        d_chunks->at(i) = this->allocator.allocate(chunk_size);
    }
    for (i = 0; i < n_c; i++) {
        c_chunks->at(i) = this->allocator.allocate(chunk_size);
    }
    for (i = 0; i < k; i++) {
        r_chunks->at(i) = this->allocator.allocate(chunk_size);
    }

    // Allocate memory for iostreambufs
    d_istreambufs = new std::vector<istreambuf<char>*>(k);
    c_istreambufs = new std::vector<istreambuf<char>*>(n_c);
    c_ostreambufs = new std::vector<ostreambuf<char>*>(n_c);
    r_ostreambufs = new std::vector<ostreambuf<char>*>(k);
    for (i = 0; i < k; i++) {
        d_istreambufs->at(i) = new istreambuf<char>(
            reinterpret_cast<char*>(d_chunks->at(i)), chunk_size);
    }
    for (i = 0; i < n_c; i++) {
        c_istreambufs->at(i) = new istreambuf<char>(
            reinterpret_cast<char*>(c_chunks->at(i)), chunk_size);
        c_ostreambufs->at(i) = new ostreambuf<char>(
            reinterpret_cast<char*>(c_chunks->at(i)), chunk_size);
    }
    for (i = 0; i < k; i++) {
        r_ostreambufs->at(i) = new ostreambuf<char>(
            reinterpret_cast<char*>(r_chunks->at(i)), chunk_size);
    }

    // Allocate memory for streams
    d_streams = new std::vector<std::istream*>(k);
    c_streams = new std::vector<std::ostream*>(n_c);
    a_streams = new std::vector<std::istream*>(n);
    r_streams = new std::vector<std::ostream*>(k);
    c_props = std::vector<quadiron::Properties>(n_c);

    for (i = 0; i < k; i++) {
        d_streams->at(i) = new std::istream(d_istreambufs->at(i));
    }

    for (i = 0; i < n_c; i++) {
        c_streams->at(i) = new std::ostream(c_ostreambufs->at(i));
    }

    if (systematic_ec) {
        for (i = 0; i < k; i++) {
            a_streams->at(i) = new std::istream(d_istreambufs->at(i));
        }
        for (; i < n; i++) {
            a_streams->at(i) = new std::istream(c_istreambufs->at(i - k));
        }
    } else {
        for (i = 0; i < n; i++) {
            a_streams->at(i) = new std::istream(c_istreambufs->at(i));
        }
    }

    for (i = 0; i < k; i++) {
        r_streams->at(i) = new std::ostream(r_ostreambufs->at(i));
    }

    // init index of chunks
    c_chunks_id = new std::vector<int>();
    for (i = 0; i < n; i++) {
        c_chunks_id->push_back(i);
    }

    this->enc_stats = new Stats_t("Encode", chunk_size * n_c);
    this->dec_stats = new Stats_t("Decode", chunk_size * k);
    return 1;
}

template <typename T>
int Benchmark<T>::check_params()
{
    if (word_size <= 0) {
        return ERR_WORD_SIZE;
    }

    if (fec_type == EC_TYPE_RS_NF4) {
        if (word_size < 2) {
            return ERR_WORD_SIZE;
        }
    }

    if (sizeof(T) < word_size) {
        return ERR_COMPT_WORD_SIZE_T;
    }
    if (fec_type == EC_TYPE_RS_FNT || fec_type == EC_TYPE_RS_FNT_SYS
        || fec_type == EC_TYPE_RS_GFP_FFT) {
        if (sizeof(T) <= word_size) {
            return ERR_COMPT_WORD_SIZE_T;
        }
    }
    if (fec_type == EC_TYPE_RS_FNT || fec_type == EC_TYPE_RS_FNT_SYS) {
        if (word_size > 2)
            return ERR_WORD_SIZE;
    }
    if (fec_type == EC_TYPE_RS_NF4) {
        if (sizeof(T) < 2 * word_size) {
            return ERR_COMPT_WORD_SIZE_T;
        }
    }

    size_t wordsize_limit = quadiron::arith::log2<T>(n) + 1;
    if (wordsize_limit > 8 * word_size) {
        return ERR_COMPT_CODE_LEN_T;
    }

    // adjust chunk_size
    if (pkt_size > 0 && chunk_size % (pkt_size * word_size) > 0)
        chunk_size = ((chunk_size - 1) / (pkt_size * word_size) + 1)
                     * (pkt_size * word_size);

    return 1;
}

template <typename T>
void Benchmark<T>::gen_data()
{
    for (int i = 0; i < k; i++) {
        prng->gen_chunk(d_chunks->at(i), chunk_size);
    }
    if (!check(d_chunks)) {
        std::cerr << errors_desc.at(ERR_FAILED_CHUNK) << std::endl;
        // NOLINTNEXTLINE(cert-err60-cpp)
        throw std::runtime_error(errors_desc.at(ERR_FAILED_CHUNK));
    }
}

template <typename T>
bool Benchmark<T>::check(std::vector<uint8_t*>* chunks)
{
    for (std::vector<int>::size_type i = 0; i < chunks->size(); i++) {
        if (!prng->check_chunk(chunks->at(i), chunk_size)) {
            // dump_chunk("failed check chunk", chunks->at(i));
            return false;
        }
    }
    return true;
}

template <typename T>
bool Benchmark<T>::compare(
    std::vector<uint8_t*>* arr1,
    std::vector<uint8_t*>* arr2)
{
    if (arr1->size() != arr2->size()) {
        std::cout << "Sizes are different\n";
        return false;
    }
    for (std::vector<int>::size_type i = 0; i < arr1->size(); i++) {
        if (prng->get_crc(arr1->at(i), chunk_size)
            != prng->get_crc(arr2->at(i), chunk_size)) {
            // dump_chunk("chunk1", arr1->at(i));
            // dump_chunk("chunk2", arr2->at(i));
            std::cout << "CRCs are different\n";
            return false;
        }
    }
    if (!check(arr1) || !check(arr2)) {
        std::cout << "Contents are different\n";
        return false;
    }
    return true;
}

template <typename T>
void Benchmark<T>::dump(const char* name, std::vector<uint8_t*>* chunks)
{
    std::cout << name << ": ";
    for (std::vector<int>::size_type i = 0; i < chunks->size(); i++) {
        uint8_t* chunk = chunks->at(i);
        for (int j = 0; j < chunk_size; j++) {
            std::cout << unsigned(chunk[j]) << " ";
        }
        std::cout << "|";
    }
    std::cout << std::endl;
}

template <typename T>
void Benchmark<T>::dump_chunk(const char* name, uint8_t* chunk)
{
    std::cout << name << ": ";
    for (int j = 0; j < chunk_size; j++) {
        std::cout << unsigned(chunk[j]) << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void Benchmark<T>::reset_d_streams()
{
    for (int i = 0; i < k; i++) {
        d_streams->at(i)->clear();
        d_streams->at(i)->rdbuf()->pubseekpos(0);
    }
}

template <typename T>
void Benchmark<T>::reset_c_streams()
{
    for (int i = 0; i < n_c; i++) {
        c_streams->at(i)->clear();
        c_streams->at(i)->rdbuf()->pubseekpos(0);
        c_props.at(i).clear();
    }
}

template <typename T>
void Benchmark<T>::reset_a_streams()
{
    int i;
    for (i = 0; i < n; i++) {
        a_streams->at(i)->clear();
        a_streams->at(i)->rdbuf()->pubseekpos(0);
    }
}

template <typename T>
void Benchmark<T>::reset_r_streams()
{
    for (int i = 0; i < k; i++) {
        r_streams->at(i)->clear();
        r_streams->at(i)->rdbuf()->pubseekpos(0);
    }
}

template <typename T>
void Benchmark<T>::get_avail_chunks(
    std::vector<std::istream*>* avail_d_chunks,
    std::vector<std::istream*>* avail_c_chunks,
    std::vector<quadiron::Properties>& avail_c_props)
{
    std::random_shuffle(c_chunks_id->begin(), c_chunks_id->end());

    int i;
    for (i = 0; i < k; i++) {
        avail_d_chunks->at(i) = nullptr;
    }
    for (i = 0; i < n_c; i++) {
        avail_c_chunks->at(i) = nullptr;
        avail_c_props.at(i) = quadiron::Properties();
    }
    if (systematic_ec) {
        int avail_d_chunks_nb = 0;
        for (i = 0; i < k; i++) {
            int j = c_chunks_id->at(i);
            if (j < k) {
                avail_d_chunks->at(j) = a_streams->at(j);
                avail_d_chunks_nb++;
            } else {
                avail_c_chunks->at(j - k) = a_streams->at(j);
                avail_c_props.at(j - k) = c_props.at(j - k);
            }
        }
        // shuffle again if all data are available
        if (avail_d_chunks_nb == k)
            get_avail_chunks(avail_d_chunks, avail_c_chunks, avail_c_props);
    } else {
        for (i = 0; i < k; i++) {
            int j = c_chunks_id->at(i);
            avail_c_chunks->at(j) = a_streams->at(j);
            avail_c_props.at(j) = c_props.at(j);
        }
    }
}

template <typename T>
bool Benchmark<T>::encode()
{
    // this operation is done per trail
    reset_d_streams();
    reset_c_streams();

    if (operation_on_packet)
        fec->encode_packet(*d_streams, *c_streams, c_props);
    else
        fec->encode_bufs(*d_streams, *c_streams, c_props);

    // update stats
    enc_stats->add(fec->total_enc_usec);

    // dump("d_chunks", d_chunks);
    // dump("c_chunks", c_chunks);

    return true;
}

template <typename T>
bool Benchmark<T>::decode()
{
    std::vector<std::istream*> d_streams_shuffled(k, nullptr);
    std::vector<std::istream*> c_streams_shuffled(n_c, nullptr);
    std::vector<quadiron::Properties> c_props_shuffled(n_c);

    get_avail_chunks(
        &d_streams_shuffled, &c_streams_shuffled, c_props_shuffled);

    // this operation is done per trail
    reset_a_streams();
    reset_r_streams();

    if (operation_on_packet) {
        if (!fec->decode_packet(
                d_streams_shuffled,
                c_streams_shuffled,
                c_props_shuffled,
                *r_streams))
            return false;
    } else {
        if (!fec->decode_bufs(
                d_streams_shuffled,
                c_streams_shuffled,
                c_props_shuffled,
                *r_streams))
            return false;
    }

    if (!compare(d_chunks, r_chunks)) {
        std::cerr << errors_desc.at(ERR_FAILED_REPAIR_CHUNK) << std::endl;
        return false;
    }

    // dump("r_chunks", r_chunks);

    if (!check(r_chunks)) {
        std::cerr << errors_desc.at(ERR_FAILED_CHUNK) << std::endl;
        return false;
    }

    // update stats
    dec_stats->add(fec->total_dec_usec);

    return true;
}

template <typename T>
void Benchmark<T>::show(Stats_t* stats)
{
    if (params->compact_print == 2) {
        params->show_ec_desc();
        stats->show();
    } else {
        params->show_params();
        params->show_speed(
            stats->get_avg(), stats->get_std_dev(), stats->get_thrpt());
        params->show_end();
    }
}

template <typename T>
void Benchmark<T>::show(Stats_t* stats1, Stats_t* stats2)
{
    if (params->compact_print == 2) {
        params->show_ec_desc();
        stats1->show();
        stats2->show();
    } else {
        params->show_params();
        params->show_speed(
            stats1->get_avg(), stats1->get_std_dev(), stats1->get_thrpt());
        params->show_speed(
            stats2->get_avg(), stats2->get_std_dev(), stats2->get_thrpt());
        params->show_end();
    }
}

template <typename T>
bool Benchmark<T>::enc_only()
{
    enc_stats->begin();
    // this operation is done once per benchmark
    gen_data();

    for (uint32_t i = 0; i < samples_nb; i++) {
        if (!encode())
            return false;
    }

    enc_stats->end();
    show(enc_stats);

    return true;
}

template <typename T>
bool Benchmark<T>::dec_only()
{
    enc_stats->begin();
    dec_stats->begin();

    // this operation is done once per benchmark
    gen_data();

    if (!encode())
        return false;

    for (uint32_t i = 0; i < samples_nb; i++) {
        if (!decode())
            return false;
    }

    enc_stats->end();
    dec_stats->end();
    show(enc_stats, dec_stats);

    return true;
}

template <typename T>
bool Benchmark<T>::enc_dec()
{
    enc_stats->begin();
    dec_stats->begin();

    // this operation is done once per benchmark
    gen_data();

    for (uint32_t i = 0; i < samples_nb; i++) {
        if (!encode())
            return false;
        if (!decode())
            return false;
    }

    enc_stats->end();
    dec_stats->end();
    show(enc_stats, dec_stats);
    return true;
}

[[noreturn]] void xusage()
{
    std::cerr << "Usage: benchmark [options]\n"
              << "Options:\n"
              << "\t-e \tType of Reed-Solomon codes, either\n"
              << "\t\t\trs-gf2n-v: " << ec_desc.at(EC_TYPE_RS_GF2N_V) << '\n'
              << "\t\t\trs-gf2n-c: " << ec_desc.at(EC_TYPE_RS_GF2N_C) << '\n'
              << "\t\t\trs-gf2n-fft: " << ec_desc.at(EC_TYPE_RS_GF2N_FFT)
              << '\n'
              << "\t\t\trs-gf2n-fft-add: "
              << ec_desc.at(EC_TYPE_RS_GF2N_FFT_ADD) << '\n'
              << "\t\t\trs-gfp-fft: " << ec_desc.at(EC_TYPE_RS_GFP_FFT) << '\n'
              << "\t\t\trs-fnt: " << ec_desc.at(EC_TYPE_RS_FNT) << '\n'
              << "\t\t\trs-fnt-sys: " << ec_desc.at(EC_TYPE_RS_FNT_SYS) << '\n'
              << "\t\t\trs-nf4: " << ec_desc.at(EC_TYPE_RS_NF4) << '\n'
              << "\t\t\tall: All available Reed-solomon codes\n"
              << "\t-s \tScenario for benchmark, either\n"
              << "\t\t\tenc_only: Only encodings\n"
              << "\t\t\tdec_only: Only decodings\n"
              << "\t\t\tenc_dec: Encodings and decodings\n"
              << "\t-w \tWord size (bytes)\n"
              << "\t-k \tNumber of data chunks\n"
              << "\t-m \tNumber of parity chunks\n"
              << "\t-c \tChunk size (bytes)\n"
              << "\t-n \tNumber of samples per operation\n"
              << "\t-t \tSize of used integer type, either "
              << "2, 4, 8, 16 for uint16_t, uint32_t, uint64_t, __uint128_t\n"
              << "\t-g \tNumber of threads\n"
              << "\t-x \tExtra parameter\n\n";
    std::exit(EXIT_FAILURE);
}

template <typename T>
void run(Benchmark<T>* bench, Params_t* params)
{
    switch (params->sce_type) {
    case ENC_ONLY:
        bench->enc_only();
        break;
    case DEC_ONLY:
        bench->dec_only();
        break;
    case ENC_DEC:
        bench->enc_dec();
        break;
    }
}

template <typename T>
void init_run_bench(Params_t* params)
{
    try {
        Benchmark<T> bench(params);
        run<T>(&bench, params);
    } catch (const std::exception& e) {
        return;
    }
}

void run_scenario(Params_t* params)
{
    // get sizeof_T if necessary
    params->get_sizeof_T();

    switch (params->sizeof_T) {
    case 2:
        init_run_bench<uint16_t>(params);
        break;
    case 4:
        init_run_bench<uint32_t>(params);
        break;
    case 8:
        init_run_bench<uint64_t>(params);
        break;
    case 16:
        init_run_bench<__uint128_t>(params);
        break;
    default:
        std::cerr << errors_desc.at(ERR_T_NOT_SUPPORTED)
                  << " T: " << params->sizeof_T << std::endl;
        exit(0);
    }
}

void run_benchmark(Params_t* params)
{
    if (params->fec_type == EC_TYPE_ALL) {
        for (int type = EC_TYPE_ALL + 1; type < EC_TYPE_END; type++) {
            params->fec_type = static_cast<ec_type>(type);
            run_scenario(params);
        }
    } else {
        run_scenario(params);
    }
}

void do_join(std::thread& t)
{
    t.join();
}

int main(int argc, char** argv)
{
    PRNG prng;
    Params_t* params;
    int opt;

    params = new Params_t();
    while ((opt = getopt(argc, argv, "t:e:w:k:m:c:n:s:x:g:p:f:")) != -1) {
        switch (opt) {
        case 't':
            params->sizeof_T = std::stoi(optarg);
            break;
        case 'e':
            if (fec_type_map.find(optarg) == fec_type_map.end()) {
                std::cerr << errors_desc.at(ERR_FEC_TYPE_NOT_SUPPORTED)
                          << std::endl;
                xusage();
            }
            params->fec_type = fec_type_map.at(optarg);
            break;
        case 's':
            if (sce_type_map.find(optarg) == sce_type_map.end()) {
                std::cerr << errors_desc.at(ERR_SCENARIO_TYPE_NOT_SUPPORTED)
                          << std::endl;
                xusage();
            }
            params->sce_type = sce_type_map.at(optarg);
            break;
        case 'w':
            params->word_size = std::stoi(optarg);
            break;
        case 'k':
            params->k = std::stoi(optarg);
            break;
        case 'm':
            params->m = std::stoi(optarg);
            break;
        case 'c':
            params->chunk_size = std::stoi(optarg);
            break;
        case 'p':
            params->pkt_size = std::stoi(optarg);
            break;
        case 'n':
            params->samples_nb = std::stoi(optarg);
            break;
        case 'x':
            params->extra_param = std::stoi(optarg);
            break;
        case 'g':
            params->threads_nb = std::stoi(optarg);
            break;
        case 'f':
            params->compact_print = std::stoi(optarg);
            break;
        default:
            xusage();
        }
    }

    // Currently support operating on packet:RS_FNT
    if (params->fec_type != EC_TYPE_RS_FNT
        && params->fec_type != EC_TYPE_RS_FNT_SYS
        && params->fec_type != EC_TYPE_RS_NF4) {
        params->operation_on_packet = false;
    }

    if (params->pkt_size <= 0) {
        params->operation_on_packet = false;
    }
    if (params->compact_print == 2)
        params->print();
    else if (params->compact_print == 1)
        params->show_header();

    std::vector<std::thread> threads;
    for (uint32_t thr = 0; thr < params->threads_nb; thr++) {
        threads.push_back(std::thread(run_benchmark, params));
    }
    std::for_each(threads.begin(), threads.end(), do_join);

    delete params;
    return 0;
}
