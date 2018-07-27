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
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "quadiron.h"

int vflag = 0;
int tflag = 0;
char* prefix = nullptr;

void xusage()
{
    std::cerr << std::string("Usage: ") +
    "ec [-e rs-gf2n-v|rs-gf2n-c|rs-gf2n-fft|rs-gf2n-fft-add|rs-gfp-fft|rs-fnt|rs-nf4]" +
    "[-w word_size][-n n_data][-m n_parities][-p prefix][-v (verbose)]" +
    " -c (encode) | -r (repair)\n";
    exit(1);
}

void xperror(const char* str)
{
    std::cerr << str << "\n";
    exit(1);
}

char* xstrdup(const char* str)
{
    char* n;

    if (nullptr == (n = strdup(str)))
        xperror("malloc");
    return n;
}

/**
 * (re-)create missing prefix.c1 ... cm files
 *
 */
template <typename T>
void create_coding_files(
    quadiron::fec::FecCode<T>* fec,
    bool operation_on_packet = false)
{
    char filename[1024];
    std::vector<std::istream*> d_files(fec->n_data, nullptr);
    std::vector<std::ostream*> c_files(fec->n_outputs, nullptr);
    std::vector<std::ostream*> c_props_files(fec->n_outputs, nullptr);
    std::vector<quadiron::Properties> c_props(fec->n_outputs);

    for (unsigned i = 0; i < fec->n_data; i++) {
        snprintf(filename, sizeof(filename), "%s.d%d", prefix, i);
        if (vflag)
            std::cerr << "create: opening data " << filename << "\n";
        d_files[i] = new std::ifstream(filename);
    }

    for (unsigned i = 0; i < fec->n_outputs; i++) {
        snprintf(filename, sizeof(filename), "%s.c%d", prefix, i);
        if (vflag)
            std::cerr << "create: opening coding for writing " << filename
                      << "\n";
        c_files[i] = new std::ofstream(filename);
        snprintf(filename, sizeof(filename), "%s.c%d.props", prefix, i);
        if (vflag)
            std::cerr << "create: opening coding props for writing " << filename
                      << "\n";
        c_props_files[i] = new std::ofstream(filename);
    }

    if (operation_on_packet)
        fec->encode_packet(d_files, c_files, c_props);
    else
        fec->encode_bufs(d_files, c_files, c_props);

    for (unsigned i = 0; i < fec->n_data; i++) {
        (static_cast<std::ifstream*>(d_files[i]))->close();
        delete d_files[i];
    }

    for (unsigned i = 0; i < fec->n_outputs; i++) {
        *(c_props_files[i]) << c_props[i];

        (static_cast<std::ofstream*>(c_props_files[i]))->close();
        delete c_props_files[i];

        (static_cast<std::ofstream*>(c_files[i]))->close();
        delete c_files[i];
    }
}

/**
 * repair data files
 *
 */
template <typename T>
bool repair_data_files(quadiron::fec::FecCode<T>* fec)
{
    char filename[1024];
    std::vector<std::istream*> d_files(fec->n_data, nullptr);
    std::vector<std::istream*> c_files(fec->n_outputs, nullptr);
    std::vector<std::istream*> c_props_files(fec->n_outputs, nullptr);
    std::vector<quadiron::Properties> c_props(fec->n_outputs);
    std::vector<std::ostream*> r_files(fec->n_data, nullptr);

    // re-read data
    for (unsigned i = 0; i < fec->n_data; i++) {
        snprintf(filename, sizeof(filename), "%s.d%d", prefix, i);
        if (vflag)
            std::cerr << "repair: checking data " << filename << "\n";
        if (-1 == access(filename, F_OK)) {
            if (vflag)
                std::cerr << filename << " is missing\n";
            d_files[i] = nullptr;
            r_files[i] = new std::ofstream(filename);
        } else {
            r_files[i] = nullptr;
            d_files[i] = new std::ifstream(filename);
        }
    }

    for (unsigned i = 0; i < fec->n_outputs; i++) {
        snprintf(filename, sizeof(filename), "%s.c%d", prefix, i);
        if (vflag)
            std::cerr << "repair: checking coding " << filename << "\n";
        if (access(filename, F_OK)) {
            if (vflag)
                std::cerr << filename << " is missing\n";
            c_files[i] = nullptr;
        } else {
            c_files[i] = new std::ifstream(filename);
        }

        snprintf(filename, sizeof(filename), "%s.c%d.props", prefix, i);
        if (vflag)
            std::cerr << "repair: checking coding props " << filename << "\n";
        if (access(filename, F_OK)) {
            c_props_files[i] = nullptr;
        } else {
            c_props_files[i] = new std::ifstream(filename);
            *(c_props_files[i]) >> c_props[i];
        }
    }

    fec->decode_bufs(d_files, c_files, c_props, r_files);

    for (unsigned i = 0; i < fec->n_data; i++) {
        if (nullptr != d_files[i]) {
            (static_cast<std::ifstream*>(d_files[i]))->close();
            delete d_files[i];
        }
    }

    for (unsigned i = 0; i < fec->n_outputs; i++) {
        if (nullptr != c_props_files[i]) {
            (static_cast<std::ifstream*>(c_props_files[i]))->close();
            delete c_props_files[i];
        }

        if (nullptr != c_files[i]) {
            (static_cast<std::ifstream*>(c_files[i]))->close();
            delete c_files[i];
        }
    }

    for (unsigned i = 0; i < fec->n_data; i++) {
        if (nullptr != r_files[i]) {
            (static_cast<std::ofstream*>(r_files[i]))->close();
            delete r_files[i];
        }
    }

    return 0;
}

template <typename T>
void print_stats(quadiron::fec::FecCode<T>* fec)
{
    std::cerr << "enc,"
              << (fec->n_encode_ops != 0
                      ? fec->total_encode_cycles / fec->n_encode_ops
                      : 0)
              << ",";
    std::cerr << "dec,"
              << (fec->n_decode_ops != 0
                      ? fec->total_decode_cycles / fec->n_decode_ops
                      : 0)
              << ",";
}

template <typename T>
void print_fec_type(quadiron::fec::FecCode<T>* fec)
{
    switch (fec->type) {
    case quadiron::fec::FecType::SYSTEMATIC:
        std::cout << "SYSTEMATIC\n";
        break;
    case quadiron::fec::FecType::NON_SYSTEMATIC:
        std::cout << "NON_SYSTEMATIC\n";
        break;
    default:
        std::cout << "unknown\n";
        break;
    }
}

template <typename T>
void run_fec_rs_gf2n(
    int word_size,
    int n_data,
    int n_parities,
    quadiron::fec::RsMatrixType mflag,
    int rflag)
{
    quadiron::fec::RsGf2n<T>* fec;
    typename quadiron::fec::RsMatrixType gf2nrs_type;
    ;
    if (mflag == quadiron::fec::RsMatrixType::VANDERMONDE) {
        gf2nrs_type = quadiron::fec::RsMatrixType::VANDERMONDE;
    } else {
        gf2nrs_type = quadiron::fec::RsMatrixType::CAUCHY;
    }
    fec = new quadiron::fec::RsGf2n<T>(word_size, n_data, n_parities, gf2nrs_type);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec);
    print_stats<T>(fec);
    delete fec;
}

template <typename T>
void run_fec_rs_gf2n_fft(int word_size, int n_data, int n_parities, int rflag)
{
    quadiron::fec::RsGf2nFft<T>* fec;
    fec = new quadiron::fec::RsGf2nFft<T>(word_size, n_data, n_parities);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec);
    print_stats<T>(fec);
    delete fec;
}

template <typename T>
void run_fec_rs_gf2n_fft_add(
    int word_size,
    int n_data,
    int n_parities,
    int rflag)
{
    quadiron::fec::RsGf2nFftAdd<T>* fec;
    fec = new quadiron::fec::RsGf2nFftAdd<T>(word_size, n_data, n_parities);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec);
    print_stats<T>(fec);
    delete fec;
}

template <typename T>
void run_fec_rs_gfp_fft(
    unsigned word_size,
    int n_data,
    int n_parities,
    int rflag)
{
    assert(sizeof(T) > word_size);

    quadiron::fec::RsGfpFft<T>* fec;
    fec = new quadiron::fec::RsGfpFft<T>(word_size, n_data, n_parities);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec);
    print_stats<T>(fec);
    delete fec;
}

template <typename T>
void run_fec_rs_fnt(int word_size, int n_data, int n_parities, int rflag)
{
    quadiron::fec::RsFnt<T>* fec;
    size_t pkt_size = 1024;
    fec = new quadiron::fec::RsFnt<T>(word_size, n_data, n_parities, pkt_size);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec, true);
    print_stats<T>(fec);
    delete fec;
}

template <typename T>
void run_fec_rs_nf4(int word_size, int n_data, int n_parities, int rflag)
{
    quadiron::fec::RsNf4<T>* fec;
    fec = new quadiron::fec::RsNf4<T>(word_size, n_data, n_parities);

    if (tflag) {
        print_fec_type<T>(fec);
        exit(1);
    }
    if (rflag) {
        if (0 != repair_data_files<T>(fec)) {
            exit(1);
        }
    }
    create_coding_files<T>(fec);
    print_stats<T>(fec);
    delete fec;
}

enum ec_type {
    EC_TYPE_UNDEF = 0,
    EC_TYPE_RS_GF2N,
    EC_TYPE_RS_GF2N_FFT,
    EC_TYPE_RS_GF2N_FFT_ADD,
    EC_TYPE_RS_GFP_FFT,
    EC_TYPE_RS_FNT,
    EC_TYPE_RS_NF4,
};

bool check(int n, int word_size, ec_type eflag)
{
    // we suppose that code length is not too long, i.e. > 2^32
    if (word_size >= 4)
        return true;
    if (eflag == EC_TYPE_RS_FNT) {
        return (n <= (1LL << (8 * word_size)) + 1);
    } else {
        return (n <= (1LL << (8 * word_size)));
    }
}

int main(int argc, char** argv)
{
    int n_data, n_parities, opt;
    int cflag = 0;
    int rflag = 0;
    int uflag = 0;
    ec_type eflag = EC_TYPE_UNDEF;
    quadiron::fec::RsMatrixType mflag = quadiron::fec::RsMatrixType::VANDERMONDE;
    unsigned word_size = 0;

    n_data = n_parities = -1;
    while ((opt = getopt(argc, argv, "n:m:p:cruve:w:t")) != -1) {
        switch (opt) {
        case 'e':
            if (!strcmp(optarg, "rs-gf2n-v")) {
                eflag = EC_TYPE_RS_GF2N;
                mflag = quadiron::fec::RsMatrixType::VANDERMONDE;
            } else if (!strcmp(optarg, "rs-gf2n-c")) {
                eflag = EC_TYPE_RS_GF2N;
                mflag = quadiron::fec::RsMatrixType::CAUCHY;
            } else if (!strcmp(optarg, "rs-gf2n-fft")) {
                eflag = EC_TYPE_RS_GF2N_FFT;
            } else if (!strcmp(optarg, "rs-gf2n-fft-add")) {
                eflag = EC_TYPE_RS_GF2N_FFT_ADD;
            } else if (!strcmp(optarg, "rs-gfp-fft")) {
                eflag = EC_TYPE_RS_GFP_FFT;
            } else if (!strcmp(optarg, "rs-nf4")) {
                eflag = EC_TYPE_RS_NF4;
            } else if (!strcmp(optarg, "rs-fnt")) {
                eflag = EC_TYPE_RS_FNT;
            } else {
                xusage();
            }
            break;
        case 'w':
            word_size = std::stoi(optarg);
            break;
        case 'v':
            vflag = 1;
            break;
        case 'u':
            uflag = 1;
            break;
        case 'c':
            cflag = 1;
            break;
        case 'r':
            rflag = 1;
            break;
        case 'n':
            n_data = std::stoi(optarg);
            break;
        case 'm':
            n_parities = std::stoi(optarg);
            break;
        case 'p':
            prefix = xstrdup(optarg);
            break;
        case 't':
            tflag = 1;
            break;
        default: /* '?' */
            xusage();
        }
    }

    if (!eflag)
        xusage();

    if (!(uflag || cflag || rflag))
        xusage();

    if (0 == check(n_data + n_parities, word_size, eflag)) {
        std::cerr
            << "Number of fragments is too big compared to Galois field size\n";
        exit(1);
    }

    if (-1 == n_data || -1 == n_parities || nullptr == prefix)
        xusage();

    if (eflag == EC_TYPE_RS_FNT) {
        if (word_size <= 4) {
            run_fec_rs_fnt<uint32_t>(word_size, n_data, n_parities, rflag);
        } else if (word_size <= 8) {
            run_fec_rs_fnt<uint64_t>(word_size, n_data, n_parities, rflag);
        }
        // else if (word_size <= 16)
        // run_fec_rs_fnt<__uint128_t>(word_size, n_data, n_parities, rflag);
    } else if (eflag == EC_TYPE_RS_NF4) {
        if (word_size <= 2) {
            run_fec_rs_nf4<uint32_t>(word_size, n_data, n_parities, rflag);
        } else if (word_size <= 4) {
            run_fec_rs_nf4<uint64_t>(word_size, n_data, n_parities, rflag);
        } else if (word_size <= 8) {
            run_fec_rs_nf4<__uint128_t>(word_size, n_data, n_parities, rflag);
        }
    } else if (eflag == EC_TYPE_RS_GF2N) {
        if (word_size <= 4) {
            run_fec_rs_gf2n<uint32_t>(
                word_size, n_data, n_parities, mflag, rflag);
        } else if (word_size <= 8) {
            run_fec_rs_gf2n<uint64_t>(
                word_size, n_data, n_parities, mflag, rflag);
        } else if (word_size <= 16) {
            run_fec_rs_gf2n<__uint128_t>(
                word_size, n_data, n_parities, mflag, rflag);
        }
    } else if (eflag == EC_TYPE_RS_GFP_FFT) {
        if (word_size <= 7) {
            run_fec_rs_gfp_fft<uint64_t>(word_size, n_data, n_parities, rflag);
        } else if (word_size <= 15) {
            run_fec_rs_gfp_fft<__uint128_t>(
                word_size, n_data, n_parities, rflag);
        }
    } else if (eflag == EC_TYPE_RS_GF2N_FFT) {
        if (word_size <= 4) {
            run_fec_rs_gf2n_fft<uint32_t>(word_size, n_data, n_parities, rflag);
        } else if (word_size <= 8) {
            run_fec_rs_gf2n_fft<uint64_t>(word_size, n_data, n_parities, rflag);
        } else if (word_size <= 16) {
            run_fec_rs_gf2n_fft<__uint128_t>(
                word_size, n_data, n_parities, rflag);
        }
    } else if (eflag == EC_TYPE_RS_GF2N_FFT_ADD) {
        if (word_size <= 4) {
            run_fec_rs_gf2n_fft_add<uint32_t>(
                word_size, n_data, n_parities, rflag);
        } else if (word_size <= 8) {
            run_fec_rs_gf2n_fft_add<uint64_t>(
                word_size, n_data, n_parities, rflag);
        } else if (word_size <= 16) {
            run_fec_rs_gf2n_fft_add<__uint128_t>(
                word_size, n_data, n_parities, rflag);
        }
    }

    free(prefix);
    return 0;
}
