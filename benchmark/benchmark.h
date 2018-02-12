/* -*- mode: c++ -*- */
#ifndef __NTTEC_BENCH_BENCHMARK_H__
#define __NTTEC_BENCH_BENCHMARK_H__

#include <iomanip>
#include <thread>
#include <sys/stat.h>
#include <sys/types.h>

#include "nttec.h"

#include "iostreambuf.h"
#include "prng.h"
#include "stats.h"

enum ec_type {
    EC_TYPE_ALL = 0,
    EC_TYPE_FNTRS,
    EC_TYPE_NGFF4RS,
    EC_TYPE_GFPFFTRS,
    EC_TYPE_GF2NFFTADDRS,
    EC_TYPE_GF2NRSV,
    EC_TYPE_GF2NRSC,
    EC_TYPE_GF2NFFTRS,
    EC_TYPE_END,
};

std::map<int, std::string> ec_desc = {
    {EC_TYPE_ALL, "All available Reed-solomon codes"},
    {EC_TYPE_GF2NRSV, "Classical Vandermonde Reed-solomon codes over GF(2^n)"},
    {EC_TYPE_GF2NRSC, "Classical Cauchy Reed-solomon codes over GF(2^n)"},
    {EC_TYPE_GF2NFFTRS, "Reed-solomon codes over GF(2^n) using FFT"},
    {EC_TYPE_GF2NFFTADDRS,
     "Reed-solomon codes over GF(2^n) using additive FFT"},
    {EC_TYPE_GFPFFTRS, "Reed-solomon codes over GF(p) using FFT"},
    {EC_TYPE_FNTRS, "Reed-solomon codes over GF(p = Fermat number) using FFT"},
    {EC_TYPE_NGFF4RS,
     "Reed-solomon codes over GF(65537) using FFT on pack of codewords"},
};

std::map<int, std::string> ec_desc_short = {
    {EC_TYPE_ALL, "all"},
    {EC_TYPE_GF2NRSV, "gf2nrsv"},
    {EC_TYPE_GF2NRSC, "gf2nrsc"},
    {EC_TYPE_GF2NFFTRS, "gf2nfftrs"},
    {EC_TYPE_GF2NFFTADDRS, "gf2nfftaddrs"},
    {EC_TYPE_GFPFFTRS, "gfpfftrs"},
    {EC_TYPE_FNTRS, "fntrs"},
    {EC_TYPE_NGFF4RS, "ngff4rs"},
};

enum gf2nrs_type {
    VANDERMONDE = 0,
    CAUCHY,
};

enum errors {
    ERR_COMPT_WORD_SIZE_T = -4,
    ERR_WORD_SIZE,
    ERR_COMPT_CODE_LEN_T,
    ERR_FEC_TYPE_NOT_SUPPORTED,
    ERR_SCENARIO_TYPE_NOT_SUPPORTED,
    ERR_FAILED_CHUNK,
    ERR_FAILED_REPAIR_CHUNK,
    ERR_T_NOT_SUPPORTED,
};

std::map<int, std::string> errors_desc = {
    {ERR_COMPT_WORD_SIZE_T, "Word size and type T is not compatible"},
    {ERR_WORD_SIZE, "Word size is incorrect"},
    {ERR_COMPT_CODE_LEN_T, "Code length is too long vs. type T"},
    {ERR_FEC_TYPE_NOT_SUPPORTED, "Fec type is not recognised"},
    {ERR_SCENARIO_TYPE_NOT_SUPPORTED, "Scenario type is not recognised"},
    {ERR_FAILED_CHUNK, "ERROR: Chunks are incorrect"},
    {ERR_FAILED_REPAIR_CHUNK, "ERROR: Repaired chunks are incorrect"},
    {ERR_T_NOT_SUPPORTED, "Type T is not supported"},
};

std::map<std::string, ec_type> fec_type_map = {
    {"all", EC_TYPE_ALL},
    {"gf2nrsv", EC_TYPE_GF2NRSV},
    {"gf2nrsc", EC_TYPE_GF2NRSC},
    {"gf2nfftrs", EC_TYPE_GF2NFFTRS},
    {"gf2nfftaddrs", EC_TYPE_GF2NFFTADDRS},
    {"gfpfftrs", EC_TYPE_GFPFFTRS},
    {"fntrs", EC_TYPE_FNTRS},
    {"ngff4rs", EC_TYPE_NGFF4RS},
};

enum scenario_type {
    ENC_ONLY = 0,
    DEC_ONLY,
    ENC_DEC,
};

std::map<std::string, scenario_type> sce_type_map = {
    {"enc_only", ENC_ONLY},
    {"dec_only", DEC_ONLY},
    {"enc_dec", ENC_DEC},
};

std::map<int, std::string> sce_desc = {
    {ENC_ONLY, "Only encodings"},
    {DEC_ONLY, "Encode once and many decodings"},
    {ENC_DEC, "Encodings and decodings"},
};

std::map<int, std::string> sce_desc_short = {
    {ENC_ONLY, "enc"},
    {DEC_ONLY, "dec"},
    {ENC_DEC, "enc & dec"},
};

struct Params_t {
    ec_type fec_type = EC_TYPE_ALL;
    size_t word_size = 2;
    int k = 3;
    int m = 2;
    bool operation_on_packet = true;
    size_t pkt_size = 1024;
    size_t chunk_size = 512;
    uint32_t samples_nb = 100;
    int extra_param = -1;
    int sizeof_T = -1;
    scenario_type sce_type = ENC_DEC;
    uint32_t threads_nb = 4;
    // 0: show only params + speed
    // 1: show header + params + speed
    // 2: full show
    int compact_print = 1;
    std::string* header = nullptr;

    void print()
    {
        std::cout << "\n--------------- Parameters ----------------\n";
        std::cout << "FEC type:             " << ec_desc[fec_type] << std::endl;
        std::cout << "Scenario benchmark:   " << sce_desc[sce_type]
                  << std::endl;
        std::cout << "Word size:            " << word_size << std::endl;
        std::cout << "Number of data:       " << k << std::endl;
        std::cout << "Number of parity:     " << m << std::endl;
        std::cout << "Packet size           " << pkt_size << std::endl;
        std::cout << "Chunk size:           " << chunk_size << std::endl;
        std::cout << "Operations on packet: " << operation_on_packet
                  << std::endl;
        std::cout << "Number of samples:    " << samples_nb << std::endl;
        std::cout << "Number of threads:    " << threads_nb << std::endl;
        if (sizeof_T > -1)
            std::cout << "Size of integer type: " << sizeof_T << std::endl;
        if (extra_param > -1)
            std::cout << "Extra parameter:      " << extra_param << std::endl;
        std::cout << "-------------------------------------------\n";
    }

    void show_ec_desc()
    {
        std::cout << "\nBenchmarking for " << ec_desc[fec_type] << std::endl;
    }

    std::string* get_header()
    {
        if (header != nullptr)
            return header;

        std::string begin("| ");
        std::string end(" |\n");

        std::string content =
            begin + "         FEC" + "| scenario  " + "|       k" + "|       m"
            + "|       w" + "| T size " + "| chunk size  " + "| packet size "
            + "| #samples " + "| #threads "
            + "| Encoding   lat (us)       throuput (MB/s) "
            + "| Decoding   lat (us)       throuput (MB/s) " + end;

        std::string bare(content.length() - begin.length() - end.length(), '-');
        bare = begin + bare + end;

        header = new std::string();
        header->append(bare + content + bare);

        return header;
    }

    void show_header()
    {
        std::cout << get_header()->c_str();
    }

    void show_params()
    {
        std::cout << "  " << std::setw(12) << ec_desc_short[fec_type];
        std::cout << "  " << std::setw(10) << sce_desc_short[sce_type];
        std::cout << "  " << std::setw(7) << k;
        std::cout << "  " << std::setw(7) << m;
        std::cout << "  " << std::setw(7) << word_size;
        std::cout << "  " << std::setw(7) << sizeof_T;
        std::cout << "  " << std::setw(12) << chunk_size;
        if (operation_on_packet)
            std::cout << "  " << std::setw(12) << pkt_size;
        else
            std::cout << "  " << std::setw(12) << "-";
        std::cout << "  " << std::setw(9) << samples_nb;
        std::cout << "  " << std::setw(9) << threads_nb;
    }

    void show_speed(double avg, double std_dev, double thrput)
    {
        std::cout << "  " << std::setw(8) << avg << "+/-" << std::setw(8)
                  << std_dev << std::setw(20) << thrput;
    }

    void show_end()
    {
        std::cout << "\n";
    }

    void get_sizeof_T()
    {
        if (sizeof_T == -1) {
            if (fec_type == EC_TYPE_NGFF4RS) {
                sizeof_T = (((word_size - 1) / 2) + 1) * 4;
            } else if (fec_type == EC_TYPE_GFPFFTRS && word_size == 4) {
                sizeof_T = 8;
            } else {
                sizeof_T = (((word_size - 1) / 4) + 1) * 4;
            }
        }
    }
};

template <typename T>
class Benchmark {
  public:
    Benchmark(Params_t* params);
    ~Benchmark();
    bool enc_only();
    bool dec_only();
    bool enc_dec();

  private:
    int k;
    int m;
    int n;
    int n_c;
    size_t word_size;
    ec_type fec_type;
    bool operation_on_packet;
    int extra_param;
    size_t pkt_size;
    size_t chunk_size;
    uint32_t samples_nb;
    PRNG* prng = nullptr;
    nttec::fec::FEC<T>* fec = nullptr;
    Params_t* params = nullptr;

    bool systematic_ec = false;

    std::vector<int>* c_chunks_id = nullptr;

    std::vector<uint8_t*>* d_chunks = nullptr;
    std::vector<uint8_t*>* c_chunks = nullptr;
    std::vector<uint8_t*>* r_chunks = nullptr;

    std::vector<istreambuf<char>*>* d_istreambufs = nullptr;
    std::vector<istreambuf<char>*>* c_istreambufs = nullptr;
    std::vector<ostreambuf<char>*>* c_ostreambufs = nullptr;
    std::vector<ostreambuf<char>*>* r_ostreambufs = nullptr;

    Stats_t* enc_stats = nullptr;
    Stats_t* dec_stats = nullptr;

    // streams of data chunks
    std::vector<std::istream*>* d_streams = nullptr;
    // streams of coded chunks
    std::vector<std::ostream*>* c_streams = nullptr;
    // streams of available chunks
    std::vector<std::istream*>* a_streams = nullptr;
    // streams of repair chunks
    std::vector<std::ostream*>* r_streams = nullptr;
    // propos vector
    std::vector<nttec::KeyValue*>* c_propos = nullptr;

    int init();
    int check_params();
    void gen_data();
    bool check(std::vector<uint8_t*>* chunks);
    bool compare(std::vector<uint8_t*>* arr1, std::vector<uint8_t*>* arr2);
    void dump(const char* name, std::vector<uint8_t*>* chunks);
    void dump_chunk(const char* name, uint8_t* chunk);
    void reset_d_streams();
    void reset_c_streams();
    void reset_a_streams();
    void reset_r_streams();
    void get_avail_chunks(
        std::vector<std::istream*>* avail_d_chunks,
        std::vector<std::istream*>* avail_c_chunks,
        std::vector<nttec::KeyValue*>* avail_c_props);
    bool encode();
    bool decode();
    void show(Stats_t* stats);
    void show(Stats_t* stats1, Stats_t* stats2);
};

#endif
