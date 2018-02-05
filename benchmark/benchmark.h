/* -*- mode: c++ -*- */
// #pragma once

#include <sys/types.h>
#include <sys/stat.h>
#include <thread>

#include "ntl.h"

#include "iostreambuf.h"
#include "prng.h"
#include "stats.h"

template class FECGF2NRS<uint32_t>;
template class FECGF2NRS<uint64_t>;
template class FECGF2NRS<__uint128_t>;

template class FECGF2NFFTRS<uint32_t>;
template class FECGF2NFFTRS<uint64_t>;
template class FECGF2NFFTRS<__uint128_t>;

template class FECGF2NFFTADDRS<uint32_t>;
template class FECGF2NFFTADDRS<uint64_t>;
template class FECGF2NFFTADDRS<__uint128_t>;

template class FECGFPFFTRS<uint32_t>;
template class FECGFPFFTRS<uint64_t>;
template class FECGFPFFTRS<__uint128_t>;

template class FECFNTRS<uint32_t>;
template class FECFNTRS<uint64_t>;
template class FECFNTRS<__uint128_t>;

template class FECNGFF4RS<uint32_t>;
template class FECNGFF4RS<uint64_t>;
template class FECNGFF4RS<__uint128_t>;

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
  { EC_TYPE_ALL, "All available Reed-solomon codes" },
  { EC_TYPE_GF2NRSV, "Classical Vandermonde Reed-solomon codes over GF(2^n)" },
  { EC_TYPE_GF2NRSC, "Classical Cauchy Reed-solomon codes over GF(2^n)" },
  { EC_TYPE_GF2NFFTRS, "Reed-solomon codes over GF(2^n) using FFT" },
  { EC_TYPE_GF2NFFTADDRS, "Reed-solomon codes over GF(2^n) using additive FFT"},
  { EC_TYPE_GFPFFTRS, "Reed-solomon codes over GF(p) using FFT" },
  { EC_TYPE_FNTRS, "Reed-solomon codes over GF(p = Fermat number) using FFT" },
  { EC_TYPE_NGFF4RS,
    "Reed-solomon codes over GF(65537) using FFT on pack of codewords" },
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
  { ERR_COMPT_WORD_SIZE_T, "Word size and type T is not compatible" },
  { ERR_WORD_SIZE, "Word size is incorrect" },
  { ERR_COMPT_CODE_LEN_T, "Code length is too long vs. type T" },
  { ERR_FEC_TYPE_NOT_SUPPORTED, "Fec type is not recognised" },
  { ERR_SCENARIO_TYPE_NOT_SUPPORTED, "Scenario type is not recognised" },
  { ERR_FAILED_CHUNK, "ERROR: Chunks are incorrect" },
  { ERR_FAILED_REPAIR_CHUNK, "ERROR: Repaired chunks are incorrect" },
  { ERR_T_NOT_SUPPORTED, "Type T is not supported" },
};

std::map<std::string, ec_type> fec_type_map = {
  { "all", EC_TYPE_ALL },
  { "gf2nrsv", EC_TYPE_GF2NRSV },
  { "gf2nrsc", EC_TYPE_GF2NRSC },
  { "gf2nfftrs", EC_TYPE_GF2NFFTRS },
  { "gf2nfftaddrs", EC_TYPE_GF2NFFTADDRS },
  { "gfpfftrs", EC_TYPE_GFPFFTRS },
  { "fntrs", EC_TYPE_FNTRS },
  { "ngff4rs", EC_TYPE_NGFF4RS },
};

enum scenario_type {
  ENC_ONLY = 0,
  DEC_ONLY,
  ENC_DEC,
};

std::map<std::string, scenario_type> sce_type_map = {
  { "enc_only", ENC_ONLY },
  { "dec_only", DEC_ONLY },
  { "enc_dec", ENC_DEC },
};

std::map<int, std::string> sce_desc = {
  { ENC_ONLY, "Only encodings" },
  { DEC_ONLY, "Encode once and many decodings" },
  { ENC_DEC, "Encodings and decodings" },
};

struct Params_t
{
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

  void print() {
    std::cout << "\n--------------- Parameters ----------------\n";
    std::cout << "FEC type:             " << ec_desc[fec_type] << std::endl;
    std::cout << "Scenario benchmark:   " << sce_desc[sce_type] << std::endl;
    std::cout << "Word size:            " << word_size << std::endl;
    std::cout << "Number of data:       " << k << std::endl;
    std::cout << "Number of parity:     " << m << std::endl;
    std::cout << "Packet size           " << pkt_size << std::endl;
    std::cout << "Chunk size:           " << chunk_size << std::endl;
    std::cout << "Operations on packet: " << operation_on_packet << std::endl;
    std::cout << "Number of samples:    " << samples_nb << std::endl;
    std::cout << "Number of threads:    " << threads_nb << std::endl;
    if (sizeof_T > -1)
      std::cout << "Size of integer type: " << sizeof_T << std::endl;
    if (extra_param > -1)
      std::cout << "Extra parameter:      " << extra_param << std::endl;
    std::cout << "-------------------------------------------\n";
  }

  void get_sizeof_T() {
    if (sizeof_T == -1) {
      if (fec_type == EC_TYPE_NGFF4RS) {
        sizeof_T = (((word_size-1) / 2) + 1) * 4;
      } else if (fec_type == EC_TYPE_GFPFFTRS && word_size == 4) {
        sizeof_T = 8;
      } else {
        sizeof_T = (((word_size-1) / 4) + 1) * 4;
      }
      std::cout << "Size of integer type: " << sizeof_T << std::endl;
    }
  }
};

template<typename T>
class Benchmark
{
 public:
  Benchmark(Params_t params);
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
  PRNG *prng = nullptr;
  FEC<T> *fec = nullptr;

  bool systematic_ec = false;

  std::vector<int> *c_chunks_id = nullptr;

  std::vector<uint8_t*> *d_chunks = nullptr;
  std::vector<uint8_t*> *c_chunks = nullptr;
  std::vector<uint8_t*> *r_chunks = nullptr;

  std::vector<istreambuf<char>*> *d_istreambufs = nullptr;
  std::vector<istreambuf<char>*> *c_istreambufs = nullptr;
  std::vector<ostreambuf<char>*> *c_ostreambufs = nullptr;
  std::vector<ostreambuf<char>*> *r_ostreambufs = nullptr;

  Stats_t *enc_stats = nullptr;
  Stats_t *dec_stats = nullptr;

  // streams of data chunks
  std::vector<std::istream*> *d_streams = nullptr;
  // streams of coded chunks
  std::vector<std::ostream*> *c_streams = nullptr;
  // streams of available chunks
  std::vector<std::istream*> *a_streams = nullptr;
  // streams of repair chunks
  std::vector<std::ostream*> *r_streams = nullptr;
  // propos vector
  std::vector<KeyValue*> *c_propos = nullptr;

  int init();
  int check_params();
  void gen_data();
  bool check(std::vector<uint8_t*> *chunks);
  bool compare(std::vector<uint8_t*> *arr1, std::vector<uint8_t*> *arr2);
  void dump(const char* name, std::vector<uint8_t*> *chunks);
  void dump_chunk(const char* name, uint8_t *chunk);
  void reset_d_streams();
  void reset_c_streams();
  void reset_a_streams();
  void reset_r_streams();
  void get_avail_chunks(std::vector<std::istream*> *avail_d_chunks,
                       std::vector<std::istream*> *avail_c_chunks,
                       std::vector<KeyValue*> *avail_c_props);
  bool encode();
  bool decode();
};
