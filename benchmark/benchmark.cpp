/* -*- mode: c++ -*- */

#include "ntl.h"
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

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

// CRC-32C (iSCSI) polynomial in reversed bit order
#define POLY 0x82f63b78

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
  ERR_T_NOT_SUPPORTED,
};

std::map<int, std::string> errors_desc = {
  { ERR_COMPT_WORD_SIZE_T, "Word size and type T is not compatible" },
  { ERR_WORD_SIZE, "Word size is incorrect" },
  { ERR_COMPT_CODE_LEN_T, "Code length is too long vs. type T" },
  { ERR_FEC_TYPE_NOT_SUPPORTED, "Fec type is not recognised" },
  { ERR_SCENARIO_TYPE_NOT_SUPPORTED, "Scenario type is not recognised" },
  { ERR_FAILED_CHUNK, "ERROR: Repaired chunks are incorrect" },
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
  int word_size = 2;
  int k = 3;
  int m = 2;
  size_t chunk_size = 512;
  uint32_t samples_nb = 100;
  int extra_param = -1;
  int sizeof_T = -1;
  scenario_type sce_type = ENC_DEC;

  void print() {
    std::cout << "\n--------------- Parameters ----------------\n";
    std::cout << "FEC type:             " << ec_desc[fec_type] << std::endl;
    std::cout << "Scenario benchmark:   " << sce_desc[sce_type] << std::endl;
    std::cout << "Word size:            " << word_size << std::endl;
    std::cout << "Number of data:       " << k << std::endl;
    std::cout << "Number of parity:     " << m << std::endl;
    std::cout << "Chunk size:           " << chunk_size << std::endl;
    std::cout << "Number of samples:    " << samples_nb << std::endl;
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
      } else {
        sizeof_T = (((word_size-1) / 4) + 1) * 4;
      }
      std::cout << "Size of integer type: " << sizeof_T << std::endl;
    }
  }
};

class Stats_t {
public:
  Stats_t(const std::string& str, size_t work_load) {
    this->name = str;
    this->work_load = work_load;

    this->nb = 0;
    this->sum = 0;
    this->sum_2 = 0;
  }

  void begin() {
    nb = 0;
    sum = 0;
    sum_2 = 0;
  }

  void add(uint64_t val) {
    nb++;
    sum += val;
    sum_2 += (val * val);
  }

  void end() {
    avg = (double)sum / (double)nb;
    std_dev = sqrt((double)sum_2 / (double)nb - avg * avg);
  }

  void show() {
    std::cout << name << ":\tLatency(us) " << avg << " +/- " << std_dev;
    std::cout << "\t\tThroughput " << work_load/avg << " (MB/s)" << std::endl;
  }

private:
  uint64_t nb;
  uint64_t sum;
  uint64_t sum_2;
  double avg;
  double std_dev;
  size_t work_load;
  std::string name;
};

template <typename CharT>
class ostreambuf : public std::basic_streambuf<CharT>
{
public:
  ostreambuf(CharT* buffer, size_t buf_len)
  {
    std::basic_streambuf<CharT>::setp(buffer, buffer + buf_len);
  }
};

template <typename CharT>
class istreambuf : public std::basic_streambuf<CharT>
{
public:
  istreambuf(CharT* buffer, size_t buf_len)
  {
    std::basic_streambuf<CharT>::setg(buffer, buffer, buffer + buf_len);
  }
};

/*
 *  From http://www.pcg-random.org/
 *
 * *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
 * Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
 */
class PRNG {

private:
  uint64_t state;
  uint64_t inc;

public:
  PRNG(uint64_t initstate = 0, uint64_t initseq = 0) {
    this->state = initstate;
    this->inc = initseq;
  }

  void srandom(uint64_t initstate, uint64_t initseq) {
    state = 0U;
    inc = (initseq << 1u) | 1u;
    rand();
    state += initstate;
    rand();
  }

  uint32_t _rand() {
    uint64_t oldstate = state;
    state = oldstate * 6364136223846793005ULL + inc;
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
  }

  void gen_chunk(void* chunk, size_t size) {
    uint8_t* buffer = (uint8_t*)chunk;
    if (size < 4) {
      buffer[0] = (uint8_t)_rand();
      for (size_t i = 1; i < size; i++) {
        buffer[i] = buffer[0];
      }
    } else {
      uint32_t crc = 0;
      crc = ~crc;
      for (size_t i = 4; i < size; i++) {
        uint8_t val = (uint8_t)_rand();
        buffer[i] = val;
        crc ^= val;
        for (int k = 0; k < 8; k++)
          crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
      }
      *(uint32_t*)buffer = ~crc;
    }
  }

  bool check_chunk(void* chunk, size_t size) {
    uint8_t* buffer = (uint8_t*)chunk;
    if (size < 4) {
      uint8_t val = buffer[0];
      for (size_t i = 1; i < size; i++) {
        if (buffer[i] != val)
          return false;
      }
    } else {
      uint32_t crc = 0;
      crc = ~crc;
      for (size_t i = 4; i < size; i++) {
        crc ^= buffer[i];
        for (int k = 0; k < 8; k++)
          crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
      }
      uint32_t _crc = *(uint32_t*)buffer;
      if (_crc != ~crc)
        return false;
    }
    return true;
  }
};

template<typename T>
class Benchmark
{
private:
  int k;
  int m;
  int n;
  int n_c;
  int word_size;
  ec_type fec_type;
  int extra_param;
  size_t chunk_size;
  uint32_t samples_nb;
  PRNG *prng = nullptr;
  FEC<T> *fec = nullptr;

  bool systematic_ec = false;

  std::vector<int> *c_chunks_id;

  std::vector<uint8_t*> *d_chunks = nullptr;
  std::vector<uint8_t*> *c_chunks = nullptr;
  std::vector<uint8_t*> *r_chunks = nullptr;

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

public:
  Benchmark(Params_t params) {
    this->fec_type = params.fec_type;
    this->word_size = params.word_size;
    this->k = params.k;
    this->m = params.m;
    this->n = params.k + params.m;
    this->n_c = this->n;  // if systematic, it will be updated in init()
    this->chunk_size = params.chunk_size;
    this->samples_nb = params.samples_nb;
    if (params.extra_param > -1) {
      this->extra_param = params.extra_param;
    }

    int error;
    if ((error = this->check_params()) < 0) {
      std::cerr << errors_desc[error] << std::endl;
      throw std::invalid_argument(errors_desc[error]);
    }

    if ((error = this->init()) < 0) {
      std::cerr << errors_desc[error] << std::endl;
      throw std::invalid_argument(errors_desc[error]);
    }
  }

  ~Benchmark() {
    if (fec != nullptr) delete fec;
    if (prng != nullptr) delete prng;

    if (d_streams != nullptr) {
      for (int i = 0; i < k; i++) {
        delete d_streams->at(i);
      }
      delete d_streams;
    }
    if (c_streams != nullptr) {
      for (int i = 0; i < n_c; i++) {
        delete c_streams->at(i);
        delete c_propos->at(i);
      }
      delete c_streams;
      delete c_propos;
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

    if (d_chunks != nullptr) {
      for (int i = 0; i < k; i++) {
        delete[] d_chunks->at(i);
      }
      delete d_chunks;
    }
    if (c_chunks != nullptr) {
      for (int i = 0; i < n_c; i++) {
        delete[] c_chunks->at(i);
      }
      delete c_chunks;
    }
    if (r_chunks != nullptr) {
      for (int i = 0; i < k; i++) {
        delete[] r_chunks->at(i);
      }
      delete r_chunks;
    }

    if (enc_stats != nullptr)
      delete enc_stats;
    if (dec_stats != nullptr)
      delete dec_stats;
  }

private:
  int init() {
    switch (fec_type) {
      case EC_TYPE_GF2NRSV:
        fec = new FECGF2NRS<T>(word_size, k, m, FECGF2NRS<T>::VANDERMONDE);
        break;
      case EC_TYPE_GF2NRSC:
        fec = new FECGF2NRS<T>(word_size, k, m, FECGF2NRS<T>::CAUCHY);
        break;
      case EC_TYPE_GF2NFFTRS:
        fec = new FECGF2NFFTRS<T>(word_size, k, m);
        break;
      case EC_TYPE_GF2NFFTADDRS:
        fec = new FECGF2NFFTADDRS<T>(word_size, k, m);
        break;
      case EC_TYPE_GFPFFTRS:
        fec = new FECGFPFFTRS<T>(word_size, k, m);
        break;
      case EC_TYPE_NGFF4RS:
        fec = new FECNGFF4RS<T>(word_size, k, m);
        break;
      case EC_TYPE_FNTRS:
        fec = new FECFNTRS<T>(word_size, k, m);
        break;
      default:
        return ERR_FEC_TYPE_NOT_SUPPORTED;
    }

    this->systematic_ec = (fec->type == FEC<T>::TYPE_1);
    if (this->systematic_ec) {
      this->n_c = this->m;
    }

    this->prng = new PRNG(time(NULL));

    // Allocate memory for data
    int i;
    d_chunks = new std::vector<uint8_t*>(k);
    c_chunks = new std::vector<uint8_t*>(n_c);
    r_chunks = new std::vector<uint8_t*>(k);

    for (i = 0; i < k; i++) {
      d_chunks->at(i) = new uint8_t[chunk_size];
    }
    for (i = 0; i < n_c; i++) {
      c_chunks->at(i) = new uint8_t[chunk_size];
    }
    for (i = 0; i < k; i++) {
      r_chunks->at(i) = new uint8_t[chunk_size];
    }

    // Allocate memory for streams
    d_streams = new std::vector<std::istream*>(k);
    c_streams = new std::vector<std::ostream*>(n_c);
    a_streams = new std::vector<std::istream*>(n);
    r_streams = new std::vector<std::ostream*>(k);
    c_propos = new std::vector<KeyValue*> (n_c);

    for (i = 0; i < k; i++) {
      d_streams->at(i) = new std::istream(
        new istreambuf<char>((char*)d_chunks->at(i), chunk_size));
    }

    for (i = 0; i < n_c; i++) {
      c_streams->at(i) = new std::ostream(
        new ostreambuf<char>((char*)c_chunks->at(i), chunk_size));
      c_propos->at(i) = new KeyValue();
    }

    if (systematic_ec) {
      for (i = 0; i < k; i++) {
        a_streams->at(i) = new std::istream(
          new istreambuf<char>((char*)d_chunks->at(i), chunk_size));
      }
      for (; i < n; i++) {
        a_streams->at(i) = new std::istream(
          new istreambuf<char>((char*)c_chunks->at(i-k), chunk_size));
      }
    } else {
      for (i = 0; i < n; i++) {
        a_streams->at(i) = new std::istream(
          new istreambuf<char>((char*)c_chunks->at(i), chunk_size));
      }
    }

    for (i = 0; i < k; i++) {
      r_streams->at(i) = new std::ostream(
        new ostreambuf<char>((char*)r_chunks->at(i), chunk_size));
    }

    // init index of chunks
    c_chunks_id = new std::vector<int>();
    for (i = 0; i < n; i++){
      c_chunks_id->push_back(i);
    }

    this->enc_stats = new Stats_t("Encode", chunk_size*k);
    this->dec_stats = new Stats_t("Decode", chunk_size*k);
    return 1;
  }

  int check_params() {
    if (word_size <= 0) {
      return ERR_WORD_SIZE;
    }

    if (fec_type == EC_TYPE_NGFF4RS) {
      if (word_size < 2) {
        return ERR_WORD_SIZE;
      }
    }

    if (sizeof(T) < word_size) {
      return ERR_COMPT_WORD_SIZE_T;
    }
    if (fec_type == EC_TYPE_FNTRS ||
        fec_type == EC_TYPE_GFPFFTRS) {
      if (sizeof(T) <= word_size) {
        return ERR_COMPT_WORD_SIZE_T;
      }
    }
    if (fec_type == EC_TYPE_NGFF4RS) {
      if (sizeof(T) < 2 * word_size) {
        return ERR_COMPT_WORD_SIZE_T;
      }
    }

    uint64_t n_max;
    if (fec_type == EC_TYPE_FNTRS) {
      n_max = (1ULL << (8*word_size)) + 1;
    } else if (fec_type == EC_TYPE_NGFF4RS) {
      n_max = 65536;
    } else {
      n_max = (1ULL << (8*word_size));
    }
    if (n > n_max) {
      return ERR_COMPT_CODE_LEN_T;
    }

    // adjust chunk_size
    if (chunk_size % word_size > 0)
      chunk_size = ((chunk_size - 1)/word_size + 1)*word_size;

    return 1;
  }

  void gen_data() {
    for (int i = 0; i < k; i++) {
      prng->gen_chunk(d_chunks->at(i), chunk_size);
    }
    if (!check(d_chunks)) {
      std::cerr << errors_desc[ERR_FAILED_CHUNK] << std::endl;
      throw errors_desc[ERR_FAILED_CHUNK];
    }
  }

  bool check(std::vector<uint8_t*> *chunks) {
    for (std::vector<int>::size_type i = 0; i < chunks->size(); i++) {
      if (!prng->check_chunk(chunks->at(i), chunk_size))
        return false;
    }
    return true;
  }

  void dump(const char* name, std::vector<uint8_t*> *chunks) {
    std::cout << name << ": ";
    for (std::vector<int>::size_type i = 0; i < chunks->size(); i++) {
      uint8_t *chunk = chunks->at(i);
      for (int j = 0; j < chunk_size; j++) {
        std::cout << unsigned(chunk[j]) << " ";
      }
      std::cout << "|";
    }
    std::cout << std::endl;
  }

  void dump_chunk(const char* name, uint8_t *chunk) {
    std::cout << name << ": ";
    for (int j = 0; j < chunk_size; j++) {
      std::cout << unsigned(chunk[j]) << " ";
    }
    std::cout << std::endl;
  }

  void reset_d_streams() {
    for (int i = 0; i < k; i++) {
      d_streams->at(i)->rdbuf(
        new istreambuf<char>((char*)d_chunks->at(i), chunk_size));
    }
  }

  void reset_c_streams() {
    for (int i = 0; i < n_c; i++) {
      c_streams->at(i)->rdbuf(
        new ostreambuf<char>((char*)c_chunks->at(i), chunk_size));
      c_propos->at(i)->clear();
    }
  }

  void reset_a_streams() {
    int i;
    if (systematic_ec) {
      for (i = 0; i < k; i++) {
        a_streams->at(i)->rdbuf(
          new istreambuf<char>((char*)d_chunks->at(i), chunk_size));
      }
      for (; i < n; i++) {
        a_streams->at(i)->rdbuf(
          new istreambuf<char>((char*)c_chunks->at(i-k), chunk_size));
      }
    } else {
      for (i = 0; i < n; i++) {
        a_streams->at(i)->rdbuf(
          new istreambuf<char>((char*)c_chunks->at(i), chunk_size));
      }
    }
  }

  void reset_r_streams() {
    for (int i = 0; i < k; i++) {
      r_streams->at(i)->rdbuf(
        new ostreambuf<char>((char*)r_chunks->at(i), chunk_size));
    }
  }

  void get_avail_chunks(std::vector<std::istream*> *avail_d_chunks,
                        std::vector<std::istream*> *avail_c_chunks,
                        std::vector<KeyValue*> *avail_c_props) {

    std::random_shuffle(c_chunks_id->begin(), c_chunks_id->end());

    int i;
    for (i = 0; i < k; i++) {
      avail_d_chunks->at(i) = nullptr;
    }
    for (i = 0; i < n_c; i++) {
      avail_c_chunks->at(i) = nullptr;
      avail_c_props->at(i) = nullptr;
    }
    if (systematic_ec) {
      for (i = 0; i < k; i++) {
        int j = c_chunks_id->at(i);
        if (j < k) {
          avail_d_chunks->at(j) = a_streams->at(j);
        } else {
          avail_c_chunks->at(j-k) = a_streams->at(j);
          avail_c_props->at(j-k) = c_propos->at(j-k);
        }
      }
    } else {
      for (i = 0; i < k; i++) {
        int j = c_chunks_id->at(i);
        avail_c_chunks->at(j) = a_streams->at(j);
        avail_c_props->at(j) = c_propos->at(j);
      }
    }
  }

  bool encode() {
    // this operation is done per trail
    reset_d_streams();
    fec->encode_bufs(*d_streams, *c_streams, *c_propos);

    // update stats
    enc_stats->add(fec->total_enc_usec);

    // dump("d_chunks", d_chunks);
    // dump("c_chunks", c_chunks);

    return true;
  }

  bool decode() {
    std::vector<std::istream*> d_streams_shuffled(k, nullptr);
    std::vector<std::istream*> c_streams_shuffled(n_c, nullptr);
    std::vector<KeyValue*> c_props_shuffled(n_c, nullptr);

    get_avail_chunks(&d_streams_shuffled, &c_streams_shuffled,
      &c_props_shuffled);

    // this operation is done per trail
    reset_a_streams();
    reset_r_streams();

    fec->decode_bufs(d_streams_shuffled, c_streams_shuffled, c_props_shuffled,
      *r_streams);

    // dump("r_chunks", r_chunks);

    if (!check(r_chunks)) {
      std::cerr << errors_desc[ERR_FAILED_CHUNK] << std::endl;
      return false;
    }

    // update stats
    dec_stats->add(fec->total_dec_usec);

    return true;
  }

public:
  bool enc_only() {
    enc_stats->begin();
    // this operation is done once per benchmark
    gen_data();

    for (int i = 0; i < samples_nb; i++) {
      if (!encode())
        return false;
    }

    enc_stats->end();
    enc_stats->show();

    return true;
  }

  bool dec_only() {
    enc_stats->begin();
    dec_stats->begin();

    // this operation is done once per benchmark
    gen_data();

    if (!encode())
      return false;

    for (int i = 0; i < samples_nb; i++) {
      if (!decode())
        return false;
    }

    enc_stats->end();
    enc_stats->show();
    dec_stats->end();
    dec_stats->show();

    return true;
  }

  bool enc_dec() {
    enc_stats->begin();
    dec_stats->begin();

    // this operation is done once per benchmark
    gen_data();

    for (int i = 0; i < samples_nb; i++) {
      if (!encode())
        return false;
      if (!decode())
        return false;
    }

    enc_stats->end();
    enc_stats->show();
    dec_stats->end();
    dec_stats->show();
    return true;
  }
};

void xusage()
{
  std::cerr << std::string("Usage: benchmark [options]\n") +
    "Options:\n" +
    "\t-e \tType of Reed-Solomon codes, either\n" +
    "\t\t\tgf2nrsv: " + ec_desc[EC_TYPE_GF2NRSV] + "\n" +
    "\t\t\tgf2nrsc: " + ec_desc[EC_TYPE_GF2NRSC] + "\n" +
    "\t\t\tgf2nfftrs: " + ec_desc[EC_TYPE_GF2NFFTRS] + "\n" +
    "\t\t\tgf2nfftaddrs: " + ec_desc[EC_TYPE_GF2NFFTADDRS] + "\n" +
    "\t\t\tgfpfftrs: " + ec_desc[EC_TYPE_GFPFFTRS] + "\n" +
    "\t\t\tfntrs: " + ec_desc[EC_TYPE_FNTRS] + "\n" +
    "\t\t\tngff4rs: " + ec_desc[EC_TYPE_NGFF4RS] + "\n" +
    "\t\t\tall: All available Reed-solomon codes\n" +
    "\t-s \tScenario for benchmark, either\n" +
    "\t\t\tenc_only: Only encodings\n" +
    "\t\t\tdec_only: Only decodings\n" +
    "\t\t\tenc_dec: Encodings and decodings\n" +
    "\t-w \tWord size (bytes)\n" +
    "\t-k \tNumber of data chunks\n" +
    "\t-m \tNumber of parity chunks\n" +
    "\t-c \tChunk size (bytes)\n" +
    "\t-n \tNumber of samples per operation\n" +
    "\t-t \tSize of used integer type, either " +
      "4, 8, 16 for uint32_t, uint64_t, __uint128_t\n" +
    "\t-x \tExtra parameter\n\n";
  exit(1);
}

template <typename T>
void run(Benchmark<T> *bench, Params_t params) {
  switch (params.sce_type) {
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

void run_scenario(Params_t params) {
  // get sizeof_T if necessary
  params.get_sizeof_T();

  switch (params.sizeof_T) {
    case 4: {
      try {
        Benchmark<uint32_t> bench(params);
        run(&bench, params);
      } catch(const std::exception& e) {
        return;
      }
      break;
    }
    case 8: {
      try {
        Benchmark<uint64_t> bench(params);
        run(&bench, params);
      } catch(const std::exception& e) {
        return;
      }
      break;
    }
    case 16: {
      try {
        Benchmark<__uint128_t> bench(params);
        run(&bench, params);
      } catch(const std::exception& e) {
        return;
      }
      break;
    }
    default:
      std::cerr << errors_desc[ERR_T_NOT_SUPPORTED] <<
        " T: " << params.sizeof_T << std::endl;
      exit(0);
  }
}

void run_benchmark(Params_t params) {
  if (params.fec_type == EC_TYPE_ALL) {
    for (int type = EC_TYPE_ALL + 1; type < EC_TYPE_END; type++) {
      params.fec_type = (ec_type)type;
      std::cout << "\nBenchmarking for " <<
        ec_desc[params.fec_type] << std::endl;
      run_scenario(params);
    }
  } else {
    run_scenario(params);
  }
}

int main(int argc, char **argv)
{
  PRNG prng;
  Params_t params;
  int opt;

  while ((opt = getopt(argc, argv, "t:e:w:k:m:c:n:s:x:")) != -1) {
    switch (opt) {
    case 't':
      params.sizeof_T = atoi(optarg);
      break;
    case 'e':
      if (fec_type_map.find(optarg) == fec_type_map.end()) {
        std::cerr << errors_desc[ERR_FEC_TYPE_NOT_SUPPORTED] << std::endl;
        xusage();
      }
      params.fec_type = fec_type_map[optarg];
      break;
    case 's':
      if (sce_type_map.find(optarg) == sce_type_map.end()) {
        std::cerr << errors_desc[ERR_SCENARIO_TYPE_NOT_SUPPORTED] << std::endl;
        xusage();
      }
      params.sce_type = sce_type_map[optarg];
      break;
    case 'w':
      params.word_size = atoi(optarg);
      break;
    case 'k':
      params.k = atoi(optarg);
      break;
    case 'm':
      params.m = atoi(optarg);
      break;
    case 'c':
      params.chunk_size = atoi(optarg);
      break;
    case 'n':
      params.samples_nb = atoi(optarg);
      break;
    case 'x':
      params.extra_param = atoi(optarg);
      break;
    default:
      xusage();
    }
  }

  params.print();

  run_benchmark(params);

  return 0;
}
