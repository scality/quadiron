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
  EC_TYPE_UNDEF = 0,
  EC_TYPE_GF2NRSV,
  EC_TYPE_GF2NRSC,
  EC_TYPE_GF2NFFTRS,
  EC_TYPE_GF2NFFTADDRS,
  EC_TYPE_GFPFFTRS,
  EC_TYPE_FNTRS,
  EC_TYPE_NGFF4RS,
};

std::map<int, std::string> ec_desc = {
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
};

std::map<int, std::string> errors_desc = {
  { ERR_COMPT_WORD_SIZE_T, "Word size and type T is not compatible" },
  { ERR_WORD_SIZE, "Word size is incorrect" },
  { ERR_COMPT_CODE_LEN_T, "Code length is too long vs. type T" },
  { ERR_FEC_TYPE_NOT_SUPPORTED, "Fec type is not recognised" },
};

std::map<std::string, ec_type> fec_type_map = {
  { "gf2nrsv", EC_TYPE_GF2NRSV },
  { "gf2nrsc", EC_TYPE_GF2NRSC },
  { "gf2nfftrs", EC_TYPE_GF2NFFTRS },
  { "gf2nfftaddrs", EC_TYPE_GF2NFFTADDRS },
  { "gfpfftrs", EC_TYPE_GFPFFTRS },
  { "fntrs", EC_TYPE_FNTRS },
  { "ngff4rs", EC_TYPE_NGFF4RS },
};

struct Params_t
{
  ec_type fec_type = EC_TYPE_FNTRS;
  int word_size = 2;
  int k = 3;
  int m = 2;
  size_t chunk_size = 8;
  uint32_t samples_nb = 3;
  int extra_param = -1;
  int sizeof_T = 4;
  void print() {
    std::cout << "----- Parameters -----\n";
    std::cout << "fec_type: " << ec_desc[fec_type] << std::endl;
    std::cout << "word_size: " << word_size << std::endl;
    std::cout << "k: " << k << std::endl;
    std::cout << "m: " << m << std::endl;
    std::cout << "chunk_size: " << chunk_size << std::endl;
    std::cout << "samples_nb: " << samples_nb << std::endl;
    std::cout << "extra_param: " << extra_param << std::endl;
    std::cout << "sizeof_T: " << sizeof_T << std::endl;
    std::cout << "----------------------\n";
  }
};

struct Stats_t {
  uint64_t nb = 0;
  uint64_t sum = 0;
  uint64_t sum_2 = 0;
  double avg;
  double std_dev;

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

  void show(const char* name) {
    std::cout << name << ": " << avg << " +/- " << std_dev << std::endl;
  }
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
  int word_size;
  ec_type fec_type;
  int extra_param;
  size_t chunk_size;
  uint32_t samples_nb;
  PRNG *prng = nullptr;
  FEC<T> *fec = nullptr;

  std::vector<int> *c_chunks_id;

  std::vector<uint8_t*> *d_chunks = nullptr;
  std::vector<uint8_t*> *c_chunks = nullptr;
  std::vector<uint8_t*> *r_chunks = nullptr;

  Stats_t enc_stats;
  Stats_t dec_stats;

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
    this->chunk_size = params.chunk_size;
    this->samples_nb = params.samples_nb;
    if (params.extra_param > -1) {
      this->extra_param = params.extra_param;
    }

    int error;
    if ((error = this->check_params()) < 0) {
      std::cerr << errors_desc[error] << std::endl;
      exit(1);
    }

    if ((error = this->init()) < 0) {
      std::cerr << errors_desc[error] << std::endl;
      exit(1);
    }

    std::cout << "Bench init done" << std::endl;
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
      for (int i = 0; i < n; i++) {
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
      for (int i = 0; i < n; i++) {
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

    this->prng = new PRNG(time(NULL));

    // Allocate memory for data
    int i;
    d_chunks = new std::vector<uint8_t*>(k);
    c_chunks = new std::vector<uint8_t*>(n);
    r_chunks = new std::vector<uint8_t*>(k);

    for (i = 0; i < k; i++) {
      d_chunks->at(i) = new uint8_t[chunk_size];
    }
    for (i = 0; i < n; i++) {
      c_chunks->at(i) = new uint8_t[chunk_size];
    }
    for (i = 0; i < k; i++) {
      r_chunks->at(i) = new uint8_t[chunk_size];
    }

    // Allocate memory for streams
    d_streams = new std::vector<std::istream*>(k);
    c_streams = new std::vector<std::ostream*>(n);
    a_streams = new std::vector<std::istream*>(n);
    r_streams = new std::vector<std::ostream*>(k);
    c_propos = new std::vector<KeyValue*> (n);

    for (i = 0; i < k; i++) {
      d_streams->at(i) = new std::istream(
        new istreambuf<char>((char*)d_chunks->at(i), chunk_size));
    }

    for (i = 0; i < n; i++) {
      c_streams->at(i) = new std::ostream(
        new ostreambuf<char>((char*)c_chunks->at(i), chunk_size));
      c_propos->at(i) = new KeyValue();
    }

    for (i = 0; i < n; i++) {
      a_streams->at(i) = new std::istream(
        new istreambuf<char>((char*)c_chunks->at(i), chunk_size));
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

    return 1;
  }

  int check_params() {
    if (word_size <= 0) {
      return ERR_WORD_SIZE;
    }

    if (sizeof(T) < word_size) {
      return ERR_COMPT_WORD_SIZE_T;
    }

    uint64_t n_max;
    if (fec_type == EC_TYPE_FNTRS) {
      if (sizeof(T) <= word_size) {
        return ERR_COMPT_WORD_SIZE_T;
      }
      n_max = (1ULL << (8*word_size)) + 1;
    } else if (fec_type == EC_TYPE_NGFF4RS) {
      n_max = 65536;
    } else {
      n_max = (1ULL << (8*word_size));
    }
    if (n > n_max) {
      return ERR_COMPT_CODE_LEN_T;
    }

    return 1;
  }

  void gen_data() {
    for (int i = 0; i < k; i++) {
      prng->gen_chunk(d_chunks->at(i), chunk_size);
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

  void reset_d_streams() {
    for (int i = 0; i < k; i++) {
      d_streams->at(i)->rdbuf(
        new istreambuf<char>((char*)d_chunks->at(i), chunk_size));
    }
  }

  void reset_c_streams() {
    for (int i = 0; i < n; i++) {
      c_streams->at(i)->rdbuf(
        new ostreambuf<char>((char*)c_chunks->at(i), chunk_size));
      c_propos->at(i)->clear();
    }
  }

  void reset_a_streams() {
    for (int i = 0; i < n; i++) {
      a_streams->at(i)->rdbuf(
        new istreambuf<char>((char*)c_chunks->at(i), chunk_size));
    }
  }

  void reset_r_streams() {
    for (int i = 0; i < k; i++) {
      r_streams->at(i)->rdbuf(
        new ostreambuf<char>((char*)r_chunks->at(i), chunk_size));
    }
  }

  void get_avail_chunks(std::vector<std::istream*> *avail_c_chunks,
                        std::vector<KeyValue*> *avail_c_props) {

    std::random_shuffle(c_chunks_id->begin(), c_chunks_id->end());

    int i;
    for (i = 0; i < k; i++) {
      int j = c_chunks_id->at(i);
      avail_c_chunks->at(j) = a_streams->at(j);
      avail_c_props->at(j) = c_propos->at(j);
    }
    for (; i < n; i++) {
      int j = c_chunks_id->at(i);
      avail_c_chunks->at(j) = nullptr;
      avail_c_props->at(j) = nullptr;
    }
  }

  bool encode() {
    // this operation is done per trail
    reset_d_streams();
    fec->encode_bufs(*d_streams, *c_streams, *c_propos);

    if (!check(d_chunks)) {
      std::cout << "check d_chunks failed\n";
      return false;
    }

    // update stats
    enc_stats.add(fec->total_enc_usec);
    fec->reset_stats_enc();

    // dump("d_chunks", d_chunks);
    // dump("c_chunks", c_chunks);

    return true;
  }

  bool decode() {
    std::vector<std::istream*> c_streams_shuffled(n, nullptr);
    std::vector<KeyValue*> c_props_shuffled(n, nullptr);

    get_avail_chunks(&c_streams_shuffled, &c_props_shuffled);

    // this operation is done per trail
    reset_a_streams();
    reset_r_streams();

    fec->decode_bufs(*d_streams, c_streams_shuffled, c_props_shuffled,
      *r_streams);

    // dump("r_chunks", r_chunks);

    if (!check(r_chunks)) {
      std::cout << "check r_chunks failed\n";
      return false;
    }

    // update stats
    dec_stats.add(fec->total_dec_usec);
    fec->reset_stats_dec();

    return true;
  }

public:
  bool enc_only() {
    enc_stats.begin();
    std::cout << samples_nb << " encodings..";
    // this operation is done once per benchmark
    gen_data();

    for (int i = 0; i < samples_nb; i++) {
      if (!encode())
        return false;
    }
    std::cout << "done\n";

    enc_stats.end();
    enc_stats.show("Enc (usec)");

    return true;
  }

  bool dec_only() {
    enc_stats.begin();
    dec_stats.begin();

    std::cout << "encode once then " << samples_nb << " decodings..";
    // this operation is done once per benchmark
    gen_data();

    if (!encode())
      return false;

    for (int i = 0; i < samples_nb; i++) {
      if (!decode())
        return false;
    }
    std::cout << "done\n";

    enc_stats.end();
    enc_stats.show("Enc (usec)");
    dec_stats.end();
    dec_stats.show("Dec (usec)");

    return true;
  }

  bool enc_dec() {
    enc_stats.begin();
    dec_stats.begin();

    std::cout << samples_nb << " encodings and decodings..";
    // this operation is done once per benchmark
    gen_data();

    for (int i = 0; i < samples_nb; i++) {
      if (!encode())
        return false;
      if (!decode())
        return false;
    }
    std::cout << "done\n";

    enc_stats.end();
    enc_stats.show("Enc (usec)");
    dec_stats.end();
    dec_stats.show("Dec (usec)");
    return true;
  }
};

void xusage()
{
  std::cerr << std::string("Usage: benchmark") +
    "[-e gf2nrsv|gf2nrsc|gf2nfftrs|gf2nfftaddrs|gfpfftrs|fntrs|ngff4rs]" +
    "[-w word_size][-k k][-m k][-c chunk_size][-s samples_nb]" +
    "[-x extra_param]\n";
  exit(1);
}

int main(int argc, char **argv)
{
  PRNG prng;
  Params_t params;
  int opt;

  while ((opt = getopt(argc, argv, "t:e:w:k:m:c:s:x:")) != -1) {
    switch (opt) {
    case 't':
      params.sizeof_T = atoi(optarg);
      break;
    case 'e':
      if (fec_type_map.find(optarg) == fec_type_map.end()) {
        assert("Input fec type is incorrect");
        xusage();
      }
      params.fec_type = fec_type_map[optarg];
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
    case 's':
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

  switch (params.sizeof_T) {
    case 4: {
      Benchmark<uint32_t> bench(params);
      bench.enc_dec();
      break;
    }
    case 8: {
      Benchmark<uint64_t> bench(params);
      bench.enc_dec();
      break;
    }
    case 16: {
      Benchmark<__uint128_t> bench(params);
      bench.enc_dec();
      break;
    }
    default:
      assert("T type is not supported");
  }

  return 0;
}
