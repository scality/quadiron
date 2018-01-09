/* -*- mode: c++ -*- */

#include "benchmark.h"

template<typename T>
Benchmark<T>::Benchmark(Params_t params) {
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

template<typename T>
Benchmark<T>::~Benchmark() {
  if (fec != nullptr) delete fec;
  if (prng != nullptr) delete prng;
  if (c_chunks_id != nullptr) delete c_chunks_id;

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

template<typename T>
int Benchmark<T>::init() {
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

  // Allocate memory for iostreambufs
  d_istreambufs = new std::vector<istreambuf<char>*>(k);
  c_istreambufs = new std::vector<istreambuf<char>*>(n_c);
  c_ostreambufs = new std::vector<ostreambuf<char>*>(n_c);
  r_ostreambufs = new std::vector<ostreambuf<char>*>(k);
  for (i = 0; i < k; i++) {
    d_istreambufs->at(i) =
      new istreambuf<char>((char*)d_chunks->at(i), chunk_size);
  }
  for (i = 0; i < n_c; i++) {
    c_istreambufs->at(i) =
      new istreambuf<char>((char*)c_chunks->at(i), chunk_size);
    c_ostreambufs->at(i) =
      new ostreambuf<char>((char*)c_chunks->at(i), chunk_size);
  }
  for (i = 0; i < k; i++) {
    r_ostreambufs->at(i) =
      new ostreambuf<char>((char*)r_chunks->at(i), chunk_size);
  }

  // Allocate memory for streams
  d_streams = new std::vector<std::istream*>(k);
  c_streams = new std::vector<std::ostream*>(n_c);
  a_streams = new std::vector<std::istream*>(n);
  r_streams = new std::vector<std::ostream*>(k);
  c_propos = new std::vector<KeyValue*>(n_c);

  for (i = 0; i < k; i++) {
    d_streams->at(i) = new std::istream(d_istreambufs->at(i));
  }

  for (i = 0; i < n_c; i++) {
    c_streams->at(i) = new std::ostream(c_ostreambufs->at(i));
    c_propos->at(i) = new KeyValue();
  }

  if (systematic_ec) {
    for (i = 0; i < k; i++) {
      a_streams->at(i) = new std::istream(d_istreambufs->at(i));
    }
    for (; i < n; i++) {
      a_streams->at(i) = new std::istream(c_istreambufs->at(i-k));
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
  for (i = 0; i < n; i++){
    c_chunks_id->push_back(i);
  }

  this->enc_stats = new Stats_t("Encode", chunk_size*n_c);
  this->dec_stats = new Stats_t("Decode", chunk_size*k);
  return 1;
}

template<typename T>
int Benchmark<T>::check_params() {
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

  int wordsize_limit = _log2<T>(n) + 1;
  if (wordsize_limit > 8 * word_size) {
    return ERR_COMPT_CODE_LEN_T;
  }

  // adjust chunk_size
  if (chunk_size % word_size > 0)
    chunk_size = ((chunk_size - 1)/word_size + 1)*word_size;

  return 1;
}

template<typename T>
void Benchmark<T>::gen_data() {
  for (int i = 0; i < k; i++) {
    prng->gen_chunk(d_chunks->at(i), chunk_size);
  }
  if (!check(d_chunks)) {
    std::cerr << errors_desc[ERR_FAILED_CHUNK] << std::endl;
    throw errors_desc[ERR_FAILED_CHUNK];
  }
}

template<typename T>
bool Benchmark<T>::check(std::vector<uint8_t*> *chunks) {
  for (std::vector<int>::size_type i = 0; i < chunks->size(); i++) {
    if (!prng->check_chunk(chunks->at(i), chunk_size))
      return false;
  }
  return true;
}

template<typename T>
bool Benchmark<T>::compare(std::vector<uint8_t*> *arr1,
                           std::vector<uint8_t*> *arr2) {
  if (arr1->size() != arr2->size()) {
    std::cout << "Sizes are different\n";
    return false;
  }
  for (std::vector<int>::size_type i = 0; i < arr1->size(); i++) {
    if (prng->get_crc(arr1->at(i), chunk_size) !=
        prng->get_crc(arr2->at(i), chunk_size)) {
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

template<typename T>
void Benchmark<T>::dump(const char* name, std::vector<uint8_t*> *chunks) {
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

template<typename T>
void Benchmark<T>::dump_chunk(const char* name, uint8_t *chunk) {
  std::cout << name << ": ";
  for (int j = 0; j < chunk_size; j++) {
    std::cout << unsigned(chunk[j]) << " ";
  }
  std::cout << std::endl;
}

template<typename T>
void Benchmark<T>::reset_d_streams() {
  for (int i = 0; i < k; i++) {
    d_streams->at(i)->clear();
    d_streams->at(i)->rdbuf()->pubseekpos(0);
  }
}

template<typename T>
void Benchmark<T>::reset_c_streams() {
  for (int i = 0; i < n_c; i++) {
    c_streams->at(i)->clear();
    c_streams->at(i)->rdbuf()->pubseekpos(0);
    c_propos->at(i)->clear();
  }
}

template<typename T>
void Benchmark<T>::reset_a_streams() {
  int i;
  for (i = 0; i < n; i++) {
    a_streams->at(i)->clear();
    a_streams->at(i)->rdbuf()->pubseekpos(0);
  }
}

template<typename T>
void Benchmark<T>::reset_r_streams() {
  for (int i = 0; i < k; i++) {
    r_streams->at(i)->clear();
    r_streams->at(i)->rdbuf()->pubseekpos(0);
  }
}

template<typename T>
void Benchmark<T>::get_avail_chunks(std::vector<std::istream*> *avail_d_chunks,
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
    int avail_d_chunks_nb = 0;
    for (i = 0; i < k; i++) {
      int j = c_chunks_id->at(i);
      if (j < k) {
        avail_d_chunks->at(j) = a_streams->at(j);
        avail_d_chunks_nb++;
      } else {
        avail_c_chunks->at(j-k) = a_streams->at(j);
        avail_c_props->at(j-k) = c_propos->at(j-k);
      }
    }
    // shuffle again if all data are available
    if (avail_d_chunks_nb == k)
      get_avail_chunks(avail_d_chunks, avail_c_chunks, avail_c_props);
  } else {
    for (i = 0; i < k; i++) {
      int j = c_chunks_id->at(i);
      avail_c_chunks->at(j) = a_streams->at(j);
      avail_c_props->at(j) = c_propos->at(j);
    }
  }
}

template<typename T>
bool Benchmark<T>::encode() {
  // this operation is done per trail
  reset_d_streams();
  fec->encode_bufs(*d_streams, *c_streams, *c_propos);

  // update stats
  enc_stats->add(fec->total_enc_usec);

  // dump("d_chunks", d_chunks);
  // dump("c_chunks", c_chunks);

  return true;
}

template<typename T>
bool Benchmark<T>::decode() {
  std::vector<std::istream*> d_streams_shuffled(k, nullptr);
  std::vector<std::istream*> c_streams_shuffled(n_c, nullptr);
  std::vector<KeyValue*> c_props_shuffled(n_c, nullptr);

  get_avail_chunks(&d_streams_shuffled, &c_streams_shuffled,
    &c_props_shuffled);

  // this operation is done per trail
  reset_a_streams();
  reset_r_streams();

  if (!fec->decode_bufs(d_streams_shuffled, c_streams_shuffled,
    c_props_shuffled, *r_streams))
    return false;

  if (!compare(d_chunks, r_chunks)) {
    std::cerr << errors_desc[ERR_FAILED_REPAIR_CHUNK] << std::endl;
    return false;
  }

  // dump("r_chunks", r_chunks);

  if (!check(r_chunks)) {
    std::cerr << errors_desc[ERR_FAILED_CHUNK] << std::endl;
    return false;
  }

  // update stats
  dec_stats->add(fec->total_dec_usec);

  return true;
}

template<typename T>
bool Benchmark<T>::enc_only() {
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

template<typename T>
bool Benchmark<T>::dec_only() {
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

template<typename T>
bool Benchmark<T>::enc_dec() {
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
    "\t-g \tNumber of threads\n" +
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
        run<uint32_t>(&bench, params);
      } catch(const std::exception& e) {
        return;
      }
      break;
    }
    case 8: {
      try {
        Benchmark<uint64_t> bench(params);
        run<uint64_t>(&bench, params);
      } catch(const std::exception& e) {
        return;
      }
      break;
    }
    case 16: {
      try {
        Benchmark<__uint128_t> bench(params);
        run<__uint128_t>(&bench, params);
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

void do_join(std::thread& t)
{
    t.join();
}

int main(int argc, char **argv)
{
  PRNG prng;
  Params_t params;
  int opt;

  while ((opt = getopt(argc, argv, "t:e:w:k:m:c:n:s:x:g:")) != -1) {
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
    case 'g':
      params.threads_nb = atoi(optarg);
      break;
    default:
      xusage();
    }
  }

  params.print();

  std::vector<std::thread> threads;
  for (uint32_t thr = 0; thr < params.threads_nb; thr++) {
    threads.push_back(std::thread(run_benchmark, params));
  }
  std::for_each(threads.begin(), threads.end(), do_join);

  return 0;
}
