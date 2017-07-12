
#include "ntl.h"
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;
template class GF2N<uint64_t>;
template class Mat<uint32_t>;
template class Vec<uint32_t>;
template class FEC<uint32_t>;
template class FECGF2NRS<uint32_t>;
template class FECGF2NRS<uint64_t>;
template class FECFNTRS<uint32_t>;

int vflag = 0;
char *prefix = NULL;

void xusage()
{
  std::cerr << "Usage: ec [-e gf2nrsv|gf2nrsc|gf2nrsv-bign|gf2nrsc-bign|fntrs][-w word_size][-n n_data][-m n_parities][-p prefix][-v (verbose)] -c (encode) | -r (repair)\n";
  exit(1);
}

void xperror(const char *str)
{
  std::cerr << str << "\n";
  exit(1);
}

char *xstrdup(const char *str)
{
  char *n;

  if (NULL == (n = strdup(str)))
    xperror("malloc");
  return n;
}

/**
 * (re-)create missing prefix.c1 ... cm files
 *
 */
template <typename T>
void create_coding_files(FEC<T> *fec)
{
  char filename[1024];
  std::vector<std::istream*> d_files(fec->n_data, nullptr);
  std::vector<std::ostream*> c_files(fec->get_n_outputs(), nullptr);
  std::vector<std::ostream*> c_props_files(fec->get_n_outputs(), nullptr);
  std::vector<KeyValue*> c_props(fec->get_n_outputs(), nullptr);

  for (int i = 0;i < fec->n_data;i++) {
    snprintf(filename, sizeof (filename), "%s.d%d", prefix, i);
    if (vflag)
      std::cerr << "create: opening data " << filename << "\n";
    d_files[i] = new std::ifstream(filename);
  }

  for (int i = 0;i < fec->get_n_outputs();i++) {
    snprintf(filename, sizeof (filename), "%s.c%d", prefix, i);
    if (vflag)
      std::cerr<< "create: opening coding for writing " << filename << "\n";
    c_files[i] = new std::ofstream(filename);
    snprintf(filename, sizeof (filename), "%s.c%d.props", prefix, i);
    if (vflag)
      std::cerr<< "create: opening coding props for writing " << filename << "\n";
    c_props_files[i] = new std::ofstream(filename);
    c_props[i] = new KeyValue();
  }

  fec->encode_bufs(d_files, c_files, c_props);

  for (int i = 0;i < fec->n_data;i++) {
    (static_cast<std::ifstream*>(d_files[i]))->close();
    delete d_files[i];
  }

  for (int i = 0;i < fec->get_n_outputs();i++) {

    *(c_props_files[i]) << *(c_props[i]);

    (static_cast<std::ofstream*>(c_props_files[i]))->close();
    delete c_props_files[i];
    delete c_props[i];

    (static_cast<std::ofstream*>(c_files[i]))->close();
    delete c_files[i];
  }
}

/**
 * repair data files
 *
 */
template <typename T>
bool repair_data_files(FEC<T> *fec)
{
  char filename[1024];
  std::vector<std::istream*> d_files(fec->n_data, nullptr);
  std::vector<std::istream*> c_files(fec->get_n_outputs(), nullptr);
  std::vector<std::istream*> c_props_files(fec->get_n_outputs(), nullptr);
  std::vector<KeyValue*> c_props(fec->get_n_outputs(), nullptr);
  std::vector<std::ostream*> r_files(fec->n_data, nullptr);

  //re-read data
  for (int i = 0;i < fec->n_data;i++) {
    snprintf(filename, sizeof (filename), "%s.d%d", prefix, i);
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

  for (int i = 0;i < fec->get_n_outputs();i++) {
    snprintf(filename, sizeof (filename), "%s.c%d", prefix, i);
    if (vflag)
      std::cerr << "repair: checking coding " << filename << "\n";
    if (access(filename, F_OK)) {
      if (vflag)
        std::cerr << filename << " is missing\n";
      c_files[i] = nullptr;
    } else {
      c_files[i] = new std::ifstream(filename);
    }

    snprintf(filename, sizeof (filename), "%s.c%d.props", prefix, i);
    if (vflag)
      std::cerr << "repair: checking coding props " << filename << "\n";
    if (access(filename, F_OK)) {
      c_props_files[i] = nullptr;
      c_props[i] = nullptr;
    } else {
      c_props_files[i] = new std::ifstream(filename);
      c_props[i] = new KeyValue();
      *(c_props_files[i]) >> *(c_props[i]);
    }
  }

  fec->decode_bufs(d_files, c_files, c_props, r_files);

  for (int i = 0;i < fec->n_data;i++) {
    if (nullptr != d_files[i]) {
      (static_cast<std::ifstream*>(d_files[i]))->close();
      delete d_files[i];
    }
  }

  for (int i = 0;i < fec->get_n_outputs();i++) {
    if (nullptr != c_props_files[i]) {
      (static_cast<std::ifstream*>(c_props_files[i]))->close();
      delete c_props_files[i];
    }

    if (nullptr != c_props[i])
      delete c_props[i];

    if (nullptr != c_files[i]) {
      (static_cast<std::ifstream*>(c_files[i]))->close();
      delete c_files[i];
    }
  }

  for (int i = 0;i < fec->n_data;i++) {
    if (nullptr != r_files[i]) {
      (static_cast<std::ofstream*>(r_files[i]))->close();
      delete r_files[i];
    }
  }

  return 0;
}

template <typename T>
void print_stats(FEC<T> *fec)
{
  std::cerr << "enc," << (fec->n_encode_ops != 0 ? fec->total_encode_cycles / fec->n_encode_ops : 0) << ",";
  std::cerr << "dec," << (fec->n_decode_ops != 0 ? fec->total_decode_cycles / fec->n_decode_ops : 0) << ",";
}

enum ec_type
  {
    EC_TYPE_UNDEF = 0,
    EC_TYPE_GF2NRS,
    EC_TYPE_GF2NRS_BIGN,
    EC_TYPE_FNTRS,
  };

int main(int argc, char **argv)
{
  int n_data, n_parities, opt;
  int cflag = 0;
  int rflag = 0;
  int uflag = 0;
  ec_type eflag = EC_TYPE_UNDEF;
  FECGF2NRS<uint32_t>::FECGF2NRSType gf2nrs_type;
  FECGF2NRS<uint64_t>::FECGF2NRSType gf2nrs_bign_type;
  int word_size = 0;

  n_data = n_parities = -1;
  while ((opt = getopt(argc, argv, "n:m:p:cruve:w:")) != -1) {
    switch (opt) {
    case 'e':
      if (!strcmp(optarg, "gf2nrsv")) {
        eflag = EC_TYPE_GF2NRS;
        gf2nrs_type = FECGF2NRS<uint32_t>::VANDERMONDE;
      } else if (!strcmp(optarg, "gf2nrsc")) {
        eflag = EC_TYPE_GF2NRS;
        gf2nrs_type = FECGF2NRS<uint32_t>::CAUCHY;
      } else if (!strcmp(optarg, "gf2nrsv-bign")) {
        eflag = EC_TYPE_GF2NRS_BIGN;
        gf2nrs_bign_type = FECGF2NRS<uint64_t>::VANDERMONDE;
      } else if (!strcmp(optarg, "gf2nrsc-bign")) {
        eflag = EC_TYPE_GF2NRS_BIGN;
        gf2nrs_bign_type = FECGF2NRS<uint64_t>::CAUCHY;
      } else if (!strcmp(optarg, "fntrs"))
        eflag = EC_TYPE_FNTRS;
      else
        xusage();
      break ;
    case 'w':
      word_size = atoi(optarg);
      break ;
    case 'v':
      vflag = 1;
      break ;
    case 'u':
      uflag = 1;
      break ;
    case 'c':
      cflag = 1;
      break ;
    case 'r':
      rflag = 1;
      break ;
    case 'n':
      n_data = atoi(optarg);
      break;
    case 'm':
      n_parities = atoi(optarg);
      break;
    case 'p':
      prefix = xstrdup(optarg);
      break;
    default: /* '?' */
      xusage();
    }
  }

  if (!eflag)
    xusage();

  if (!(uflag || cflag || rflag))
    xusage();

#if 0
  //XXX TODO
  if (0 != check(n_data + n_parities)) {
    std::cerr << "Number of fragments is too big compared to Galois field size\n";
    exit(1);
  }
#endif

  if (-1 == n_data || -1 == n_parities || NULL == prefix)
    xusage();

  if (eflag == EC_TYPE_FNTRS) {
    FECFNTRS<uint32_t> *fec;
    if (word_size != 2) {
      std::cerr << "only supports -w 2 for now\n";
      exit(1);
    }
    GFP<uint32_t> *gf = new GFP<uint32_t>(__gf32._exp(2, word_size * 8)+1); //2^2^4+1
    fec = new FECFNTRS<uint32_t>(gf, word_size, n_data, n_parities);

    if (rflag) {
      if (0 != repair_data_files<uint32_t>(fec)) {
        exit(1);
      }
    }
    create_coding_files<uint32_t>(fec);
    print_stats<uint32_t>(fec);
    delete fec;
    delete gf;
  } else if (eflag == EC_TYPE_GF2NRS) {
    FECGF2NRS<uint32_t> *fec;
    GF2N<uint32_t> *gf = new GF2N<uint32_t>(word_size * 8);
    fec = new FECGF2NRS<uint32_t>(gf, word_size, n_data, n_parities, gf2nrs_type);

    if (rflag) {
      if (0 != repair_data_files<uint32_t>(fec)) {
        exit(1);
      }
    }
    create_coding_files<uint32_t>(fec);
    print_stats<uint32_t>(fec);
    delete fec;
    delete gf;
  } else if (eflag == EC_TYPE_GF2NRS_BIGN) {
    FECGF2NRS<uint64_t> *fec;
    GF2N<uint64_t> *gf = new GF2N<uint64_t>(word_size * 8);
    fec = new FECGF2NRS<uint64_t>(gf, word_size, n_data, n_parities,
      gf2nrs_bign_type);

    if (rflag) {
      if (0 != repair_data_files<uint64_t>(fec)) {
        exit(1);
      }
    }
    create_coding_files<uint64_t>(fec);
    print_stats<uint64_t>(fec);
    delete fec;
    delete gf;
  }

 end:
  free(prefix);
  return 0;
}
