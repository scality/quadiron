
#include "ntl.h"
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;
template class Mat<uint32_t>;
template class Vec<uint32_t>;
template class FEC<uint32_t>;
template class FECGF2NRS<uint32_t>;
template class FECFNTRS<uint32_t>;

int vflag = 0;
char *prefix = NULL;

void xusage()
{
  std::cerr << "Usage: ec [-e 8|16|65537 (fermat fields)][-n n_data][-m n_parities][-s (use cauchy instead of vandermonde)][-p prefix][-v (verbose)] -c (encode) | -r (repair)\n";
  exit(1);
}

void xperror(const char *str)
{
  std::cerr << str << "\n";
  exit(1);
}

void xerrormsg(const char *str1, const char *str2)
{
  std::cerr << str1 << " " << str2 << ": " << strerror(errno) << "\n";
  exit(1);
}

void xmsg(const char *str1, const char *str2)
{
  std::cerr << str1 << " " << str2 << "\n";
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
  }

  fec->encode_bufs(d_files, c_files);

  //XXX close
#if 0
  for (int i = 0;i < fec->n_data;i++) {
    d_files[i]->close();
  }
#endif
  
  //XXX close
  for (int i = 0;i < fec->get_n_outputs();i++) {
    c_files[i]->flush();
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
  std::vector<std::ostream*> r_files(fec->n_data, nullptr);

  if (fec->fec_type == FEC<T>::TYPE_1) {
    //re-read data
    for (int i = 0;i < fec->n_data;i++) {
      snprintf(filename, sizeof (filename), "%s.d%d", prefix, i);
      if (vflag)
        std::cerr << "repair: stating data " << filename << "\n";
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
  }
  
  for (int i = 0;i < fec->get_n_outputs();i++) {
    snprintf(filename, sizeof (filename), "%s.c%d", prefix, i);
    if (vflag)
      std::cerr << "repair: stating coding " << filename << "\n";
    if (access(filename, F_OK)) {
      if (vflag)
        std::cerr << filename << " is missing\n";
      c_files[i] = nullptr;
    } else {
      c_files[i] = new std::ifstream(filename);
    }
  }

  fec->decode_bufs(d_files, c_files, r_files);

  //XXX close
#if 0
  for (int i = 0;i < fec->n_data;i++) { 
    if (nullptr != d_files[i])                         
      d_files[i]->flush();
  }
#endif

  //XXX close
  for (int i = 0;i < fec->n_data;i++) {                
    if (nullptr != r_files[i])                         
      r_files[i]->flush();
  }                       
             
  //XXX close
#if 0
  for (int i = 0;i < fec->get_n_outputs();i++) {       
    if (nullptr != c_files[i])                         
      c_files[i]->flush();
  }
#endif
  
  return 0;
}


int main(int argc, char **argv)
{
  int n_data, n_parities, opt;
  int cflag = 0;
  int rflag = 0;
  int uflag = 0;
  int sflag = 0;
  int eflag = 0;

  n_data = n_parities = -1;
  while ((opt = getopt(argc, argv, "n:m:p:scruve:")) != -1) {
    switch (opt) {
    case 'e':
      if (!strcmp(optarg, "8"))
        eflag = 8;
      else if (!strcmp(optarg, "16"))
        eflag = 16;
      else if (!strcmp(optarg, "65537"))
        eflag = 65537;
      else
        xusage();
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
    case 's':
      sflag = 1;
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

  if (eflag == 65537) {
    FECFNTRS<uint32_t> *fec;
    GFP<uint32_t> *gf = new GFP<uint32_t>(65537); //2^2^4+1
    fec = new FECFNTRS<uint32_t>(gf, 2, n_data, n_parities);
    
    if (rflag) {
      if (0 != repair_data_files<uint32_t>(fec)) {
        exit(1);
      }
    }
    create_coding_files<uint32_t>(fec);
    delete fec;
  } else {
    FECGF2NRS<uint32_t> *fec;
    GF2N<uint32_t> *gf = new GF2N<uint32_t>(eflag);
    fec = new FECGF2NRS<uint32_t>(gf, eflag / 8, n_data, n_parities, sflag ? FECGF2NRS<uint32_t>::CAUCHY : FECGF2NRS<uint32_t>::VANDERMONDE);
    
    if (rflag) {
      if (0 != repair_data_files<uint32_t>(fec)) {
        exit(1);
      }
    }
    create_coding_files<uint32_t>(fec);
    delete fec;
  }

 end:
  free(prefix);
  return 0;
}
