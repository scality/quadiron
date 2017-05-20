
#include "ntl.h"
#include "gf.cpp"
#include "gfp.cpp"
#include "gf2n.cpp"
#include "mat.cpp"
#include "vec.cpp"
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;
template class Mat<uint32_t>;
template class Vec<uint32_t>;

int vflag = 0;
int sflag = 0;

void xusage()
{
  std::cerr << "Usage: ec [-e 8|16][-n n_data][-m n_coding][-s (use cauchy instead of vandermonde)][-p prefix][-v (verbose)] -c (encode) | -r (repair) | -u (utest)\n";
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

template<typename T>
class EC
{
 public:
  GF<T> *gf;
  Mat<T> *mat;
  const char *prefix;

  EC(GF<T> *gf, int n_coding, int n_data, const char *prefix) {
    this->gf = gf;
    this->mat = new Mat<T>(gf, n_coding, n_data);
    if (sflag) {
      mat->cauchy();
    } else {
      mat->vandermonde_suitable_for_ec();
    }
    this->prefix = xstrdup(prefix);
  }

  u_int log2(u_int v)
  {
    u_int r = 0;

    if (v == 0) 
      return (u_int) -1;
    if (v == 1) 
      return 0;

    while (v > 1) {
      v >>= 1;
      r++;
    }
    return r;
  }

  int word_size()
  {
    return log2(gf->card()) / 8;
  }
  
  size_t sizew(size_t size)
  {
    return size / word_size();
  }
  
  size_t freadw(void *ptr, FILE *stream)
  {
    return fread(ptr, word_size(), 1, stream);
  }

  size_t fwritew(const void *ptr, FILE *stream)
  {
    return fwrite(ptr, word_size(), 1, stream);
  }
  
  /** 
   * (re-)create missing prefix.c1 ... cm files acc/to Vandermonde matrix
   *
   */
  void create_coding_files()
  {
    int i, j;
    char filename[1024];
    struct stat stbuf;
    size_t size = -1, wc;
    
    FILE *d_files[mat->n_cols];
    FILE *c_files[mat->n_rows];
    
    for (i = 0;i < mat->n_cols;i++) {
      snprintf(filename, sizeof (filename), "%s.d%d", prefix, i);
      if (NULL == (d_files[i] = fopen(filename, "r")))
        xerrormsg("error opening", filename);
      if (-1 == fstat(fileno(d_files[i]), &stbuf))
        xerrormsg("error stating", filename);
      if (-1 == size)
        size = stbuf.st_size;
      else if (size != stbuf.st_size)
        xmsg("bad size", filename);
    }
    
    for (i = 0;i < mat->n_rows;i++) {
      snprintf(filename, sizeof (filename), "%s.c%d", prefix, i);
      if (NULL == (c_files[i] = fopen(filename, "w")))
        xerrormsg("error opening", filename);
    }
    
    Vec<T> words = Vec<T>(mat->gf, mat->n_cols);
    Vec<T> output = Vec<T>(mat->gf, mat->n_cols);
    
    for (i = 0;i < sizew(size);i++) {
      words.zero();
      for (j = 0;j < mat->n_cols;j++) {
        wc = freadw(&VEC_ITEM(&words, j), d_files[j]);
        if (1 != wc)
          xperror("short read data");
      }
      mat->mult(&output, mat, &words);
      for (j = 0;j < mat->n_rows;j++) {
        wc = fwritew(& VEC_ITEM(&output, j), c_files[j]);
        if (1 != wc)
          xperror("short write coding");
      }
    } 
    
    for (i = 0;i < mat->n_cols;i++) {
      fclose(d_files[i]);
    }
    
    for (i = 0;i < mat->n_rows;i++) {
      fclose(c_files[i]);
    }
  }
  
  /** 
   * repair data files
   * 
   */
  int repair_data_files()
  {
    int i, j, k;
    char filename[1024];
    struct stat stbuf;
    size_t size = -1, wc;
    u_int n_data_ok = 0;
    u_int n_coding_ok = 0;
    u_int n_total_ok;
    
    FILE *d_files[mat->n_cols];
    FILE *r_files[mat->n_cols];
    FILE *c_files[mat->n_rows];
    
#define CLEANUP()                               \
    for (i = 0;i < mat->n_cols;i++) {           \
      if (NULL != d_files[i])                   \
        fclose(d_files[i]);                     \
      if (NULL != r_files[i])                   \
        fclose(r_files[i]);                     \
    }                                           \
    for (i = 0;i < mat->n_rows;i++) {           \
      if (NULL != c_files[i])                   \
        fclose(c_files[i]);                     \
    }
    
    for (i = 0;i < mat->n_cols;i++) {
      snprintf(filename, sizeof (filename), "%s.d%d", prefix, i);
      if (-1 == access(filename, F_OK)) {
        if (vflag)
          fprintf(stderr, "%s is missing\n", filename);
        d_files[i] = NULL;
        if (NULL == (r_files[i] = fopen(filename, "w")))
          xerrormsg("error opening", filename);
      } else {
        r_files[i] = NULL;
        if (NULL == (d_files[i] = fopen(filename, "r")))
          xerrormsg("error opening", filename);
        if (-1 == fstat(fileno(d_files[i]), &stbuf))
          xerrormsg("error stating", filename);
        if (-1 == size)
          size = stbuf.st_size;
        else if (size != stbuf.st_size)
          xmsg("bad size", filename);
        n_data_ok++;
      }
    }
    
    for (i = 0;i < mat->n_rows;i++) {
      snprintf(filename, sizeof (filename), "%s.c%d", prefix, i);
      if (access(filename, F_OK)) {
        if (vflag)
          fprintf(stderr, "%s is missing\n", filename);
        c_files[i] = NULL;
      } else {
        if (NULL == (c_files[i] = fopen(filename, "r")))
          xerrormsg("error opening", filename);
        n_coding_ok++;
      }
    }
    
    if (n_data_ok == mat->n_cols) {
      CLEANUP();
      return 0;
    }
    
    if (n_coding_ok < (mat->n_cols-n_data_ok)) {
      fprintf(stderr, "too many losses\n");
      CLEANUP();
      return -1;
    }
    
    n_total_ok = n_data_ok + n_coding_ok;
    
    if (vflag)
      fprintf(stderr, "n_data_ok=%d n_coding_ok=%d\n", n_data_ok, n_coding_ok);
    
    //generate a_prime
    Mat<T> a_prime(gf, n_total_ok, mat->n_cols);
    
    //for each data available generate the corresponding identity
    k = 0;
    for (i = 0;i < mat->n_cols;i++) {
      if (NULL != d_files[i]) {
        for (j = 0;j < mat->n_cols;j++) {
          if (i == j)
            MAT_ITEM(&a_prime, k, j) = 1;
          else
            MAT_ITEM(&a_prime, k, j) = 0;
        }
        k++;
      }
    }
    //finish the matrix with every coding available
    for (i = 0;i < mat->n_rows;i++) {
      if (NULL != c_files[i]) {
        //copy corresponding row in vandermonde matrix
        for (j = 0;j < mat->n_cols;j++) {
          MAT_ITEM(&a_prime, k, j) = MAT_ITEM(mat, i, j);
        }
        k++;
        //stop when we have enough codings
        if (mat->n_cols == k)
          break ;
      }
    }
    
    if (vflag) {
      fprintf(stderr, "rebuild matrix:\n");
      a_prime.dump();
    }
    
    a_prime.inv();
    
    //read-and-repair
    Vec<T> words(mat->gf, mat->n_cols);
    Vec<T> output(mat->gf, mat->n_cols);
    
    for (i = 0;i < sizew(size);i++) {
      words.zero();
      k = 0;
      for (j = 0;j < mat->n_cols;j++) {
        if (NULL != d_files[j]) {
          wc = freadw(&VEC_ITEM(&words, k), d_files[j]);
          if (1 != wc)
            xperror("short read data");
          k++;
        }
      }
      for (j = 0;j < mat->n_rows;j++) {
        if (NULL != c_files[j]) {
          wc = freadw(&VEC_ITEM(&words, k), c_files[j]);
          if (1 != wc)
            xperror("short read coding");
          k++;
          //stop when we have enough codings
          if (mat->n_cols == k)
            break ;
        }
      }
      
      mat->mult(&output, &a_prime, &words);
      
      for (j = 0;j < mat->n_cols;j++) {
        if (NULL != r_files[j]) {
          wc = fwritew(& VEC_ITEM(&output, j), r_files[j]);
          if (1 != wc)
            xperror("short write coding");
        }
      }
    } 
    
    return 0;
  }

};

template class EC<uint32_t>;

int main(int argc, char **argv)
{
  int n_data, n_coding, opt;
  char *prefix = NULL;
  int cflag = 0;
  int rflag = 0;
  int uflag = 0;
  int sflag = 0;
  int eflag = 0;

  n_data = n_coding = -1;
  prefix = NULL;
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
      n_coding = atoi(optarg);
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
  if (0 != check(n_data + n_coding)) {
    std::cerr << "Number of fragments is too big compared to Galois field size\n";
    exit(1);
  }
#endif

  if (-1 == n_data || -1 == n_coding || NULL == prefix)
    xusage();

  EC<uint32_t> *ec;

  if (eflag == 65537) {
    GFP<uint32_t> *gf = new GFP<uint32_t>(eflag);
    ec = new EC<uint32_t> (gf, n_coding, n_data, prefix);
  } else {
    GF2N<uint32_t> *gf = new GF2N<uint32_t>(eflag);
    ec = new EC<uint32_t> (gf, n_coding, n_data, prefix);
  }

  if (rflag) {
    if (0 != ec->repair_data_files()) {
      exit(1);
    }
  }
  
  ec->create_coding_files();

 end:
  free(prefix);
  return 0;
}
