
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

int vflag = 0;
int sflag = 0;

void xusage()
{
  std::cerr << "Usage: ec [-e 8|16|65537 (fermat fields)][-n n_data][-m n_coding][-s (use cauchy instead of vandermonde)][-p prefix][-v (verbose)] -c (encode) | -r (repair)\n";
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

enum tec_type
  {
    EC_TYPE_1, //take n_data input, generate n_coding outputs
    EC_TYPE_2  //take n_data input, generate n_data+n_coding outputs
  };

template<typename T>
class EC
{
protected:
  GF<T> *gf;
  tec_type ec_type;
  u_int word_size;
  char *prefix;
  u_int n_data;
  u_int n_coding;

 public:

  /** 
   * Create an encode
   * 
   * @param gf 
   * @param word_size in bytes
   * @param prefix fragment files prefix
   * @param n_data 
   * @param n_coding 
   */
  EC(GF<T> *gf, tec_type type, u_int word_size, const char *prefix, u_int n_data, u_int n_coding)
  {
    assert(ec_type == EC_TYPE_1 || ec_type == EC_TYPE_2);

    this->gf = gf;
    this->ec_type = ec_type;
    this->word_size = word_size;
    this->prefix = xstrdup(prefix);
    this->n_data = n_data;
    this->n_coding = n_coding;
  }

  ~EC()
  {
    free(prefix);
  }

public:

  virtual u_int get_n_outputs() = 0;
  virtual void encode(Vec<T> *output, Vec<T> *words) = 0;
  virtual void repair_init(void) = 0;
  virtual void repair_add_data(int k, int i) = 0;
  virtual void repair_add_coding(int k, int i) = 0;
  virtual void repair_build(void) = 0;
  virtual void repair(Vec<T> *output, Vec<T> *words) = 0;

private:

  size_t sizew(size_t size)
  {
    return size / word_size;
  }
  
  size_t freadw(T *ptr, FILE *stream)
  {
    size_t s;

    if (word_size == 1) {
      u_char c;
      if ((s = fread(&c, 1, 1, stream)) != 1)
        return s;
      *ptr = c;
      return s;
    } else if (word_size == 2) {
      u_short c;
      if ((s = fread(&c, 2, 1, stream)) != 1)
        return s;
      *ptr = c;
      return s;
    } else {
      assert(false && "no such size");
      return 0;
    }
  }

  size_t fwritew(T val, FILE *stream)
  {
    size_t s;

    if (word_size == 1) {
      u_char c = val;
      if ((s = fwrite(&c, 1, 1, stream)) != 1)
        return s;
      return s;
    } else if (word_size == 2) {
      u_short c = val;
      if ((s = fwrite(&c, 2, 1, stream)) != 1)
        return s;
      return s;
    } else {
      assert(false && "no such size");
      return 0;
    }
  }

public:

  /** 
   * (re-)create missing prefix.c1 ... cm files
   *
   */
  void create_coding_files()
  {
    int i, j;
    char filename[1024];
    struct stat stbuf;
    size_t size = -1, wc;
    FILE *d_files[n_data];
    FILE *c_files[get_n_outputs()];
    
    memset(d_files, 0, n_data * sizeof (FILE *));
    memset(c_files, 0, get_n_outputs() * sizeof (FILE *));

    for (i = 0;i < n_data;i++) {
      snprintf(filename, sizeof (filename), "%s.d%d", prefix, i);
      if (vflag)
        std::cerr << "create: opening data " << filename << "\n";
      if (NULL == (d_files[i] = fopen(filename, "r")))
        xerrormsg("error opening", filename);
      if (-1 == fstat(fileno(d_files[i]), &stbuf))
        xerrormsg("error stating", filename);
      if (-1 == size)
        size = stbuf.st_size; 
      else if (size != stbuf.st_size)
        xmsg("bad size", filename);
    }
    
    for (i = 0;i < get_n_outputs();i++) {
      snprintf(filename, sizeof (filename), "%s.c%d", prefix, i);
      if (vflag)
        std::cerr<< "create: opening coding for writing " << filename << "\n";
      if (NULL == (c_files[i] = fopen(filename, "w")))
        xerrormsg("error opening", filename);
    }
    
    Vec<T> words = Vec<T>(gf, n_data);
    Vec<T> output = Vec<T>(gf, get_n_outputs());
    
    for (i = 0;i < sizew(size);i++) {
      words.zero_fill();
      for (j = 0;j < n_data;j++) {
        T tmp;
        wc = freadw(&tmp, d_files[j]);
        if (1 != wc)
          xperror("short read data");
        words.set(j, tmp);
      }
      encode(&output, &words);
      for (j = 0;j < get_n_outputs();j++) {
        T tmp = output.get(j);
        wc = fwritew(tmp, c_files[j]);
        if (1 != wc)
          xperror("short write coding");
      }
    } 
    
    for (i = 0;i < n_data;i++) {
      fclose(d_files[i]);
    }
    
    for (i = 0;i < get_n_outputs();i++) {
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
    FILE *d_files[n_data];
    FILE *r_files[n_data];
    FILE *c_files[get_n_outputs()];

    memset(d_files, 0, n_data * sizeof (FILE *));
    memset(r_files, 0, n_data * sizeof (FILE *));
    memset(c_files, 0, get_n_outputs() * sizeof (FILE *));
    
#define CLEANUP()                               \
    for (i = 0;i < n_data;i++) {                \
      if (NULL != d_files[i])                   \
        fclose(d_files[i]);                     \
      if (NULL != r_files[i])                   \
        fclose(r_files[i]);                     \
    }                                           \
    for (i = 0;i < get_n_outputs();i++) {       \
      if (NULL != c_files[i])                   \
        fclose(c_files[i]);                     \
    }
    
    if (ec_type == EC_TYPE_1) {
      //re-read data
      for (i = 0;i < n_data;i++) {
        snprintf(filename, sizeof (filename), "%s.d%d", prefix, i);
        if (vflag)
          std::cerr << "repair: stating data " << filename << "\n";
        if (-1 == access(filename, F_OK)) {
          if (vflag)
            std::cerr << filename << " is missing\n";
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
    }
    
    for (i = 0;i < get_n_outputs();i++) {
      snprintf(filename, sizeof (filename), "%s.c%d", prefix, i);
      if (vflag)
        std::cerr << "repair: stating coding " << filename << "\n";
      if (access(filename, F_OK)) {
        if (vflag)
          std::cerr << filename << " is missing\n";
        c_files[i] = NULL;
      } else {
        if (NULL == (c_files[i] = fopen(filename, "r")))
          xerrormsg("error opening", filename);
        n_coding_ok++;
      }
    }
    
    if (n_data_ok == n_data) {
      CLEANUP();
      return 0;
    }
    
    if (n_coding_ok < (n_data - n_data_ok)) {
      std::cerr << "too many losses\n";
      CLEANUP();
      return -1;
    }
    
    if (vflag)
      std::cerr << "n_data_ok=" << n_data_ok << 
        " n_coding_ok=" << n_coding_ok << "\n";
    
    repair_init();

    k = 0;
    for (i = 0;i < n_data;i++) {
      if (NULL != d_files[i]) {
        repair_add_data(k, i);
        k++;
      }
    }

    //finish with codings available
    for (i = 0;i < get_n_outputs();i++) {
      if (NULL != c_files[i]) {
        repair_add_coding(k, i);
        k++;
        //stop when we have enough codings - XXX is it better to have more ?
        if (n_data == k)
          break ;
      }
    }

    repair_build();
    
    //read-and-repair
    Vec<T> words(gf, n_data);
    Vec<T> output(gf, n_data);
    
    for (i = 0;i < sizew(size);i++) {
      words.zero_fill();
      k = 0;
      for (j = 0;j < n_data;j++) {
        if (NULL != d_files[j]) {
          T tmp;
          wc = freadw(&tmp, d_files[j]);
          if (1 != wc)
            xperror("short read data");
          words.set(k, tmp);
          k++;
        }
      }
      for (j = 0;j < get_n_outputs();j++) {
        if (NULL != c_files[j]) {
          T tmp;
          wc = freadw(&tmp, c_files[j]);
          if (1 != wc)
            xperror("short read coding");
          words.set(k, tmp);
          k++;
          //stop when we have enough codings
          if (n_data == k)
            break ;
        }
      }

      repair(&output, &words);
      
      for (j = 0;j < n_data;j++) {
        if (NULL != r_files[j]) {
          T tmp = output.get(j);
          wc = fwritew(tmp, r_files[j]);
          if (1 != wc)
            xperror("short write coding");
        }
      }
    } 

    CLEANUP();
    return 0;
  }

};

/** 
 * GF_2^n based RS (Cauchy or Vandermonde)
 */
template<typename T>
class ECGF2NRS : public EC<T>
{
private:
  Mat<T> *mat = NULL;
  Mat<T> *repair_mat = NULL;

public:

  ECGF2NRS(GF<T> *gf, u_int word_size, const char *prefix, u_int n_data, u_int n_coding) : 
    EC<T>(gf, EC_TYPE_1, word_size, prefix, n_data, n_coding)
  {
    this->mat = new Mat<T>(gf, n_coding, n_data);
    if (sflag) {
      mat->cauchy();
    } else {
      mat->vandermonde_suitable_for_ec();
    }

    if (vflag) {
      std::cerr << "encode matrix:\n";
      mat->dump();
    }
  }

  ~ECGF2NRS()
  {
    delete mat;
    delete repair_mat;
  }

  u_int get_n_outputs()
  {
    return this->n_coding;
  }

  void encode(Vec<T> *output, Vec<T> *words)
  {
    mat->mult(output, words);
  }
  
  void repair_init(void)
  {
    //has to be a n_data*n_data invertible square matrix
    repair_mat = new Mat<T>(this->gf, mat->n_cols, mat->n_cols);
  }

  void repair_add_data(int k, int i)
  {
    //for each data available generate the corresponding identity
    for (int j = 0;j < mat->n_cols;j++) {
      if (i == j)
        repair_mat->set(k, j, 1);
      else
        repair_mat->set(k, j, 0);
    }
  }

  void repair_add_coding(int k, int i)
  {
    //copy corresponding row in vandermonde matrix
    for (int j = 0;j < mat->n_cols;j++) {
      repair_mat->set(k, j, mat->get(i, j));
    }
  }

  void repair_build()
  {
    if (vflag) {
      std::cerr << "rebuild matrix:\n";
      repair_mat->dump();
    }
    
    repair_mat->inv();
  }

  void repair(Vec<T> *output, Vec<T> *words)
  {
    repair_mat->mult(output, words);
  }
};

/** 
 * GF_2^2^k+1 based RS (Fermat Number Transform)
 * As suggested by the paper:
 * FNT-based Reed-Solomon Erasure Codes
 * by Alexandre Soro and Jerome Lacan
 */
template<typename T>
class ECFNTRS : public EC<T>
{
private:
  FFT<T> *fft = NULL;

public:
  u_int n;
  u_int N;
  u_int r;

  ECFNTRS(GF<T> *gf, u_int word_size, const char *prefix, u_int n_data, u_int n_coding) : 
    EC<T>(gf, EC_TYPE_2, word_size, prefix, n_data, n_coding)
  {
    u_int q = 65537;
    //order of a number is the lowest power of the number equals 1.
    //if q=65537 in GF_q 3 is a primitive root because is order is q-1: 
    //3^(q-1) = 1 mod q, and for each i such as 0 < i < q-1, 3^i mod q != 1
    assert(gf->_jacobi(3, q) == -1);

    //with this encoder we cannot exactly satisfy users request, we need to pad
    n = __gf._log2(n_coding + n_data) + 1;
    std::cerr << "n=" << n << "\n";
    N = __gf._exp(2, n);
    std::cerr << "N=" << N << "\n";

    //find r such as r^(N-1)=1 mod q
    bool found = false;
    for (r = 2;r < q;r++) {
      if ((__gf._exp(r, N-1) % q) == 1) {
        std::cerr << "r=" << r << "\n";
        found = true;
        break ;
      }
    } 
    assert(found);

    this->fft = new FFT<T>(gf, n, r);
  }

  ~ECFNTRS()
  {
    delete fft;
  }

  u_int get_n_outputs()
  {
    return this->N;
  }

  void encode(Vec<T> *output, Vec<T> *words)
  {
    VVec<T> vwords(words, N);
    fft->fft(output, &vwords);
  }
  
  void repair_init(void)
  {
  }

  void repair_add_data(int k, int i)
  {
    //not applicable
  }

  void repair_add_coding(int k, int i)
  {
  }

  void repair_build()
  {
  }

  void repair(Vec<T> *output, Vec<T> *words)
  {
  }
};

template class EC<uint32_t>;
template class ECGF2NRS<uint32_t>;
template class ECFNTRS<uint32_t>;

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
  //XXX TODO
  if (0 != check(n_data + n_coding)) {
    std::cerr << "Number of fragments is too big compared to Galois field size\n";
    exit(1);
  }
#endif

  if (-1 == n_data || -1 == n_coding || NULL == prefix)
    xusage();

  if (eflag == 65537) {
    ECFNTRS<uint32_t> *ec;
    GFP<uint32_t> *gf = new GFP<uint32_t>(65537); //2^2^4+1
    ec = new ECFNTRS<uint32_t>(gf, 2, prefix, n_data, n_coding);
    
    if (rflag) {
      if (0 != ec->repair_data_files()) {
        exit(1);
      }
    }
    ec->create_coding_files();
    delete ec;
  } else {
    ECGF2NRS<uint32_t> *ec;
    GF2N<uint32_t> *gf = new GF2N<uint32_t>(eflag);
    ec = new ECGF2NRS<uint32_t>(gf, eflag / 8, prefix, n_data, n_coding);
    
    if (rflag) {
      if (0 != ec->repair_data_files()) {
        exit(1);
      }
    }
    ec->create_coding_files();
    delete ec;
  }

 end:
  free(prefix);
  return 0;
}
