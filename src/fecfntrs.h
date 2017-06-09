/* -*- mode: c++ -*- */
#pragma once

/** 
 * GF_2^2^k+1 based RS (Fermat Number Transform)
 * As suggested by the paper:
 * FNT-based Reed-Solomon Erasure Codes
 * by Alexandre Soro and Jerome Lacan
 */
template<typename T>
class FECFNTRS : public FEC<T>
{
private:
  FFT<T> *fft = NULL;

public:
  u_int n;
  u_int N;
  u_int r;
  T inv_N;

  FECFNTRS(GF<T> *gf, u_int word_size, u_int n_data, u_int n_parities) : 
    FEC<T>(gf, FEC<T>::TYPE_2, word_size, n_data, n_parities)
  {
    u_int q = 65537;
    //order of a number is the lowest power of the number equals 1.
    //if q=65537 in GF_q 3 is a primitive root because is order is q-1: 
    //3^(q-1) = 1 mod q, and for each i such as 0 < i < q-1, 3^i mod q != 1
    assert(gf->_jacobi(3, q) == -1);

#if 1
    n = 3;
    N = __gf._exp(2, n);
    //r = gf->_exp(3, gf->_exp(2, 16-n)) % gf->p;
    r = 4096;
#else
    //with this encoder we cannot exactly satisfy users request, we need to pad
    n = __gf._log2(n_parities + n_data) + 1;
    N = __gf._exp(2, n);
    //find r such as r^(N-1)=1 mod q
    bool found = false;
    for (r = 2;r < q;r++) {
      if (gf->exp(r, N-1) == 1) {
        std::cerr << "r=" << r << "\n";
        found = true;
        break ;
      }
    } 
    assert(found);
#endif

    inv_N = gf->inv(N);

    std::cerr << "n=" << n << "\n";
    std::cerr << "N=" << N << "\n";
    std::cerr << "r=" << r << "\n";

    //std::cerr << gf->exp(r, N-1) << "\n";
    //assert(gf->exp(r, N-1) == 1);

    this->fft = new FFT<T>(gf, n, r);
  }

  ~FECFNTRS()
  {
    delete fft;
  }

  int get_n_fragments_required()
  {
    //XXX temp
    return this->N;
  }

  int get_n_inputs()
  {
    return this->N;
  }

  int get_n_outputs()
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

  void repair_add_data(int fragment_index, int row)
  {
    //not applicable
    assert(false);
  }

  void repair_add_parities(int fragment_index, int row)
  {
    std::cerr << "fragment_index=" << fragment_index << " row=" << row << "\n";
  }

  void repair_build()
  {
  }

  void repair(Vec<T> *output, Vec<T> *words)
  {
    VVec<T> vwords(words, N);
    fft->ifft(output, &vwords);
  }
};
