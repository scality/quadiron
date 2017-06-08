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

  FECFNTRS(GF<T> *gf, u_int word_size, u_int n_data, u_int n_parities) : 
    FEC<T>(gf, FEC<T>::TYPE_2, word_size, n_data, n_parities)
  {
    u_int q = 65537;
    //order of a number is the lowest power of the number equals 1.
    //if q=65537 in GF_q 3 is a primitive root because is order is q-1: 
    //3^(q-1) = 1 mod q, and for each i such as 0 < i < q-1, 3^i mod q != 1
    assert(gf->_jacobi(3, q) == -1);

    //with this encoder we cannot exactly satisfy users request, we need to pad
    n = __gf._log2(n_parities + n_data) + 1;
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

  ~FECFNTRS()
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

  void repair_add_parities(int k, int i)
  {
  }

  void repair_build()
  {
  }

  void repair(Vec<T> *output, Vec<T> *words)
  {
  }
};
