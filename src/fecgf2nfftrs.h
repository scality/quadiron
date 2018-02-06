/* -*- mode: c++ -*- */
#ifndef __NTL_FECGF2NFFTRS_H__
#define __NTL_FECGF2NFFTRS_H__

#include "fec.h"
#include "fftct.h"
#include "gf2n.h"
#include "poly.h"
#include "vec.h"
#include "vvec.h"

/**
 * GF_2^n based RS using FFT transformation
 */
template<typename T>
class FECGF2NFFTRS : public FEC<T>
{
 private:
  FFTCT<T> *fft = NULL;

 public:
  T n;
  T r;
  // NOTE: only type2 is supported now
  FECGF2NFFTRS(unsigned word_size, unsigned n_data, unsigned n_parities) :
    FEC<T>(FEC<T>::TYPE_2, word_size, n_data, n_parities)
  {
    if (word_size > 16)
      assert(false);  // not support yet
    unsigned gf_n = 8*word_size;
    this->gf = new GF2N<T>(gf_n);

    // with this encoder we cannot exactly satisfy users request, we need to pad
    // n = minimal divisor of (q-1) that is at least (n_parities + n_data)
    n = this->gf->get_code_len_high_compo(n_parities + n_data);

    // compute root of order n such as r^n == 1
    this->r = this->gf->get_nth_root(n);

    // std::cerr << "n_parities=" << n_parities << "\n";
    // std::cerr << "n_data=" << n_data << "\n";
    // std::cerr << "n=" << n << "\n";
    // std::cerr << "r=" << r << "\n";

    this->fft = new FFTCT<T>(this->gf, n);
  }

  ~FECGF2NFFTRS()
  {
    delete fft;
    delete this->gf;
  }

  int get_n_outputs()
  {
    return this->n;
  }

  /**
   * Encode vector
   *
   * @param output must be n
   * @param props special values dictionary must be exactly n_data
   * @param offset used to locate special values
   * @param words must be n_data
   */
  void encode(Vec<T> *output, std::vector<KeyValue*> props, off_t offset,
    Vec<T> *words)
  {
    VVec<T> vwords(words, n);
    fft->fft(output, &vwords);
  }

  void decode_add_data(int fragment_index, int row)
  {
    // not applicable
    assert(false);
  }

  void decode_add_parities(int fragment_index, int row)
  {
    // we can't anticipate here
  }

  void decode_build()
  {
    // nothing to do
  }

  /**
   * Perform a Lagrange interpolation to find the coefficients of the polynomial
   *
   * @note If all fragments are available ifft(words) is enough
   *
   * @param output must be exactly n_data
   * @param props special values dictionary must be exactly n_data
   * @param offset used to locate special values
   * @param fragments_ids unused
   * @param words v=(v_0, v_1, ..., v_k-1) k must be exactly n_data
   */
  void decode(Vec<T> *output, std::vector<KeyValue*> props, off_t offset,
    Vec<T> *fragments_ids, Vec<T> *words)
  {
    int k = this->n_data;  // number of fragments received
    // vector x=(x_0, x_1, ..., x_k-1)
    Vec<T> vx(this->gf, k);
    for (int i = 0; i < k; i++) {
      vx.set(i, this->gf->exp(r, fragments_ids->get(i)));
    }

    // Lagrange interpolation
    Poly<T> A(this->gf), _A(this->gf);

    // compute A(x) = prod_j(x-x_j)
    A.set(0, 1);
    for (int i = 0; i < k; i++) {
      Poly<T> _t(this->gf);
      _t.set(1, 1);
      _t.set(0, vx.get(i));
      // _t.dump();
      A.mul(&_t);
    }
    // std::cout << "A(x)="; A.dump();

    // compute A'(x) since A_i(x_i) = A'_i(x_i)
    _A.copy(&A);
    _A.derivative();
    // std::cout << "A'(x)="; _A.dump();

    // evaluate n_i=v_i/A'_i(x_i)
    Vec<T> _n(this->gf, k);
    for (int i = 0; i < k; i++) {
      _n.set(i,
             this->gf->div(words->get(i),
                           _A.eval(vx.get(i))));
    }

    // compute N'(x) = sum_i{n_i * x^z_i}
    Poly<T> N_p(this->gf);
    for (int i = 0; i <= k-1; i++) {
      N_p.set(fragments_ids->get(i), _n.get(i));
    }

    // We have to find the numerator of the following expression:
    // P(x)/A(x) = sum_i=0_k-1(n_i/(x-x_i)) mod x^n
    // using Taylor series we rewrite the expression into
    // P(x)/A(x) = -sum_i=0_k-1(sum_j=0_n-1(n_i*x_i^(-j-1)*x^j))

    Poly<T> S(this->gf);
    for (T i = 0; i <= n-1; i++) {
      T val = this->gf->inv(this->gf->exp(r, i+1));
      S.set(i, N_p.eval(val));
    }
    S.mul(&A);

    // output is n_data length
    for (unsigned i = 0; i < this->n_data; i++)
      output->set(i, S.get(i));
  }
};

#endif
