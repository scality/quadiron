/* -*- mode: c++ -*- */
#pragma once

/**
 * GF_2^n based RS using additive FFT transformation
 */
template<typename T>
class FECGF2NFFTADDRS : public FEC<T>
{
 private:
  FFTADD<T> *fft = NULL;
  Vec<T> *betas = NULL;

 public:
  T n;
  T m;
  // NOTE: only type2 is supported now
  FECGF2NFFTADDRS(u_int word_size, u_int n_data, u_int n_parities) :
    FEC<T>(FEC<T>::TYPE_2, word_size, n_data, n_parities)
  {
    if (word_size > 16)
      assert(false);  // not support yet
    u_int gf_n = 8*word_size;
    this->gf = new GF2N<T>(gf_n);

    // with this encoder we cannot exactly satisfy users request, we need to pad
    // n = smallest power of 2 and at least (n_parities + n_data)
    n = _get_smallest_power_of_2<T>(n_data + n_parities);
    m = _log2<T>(n);

    // std::cerr << "n_parities=" << n_parities << "\n";
    // std::cerr << "n_data=" << n_data << "\n";
    // std::cerr << "n=" << n << "\n";
    // std::cerr << "m=" << m << "\n";

    this->fft = new FFTADD<T>(this->gf, m);

    // subspace spanned by <beta_i>
    this->betas = new Vec<T>(this->gf, n);
    this->fft->compute_B(this->betas);
  }

  ~FECGF2NFFTADDRS()
  {
    delete fft;
    delete this->gf;
    if (betas) delete betas;
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
   * @note FFT(f, m, B) = ( f(B[0]), f(B[1]), .., f(B[n-1]) )
   *  where n=2^m, and B is a subspace spanned by m linearly independent
   *    beta_1, .., beta_m over GF(2).
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
    int i, j;
    int vx_zero = -1;
    int k = this->n_data;  // number of fragments received
    // vector x=(x_0, x_1, ..., x_k-1)
    Vec<T> vx(this->gf, k);
    for (i = 0; i < k; i++) {
      int _vx = this->betas->get(fragments_ids->get(i));
      vx.set(i, _vx);
      if (_vx == 0)
        vx_zero = i;
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
    for (i = 0; i < k; i++) {
      _n.set(i,
             this->gf->div(words->get(i),
                           _A.eval(vx.get(i))));
    }

    // We have to find the numerator of the following expression:
    // P(x)/A(x) = sum_i=0_k-1(n_i/(x-x_i)) mod x^n
    // using Taylor series we rewrite the expression into
    // P(x)/A(x) = -sum_i=0_k-1(sum_j=0_n-1(n_i*x_i^(-j-1)*x^j))
    Poly<T> S(this->gf);
    for (j = 0; j <= n-1; j++) {
      T val = 0;
      for (i = 0; i <= k-1; i++) {
        if (i == vx_zero)
          continue;
        // perform Taylor series at 0
        T xi_j_1 = this->gf->inv(this->gf->exp(vx.get(i), j+1));
        val = this->gf->add(val, this->gf->mul(_n.get(i), xi_j_1));
      }
      S.set(j, val);
    }
    // S.dump();
    S.mul(&A);
    if (vx_zero > -1) {
      assert(A.get(0) == 0);
      // P(x) = A(x)*S(x) + _n[vx_zero] * A(x) / x
      //  as S(x) does not include the term of vx_zero
      int deg = A.degree();
      for (i = 1; i <= deg; i++)
        S.set(i-1, this->gf->add(S.get(i-1),
                                 this->gf->mul(_n.get(vx_zero), A.get(i))));
    }
    // S.dump();

    // output is n_data length
    for (int i = 0; i < this->n_data; i++)
      output->set(i, S.get(i));
  }
};
