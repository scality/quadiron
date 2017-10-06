
/* -*- mode: c++ -*- */

#pragma once

/**
 * Additive FFT algorithm of length n = 2^m (arbitrary m)
 *  Algorithm 2 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 * @param output
 * @param input
 */
template<typename T>
class FFTADD : public FFT<T>
{
 private:
  bool create_betas;
  T m;
  T k, deg0, deg1, deg2;
  T beta_1, inv_beta_1;
  T beta_m, inv_beta_m;
  Vec<T> *betas = NULL;
  Vec<T> *gammas = NULL;
  Vec<T> *deltas = NULL;
  Vec<T> *beta_m_powers = NULL;
  Vec<T> *G = NULL;
  Vec<T> *g0 = NULL;
  Vec<T> *g1 = NULL;
  Vec<T> *u = NULL;
  Vec<T> *v = NULL;
  Vec<T> *mem = NULL;
  FFTADD<T> *fft_add = NULL;
 public:
  FFTADD(GF<T> *gf, T m, Vec<T>* betas = NULL);
  ~FFTADD();
  void compute_basis();
  void compute_beta_m_powers();
  void compute_G();
  void compute_B(Vec<T> *B);
  void compute_subspace(Vec<T> *basis, Vec<T> *subspace);
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
  void fft_inv(Vec<T> *output, Vec<T> *input);
  void taylor_expand(Vec<T> *output, Vec<T> *input, int n, int t,
    bool do_copy=false);
  void taylor_expand_t2(Vec<T> *g0, Vec<T> *g1, Vec<T> *input, int n,
    bool do_copy=false);
  void inv_taylor_expand(Vec<T> *output, Vec<T> *input, int t);
  void inv_taylor_expand_t2(Vec<T> *output, Vec<T> *g0, Vec<T> *g1);
 private:
  int find_k(int n, int t);
  void compute_g0_g1(Vec<T> *vec, int deg0, int deg2, Vec<T> *g0, Vec<T> *g1);
  void _taylor_expand(Vec<T> *output, Vec<T> *input, int n, int t);
  void _taylor_expand_t2(Vec<T> *g0, Vec<T> *g1, Vec<T> *input, int n, int k,
    int start);
  void mul_xt_x(Vec<T> *vec, int t);
  void _fft(Vec<T> *output, Vec<T> *input);
  void _ifft(Vec<T> *output, Vec<T> *input);
};

/**
 * @param gf
 * @param n
 * @param id: index in the list of factors of n
 *
 * @return
 */
template <typename T>
FFTADD<T>::FFTADD(GF<T> *gf, T m, Vec<T>* betas) :
  FFT<T>(gf, gf->arith->exp2(m))
{
  assert(m >= 1);
  this->m = m;
  create_betas = false;
  if (betas == NULL) {
    // it supports only GF2N
    assert(gf->get_sub_field()->card() == 2);

    create_betas = true;
    betas = new Vec<T>(gf, m);
    T beta = this->gf->get_prime_root();
    betas->set(0, 1);
    betas->set(1, beta);
    for (T i = 2; i < m; i++) {
      betas->set(i, this->gf->exp(beta, i));
    }
}
  this->betas = betas;
  this->beta_1 = betas->get(0);
  this->inv_beta_1 = gf->inv(this->beta_1);
  this->beta_m = betas->get(m-1);
  this->inv_beta_m = gf->inv(this->beta_m);

  this->beta_m_powers = new Vec<T>(gf, this->n);
  this->compute_beta_m_powers();
  if (m > 1) {
    this->k = this->gf->arith->exp2(m - 1);

    // compute gammas and deltas
    this->gammas = new Vec<T>(gf, m-1);
    this->deltas = new Vec<T>(gf, m-1);
    this->G = new Vec<T>(gf, this->k);
    this->compute_basis();

    this->u = new Vec<T>(gf, this->k);
    this->v = new Vec<T>(gf, this->k);
    this->g0 = new Vec<T>(gf, this->k);
    this->g1 = new Vec<T>(gf, this->k);

    this->mem = new Vec<T>(gf, this->n);
    this->fft_add = new FFTADD(gf, m-1, this->deltas);
  }
}

template <typename T>
FFTADD<T>::~FFTADD()
{
  if (create_betas && betas) delete betas;
  if (gammas) delete gammas;
  if (deltas) delete deltas;
  if (beta_m_powers) delete beta_m_powers;
  if (G) delete G;
  if (u) delete u;
  if (v) delete v;
  if (g0) delete g0;
  if (g1) delete g1;
  if (mem) delete mem;
  if (fft_add) delete fft_add;
}

template <typename T>
void FFTADD<T>::compute_beta_m_powers()
{
  int i;
  // compute beta_m^i
  T beta = beta_m;
  this->beta_m_powers->set(0, 1);
  this->beta_m_powers->set(1, beta);
  for (i = 2; i < this->n; i++) {
    beta = this->gf->mul(beta, beta_m);
    this->beta_m_powers->set(i, beta);
  }
}

template <typename T>
void FFTADD<T>::compute_basis()
{
  int i;
  for (i = 0; i < m-1; i++) {
    T gamma = this->gf->mul(inv_beta_m, betas->get(i));
    this->gammas->set(i, gamma);
    this->deltas->set(i, this->gf->add(this->gf->exp(gamma, 2), gamma));
  }
  compute_G();
}

/*
 * Compute subspace spanned by gammas
 */
template <typename T>
void FFTADD<T>::compute_G()
{
  compute_subspace(gammas, G);
}

/*
 * Compute subspace spanned by betas
 */
template <typename T>
void FFTADD<T>::compute_B(Vec<T> *B)
{
  compute_subspace(betas, B);
}

/*
 * Compute subspace spanned by basis
 *  We use balanced Gray sequence to decrease computations number
 *  Intuitively, binary image of two consecutive number of the Gray sequence
 *    differ at 1 position
 *
 * @param basis: vector of `dim` linear independent elements
 * @param subspace: output, a vector of 2^dim length
 */
template <typename T>
void FFTADD<T>::compute_subspace(Vec<T> *basis, Vec<T> *subspace)
{
  int i, j;
  int size = 1;
  int dim = basis->get_n();
  T val = 0;
  assert(subspace->get_n() == this->gf->arith->exp2(dim));
  // index of element in G
  std::vector<int> ids;
  // create sequence 2^i
  Vec<T> powers2(this->gf, dim);
  for (i = 0; i < dim; i++)
    powers2.set(i, this->gf->arith->exp2(i));

  ids.push_back(0);
  subspace->set(0, 0);
  for (i = 0; i < dim; i++) {
    for (j = size - 1; j >= 0; j--) {
      int id = this->gf->add(ids.at(j), powers2.get(i));
      int diff = this->gf->arith->log2(this->gf->add(id, ids.back()));
      val = this->gf->add(val, basis->get(diff));
      subspace->set(id, val);
      // for next i
      ids.push_back(id);
    }
    // for next i;
    size *= 2;
  }
}

template <typename T>
void FFTADD<T>::_fft(Vec<T> *output, Vec<T> *input)
{
  int i;
  int deg = input->get_n();
  mem->set(0, input->get(0));
  for (i = 1; i < deg; i++) {
    mem->set(i, this->gf->mul(beta_m_powers->get(i), input->get(i)));
  }

  // compute taylor expansion of g(x) at (x^2 - x)
  taylor_expand_t2(g0, g1, mem, this->n);

  this->fft_add->fft(u, g0);
  this->fft_add->fft(v, g1);

  for (i = 0; i < k; i++) {
    T w_i = this->gf->add(u->get(i), this->gf->mul(G->get(i), v->get(i)));
    output->set(i, w_i);
    output->set(k + i, this->gf->add(w_i, v->get(i)));
  }
}

template <typename T>
void FFTADD<T>::_ifft(Vec<T> *output, Vec<T> *input)
{
  int i;
  output->zero_fill();
  for (i = 0; i < k; i++) {
    v->set(i, this->gf->add(input->get(k+i), input->get(i)));
    u->set(i, this->gf->add(input->get(i),
                            this->gf->mul(G->get(i), v->get(i))));
  }

  this->fft_add->ifft(g0, u);
  this->fft_add->ifft(g1, v);

  inv_taylor_expand_t2(output, g0, g1);

  T beta = inv_beta_m;
  T deg = output->get_n();
  for (i = 1; i < deg; i++) {
    output->set(i, this->gf->mul(beta, output->get(i)));
    beta = this->gf->mul(beta, inv_beta_m);
  }
}

/**
 * Additive FFT
 *  Algorithm 2 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 * Input: (f, m, B) where m >= 1, f(x) of degree < n = 2
 *    B = <beta_i> where beta_i are linearly independent over F2
 * Output: FFT(f, m, B) = ( f(B[0]), f(B[1]), .., f(B[n-1]) )
 *
 * Step1: if m = 1, return ( f(0), f(1) )
 * Step2: Compute g(x) = f(beta_m * x)
 * Step3: Compute Taylor expansion of g(x) over (x^2 - x), i.e. with k = 2^(m-1)
 *         g(x) = sum_{i=0, k-1}(g_i0 + g_i1 * x) * (x^2 - x)^i
 *     Define g0(x), g1(x) as below:
 *         g0(x) = sum_{i=0, k-1} g_i0 * x^i
 *         g1(x) = sum_{i=0, k-1} g_i1 * x^i
 * Step4: For 1 <= i <= m-1, compute
 *                 gamma_i = beta_i * beta_m^(-1)
 *                 delta_i = gamma_i^2 - gamma_i
 *     Define G = <gamma_i>, D = <delta_i>
 * Step5: Compute
 *         FFT(g0, m-1, D) = ( u_0, u_1, .., u_{k-1} )
 *         FFT(g1, m-1, D) = ( v_0, v_1, .., v_{k-1} )
 * Step6: for i = 0, .., k, set
 *         w_i = u_i + G[i] * v_i
 *         w_{k+1} = w_i + v_i
 * Step7: return (w_0, w_1, .., w_{n-1})
 *
 * @param output: output of hi(x) polynomials
 * @param input: input polynomial f(x)
 */
template <typename T>
void FFTADD<T>::fft(Vec<T> *output, Vec<T> *input)
{
  if (m > 1)
    return _fft(output, input);

  // m == 1 -> output = (f(0), f(beta_1))
  output->set(0, input->get(0));
  output->set(1, this->gf->add(input->get(0),
                               this->gf->mul(input->get(1), this->beta_1)));
}

template <typename T>
void FFTADD<T>::fft_inv(Vec<T> *output, Vec<T> *input)
{
  if (m > 1)
    return _ifft(output, input);

  // m == 1 -> return ( input[0], (input[1] - input[0])*beta_1^-1 )
  //  because input[0] = f(0) = f0,
  //          input[1] = f(beta_1) = f0 + beta_1 * f1
  output->set(0, input->get(0));
  output->set(1, this->gf->mul(this->inv_beta_1,
                               this->gf->add(input->get(0), input->get(1))));
}

/**
 * Additive inverse FFT
 *  Algorithm 2 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 * Input: (W, m, B) where m >= 1, where W = (w_0, w_1, .., w_{n-1})
 *    B = <beta_i> where beta_i are linearly independent over F2
 * Output: IFFT(W, m, B) = f(x)
 *
 * Step1: if m = 1, return ( w_0, (w_1 - w_0)*beta_1^-1 )
 *     since w_0 = f(0), w_1 = f(1) = f0 + beta_1 * f1 where f(x) = f0 + f1 * x
 * Step2: for i = 0, .., k, set
 *         v_i = w_{k+i} - w_i
 *         u_i = w_i - G[i] * v_i
 * Step3: If D and G are note computed, compute gamma_i = beta_i * beta_m^(-1)
 *                delta_i = gamma_i^2 - gamma_i
 *     Define G = <gamma_i>, D = <delta_i>
 * Step4: Compute
 *         iFFT(u, m-1, D) = g0, where u = ( u_0, u_1, .., u_{k-1} )
 *         iFFT(v, m-1, D) = g1, where v = ( v_0, v_1, .., v_{k-1} )
 * Step5: Compute
 *         g(x) = sum_{i=0, k-1}(g_i0 + g_i1 * x) * (x^2 - x)^i
 *     where g0(x), g1(x) are computed as above, then:
 *         g0(x) = sum_{i=0, k-1} g_i0 * x^i
 *         g1(x) = sum_{i=0, k-1} g_i1 * x^i
 * Step6: Compute f(x) = g(beta_m^(-1) * x)
 *
 * @param output: output of hi(x) polynomials
 * @param input: input polynomial f(x)
 */
template <typename T>
void FFTADD<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  fft_inv(output, input);
}

/**
 * Taylor expansion at (x^t - x)
 *  Algorithm 1 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 *  f(x) = sum_i h_i(x) * (x^t - x)^i
 *
 * @param output: output of hi(x) polynomials
 * @param input: input polynomial f(x)
 * @param n
 * @param t
 * @param do_copy: flag to do a copy of input or not since this function will
 *  modify input vector
 */
template <typename T>
void FFTADD<T>::taylor_expand(Vec<T> *output, Vec<T> *input, int n, int t,
  bool do_copy)
{
  assert(n >= 1);
  assert(t > 1);
  assert(input->get_n() <= n);
  output->zero_fill();

  if (!do_copy) {
    _taylor_expand(output, input, n, t);
    return;
  }

  Vec<T> *_input = new Vec<T>(this->gf, input->get_n());
  _input->copy(input);
  _taylor_expand(output, _input, n, t);

  delete _input;
}

template <typename T>
void FFTADD<T>::_taylor_expand(Vec<T> *output, Vec<T> *input, int n, int t)
{
  if (n <= t) {
    // TODO: copy_mutual has not worked yet for VmVec due to copying of `mem`
    output->add_mutual(input);
    return;
  }
  // find k s.t. t2^k < n <= 2 *t2^k
  int k = find_k(n, t);
  int deg2 = this->gf->arith->exp2(k);
  int deg0 = t*deg2;

  VmVec<T> g0(input, deg0);
  VmVec<T> g1(input, n - deg0, deg0);
  VmVec<T> V1(output, deg0);
  VmVec<T> V2(output, n - deg0, deg0);

  compute_g0_g1(input, deg0, deg2, &g0, &g1);
  _taylor_expand(&V1, &g0, deg0, t);
  _taylor_expand(&V2, &g1, n - deg0, t);
}


/**
 * Taylor expansion at (x^2 - x)
 *  Algorithm 1 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 * This function is used for FFT over n=2^m, hence n must be power of 2
 *
 * @param G0: vector for poly of gi0 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most t*2^k
 * @param G1: vector poly of gi1 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most (n-t*2^k)
 * @param n: a power of 2
 * @param do_copy: flag to do a copy of input or not since this function will
 *  modify input vector
 */
template <typename T>
void FFTADD<T>::taylor_expand_t2(Vec<T> *G0, Vec<T> *G1, Vec<T> *input, int n,
  bool do_copy)
{
  assert(n >= 1);
  assert(input->get_n() <= n);

  // find k s.t. t2^k < n <= 2 *t2^k
  int _k = find_k(n, 2);

  if (!do_copy) {
    _taylor_expand_t2(G0, G1, input, n, _k, 0);
    return;
  }

  Vec<T> *_input = new Vec<T>(this->gf, input->get_n());
  _input->copy(input);
  _taylor_expand_t2(G0, G1, _input, n, _k, 0);

  delete _input;
}

/**
 * Taylor expansion at (x^2 - x)
 *  Algorithm 1 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 * This function is used for FFT over n=2^m, hence n must be power of 2
 *
 * @param G0: vector for poly of gi0 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most t*2^k
 * @param G1: vector poly of gi1 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most (n-t*2^k)
 *
 * @param n: a power of 2
 * @param k: a number st t*2^k < n <= 2*t*2^k
 * @param s_deg: start degree of G0 and G1
 */
template <typename T>
void FFTADD<T>::_taylor_expand_t2(Vec<T> *G0, Vec<T> *G1, Vec<T> *input, int n,
  int _k, int s_deg)
{
  if (n <= 2) {
    G0->set(s_deg, input->get(0));
    if (n == 2)
      G1->set(s_deg, input->get(1));
    else
      G1->set(s_deg, 0);
    return;
  }

  int deg2 = this->gf->arith->exp2(_k);
  int deg0 = 2*deg2;

  VmVec<T> g0(input, deg0);
  VmVec<T> g1(input, n - deg0, deg0);
  compute_g0_g1(input, deg0, deg2, &g0, &g1);

  _taylor_expand_t2(G0, G1, &g0, deg0, _k-1, s_deg);
  _taylor_expand_t2(G0, G1, &g1, n - deg0, _k-1, s_deg + deg2);
}

template <typename T>
void FFTADD<T>::compute_g0_g1(Vec<T> *vec, int deg0, int deg2,
  Vec<T> *g0, Vec<T> *g1) {
  // find k s.t. t*2^k < n <= 2*t*2^k
  int deg1 = deg0 - deg2;

  // find f0, f1, f2 s.t.
  // f = f0 + x^(t2^k)( f1 + x^( (t-1)2^k) f2 )
  // NOTE: f0 === g0
  VmVec<T> f1(vec, deg1, deg0, 1);
  VmVec<T> f2(vec, deg2, deg0+deg1, 1);

  // since deg1 >= deg2, add f2 into f1 to perform h === f1
  // that is also actually g1
  f1.add_mutual(&f2);
  // add h (=== f1) into g0 with offset deg2=2^k
  g0->add_mutual(&f1, deg2);
  // since g1 = (f1+f2, f2) so nothing to do
}

/*
 * This function compute f(x) from its taylor expansion
 *
 *  f(x) = sum_i res_i(x) * (x^t - x)^i
 * We perform it in the following way
 *  f(x) = (..(( h_k*y + h_{k-1} )*y + h_{k-2})*y ...)
 *  where y = x^t - x
 *        h_i := ith polynomial of taylor expansion output
 *
 * @param result: vector of hi(x) polynomials
 * @param t
 */
template <typename T>
void FFTADD<T>::inv_taylor_expand(Vec<T> *output, Vec<T> *input, int t) {
  output->zero_fill();

  int i = input->get_n() - t;
  VmVec<T> hi(input, t, i);
  output->add_mutual(&hi);
  while(i >= t) {
    i -= t;
    // multiply output to (x^t - x)
    mul_xt_x(output, t);
    hi.set_map(i, 1);
    output->add_mutual(&hi);
  }
}

/*
 * This function compute f(x) from its taylor expansion
 *  f(x) = sum_i (gi0 + gi1 * x) * (x^2 - x)^i
 *
 * We perform it in the following way
 *  f(x) = (..(( h_k*y + h_{k-1} )*y + h_{k-2})*y ...)
 *  where y = x^2 - x
 *        h_i = gi0 + gi1*x
 *
 * @param G0: vector for poly of gi0 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most t*2^k
 * @param G1: vector poly of gi1 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most (n-t*2^k)
 */
template <typename T>
void FFTADD<T>::inv_taylor_expand_t2(Vec<T> *output, Vec<T> *G0, Vec<T> *G1) {
  assert(G0->get_n() == G1->get_n());
  output->zero_fill();

  int i = G0->get_n() - 1;
  output->set(0, G0->get(i));
  output->set(1, G1->get(i));
  while(--i >= 0) {
    // multiply output to (x^2 - x)
    mul_xt_x(output, 2);
    output->set(0, G0->get(i));
    if (G1->get(i) > 0)
      output->set(1, this->gf->add(output->get(1), G1->get(i)));
  }
}

// multiply vec as f(x) to (x^t + x)
template <typename T>
void FFTADD<T>::mul_xt_x(Vec<T> *vec, int t) {
  int i;
  // f_i = f_{i-1} + f_{i-t} for i >= t
  for (i = vec->get_n() - 1; i >= t; i--) {
    vec->set(i, this->gf->add(vec->get(i-t), vec->get(i-1)));
  }
  // f_i = f_{i-1} for i >= 1
  for (; i >= 1; i--) {
    vec->set(i, vec->get(i-1));
  }
  // f_0 = 0
  vec->set(0, 0);
}

/**
 * find k s.t. t2^k < n <= 2 *t2^k
 *
 * @param n: n
 * @param t: t
 * @return k: such that t2^k < n <= 2 *t2^k
 */
template <typename T>
inline int FFTADD<T>::find_k(int n, int t)
{
  int k;
  int t2k = t;  // init for t*2^k with k = 0
  for (k = 0; k < n; k++) {
    if (t2k < n && 2*t2k >= n)
      return k;
    // next t2k
    t2k *= 2;
  }
  return 0;
}
