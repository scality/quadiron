/* -*- mode: c++ -*- */
#ifndef __NTL_FFTADD_H__
#define __NTL_FFTADD_H__

#include "dft.h"
#include "vec.h"
#include "vmvec.h"

/**
 * Additive FFT algorithm of length n = 2^m (arbitrary m)
 *  Algorithm 2 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 * @param output
 * @param input
 */
template<typename T>
class FFTADD : public DFT<T>
{
 private:
  bool create_betas;
  T m;
  T k, deg0, deg1, deg2;
  T beta_1, inv_beta_1;
  T beta_m, inv_beta_m;
  Vec<T> *betas = nullptr;
  Vec<T> *gammas = nullptr;
  Vec<T> *deltas = nullptr;
  Vec<T> *beta_m_powers = nullptr;
  Vec<T> *G = nullptr;
  Vec<T> *g0 = nullptr;
  Vec<T> *g1 = nullptr;
  Vec<T> *u = nullptr;
  Vec<T> *v = nullptr;
  Vec<T> *mem = nullptr;
  FFTADD<T> *fft_add = nullptr;
 public:
  FFTADD(GF<T> *gf, T m, Vec<T>* betas = nullptr);
  ~FFTADD();
  void compute_basis();
  void compute_beta_m_powers();
  void compute_G();
  void compute_B(Vec<T> *B);
  void compute_subspace(Vec<T> *basis, Vec<T> *subspace);
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
  void fft_inv(Vec<T> *output, Vec<T> *input);
  void taylor_expand_t2(Vec<T> *input, int n, bool do_copy=false);
  void taylor_expand(Vec<T> *output, Vec<T> *input, int n, int t);
  void inv_taylor_expand_t2(Vec<T> *output);
  void inv_taylor_expand(Vec<T> *output, Vec<T> *input, int t);
 private:
  int find_k(int n, int t);
  void _taylor_expand_t2(Vec<T> *input, int n, int k, int start);
  void _taylor_expand(Vec<T> *input, int n, int t);
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
  DFT<T>(gf, _exp2<T>(m))
{
  assert(m >= 1);
  this->m = m;
  create_betas = false;
  if (betas == nullptr) {
    // it supports only GF2N
    assert(gf->get_sub_field()->card() == 2);

    create_betas = true;
    betas = new Vec<T>(gf, m);
    T beta = this->gf->get_prime_root();
    betas->set(0, beta);
    for (T i = 1; i < m-1; i++) {
      betas->set(i, this->gf->exp(beta, i+1));
    }
    betas->set(m-1, 1);
  }
  this->betas = betas;
  this->beta_1 = betas->get(0);
  this->inv_beta_1 = gf->inv(this->beta_1);
  this->beta_m = betas->get(m-1);
  this->inv_beta_m = gf->inv(this->beta_m);

  this->beta_m_powers = new Vec<T>(gf, this->n);
  this->compute_beta_m_powers();
  if (m > 1) {
    this->k = _exp2<T>(m - 1);

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
  for (T i = 0; i < m-1; i++) {
    T gamma = betas->get(i);
    if (beta_m > 1)
      gamma = this->gf->mul(inv_beta_m, betas->get(i));
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
  assert(subspace->get_n() == _exp2<T>(dim));
  // index of element in G
  std::vector<int> ids;
  // create sequence 2^i
  Vec<T> powers2(this->gf, dim);
  for (i = 0; i < dim; i++)
    powers2.set(i, _exp2<T>(i));

  ids.push_back(0);
  subspace->set(0, 0);
  for (i = 0; i < dim; i++) {
    for (j = size - 1; j >= 0; j--) {
      int id = this->gf->add(ids.at(j), powers2.get(i));
      int diff = _log2<T>(this->gf->add(id, ids.back()));
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
  mem->copy(input, this->n);
  if (beta_m > 1)
    mem->hadamard_mul(beta_m_powers);

  // compute taylor expansion of g(x) at (x^2 - x)
  // outputs are this->g0 and this->g1
  taylor_expand_t2(mem, this->n);

  this->fft_add->fft(u, g0);
  this->fft_add->fft(v, g1);

  // copy output = (undefined, v)
  output->copy(v, k, k);
  // perform G * v hadamard multiplication
  v->hadamard_mul(G);
  // v += u
  v->add(u);
  // output = (u + G*v, v)
  output->copy(v, k);
  // output = (u + G*v, (u + G*v) + v)
  output->add(v, k);
}

template <typename T>
void FFTADD<T>::_ifft(Vec<T> *output, Vec<T> *input)
{
  output->zero_fill();
  /*
   * input = (w0, w1)
   * calculate u, v s.t. v = w1 - w0, u = w0 - G * v
   */
  VmVec<T> w0(input, k);
  // v = w_1
  v->copy(input, k, 0, k);
  // v = w0 + w1
  v->add(&w0);
  // u = v
  u->copy(v);
  // u = G * v
  u->hadamard_mul(G);
  // u = w0 + G * v;
  u->add(&w0);

  this->fft_add->ifft(g0, u);
  this->fft_add->ifft(g1, v);

  inv_taylor_expand_t2(output);

  if (beta_m > 1) {
    output->mul_beta(inv_beta_m);
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
     _fft(output, input);
  else {
    // m == 1 -> output = (f(0), f(beta_1))
    output->set(0, input->get(0));
    output->set(1, this->gf->add(input->get(0),
                                 this->gf->mul(input->get(1), this->beta_1)));
  }
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
 * Taylor expansion at (x^2 - x)
 *  Algorithm 1 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 * This function is used for FFT over n=2^m, hence n must be power of 2
 *
 * @param g0: vector for poly of gi0 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most t*2^k
 * @param g1: vector poly of gi1 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most (n-t*2^k)
 * @param n: a power of 2
 * @param do_copy: flag to do a copy of input or not since this function will
 *  modify input vector
 */
template <typename T>
void FFTADD<T>::taylor_expand_t2(Vec<T> *input, int n, bool do_copy)
{
  assert(n >= 1);
  assert(input->get_n() <= n);

  // find k s.t. t2^k < n <= 2 *t2^k
  int _k = find_k(n, 2);

  Vec<T> *_input;
  if (do_copy) {
    _input = new Vec<T>(this->gf, input->get_n());
    _input->copy(input, input->get_n());
  } else
    _input = input;

  _taylor_expand_t2(_input, n, _k, 0);

  // get g0, g1 from mem
  for (int i = 0; i < this->n; i += 2) {
    g0->set(i/2, _input->get(i));
    g1->set(i/2, _input->get(i+1));
  }

  if (do_copy)
    delete _input;
}

/**
 * Taylor expansion at (x^2 - x)
 *  Algorithm 1 in the paper of Shuhong Gao and Todd Mateer:
 *    "Additive Fast Fourier Transforms Over Finite Fields"
 *
 * This function is used for FFT over n=2^m, hence n must be power of 2
 *
 * @param g0: vector for poly of gi0 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most t*2^k
 * @param g1: vector poly of gi1 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most (n-t*2^k)
 *
 * @param n: a power of 2
 * @param k: a number st t*2^k < n <= 2*t*2^k
 * @param s_deg: start degree of g0 and g1
 */
template <typename T>
void FFTADD<T>::_taylor_expand_t2(Vec<T> *input, int n, int _k, int s_deg)
{
  int deg2 = _exp2<T>(_k);
  int deg0 = 2*deg2;
  int deg1 = deg0 - deg2;

  VmVec<T> _g0(input, deg0);
  VmVec<T> _g1(input, n - deg0, deg0);

  // find f0, f1, f2 s.t.
  // f = f0 + x^(t2^k)( f1 + x^( (t-1)2^k) f2 )
  // NOTE: f0 === g0
  VmVec<T> f2(input, deg2, deg0+deg1);

  // since deg1 >= deg2, add f2 into f1 to perform h === f1
  // that is also actually g1
  _g1.add_mutual(&f2);

  // add h (=== f1) into g0 with offset deg2=2^k
  _g0.add_mutual(&_g1, deg2, deg1);

  if (deg0 > 2) {
    _taylor_expand_t2(&_g0, deg0, _k-1, s_deg);
    _taylor_expand_t2(&_g1, n - deg0, _k-1, s_deg + deg2);
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
 * @param g0: vector for poly of gi0 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most t*2^k
 * @param g1: vector poly of gi1 of hi(x) polynomials = gi0 + gi1 * x
 *            The polynomial is of degree at most (n-t*2^k)
 */
template <typename T>
void FFTADD<T>::inv_taylor_expand_t2(Vec<T> *output) {
  assert(g0->get_n() == g1->get_n());
  output->zero_fill();

  int i = output->get_n() / 2 - 1;
  output->set(0, g0->get(i));
  output->set(1, g1->get(i));
  while(--i >= 0) {
    // multiply output to (x^2 - x)
    mul_xt_x(output, 2);
    output->set(0, g0->get(i));
    if (g1->get(i) > 0)
      output->set(1, this->gf->add(output->get(1), g1->get(i)));
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
 */
template <typename T>
void FFTADD<T>::taylor_expand(Vec<T> *output, Vec<T> *input, int n, int t)
{
  assert(n >= 1);
  assert(t > 1);
  assert(input->get_n() <= n);

  if (output->get_n() > input->get_n())
    output->zero_fill();

  // set output as input
  output->copy(input);
  _taylor_expand(output, n, t);
}

template <typename T>
void FFTADD<T>::_taylor_expand(Vec<T> *input, int n, int t)
{
  // find k s.t. t2^k < n <= 2 *t2^k
  int k = find_k(n, t);
  int deg2 = _exp2<T>(k);
  int deg0 = t*deg2;
  int deg1 = deg0 - deg2;
  int g1deg = n - deg0;
  int hdeg = deg1;
  if (g1deg < hdeg)
    hdeg = g1deg;
  VmVec<T> _g0(input, deg0);
  VmVec<T> _g1(input, g1deg, deg0);

  // find f0, f1, f2 s.t.
  // f = f0 + x^(t2^k)( f1 + x^( (t-1)2^k) f2 )
  // NOTE: f0 is stored by _g0, f1 is stored by _g1
  VmVec<T> f2(input, deg2, deg0+deg1);

  // since deg1 >= deg2, add f2 into f1 to perform h == f1 + f2
  // that is also actually g1
  _g1.add_mutual(&f2);

  // add h (=== f1) into g0 with offset deg2=2^k
  _g0.add_mutual(&_g1, deg2, hdeg);

  if (deg0 > t)
    _taylor_expand(&_g0, deg0, t);
  if (g1deg > t)
    _taylor_expand(&_g1, g1deg, t);
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
    hi.set_map(i);
    output->add_mutual(&hi);
  }
}

#endif
