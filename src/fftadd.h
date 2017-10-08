
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
  Vec<T> *G = NULL;
  Vec<T> *g0 = NULL;
  Vec<T> *g1 = NULL;
  Vec<T> *u = NULL;
  Vec<T> *v = NULL;
  FFTADD<T> *fft_add = NULL;
 public:
  FFTADD(GF<T> *gf, T m, Vec<T>* betas = NULL);
  ~FFTADD();
  void compute_GD();
  void compute_G();
  void compute_B(Vec<T> *B);
  void compute_subspace(Vec<T> *basis, Vec<T> *subspace);
  void fft(Vec<T> *output, Vec<T> *input);
  void ifft(Vec<T> *output, Vec<T> *input);
  void fft_inv(Vec<T> *output, Vec<T> *input);
 private:
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
  if (m > 1) {
    this->k = this->gf->arith->exp2(m - 1);

    // compute gammas and deltas
    this->gammas = new Vec<T>(gf, m-1);
    this->deltas = new Vec<T>(gf, m-1);
    this->G = new Vec<T>(gf, this->k);
    this->compute_GD();

    this->u = new Vec<T>(gf, this->k);
    this->v = new Vec<T>(gf, this->k);
    this->g0 = new Vec<T>(gf, this->k);
    this->g1 = new Vec<T>(gf, this->k);
    this->fft_add = new FFTADD(gf, m-1, this->deltas);
  }
}

template <typename T>
FFTADD<T>::~FFTADD()
{
  if (create_betas && betas) delete betas;
  if (gammas) delete gammas;
  if (deltas) delete deltas;
  if (G) delete G;
  if (u) delete u;
  if (v) delete v;
  if (g0) delete g0;
  if (g1) delete g1;
  if (fft_add) delete fft_add;
}

template <typename T>
void FFTADD<T>::compute_GD()
{
  for (int i = 0; i < m-1; i++) {
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
  Poly<T> g(this->gf);
  T beta = beta_m;
  g.set(0, input->get(0));
  for (i = 1; i < deg; i++) {
    g.set(i, this->gf->mul(beta, input->get(i)));
    beta = this->gf->mul(beta, beta_m);
  }

  // compute taylor expansion of g(x) at (x^2 - x)
  g.taylor_expand_t2(g0, g1, this->get_n());

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
  for (i = 0; i < k; i++) {
    v->set(i, this->gf->add(input->get(k+i), input->get(i)));
    u->set(i, this->gf->add(input->get(i),
                            this->gf->mul(G->get(i), v->get(i))));
  }

  this->fft_add->ifft(g0, u);
  this->fft_add->ifft(g1, v);

  Poly<T> g(this->gf);
  g.inv_taylor_expand_t2(g0, g1);

  T beta = inv_beta_m;
  T deg = g.degree();
  for (i = 1; i <= deg; i++) {
    g.set(i, this->gf->mul(beta, g.get(i)));
    beta = this->gf->mul(beta, inv_beta_m);
  }

  g.to_vec(output);
}

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

template <typename T>
void FFTADD<T>::ifft(Vec<T> *output, Vec<T> *input)
{
  fft_inv(output, input);
}
