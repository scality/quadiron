/* -*- mode: c++ -*- */
#pragma once

/**
 * GF_2^n based RS (Cauchy or Vandermonde)
 */
template<typename T>
class FECGF2NRS : public FEC<T>
{
 public:
  enum FECGF2NRSType
  {
    VANDERMONDE,
    CAUCHY,
  };

  FECGF2NRSType type;

 private:
  Mat<T> *mat = NULL;
  Mat<T> *decode_mat = NULL;

 public:
  FECGF2NRS(u_int word_size, u_int n_data, u_int n_parities,
    FECGF2NRSType type) :
    FEC<T>(FEC<T>::TYPE_1, word_size, n_data, n_parities)
  {
    assert(type == VANDERMONDE || type == CAUCHY);
    this->type = type;

    if (word_size > 16)
      assert(false);  // not support yet
    u_int gf_n = 8*word_size;
    this->gf = new GF2N<T>(gf_n);

    this->mat = new Mat<T>(this->gf, n_parities, n_data);
    if (type == CAUCHY) {
      mat->cauchy();
    } else if (type == VANDERMONDE) {
      mat->vandermonde_suitable_for_ec();
    }

    // has to be a n_data*n_data invertible square matrix
    decode_mat = new Mat<T>(this->gf, mat->n_cols, mat->n_cols);
  }

  ~FECGF2NRS()
  {
    delete mat;
    delete decode_mat;
    delete this->gf;
  }

  int get_n_outputs()
  {
    return this->n_parities;
  }

  void encode(Vec<T> *output, std::vector<KeyValue*> props, off_t offset,
    Vec<T> *words)
  {
    mat->mul(output, words);
  }

  void decode_add_data(int fragment_index, int row)
  {
    // for each data available generate the corresponding identity
    for (int j = 0; j < mat->n_cols; j++) {
      if (row == j)
        decode_mat->set(fragment_index, j, 1);
      else
        decode_mat->set(fragment_index, j, 0);
    }
  }

  void decode_add_parities(int fragment_index, int row)
  {
    // copy corresponding row in vandermonde matrix
    for (int j = 0; j < mat->n_cols; j++) {
      decode_mat->set(fragment_index, j, mat->get(row, j));
    }
  }

  void decode_build()
  {
    decode_mat->inv();
  }

  void decode(Vec<T> *output, std::vector<KeyValue*> props, off_t offset,
    Vec<T> *fragments_ids, Vec<T> *words)
  {
    decode_mat->mul(output, words);
  }
};
