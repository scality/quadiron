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

  FECGF2NRS(GF<T> *gf, u_int word_size, u_int n_data, u_int n_parities, FECGF2NRSType type) : 
    FEC<T>(gf, FEC<T>::TYPE_1, word_size, n_data, n_parities)
  {
    assert(type == VANDERMONDE || type == CAUCHY);
    this->type = type;

    this->mat = new Mat<T>(gf, n_parities, n_data);
    if (type == CAUCHY) {
      mat->cauchy();
    } else if (type == VANDERMONDE) {
      mat->vandermonde_suitable_for_ec();
    }

    //has to be a n_data*n_data invertible square matrix
    decode_mat = new Mat<T>(this->gf, mat->n_cols, mat->n_cols);
  }

  ~FECGF2NRS()
  {
    delete mat;
    delete decode_mat;
  }

  int get_n_fragments_required()
  {
    return this->n_data;
  } 

  int get_n_inputs()
  {
    return this->n_data;
  }

  int get_n_outputs()
  {
    return this->n_parities;
  }

  void encode(std::vector<KeyValue*> props, off_t offset, Vec<T> *output, Vec<T> *words)
  {
    mat->mul(output, words);
  }

  void decode_add_data(int fragment_index, int row)
  {
    //for each data available generate the corresponding identity
    for (int j = 0;j < mat->n_cols;j++) {
      if (row == j)
        decode_mat->set(fragment_index, j, 1);
      else
        decode_mat->set(fragment_index, j, 0);
    }
  }

  void decode_add_parities(int fragment_index, int row)
  {
    //copy corresponding row in vandermonde matrix
    for (int j = 0;j < mat->n_cols;j++) {
      decode_mat->set(fragment_index, j, mat->get(row, j));
    }
  }

  void decode_build()
  {
    decode_mat->inv();
  }

  void decode(std::vector<KeyValue*> props, off_t offset, Vec<T> *output, Vec<T> *words)
  {
    decode_mat->mul(output, words);
  }
};
