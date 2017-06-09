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

private:
  Mat<T> *mat = NULL;
  Mat<T> *repair_mat = NULL;

public:

  FECGF2NRS(GF<T> *gf, u_int word_size, u_int n_data, u_int n_parities, FECGF2NRSType type) : 
    FEC<T>(gf, FEC<T>::TYPE_1, word_size, n_data, n_parities)
  {
    this->mat = new Mat<T>(gf, n_parities, n_data);
    if (type == CAUCHY) {
      mat->cauchy();
    } else if (type == VANDERMONDE) {
      mat->vandermonde_suitable_for_ec();
    } else {
      throw NTL_EX_INVAL;
    }
  }

  ~FECGF2NRS()
  {
    delete mat;
    delete repair_mat;
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

  void encode(Vec<T> *output, Vec<T> *words)
  {
    mat->mult(output, words);
  }
  
  void repair_init(void)
  {
    //has to be a n_data*n_data invertible square matrix
    repair_mat = new Mat<T>(this->gf, mat->n_cols, mat->n_cols);
  }

  void repair_add_data(int fragment_index, int row)
  {
    //for each data available generate the corresponding identity
    for (int j = 0;j < mat->n_cols;j++) {
      if (row == j)
        repair_mat->set(fragment_index, j, 1);
      else
        repair_mat->set(fragment_index, j, 0);
    }
  }

  void repair_add_parities(int fragment_index, int row)
  {
    //copy corresponding row in vandermonde matrix
    for (int j = 0;j < mat->n_cols;j++) {
      repair_mat->set(fragment_index, j, mat->get(row, j));
    }
  }

  void repair_build()
  {
    repair_mat->inv();
  }

  void repair(Vec<T> *output, Vec<T> *words)
  {
    repair_mat->mult(output, words);
  }
};
