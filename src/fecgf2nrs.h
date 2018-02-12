/* -*- mode: c++ -*- */
#ifndef __NTTEC_FECGF2NRS_H__
#define __NTTEC_FECGF2NRS_H__

#include "fec.h"
#include "gf2n.h"
#include "mat.h"

namespace nttec {
namespace fec {

/**
 * GF_2^n based RS (Cauchy or Vandermonde)
 */
template <typename T>
class FECGF2NRS : public FEC<T> {
  public:
    enum FECGF2NRSType {
        VANDERMONDE,
        CAUCHY,
    };

    FECGF2NRSType type;

  private:
    Mat<T>* mat = nullptr;
    Mat<T>* decode_mat = nullptr;

  public:
    FECGF2NRS(
        unsigned word_size,
        unsigned n_data,
        unsigned n_parities,
        FECGF2NRSType type)
        : FEC<T>(FEC<T>::TYPE_1, word_size, n_data, n_parities)
    {
        assert(type == VANDERMONDE || type == CAUCHY);

        if (word_size > 16)
            assert(false); // not support yet
        unsigned gf_n = 8 * word_size;
        this->gf = new gf::GF2N<T>(gf_n);

        this->mat = new Mat<T>(this->gf, n_parities, n_data);
        if (type == CAUCHY) {
            mat->cauchy();
        } else if (type == VANDERMONDE) {
            mat->vandermonde_suitable_for_ec();
        }

        // has to be a n_data*n_data invertible square matrix
        decode_mat = new Mat<T>(this->gf, mat->get_n_cols(), mat->get_n_cols());
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

    void encode(
        vec::Vec<T>* output,
        std::vector<KeyValue*> props,
        off_t offset,
        vec::Vec<T>* words)
    {
        mat->mul(output, words);
    }

    void decode_add_data(int fragment_index, int row)
    {
        // for each data available generate the corresponding identity
        for (int j = 0; j < mat->get_n_cols(); j++) {
            if (row == j)
                decode_mat->set(fragment_index, j, 1);
            else
                decode_mat->set(fragment_index, j, 0);
        }
    }

    void decode_add_parities(int fragment_index, int row)
    {
        // copy corresponding row in vandermonde matrix
        for (int j = 0; j < mat->get_n_cols(); j++) {
            decode_mat->set(fragment_index, j, mat->get(row, j));
        }
    }

    void decode_build()
    {
        decode_mat->inv();
    }

    void decode(
        vec::Vec<T>* output,
        std::vector<KeyValue*> props,
        off_t offset,
        vec::Vec<T>* fragments_ids,
        vec::Vec<T>* words)
    {
        decode_mat->mul(output, words);
    }
};

} // namespace fec
} // namespace nttec

#endif
