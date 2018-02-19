/* -*- mode: c++ -*- */
#ifndef __NTTEC_FEC_RS_GF2N_H__
#define __NTTEC_FEC_RS_GF2N_H__

#include "fec_base.h"
#include "gf_bin_ext.h"
#include "matrix.h"
#include "vec_vector.h"

namespace nttec {
namespace fec {

enum class RsMatrixType { VANDERMONDE, CAUCHY };

/** Reed-Solomon (RS) Erasure code over GF(2<sup>n</sup>) (Cauchy or
 *  Vandermonde).
 */
template <typename T>
class RsGf2n : public FecCode<T> {
  public:
    RsMatrixType type;

    RsGf2n(
        unsigned word_size,
        unsigned n_data,
        unsigned n_parities,
        RsMatrixType type)
        : FecCode<T>(FecType::SYSTEMATIC, word_size, n_data, n_parities)
    {
        assert(
            type == RsMatrixType::VANDERMONDE || type == RsMatrixType::CAUCHY);

        if (word_size > 16)
            assert(false); // not support yet
        unsigned gf_n = 8 * word_size;
        this->gf = new gf::BinExtension<T>(gf_n);

        this->mat = new Matrix<T>(this->gf, n_parities, n_data);
        if (type == RsMatrixType::CAUCHY) {
            mat->cauchy();
        } else if (type == RsMatrixType::VANDERMONDE) {
            mat->vandermonde_suitable_for_ec();
        }

        // has to be a n_data*n_data invertible square matrix
        decode_mat =
            new Matrix<T>(this->gf, mat->get_n_cols(), mat->get_n_cols());
    }

    ~RsGf2n()
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
        vec::Vector<T>* output,
        std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* words)
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
        vec::Vector<T>* output,
        const std::vector<Properties>& props,
        off_t offset,
        vec::Vector<T>* fragments_ids,
        vec::Vector<T>* words)
    {
        decode_mat->mul(output, words);
    }

  private:
    Matrix<T>* mat = nullptr;
    Matrix<T>* decode_mat = nullptr;
};

} // namespace fec
} // namespace nttec

#endif
