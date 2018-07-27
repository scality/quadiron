/* -*- mode: c++ -*- */
/*
 * Copyright 2017-2018 Scality
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef __QUAD_VEC_MATRIX_H__
#define __QUAD_VEC_MATRIX_H__

#include <iostream>

#include "gf_ring.h"
#include "vec_vector.h"

namespace quadiron {
namespace vec {

/** A 2D matrix
 *
 * A 2D matrix such as:
 *
 * \f[
 * Matrix(rows, cols) =
 * \begin{pmatrix}
 *   matrix_{0, 0}        & matrix_{0, 1}   & \cdots & matrix_{0, cols - 1}   \\
 *   matrix_{1, 0}        & matrix_{1, 1}   & \cdots & matrix_{1, cols - 1}   \\
 *   \vdots               & \vdots          & \ddots & \vdots                 \\
 *   matrix_{rows - 1 ,0} & a_{rows - 1, 2} & \cdots & a_{rows - 1, cols - 1}
 * \end{pmatrix}
 * \f]
 */
template <typename T>
class Matrix {
  public:
    Matrix(const gf::RingModN<T>& rn, int n_rows, int n_cols);
    virtual ~Matrix();
    virtual int get_n_rows();
    virtual int get_n_cols();
    void zero_fill(void);
    void set(int i, int j, T val);
    virtual T get(int i, int j);
    void inv(void);
    void mul(vec::Vector<T>* output, vec::Vector<T>* v);
    void vandermonde(void);
    void vandermonde_suitable_for_ec(void);
    void cauchy(void);
    void dump_row(int row);
    void dump(void);

  private:
    const gf::RingModN<T>* rn;
    int n_rows;
    int n_cols;
    T* mem;
    void swap_rows(int row0, int row1);
    void mul_row(int row, T factor);
    void add_rows(int src_row, int dst_row, T factor);
    void reduced_row_echelon_form(void);
    void ec_transform1(int i);
    void ec_transform2(int i, int j);
    bool row_is_identity(int row);
};

template <typename T>
Matrix<T>::Matrix(const gf::RingModN<T>& rn, int n_rows, int n_cols)
{
    this->rn = &rn;
    this->n_rows = n_rows;
    this->n_cols = n_cols;
    this->mem = new T[n_rows * n_cols];
}

template <typename T>
Matrix<T>::~Matrix()
{
    delete[] this->mem;
}

template <typename T>
int Matrix<T>::get_n_rows(void)
{
    return n_rows;
}

template <typename T>
int Matrix<T>::get_n_cols(void)
{
    return n_cols;
}

template <typename T>
void Matrix<T>::zero_fill(void)
{
    int i, j;

    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++) {
            set(i, j, 0);
        }
    }
}

template <typename T>
void Matrix<T>::set(int i, int j, T val)
{
    assert(i >= 0 && i < n_rows);
    assert(j >= 0 && j < n_cols);

    mem[i * n_cols + j] = val;
}

template <typename T>
T Matrix<T>::get(int i, int j)
{
    assert(i >= 0 && i < n_rows);
    assert(j >= 0 && j < n_cols);

    return mem[i * n_cols + j];
}

template <typename T>
void Matrix<T>::swap_rows(int row0, int row1)
{
    int j;
    T tmp;

    assert(row0 >= 0 && row0 < n_rows);
    assert(row1 >= 0 && row1 < n_rows);

    for (j = 0; j < n_cols; j++) {
        tmp = get(row0, j);
        set(row0, j, get(row1, j));
        set(row1, j, tmp);
    }
}

template <typename T>
void Matrix<T>::mul_row(int row, T factor)
{
    int j;

    assert(row >= 0 && row < n_rows);

    for (j = 0; j < n_cols; j++) {
        set(row, j, rn->mul(get(row, j), factor));
    }
}

template <typename T>
void Matrix<T>::add_rows(int src_row, int dst_row, T factor)
{
    int j;

    assert(src_row >= 0 && src_row < n_rows);
    assert(dst_row >= 0 && dst_row < n_rows);

    for (j = 0; j < n_cols; j++) {
        set(dst_row,
            j,
            rn->add(get(dst_row, j), rn->mul(get(src_row, j), factor)));
    }
}

template <typename T>
void Matrix<T>::reduced_row_echelon_form(void)
{
    int i, j;

    // Compute row echelon form (REF)
    int num_pivots = 0;
    for (j = 0; j < n_cols && num_pivots < n_rows; j++) { // For each column
        // Find a pivot row for this column
        int pivot_row = num_pivots;
        while (pivot_row < n_rows && get(pivot_row, j) == 0)
            pivot_row++;
        if (pivot_row == n_rows)
            continue; // Cannot eliminate on this column
        swap_rows(num_pivots, pivot_row);
        pivot_row = num_pivots;
        num_pivots++;

        // Simplify the pivot row
        mul_row(pivot_row, rn->div(1, get(pivot_row, j)));

        // Eliminate rows below
        for (i = pivot_row + 1; i < n_rows; i++)
            add_rows(pivot_row, i, rn->neg(get(i, j)));
    }

    // Compute reduced row echelon form (RREF)
    for (int i = num_pivots - 1; i >= 0; i--) {
        // Find pivot
        int pivot_col = 0;
        while (pivot_col < n_cols && get(i, pivot_col) == 0)
            pivot_col++;
        if (pivot_col == n_cols)
            continue; // Skip this all-zero row

        // Eliminate rows above
        for (int j = i - 1; j >= 0; j--)
            add_rows(i, j, rn->neg(get(j, pivot_col)));
    }
}

/**
 * taken from
 * https://www.nayuki.io/res/gauss-jordan-elimination-over-any-field/Matrix.java
 */
template <typename T>
void Matrix<T>::inv(void)
{
    int i, j;

    assert(n_rows == n_cols);

    // Build augmented matrix: [this | identity]
    Matrix<T> tmp(*rn, n_rows, n_cols * 2);
    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++) {
            tmp.set(i, j, get(i, j));
            tmp.set(i, j + n_cols, (i == j) ? 1 : 0);
        }
    }

    // Do the main calculation
    tmp.reduced_row_echelon_form();

    // Check that the left half is the identity matrix
    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++) {
            if (tmp.get(i, j) != (i == j) ? 1 : 0)
                throw DomainError("non invertible matrix");
        }
    }

    // Extract inverse matrix from: [identity | inverse]
    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++)
            set(i, j, tmp.get(i, j + n_cols));
    }
}

template <typename T>
bool Matrix<T>::row_is_identity(int row)
{
    int j;

    assert(row >= 0 && row < n_rows);

    for (j = 0; j < n_cols; j++) {
        if (get(row, j) != ((j == row) ? 1 : 0))
            return false;
    }

    return true;
}

template <typename T>
void Matrix<T>::vandermonde(void)
{
    int i, j;

    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++) {
            set(i, j, rn->exp(i, j));
        }
    }
}

/**
 * transform c_i into f_i_i_minus_1 * c_i
 */
template <typename T>
void Matrix<T>::ec_transform1(int i)
{
    int k;

    assert(i >= 0 && i < n_rows);
    assert(i >= 0 && i < n_cols);

    T f_minus_1 = rn->inv(get(i, i));

    for (k = 0; k < n_rows; k++) {
        set(k, i, rn->mul(f_minus_1, get(k, i)));
    }
}

/**
 * transform c_j into c_j - f_i_j * c_i
 */
template <typename T>
void Matrix<T>::ec_transform2(int i, int j)
{
    int k;

    assert(i >= 0 && i < n_rows);
    assert(j >= 0 && i < n_cols);

    T f_i_j = get(i, j);

    for (k = 0; k < n_rows; k++) {
        set(k, j, rn->sub(get(k, j), rn->mul(f_i_j, get(k, i))));
    }
}

/**
 * see http://web.eecs.utk.edu/~plank/plank/papers/CS-03-504.html
 */
template <typename T>
void Matrix<T>::vandermonde_suitable_for_ec(void)
{
    int i, j, dim;

    dim = n_rows + n_cols;

    Matrix<T> tmp(*rn, dim, n_cols);

    tmp.vandermonde();

    /* perform transformations to get the identity matrix on the top rows */
    i = 0;
    while (i < n_cols) {
        if (tmp.row_is_identity(i)) {
            i++;
            continue;
        }

        // this case is mentionned in the paper but cannot happen
        /* if (0 == tmp.get(i, i)) {
           for (j = i + 1; j < &tmp->n_cols; j++) {
           if (0 != tmp.get(i, j)) {
           tmp.swap_cols(i, j);
           continue ;
           }
           } */

        // check if f_i_i == 1
        if (1 != tmp.get(i, i)) {
            // check for inverse since f_i_i != 0
            tmp.ec_transform1(i);
        }

        // now f_i_i == 1
        for (j = 0; j < tmp.n_cols; j++) {
            if (i != j) {
                if (0 != tmp.get(i, j)) {
                    tmp.ec_transform2(i, j);
                }
            }
        }

        i++;
    }

    // copy last n_rows rows of tmp into mat
    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++) {
            set(i, j, tmp.get(n_cols + i, j));
        }
    }
}

template <typename T>
void Matrix<T>::mul(vec::Vector<T>* output, vec::Vector<T>* v)
{
    int i, j;

    assert(get_n_cols() == v->get_n());
    assert(get_n_rows() == output->get_n());

    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++) {
            T x = rn->mul(get(i, j), v->get(j));
            if (0 == j)
                output->set(i, x);
            else
                output->set(i, rn->add(output->get(i), x));
        }
    }
}

template <typename T>
void Matrix<T>::cauchy()
{
    int i, j;

    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++) {
            set(i, j, rn->inv(rn->add(i, (j + n_rows))));
        }
    }

    /* do optimise */
    // convert 1st row to all 1s
    for (j = 0; j < n_cols; j++) {
        for (i = 0; i < n_rows; i++) {
            set(i, j, rn->div(get(i, j), get(0, j)));
        }
    }
    // convert 1st element of each row to 1
    for (i = 1; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++) {
            set(i, j, rn->div(get(i, j), get(i, 0)));
        }
    }
}

template <typename T>
void Matrix<T>::dump_row(int row)
{
    assert(row >= 0 && row < n_rows);

    for (int j = 0; j < n_cols; j++) {
        std::cout << " " << get(row, j);
    }
    std::cout << "\n";
}

template <typename T>
void Matrix<T>::dump(void)
{
    int i, j;

    std::cout << "--\n";
    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++) {
            std::cout << " " << get(i, j);
        }
        std::cout << "\n";
    }
}

} // namespace vec
} // namespace quadiron

#endif
