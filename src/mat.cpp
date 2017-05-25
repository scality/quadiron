
#include "ntl.h"

template <typename T>
Mat<T>::Mat(GF<T> *gf, int n_rows, int n_cols)
{
  this->gf = gf;
  this->n_rows = n_rows;
  this->n_cols = n_cols;
  this->mem = new T[n_rows * n_cols];
}

template <typename T>
void Mat<T>::zero_fill(void)
{
  int i, j;

  for (i = 0;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(this, i, j) = 0;
    }
  }
}

template <typename T>
void Mat<T>::swap_rows(int row0, int row1)
{
  int j;
  T tmp;
  
  assert(row0 >= 0 && row0 < n_rows);
  assert(row1 >= 0 && row1 < n_rows);

  for (j = 0;j < n_cols;j++) {
    tmp = MAT_ITEM(this, row0, j);
    MAT_ITEM(this, row0, j) = MAT_ITEM(this, row1, j);
    MAT_ITEM(this, row1, j) = tmp;
  }
}

template <typename T>
void Mat<T>::mul_row(int row, T factor)
{
  int j;

  assert(row >= 0 && row < n_rows);

  for (j = 0; j < n_cols; j++) {
    MAT_ITEM(this, row, j) = gf->mul(MAT_ITEM(this, row, j), factor);
  }
}

template <typename T>
void Mat<T>::add_rows(int src_row, int dst_row, T factor)
{
  int j;

  assert(src_row >= 0 && src_row < n_rows);
  assert(dst_row >= 0 && dst_row < n_rows);

  for (j = 0; j < n_cols; j++) {
    MAT_ITEM(this, dst_row, j) = 
      gf->add(MAT_ITEM(this, dst_row, j),
              gf->mul(MAT_ITEM(this, src_row, j), factor));
  }
}

template <typename T>
void Mat<T>::reduced_row_echelon_form(void)
{
  int i, j;

  // Compute row echelon form (REF)
  int num_pivots = 0;
  for (j = 0; j < n_cols && num_pivots < n_rows; j++) { // For each column
    // Find a pivot row for this column
    int pivot_row = num_pivots;
    while (pivot_row < n_rows && 
           MAT_ITEM(this, pivot_row, j) == 0)
      pivot_row++;
    if (pivot_row == n_rows)
      continue;  // Cannot eliminate on this column
    swap_rows(num_pivots, pivot_row);
    pivot_row = num_pivots;
    num_pivots++;
    
    // Simplify the pivot row
    mul_row(pivot_row, gf->div(1, MAT_ITEM(this, pivot_row, j)));
    
    // Eliminate rows below
    for (i = pivot_row + 1; i < n_rows; i++)
      add_rows(pivot_row, i, gf->neg(MAT_ITEM(this, i, j)));
  }

  // Compute reduced row echelon form (RREF)
  for (int i = num_pivots - 1; i >= 0; i--) {
    // Find pivot
    int pivot_col = 0;
    while (pivot_col < n_cols && 
           MAT_ITEM(this, i, pivot_col) == 0)
      pivot_col++;
    if (pivot_col == n_cols)
      continue;  // Skip this all-zero row
    
    // Eliminate rows above
    for (int j = i - 1; j >= 0; j--)
      add_rows(i, j, gf->neg(MAT_ITEM(this, j, pivot_col)));
  }
}

/** 
 * taken from 
 * https://www.nayuki.io/res/gauss-jordan-elimination-over-any-field/Matrix.java
 */
template <typename T>
void Mat<T>::inv(void)
{
  int i, j;

  assert(n_rows == n_cols);

  // Build augmented matrix: [this | identity]
  Mat<T> tmp(this->gf, n_rows, n_cols * 2);
  for (i = 0; i < n_rows; i++) {
    for (j = 0; j < n_cols; j++) {
      MAT_ITEM(&tmp, i, j) = MAT_ITEM(this, i, j);
      MAT_ITEM(&tmp, i, j + n_cols) = (i == j) ? 1 : 0;
    }
  }

  // Do the main calculation
  tmp.reduced_row_echelon_form();
  
  // Check that the left half is the identity matrix
  for (i = 0; i < n_rows; i++) {
    for (j = 0; j < n_cols; j++) {
      if (MAT_ITEM(&tmp, i, j) != (i == j) ? 1 : 0)
        throw NTL_EX_MAT_NOT_INVERTIBLE;
    }
  }
  
  // Extract inverse matrix from: [identity | inverse]
  for (i = 0; i < n_rows; i++) {
    for (j = 0; j < n_cols; j++)
      MAT_ITEM(this, i, j) = MAT_ITEM(&tmp, i, j + n_cols);
  }
}

template <typename T>
bool Mat<T>::row_is_identity(int row)
{
  int j;

  assert(row >= 0 && row < n_rows);

  for (j = 0;j < n_cols;j++) {
    if (MAT_ITEM(this, row, j) != ((j == row) ? 1 : 0))
      return false;
  }

  return true;
}

template <typename T>
void Mat<T>::vandermonde(void)
{
  int i, j;

  for (i = 0;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(this, i, j) = gf->exp(i, j); 
    }
  }
}

/**
 * transform c_i into f_i_i_minus_1 * c_i 
 */
template <typename T>
void Mat<T>::ec_transform1(int i)
{
  int k;

  assert(i >= 0 && i < n_rows);
  assert(i >= 0 && i < n_cols);

  T f_minus_1 = gf->inv(MAT_ITEM(this, i, i));

  for (k = 0;k < n_rows;k++) {
    MAT_ITEM(this, k, i) = 
      gf->mul(f_minus_1, 
              MAT_ITEM(this, k, i));
  }
}

/**
 * transform c_j into c_j - f_i_j * c_i 
 */
template <typename T>
void Mat<T>::ec_transform2(int i, int j)
{
  int k;

  assert(i >= 0 && i < n_rows);
  assert(j >= 0 && i < n_cols);

  T f_i_j = MAT_ITEM(this, i, j);

  for (k = 0;k < n_rows;k++) {
    MAT_ITEM(this, k, j) = 
      gf->sub(MAT_ITEM(this, k, j),
              gf->mul(f_i_j, MAT_ITEM(this, k, i)));
  }
}

/** 
 * see http://web.eecs.utk.edu/~plank/plank/papers/CS-03-504.html
 */
template <typename T>
void Mat<T>::vandermonde_suitable_for_ec(void)
{
  int i, j, dim;
  
  dim = n_rows + n_cols;

  Mat<T> tmp(gf, dim, n_cols);

  tmp.vandermonde();

  /* perform transformations to get the identity matrix on the top rows */
  i = 0;
  while (i < n_cols) {
    
    if (tmp.row_is_identity(i)) {
      i++;
      continue ;
    }

    //this case is mentionned in the paper but cannot happen
    /* if (0 == MAT_ITEM(&tmp, i, i)) {
       for (j = i + 1;j < &tmp->n_cols;j++) {
       if (0 != MAT_ITEM(&tmp, i, j)) {
       tmp.swap_cols(i, j);
       continue ;
       }
       } */

    //check if f_i_i == 1
    if (1 != MAT_ITEM(&tmp, i, i)) {
      //check for inverse since f_i_i != 0
      tmp.ec_transform1(i);
    }

    //now f_i_i == 1
    for (j = 0;j < tmp.n_cols;j++) {
      if (i != j) {
        if (0 != MAT_ITEM(&tmp, i, j)) {
          tmp.ec_transform2(i, j);
        }
      }
    }

    i++;
  }

  //copy last n_rows rows of tmp into mat
  for (i = 0;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(this, i, j) = MAT_ITEM(&tmp, n_cols + i, j);
    }
  }
}

template <typename T>
void Mat<T>::mult(Vec<T> *output, Vec<T> *v)
{
  int i, j;

  assert(output->n == v->n);
  assert(output->n == n_cols);
  for (i = 0;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      T x = gf->mul(MAT_ITEM(this, i, j), VEC_ITEM(v, j));
      if (0 == j)
        VEC_ITEM(output, i) = x;
      else
        VEC_ITEM(output, i) = gf->add(VEC_ITEM(output, i), x);
    }
  }
}

template <typename T>
void Mat<T>::cauchy()
{
  int i, j;

  for (i = 0;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(this, i, j) = gf->inv(gf->add(i, (j + n_rows)));
    }
  }

  /* do optimise */
  // convert 1st row to all 1s
  for (j = 0;j < n_cols;j++) {
    for (i = 0;i < n_rows;i++) {
      MAT_ITEM(this, i, j) = gf->div(MAT_ITEM(this, i, j), MAT_ITEM(this, 0, j));
    }
  }
  // convert 1st element of each row to 1
  for (i = 1;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(this, i, j) = gf->div(MAT_ITEM(this, i, j), MAT_ITEM(this, i, 0));
    }
  }
}

template <typename T>
void Mat<T>::dump(void)
{
  int i, j;
  
  std::cout << "--\n";
  for (i = 0;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      std::cout << " " << MAT_ITEM(this, i, j);
    }
    std::cout << "\n";
  }
}
