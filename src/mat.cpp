
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
void Mat<T>::inv()
{
  int dim, i, j, k, tpos, tval, r;

  assert(n_rows == n_cols);
  dim = n_rows;

  Mat<T> aug(this->gf, dim, dim * 2);

  for (i = 0;i < dim;i++) {
    for (j = 0;j < dim;j++) {
      MAT_ITEM(&aug, i, j) = MAT_ITEM(this, i, j);
    }
  }

  for (i = 0;i < dim;i++) {
    for (j = dim;j < 2*dim;j++) {
      if (i == (j % dim))
        MAT_ITEM(&aug, i, j) = 1;
      else
        MAT_ITEM(&aug, i, j) = 0;
    }
  }

  /* using gauss-jordan elimination */
  for (j = 0;j < dim;j++) {
    tpos = j;
    
    /* finding maximum jth column element in last (dimension-j) rows */
    for (i = j+1;i < dim;i++) {
      if (MAT_ITEM(&aug, i, j) > MAT_ITEM(&aug, tpos, j))
        tpos = i;
    }
    
    /* swapping row which has maximum jth column element */
    if (tpos != j) {
      for (k=0;k < 2*dim;k++) {
        tval = MAT_ITEM(&aug, j, k);
        MAT_ITEM(&aug, j, k) = MAT_ITEM(&aug, tpos, k);
        MAT_ITEM(&aug, tpos, k) = tval;
      }
    }
    
    /* performing row operations to form required identity matrix out 
       of the input matrix */
    for (i = 0;i < dim;i++) {
      if (i != j) {
        r = MAT_ITEM(&aug, i, j);
        for (k = 0;k < 2*dim;k++) {
          MAT_ITEM(&aug, i, k) =
            gf->add(MAT_ITEM(&aug, i, k), 
                    gf->mul(gf->div(MAT_ITEM(&aug, j, k), 
                                    MAT_ITEM(&aug, j, j)),
                            r));
        }
      } else {
        r = MAT_ITEM(&aug, i, j);
        for (k = 0;k < 2*dim;k++) {
          MAT_ITEM(&aug, i, k) = gf->div(MAT_ITEM(&aug, i, k), r);
        }
      }
    }
  }

  for (i = 0;i < dim;i++) {
    for (j = 0;j < dim;j++) {
      MAT_ITEM(this, i, j) = MAT_ITEM(&aug, i, dim + j);
    }
  }
}

template <typename T>
bool Mat<T>::row_is_identity(Mat<T> *tmp, int row)
{
  int j;

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
      MAT_ITEM(this, i, j) = gf->pow(j + 1, i); 
    }
  }
}

/**
 * transform c_i into f_i_i_minus_1 * c_i 
 */
template <typename T>
void Mat<T>::ec_transform1(Mat<T> *tmp, int i)
{
  int k;
  int f_minus_1 = gf->div(1, MAT_ITEM(tmp, i, i));

  for (k = 0;k < tmp->n_rows;k++) {
    MAT_ITEM(tmp, k, i) = gf->mul(f_minus_1, MAT_ITEM(tmp, k, i));
  }
}

/**
 * transform c_j into c_j - f_i_j * c_i 
 */
template <typename T>
void Mat<T>::ec_transform2(Mat<T> *tmp, int i, int j)
{
  int k;
  int f_i_j = MAT_ITEM(tmp, i, j);

  for (k = 0;k < tmp->n_rows;k++) {
    MAT_ITEM(tmp, k, j) = 
      gf->add(MAT_ITEM(tmp, k, j),
              gf->mul(f_i_j, MAT_ITEM(tmp, k, i)));
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

  for (i = 0;i < dim;i++) {
    for (j = 0;j < n_cols;j++) {
      MAT_ITEM(&tmp, i, j) = gf->pow(i, j); 
    }
  }

  /* perform transformations to get the identity matrix on the top rows */
  i = 0;
  while (i < n_cols) {
    
    if (row_is_identity(&tmp, i)) {
      i++;
      continue ;
    }
    
    //this case is mentionned in the paper but cannot happen
    /* if (0 == MAT_ITEM(&tmp, i, i)) {
       for (j = i + 1;j < &tmp->n_cols;j++) {
       if (0 != MAT_ITEM(&tmp, i, j)) {
       mat_swap_cols(&tmp, i, j);
       continue ;
       }
       } */
    
    //check if f_i_i == 1
    if (1 != MAT_ITEM(&tmp, i, i)) {
      //check for inverse since f_i_i != 0
      ec_transform1(&tmp, i);
    }
    
    //now f_i_i == 1
    for (j = 0;j < tmp.n_cols;j++) {
      if (i != j) {
        if (0 != MAT_ITEM(&tmp, i, j)) {
          ec_transform2(&tmp, i, j);
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
void Mat<T>::mult(Vec<T> *output, Mat<T> *a, Vec<T> *b)
{
  int i, j;

  assert(output->n == b->n);
  assert(output->n == a->n_cols);
  for (i = 0;i < a->n_rows;i++) {
    for (j = 0;j < a->n_cols;j++) {
      int x = gf->mul(MAT_ITEM(a, i, j), VEC_ITEM(b, j));
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
      MAT_ITEM(this, i, j) = gf->div(1, (gf->add(i, (j + n_rows))));
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
  
  for (i = 0;i < n_rows;i++) {
    for (j = 0;j < n_cols;j++) {
      std::cout << " " << MAT_ITEM(this, i, j);
    }
    std::cout << "\n";
  }
}
