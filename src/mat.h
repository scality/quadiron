
#ifndef __MAT_H__
#define __MAT_H__ 1

/*
 * mat[n_rows][n_cols]:
 * mat[0][0] mat[0][1] ... mat[0][n_cols]
 * mat[1][0] ...
 * mat[n_rows-1] ...       mat[n_rows][n_cols]
 */
template<typename T>
class Mat
{
 public:
  GF<T> *gf;
  int n_rows;
  int n_cols;
  T *mem;
#define MAT_ITEM(mat, i, j) ((mat)->mem[(i) * (mat)->n_cols + (j)])
  Mat(GF<T> *gf, int n_rows, int n_cols);
  void zero(void);
  void inv(void);
  void mult(Vec<T> *output, Mat<T> *a, Mat<T> *b);
  void vandermonde();
  void vandermonde_transformed();
  void cauchy();
};

#endif
