
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
  void mult(Vec<T> *output, Mat<T> *a, Vec<T> *b);
  void vandermonde(void);
  void vandermonde_suitable_for_ec(void);
  void cauchy(void);
  void dump(void);
 private:
  void ec_transform1(Mat<T> *tmp, int i);
  void ec_transform2(Mat<T> *tmp, int i, int j);
  bool row_is_identity(Mat<T> *tmp, int row);
};

#endif
