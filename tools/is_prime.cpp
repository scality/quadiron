
#include "ntl.h"
#include "gf.cpp"
#include "gfp.cpp"
#include "gf2n.cpp"

template class GF<mpz_class>;
template class GFP<mpz_class>;

bool _solovay_strassen(GF<mpz_class> *gf, mpz_class a, mpz_class n)
{
  mpz_class _n = (n - 1) / 2;
  mpz_class _j = gf->_mod_exp(a, _n, n);
  int j = gf->_jacobi(a, n);
  
  return j == _j;
}

int main(int argc, char **argv)
{
  GFP<mpz_class> gfp(3);

  if (argc != 2) {
    std::cerr << "usage: solovay-strassen n\n";
    exit(1);
  }

  mpz_class n(argv[1]);

  if (n % 2 == 0) {
    std::cerr << "please choose an odd number\n";
    exit(1);
  }

  int ok = 0;
  for (int a = 0;a < 100;a++) {
    if (_solovay_strassen(&gfp, a, n))
      ok++;
  }
  
  std::cout << ok << "\n";
  exit(0);
}
