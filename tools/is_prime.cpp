
#include "ntl.h"
#include "gf.cpp"
#include "gfp.cpp"
#include "gf2n.cpp"

template class GF<mpz_class>;
template class GFP<mpz_class>;

int main(int argc, char **argv)
{
  GFP<mpz_class> gfp(3);

  if (argc != 2) {
    std::cerr << "usage: is_prime n\n";
    exit(1);
  }

  mpz_class n(argv[1]);

  if (n % 2 == 0) {
    std::cerr << "please choose an odd number\n";
    exit(1);
  }

  bool result = gfp._solovay_strassen(n); 
  std::cerr << result << "\n";
  exit(result == true ? 0 : 1);
}
