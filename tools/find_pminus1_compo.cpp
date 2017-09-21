
#include "ntl.h"

template class GF<mpz_class>;
template class Arith<mpz_class>;

int main(int argc, char **argv)
{
  Arith<mpz_class> arith;

  if (argc != 3) {
    std::cerr << "usage: find_pminus1_compo start_n min_n_factors\n";
    exit(1);
  }

  mpz_class n(argv[1]);
  int min_n_factors(atoi(argv[2]));

  while (true) {
    if (n % 2 == 0) {
      n++;
      continue;
    }
    if (arith.solovay_strassen(n)) {
      // number is probably prime: double-check
      if (arith.is_prime(n)) {
        std::vector<mpz_class> primes;
        std::vector<mpz_class> exponents;
        arith.factor_prime(n-1, &primes, &exponents);
        typename std::vector<mpz_class>::size_type i;
        mpz_class n_factors = 0;
        for (i = 0; i != primes.size(); i++) {
          n_factors += exponents.at(i);
        }
        if (n_factors > min_n_factors) {
          std::cout << n << "-1=";
          for (i = 0; i != primes.size(); i++) {
            std::cout << primes.at(i) << "^" << exponents.at(i) << " ";
          }
          std::cout << std::endl;
        }
      }
    }
    n++;
  }
}
