
#include "ntl.h"
#include <signal.h>

template class GF<mpz_class>;
template class GFP<mpz_class>;

double tot = 0;
double ok = 0;

void sigint(int foo)
{
  std::cerr << "total=" << tot << " ok=" << ok << " succ_rate=" << ok / tot * 100 << "%\n";
  exit(1);
}

int main(int argc, char **argv)
{
  GFP<mpz_class> gfp(3);

  if (argc != 2) {
    std::cerr << "usage: find_primes n_bits\n";
    exit(1);
  }

  int n_bits = atoi(argv[1]);
  gmp_randclass r(gmp_randinit_default);

  signal(SIGINT, sigint);

  for (int i = 0;i < 100000;i++) {
    mpz_class n = r.get_z_bits(n_bits);
    if (n % 2 == 0)
      continue ;
    bool result = gfp._solovay_strassen(n);
    tot++;
    if (result) {
      ok++;
      std::cout << n << "\n";
    }
  }
  exit(0);
}
