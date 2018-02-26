/*
 * Copyright 2017-2018 the NTTEC authors
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
#include "nttec.h"

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
