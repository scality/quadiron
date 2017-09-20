
#include "ntl.h"

Arith<uint32_t> __arith32;
Arith<uint64_t> __arith64;
Arith<mpz_class> __arithmpz;

std::ostream&
operator<<(std::ostream& dest, __uint128_t value)
{
  std::ostream::sentry s(dest);
  if (s) {
    __uint128_t tmp = value;
    char buffer[ 128 ];
    char* d = std::end(buffer);
    do
      {
        --d;
        *d = "0123456789"[ tmp % 10 ];
        tmp /= 10;
      } while (tmp != 0);

    int len = std::end(buffer) - d;
    if (dest.rdbuf()->sputn(d, len) != len) {
      dest.setstate(std::ios_base::badbit);
    }
  }
  return dest;
}

std::ostream&
operator<<(std::ostream& dest, __int128_t value)
{
  std::ostream::sentry s(dest);
  if (s) {
    __uint128_t tmp = value < 0 ? -value : value;
    char buffer[ 128 ];
    char* d = std::end(buffer);
    do
      {
        --d;
        *d = "0123456789"[ tmp % 10 ];
        tmp /= 10;
      } while (tmp != 0);
    if (value < 0) {
      --d;
      *d = '-';
    }
    int len = std::end(buffer) - d;
    if (dest.rdbuf()->sputn(d, len) != len) {
      dest.setstate(std::ios_base::badbit);
    }
  }
  return dest;
}
