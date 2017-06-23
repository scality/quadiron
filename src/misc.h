/* -*- mode: c++ -*- */
#pragma once

extern GFP<uint32_t> __gf32;
extern GFP<uint64_t> __gf64;
extern GFP<mpz_class> __gfmpz;

#if defined(__i386__)
static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int x;
  __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
  return x;
}
#elif defined(__x86_64__)
static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#endif

std::ostream& operator<<( std::ostream& dest, __uint128_t value );
std::ostream& operator<<( std::ostream& dest, __int128_t value );
