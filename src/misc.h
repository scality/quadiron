/* -*- mode: c++ -*- */
#pragma once

extern GFP<uint32_t> __gf32;
extern GFP<uint64_t> __gf64;
extern GFP<mpz_class> __gfmpz;

std::ostream& operator<<( std::ostream& dest, __uint128_t value );
std::ostream& operator<<( std::ostream& dest, __int128_t value );
