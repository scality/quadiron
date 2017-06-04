/* -*- mode: c++ -*- */
#pragma once

extern GFP<uint64_t> __gf;

std::ostream& operator<<( std::ostream& dest, __uint128_t value );
std::ostream& operator<<( std::ostream& dest, __int128_t value );
