/* -*- mode: c++ -*- */
#ifndef __NTTEC_BIG_INT_H__
#define __NTTEC_BIG_INT_H__

#include <iostream>

namespace std {

std::ostream& operator<<(std::ostream& dest, __uint128_t value);
std::ostream& operator<<(std::ostream& dest, __int128_t value);

} // namespace std

namespace nttec {

struct UInt256 {
    __uint128_t lo;
    __uint128_t hi;
    UInt256() : lo(), hi() {}
    UInt256(__uint128_t x) : lo(x), hi(0) {}

    UInt256 operator*(__uint128_t x)
    {
        return lo * x;
    }
    UInt256 operator%(__uint128_t x)
    {
        return lo % x;
    }
    operator __uint128_t()
    {
        return lo;
    }
};

struct Int256 {
    __uint128_t lo;
    __uint128_t hi;
    Int256() : lo(), hi() {}
    Int256(__uint128_t x) : lo(x), hi(0) {}

    bool operator<(__uint128_t x)
    {
        return (lo < x);
    }
    bool operator<(int x)
    {
        return (lo < (__uint128_t)x);
    }
    bool operator!=(__uint128_t x)
    {
        return (lo != x);
    }
    bool operator!=(int x)
    {
        return (lo != (__uint128_t)x);
    }
    Int256 operator+(__uint128_t x)
    {
        return lo + x;
    }
    Int256 operator-(Int256 x)
    {
        return lo - x.lo;
    }
    Int256 operator*(Int256 x)
    {
        return lo * x.lo;
    }
    Int256 operator/(Int256 x)
    {
        return lo / x.lo;
    }
    Int256 operator/(int x)
    {
        return lo / x;
    }
    operator __uint128_t()
    {
        return lo;
    }
};

} // namespace nttec

#endif
