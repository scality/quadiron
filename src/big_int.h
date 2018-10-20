/* -*- mode: c++ -*- */
/*
 * Copyright 2017-2018 Scality
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
#ifndef __QUAD_BIG_INT_H__
#define __QUAD_BIG_INT_H__

#include <iostream>

namespace std {

std::ostream& operator<<(std::ostream& dest, __uint128_t value);
std::ostream& operator<<(std::ostream& dest, __int128_t value);

} // namespace std

namespace quadiron {

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
        return lo < x;
    }

    bool operator<(int x)
    {
        return lo < static_cast<__uint128_t>(x);
    }

    bool operator!=(__uint128_t x)
    {
        return lo != x;
    }

    bool operator!=(int x)
    {
        return lo != static_cast<__uint128_t>(x);
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

} // namespace quadiron

#endif
