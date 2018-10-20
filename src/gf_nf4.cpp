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

#include "gf_nf4.h"

#ifdef QUADIRON_USE_SIMD

#include "simd.h"

namespace quadiron {
namespace gf {

template <>
__uint128_t NF4<__uint128_t>::expand16(uint16_t* arr) const
{
    return simd::expand16(arr, this->n);
}

template <>
__uint128_t NF4<__uint128_t>::expand32(uint32_t* arr) const
{
    return simd::expand32(arr, this->n);
}

template <>
__uint128_t NF4<__uint128_t>::add(__uint128_t a, __uint128_t b) const
{
    __uint128_t c = simd::add(a, b);
    return c;
}

template <>
__uint128_t NF4<__uint128_t>::sub(__uint128_t a, __uint128_t b) const
{
    return simd::sub(a, b);
}

template <>
__uint128_t NF4<__uint128_t>::mul(__uint128_t a, __uint128_t b) const
{
    return simd::mul(a, b);
}

template <>
void NF4<__uint128_t>::hadamard_mul(int n, __uint128_t* x, __uint128_t* y) const
{
    simd::hadamard_mul(n, x, y);
}

template <>
GroupedValues<__uint128_t> NF4<__uint128_t>::unpack(__uint128_t a) const
{
    return simd::unpack(a);
}

template <>
void NF4<__uint128_t>::unpack(__uint128_t a, GroupedValues<__uint128_t>& b)
    const
{
    simd::unpack(a, b);
}

template <>
__uint128_t NF4<__uint128_t>::pack(__uint128_t a) const
{
    return simd::pack(a);
}

template <>
__uint128_t NF4<__uint128_t>::pack(__uint128_t a, uint32_t flag) const
{
    return simd::pack(a, flag);
}

} // namespace gf
} // namespace quadiron

#endif // #ifdef QUADIRON_USE_SIMD
