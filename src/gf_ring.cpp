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

#include "gf_ring.h"

#ifdef NTTEC_USE_SIMD
#include "simd.h"

namespace nttec {
namespace gf {

template <>
void RingModN<uint32_t>::mul_coef_to_buf(
    uint32_t a,
    uint32_t* src,
    uint32_t* dest,
    size_t len)
{
    simd::mul_coef_to_buf(a, src, dest, len, this->_card);
}

template <>
void RingModN<uint32_t>::add_two_bufs(uint32_t* src, uint32_t* dest, size_t len)
{
    simd::add_two_bufs(src, dest, len, this->_card);
}

template <>
void RingModN<uint32_t>::sub_two_bufs(
    uint32_t* bufa,
    uint32_t* bufb,
    uint32_t* res,
    size_t len)
{
    simd::sub_two_bufs(bufa, bufb, res, len, this->_card);
}

template <>
void RingModN<uint16_t>::add_doubled(int n, uint16_t* x_u16, uint16_t* y_u16)
{
    const int half = n / 2;
    uint16_t* x_next = x_u16 + half;
    simd::add_two_bufs(y_u16, x_u16, half, this->_card);
    simd::add_two_bufs(y_u16, x_next, half, this->_card);
}

template <>
void RingModN<uint32_t>::add_doubled(int n, uint32_t* x_u32, uint32_t* y_u32)
{
    const int half = n / 2;
    uint32_t* x_next = x_u32 + half;
    simd::add_two_bufs(y_u32, x_u32, half, this->_card);
    simd::add_two_bufs(y_u32, x_next, half, this->_card);
}

// TODO
// template <>
// void RingModN<uint64_t>::add_doubled(int n, uint64_t* x_u64, uint64_t* y_u64)
// {
// }
// template <>
// void RingModN<__uint128_t>::add_doubled(int n, __uint128_t* x_u128,
// __uint128_t* y_u128)
// {
// }

template <>
void RingModN<uint16_t>::hadamard_mul(int n, uint16_t* x_u16, uint16_t* y_u16)
{
    simd::mul_two_bufs(y_u16, x_u16, n, this->_card);
}

template <>
void RingModN<uint32_t>::hadamard_mul(int n, uint32_t* x_u32, uint32_t* y_u32)
{
    simd::mul_two_bufs(y_u32, x_u32, n, this->_card);
}

template <>
void RingModN<uint16_t>::hadamard_mul_doubled(
    int n,
    uint16_t* x_u16,
    uint16_t* y_u16)
{
    const int half = n / 2;
    uint16_t* x_next = x_u16 + half;
    simd::mul_two_bufs(y_u16, x_u16, half, this->_card);
    simd::mul_two_bufs(y_u16, x_next, half, this->_card);
}

template <>
void RingModN<uint32_t>::hadamard_mul_doubled(
    int n,
    uint32_t* x_u32,
    uint32_t* y_u32)
{
    const int half = n / 2;
    uint32_t* x_next = x_u32 + half;
    simd::mul_two_bufs(y_u32, x_u32, half, this->_card);
    simd::mul_two_bufs(y_u32, x_next, half, this->_card);
}

// TODO
// template <>
// void RingModN<uint64_t>::hadamard_mul_doubled(int n, uint64_t* x_u64,
// uint64_t* y_u64)
// {
// }
// template <>
// void RingModN<__uint128_t>::hadamard_mul_doubled(int n, __uint128_t* x_u128,
// __uint128_t* y_u128)
// {
// }

} // namespace gf
} // namespace nttec

#endif // #ifdef NTTEC_USE_SIMD
