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

#ifdef QUADIRON_USE_SIMD
#include "simd.h"
#include "simd/simd.h"

namespace quadiron {
namespace gf {

template <>
void RingModN<uint16_t>::neg(size_t n, uint16_t* x) const
{
    simd::neg(n, x, this->_card);
}

template <>
void RingModN<uint32_t>::neg(size_t n, uint32_t* x) const
{
    simd::neg(n, x, this->_card);
}

template <>
void RingModN<uint32_t>::mul_coef_to_buf(
    uint32_t a,
    uint32_t* src,
    uint32_t* dest,
    size_t len) const
{
    simd::mul_coef_to_buf(a, src, dest, len, this->_card);
}

template <>
void RingModN<uint16_t>::butterfly_ct(
    uint16_t coef,
    uint16_t* buf1,
    uint16_t* buf2,
    size_t len) const
{
    unsigned ratio = simd::countof<uint16_t>();
    size_t vec_len = len / ratio;
    size_t last_len = len - vec_len * ratio;

    // perform vector operations
    simd::butterfly_ct(coef, buf1, buf2, vec_len, this->_card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        for (size_t i = vec_len * ratio; i < len; ++i) {
            uint16_t a = buf1[i];
            uint16_t b = mul(coef, buf2[i]);
            buf1[i] = add(a, b);
            buf2[i] = sub(a, b);
        }
    }
}

template <>
void RingModN<uint32_t>::butterfly_ct(
    uint32_t coef,
    uint32_t* buf1,
    uint32_t* buf2,
    size_t len) const
{
    unsigned ratio = simd::countof<uint32_t>();
    size_t vec_len = len / ratio;
    size_t last_len = len - vec_len * ratio;

    // perform vector operations
    simd::butterfly_ct(coef, buf1, buf2, vec_len, this->_card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        for (size_t i = vec_len * ratio; i < len; ++i) {
            uint32_t a = buf1[i];
            uint32_t b = mul(coef, buf2[i]);
            buf1[i] = add(a, b);
            buf2[i] = sub(a, b);
        }
    }
}

template <>
void RingModN<uint16_t>::butterfly_gs(
    uint16_t coef,
    uint16_t* buf1,
    uint16_t* buf2,
    size_t len) const
{
    unsigned ratio = simd::countof<uint16_t>();
    size_t vec_len = len / ratio;
    size_t last_len = len - vec_len * ratio;

    // perform vector operations
    simd::butterfly_gs(coef, buf1, buf2, vec_len, this->_card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        for (size_t i = vec_len * ratio; i < len; ++i) {
            uint16_t a = buf1[i];
            uint16_t b = buf2[i];
            uint16_t c = sub(a, b);
            buf1[i] = add(a, b);
            buf2[i] = mul(coef, c);
        }
    }
}

template <>
void RingModN<uint32_t>::butterfly_gs(
    uint32_t coef,
    uint32_t* buf1,
    uint32_t* buf2,
    size_t len) const
{
    unsigned ratio = simd::countof<uint32_t>();
    size_t vec_len = len / ratio;
    size_t last_len = len - vec_len * ratio;

    // perform vector operations
    simd::butterfly_gs(coef, buf1, buf2, vec_len, this->_card);

    // for last elements, perform as non-SIMD method
    if (last_len > 0) {
        for (size_t i = vec_len * ratio; i < len; ++i) {
            uint32_t a = buf1[i];
            uint32_t b = buf2[i];
            uint32_t c = sub(a, b);
            buf1[i] = add(a, b);
            buf2[i] = mul(coef, c);
        }
    }
}

template <>
void RingModN<uint32_t>::add_two_bufs(uint32_t* src, uint32_t* dest, size_t len)
    const
{
    simd::add_two_bufs(src, dest, len, this->_card);
}

template <>
void RingModN<uint32_t>::sub_two_bufs(
    uint32_t* bufa,
    uint32_t* bufb,
    uint32_t* res,
    size_t len) const
{
    simd::sub_two_bufs(bufa, bufb, res, len, this->_card);
}

template <>
void RingModN<uint32_t>::add_doubled(int n, uint32_t* x_u32, uint32_t* y_u32)
    const
{
    const int half = n / 2;
    uint32_t* x_next = x_u32 + half;
    simd::add_two_bufs(y_u32, x_u32, half, this->_card);
    simd::add_two_bufs(y_u32, x_next, half, this->_card);
}

template <>
void RingModN<uint16_t>::mul_coef_to_buf(
    uint16_t a,
    uint16_t* src,
    uint16_t* dest,
    size_t len) const
{
    simd::mul_coef_to_buf(a, src, dest, len, this->_card);
}

template <>
void RingModN<uint16_t>::add_two_bufs(uint16_t* src, uint16_t* dest, size_t len)
    const
{
    simd::add_two_bufs(src, dest, len, this->_card);
}

template <>
void RingModN<uint16_t>::sub_two_bufs(
    uint16_t* bufa,
    uint16_t* bufb,
    uint16_t* res,
    size_t len) const
{
    simd::sub_two_bufs(bufa, bufb, res, len, this->_card);
}

template <>
void RingModN<uint16_t>::add_doubled(int n, uint16_t* x_u16, uint16_t* y_u16)
    const
{
    const int half = n / 2;
    uint16_t* x_next = x_u16 + half;
    simd::add_two_bufs(y_u16, x_u16, half, this->_card);
    simd::add_two_bufs(y_u16, x_next, half, this->_card);
}

// TODO
// template <>
// void RingModN<uint64_t>::add_doubled(int n, uint64_t* x_u64, uint64_t* y_u64)
// const
// {
// }
// template <>
// void RingModN<__uint128_t>::add_doubled(int n, __uint128_t* x_u128,
// __uint128_t* y_u128) const
// {
// }

template <>
void RingModN<uint16_t>::hadamard_mul(int n, uint16_t* x_u16, uint16_t* y_u16)
    const
{
    simd::mul_two_bufs(y_u16, x_u16, n, this->_card);
}

template <>
void RingModN<uint32_t>::hadamard_mul(int n, uint32_t* x_u32, uint32_t* y_u32)
    const
{
    simd::mul_two_bufs(y_u32, x_u32, n, this->_card);
}

template <>
void RingModN<uint16_t>::hadamard_mul_doubled(
    int n,
    uint16_t* x_u16,
    uint16_t* y_u16) const
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
    uint32_t* y_u32) const
{
    const int half = n / 2;
    uint32_t* x_next = x_u32 + half;
    simd::mul_two_bufs(y_u32, x_u32, half, this->_card);
    simd::mul_two_bufs(y_u32, x_next, half, this->_card);
}

// TODO
// template <>
// void RingModN<uint64_t>::hadamard_mul_doubled(int n, uint64_t* x_u64,
// uint64_t* y_u64) const
// {
// }
// template <>
// void RingModN<__uint128_t>::hadamard_mul_doubled(int n, __uint128_t* x_u128,
// __uint128_t* y_u128) const
// {
// }

} // namespace gf
} // namespace quadiron

#endif // #ifdef QUADIRON_USE_SIMD
