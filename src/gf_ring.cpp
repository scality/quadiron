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

} // namespace gf
} // namespace quadiron

#endif // #ifdef QUADIRON_USE_SIMD
