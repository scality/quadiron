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

#ifndef __QUAD_SIMD_128_H__
#define __QUAD_SIMD_128_H__

#include <x86intrin.h>

namespace quadiron {
namespace simd {

typedef __m128i VecType;

/* ============= Constant variable  ============ */

#define F4_U32 _mm_set1_epi32(65537)
#define F4_MINUS_ONE_U32 _mm_set1_epi32(65536)
#define F3_U32 _mm_set1_epi32(257)
#define F3_MINUS_ONE_U32 _mm_set1_epi32(256)

#define F3_U16 _mm_set1_epi16(257)
#define F3_MINUS_ONE_U16 _mm_set1_epi16(256)

#define ZERO (_mm_setzero_si128())
#define ONE_U16 (_mm_set1_epi16(1))
#define ONE_U32 (_mm_set1_epi32(1))

#define MASK8_LO (_mm_set1_epi16(0x80))

/* ============= Essential Operations for SSE w/ both u16 & u32 ============ */

inline VecType load_to_reg(VecType* address)
{
    return _mm_load_si128(address);
}
inline void store_to_mem(VecType* address, VecType reg)
{
    _mm_store_si128(address, reg);
}

inline VecType bit_and(VecType x, VecType y)
{
    return _mm_and_si128(x, y);
}
inline VecType bit_xor(VecType x, VecType y)
{
    return _mm_xor_si128(x, y);
}
inline uint16_t msb8_mask(VecType x)
{
    return _mm_movemask_epi8(x);
}
inline bool and_is_zero(VecType x, VecType y)
{
    return _mm_testz_si128(x, y);
}
inline bool is_zero(VecType x)
{
    return _mm_testc_si128(ZERO, x);
}

#define SHIFTR(x, imm8) (_mm_srli_si128(x, imm8))
#define BLEND8(x, y, mask) (_mm_blendv_epi8(x, y, mask))
#define BLEND16(x, y, imm8) (_mm_blend_epi16(x, y, imm8))

/* ================= Essential Operations for SSE ================= */

template <typename T>
inline VecType set_one(T val);
template <>
inline VecType set_one(uint32_t val)
{
    return _mm_set1_epi32(val);
}
template <>
inline VecType set_one(uint16_t val)
{
    return _mm_set1_epi16(val);
}

template <typename T>
inline VecType add(VecType x, VecType y);
template <>
inline VecType add<uint32_t>(VecType x, VecType y)
{
    return _mm_add_epi32(x, y);
}
template <>
inline VecType add<uint16_t>(VecType x, VecType y)
{
    return _mm_add_epi16(x, y);
}

template <typename T>
inline VecType sub(VecType x, VecType y);
template <>
inline VecType sub<uint32_t>(VecType x, VecType y)
{
    return _mm_sub_epi32(x, y);
}
template <>
inline VecType sub<uint16_t>(VecType x, VecType y)
{
    return _mm_sub_epi16(x, y);
}

template <typename T>
inline VecType mul(VecType x, VecType y);
template <>
inline VecType mul<uint32_t>(VecType x, VecType y)
{
    return _mm_mullo_epi32(x, y);
}
template <>
inline VecType mul<uint16_t>(VecType x, VecType y)
{
    return _mm_mullo_epi16(x, y);
}

template <typename T>
inline VecType compare_eq(VecType x, VecType y);
template <>
inline VecType compare_eq<uint32_t>(VecType x, VecType y)
{
    return _mm_cmpeq_epi32(x, y);
}
template <>
inline VecType compare_eq<uint16_t>(VecType x, VecType y)
{
    return _mm_cmpeq_epi16(x, y);
}

template <typename T>
inline VecType min(VecType x, VecType y);
template <>
inline VecType min<uint32_t>(VecType x, VecType y)
{
    return _mm_min_epu32(x, y);
}
template <>
inline VecType min<uint16_t>(VecType x, VecType y)
{
    return _mm_min_epu16(x, y);
}

} // namespace simd
} // namespace quadiron

#endif
