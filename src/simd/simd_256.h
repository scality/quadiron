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

#ifndef __QUAD_SIMD_256_H__
#define __QUAD_SIMD_256_H__

/* GCC doesn't include the split store intrinsics so define them here. */
#if defined(__GNUC__) && !defined(__clang__)

static inline void __attribute__((__always_inline__))
_mm256_storeu2_m128i(__m128i* const hi, __m128i* const lo, const __m256i a)
{
    _mm_storeu_si128(lo, _mm256_castsi256_si128(a));
    _mm_storeu_si128(hi, _mm256_extracti128_si256(a, 1));
}

#endif /* defined(__GNUC__) */

namespace quadiron {
namespace simd {

typedef __m256i VecType;
typedef __m128i HalfVecType;

/* ============= Constant variable  ============ */

#define F4_U32 _mm256_set1_epi32(65537)
#define F4_MINUS_ONE_U32 _mm256_set1_epi32(65536)
#define F3_U32 _mm256_set1_epi32(257)
#define F3_MINUS_ONE_U32 _mm256_set1_epi32(256)

#define F3_U16 _mm256_set1_epi16(257)
#define F3_MINUS_ONE_U16 _mm256_set1_epi16(256)

#define ZERO (_mm256_setzero_si256())
#define ONE_U16 (_mm256_set1_epi16(1))
#define ONE_U32 (_mm256_set1_epi32(1))

#define MASK8_LO (_mm256_set1_epi16(0x80))

/* ============= Essential Operations for AVX2 w/ both u16 & u32 ============ */

inline VecType LoadToReg(VecType* address)
{
    return _mm256_load_si256(address);
}
inline void StoreToMem(VecType* address, VecType reg)
{
    _mm256_store_si256(address, reg);
}

inline VecType And(VecType x, VecType y)
{
    return _mm256_and_si256(x, y);
}
inline VecType Xor(VecType x, VecType y)
{
    return _mm256_xor_si256(x, y);
}
inline uint32_t Msb8Mask(VecType x)
{
    return _mm256_movemask_epi8(x);
}
inline uint32_t AndIsZero(VecType x, VecType y)
{
    return _mm256_testz_si256(x, y);
}
inline int IsZero(VecType x)
{
    return _mm256_testc_si256(ZERO, x);
}

#define SHIFTR(x, imm8) (_mm256_srli_si256(x, imm8))
#define BLEND8(x, y, mask) (_mm256_blendv_epi8(x, y, mask))
#define BLEND16(x, y, imm8) (_mm256_blend_epi16(x, y, imm8))

/* ================= Essential Operations for AVX2 ================= */

template <typename T>
inline VecType SetOne(T val);
template <>
inline VecType SetOne(uint32_t val)
{
    return _mm256_set1_epi32(val);
}
template <>
inline VecType SetOne(uint16_t val)
{
    return _mm256_set1_epi16(val);
}

template <typename T>
inline VecType Add(VecType x, VecType y);
template <>
inline VecType Add<uint32_t>(VecType x, VecType y)
{
    return _mm256_add_epi32(x, y);
}
template <>
inline VecType Add<uint16_t>(VecType x, VecType y)
{
    return _mm256_add_epi16(x, y);
}

template <typename T>
inline VecType Sub(VecType x, VecType y);
template <>
inline VecType Sub<uint32_t>(VecType x, VecType y)
{
    return _mm256_sub_epi32(x, y);
}
template <>
inline VecType Sub<uint16_t>(VecType x, VecType y)
{
    return _mm256_sub_epi16(x, y);
}

template <typename T>
inline VecType Mul(VecType x, VecType y);
template <>
inline VecType Mul<uint32_t>(VecType x, VecType y)
{
    return _mm256_mullo_epi32(x, y);
}
template <>
inline VecType Mul<uint16_t>(VecType x, VecType y)
{
    return _mm256_mullo_epi16(x, y);
}

template <typename T>
inline VecType CompareEq(VecType x, VecType y);
template <>
inline VecType CompareEq<uint32_t>(VecType x, VecType y)
{
    return _mm256_cmpeq_epi32(x, y);
}
template <>
inline VecType CompareEq<uint16_t>(VecType x, VecType y)
{
    return _mm256_cmpeq_epi16(x, y);
}

template <typename T>
inline VecType Min(VecType x, VecType y);
template <>
inline VecType Min<uint32_t>(VecType x, VecType y)
{
    return _mm256_min_epu32(x, y);
}
template <>
inline VecType Min<uint16_t>(VecType x, VecType y)
{
    return _mm256_min_epu16(x, y);
}

} // namespace simd
} // namespace quadiron

#endif
