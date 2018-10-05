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
typedef uint32_t MaskIntType;

#define F4_u32 _mm_set1_epi32(65537)
#define F4m1_u32 _mm_set1_epi32(65536)
#define F3_u32 _mm_set1_epi32(257)
#define F3m1_u32 _mm_set1_epi32(256)

#define F3_u16 _mm_set1_epi16(257)
#define F3m1_u16 _mm_set1_epi16(256)

#define CARD(q) (EITHER(q == F3, F3_u32, F4_u32))
#define CARD_M_1(q) (EITHER(q == F3, F3m1_u32, F4m1_u32))

/* ============= Essential Operations for AVX2 w/ both u16 & u32 ============ */

#define ZERO (_mm_setzero_si128())
#define ONE16 (_mm_set1_epi16(1))
#define ONE32 (_mm_set1_epi32(1))

inline VecType LOAD(VecType* address)
{
    return _mm_load_si128(address);
}
inline void STORE(VecType* address, VecType reg)
{
    _mm_store_si128(address, reg);
}

inline VecType AND(VecType x, VecType y)
{
    return _mm_and_si128(x, y);
}
inline VecType XOR(VecType x, VecType y)
{
    return _mm_xor_si128(x, y);
}
inline VecType SHIFTR_1(VecType x)
{
    return _mm_srli_si128(x, 1);
}
inline VecType SHIFTR_2(VecType x)
{
    return _mm_srli_si128(x, 2);
}
inline uint16_t MVMSK8(VecType x)
{
    return _mm_movemask_epi8(x);
}
inline uint16_t TESTZ(VecType x, VecType y)
{
    return _mm_testz_si128(x, y);
}

/* ================= Essential Operations for AVX2 w/ u32 ================= */

inline VecType SET1(uint32_t val)
{
    return _mm_set1_epi32(val);
}
inline VecType ADD32(VecType x, VecType y)
{
    return _mm_add_epi32(x, y);
}
inline VecType SUB32(VecType x, VecType y)
{
    return _mm_sub_epi32(x, y);
}
inline VecType MUL32(VecType x, VecType y)
{
    return _mm_mullo_epi32(x, y);
}

inline VecType CMPEQ32(VecType x, VecType y)
{
    return _mm_cmpeq_epi32(x, y);
}
inline VecType CMPGT32(VecType x, VecType y)
{
    return _mm_cmpgt_epi32(x, y);
}
inline VecType MINU32(VecType x, VecType y)
{
    return _mm_min_epu32(x, y);
}
#define MASK8_LO (_mm_set1_epi16(0x80))
#define BLEND8(x, y, mask) (_mm_blendv_epi8(x, y, mask))
#define BLEND16(x, y, imm8) (_mm_blend_epi16(x, y, imm8))

/* ================= Essential Operations for AVX2 w/ u16 ================= */

inline VecType SET1(uint16_t val)
{
    return _mm_set1_epi16(val);
}
inline VecType ADD16(VecType x, VecType y)
{
    return _mm_add_epi16(x, y);
}
inline VecType SUB16(VecType x, VecType y)
{
    return _mm_sub_epi16(x, y);
}
inline VecType MUL16(VecType x, VecType y)
{
    return _mm_mullo_epi16(x, y);
}

inline VecType CMPEQ16(VecType x, VecType y)
{
    return _mm_cmpeq_epi16(x, y);
}
inline VecType CMPGT16(VecType x, VecType y)
{
    return _mm_cmpgt_epi16(x, y);
}
inline VecType MINU16(VecType x, VecType y)
{
    return _mm_min_epu16(x, y);
}

} // namespace simd
} // namespace quadiron

#endif
