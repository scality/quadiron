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

#include <x86intrin.h>

typedef __m256i m256i;

/* GCC doesn't include the split store intrinsics so define them here. */
#if defined(__GNUC__) && !defined(__clang__)

static inline void __attribute__((__always_inline__))
_mm256_storeu2_m128i(__m128i* const hi, __m128i* const lo, const __m256i a)
{
    _mm_storeu_si128(lo, _mm256_castsi256_si128(a));
    _mm_storeu_si128(hi, _mm256_extracti128_si256(a, 1));
}

#endif /* defined(__GNUC__) */

    // Following functions are used for AVX2 w/ both of u16 & u32

#define ZERO (_mm256_setzero_si256())
#define ONE (_mm256_set1_epi32(1))

inline m256i LOAD(m256i* address)
{
    return _mm256_load_si256(address);
}
inline void STORE(m256i* address, m256i reg)
{
    _mm256_store_si256(address, reg);
}

inline m256i AND(m256i x, m256i y)
{
    return _mm256_and_si256(x, y);
}
inline m256i XOR(m256i x, m256i y)
{
    return _mm256_xor_si256(x, y);
}
inline m256i SHIFTR_1(m256i x)
{
    return _mm256_srli_si256(x, 1);
}
inline m256i SHIFTR_2(m256i x)
{
    return _mm256_srli_si256(x, 2);
}
inline uint32_t MVMSK8(m256i x)
{
    return _mm256_movemask_epi8(x);
}
inline uint32_t TESTZ(m256i x, m256i y)
{
    return _mm256_testz_si256(x, y);
}

#include "simd_256_u16.h"
#include "simd_256_u32.h"

#endif
