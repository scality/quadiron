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

typedef __m128i m128i;

// Disable `cert-err58-cpp` on these: AFAIK they cannot throw.
// (probably a false positive present in Clang 5 and fixed in Clang 6).
const m128i F4_m128i = _mm_set1_epi32(65537);       // NOLINT(cert-err58-cpp)
const m128i F4minus1_m128i = _mm_set1_epi32(65536); // NOLINT(cert-err58-cpp)
const m128i F4minus2_m128i = _mm_set1_epi32(65535); // NOLINT(cert-err58-cpp)
const m128i F3_m128i = _mm_set1_epi32(257);         // NOLINT(cert-err58-cpp)
const m128i F3minus1_m128i = _mm_set1_epi32(256);   // NOLINT(cert-err58-cpp)
const m128i F3minus2_m128i = _mm_set1_epi32(255);   // NOLINT(cert-err58-cpp)

const m128i F3_m128i_u16 = _mm_set1_epi16(257);       // NOLINT(cert-err58-cpp)
const m128i F3minus1_m128i_u16 = _mm_set1_epi16(256); // NOLINT(cert-err58-cpp)

#include "simd_128_u16.h"
#include "simd_128_u32.h"

#endif
