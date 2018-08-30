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

#ifndef __QUAD_SIMD_H__
#define __QUAD_SIMD_H__

#ifdef QUADIRON_USE_SIMD

const unsigned F4 = 65537;
const unsigned F3 = 257;

typedef uint8_t aint8 __attribute__((aligned(ALIGN_SIZE)));
typedef uint16_t aint16 __attribute__((aligned(ALIGN_SIZE)));
typedef uint32_t aint32 __attribute__((aligned(ALIGN_SIZE)));
typedef uint64_t aint64 __attribute__((aligned(ALIGN_SIZE)));
typedef __uint128_t aint128 __attribute__((aligned(ALIGN_SIZE)));

namespace quadiron {
/** The namespace simd contains functions for GF-NF4 that are accelerated by
 *  using SIMD operations over 128bits
 *
 *  It supports operations on 32-bit numbers
 */
namespace simd {

// Vectorized operations are implemented in appropriated headers simd*.h

} // namespace simd
} // namespace quadiron

#ifdef QUADIRON_USE_SSE4
#include "simd_128.h"
#elif defined QUADIRON_USE_AVX2
#include "simd_256.h"
#endif

#include "simd_nf4.h"

#endif // #ifdef QUADIRON_USE_SIMD

#endif
