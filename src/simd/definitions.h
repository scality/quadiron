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

/** @file definitions.h
 *
 * Provide includes and constants that are platform-specific.
 */

#ifndef __QUAD_SIMD_SIMD_DEFINITIONS_H__
#define __QUAD_SIMD_SIMD_DEFINITIONS_H__

#include <climits>
#include <cstdint>

// SIMD-specific includes {{{

#ifdef QUADIRON_USE_SIMD

// SIMD for x86 and x86_64
#if defined(__i386__) || defined(__x86_64__)
#include <x86intrin.h>
#endif

// SIMD for ARM
#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif

#endif

// }}}

namespace quadiron {
namespace simd {

/// Supported instruction set.
enum class InstructionSet {
    NONE, ///< No SIMD instruction (fallback).
    SSE,  ///< SSE4.1
    AVX,  ///< AVX2
};

// Definitions for Intel AVX-256 {{{

// We required AVX2 because we relies on some instructions (such as
// `_mm256_add_epi16` and others) that aren't available in the first version of
// AVX.
#if defined(__AVX__) && defined(__AVX2__)

using RegisterType = __m256;
using MaskType = __m256i;

static constexpr InstructionSet INSTRUCTION_SET = InstructionSet::AVX;

// }}}
// Definitions for Intel SSE {{{

#elif defined(__SSE__) && defined(__SSE4_1__)

using RegisterType = __m128;
using MaskType = __m128i;

static constexpr InstructionSet INSTRUCTION_SET = InstructionSet::SSE;

// }}}
// Definitions for scalar fallback {{{

#else

// Use The Widest Integer Type Available.
#if __LP64__ == 1
using RegisterType = uint64_t;
using MaskType = uint64_t;
#else
using RegisterType = uint32_t;
using MaskType = uint32_t;
#endif

static constexpr InstructionSet INSTRUCTION_SET = InstructionSet::NONE;

#endif

// }}}
// Portable definitions {{{

/// Alignment constraint (in bytes).
static constexpr std::size_t ALIGNMENT = alignof(RegisterType);

/// Register size (in bits).
static constexpr std::size_t REG_BITSZ = sizeof(RegisterType) * CHAR_BIT;

// }}}

} // namespace simd
} // namespace quadiron

#endif
