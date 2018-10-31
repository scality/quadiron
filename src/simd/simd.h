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

#ifndef __QUAD_SIMD_SIMD_H__
#define __QUAD_SIMD_SIMD_H__

#include "core.h"
#include "property.h"
#include "vec_buffers.h"

#include "simd/allocator.h"
#include "simd/definitions.h"

namespace quadiron {

/** SIMD abstraction layer.
 *
 * This module provides a class, Register, that allows to write code that
 * leverage SIMD acceleration regardless of the supported instruction set.
 */
namespace simd {

/// Return the number of element of type T that can fit into a SIMD register.
template <typename T>
static constexpr std::size_t countof()
{
    // Without SIMD, there is no parallelism.
    if (INSTRUCTION_SET == InstructionSet::NONE) {
        return 1;
    }
    return REG_BITSZ / (sizeof(T) * CHAR_BIT);
}

} // namespace simd
} // namespace quadiron

#ifdef QUADIRON_USE_SIMD

const unsigned F4 = 65537;
const unsigned F3 = 257;

// Include essential operations that use SIMD functions
#if defined(__AVX2__)

#include "simd_256.h"

#elif defined(__SSE4_1__)

#include "simd_128.h"

#endif

// Include basic operations
#include "simd_basic.h"

// Include accelerated operations dedicated for FNT
#include "simd_fnt.h"

// Include accelerated operations dedicated for NF4
#include "simd_nf4.h"

#endif // #ifdef QUADIRON_USE_SIMD

#endif
