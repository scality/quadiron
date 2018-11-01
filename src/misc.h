/* -*- mode: c++ -*- */
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
#ifndef __QUAD_MISC_H__
#define __QUAD_MISC_H__

#include <cassert>
#include <cstdint>
#include <iostream>

namespace std {

std::ostream& operator<<(std::ostream& dest, __uint128_t value);
std::ostream& operator<<(std::ostream& dest, __int128_t value);

} // namespace std

namespace quadiron {

static inline uint64_t hw_timer()
{
    uint64_t x;
#if defined(__i386__)
    __asm__ volatile("rdtsc" : "=A"(x));
#elif defined(__x86_64__)
    uint64_t lo, hi;
    __asm__ volatile("rdtsc" : "=a"(lo), "=d"(hi));
    x = (hi << 32) | lo;
#elif defined(__aarch64__)
    asm volatile("mrs %0, cntvct_el0" : "=r"(x));
#elif defined(__ARM_ARCH) && (__ARM_ARCH >= 6)
    uint32_t pmccntr;
    uint32_t pmuseren;
    uint32_t pmcntenset;

    // Check that unprivileged PMCCNTR reads are alowed.
    // Source:
    // http://infocenter.arm.com/help/index.jsp?topic=/com.arm.doc.ddi0460c/Bgbhfegi.html
    asm volatile("mrc p15, 0, %0, c9, c14, 0" : "=r"(pmuseren));
    assert(pmuseren & 1);

    // Check that PMCCNTR is actually counting cycles.
    // Source:
    // http://infocenter.arm.com/help/index.jsp?topic=/com.arm.doc.ddi0460c/Bgbhfegi.html
    asm volatile("mrc p15, 0, %0, c9, c12, 1" : "=r"(pmcntenset));
    assert(pmcntenset & 0x80000000ul);

    asm volatile("mrc p15, 0, %0, c9, c13, 0" : "=r"(pmccntr));
    x = static_cast<uint64_t>(pmccntr) << 6;
#else
#error hw_timer is not defined for your CPU
#endif
    return x;
}

} // namespace quadiron

#endif
