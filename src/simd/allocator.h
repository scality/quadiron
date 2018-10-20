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

/** @file allocator.h
 *
 * Provide a custom memory allocator for SIMD.
 *
 * This allocator always returns memory that is suitably aligned to be loaded
 * efficiently by a Register object.
 */

#ifndef __QUAD_SIMD_SIMD_ALLOCATOR_H__
#define __QUAD_SIMD_SIMD_ALLOCATOR_H__

#include <cassert>
#include <cstdint>
#include <limits>

#include "simd/definitions.h"

namespace quadiron {
namespace simd {

/// Check that the given address is properly aligned.
template <typename T>
inline bool addr_is_aligned(const T* addr)
{
    // Without SIMD, there is no specific alignment constraint.
    if (INSTRUCTION_SET == InstructionSet::NONE) {
        return true;
    }
    const std::uintptr_t address = reinterpret_cast<std::uintptr_t>(addr);
    return (address & (ALIGNMENT - 1)) == 0;
}

/** Custom allocator to take advantage of SIMD processing.
 *
 * This allocator always return memory that is suitably aligned for the current
 * SIMD instruction set. Thanks to this property, you can safely use the aligned
 * load from the Register class in order to increase performance.
 */
template <typename T>
class AlignedAllocator {
  public:
    using value_type = T;

    AlignedAllocator() noexcept {}
    // No state, => nothing to copy.
    template <class U>
    AlignedAllocator(AlignedAllocator<U> const& /* other */) noexcept
    {
    }

    value_type* allocate(std::size_t count)
    {
        // Guard against overflow!
        if (count > max_size()) {
            throw std::bad_alloc();
        }

        // No SIMD: default allocator is good enough!
        if (INSTRUCTION_SET == InstructionSet::NONE) {
            return static_cast<value_type*>(
                ::operator new(count * sizeof(value_type)));
        }

        // Overallocate just enough to have room for alignment adjustment.
        const std::size_t size = count * sizeof(value_type) + ALIGNMENT;
        unsigned char* ptr = static_cast<unsigned char*>(::operator new(size));

        // Align the allocated memory.
        const std::uintptr_t address = reinterpret_cast<std::uintptr_t>(ptr);
        const unsigned offset = ALIGNMENT - (address % ALIGNMENT);
        assert(offset >= 1); // We need a byte to store the offset itself.
        unsigned char* aligned_ptr = ptr + offset;

        // Store the offset just before the aligned memory.
        assert(offset <= std::numeric_limits<unsigned char>::max());
        *(aligned_ptr - 1) = static_cast<unsigned char>(offset);

        // Return the aligned pointer.
        //
        // Clang analyser think that we leak `ptr`, whereas we can re-compute
        // it from `aligned_ptr` and free it in `deallocate`.
        // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
        return reinterpret_cast<value_type*>(aligned_ptr);
    }

    void deallocate(value_type* ptr, std::size_t /* count */) noexcept
    {
        // No SIMD: default allocator is good enough!
        if (INSTRUCTION_SET == InstructionSet::NONE) {
            ::operator delete(ptr);
            return;
        }

        if (ptr == nullptr) {
            return;
        }
        // Respect strict aliasing rules: read through a character type.
        unsigned char* raw = reinterpret_cast<unsigned char*>(ptr);
        // Get the alignment offset stored just before the aligned pointer.
        const unsigned offset = *(raw - 1);
        ::operator delete(raw - offset);
    }

    std::size_t max_size() const noexcept
    {
        const std::size_t max_size = std::numeric_limits<std::size_t>::max();
        return (max_size - ALIGNMENT) / sizeof(value_type);
    }

    // Our allocator is stateless.
    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;
};

template <class T, class U>
bool operator==(AlignedAllocator<T> const&, AlignedAllocator<U> const&) noexcept
{
    // Our allocator is stateless: Any instance of our allocator can deallocate
    // the memory from another instance.
    return true;
}

template <class T, class U>
bool operator!=(
    AlignedAllocator<T> const& x,
    AlignedAllocator<U> const& y) noexcept
{
    return !(x == y);
}

} // namespace simd
} // namespace quadiron

#endif
