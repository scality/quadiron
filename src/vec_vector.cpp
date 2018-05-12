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
#include "vec_vector.h"

namespace nttec {
namespace vec {

void _vec_hadamard_mul_257(int n, uint32_t* x, uint32_t* y)
{
    int i;

    for (i = 0; i < n / 2; i++) {
        x[i] = (((uint64_t)x[i]) * y[i]) % 257;
    }

    for (int j = 0; j < n / 2; j++) {
        x[i + j] = (((uint64_t)x[i + j]) * y[j]) % 257;
    }
}

void _vec_hadamard_mul_65537(int n, uint32_t* x, uint32_t* y)
{
    int i;

    for (i = 0; i < n / 2; i++) {
        x[i] = (((uint64_t)x[i]) * y[i]) % 65537;
    }

    for (int j = 0; j < n / 2; j++) {
        x[i + j] = (((uint64_t)x[i + j]) * y[j]) % 65537;
    }
}

void _vec_add_257(int n, uint32_t* x, uint32_t* y)
{
    int i;

    for (i = 0; i < n / 2; i++) {
        x[i] = (x[i] + y[i]) % 257;
    }

    for (int j = 0; j < n / 2; j++) {
        x[i + j] = (x[i + j] + y[j]) % 257;
    }
}

void _vec_add_65537(int n, uint32_t* x, uint32_t* y)
{
    int i;

    for (i = 0; i < n / 2; i++) {
        x[i] = (x[i] + y[i]) % 65537;
    }

    for (int j = 0; j < n / 2; j++) {
        x[i + j] = (x[i + j] + y[j]) % 65537;
    }
}

template <>
void Vector<uint32_t>::hadamard_mul(Doubled<uint32_t>* v)
{
    assert(n == v->get_n());

    // typical butterfly operation
    if (rn->card() == 257) {
        uint32_t* a = get_mem();
        uint32_t* b = v->get_mem();
        _vec_hadamard_mul_257(n, a, b);
        return;
    } else if (rn->card() == 65537) {
        uint32_t* a = get_mem();
        uint32_t* b = v->get_mem();
        _vec_hadamard_mul_65537(n, a, b);
        return;
    }
    uint32_t* src = v->get_mem();
    int i;
    int j;
    for (i = 0; i < n / 2; i++)
        mem[i] = rn->mul(mem[i], src[i]);
    for (j = 0; i < n; i++, j++)
        mem[i] = rn->mul(mem[i], src[j]);
}

template <>
void Vector<uint32_t>::add(Doubled<uint32_t>* v)
{
    assert(n == v->get_n());

    // typical butterfly operation
    if (rn->card() == 257) {
        uint32_t* a = get_mem();
        uint32_t* b = v->get_mem();
        _vec_add_257(n, a, b);
        return;
    } else if (rn->card() == 65537) {
        uint32_t* a = get_mem();
        uint32_t* b = v->get_mem();
        _vec_add_65537(n, a, b);
        return;
    }
    uint32_t* src = v->get_mem();
    int i;
    int j;
    for (i = 0; i < n / 2; i++)
        mem[i] = rn->add(mem[i], src[i]);
    for (j = 0; i < n; i++, j++)
        mem[i] = rn->add(mem[i], src[j]);
}

void _vec_hadamard_mul_257(int n, uint64_t* x, uint64_t* y)
{
    int i;

    for (i = 0; i < n / 2; i++) {
        x[i] = (x[i] * y[i]) % 257;
    }

    for (int j = 0; j < n / 2; j++) {
        x[i + j] = (x[i + j] * y[j]) % 257;
    }
}

void _vec_hadamard_mul_65537(int n, uint64_t* x, uint64_t* y)
{
    int i;

    for (i = 0; i < n / 2; i++) {
        x[i] = (x[i] * y[i]) % 65537;
    }

    for (int j = 0; j < n / 2; j++) {
        x[i + j] = (x[i + j] * y[j]) % 65537;
    }
}

void _vec_add_257(int n, uint64_t* x, uint64_t* y)
{
    int i;

    for (i = 0; i < n / 2; i++) {
        x[i] = (x[i] + y[i]) % 257;
    }

    for (int j = 0; j < n / 2; j++) {
        x[i + j] = (x[i + j] + y[j]) % 257;
    }
}

void _vec_add_65537(int n, uint64_t* x, uint64_t* y)
{
    int i;

    for (i = 0; i < n / 2; i++) {
        x[i] = (x[i] + y[i]) % 65537;
    }

    for (int j = 0; j < n / 2; j++) {
        x[i + j] = (x[i + j] + y[j]) % 65537;
    }
}

template <>
void Vector<uint64_t>::hadamard_mul(Doubled<uint64_t>* v)
{
    assert(n == v->get_n());

    // typical butterfly operation
    if (rn->card() == 257) {
        uint64_t* a = get_mem();
        uint64_t* b = v->get_mem();
        _vec_hadamard_mul_257(n, a, b);
        return;
    } else if (rn->card() == 65537) {
        uint64_t* a = get_mem();
        uint64_t* b = v->get_mem();
        _vec_hadamard_mul_65537(n, a, b);
        return;
    }
    uint64_t* src = v->get_mem();
    int i;
    int j;
    for (i = 0; i < n / 2; i++)
        mem[i] = rn->mul(mem[i], src[i]);
    for (j = 0; i < n; i++, j++)
        mem[i] = rn->mul(mem[i], src[j]);
}

template <>
void Vector<uint64_t>::add(Doubled<uint64_t>* v)
{
    assert(n == v->get_n());

    // typical butterfly operation
    if (rn->card() == 257) {
        uint64_t* a = get_mem();
        uint64_t* b = v->get_mem();
        _vec_add_257(n, a, b);
        return;
    } else if (rn->card() == 65537) {
        uint64_t* a = get_mem();
        uint64_t* b = v->get_mem();
        _vec_add_65537(n, a, b);
        return;
    }
    uint64_t* src = v->get_mem();
    int i;
    int j;
    for (i = 0; i < n / 2; i++)
        mem[i] = rn->add(mem[i], src[i]);
    for (j = 0; i < n; i++, j++)
        mem[i] = rn->add(mem[i], src[j]);
}

} // namespace vec
} // namespace nttec
