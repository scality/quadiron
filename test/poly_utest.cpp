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
#include "nttec.h"

void poly_utest1()
{
    nttec::gf::Prime<uint32_t> gfp(1000);

    nttec::Polynomial<uint32_t> p0(&gfp);
    assert(0 == p0.degree());
}

/**
 * Examples taken from https://www.doc.ic.ac.uk/~mrh/330tutor/ch04s02.html
 */
void poly_utest2()
{
    std::cout << "poly_utest2\n";

    nttec::gf::Prime<uint32_t> gfp(11);
    nttec::Polynomial<uint32_t> p1(&gfp);
    p1.set(5, 1);
    p1.set(3, 3);
    p1.set(0, 4);
    // p1.dump();
    nttec::Polynomial<uint32_t> p2(&gfp);
    p2.set(6, 6);
    p2.set(3, 4);
    // p2.dump();
    nttec::Polynomial<uint32_t> p3(&gfp);
    p1._add(&p3, &p1, &p2);
    // p3.dump();
    assert(p3.degree() == 6);
    assert(p3.get(6) == 6);
    assert(p3.get(5) == 1);
    assert(p3.get(3) == 7);
    assert(p3.get(0) == 4);
}

void poly_utest3()
{
    std::cout << "poly_utest3\n";

    nttec::gf::Prime<uint32_t> gfp(11);
    nttec::Polynomial<uint32_t> p1(&gfp);
    p1.set(5, 1);
    p1.set(3, 3);
    p1.set(0, 4);
    // p1.dump();
    nttec::Polynomial<uint32_t> p2(&gfp);
    p2.set(6, 6);
    p2.set(3, 4);
    // p2.dump();
    nttec::Polynomial<uint32_t> p3(&gfp);
    p1._sub(&p3, &p1, &p2);
    assert(p3.degree() == 6);
    assert(p3.get(6) == 5);
    assert(p3.get(5) == 1);
    assert(p3.get(3) == 10);
    assert(p3.get(0) == 4);
}

void poly_utest4()
{
    std::cout << "poly_utest4\n";

    nttec::gf::Prime<uint32_t> gfp(11);
    nttec::Polynomial<uint32_t> p1(&gfp);
    p1.set(5, 1);
    p1.set(3, 3);
    p1.set(0, 4);
    // p1.dump();
    nttec::Polynomial<uint32_t> p2(&gfp);
    p2.set(6, 6);
    p2.set(3, 4);
    // p2.dump();
    nttec::Polynomial<uint32_t> p3(&gfp);
    p1._mul(&p3, &p1, &p2);
    assert(p3.degree() == 11);
    assert(p3.get(11) == 6);
    assert(p3.get(9) == 7);
    assert(p3.get(8) == 4);
    assert(p3.get(6) == 3);
    assert(p3.get(3) == 5);
}

void poly_utest5()
{
    std::cout << "poly_utest5\n";
    nttec::gf::Prime<uint32_t> gfp(11);
    nttec::Polynomial<uint32_t> p1(&gfp);
    p1.set(6, 3);
    p1.set(4, 7);
    p1.set(3, 4);
    p1.set(0, 5);
    // p1.dump();
    nttec::Polynomial<uint32_t> p2(&gfp);
    p2.set(4, 1);
    p2.set(3, 3);
    p2.set(0, 4);
    // p2.dump();
    nttec::Polynomial<uint32_t> p3(&gfp);
    nttec::Polynomial<uint32_t> p4(&gfp);
    p1._div(&p3, &p4, &p1, &p2);
    // p3.dump();
    // p4.dump();
    assert(p3.degree() == 2);
    assert(p3.get(2) == 3);
    assert(p3.get(1) == 2);
    assert(p3.get(0) == 1);
    assert(p4.degree() == 3);
    assert(p4.get(3) == 1);
    assert(p4.get(2) == 10);
    assert(p4.get(1) == 3);
    assert(p4.get(0) == 1);
}

void poly_utest6()
{
    std::cout << "poly_utest6\n";
    nttec::gf::Prime<uint32_t> gfp(11);
    nttec::Polynomial<uint32_t> p1(&gfp);
    p1.set(4, 4);
    p1.set(1, 1);
    p1.set(0, 8);
    nttec::Polynomial<uint32_t> p2(&gfp);
    p1._derivative(&p2, &p1);
    assert(p2.degree() == 3);
    assert(p2.get(3) == 5);
    assert(p2.get(0) == 1);
}

void poly_utest7()
{
    std::cout << "poly_utest7\n";
    nttec::gf::Prime<uint32_t> gfp(11);
    nttec::Polynomial<uint32_t> p1(&gfp);
    p1.set(4, 4);
    p1.set(1, 1);
    p1.set(0, 8);
    uint32_t result = p1.eval(3);
    assert(result == 5);
}

void poly_utest8()
{
    std::cout << "poly_utest8\n";
    nttec::gf::Prime<uint32_t> gfp(3);
    nttec::Polynomial<uint32_t> p1(&gfp);
    p1.set(4, 1);
    // p1.dump();
    nttec::Polynomial<uint32_t> p2(&gfp);
    p2.set(3, 1);
    p2.set(2, 2);
    p2.set(0, 1);
    // p2.dump();
    nttec::Polynomial<uint32_t> p3(&gfp);
    nttec::Polynomial<uint32_t> p4(&gfp);
    p1._div(&p3, &p4, &p1, &p2);
    // p3.dump();
    // p4.dump();
    assert(p3.degree() == 1);
    assert(p3.get(1) == 1);
    assert(p3.get(0) == 1);
    assert(p4.degree() == 2);
    assert(p4.get(2) == 1);
    assert(p4.get(1) == 2);
    assert(p4.get(0) == 2);
}

void poly_utest()
{
    std::cout << "poly_utest\n";

    poly_utest1();
    poly_utest2();
    poly_utest3();
    poly_utest4();
    poly_utest5();
    poly_utest6();
    poly_utest7();
    poly_utest8();
}
