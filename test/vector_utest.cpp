/*
 * Copyright 2017-2018 the NTTEC authors
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

template <typename T>
class VectorUtest {
  public:
    void vector_utest1()
    {
        std::cout << "vector_utest1\n";

        nttec::gf::Prime<T> gfp(65537);
        nttec::vec::Vector<T> vec1(&gfp, 16);
        nttec::vec::Vector<T> vec2(&gfp, 8);
        nttec::vec::Doubled<T> v2vec2(&vec2);

        vec1.set(0, 1);
        vec1.set(1, 64);
        vec1.set(2, 4096);
        vec1.set(3, 65533);
        vec1.set(4, 65281);
        vec1.set(5, 49153);
        vec1.set(6, 16);
        vec1.set(7, 1024);
        vec1.set(8, 65536);
        vec1.set(9, 65473);
        vec1.set(10, 61441);
        vec1.set(11, 4);
        vec1.set(12, 256);
        vec1.set(13, 16384);
        vec1.set(14, 65521);
        vec1.set(15, 64513);

        vec2.set(0, 6505);
        vec2.set(1, 23324);
        vec2.set(2, 7736);
        vec2.set(3, 43678);
        vec2.set(4, 27296);
        vec2.set(5, 50697);
        vec2.set(6, 59274);
        vec2.set(7, 48649);

        vec1.hadamard_mul(&v2vec2);

        assert(vec1.get(0) == 6505);
        assert(vec1.get(1) == 50922);
        assert(vec1.get(2) == 32285);
        assert(vec1.get(3) == 21899);
        assert(vec1.get(4) == 24683);
        assert(vec1.get(5) == 61827);
        assert(vec1.get(6) == 30866);
        assert(vec1.get(7) == 8456);
        assert(vec1.get(8) == 59032);
        assert(vec1.get(9) == 14615);
        assert(vec1.get(10) == 33252);
        assert(vec1.get(11) == 43638);
        assert(vec1.get(12) == 40854);
        assert(vec1.get(13) == 3710);
        assert(vec1.get(14) == 34671);
        assert(vec1.get(15) == 57081);
    }

    void vector_utest2()
    {
        std::cout << "vector_utest2\n";

        nttec::gf::Prime<T> gfp(65537);
        nttec::vec::Vector<T> vec1(&gfp, 8);
        nttec::vec::Vector<T> vec2(&gfp, 4);
        nttec::vec::Doubled<T> v2vec2(&vec2);

        vec1.set(0, 5459);
        vec1.set(1, 11947);
        vec1.set(2, 44310);
        vec1.set(3, 21807);
        vec1.set(4, 60078);
        vec1.set(5, 53590);
        vec1.set(6, 21227);
        vec1.set(7, 43730);

        vec2.set(0, 39466);
        vec2.set(1, 40329);
        vec2.set(2, 16012);
        vec2.set(3, 15149);

        vec1.add(&v2vec2);

        assert(vec1.get(0) == 44925);
        assert(vec1.get(1) == 52276);
        assert(vec1.get(2) == 60322);
        assert(vec1.get(3) == 36956);
        assert(vec1.get(4) == 34007);
        assert(vec1.get(5) == 28382);
        assert(vec1.get(6) == 37239);
        assert(vec1.get(7) == 58879);
    }

    void vector_utest3()
    {
        std::cout << "vector_utest3\n";

        nttec::gf::Prime<T> gfp(65537);
        int len = 20;
        nttec::vec::Vector<T> base_vec(&gfp, len);
        for (int i = 0; i < len; i++)
            base_vec.set(i, gfp.weak_rand());
        int len1 = 7;
        int len2 = 5;
        int offset1 = 2;
        int offset2 = 3;
        int len3 = len2;
        if (offset2 + len2 > len1)
            len3 = len1 - offset2;
        // vmvec1 = base_vec[offset1, .., offset1 + len1 - 1]
        nttec::vec::Slice<T> vmvec1(&base_vec, len1, offset1);
        // vmvec2 = vmvec1[offset2, .., min(offset2 + len2 - 1, len1 - 1)]
        nttec::vec::Slice<T> vmvec2(&vmvec1, len2, offset2);
        // vmvec3 = base_vec[offset1 + offset2, .., offset1 + len1 - 1]
        nttec::vec::Slice<T> vmvec3(&base_vec, len3, offset1 + offset2);
        assert(vmvec3.eq(&vmvec2));
    }

    void vector_utest()
    {
        std::cout << "vector_utest\n";

        vector_utest1();
        vector_utest2();
        vector_utest3();
    }
};

void vector_utest()
{
    VectorUtest<uint32_t> vecutest_uint32;
    vecutest_uint32.vector_utest();
    VectorUtest<uint64_t> vecutest_uint64;
    vecutest_uint64.vector_utest();
}
