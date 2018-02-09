#include "nttec.h"

template <typename T>
class VECUtest {
  public:
    void vec_utest1()
    {
        std::cout << "vec_utest1\n";

        GFP<T> gfp(65537);
        Vec<T> vec1(&gfp, 16);
        Vec<T> vec2(&gfp, 8);
        V2Vec<T> v2vec2(&vec2);

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

    void vec_utest2()
    {
        std::cout << "vec_utest2\n";

        GFP<T> gfp(65537);
        Vec<T> vec1(&gfp, 8);
        Vec<T> vec2(&gfp, 4);
        V2Vec<T> v2vec2(&vec2);

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

    void vec_utest3()
    {
        std::cout << "vec_utest3\n";

        GFP<T> gfp(65537);
        int len = 20;
        Vec<T> base_vec(&gfp, len);
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
        VmVec<T> vmvec1(&base_vec, len1, offset1);
        // vmvec2 = vmvec1[offset2, .., min(offset2 + len2 - 1, len1 - 1)]
        VmVec<T> vmvec2(&vmvec1, len2, offset2);
        // vmvec3 = base_vec[offset1 + offset2, .., offset1 + len1 - 1]
        VmVec<T> vmvec3(&base_vec, len3, offset1 + offset2);
        assert(vmvec3.eq(&vmvec2));
    }

    void vec_utest()
    {
        std::cout << "vec_utest\n";

        vec_utest1();
        vec_utest2();
        vec_utest3();
    }
};

template class GF<uint32_t>;
template class GFP<uint32_t>;
template class GF2N<uint32_t>;
template class Vec<uint32_t>;

template class GF<uint64_t>;
template class GFP<uint64_t>;
template class GF2N<uint64_t>;
template class Vec<uint64_t>;

void vec_utest()
{
    VECUtest<uint32_t> vecutest_uint32;
    vecutest_uint32.vec_utest();
    VECUtest<uint64_t> vecutest_uint64;
    vecutest_uint64.vec_utest();
}
