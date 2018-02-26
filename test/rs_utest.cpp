#include "nttec.h"

/**
 * see http://web.eecs.utk.edu/~plank/plank/papers/CS-03-504.html
 */
void rs_utest1()
{
    std::cout << "rs_utest1\n";

    nttec::gf::BinExtension<uint32_t> gf16(4);
    nttec::vec::Matrix<uint32_t> mat(&gf16, 3, 3);

    mat.vandermonde_suitable_for_ec();
    // mat.dump();

    assert(mat.get(0, 0) == 1);
    assert(mat.get(0, 1) == 1);
    assert(mat.get(0, 2) == 1);
    assert(mat.get(1, 0) == 15);
    assert(mat.get(1, 1) == 8);
    assert(mat.get(1, 2) == 6);
    assert(mat.get(2, 0) == 14);
    assert(mat.get(2, 1) == 9);
    assert(mat.get(2, 2) == 6);
}

/**
 * Example taken from
 * An Introduction to Galois Fields and Reed-Solomon Coding
 * James Westall James Martin
 *
 * Note: the Vandermonde matrix is not generated the same way as we do
 */
void rs_utest2()
{
    std::cout << "rs_utest2\n";

    nttec::gf::BinExtension<uint32_t> gf8(3);
    nttec::vec::Matrix<uint32_t> mat(&gf8, 5, 3);
    nttec::vec::Vector<uint32_t> vec(&gf8, 3);
    nttec::vec::Vector<uint32_t> output(&gf8, 5);

    mat.set(0, 0, 1);
    mat.set(0, 1, 1);
    mat.set(0, 2, 6);
    mat.set(1, 0, 4);
    mat.set(1, 1, 3);
    mat.set(1, 2, 2);
    mat.set(2, 0, 5);
    mat.set(2, 1, 2);
    mat.set(2, 2, 2);
    mat.set(3, 0, 5);
    mat.set(3, 1, 3);
    mat.set(3, 2, 4);
    mat.set(4, 0, 4);
    mat.set(4, 1, 2);
    mat.set(4, 2, 4);
    // mat.dump();
    vec.set(0, 4);
    vec.set(1, 5);
    vec.set(2, 6);
    // vec.dump();
    mat.mul(&output, &vec);
    // output.dump();
    assert(output.get(0) == 3);
    assert(output.get(1) == 5);
    assert(output.get(2) == 4);
    assert(output.get(3) == 3);
    assert(output.get(4) == 2);
}

void rs_utest3()
{
    std::cout << "rs_utest3\n";

    nttec::gf::BinExtension<uint32_t> gf256(8);

    assert(gf256.mul(3, 7) == 9);
    assert(gf256.mul(13, 10) == 114);
    assert(gf256.div(13, 10) == 40);
    assert(gf256.div(3, 7) == 211);
}

void rs_utest()
{
    std::cout << "rs_utest\n";

    rs_utest1();
    rs_utest2();
    rs_utest3();
}
