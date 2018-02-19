#include "nttec.h"

void mat_utest1()
{
    std::cout << "mat_utest1\n";

    nttec::gf::Prime<uint32_t> gf11(11);
    nttec::Matrix<uint32_t> mat(&gf11, 3, 3);

    mat.set(0, 0, 2);
    mat.set(0, 1, 1);
    mat.set(0, 2, 2);
    mat.set(1, 0, 1);
    mat.set(1, 1, 2);
    mat.set(1, 2, 9);
    mat.set(2, 0, 1);
    mat.set(2, 1, 2);
    mat.set(2, 2, 7);

    mat.inv();

    assert(mat.get(0, 0) == 8);
    assert(mat.get(0, 1) == 6);
    assert(mat.get(0, 2) == 1);
    assert(mat.get(1, 0) == 7);
    assert(mat.get(1, 1) == 9);
    assert(mat.get(1, 2) == 10);
    assert(mat.get(2, 0) == 0);
    assert(mat.get(2, 1) == 6);
    assert(mat.get(2, 2) == 5);
}

void mat_utest2()
{
    std::cout << "mat_utest2\n";

    nttec::gf::Prime<uint32_t> gf29(29);
    nttec::Matrix<uint32_t> mat(&gf29, 3, 3);

    mat.set(0, 0, 22);
    mat.set(0, 1, 27);
    mat.set(0, 2, 18);
    mat.set(1, 0, 18);
    mat.set(1, 1, 28);
    mat.set(1, 2, 5);
    mat.set(2, 0, 4);
    mat.set(2, 1, 17);
    mat.set(2, 2, 1);

    mat.inv();

    assert(mat.get(0, 0) == 1);
    assert(mat.get(0, 1) == 18);
    assert(mat.get(0, 2) == 8);
    assert(mat.get(1, 0) == 2);
    assert(mat.get(1, 1) == 8);
    assert(mat.get(1, 2) == 11);
    assert(mat.get(2, 0) == 20);
    assert(mat.get(2, 1) == 24);
    assert(mat.get(2, 2) == 14);
}

void mat_utest()
{
    std::cout << "mat_utest\n";

    mat_utest1();
    mat_utest2();
}
