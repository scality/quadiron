#include "nttec.h"

extern void arith_utest();
extern void gf_utest();
extern void rs_utest();
extern void vector_utest();
extern void buffers_utest();
extern void mat_utest();
extern void fft_utest();
extern void poly_utest();
extern void fec_utest();

struct test {
    const char* name;
    void (*func)();
} tests[] = {{"arith", arith_utest},
             {"gf", gf_utest},
             {"vec", vector_utest},
             {"vecp", buffers_utest},
             {"mat", mat_utest},
             {"rs", rs_utest},
             {"poly", poly_utest},
             {"fft", fft_utest},
             {"fec", fec_utest},
             {nullptr, nullptr}};

int main(int argc, char** argv)
{
    struct test* p;

    srand(0);

    if (2 == argc) {
        std::stringstream ss(argv[1]);
        std::vector<int> output;

        std::string token;

        while (std::getline(ss, token, ',')) {
            for (p = tests; nullptr != p->name; p++) {
                std::string s(p->name);
                if (s == token) {
                    std::cout << p->name << "\n";
                    p->func();
                    break;
                }
            }
        }
    } else {
        for (p = tests; nullptr != p->name; p++) {
            std::cout << p->name << "\n";
            p->func();
        }
    }
}
