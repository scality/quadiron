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

template <typename T>
class BuffersUtest {
  private:
    nttec::gf::Prime<T>* gfp;
    T max_val;

  public:
    BuffersUtest()
    {
        this->gfp = new nttec::gf::Prime<T>(65537);
        this->max_val = 65537;
    }
    ~BuffersUtest()
    {
        delete this->gfp;
    }

    nttec::vec::Buffers<T>* gen_buffers_rand_data(int n, int size, int _max = 0)
    {
        int max = (_max == 0) ? max_val : _max;
        nttec::vec::Buffers<T>* vec = new nttec::vec::Buffers<T>(n, size);
        for (int i = 0; i < n; i++) {
            T* buf = new T[size];
            for (int j = 0; j < size; j++) {
                buf[j] = rand() % max;
            }
            vec->set(i, buf);
        }
        // vec.dump();

        return vec;
    }

    void buffers_utest1()
    {
        std::cout << "buffers_utest1\n";

        int n = 16;
        int begin = 5;
        int end = 12;
        int size = 32;
        int i, j;

        nttec::vec::Buffers<T>* vec1 = gen_buffers_rand_data(n, size);
        std::vector<T*>* mem1 = vec1->get_mem();
        // vec1->dump();

        nttec::vec::Buffers<T> vec2(vec1, begin, end);
        std::vector<T*>* mem2 = vec2.get_mem();
        assert(vec2.get_n() == end - begin);
        assert(vec2.get_size() == vec1->get_size());
        for (i = 0; i < end - begin; i++) {
            for (j = 0; j < size; j++) {
                mem2->at(i)[j] = mem1->at(i + begin)[j];
            }
        }
        // vec2.dump();

        std::vector<T*> mem3(end - begin, nullptr);
        for (int i = 0; i < end - begin; i++) {
            mem3[i] = mem1->at(i + begin);
        }
        nttec::vec::Buffers<T> vec3(end - begin, size, &mem3);
        // vec3.dump();

        assert(vec2.eq(&vec3));

        delete vec1;
    }

    void buffers_utest2()
    {
        std::cout << "buffers_utest2\n";

        int n = 8;
        int size = 32;
        int i, j;
        int half = n / 2;

        nttec::vec::Buffers<T>* vec1 = gen_buffers_rand_data(n, size);
        nttec::vec::Buffers<T> vec2(n, size);
        vec2.copy(vec1);

        std::vector<T*>* even_mem = new std::vector<T*>(half, nullptr);
        std::vector<T*>* odd_mem = new std::vector<T*>(half, nullptr);
        nttec::vec::Buffers<T>* i_even =
            new nttec::vec::Buffers<T>(half, size, even_mem);
        nttec::vec::Buffers<T>* i_odd =
            new nttec::vec::Buffers<T>(half, size, odd_mem);
        vec1->separate_even_odd(i_even, i_odd);

        // vec1->dump();
        vec1->separate_even_odd();
        // vec1->dump();

        nttec::vec::Buffers<T> _i_even(vec1, 0, half);
        nttec::vec::Buffers<T> _i_odd(vec1, half, n);
        assert(i_even->eq(&_i_even));
        assert(i_odd->eq(&_i_odd));

        // vec2.dump();

        std::vector<T*>* mem1 = vec1->get_mem();
        std::vector<T*>* mem2 = vec2.get_mem();

        bool ok = true;
        for (i = 0; i < n / 2; i += 2) {
            T* even1 = mem1->at(i);
            T* even2 = mem2->at(i * 2);
            T* odd1 = mem1->at(i + n / 2);
            T* odd2 = mem2->at(i * 2 + 1);
            for (j = 0; j < size; j++) {
                if (even1[j] != even2[j] || odd1[j] != odd2[j]) {
                    ok = false;
                    i = n;
                    break;
                }
            }
        }
        assert(ok);

        delete vec1;
    }

    void buffers_utest3()
    {
        std::cout << "buffers_utest3\n";

        int n = 8;
        int size = 32;
        int n1 = 4;
        int n2 = 10;

        nttec::vec::Buffers<T>* vec = gen_buffers_rand_data(n, size);
        nttec::vec::Buffers<T> vec1(vec, n1);
        nttec::vec::Buffers<T> vec2(vec, n2);

        nttec::vec::Buffers<T> _vec1(vec, 0, n1);
        nttec::vec::Buffers<T>* _vec2 =
            new nttec::vec::BuffersZeroExtended<T>(vec, n2);

        assert(vec1.eq(&_vec1));
        assert(vec2.eq(_vec2));

        delete vec;

        nttec::vec::Buffers<T> vec3(&vec2, n1);
        assert(vec3.eq(&vec1));

        delete _vec2;
    }

    void buffers_utest4()
    {
        std::cout << "buffers_utest4\n";
        int i;
        int word_size;

        for (i = 0; i <= nttec::arith::log2<T>(sizeof(T)); i++) {
            word_size = nttec::arith::exp2<T>(i);
            pack_unpack(word_size);
        }
    }

    void pack_unpack(int word_size)
    {
        std::cout << "pack_unpack with word_size=" << word_size << "\n";
        int n = 8;
        int size = 32;
        int bytes_size = size * word_size;
        int i, j, t, u;

        T symb;
        T max = ((T)1 << word_size) + 1;

        nttec::vec::Buffers<T>* words = gen_buffers_rand_data(n, size, max);
        std::vector<T*>* mem_T = words->get_mem();
        // std::cout << "words:"; words->dump();

        // pack manually from T to uint8_t
        nttec::vec::Buffers<uint8_t> vec_char(n, bytes_size);
        std::vector<uint8_t*>* mem_char = vec_char.get_mem();
        for (i = 0; i < n; i++) {
            t = 0;
            T* buf_T = mem_T->at(i);
            uint8_t* buf_char = mem_char->at(i);
            for (j = 0; j < size; j++) {
                symb = buf_T[j];
                buf_char[t] = (uint8_t)(symb & 0xFF);
                t++;
                for (u = 1; u < word_size; u++) {
                    symb >>= 8;
                    buf_char[t] = (uint8_t)(symb & 0xFF);
                    t++;
                }
            }
        }

        /*
         * pack bufs of type uint8_t to bufs of type T
         */
        // tmp vectors to store results
        nttec::vec::Buffers<T> vec_T_tmp(n, size);
        std::vector<T*>* mem_T_tmp = vec_T_tmp.get_mem();
        nttec::vec::pack<uint8_t, T>(mem_char, mem_T_tmp, n, size, word_size);
        // std::cout << "vec_char:"; vec_char.dump();
        // std::cout << "vec_T_tmp:"; vec_T_tmp.dump();
        // check
        assert(vec_T_tmp.eq(words));

        /*
         * unpack bufs of type T to bufs of type uint8_t
         */
        // tmp vectors to store results
        nttec::vec::Buffers<uint8_t> vec_char_tmp(n, bytes_size);
        std::vector<uint8_t*>* mem_char_tmp = vec_char_tmp.get_mem();
        nttec::vec::unpack<T, uint8_t>(
            mem_T_tmp, mem_char_tmp, n, size, word_size);
        // std::cout << "vec_T_tmp:"; vec_T_tmp.dump();
        // std::cout << "vec_char_tmp:"; vec_char_tmp.dump();
        // check
        assert(vec_char_tmp.eq(&vec_char));

        delete words;
    }

    void buffers_utest()
    {
        std::cout << "buffers_utest with sizeof(T)=" << sizeof(T) << "\n";

        buffers_utest1();
        buffers_utest2();
        buffers_utest3();
        buffers_utest4();
    }
};

void buffers_utest()
{
    BuffersUtest<uint32_t> buffersutest_uint32;
    buffersutest_uint32.buffers_utest();
    BuffersUtest<uint64_t> buffersutest_uint64;
    buffersutest_uint64.buffers_utest();
}
