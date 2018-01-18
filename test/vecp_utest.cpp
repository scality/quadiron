
#include "ntl.h"

template<typename T>
class VECPUtest
{
 private:
   GFP<T> *gfp;

 public:
  VECPUtest() {
    this->gfp = new GFP<T>(65537);
  }
  ~VECPUtest() {
    delete this->gfp;
  }

  Vecp<T>* gen_vecp_rand_data(int n, int size) {
    Vecp<T> *vec = new Vecp<T>(n, size);
    for (int i = 0; i < n; i++) {
      T *buf = new T[size];
      for (int j = 0; j < size; j++) {
        buf[j] = this->gfp->weak_rand();
      }
      vec->set(i, buf);
    }
    // vec.dump();

    return vec;
  }

  void vecp_utest1()
  {
    std::cout << "vecp_utest1\n";

    int n = 16;
    int begin = 5;
    int end = 12;
    int size = 32;
    int i, j;

    Vecp<T> *vec1 = gen_vecp_rand_data(n, size);
    std::vector<T*> *mem1 = vec1->get_mem();
    // vec1->dump();

    Vecp<T> *vec2 = vec1->slice(begin, end);
    std::vector<T*> *mem2 = vec2->get_mem();
    assert(vec2->get_n() == end - begin);
    assert(vec2->get_size() == vec1->get_size());
    for (i = 0; i < end - begin; i++) {
      for (j = 0; j < size; j++) {
        mem2->at(i)[j] = mem1->at(i + begin)[j];
      }
    }
    // vec2->dump();

    std::vector<T*> mem3(end - begin, nullptr);
    for (int i = 0; i < end - begin; i++) {
      mem3[i] = mem1->at(i + begin);
    }
    Vecp<T> vec3(end - begin, size, &mem3);
    // vec3.dump();

    assert(vec2->eq(&vec3));

    delete vec1;
    delete vec2;
  }

  void vecp_utest2()
  {
    std::cout << "vecp_utest2\n";

    int n = 8;
    int size = 32;
    int i, j;
    int half = n / 2;

    Vecp<T> *vec1 = gen_vecp_rand_data(n, size);
    Vecp<T> vec2(n, size);
    vec2.copy(vec1);

    std::vector<T*> *even_mem = new std::vector<T*>(half, nullptr);
    std::vector<T*> *odd_mem = new std::vector<T*>(half, nullptr);
    Vecp<T> *i_even = new Vecp<T>(half, size, even_mem);
    Vecp<T> *i_odd = new Vecp<T>(half, size, odd_mem);
    vec1->separate_even_odd(i_even, i_odd);

    // vec1->dump();
    vec1->separate_even_odd();
    // vec1->dump();

    assert(i_even->eq(vec1->slice(0, half)));
    assert(i_odd->eq(vec1->slice(half, n)));

    // vec2.dump();

    std::vector<T*> *mem1 = vec1->get_mem();
    std::vector<T*> *mem2 = vec2.get_mem();

    bool ok = true;
    for (i = 0; i < n/2; i += 2) {
      T *even1 = mem1->at(i);
      T *even2 = mem2->at(i*2);
      T *odd1 = mem1->at(i+n/2);
      T *odd2 = mem2->at(i*2+1);
      for (j = 0; j < size; j++) {
        if(even1[j] != even2[j] ||
           odd1[j] != odd2[j]) {
          ok = false;
          i = n;
          break;
        }
      }
    }
    assert(ok);


    delete vec1;
  }

  void vecp_utest()
  {
    std::cout << "vecp_utest with sizeof(T)=" << sizeof(T) << "\n";

    vecp_utest1();
    vecp_utest2();
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

void vecp_utest()
{
  VECPUtest<uint32_t> vecputest_uint32;
  vecputest_uint32.vecp_utest();
  VECPUtest<uint64_t> vecputest_uint64;
  vecputest_uint64.vecp_utest();
}
