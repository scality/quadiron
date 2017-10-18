#include "./ntl.h"

void _vec_hadamard_mul_257(int n,
                           uint32_t *x,
                           uint32_t *y)
{
  int i;

  for (i = 0; i < n/2; i++) {
    x[i] = (((uint64_t) x[i]) * y[i]) % 257;
  }

  for (int j = 0; j < n/2; j++) {
    x[i + j] = (((uint64_t) x[i + j]) * y[j]) % 257;
  }
}

void _vec_hadamard_mul_65537(int n,
                             uint32_t *x,
                             uint32_t *y)
{
  int i;

  for (i = 0; i < n/2; i++) {
    x[i] = (((uint64_t) x[i]) * y[i]) % 65537;
  }

  for (int j = 0; j < n/2; j++) {
    x[i + j] = (((uint64_t) x[i + j]) * y[j]) % 65537;
  }
}

void _vec_add_257(int n,
                  uint32_t *x,
                  uint32_t *y)
{
  int i;

  for (i = 0; i < n/2; i++) {
    x[i] = (x[i] + y[i]) % 257;
  }

  for (int j = 0; j < n/2; j++) {
    x[i + j] = (x[i + j] + y[j]) % 257;
  }
}

void _vec_add_65537(int n,
                    uint32_t *x,
                    uint32_t *y)
{
  int i;

  for (i = 0; i < n/2; i++) {
    x[i] = (x[i] + y[i]) % 65537;
  }

  for (int j = 0; j < n/2; j++) {
    x[i + j] = (x[i + j] + y[j]) % 65537;
  }
}

template <>
void Vec<uint32_t>::hadamard_mul(V2Vec<uint32_t> *v)
{
  assert(n == v->get_n());

  //typical butterfly operation
  if (rn->card() == 257) {
    uint32_t *a = get_mem();
    uint32_t *b = v->get_mem();
    _vec_hadamard_mul_257(n, a, b);
    return ;
  } else if (rn->card() == 65537) {
    uint32_t *a = get_mem();
    uint32_t *b = v->get_mem();
    _vec_hadamard_mul_65537(n, a, b);
    return ;
  }
  uint32_t *src = v->get_mem();
  int i;
  int j;
  for (i = 0; i < n / 2; i++)
    mem[i] = rn->mul(mem[i], src[i]);
  for (j = 0; i < n; i++, j++)
    mem[i] = rn->mul(mem[i], src[j]);
}

template <>
void Vec<uint32_t>::add(V2Vec<uint32_t> *v)
{
  assert(n == v->get_n());

  //typical butterfly operation
  if (rn->card() == 257) {
    uint32_t *a = get_mem();
    uint32_t *b = v->get_mem();
    _vec_add_257(n, a, b);
    return ;
  } else if (rn->card() == 65537) {
    uint32_t *a = get_mem();
    uint32_t *b = v->get_mem();
    _vec_add_65537(n, a, b);
    return ;
  }
  uint32_t *src = v->get_mem();
  int i;
  int j;
  for (i = 0; i < n / 2; i++)
    mem[i] = rn->add(mem[i], src[i]);
  for (j = 0; i < n; i++, j++)
    mem[i] = rn->add(mem[i], src[j]);
}

void _vec_hadamard_mul_257(int n,
                           uint64_t *x,
                           uint64_t *y)
{
  int i;

  for (i = 0; i < n/2; i++) {
    x[i] = (x[i] * y[i]) % 257;
  }

  for (int j = 0; j < n/2; j++) {
    x[i + j] = (x[i + j] * y[j]) % 257;
  }
}

void _vec_hadamard_mul_65537(int n,
                             uint64_t *x,
                             uint64_t *y)
{
  int i;

  for (i = 0; i < n/2; i++) {
    x[i] = (x[i] * y[i]) % 65537;
  }

  for (int j = 0; j < n/2; j++) {
    x[i + j] = (x[i + j] * y[j]) % 65537;
  }
}

void _vec_add_257(int n,
                  uint64_t *x,
                  uint64_t *y)
{
  int i;

  for (i = 0; i < n/2; i++) {
    x[i] = (x[i] + y[i]) % 257;
  }

  for (int j = 0; j < n/2; j++) {
    x[i + j] = (x[i + j] + y[j]) % 257;
  }
}

void _vec_add_65537(int n,
                    uint64_t *x,
                    uint64_t *y)
{
  int i;

  for (i = 0; i < n/2; i++) {
    x[i] = (x[i] + y[i]) % 65537;
  }

  for (int j = 0; j < n/2; j++) {
    x[i + j] = (x[i + j] + y[j]) % 65537;
  }
}

template <>
void Vec<uint64_t>::hadamard_mul(V2Vec<uint64_t> *v)
{
  assert(n == v->get_n());

  //typical butterfly operation
  if (rn->card() == 257) {
    uint64_t *a = get_mem();
    uint64_t *b = v->get_mem();
    _vec_hadamard_mul_257(n, a, b);
    return ;
  } else if (rn->card() == 65537) {
    uint64_t *a = get_mem();
    uint64_t *b = v->get_mem();
    _vec_hadamard_mul_65537(n, a, b);
    return ;
  }
  uint64_t *src = v->get_mem();
  int i;
  int j;
  for (i = 0; i < n / 2; i++)
    mem[i] = rn->mul(mem[i], src[i]);
  for (j = 0; i < n; i++, j++)
    mem[i] = rn->mul(mem[i], src[j]);
}

template <>
void Vec<uint64_t>::add(V2Vec<uint64_t> *v)
{
  assert(n == v->get_n());

  //typical butterfly operation
  if (rn->card() == 257) {
    uint64_t *a = get_mem();
    uint64_t *b = v->get_mem();
    _vec_add_257(n, a, b);
    return ;
  } else if (rn->card() == 65537) {
    uint64_t *a = get_mem();
    uint64_t *b = v->get_mem();
    _vec_add_65537(n, a, b);
    return ;
  }
  uint64_t *src = v->get_mem();
  int i;
  int j;
  for (i = 0; i < n / 2; i++)
    mem[i] = rn->add(mem[i], src[i]);
  for (j = 0; i < n; i++, j++)
    mem[i] = rn->add(mem[i], src[j]);
}
