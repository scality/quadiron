
#include "ntl.h"

extern void gf_utest();
extern void rs_utest();

int main(int argc, char **argv)
{
  gf_utest();
  rs_utest();
}
