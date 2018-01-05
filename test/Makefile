
#CXX = clang++

CXXFLAGS = -std=c++11 -I../src -g -Werror -Wall
LDFLAGS = -L../src -lntl -lgmpxx -lgmp

OBJS = arith_utest.o gf_utest.o rs_utest.o vec_utest.o mat_utest.o fft_utest.o poly_utest.o fec_utest.o main.o

all: tests

tests: $(OBJS) ../src/libntl.a
	$(CXX) -o tests $(OBJS) $(LDFLAGS)

clean:
	$(RM) $(OBJS) tests
