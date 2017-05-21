
CXXFLAGS = -std=c++11 -I../src -g
LDFLAGS = -L../src -lstl

OBJS = gf_utest.o rs_utest.o mat_utest.o main.o

all: tests

tests: $(OBJS) ../src/libntl.a
	$(CXX) -o tests $(OBJS)

clean:
	$(RM) $(OBJS) tests
