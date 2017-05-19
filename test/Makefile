
CXXFLAGS = -std=c++11 -I../src
LDFLAGS = -L../src -lstl

OBJS = gf_utest.o main.o

all: tests

tests: $(OBJS)
	$(CXX) -o tests $(OBJS)

clean:
	$(RM) $(OBJS) tests
