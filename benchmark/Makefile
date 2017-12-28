
#CXX = clang++

CXXFLAGS = -std=c++11  -I../benchmark -I../src -g -O2 -DNDEBUG
LDFLAGS = -L../src -lntl -lgmpxx -lgmp

OBJS = benchmark.o

all: benchmark

benchmark: benchmark.o ../src/libntl.a
	$(CXX) -o benchmark benchmark.o $(LDFLAGS)

clean:
	$(RM) $(OBJS) benchmark
