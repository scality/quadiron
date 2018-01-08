
#CXX = clang++

CXXFLAGS = -std=c++11 -I../src -g -O2 -Werror -Wall -Wl,--no-as-needed #-DNDEBUG
LDFLAGS = -L../src -lntl -lgmpxx -lgmp -lpthread -pthread

OBJS = benchmark.o

all: benchmark

benchmark: benchmark.o ../src/libntl.a
	$(CXX) -o benchmark benchmark.o $(LDFLAGS)

clean:
	$(RM) $(OBJS) benchmark
