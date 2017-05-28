
#CXX = clang++

CXXFLAGS = -std=c++11 -I../src -g
LDFLAGS = -L../src -lntl -lgmpxx -lgmp

all: is_prime find_primes

is_prime: is_prime.o ../src/libntl.a
	$(CXX) -o is_prime is_prime.o $(LDFLAGS)

find_primes: find_primes.o ../src/libntl.a
	$(CXX) -o find_primes find_primes.o $(LDFLAGS)

clean:
	$(RM) *.o is_prime find_primes
