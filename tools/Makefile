
#CXX = clang++

CXXFLAGS = -std=c++11 -I../src -g -O2 #-DNDEBUG
LDFLAGS = -L../src -lntl -lgmpxx -lgmp

all: ec find_pminus1_compo #is_prime find_primes

ec: ec.o ../src/libntl.a
	$(CXX) -o ec ec.o $(LDFLAGS)

is_prime: is_prime.o ../src/libntl.a
	$(CXX) -o is_prime is_prime.o $(LDFLAGS)

find_primes: find_primes.o ../src/libntl.a
	$(CXX) -o find_primes find_primes.o $(LDFLAGS)

find_pminus1_compo: find_pminus1_compo.o ../src/libntl.a
	$(CXX) -o find_pminus1_compo find_pminus1_compo.o $(LDFLAGS)

clean:
	$(RM) *.o is_prime find_primes find_pminus1_compo ec
