
#CXX = clang++

CXXFLAGS = -std=c++11 -I../src -g
LDFLAGS = -L../src -lntl -lgmpxx -lgmp

OBJS = is_prime.o

all: is_prime

is_prime: $(OBJS) ../src/libntl.a
	$(CXX) -o is_prime is_prime.o $(LDFLAGS)

clean:
	$(RM) $(OBJS) is_prime
