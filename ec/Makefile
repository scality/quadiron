
#CXX = clang++

CXXFLAGS = -std=c++11 -I../src -g
LDFLAGS = -L../src -lstl

OBJS = main.o

all: ec

ec: $(OBJS) ../src/libntl.a
	$(CXX) -o ec $(OBJS)

clean:
	$(RM) $(OBJS) ec
