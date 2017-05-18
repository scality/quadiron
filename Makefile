
CXXFLAGS = -std=c++11

OBJS = main.o gf.o gfp.o gf2n.o

all: main

main: $(OBJS)
	$(CXX) -o main $(OBJS)

clean:
	$(RM) $(OBJS) main
