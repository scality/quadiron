
#CXX = clang++

AR = ar cq

CXXFLAGS = -std=c++11 -g #-DNDEBUG

OBJS = gf.o gfp.o gf2n.o vec.o mat.o

LIBNTL = libntl.a

all: $(LIBNTL)

$(LIBNTL): $(OBJS)
	rm -f $(LIBNTL)
	$(AR) $(LIBNTL) $(OBJS)

clean:
	$(RM) $(OBJS) $(LIBNTL)
