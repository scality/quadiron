
#CXX = clang++

AR = ar cq

CXXFLAGS = -std=c++11 -g -O2 #-DNDEBUG

OBJS = misc.o config.o vec.o

LIBNTL = libntl.a

all: $(LIBNTL)

$(LIBNTL): $(OBJS)
	rm -f $(LIBNTL)
	$(AR) $(LIBNTL) $(OBJS)

clean:
	$(RM) $(OBJS) $(LIBNTL)
