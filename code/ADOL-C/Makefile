CXX=/usr/local/bin/g++
CXXFLAGS = -O3

ADOLCPATH=$(HOME)/packages/adolc_edge
ADOLCLIBS=$(ADOLCPATH)/lib64
ADOLCINCL=$(ADOLCPATH)/include

all: adolc

adolc: adolc.cpp ./../obj_func.hpp
	$(CXX) $(CXXFLAGS) -I$(ADOLCINCL) adolc.cpp -o $@ -L$(ADOLCLIBS) -ladolc

clean:
	rm *.o *.tap adolc
