CXX=/usr/local/bin/g++
CXXFLAGS = -O3

ADOLCPATH=$(HOME)/packages/adolc_edge
ADOLCLIBS=$(ADOLCPATH)/lib64
ADOLCINCL=$(ADOLCPATH)/include

all: direct indirect

direct: adolc.cpp ./../obj_func.hpp
	$(CXX) -D DIRECT $(CXXFLAGS) -I$(ADOLCINCL) adolc.cpp -o $@ -L$(ADOLCLIBS) -ladolc

indirect: adolc.cpp ./../obj_func.hpp
	$(CXX) -D INDIRECT $(CXXFLAGS) -I$(ADOLCINCL) adolc.cpp -o $@ -L$(ADOLCLIBS) -ladolc

clean:
	rm *.o *.tap direct indirect
