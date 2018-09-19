ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)	

CXX           = g++
CXXFLAGS      = -Wall  -fno-exceptions -fPIC  \
	-DLINUXVERS -I$(ROOTSYS)/include -O

LIBS = $(ROOTLIBS) $(ROOTGLIBS) $(ROOTFUNCS)

all: Test

Test: Test.C
	$(CXX) -g $(CXXFLAGS) -o $@ $< -I. $(LIBS) -lMathMore -lMinuit


