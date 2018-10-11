ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)	

CXX           = g++
CXXFLAGS      = -Wall  -fno-exceptions -fPIC  \
	-DLINUXVERS -I$(ROOTSYS)/include -O

LIBS = $(ROOTLIBS) $(ROOTGLIBS) $(ROOTFUNCS)

all: Global_Fit_3He_SOG Global_Fit_3H_SOG

Global_Fit_3He_SOG: Global_Fit_3He_SOG.C
	$(CXX) -g $(CXXFLAGS) -o $@ $< -I. $(LIBS) -lMathMore -lMinuit

Global_Fit_3H_SOG: Global_Fit_3H_SOG.C
	$(CXX) -g $(CXXFLAGS) -o $@ $< -I. $(LIBS) -lMathMore -lMinuit


