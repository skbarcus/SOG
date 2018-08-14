ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -Wall  -fno-exceptions -fPIC  \
	-DLINUXVERS -I$(ROOTSYS)/include -O

LIBS = $(ROOTLIBS) $(ROOTGLIBS)

all: vdc_cross_talk_R vdc_cross_talk_L vdc_cross_talk_single_pair_R  vdc_cross_talk_single_pair_L vdc_cross_talk_single_pair_L_multicluster

vdc_cross_talk_R: vdc_cross_talk_R.C
	$(CXX) -g $(CXXFLAGS) -o $@ $< -I. $(LIBS)

vdc_cross_talk_L: vdc_cross_talk_L.C
	$(CXX) -g $(CXXFLAGS) -o $@ $< -I. $(LIBS)

vdc_cross_talk_single_pair_R: vdc_cross_talk_single_pair_R.C
	$(CXX) -g $(CXXFLAGS) -o $@ $< -I. $(LIBS)

vdc_cross_talk_single_pair_L: vdc_cross_talk_single_pair_L.C
	$(CXX) -g $(CXXFLAGS) -o $@ $< -I. $(LIBS)

vdc_cross_talk_single_pair_L_multicluster: vdc_cross_talk_single_pair_L_multicluster.C
	$(CXX) -g $(CXXFLAGS) -o $@ $< -I. $(LIBS)
