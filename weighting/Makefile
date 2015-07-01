ROOTCFLAGS  := $(shell root-config --cflags --glibs)

NEUROBAYES_PATH=$(NEUROBAYES)
NEUROBAYES_LIBS=$(NEUROBAYES_PATH)/lib
NEUROBAYESFLAGS=-I$(NEUROBAYES_PATH)/include -L$(NEUROBAYES_LIBS) -lNeuroBayesTeacherCPP  -lNeuroBayes

TOOLSDIR  = $(HOME)/work/tools

CXX         = g++
CXXFLAGS    = -g -fPIC -Wall -Wno-write-strings -O2 $(ROOTCFLAGS) -lTMVA -lRooFit -lRooStats -I$(TOOLSDIR) -I$(TOOLSDIR)/analysis -L$(TOOLSDIR)/lib $(NEUROBAYESFLAGS)
#-fopenmp

SOURCES   = $(wildcard src/*.cc)
EXE       = lib/LbLmumuAng.so lib/LbJpsiLAng.so

PRINT=@echo "<*** $@ compiled ***>"



all: $(EXE)

#Compiling programs and linking libraries
lib/%.so: src/%.cc $(SOURCES)
	$(CXX) $(CXXFLAGS) -c $< -o $@ 
	$(PRINT)


clean:
	rm -f *.out

