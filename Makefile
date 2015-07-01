ROOTCFLAGS  := $(shell root-config --cflags --glibs)

NEUROBAYES_PATH=$(NEUROBAYES)
NEUROBAYES_LIBS=$(NEUROBAYES_PATH)/lib
NEUROBAYESFLAGS=-I$(NEUROBAYES_PATH)/include -L$(NEUROBAYES_LIBS) -lNeuroBayesTeacherCPP  -lNeuroBayes

TOOLSDIR  = $(LHCBBHAM)/Tools
LBDIR     = $(HOME)/work/Lb/Lmumu

CXX         = g++
CXXFLAGS    = -g -fPIC -Wall -Wno-write-strings -O2 $(ROOTCFLAGS) -lTMVA -lRooFit -lRooStats -I$(TOOLSDIR) -I$(TOOLSDIR)/analysis -L$(TOOLSDIR)/lib $(NEUROBAYESFLAGS) -I$(LBDIR)/weighting/include -L$(LBDIR)/weighting/lib
CXXFLAGS_NONB    = -g -fPIC -Wall -Wno-write-strings -O2 $(ROOTCFLAGS) -lTMVA -lRooFit -lRooStats -I$(TOOLSDIR) -I$(TOOLSDIR)/analysis -L$(TOOLSDIR)/lib -I$(LBDIR)/weighting/include -L$(LBDIR)/weighting/lib
#-fopenmp

LIBS      = $(patsubst $(TOOLSDIR)/%.hpp, $(TOOLSDIR)/lib/%.a, $(wildcard $(TOOLSDIR)/*.hpp)) #$(LBDIR)/weighting/lib/*.so 
MODEL     = $(wildcard $(LBDIR)/weighting/src/*.cc)
SOURCES   = $(wildcard src/*.cpp) $(wildcard $(LBDIR)/weighting/src/*.cc)
EXE       = $(patsubst src/%.cpp, %.out, $(SOURCES))

PRINT=@echo "<*** $@ compiled ***>"



all: libraries $(EXE)

exe: $(EXE)

#Compiling libraries
libraries:
	$(MAKE) -C $(TOOLSDIR)


#Compiling programs and linking libraries
%.out: src/%.cpp $(LIBS) $(MODEL)
	$(CXX) $(CXXFLAGS_NONB) $^ -o $@ 
	$(PRINT)


clean:
	rm -f *.out

veryclean:
	rm -f *.out
	$(MAKE) -C $(TOOLSDIR) veryclean
