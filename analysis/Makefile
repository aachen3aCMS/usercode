 INCDIR = $(shell  root-config --incdir)
 
 CMSSW_RELEASE_BASE_NOPATCH:=$(shell echo $(CMSSW_RELEASE_BASE) | sed -e 's/-patch//' -e 's/_patch.//')
 CMSSW_INC:=-I$(CMSSW_RELEASE_BASE)/src -I$(CMSSW_BASE)/src
 CMSSW_LIBS:=-L$(CMSSW_BASE)/lib/$(SCRAM_ARCH) -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -L$(CMSSW_RELEASE_BASE)/external/$(SCRAM_ARCH)/lib -L$(CMSSW_RELEASE_BASE_NOPATCH)/lib/$(SCRAM_ARCH)
 ROOTGLIBS    := $(shell root-config --glibs)
 ROOTGLIBS   += -lX3d
 
 LIBS = $(ROOTGLIBS)
 LIBS += $(CMSSW_LIBS) -lCondFormatsJetMETObjects -lPhysicsToolsUtilities

CXXFLAGS = -g  --exceptions -I$(INCDIR) $(CMSSW_INC)

CPP    = g++

SRC= Run.cpp SUSYAna.cpp TreeContent.cpp

#=====================================================================

%.o: %.cpp
		$(CPP) -O3 $(CXXFLAGS) $(CMSSW_INC) -c -Wall -fPIC $<

all: test
test: Run
#==================
#$@ $(LDFLAGS) $^   $(LIBS)
Run.o: Run.cpp 
		$(CPP) $(CXXFLAGS) -c -O3 Run.cpp -Wall -fPIC
Run: $(SRC:.cpp=.o)
		$(CPP) $(CXXFLAGS) $(SRC:.cpp=.o) -o RunSUSY $(LIBS) -O3 -Wall -fPIC
		@echo
		@echo ' --> RunSUSY created'
		@echo

#====================================================================

clean:
		@echo "rm -f RunSUSY  *~  *.o"
		@rm -f *~ *.o
		@rm -f RunSUSY
		@echo ' --> clean'
