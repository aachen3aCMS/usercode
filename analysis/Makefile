#LIBS =  -L$(ROOTSYS)/lib -lNew  -lCint  \#
 LIBS = -L$(ROOTSYS)/lib -lNew  -lCint  \
         -lGraf -lGraf3d -lHist -lHtml -lMatrix  -lMinuit  \
         -lPostscript -lTree -lCore -lm -ldl  \
         -lGpad -lGui -lPhysics 
# -lDCache
#-lGX11 \
#         -L/usr/X11R6/lib -lm -ldl

CXXFLAGS = -g  --exceptions  -I$(ROOTSYS)/include 

CPP    = g++

SRC= Run.cpp SUSYAna.cpp TreeContent.cpp

#=====================================================================

%.o: %.cpp
		$(CPP) $(CXXFLAGS) -c -Wno-deprecated $<

all: test
test: Run
#==================

Run.o: Run.cpp 
		$(CPP) $(CXXFLAGS) -c Run.cpp -Wno-deprecated
Run: $(SRC:.cpp=.o)
		$(CPP) $(CXXFLAGS) $(SRC:.cpp=.o) $(LIBS) -o RunSUSY -Wno-deprecated
		@echo
		@echo ' --> RunSUSY created'
		@echo

#====================================================================

clean:
		@echo "rm -f RunSUSY  *~  *.o"
		@rm -f *~ *.o
		@rm -f RunSUSY
		@echo ' --> clean'
