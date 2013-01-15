C=g++
LD=g++

#O2 for optimization, g for debugging, pg for profiling
SPECIALFLAGS= -g #-g -Wall -pg#-O2
ROOTAUXCFLAGS=$(shell root-config --auxcflags)
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) #-lMath # -lMinuit
CFLAGS = $(SPECIALFLAGS) -Wall $(ROOTAUXCFLAGS)
LFLAGS = $(SPECIALFLAGS) #-lz 

RCXX=$(SPECIALFLAGS) -Wno-deprecated -Wall $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS) 

SRC=Event.cc Plot.cc NtupleSelector.cc JetResolution.cc

%.o: %.cc
		$(C) $(RCXX) -c $<

all: Plot

clean:
	@rm -f *.o 
	@rm -f Plot 


Event.o: Event.h Event.cc
	$(C) $(RCXX) -c Event.cc

NtupleSelector.o: NtupleSelector.h NtupleSelector.cc
	$(C) $(RCXX) -c NtupleSelector.cc

JetResolution.o: JetResolution.h JetResolution.cc
	$(C) $(RCXX) -c JetResolution.cc

Plot.o: Plot.cc
	$(C) $(RCXX) -c Plot.cc

Plot: $(SRC:.cc=.o) 
	$(LD) $(SRC:.cc=.o) $(RLXX) -lMathMore  -o Plot
	@echo '-> executable Plot created.'

