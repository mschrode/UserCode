# Make object files of util classes

ROOTLIBS   = $(shell root-config --libs)
ROOTCFLAGS = $(shell root-config --cflags)

CFLAG      = -I$(DIR) $(ROOTCFLAGS) -O2
LFLAG      = $(ROOTLIBS)

OBJ        = HistOps.o


HistOps.o: HistOps.h HistOps.cc
	g++ $(CFLAG) -c HistOps.cc
