# $Id: Makefile,v 1.2 2009/05/04 14:35:04 mschrode Exp $

UTILDIR    = ./../util/


ROOTLIBS   = $(shell root-config --libs) -lMinuit
ROOTCFLAGS = $(shell root-config --cflags)

CFLAG      = -I$(UTILDIR) $(ROOTCFLAGS)
LFLAG      = $(ROOTLIBS)

OBJ        = lbfgs.o \
             Jet.o \
             Event.o \
             NJetEvent.o \
             PhotonJetEvent.o \
             EventGenerator.o \
             ControlPlots.o \
             Fitter.o \
             $(UTILDIR)/HistOps.o \
             run.o


run: $(OBJ)
	g++ $(OBJ) $(LFLAG) -lz -lg2c -o run
	@echo 'Created executable "run".'

lbfgs.o: lbfgs.F
	g77 -fno-automatic -fno-backslash -O -c lbfgs.F

Jet.o: Jet.h Jet.cc
	g++ $(CFLAG) -c Jet.cc

Event.o: Event.h Event.cc
	g++ $(CFLAG) -c Event.cc

NJetEvent.o: NJetEvent.h NJetEvent.cc
	g++ $(CFLAG) -c NJetEvent.cc

PhotonJetEvent.o: PhotonJetEvent.h PhotonJetEvent.cc
	g++ $(CFLAG) -c PhotonJetEvent.cc

EventGenerator.o: EventGenerator.h EventGenerator.cc
	g++ $(CFLAG) -c EventGenerator.cc

ControlPlots.o: ControlPlots.h ControlPlots.cc
	g++ $(CFLAG) -c ControlPlots.cc

Fitter.o: Fitter.h Fitter.cc
	g++ $(CFLAG) -c Fitter.cc

run.o: run.cc
	g++ $(CFLAG) -c run.cc



clean:
	@rm run
	@rm *~
	@rm *.o 
	@rm *#
	@rm *.ps
	@rm *.eps
	@rm -f .#*
