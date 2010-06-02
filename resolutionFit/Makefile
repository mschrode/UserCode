# $Id: Makefile,v 1.7 2010/05/28 18:52:29 mschrode Exp $

ROOTLIBS   = $(shell root-config --libs) -lMinuit
ROOTCFLAGS = $(shell root-config --cflags)

CFLAG      = -I $(ROOTCFLAGS)
LFLAG      = $(ROOTLIBS)

UTILDIR = ./../util
UTIL    = $(UTILDIR)/LabelFactory.h $(UTILDIR)/StyleSettings.h

OBJ     = CutVariation.o FittedResolution.o KalibriFileParser.o PtBin.o Parameters.o ResponseFunction.o Uncertainty.o run.o

run: $(OBJ)	
	g++ $(OBJ) $(LFLAG) -o run
	@echo 'Created executable "run".'

CutVariation.o: CutVariation.h CutVariation.cc Uncertainty.h KalibriFileParser.h Parameters.h
	g++ $(CFLAG) -c CutVariation.cc

FittedResolution.o: FittedResolution.h FittedResolution.cc PtBin.h Parameters.h ResponseFunction.h $(UTIL)
	g++ $(CFLAG) -c FittedResolution.cc

KalibriFileParser.o: KalibriFileParser.h KalibriFileParser.cc
	g++ $(CFLAG) -c KalibriFileParser.cc

PtBin.o: PtBin.h PtBin.cc Uncertainty.h CutVariation.h KalibriFileParser.h Parameters.h
	g++ $(CFLAG) -c PtBin.cc

Parameters.o: Parameters.h Parameters.cc ResponseFunction.h
	g++ $(CFLAG) -c Parameters.cc

ResponseFunction.o: ResponseFunction.h ResponseFunction.cc
	g++ $(CFLAG) -c ResponseFunction.cc

Uncertainty.o: Uncertainty.h Uncertainty.cc
	g++ $(CFLAG) -c Uncertainty.cc

run.o: run.cc PtBin.h FittedResolution.h Parameters.h ResponseFunction.h $(UTIL)
	g++ $(CFLAG) -c run.cc



clean:
	@rm -f *.o 
	@rm -f run
	@rm -f *~
	@rm -f *#
	@rm -f *.ps
	@rm -f *.eps
	@rm -f .#*
	@rm -f *.tex
