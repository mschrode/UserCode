ROOTLIBS   = $(shell root-config --libs) -lMinuit
ROOTCFLAGS = $(shell root-config --cflags)

CFLAG      = -I $(ROOTCFLAGS)
LFLAG      = $(ROOTLIBS)

UTILDIR = ./../util
UTIL    = $(UTILDIR)/LabelFactory.h $(UTILDIR)/StyleSettings.h

OBJ     = CutVariation.o FittedResolution.o KalibriFileParser.o PtBin.o Uncertainty.o run.o

run: $(OBJ)	
	g++ $(OBJ) $(LFLAG) -o run
	@echo 'Created executable "run".'

CutVariation.o: CutVariation.h CutVariation.cc Uncertainty.h
	g++ $(CFLAG) -c CutVariation.cc

FittedResolution.o: FittedResolution.h FittedResolution.cc PtBin.h $(UTIL)
	g++ $(CFLAG) -c FittedResolution.cc

KalibriFileParser.o: KalibriFileParser.h KalibriFileParser.cc
	g++ $(CFLAG) -c KalibriFileParser.cc

PtBin.o: PtBin.h PtBin.cc Uncertainty.h CutVariation.h KalibriFileParser.h
	g++ $(CFLAG) -c PtBin.cc

Uncertainty.o: Uncertainty.h Uncertainty.cc Uncertainty.h
	g++ $(CFLAG) -c Uncertainty.cc

run.o: run.cc PtBin.h FittedResolution.h
	g++ $(CFLAG) -c run.cc



clean:
	@rm *.o 
	@rm run
	@rm *~
	@rm *#
	@rm *.ps
	@rm *.eps
	@rm -f .#*
