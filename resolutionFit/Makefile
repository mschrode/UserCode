# $Id: Makefile,v 1.12 2011/02/15 18:16:55 mschrode Exp $

ROOTLIBS   = $(shell root-config --libs) -lMinuit
ROOTCFLAGS = $(shell root-config --cflags)

CFLAG      = -I $(ROOTCFLAGS)
LFLAG      = $(ROOTLIBS)

UTILDIR = ./../util
UTIL    = $(UTILDIR)/ConfigParser.h $(UTILDIR)/LabelFactory.h $(UTILDIR)/HistOps.h $(UTILDIR)/StyleSettings.h $(UTILDIR)/utils.h
TOOLS =  Binning.h BinningAdmin.h

OBJ     = Binning.o BinningAdmin.o CommanderCool.o Extrapolation.o EtaBin.o FitResult.o Measurement.o OutputManager.o Parameters.o PlotMaker.o PtBin.o ResolutionFunction.o Sample.o run.o utils.o


run: $(OBJ)	
	g++ $(OBJ) $(LFLAG) -o run
	@echo 'Done. Type "./run" to plot everything about jet energy resolutions.'

utils.o: $(UTIL) $(UTILDIR)/utils.cc
	g++ $(CFLAG) -c  $(UTILDIR)/utils.cc

Binning.o: Binning.h Binning.cc
	g++ $(CFLAG) -c Binning.cc

BinningAdmin.o: BinningAdmin.h BinningAdmin.cc Binning.h $(UTIL)
	g++ $(CFLAG) -c BinningAdmin.cc

CommanderCool.o: CommanderCool.h CommanderCool.cc EtaBin.h FitResult.h ResolutionFunction.h Parameters.h PlotMaker.h PtBin.h
	g++ $(CFLAG) -c CommanderCool.cc

Extrapolation.o: Extrapolation.h Extrapolation.cc
	g++ $(CFLAG) -c Extrapolation.cc	

EtaBin.o: EtaBin.h EtaBin.cc FitResult.h Parameters.h PtBin.h ResolutionFunction.h Sample.h
	g++ $(CFLAG) -c EtaBin.cc	

FitResult.o: FitResult.h FitResult.cc Extrapolation.h Measurement.h
	g++ $(CFLAG) -c FitResult.cc	

ResolutionFunction.o: ResolutionFunction.h ResolutionFunction.cc
	g++ $(CFLAG) -c ResolutionFunction.cc

Measurement.o: Measurement.h Measurement.cc
	g++ $(CFLAG) -c Measurement.cc	

OutputManager.o: OutputManager.h OutputManager.cc $(UTIL)
	g++ $(CFLAG) -c OutputManager.cc	

Parameters.o: Parameters.h Parameters.cc OutputManager.h $(UTIL) $(TOOLS)
	g++ $(CFLAG) -c Parameters.cc	

PtBin.o: PtBin.h PtBin.cc FitResult.h Parameters.h Sample.h
	g++ $(CFLAG) -c PtBin.cc	

PlotMaker.o: PlotMaker.h PlotMaker.cc EtaBin.h FitResult.h OutputManager.h Parameters.h PtBin.h $(UTIL)
	g++ $(CFLAG) -c PlotMaker.cc	

Sample.o: Sample.h Sample.cc FitResult.h Measurement.h
	g++ $(CFLAG) -c Sample.cc	

run.o: run.cc CommanderCool.h FitResult.h Parameters.h ResolutionFunction.h $(UTIL)
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
