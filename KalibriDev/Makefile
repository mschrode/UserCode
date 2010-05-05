C=g++
LD=g++
GCC_VERSION=$(shell g++ --version | head -1 | cut -f 3 -d' ' | cut -f 1 -d'.')
ifeq ($(GCC_VERSION), 3)
 F77=g77
 F77LDFLAGS=-lg2c
else
 F77=gfortran
 F77LDFLAGS=-lgfortran
endif


#O2 for optimization, g for debugging, pg for profiling
SPECIALFLAGS= -fpic -g -O2 -pg# -O2
ROOTAUXCFLAGS=$(shell root-config --auxcflags)
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) -lMinuit
#-I. -I./include -I$(SRT_PUBLIC_CONTEXT)/include 
CFLAGS= $(SPECIALFLAGS) -Wall $(ROOTAUXCFLAGS)
#-L../../lib/$(SRT_SUBDIR)/
LFLAGS= $(SPECIALFLAGS) -lz $(F77LDFLAGS) -lgsl -lgslcblas -lm
BOOSTFLAGS=-I/usr/include/boost
BOOSTLINKFLAGS=-lboost_thread -lpthread
# change path for MacPort or fink@MacOS
 ifeq (exists, $(shell [ -d /opt/local/include/boost ] && echo exists)) 
 BOOSTFLAGS=-I/opt/local/include/boost
 SPECIALFLAGS += -I/opt/local/include
 LFLAGS += -L/opt/local/lib
 BOOSTLINKFLAGS=-lboost_thread-mt -lpthread
else ifeq (exists, $(shell [ -d /sw/include/boost ] && echo exists)) 
 BOOSTFLAGS=-I/sw/include/boost
 SPECIALFLAGS += -arch i386 -I/sw/include
 LFLAGS += -L/sw/lib
endif

RCXX=$(SPECIALFLAGS) -Wall $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS) $(BOOSTLINKFLAGS)  #-lrt -lpthread # -lposix4
ROOTSYS=$(shell root-config --prefix)

SRC=Kalibri.cc GammaJetSel.cc ZJetSel.cc NJetSel.cc TopSel.cc ConfigFile.cc CalibData.cc Parametrization.cc Parameters.cc ControlPlots.cc ControlPlotsProfile.cc ControlPlotsFunction.cc ControlPlotsConfig.cc ControlPlotsJetSmearing.cc ToyMC.cc EventReader.cc PhotonJetReader.cc DiJetReader.cc TriJetReader.cc ZJetReader.cc TopReader.cc ParameterLimitsReader.cc EventProcessor.cc EventWeightProcessor.cc Jet.cc JetTruthEvent.cc JetWithTowers.cc TwoJetsInvMassEvent.cc TwoJetsPtBalanceEvent.cc JetWithTracks.cc SmearData.cc SmearDiJet.cc SmearPhotonJet.cc JetConstraintEvent.cc CorFactorsFactory.cc

%.o: %.cc
		$(C) $(RCXX) -c $<

all: dirs lib bin

dirs:
	@mkdir -p bin
	@mkdir -p lib
	@mkdir -p tmp

lbfgs.o: lbfgs.F
	$(F77) $(RCXX) -fno-automatic -fno-backslash -O -c lbfgs.F

ConfigFile.o: ConfigFile.cc ConfigFile.h
	$(C) $(CFLAGS) -c ConfigFile.cc

GammaJetSel.o: GammaJetSel.cc GammaJetSel.h
	$(C) $(RCXX) -c GammaJetSel.cc

ZJetSel.o: ZJetSel.cc ZJetSel.h
	$(C) $(RCXX) -c ZJetSel.cc

TopSel.o: TopSel.cc TopSel.h
	$(C) $(RCXX) -c TopSel.cc

NJetSel.o: NJetSel.cc NJetSel.h
	$(C) $(RCXX) -c NJetSel.cc

CalibData.o: CalibData.cc CalibData.h Parametrization.h Parameters.h
	$(C) $(CFLAGS) -c CalibData.cc

SmearData.o: SmearData.cc SmearData.h CalibData.h SmearFunction.h
	$(C) $(CFLAGS) -c SmearData.cc

SmearDiJet.o: SmearDiJet.cc SmearDiJet.h SmearFunction.h Jet.h SmearData.h
	$(C) $(CFLAGS) -c SmearDiJet.cc

SmearPhotonJet.o: SmearPhotonJet.cc SmearPhotonJet.h Jet.h SmearData.h SmearFunction.h
	$(C) $(CFLAGS) -c SmearPhotonJet.cc

Parametrization.o: Parametrization.h Parametrization.cc
	$(C) $(RCXX) -c Parametrization.cc

Parameters.o: Parameters.cc Parameters.h Parametrization.h Function.h SmearFunction.h ConfigFile.h
	$(C) $(RCXX) -c Parameters.cc

ControlPlots.o: ControlPlots.cc ControlPlots.h ControlPlotsConfig.h ControlPlotsProfile.h ControlPlotsFunction.h CalibData.h Function.h Jet.h JetTruthEvent.h
	$(C) $(RCXX) -c ControlPlots.cc

ControlPlotsProfile.o: ControlPlotsProfile.cc ControlPlotsProfile.h ControlPlotsConfig.h ControlPlotsFunction.h ConfigFile.h CalibData.h
	$(C) $(RCXX) -c ControlPlotsProfile.cc

ControlPlotsFunction.o: ControlPlotsFunction.cc ControlPlotsFunction.h ControlPlotsConfig.h Jet.h JetTruthEvent.h CorFactors.h
	$(C) $(RCXX) -c ControlPlotsFunction.cc

ControlPlotsConfig.o: ControlPlotsConfig.cc ControlPlotsConfig.h ConfigFile.h
	$(C) $(RCXX) -c ControlPlotsConfig.cc

ControlPlotsJetSmearing.o: ControlPlotsJetSmearing.cc ControlPlotsJetSmearing.h CalibData.h ConfigFile.h Parameters.h SmearData.h SmearDiJet.h SmearPhotonJet.h Jet.h 	
	$(C) $(RCXX) -c ControlPlotsJetSmearing.cc

EventReader.o: EventReader.h EventReader.cc Parameters.h Parametrization.h ConfigFile.h CorFactorsFactory.h CorFactors.h ToyMC.h  JetConstraintEvent.h
	$(C) $(RCXX) -c EventReader.cc

PhotonJetReader.o: EventReader.h PhotonJetReader.h PhotonJetReader.cc  GammaJetSel.h ToyMC.h Parameters.h ConfigFile.h Jet.h JetTruthEvent.h JetWithTowers.h Function.h CorFactors.h CorFactorsFactory.h
	$(C) $(RCXX) -c PhotonJetReader.cc

DiJetReader.o: EventReader.h DiJetReader.h DiJetReader.cc NJetSel.h Parameters.h ConfigFile.h Jet.h JetTruthEvent.h TwoJetsPtBalanceEvent.h JetWithTowers.h Function.h SmearFunction.h CorFactors.h CorFactorsFactory.h JetConstraintEvent.h
	$(C) $(RCXX) -c DiJetReader.cc

TriJetReader.o: EventReader.h TriJetReader.h TriJetReader.cc NJetSel.h Parameters.h ConfigFile.h CorFactors.h CorFactorsFactory.h
	$(C) $(RCXX) -c TriJetReader.cc

ZJetReader.o: EventReader.h ZJetReader.h ZJetReader.cc ZJetSel.h Parameters.h ConfigFile.h Jet.h JetTruthEvent.h JetWithTowers.h JetWithTracks.h Function.h CorFactors.h CorFactorsFactory.h
	$(C) $(RCXX) -c ZJetReader.cc

TopReader.o: EventReader.h TopReader.h TopReader.cc TopSel.h Parameters.h ConfigFile.h CorFactors.h CorFactorsFactory.h Jet.h
	$(C) $(RCXX) -c TopReader.cc

ParameterLimitsReader.o: EventReader.h ParameterLimitsReader.h ParameterLimitsReader.cc Parameters.h ConfigFile.h CorFactorsFactory.h
	$(C) $(CFLAGS) -c ParameterLimitsReader.cc

EventProcessor.o: CalibData.h ConfigFile.h Parameters.h EventProcessor.h EventProcessor.cc
	$(C) $(RCXX) -c EventProcessor.cc

EventWeightProcessor.o: CalibData.h ConfigFile.h EventProcessor.h Parameters.h EventWeightProcessor.cc
	$(C) $(RCXX) -c EventWeightProcessor.cc

Jet.o: CalibData.h Jet.h Jet.cc Parametrization.h Function.h CorFactors.h
	$(C) $(CFLAGS) -c Jet.cc	

JetTruthEvent.o: CalibData.h Jet.h JetTruthEvent.h JetTruthEvent.cc
	$(C) $(CFLAGS) -c JetTruthEvent.cc

TwoJetsInvMassEvent.o: CalibData.h Jet.h TwoJetsInvMassEvent.h TwoJetsInvMassEvent.cc
	$(C) $(RCXX) -c TwoJetsInvMassEvent.cc

TwoJetsPtBalanceEvent.o: CalibData.h Jet.h TwoJetsPtBalanceEvent.h TwoJetsPtBalanceEvent.cc CorFactors.h
	$(C) $(RCXX) -c TwoJetsPtBalanceEvent.cc

JetConstraintEvent.o: JetConstraintEvent.h CalibData.h Jet.h JetConstraintEvent.cc Function.h
	$(C) $(CFLAGS) -c JetConstraintEvent.cc

JetWithTowers.o: CalibData.h Jet.h JetWithTowers.h Function.h JetWithTowers.cc Parametrization.h
	$(C) $(RCXX) -c JetWithTowers.cc

JetWithTracks.o: CalibData.h Jet.h JetWithTracks.h Function.h JetWithTracks.cc Parametrization.h
	$(C) $(RCXX) -c JetWithTracks.cc

CorFactorsFactory.o: CorFactors.h CorFactorsFactory.h CorFactorsFactory.cc
	$(C) $(CFLAGS) -c CorFactorsFactory.cc



lib/libKalibri.so: $(SRC:.cc=.o) lbfgs.o
	$(LD) $(RCXX) -shared $^ $(RLXX) -o lib/libKalibri.so
	@echo '-> Kalibri library created.'

Kalibri.o: Kalibri.cc Kalibri.h CalibMath.h external.h ConfigFile.h CalibData.h Parameters.h ControlPlots.h EventReader.h DiJetReader.h TriJetReader.h ZJetReader.h TopReader.h ParameterLimitsReader.h  EventProcessor.h  EventWeightProcessor.h Jet.h TwoJetsInvMassEvent.h TwoJetsPtBalanceEvent.h 
	$(C) $(CFLAGS) $(RCXX) $(BOOSTFLAGS) -c Kalibri.cc 

caliber.o: caliber.cc Kalibri.h JetTruthEvent.h Jet.h
	$(C) $(RCXX) -c caliber.cc


lib: dirs lib/libKalibri.so

bin: dirs junk caliber

junk: $(SRC:.cc=.o) lbfgs.o caliber.o
	$(LD) $^ $(RLXX) -o bin/junk
	@ln -s -f bin/junk
	@echo '-> static Kalibri executable created.'

caliber: caliber.o lib/libKalibri.so
	$(LD) caliber.o $(RLXX) -Llib -lKalibri -o bin/caliber
	@ln -f -s bin/caliber
	@echo '-> shared Kalibri executable created.'

clean:
	@rm -rf ti_files
	@rm -f *~
	@rm -f *.o 
	@rm -f *#
	@rm -f .nfs*
	@rm -f *.bkp
	@rm -f junk
	@rm -f caliber
	@rm -f libKalibri.so
	@rm -rf tmp
	@rm -rf lib
	@rm -rf bin
	@rm -f *.cfi
	@rm -f fort.*
	@rm -f .#*



ToyMC.o: ToyMC.h ToyMC.cc ConfigFile.h
	$(C) $(RCXX) -c ToyMC.cc

toy:	ToyMC.o toy.o ConfigFile.o
	$(LD) ToyMC.o toy.o ConfigFile.o $(RLXX) -o toy
	@echo '-> toy MC executable created.'

ControlPlotsComparison.o: ControlPlotsComparison.cc ControlPlotsComparison.h
	$(C) $(RCXX) -c ControlPlotsComparison.cc

comp: 	ControlPlotsComparison.o compareControlPlots.o
	$(LD) ControlPlotsComparison.o compareControlPlots.o $(RLXX) -o compControlPlots
	@echo '-> Comparison executable created. Type "compControlPlots" to compare control plots.'

lib/libJetMETObjects.so: dirs JetMETObjects
	@env STANDALONE_DIR=${PWD} ROOTSYS=${ROOTSYS}  CXXFLAGS='${RCXX} -DSTANDALONE -I.'  /bin/sh -c 'make -e -C JetMETObjects'

JetMETObjects:
	@cvs -d :pserver:anonymous@cmscvs.cern.ch:/cvs_server/repositories/CMSSW co -r V01-09-01-06 -d JetMETObjects CMSSW/CondFormats/JetMETObjects
	patch -d JetMETObjects/src < JetMETObjects.patch

plugins: dirs lib/libJetMETCor.so

lib/libJetMETCor.so: CorFactors.h lib/libJetMETObjects.so JetMETCorFactorsFactory.h JetMETCorFactorsFactory.cc CorFactorsFactory.h Jet.h
	$(C) $(CFLAGS) -c JetMETCorFactorsFactory.cc
	$(LD) $(CFLAGS) -shared JetMETCorFactorsFactory.o lib/libJetMETObjects.so lib/libKalibri.so $(RLXX) -o lib/libJetMETCor.so
	@echo '-> JetMETCor plugin created.'
