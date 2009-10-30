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
SPECIALFLAGS= -g -O4  #-pg#-O2
ROOTAUXCFLAGS=$(shell root-config --auxcflags)
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) -lMinuit
#-I. -I./include -I$(SRT_PUBLIC_CONTEXT)/include 
CFLAGS= $(SPECIALFLAGS) -Wall $(ROOTAUXCFLAGS)
#-L../../lib/$(SRT_SUBDIR)/
LFLAGS= $(SPECIALFLAGS) -lz $(F77LDFLAGS) -lgsl  -lgslcblas -lm

RCXX=$(SPECIALFLAGS) -Wall $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS) -lboost_thread -lpthread  #-lrt -lpthread # -lposix4


SRC=caliber.cc GammaJetSel.cc ZJetSel.cc TrackClusterSel.cc NJetSel.cc TopSel.cc ConfigFile.cc CalibData.cc Parameters.cc ControlPlots.cc ControlPlotsJetSmearing.cc ToyMC.cc EventReader.cc PhotonJetReader.cc DiJetReader.cc TriJetReader.cc ZJetReader.cc TopReader.cc ParameterLimitsReader.cc TowerConstraintsReader.cc JetConstraintsReader.cc TrackClusterReader.cc EventProcessor.cc EventWeightProcessor.cc Jet.cc JetTruthEvent.cc JetWithTowers.cc TwoJetsInvMassEvent.cc TwoJetsPtBalanceEvent.cc JetWithTracks.cc SmearData.cc SmearDiJet.cc SmearPhotonJet.cc JetConstraintEvent.cc

%.o: %.cc
		$(C) $(RCXX) -c $<

all: runjunk

lbfgs.o: lbfgs.F
	$(F77) -fno-automatic -fno-backslash -O -c lbfgs.F

ConfigFile.o: ConfigFile.cc ConfigFile.h
	$(C) $(CFLAGS) -c ConfigFile.cc

GammaJetSel.o: GammaJetSel.cc GammaJetSel.h
	$(C) $(RCXX) -c GammaJetSel.cc

ZJetSel.o: ZJetSel.cc ZJetSel.h
	$(C) $(RCXX) -c ZJetSel.cc

TopSel.o: TopSel.cc TopSel.h
	$(C) $(RCXX) -c TopSel.cc

TrackTowerSel.o: TrackTowerSel.cc TrackTowerSel.h
	$(C) $(RCXX) -c TrackTowerSel.cc

TrackClusterSel.o: TrackClusterSel.cc TrackClusterSel.h
	$(C) $(RCXX) -c TrackClusterSel.cc

NJetSel.o: NJetSel.cc NJetSel.h
	$(C) $(RCXX) -c NJetSel.cc

CalibData.o: CalibData.cc CalibData.h Parametrization.h Parameters.h
	$(C) $(RCXX) -c CalibData.cc

SmearData.o: SmearData.cc SmearData.h Parametrization.h Parameters.h CalibData.h
	$(C) $(RCXX) -c SmearData.cc

SmearDiJet.o: SmearDiJet.cc SmearDiJet.h Parametrization.h Parameters.h CalibData.h SmearData.h
	$(C) $(RCXX) -c SmearDiJet.cc

SmearPhotonJet.o: SmearPhotonJet.cc SmearPhotonJet.h Parametrization.h Parameters.h CalibData.h SmearData.h
	$(C) $(RCXX) -c SmearPhotonJet.cc

Parameters.o: Parameters.cc Parameters.h Parametrization.h Function.h ConfigFile.h
	$(C) $(RCXX) -c Parameters.cc

ControlPlots.o: ControlPlots.cc ControlPlots.h CalibData.h CalibMath.h ConfigFile.h TwoJetsInvMassEvent.h TwoJetsPtBalanceEvent.h Parameters.h
	$(C) $(RCXX) -c ControlPlots.cc

ControlPlotsJetSmearing.o: ControlPlotsJetSmearing.cc ControlPlotsJetSmearing.h CalibData.h ConfigFile.h Parameters.h SmearData.h SmearDiJet.h SmearPhotonJet.h
	$(C) $(RCXX) -c ControlPlotsJetSmearing.cc

EventReader.o: EventReader.h EventReader.cc Parameters.h ConfigFile.h 
	$(C) $(RCXX) -c EventReader.cc

PhotonJetReader.o: EventReader.h PhotonJetReader.h PhotonJetReader.cc  GammaJetSel.h ToyMC.h Parameters.h ConfigFile.h Jet.h JetTruthEvent.h JetWithTowers.h Function.h
	$(C) $(RCXX) -c PhotonJetReader.cc

DiJetReader.o: EventReader.h DiJetReader.h DiJetReader.cc NJetSel.h ToyMC.h Parameters.h ConfigFile.h Jet.h JetTruthEvent.h TwoJetsPtBalanceEvent.h JetWithTowers.h Function.h
	$(C) $(RCXX) -c DiJetReader.cc

TriJetReader.o: EventReader.h TriJetReader.h TriJetReader.cc NJetSel.h Parameters.h ConfigFile.h
	$(C) $(RCXX) -c TriJetReader.cc

ZJetReader.o: EventReader.h ZJetReader.h ZJetReader.cc ZJetSel.h Parameters.h ConfigFile.h Jet.h JetTruthEvent.h JetWithTowers.h Function.h
	$(C) $(RCXX) -c ZJetReader.cc

TopReader.o: EventReader.h TopReader.h TopReader.cc TopSel.h Parameters.h ConfigFile.h
	$(C) $(RCXX) -c TopReader.cc

ParameterLimitsReader.o: EventReader.h ParameterLimitsReader.h ParameterLimitsReader.cc Parameters.h ConfigFile.h
	$(C) $(RCXX) -c ParameterLimitsReader.cc

TowerConstraintsReader.o:  EventReader.h TowerConstraintsReader.h TowerConstraintsReader.cc Parameters.h ConfigFile.h
	$(C) $(RCXX) -c TowerConstraintsReader.cc

JetConstraintsReader.o:  EventReader.h JetConstraintsReader.h JetConstraintsReader.cc Parameters.h ConfigFile.h JetConstraintEvent.h JetConstraintsReader.cc
	$(C) $(RCXX) -c JetConstraintsReader.cc

TrackClusterReader.o: EventReader.h TrackClusterReader.h TrackClusterReader.cc TrackClusterSel.h Parameters.h ConfigFile.h
	$(C) $(RCXX) -c TrackClusterReader.cc

EventProcessor.o: CalibData.h ConfigFile.h Parameters.h EventProcessor.h EventProcessor.cc
	$(C) $(RCXX) -c EventProcessor.cc

EventWeightProcessor.o: CalibData.h ConfigFile.h EventProcessor.h Parameters.h EventWeightProcessor.cc
	$(C) $(RCXX) -c EventWeightProcessor.cc

Jet.o: CalibData.h Jet.h Jet.cc Parametrization.h Function.h
	$(C) $(RCXX) -c Jet.cc	

JetTruthEvent.o: CalibData.h Jet.h JetTruthEvent.h JetTruthEvent.cc
	$(C) $(CFLAGS) -c JetTruthEvent.cc

TwoJetsInvMassEvent.o: CalibData.h Jet.h TwoJetsInvMassEvent.h TwoJetsInvMassEvent.cc
	$(C) $(RCXX) -c TwoJetsInvMassEvent.cc

TwoJetsPtBalanceEvent.o: CalibData.h Jet.h TwoJetsPtBalanceEvent.h TwoJetsPtBalanceEvent.cc
	$(C) $(RCXX) -c TwoJetsPtBalanceEvent.cc

JetConstraintEvent.o: JetConstraintEvent.h CalibData.h Jet.h JetConstraintEvent.cc Function.h
	$(C) $(CFLAGS) -c JetConstraintEvent.cc

JetWithTowers.o: CalibData.h Jet.h JetWithTowers.h Function.h JetWithTowers.cc Parametrization.h
	$(C) $(RCXX) -c JetWithTowers.cc

JetWithTracks.o: CalibData.h Jet.h JetWithTracks.h Function.h JetWithTracks.cc Parametrization.h
	$(C) $(RCXX) -c JetWithTracks.cc

caliber.o: caliber.cc caliber.h CalibMath.h external.h ConfigFile.h CalibData.h Parameters.h ControlPlots.h EventReader.h DiJetReader.h TriJetReader.h ZJetReader.h TopReader.h ParameterLimitsReader.h TowerConstraintsReader.h JetConstraintsReader.h TrackClusterReader.h EventProcessor.h  EventWeightProcessor.h Jet.h TwoJetsInvMassEvent.h TwoJetsPtBalanceEvent.h 
	$(C) $(RCXX)  -I/usr/include/boost -c caliber.cc 

runjunk: $(SRC:.cc=.o) lbfgs.o
	$(LD) $(SRC:.cc=.o) lbfgs.o $(RLXX) $(JCORR) -o junk
	@echo '-> Calibration executable created.'

clean:
	@rm -rf ti_files
	@rm -f *~
	@rm -f *.o 
	@rm -f *#
	@rm -f .nfs*
	@rm -f *.bkp
	@rm -f junk
	@rm -f *.ps
	@rm -f *.eps
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
