# $Id: $

import FWCore.ParameterSet.Config as cms

def getSampleID(dataSetName):
    # default value
    id = -1

    if "WJetsToLNu" in dataSetName:
        id = 24
    elif "ZJetsToNuNu" in dataSetName:
        id = 23
    elif "TTJets" in dataSetName:
        id = 6
    elif "GJets" in dataSetName:
        id = 22

    if id > -1:
        print "***** SAMPLE ID = "+str(id)+" *******************************"
    else:
        print "***** UNKNOWN SAMPLE '"+dataSetName+"' --> setting SAMPLE ID = "+str(id)+" *******************************"


    return id
