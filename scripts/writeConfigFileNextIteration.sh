#!/bin/bash

for i in `ls -d testDir*/`; do
    cd $i
    CFG=(`ls -1 *.cfg`)
    /afs/naf.desy.de/user/m/mschrode/Kalibri/scripts/writeConfigFileNextIteration ${CFG} jsResponse.root ..
    cd ..
done