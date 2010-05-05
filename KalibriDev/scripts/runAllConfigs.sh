#!/bin/bash

# Read config file names
echo "%%% Enter directory name: "
read DIR

# Loop over config files
for i in `ls -1 ${DIR}/*.cfg`; do
  echo "%%% Executing file '${i}'"
  ./junk ${i}
  
  RES=(`basename ${i} .cfg`)
  echo "%%% Copying results to 'controlPlots/${RES}'"
  cd controlPlots
  mkdir ${RES}
  cp *.ps *.root ../JetSmearPar.tex ../${i} ${RES}
  cd ..
done