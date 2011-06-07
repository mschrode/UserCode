#!/bin/bash


echo "Specify directory: "
read DIR

for i in `seq 0 10`; do
    mkdir ${DIR}_${i};
done

IDX=0
for i in `seq 0 15`; do 
    cd  ${DIR}_${IDX};
    mkdir done
    cp ../${DIR}/*_Pt${i}_*.cfg .

    if   [[ i -eq 2 ]]; then IDX=$((${IDX}+1));
    elif [[ i -eq 5 ]]; then IDX=$((${IDX}+1));
    elif [[ i -eq 7 ]]; then IDX=$((${IDX}+1));
    elif [[ i -ge 9  ]]; then IDX=$((${IDX}+1));
    fi

    cd ..
done