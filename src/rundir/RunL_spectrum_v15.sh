#!/bin/bash
##RunL: run
#########################################################
export OMP_NUM_THREADS="10"
    for L0 in 0 1 2 3
    do
        FOLDER=L$L0;
        cd $FOLDER; 
        ./RunScuff15.sh
        cd ..
    done
