#!/bin/bash
##VaryL: create rundirectory
#########################################################
    for L0 in 0 1 2 3
    do
        NEWFOLDER=L$L0;
        cp -r Sample $NEWFOLDER; 
        sed -i 's/--VL \([-0-9]\{1,4\}\) --aIn \([-.0-9]\{1,10\}\)/--VL '"${L0}"' --aIn 0.5/g' $NEWFOLDER/RunScuff15_1freq.sh
        sed -i 's/--VL \([-0-9]\{1,4\}\) --aIn \([-.0-9]\{1,10\}\)/--VL '"${L0}"' --aIn 0.5/g' $NEWFOLDER/RunScuff15.sh

    done
