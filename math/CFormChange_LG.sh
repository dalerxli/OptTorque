#!/bin/bash
echo "start CFormChange.sh"
#### 2015.03.29
############################################################
for FILENAME in CFormList*
do
  cp ${FILENAME} ${FILENAME}.bak
  ##########################################################
  sed -i 's/Power(E,/exp(/g' $FILENAME; ## exp(
  sed -i 's/(I\*/(II\*/g' $FILENAME; ## I->II
  sed -i 's/Power(/pow(/g' $FILENAME; ## Power( to pow(
  ##########################################################

  sed -i 's/pow(w,2)/w2/g' $FILENAME; ## w2
  sed -i 's/pow(x,2)/x2/g' $FILENAME; ## x2
  sed -i 's/pow(y,2)/y2/g' $FILENAME; ## y2
  sed -i 's/pow(z,2)/z2/g' $FILENAME; ## z2
  ##########################################################
  sed -i 's/Complex(0,1)/II/g' $FILENAME;   ## II 
  sed -i 's/Complex(0,-1)/-II/g' $FILENAME; ## -II 
  ##########################################################
  ## integer to double for 1 
  sed -i 's/(1 +/(1.0 +/g' $FILENAME;      
  sed -i 's/(-1 +/(-1.0 +/g' $FILENAME; 
  sed -i 's/(1 -/(1.0 -/g' $FILENAME; 
  sed -i 's/(-1 -/(-1.0 -/g' $FILENAME; 

  ## integer to double for other integers 
  ## in this file, only simple integers exist. 
  for N in 2 3 4 6 8 16
  do 
      sed -i 's/Complex(0,'"${N}"')/II*('"${N}"'.0)/g' $FILENAME;   ## Complex(0,N)-> II*N.0
      sed -i 's/Complex(0,-'"${N}"')/II*(-'"${N}"'.0)/g' $FILENAME; ## Same for -N

      sed -i 's/('"${N}"'\*/('"${N}"'.0\*/g' $FILENAME;   ## (N* -> (N.0*
      sed -i 's/(-'"${N}"'\*/(-'"${N}"'.0\*/g' $FILENAME; ## Same for -N

      sed -i 's/ '"${N}"'\*/ '"${N}"'.0\*/g' $FILENAME;   ## " N*"-> " N.0*"
      sed -i 's/ -'"${N}"'\*/ -'"${N}"'.0\*/g' $FILENAME; ## Same for -N

      sed -i 's/('"${N}"' +/('"${N}"'.0 +/g' $FILENAME;   ## N + -> N.0 +
      sed -i 's/(-'"${N}"' +/(-'"${N}"'.0 +/g' $FILENAME; ## Same for -N
      sed -i 's/('"${N}"' -/('"${N}"'.0 -/g' $FILENAME;   ## N - -> N.0 -
      sed -i 's/(-'"${N}"' -/(-'"${N}"'.0 -/g' $FILENAME; ## Same for -N
  done

  ## II*0.5
  sed -i 's/Complex(0,0.5)/II*0.5/g' $FILENAME; 
  sed -i 's/Complex(0,-0.5)/II*(-0.5)/g' $FILENAME; 
  ## Complex(0,N) -> II*(double)(N) for remaining terms
  sed -i 's/Complex(0,/II*(double)(/g' $FILENAME; 

  sed -i 's/Complex(-2,2)/(-2.0+II*2.0)/g' $FILENAME;
  sed -i 's/Complex(2,-2)/(2.0-II*2.0)/g' $FILENAME;



  cat ${FILENAME}; 
done
echo "end CFormChange.sh"
