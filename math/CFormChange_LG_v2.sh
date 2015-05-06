#!/bin/bash
echo "start CFormChangev2.sh"
#### 2015.03.29
#### 2015.04.13 //regexp 
############################################################
for FILENAME in CFormList*.txt
do
  cp  ${FILENAME} ${FILENAME}.bak
  ##########################################################
  sed -i 's/Power(E,/exp(/g' $FILENAME; ## exp(
  sed -i 's/(I\*/(II\*/g' $FILENAME;    ## I->II
  sed -i 's/Power(/pow(/g' $FILENAME;   ## Power( to pow(
  sed -i 's/Sqrt(/sqrt(/g' $FILENAME;   ## Sqrt(x) -> sqrt(X)
  sed -i 's/sqrt(2)/rt2/g' $FILENAME;   ## sqrt(2) -> rt2
  ##########################################################
  ### pow(x,y) -> xpowy 
  sed -i  's/pow(\([A-Za-z0-9_-]\{1,30\}\),\([0-9]\{1,30\}\))/\1pow\2/g' $FILENAME;
  ##########################################################
  ### Complex() conversion for integers  
  ### Complex(x,y) -> (x.0+II*y.0)
  sed -i 's/Complex(0,1)/II/g' $FILENAME;   ## II 
  sed -i 's/Complex(0,-1)/-II/g' $FILENAME; ## -II 
  sed -i 's/Complex(0,\([0-9-]\{1,30\}\))/II\*(\1.0)/g' $FILENAME;
  sed -i 's/Complex(\([0-9-]\{1,30\}\),\([0-9-]\{1,30\}\))/(\1.0+II*(\2.0))/g' $FILENAME;
  ### Complex() conversion for doubles (should be performed after integer conversion) 
  sed -i 's/Complex(0,\([-.0-9]\{1,30\}\))/II*(\1)/g' $FILENAME; 
  sed -i 's/Complex(\([-.0-9]\{1,30\}\),\([-.0-9]\{1,30\}\))/(\1+II*(\2))/g' $FILENAME; 
  ##########################################################
  ### other remaining integers that should be converted to double
  sed -i 's/\([,( *-]\)\([0-9]\{1,30\}\)\([ *]\)/\1\2.0\3/g' $FILENAME;
  ##########################################################
  cat ${FILENAME}; 
  ##########################################################
  ### Create temporary file that will be prepended 
  echo "  ////////////////////////////////////////	" > CFormHeader.txt ; 
  echo "  double rt2=sqrt(2);				" >> CFormHeader.txt ; 
  echo "  cdouble uG, M[3], N[3];			" >> CFormHeader.txt ; 
  echo "  ////////////////////////////////////////	" >> CFormHeader.txt ; 
  grep -io '[a-zA-Z0-9][a-zA-Z0-9]*pow[0-9][0-9]*' ${FILENAME} | sort -u | sed 's/\([A-Za-z0-9_-]\{1,30\}\)pow\([0-9]\{1,30\}\)/  double \1pow\2 = pow(\1,\2);/g' >> CFormHeader.txt ;   
  sed -i 's/double Rhopow/cdouble Rhopow/g' CFormHeader.txt ; 
  echo "  ////////////////////////////////////////	" >> CFormHeader.txt ; 

  cat ${FILENAME} >> CFormHeader.txt; 

  mv CFormHeader.txt ${FILENAME}; 
done
############################################################

echo "end CFormChange.sh"
