#!/bin/bash
echo "start CFormChangev3.sh"
#### 2015.03.29
#### 2015.04.13 //regexp 
############################################################
for FILENAME in CFormList*.txt
do
  FILEBASE=${FILENAME%.*}
  FILENEW=${FILEBASE}_mod.txt
  cp  ${FILENAME} ${FILEBASE}.bak
  cp  ${FILENAME} ${FILENEW}
  ##########################################################
  sed -i 's/Power(E,/exp(/g' $FILENEW; ## exp(
  sed -i 's/(I\*/(II\*/g' $FILENEW;    ## I->II
  sed -i 's/Power(/pow(/g' $FILENEW;   ## Power( to pow(
  sed -i 's/Sqrt(/sqrt(/g' $FILENEW;   ## Sqrt(x) -> sqrt(X)
  sed -i 's/sqrt(2)/rt2/g' $FILENEW;   ## sqrt(2) -> rt2
  ##########################################################
  ### pow(x,y) -> xpowy 
  sed -i  's/pow(\([A-Za-z0-9_-]\{1,30\}\),\([0-9]\{1,30\}\))/\1pow\2/g' $FILENEW;
  ##########################################################
  ### Complex() conversion for integers  
  ### Complex(x,y) -> (x.0+II*y.0)
  sed -i 's/Complex(0,1)/II/g' $FILENEW;   ## II 
  sed -i 's/Complex(0,-1)/-II/g' $FILENEW; ## -II 
  sed -i 's/Complex(0,\([0-9-]\{1,30\}\))/II\*(\1.0)/g' $FILENEW;
  sed -i 's/Complex(\([0-9-]\{1,30\}\),\([0-9-]\{1,30\}\))/(\1.0+II*(\2.0))/g' $FILENEW;
  ### Complex() conversion for doubles (should be performed after integer conversion) 
  sed -i 's/Complex(0,\([-.0-9]\{1,30\}\))/II*(\1)/g' $FILENEW; 
  sed -i 's/Complex(\([-.0-9]\{1,30\}\),\([-.0-9]\{1,30\}\))/(\1+II*(\2))/g' $FILENEW; 
  ##########################################################
  ### other remaining integers that should be converted to double
  sed -i 's/\([,( =-]\)\([0-9]\{1,30\}\)\([ *;]\)/\1\2.0\3/g' $FILENEW;
  ##########################################################
  cat $FILENEW; 
  ##########################################################
  ### Create temporary file that will be prepended 
  echo "  ////////////////////////////////////////	" > CFormHeader.txt ; 
  echo "  double rt2=sqrt(2);				" >> CFormHeader.txt ; 
  echo "  cdouble M[3], N[3];	  		        " >> CFormHeader.txt ; 
  echo "  ////////////////////////////////////////	" >> CFormHeader.txt ; 
  grep -io '[a-zA-Z0-9][a-zA-Z0-9]*pow[0-9][0-9]*' ${FILENEW} | sort -u | sed 's/\([A-Za-z0-9_-]\{1,30\}\)pow\([0-9]\{1,30\}\)/  double \1pow\2 = pow(\1,\2);/g' >> CFormHeader.txt ;   
  sed -i 's/double Rhopow/cdouble Rhopow/g' CFormHeader.txt ; 
  echo "  ////////////////////////////////////////	" >> CFormHeader.txt ; 

  cat ${FILENEW} >> CFormHeader.txt; 

  mv CFormHeader.txt ${FILENEW}; 
done
############################################################

echo "end CFormChange.sh"
