#!/bin/bash
####
## ImportSTL.sh : 2015.03.02
## Meshes created from COMSOL are stored in STL format
## and then remeshed using GMSH as a whole surface. 
####
echo "Start ImportSTL.sh" 
curdir=${PWD}

for FILENAME in N*_400nm_MICRON_fine.stl
do
    FILEBASE=${FILENAME%_MICRON_fine.stl} 

    # GEOFILENAME="${FILEBASE}_Mesh20nm.geo"
    # cp RemeshSTL.geo $GEOFILENAME
    # sed -i 's/Merge "N0_400nm_MICRON.stl";/Merge "'"${FILENAME}"'";/g' $GEOFILENAME
    # sed -i 's/Mesh.CharacteristicLengthFactor=0.5;/Mesh.CharacteristicLengthFactor=0.02;/g' $GEOFILENAME
    # gmsh -2 $GEOFILENAME 
    # echo "Wrote output for ${GEOFILENAME}"


    GEOFILENAME="${FILEBASE}_Mesh70nm.geo"
    cp RemeshSTL.geo $GEOFILENAME
    sed -i 's/Merge "N0_400nm_MICRON.stl";/Merge "'"${FILENAME}"'";/g' $GEOFILENAME
    sed -i 's/Mesh.CharacteristicLengthFactor=0.5;/Mesh.CharacteristicLengthFactor=0.07;/g' $GEOFILENAME
    gmsh -2 $GEOFILENAME 
    echo "Wrote output for ${GEOFILENAME}"
done
echo "Finished ImportSTL.sh"
