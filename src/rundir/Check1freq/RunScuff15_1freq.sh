#!/bin/bash
echo "start RunScuff.sh"
CODE=./OptTorque15
curdir=${PWD}
FOLDERNAME=${PWD##*/}
CACHEDIR=/home/eunnie12/Documents/LinkDropbox/Data/BEM/F00_MeshData/RoundedPoly_20150302_Cache
MESHDIR=/home/eunnie12/Documents/LinkDropbox/Data/BEM/F00_MeshData/FVDSI

echo "start RunScuff15.sh in ${FOLDERNAME}"
#########################################################
for FILENAME in N*.scuffgeo
do
    #########################################################
    FILEBASE=${FILENAME%.scuffgeo} 
    ARGS=""
    ARGS="${ARGS} --geometry ${FILENAME}"
    ARGS="${ARGS} --Omega 10.471975511965978"
    echo "Omega: single valued (10.47)"
    #ARGS="${ARGS} --OmegaFile OmegaFile_100"
    #ARGS="${ARGS} --PWDirection 0 0 1"
    #ARGS="${ARGS} --PWPolarization 1.0 1.0I 0.0"
    ARGS="${ARGS} --Cache ${CACHEDIR}/${FILEBASE}.scuffcache"
    #ARGS="${ARGS} --glbWaist 2.0 --P 1 --L 1 --glbpolname azimuthal --glbI0 1.0"
    #ARGS="${ARGS} --ghbWaist 2.0 --HM 1 --HN 1"
    #ARGS="${ARGS} --VL 3 --aIn 10.0" 
    ARGS="${ARGS} --PARMMatrix VParameters"
    #ARGS="${ARGS} --OPFTFile ${FILEBASE}.OPFT"
    #ARGS="${ARGS} --EPPFTFile ${FILEBASE}.EPPFT"
    #ARGS="${ARGS} --DSIPFTFile ${FILEBASE}.DSIPFT"
    #ARGS="${ARGS} --DSIRadius 2.0"
    #ARGS="${ARGS} --DSIPoints 302"
    #ARGS="${ARGS} --DSIMesh Sphere.msh"
    #ARGS="${ARGS} --MatrixOut"
    #ARGS="${ARGS} --HDF5File ${FILEBASE}.HDF5"
    ARGS="${ARGS} --PlotPFTFlux"
    ARGS="${ARGS} --FVMesh ${MESHDIR}/FVMesh_N0_4um.msh"
    ARGS="${ARGS} --PlotSurfaceCurrents"
    ${CODE} ${ARGS}
    #########################################################
done
echo "end RunScuff.sh in ${FOLDERNAME}"
