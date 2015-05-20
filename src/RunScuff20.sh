#!/bin/bash
CODE=./OptTorque
curdir=${PWD}
FOLDERNAME=${PWD##*/}
CACHEDIR=/home/eunnie12/Work/OptTorque_Cache
echo "start RunScuff20_debug.sh in ${PWD}"
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
    ARGS="${ARGS} --WriteCache ${CACHEDIR}/${FILEBASE}.scuffcache"
    ARGS="${ARGS} --ReadCache ${CACHEDIR}/${FILEBASE}.scuffcache"
    #ARGS="${ARGS} --glbWaist 2.0 --P 1 --L 1 --glbpolname azimuthal --glbI0 1.0"
    #ARGS="${ARGS} --ghbWaist 2.0 --HM 1 --HN 1"
    #ARGS="${ARGS} --VL 3 --aIn 10.0" 
    ARGS="${ARGS} --PARMMatrix VParameters"
    ARGS="${ARGS} --OPFTFile ${FILEBASE}.OPFT"
    #ARGS="${ARGS} --EPPFTFile ${FILEBASE}.EPPFT"
    #ARGS="${ARGS} --DSIPFTFile ${FILEBASE}.DSIPFT"
    #ARGS="${ARGS} --DSIRadius 2.0"
    #ARGS="${ARGS} --DSIPoints 590"
    #ARGS="${ARGS} --DSIMesh Sphere.msh"
    #ARGS="${ARGS} --EPFile EXYZ"    
    #ARGS="${ARGS} --MatrixOut"
    ARGS="${ARGS} --HDF5File ${FILEBASE}.HDF5"
    ARGS="${ARGS} --PlotPFTFlux"
    ARGS="${ARGS} --FVMesh FVMesh_N0_8um.msh"
    ARGS="${ARGS} --PlotSurfaceCurrents"
    ${CODE} ${ARGS}
    #########################################################
done
echo "end RunScuff.sh in ${FOLDERNAME}"

#########################################################
### All Options 
# --geometry xx       (.scuffgeo file)
# --Omega xx          ((angular) frequency)
# --OmegaFile xx      (file listing angular frequencies)
# --FileBase xx       (base filename for output file)
# --gbDirection xx xx xx (gaussian beam direction)
# --gbPolarization xx xx xx (gaussian beam polarization)
# --gbCenter xx xx xx (gaussian beam center)
# --gbWaist xx        (gaussian beam waist)
# --pwDirection xx xx xx (plane wave direction)
# --pwPolarization xx xx xx (plane wave polarization)
# --PARMMatrix xx     (VParameters file)
# --VL xx             (VBeam mode for Single Radial Mode)
# --aIn xx            (VBeam aperture angle for Single Radial Mode)
# --P xx              (Laguerre P(radial) parameter)
# --L xx              (Laguerre L(azimuthal) parameter)
# --glbWaist xx       (Laguerre beam waist)
# --glbI0 xx          (Laguerre beam Intensity Factor)
# --glbpolname xx     (Laguerre beam polarization name, one of : xp, yp, rcp, lcp, radial, azimuthal )
# --glbCenter xx xx xx (Laguerre beam center)
# --glbDirection xx xx xx (Laguerre beam direction)
# --psLocation xx xx xx (point source location)
# --psStrength xx xx xx (point source strength)
# --MatrixOut         (output F, M, KN into .dat files)
# --PlotPFTFlux       (generate plots of spatially-resolved PFT flux)
# --EPFile xx         (list of evaluation points)
# --FVMesh xx         (field visualization mesh)
# --PFTFile xx        (name of PFT output file)
# --OPFTFile xx       (name of overlap PFT output file)
# --EPPFTFile xx      (name of equivalence-principle PFT output file)
# --EPFTOrder xx      (cubature order for equivalence-principle force/torque (1,4,9,13,20))
# --DSIPFTFile xx     (name of displaced surface-integral PFT output file)
# --DSIMesh xx        (mesh file for surface-integral PFT)
# --DSIRadius xx      (radius of bounding sphere for DSIPFT)
# --DSIPoints xx      (number of quadrature points for surface-integral PFT (6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810))
# --DSIFarField       (retain only far-field contributions to DSIPFT)
# --MomentFile xx     (name of dipole moment output file)
# --PSDFile xx        (name of panel source density file)
# --PlotSurfaceCurrents (generate surface current visualization files)
# --HDF5File xx       (name of HDF5 file for BEM matrix/vector export)
# --LogLevel xx       (none | terse | verbose | verbose2)
# --Cache xx          (read/write cache)
# --ReadCache xx      (read cache)
# --WriteCache xx     (write cache)
