#!/bin/bash
TOPDIR=${PWD}

declare -a allcases=('Advect' 'Laplace_Axisym' 'MMS3' 'BoundaryLayer' 'MMS1' 'Streamer_Axisym' 'He_RF_1d' 'MMS2' 'Streamer_Axisym_Photoion' 'GEC_RF_Cell')
export VIDYUT_DIR=/Users/hsitaram/projects/plasma/codes/Vidyut3d
for case in "${allcases[@]}";
do
	cd ${case}
        make realclean
        make -j
        mv *.ex $1
        cd ${TOPDIR}
done
unset VIDYUT_DIR
