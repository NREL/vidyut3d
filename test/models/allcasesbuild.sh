#!/bin/bash
TOPDIR=${PWD}

declare -a allcases=('Ar_DC_1d' 'Axisymmetric_streamer_Ar' 'ionelec_DC_1d' 'OffsetElectrodesRF' 'RingElectrodeRF' 'Streamer3d' 'Ar_RF_Sputter3d' 'ArH2_MultiStreamers3d')
export VIDYUT_DIR=${TOPDIR}/../../
for case in "${allcases[@]}";
do
	cd ${case}
        make realclean
        make -j
        mv *.ex $1
        cd ${TOPDIR}
done
unset VIDYUT_DIR
