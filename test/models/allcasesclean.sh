#!/bin/bash
TOPDIR=${PWD}

declare -a allcases=('Ar_DC_1d' 'Ar_RF_Sputter3d' 'Axisymmetric_streamer_Ar' 'ionelec_DC_1d' 'OffsetElectrodesRF' 'RingElectrodeRF' 'Streamer3d')
export VIDYUT_DIR=${TOPDIR}/../../
for case in "${allcases[@]}";
do
	cd ${case}
        make realclean
        rm -rf plt* Backtrace* core.* chk* $1
        cd ${TOPDIR}
done
unset VIDYUT_DIR
