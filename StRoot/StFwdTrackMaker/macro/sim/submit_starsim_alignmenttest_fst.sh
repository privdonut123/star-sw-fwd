#!/bin/bash
date
#./submit_starsim_alignmenttest_fst.sh 0.2 2.0 1 0 0 0.0 0.0 0.0

geomtag=dev2022m # dev2022, dev2022m    (ideal, misaligned)
njobs=500
nEvents=1000
nTracks=1
ptstart=$1
ptstop=$2
zeroB=$3
allFST=$4
misSensor=$5
deltau=$6
deltav=$7
deltagamma=$8

dir=/star/u/gwilks3/fst/ForwardTracking/checkpoint-07-06-2023
folder=AlignmentTest2_20231130_halfpito23rdpi_S${misSensor}_du${deltau}cm_dv${deltav}cm_rot${deltagamma}mrad_zeroB${zeroB}_noV_full_${ptstart}_${ptstop}_AllFSTPlanes${allFST} #whatever you want the output folder to be named
#folder=Fixed_AlignmentTest_20231102_NormalTracking_zeroB${zeroB}_withV_full_${ptstart}_${ptstop}_AllFSTPlanes${allFST} #whatever you want the output folder to be named
#dir=$(echo "`pwd`" | sed 's:/:\\/:g')

echo "$dir"

logdir=/gpfs01/star/pwg/gwilks3/${folder}/log    
outdir=/gpfs01/star/pwg/gwilks3/${folder}/out 

mkdir -p ${logdir}
mkdir -p ${outdir}
mkdir -p schedinfo

star-submit-template -template ${dir}/StRoot/StFwdTrackMaker/macro/sim/starsim_alignmenttest.xml -entities misSensor=$misSensor,deltau=$deltau,deltav=$deltav,deltagamma=$deltagamma,zeroB=$zeroB,allFST=$allFST,ptstart=$ptstart,ptstop=$ptstop,nEvents=$nEvents,nTracks=$nTracks,dir=$dir,geomtag=$geomtag,njobs=$njobs,logdir=$logdir,outdir=$outdir

