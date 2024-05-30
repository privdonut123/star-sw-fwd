#!/bin/bash
date
#./submit_starsim_alignmenttest_fst.sh 0.2 2.0 1 0 0 0.0 0.0 0.0

geomtag=dev2022m # dev2022, dev2022m    (ideal, misaligned)
njobs=500
nEvents=500
nTracks=1
ptstart=0.2
ptstop=2.0
zeroB=1
allFST=0
misSensor=0
deltau=0.0
deltav=0.0
deltagamma=0.0

dir=/star/u/gwilks3/ForwardTracking/star-sw-1
folder=AlignmentTest_20240308_ignorerstrip_halfpito23rdpi_S${misSensor}_du${deltau}cm_dv${deltav}cm_rot${deltagamma}mrad_zeroB${zeroB}_withV_${ptstart}_${ptstop}_AllFSTPlanes${allFST} #whatever you want the output folder to be named
#folder=Fixed_AlignmentTest_20231102_NormalTracking_zeroB${zeroB}_withV_full_${ptstart}_${ptstop}_AllFSTPlanes${allFST} #whatever you want the output folder to be named
#dir=$(echo "`pwd`" | sed 's:/:\\/:g')

echo "$dir"

logdir=/gpfs01/star/pwg/gwilks3/${folder}/log    
outdir=/gpfs01/star/pwg/gwilks3/${folder}/out 

mkdir -p ${logdir}
mkdir -p ${outdir}
mkdir -p schedinfo

star-submit-template -template ${dir}/StRoot/StFwdTrackMaker/macro/sim/starsim_alignmenttest.xml -entities misSensor=$misSensor,deltau=$deltau,deltav=$deltav,deltagamma=$deltagamma,zeroB=$zeroB,allFST=$allFST,ptstart=$ptstart,ptstop=$ptstop,nEvents=$nEvents,nTracks=$nTracks,dir=$dir,geomtag=$geomtag,njobs=$njobs,logdir=$logdir,outdir=$outdir

