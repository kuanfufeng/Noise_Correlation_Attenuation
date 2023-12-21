#!/bin/bash

dir=/home/kffeng/DATA_UTAH/UUCOR_H5/

for st in `cat stalst`
do
  pdir=${dir}/${st}-${st}_H5
  outdir=OUTPUT_ave30-30_2wvl-QT5
  
  ENDTFN=ENDTIME_QT5_${st}.txt
  fdir=figures/${outdir}
    
  if [ ! -d $fdir/${st} ]; then
    mkdir -p $fdir/${st}
  fi
  if [ ! -d $outdir ]; then
    mkdir -p $outdir
  fi
  gawk -F, 'NR>1{print $2}' $ENDTFN | gawk -F- '{printf"%s\n",$1}' | sort | uniq  > fblst
  echo "python 05_Qi_single_timelapse_ave_tlst.py $pdir $st  fblst $ENDTFN"
  python 05_Qi_single_timelapse_ave_tlst.py $pdir $st fblst $ENDTFN
  
  mv *_${st}_*.png $fdir/${st}/
  mv  OUTPUT_Qi_${st}.txt $outdir

done
