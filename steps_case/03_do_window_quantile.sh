#!/bin/bash

dir=/home/kffeng/DATA_UTAH/UUCOR_H5/
bdir=ave30-30_2wvl-QT5
for st in `cat stalst`
do
  pdir=${dir}/${st}-${st}_H5
  outdir=WINDOW_$bdir
  if [ ! -d $outdir ]; then
    mkdir -p $outdir
  fi
  noise_file=Noise_quantile_${st}.txt
  gawk -F, 'NR>1{print $2}' $noise_file | gawk -F- '{printf"%s\n%s\n",$1,$2}' | sort | uniq  > fblst

  echo "python 03_window_determine.py $pdir $st fblst $noise_file"
  python 03_window_determine.py $pdir $st  fblst $noise_file
  mv  NOISE_WINDOW_${st}.txt $outdir
done
