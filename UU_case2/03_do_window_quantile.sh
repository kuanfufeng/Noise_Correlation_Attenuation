#!/bin/bash

dir=/home/kffeng/DATA_UTAH/UUCOR_H5/
f0=0.8
f1=1.2
f2=1.6
f3=2

for st in `cat stalst`
do

  bdir=ave30-30_2wvl-QT5
  pdir=${dir}/${st}-${st}_H5
    if [ ! -d WINDOW_$bdir ]; then
      mkdir -p WINDOW_$bdir
    fi
  qt50=`gawk -F, '$1==st{print $3}' st=$st Noise_quantile.txt`
  qt51=`gawk -F, '$1==st{print $5}' st=$st Noise_quantile.txt`
  qt52=`gawk -F, '$1==st{print $7}' st=$st Noise_quantile.txt`
  echo "python 03_window_determine.py $pdir $st $bdir $qt50 $qt51 $qt52 $f0 $f1 $f2 $f3"
  python 03_window_determine.py $pdir $st $bdir $qt50 $qt51 $qt52 $f0 $f1 $f2 $f3 
  mv  NOISE_WINDOW_${st}_*.txt WINDOW_$bdir
done
