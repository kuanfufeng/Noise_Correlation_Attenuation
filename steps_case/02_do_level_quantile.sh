#!/bin/bash

bdir=LEVEL_ave30-30

for st in `cat stalst`
do
  gawk -F, 'NR>1{print $2}' $bdir/NOISE_LEVEL_$st.txt | sort | uniq  > fblst
  echo "python 02_noise_level_quantile.py $bdir $st fblst"
        python 02_noise_level_quantile.py $bdir $st fblst 
done

outdir=LEVEL_ave30-30_nwvl-QT
mkdir $outdir
mv Noise_quantile_*.txt $outdir/
