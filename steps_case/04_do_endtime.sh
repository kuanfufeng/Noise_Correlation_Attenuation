#!/bin/bash

bdir=WINDOW_ave30-30_2wvl-QT5/
for st in `cat stalst`
do
  noisefn=$bdir/NOISE_WINDOW_$st.txt
  gawk -F, 'NR>1{print $2}' $noisefn | sort | uniq  > fblst

  echo "python 04_window_endtime_quantile.py $noisefn $st fblst"
        python 04_window_endtime_quantile.py $noisefn $st fblst
done

