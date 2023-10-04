#!/bin/bash

bdir=WINDOW_ave30-30_2wvl-QT5/
fb0=F0.8-1.2
fb1=F1.2-1.6
fb2=F1.6-2.0
for st in `cat stalst`
do
  echo "python 04_window_endtime_quantile.py $bdir $st $fb0 $fb1 $fb2"
	python 04_window_endtime_quantile.py $bdir $st $fb0 $fb1 $fb2
done

