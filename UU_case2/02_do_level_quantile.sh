#!/bin/bash


bdir=LEVEL_ave30-30 #_HighFreq
fb0=F0.8-1.2
fb1=F1.2-1.6
fb2=F1.6-2
for st in `cat stalst`
do
  echo "python 02_noise_level_quantile.py $bdir $st $fb0 $fb1 $fb2"
	python 02_noise_level_quantile.py $bdir $st $fb0 $fb1 $fb2
done

