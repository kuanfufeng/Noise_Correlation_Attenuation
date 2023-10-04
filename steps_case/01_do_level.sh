#!/bin/bash

dir=/home/kffeng/DATA_UTAH/UUCOR_H5/

for st in `cat stalst`
do

  bdir=ave30-30
  pdir=${dir}/${st}-${st}_H5
    if [ ! -d LEVEL_$bdir ]; then
      mkdir -p LEVEL_$bdir
    fi

  echo "python 01_level_determine.py $pdir $st $bdir"
    python 01_level_determine.py $pdir $st $bdir
    mv  NOISE_LEVEL_${st}_*.txt LEVEL_$bdir
done
