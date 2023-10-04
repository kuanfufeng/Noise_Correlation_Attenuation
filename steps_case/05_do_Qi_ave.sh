#!/bin/bash

dir=/home/kffeng/DATA_UTAH/UUCOR_H5/
ENDTF=ENDTIME_QT5.txt
f0=0.8
f1=1.2
f2=1.6
f3=2.0
for st in `cat stalst`
do

    bdir=ave30-30_QT5_etmin
    fdir=../figures_uu/${bdir}_0.8-2
    pdir=${dir}/${st}-${st}_H5
    if [ ! -d $fdir/${st} ]; then
      mkdir -p $fdir/${st}
    fi
    if [ ! -d $bdir ]; then
      mkdir -p $bdir
    fi

    bt0=`gawk -F, '$1==st{print $2}' st=$st ${ENDTF}`
    et0=`gawk -F, '$1==st{print $3}' st=$st ${ENDTF}`
    bt1=`gawk -F, '$1==st{print $4}' st=$st ${ENDTF}`
    et1=`gawk -F, '$1==st{print $5}' st=$st ${ENDTF}`
    bt2=`gawk -F, '$1==st{print $6}' st=$st ${ENDTF}`
    et2=`gawk -F, '$1==st{print $7}' st=$st ${ENDTF}`
    #python Qi_single_timelapse_ave.py $pdir $st $bdir
    echo "python 05_Qi_single_timelapse_ave_tlst.py $pdir $st $bdir $bt0 $et0 $bt1 $et1 $bt2 $et2"
    python 05_Qi_single_timelapse_ave_tlst.py $pdir $st $bdir $bt0 $et0 $bt1 $et1 $bt2 $et2 $f0 $f1 $f2 $f3
    mv *_${st}_*.png $fdir/${st}/
    mv  OUTPUT_Qi_${st}.txt $bdir
done
