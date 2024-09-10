#!/bin/bash

job="CH01_CH04";
j=1;

for i in `ls /resources/OAS/OAS_heavy_healthy_novax/*_labels.csv`; do

    runa=`basename $i`;
    runb=${runa%%_*};
    out=`printf "%s-%s-%d.log" $job $runb $j`;
    j=$((j+1));
    echo $runa $runb $out;

    /wistar/kulp/software/slurmq --sbatch "/wistar/kulp/users/dwkulp/projects/OAS_Abs_DB/bin/grepOAS.pl --format='OASheavy' --regexHCDR3='^.{14}Y[YQK]GSG.{7}$' --Hcdr3len=20:40 --vhgene='IGHV3' --jhgene='IGHJ2' --outFasta=\"fullVDJ.fasta\" --file=$i > $out";


done
