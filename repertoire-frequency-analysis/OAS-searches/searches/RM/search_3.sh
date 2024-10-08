#!/bin/bash

job="V033";
j=1;

for i in `ls /resources/OAS/OAS_heavy_healthy_novax/*_labels.csv`; do

    runa=`basename $i`;
    runb=${runa%%_*};
    out=`printf "%s-%s-%d.log" $job $runb $j`;
    j=$((j+1));
    echo $runa $runb $out;
    while [ $(squeue -u ssolieva -t pending | wc -l) -gt 20 ];
        do
                sleep 1;
        done;

    /wistar/kulp/software/slurmq --sbatch "/wistar/kulp/users/dwkulp/projects/OAS_Abs_DB/bin/grepOAS.pl --format='OASheavy' --regexHCDR3='^.{2,}[EG]DDYG.{6,}$' --outFasta=\"fullVDJ.fasta\" --file=$i > $out";


done
