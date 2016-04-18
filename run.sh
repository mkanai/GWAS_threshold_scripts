#!/bin/bash

set -eu

workingdir=$(pwd)
basedir=$workingdir/analysis
scriptdir=$(cd $(dirname $0) && pwd)

maf=MAF${1##*.}
rstart=$2
rend=$3
N=${4:-'5'}
N2=`expr $N \* 2`

if [ $# -lt 3 ]
then
    echo "Invalid Arguments."
    exit 1
fi

arr=(0)
for i in `$scriptdir/utils/split.py $rstart $rend $N`
do
    arr+=($i)
done

for i in `seq 1 2 $N2`
do

rstart=${arr[$i]}
rend=${arr[`expr $i + 1`]}

$scriptdir/utils/enum.sh | parallel --col-sep ' ' \
                               qsub -N {2}{1}_${rstart}p$maf \
                                    -e $basedir/$maf/{2}/{1}/{1}.$maf.{2}.$rstart-$rend.stderr \
                                    -o $basedir/$maf/{2}/{1}/{1}.$maf.{2}.$rstart-$rend.stdout \
                                    $scriptdir/utils/job.sh $maf {1} {2} $rstart $rend &

sleep 3

$scriptdir/utils/split.py $rstart $rend 10 |
    parallel --col-sep ' '\
        qsub -N chr23_${rstart}m$maf \
             -e $basedir/$maf/ALL/chr23.$maf.ALL.$rstart-$rend.stderr \
             -o $basedir/$maf/ALL/chr23.$maf.ALL.$rstart-$rend.stdout \
             -hold_jid "*chr23_${rstart}p$maf" \
             -p -1 \
             $scriptdir/utils/job_analysis_chr23.sh $maf {1} {2} &

sleep 3

qsub -N ALL${rstart}a$maf \
     -e $basedir/$maf/ALL/genome.$maf.ALL.$rstart-$rend.stderr \
     -o $basedir/$maf/ALL/genome.$maf.ALL.$rstart-$rend.stdout \
     -hold_jid chr23_${rstart}m$maf,"*_${rstart}p$maf","ALL`expr $rstart - 10000`a$maf" \
     -p -2 \
     -l mem_free=54G \
     $scriptdir/utils/job_analysis.sh $maf $rstart $rend &


done
