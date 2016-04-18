#!/bin/bash

set -eu

maf=$1
chr=$2
pop=$3
rstart=$4
rend=$5
cmaf=MAF005

scriptdir=$(cd $(dirname $0) && cd ../ && pwd)
workingdir=$(pwd)
basedir=$workingdir/analysis
bfiledir=$workingdir/1KGP/plink/$maf
cfiledir=$workingdir/1KGP/plink/$cmaf
popfile=$workingdir/1KGP/integrated_call_samples.20130502.ALL.ped
memory=4096

sex=("males" "females")

echo -n "Started at "; date

if [ $chr != "chr23" ]
then

$scriptdir/plinklr --bfile $bfiledir/$pop/$chr.$maf.$pop \
                   --logistic randomization beta hide-covar gz csv \
                   --covar $cfiledir/$pop/genome.$cmaf.$pop.pca.evec.covar \
                   --covar-number 1-2 \
                   --perm-range $rstart-$rend --ci .95 \
                   --keep-allele-order \
                   --within $popfile --mwithin 5 \
                   --out $basedir/$maf/$pop/$chr/$chr.$maf.$pop.$rstart-$rend \
                   --memory $memory

else

echo -e "males\nfemales" | parallel $scriptdir/plinklr --bfile $bfiledir/$pop/$chr.$maf.$pop \
                                                       --filter-{} \
                                                       --logistic randomization beta hide-covar gz csv \
                                                       --covar $cfiledir/$pop/genome.$cmaf.$pop.pca.evec.covar \
                                                       --covar-number 1-2 \
                                                       --perm-range $rstart-$rend --ci .95 \
                                                       --keep-allele-order \
                                                       --within $popfile --mwithin 5 \
                                                       --out $basedir/$maf/$pop/$chr/$chr.$maf.$pop.{}.$rstart-$rend \
                                                       --memory $memory

fi

echo -n "Finished at "; date
