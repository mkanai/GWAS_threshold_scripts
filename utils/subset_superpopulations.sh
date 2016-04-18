#!/bin/bash

set -eu
scriptdir=$(cd $(dirname $0) && pwd)
plink=$scriptdir/../plinklr

maf=${1##*.}
n=$2
$scriptdir/enum_pop.sh | parallel mkdir -p plink/MAF$maf/{}

$scriptdir/enum.sh \
  | parallel -j$n --col-sep ' ' \
      $plink --bfile plink/{1} \
             --keep integrated_call_samples.20130502.{2}.keep \
             --maf .$maf --mac 2 \
             --keep-allele-order \
             --make-bed \
             --out plink/MAF$maf/{2}/{1}.MAF$maf.{2}

for pop in `$scriptdir/enum_pop.sh`
do
    ls plink/MAF$maf/$pop/*.bed | sed 's/.bed//' > plink/MAF$maf/$pop/mergelist.txt
    $plink --merge-list plink/MAF$maf/$pop/mergelist.txt \
           --keep-allele-order \
           --out plink/MAF$maf/$pop/genome.MAF$maf.$pop
done
