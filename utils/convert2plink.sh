#!/bin/bash

set -eu

workingdir=$(pwd)
scriptdir=$(cd $(dirname $0) && pwd)
plink=$scriptdir/../plinklr

i=$1

if [ $i != "X" ]
then
  vcf=ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
elif [ $i == "X" ]
then
  vcf=ALL.chr$i.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
fi

bcftools norm -Ob -m -any $vcf \
  | bcftools norm -Ob -f $workingdir/human_g1k_v37.fasta \
  | bcftools annotate -Ob -x ID -I +'chr%CHROM\_%POS\_b37_%REF\_%ALT' \
  | $plink --bcf /dev/stdin --keep-allele-order --double-id --allow-extra-chr --split-x b37 no-fail --make-bed --out plink/chr$i

if [ $i == "X" ]
then
  $plink --bfile plink/chr$i --keep-allele-order --chr 23 --make-bed --out plink/chr23
  $plink --bfile plink/chr$i --keep-allele-order --chr 25 --make-bed --out plink/chr25
fi

$scriptdir/align_fam.sh plink/chr$i.fam
