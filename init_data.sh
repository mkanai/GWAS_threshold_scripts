#!/bin/bash

export workingdir=$(pwd)

function separator() {
    for i in $(seq 1 ${COLUMNS:-80})
    do
        echo -n "="
    done
    echo ""
}

function section() {
    separator
    echo $1
    separator
}

function usage_exit () {
    echo "Usage: init_data.sh" 1>&2
    echo "       -s, --skip-download" 1>&2
    echo "       --maf [value]" 1>&2
    exit 1
}


set -eu

scriptdir=$(cd $(dirname $0) && pwd)

OPT=$(getopt -o s,n: --long skip-download,maf: -- "$@")
eval set -- "$OPT"

while true
do
    case "$1" in
        --maf)
            MAF=$2
            shift ;;
        -s | --skip-download)
            SKIP_DOWNLOAD=1 ;;
        -j)
            n=$2
            shift ;;
        --)
            shift
            break ;;
        -*)
            warning "Unrecognized Option $1"
            usage_exit ;;
        *)
            usage_exit ;;
    esac
    shift
done

SKIP_DOWNLOAD=${SKIP_DOWNLOAD:-0}
MAF=${MAF:-005}
MAF=${MAF##*.}
n=${n:-8}

if [ $SKIP_DOWNLOAD -ne 1 ]
then
    section "Fetch 1KGP vcf files"
    mkdir -p $workingdir/1KGP; cd $_
    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr*.genotypes.vcf.gz*"
    wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.*"
    bgzip -d human_g1k_v37.fasta.gz
    cd $workingdir
else
    echo "Fetching 1KGP vcf files was skipped."
    echo "Please make sure that appropriate vcfs are placed under ./1KGP"
fi

section "Convert to PLINK file format"
cd $workingdir/1KGP
mkdir -p plink
seq 22 | parallel -j$n $scriptdir/utils/convert2plink.sh {}
$scriptdir/utils/convert2plink.sh X
$scriptdir/utils/make_population_panel.sh

section "Subset ancestries & apply QC filters"
$scriptdir/utils/subset_superpopulations.sh $MAF $n

section "Create analysis folder"
mkdir $workingdir/analysis; cd $_
$scriptdir/utils/enum.sh | parallel --col-sep ' ' mkdir -p MAF$MAF/{2}/{1}
mkdir -p MAF$MAF/ALL

Rscript $scriptdir/utils/unify_SNPs.R MAF$MAF
