#!/bin/bash

set -eu

export workingdir=$(pwd)
export scriptdir=$(cd $(dirname $0) && pwd)

maf=$1
rstart=$2
rend=$3

echo -n "Started at "; date

Rscript $scriptdir/meta_analysis.R $maf $rstart $rend

echo -n "Finished at "; date
