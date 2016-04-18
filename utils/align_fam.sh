#!/bin/bash

set -eu

scriptdir=$(cd $(dirname $0) && pwd)
datadir=$scriptdir/data

file=$1

mv $file $file.old
awk 'NR==FNR{ref[$2]=$0;next}{if(ref[$2] != ""){print ref[$2]}}' $datadir/integrated_call_samples.20130502.ALL.ped $file.old | cut -f1-6 > $file
