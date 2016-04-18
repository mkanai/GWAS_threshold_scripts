#!/bin/bash

set -eu

workingdir=$(pwd)
scriptdir=$(cd $(dirname $0) && pwd)
datadir=$scriptdir/data

sed -e '1d' $datadir/20131219.populations.tsv | cut -f2-3 | while read line
do
  pop=`echo $line | cut -f1 -d" "`
  spop=`echo $line | cut -f2 -d" "`

  if [ -z "$pop" -o -z "$spop" ]
  then
      continue
  fi

  cat $datadir/integrated_call_samples.20130502.ALL.ped \
    | awk -v pop=$pop '{OFS="\t"}{if($7==pop){print $1,$2}}' > integrated_call_samples.20130502.$pop.keep

  cat integrated_call_samples.20130502.$pop.keep >> integrated_call_samples.20130502.$spop.keep
done

