#!/bin/bash

set -eu

scriptdir=$(cd $(dirname $0) && pwd)

for chr in `$scriptdir/enum_chr.sh`
do
    for pop in `$scriptdir/enum_pop.sh`
    do
        echo $chr $pop
    done
done

