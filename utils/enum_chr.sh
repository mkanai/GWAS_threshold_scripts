#!/bin/bash

set -eu

nchr=23
for i in $(seq $nchr)
do
    echo chr$i
done
echo chr25
