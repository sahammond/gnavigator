#!/bin/bash

export PATH=/home/shammond/src/gmap-20171115-build/bin:$PATH

db=$1; shift
dbname=$1; shift
threads=$1; shift
prefix=$1; shift
infile=$1

gmapl -D $db -d $dbname \
    --microexon-spliceprob=0 \
    --max-intronlength-ends=1000000 \
    --max-intronlength-middle=1000000 \
    --totallength=20000000 \
    -t $threads \
    -f 1 \
    -x 10 \
    -O \
    --split-output=$prefix \
    -n 10 \
    $infile

### EOF ###
