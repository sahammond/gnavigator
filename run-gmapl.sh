#!/bin/bash

export PATH=/home/shammond/src/gmap-20171115-build/bin:$PATH

db=$1; shift
dbname=$1; shift
threads=$1; shift
prefix=$1; shift
infile=$1; shift
logfile=$1; shift
splice=$1 # set to N by main script if doing cDNA:transcript alignments

if [ -z "${splice}" ]; then
    gmapl -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        -t $threads \
        -f 1 \
        -x 20 \
        -O \
        --split-output=$prefix \
        -n 10 \
        $infile > $logfile 2>&1
else
    gmapl -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        --nosplicing \
        -t $threads \
        -f 1 \
        -x 20 \
        -O \
        --split-output=$prefix \
        -n 10 \
        $infile > $logfile 2>&1
fi

### EOF ###
