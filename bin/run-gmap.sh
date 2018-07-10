#!/bin/bash

# use the gmap binary specified in the config file
source $(echo $(dirname "$0") | sed 's/bin//')/gmap_config.txt

db=$1; shift
dbname=$1; shift
threads=$1; shift
prefix=$1; shift
infile=$1; shift
logfile=$1; shift
splice=$1 # set to N by main script if doing cDNA:transcript alignments; blank otherwise

if [ -z "${splice}" ]; then
    gmap -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        -t $threads \
        -f 1 \
        -x 20 \
        -O \
        --split-output=$prefix \
        -n 10000 \
        $infile > $logfile 2>&1
else
    gmap -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        --nosplicing \
        -t $threads \
        -f 1 \
        -x 20 \
        -O \
        --split-output=$prefix \
        -n 10000 \
        $infile > $logfile 2>&1
fi

### EOF ###
