#!/bin/bash

export PATH=/home/shammond/src/gmap-20171115-build/bin:$PATH

db=$1; shift
dbname=$1; shift
infile=$1

gmap_build -D $db -d $dbname --sort=none $infile

### EOF ###
