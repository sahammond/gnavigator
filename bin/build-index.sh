#!/bin/bash

# use the gmap binary specified in the config file
source $(echo $(dirname "$0") | sed 's/bin//')/gmap_config.txt

db=$1; shift
dbname=$1; shift
infile=$1; shift
logfile=$1

gmap_build -D $db -d $dbname --sort=none $infile > $logfile 2>&1

### EOF ###
