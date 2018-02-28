#!/bin/bash

conf_out="me_*"
#rename floating point numbers in file names because maxent can't deal with them
for f in $conf_out
do
    mv $f `echo $f | sed -e 's/\(^me_.*\)\([0-9][0-9]*\)\(\.\)\([0-9][0-9]*\)/\1\2_\4/g'`
done

for f in $conf_out
do
    if [[ "$f" =~ $re_conf ]]
    then
        maxent $f 
    fi
done
