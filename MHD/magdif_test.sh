#!/bin/sh
if [ -z "$1" ]
then
    config="magdif.in"
else
    config="$1"
fi
./fem_config.py $config
../BUILD/bin/magdif_test.x $config | tee ${config%.*}.log
