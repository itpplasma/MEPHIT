#!/bin/sh
if [ -z "$1" ]
then
    config="magdif.in"
else
    config="$1"
fi
./fem_config.py $config
rm $config.log
../BUILD/bin/magdif_test.x $config & { pid=$! ; sleep 1 ; tail --pid=$pid -f $config.log ; }
