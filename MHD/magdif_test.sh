#!/bin/sh
if [ -z "$1" ]
then
    config="magdif.in"
else
    config="$1"
fi
../BUILD/bin/magdif_test.x $config
