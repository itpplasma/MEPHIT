#!/bin/sh
if [ -z "$1" ]
then
    config="magdif.in"
else
    config="$1"
fi
rm -f $config.log convergence.dat
$(dirname "$0")/../BUILD/bin/magdif_test.x $config $(dirname "$0")
