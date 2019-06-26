#!/bin/bash
target=/temp/lainer_p/NEO-EQ
shopt -s nullglob
for testcase in $*
do
    mkdir -p $target/${testcase%.*}
    ./magdif_test.sh $testcase
    make plots config="$testcase"
    mv -t $target/${testcase%.*} *.pdf *.dat
    cp -t $target/${testcase%.*} $testcase
done
