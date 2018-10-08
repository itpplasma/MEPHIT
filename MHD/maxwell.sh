#!/bin/bash
PATH=/temp/ert/local/bin:$PATH
currn_file=`realpath $1`
pushd ../FEM >/dev/null
FreeFem++ maxwell.edp 2>/tmp/freefem.out 2>&1
exitstat=$?
popd >/dev/null
exit $exitstat
