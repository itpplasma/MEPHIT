#!/bin/bash
pushd ../FEM >/dev/null
FreeFem++-nw L2int.edp 1>>/tmp/freefem.out 2>&1
exitstat=$?
popd >/dev/null
exit $exitstat
