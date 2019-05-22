#!/bin/bash
pushd ../FEM >/dev/null
FreeFem++ maxwell.edp 1>/tmp/freefem.out 2>&1
exitstat=$?
popd >/dev/null
exit $exitstat
