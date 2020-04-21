#!/bin/sh

build/bin/magdif init -c test/magdif.inp -g data/g000001.0001_TCFP_hip -w data/convexwall.kilca _test
build/bin/magdif prepare _test/
pushd _test/
pytest ../test
popd
