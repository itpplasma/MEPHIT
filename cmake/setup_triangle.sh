#!/bin/bash
set -e

cd "$1"

CC=${CC:-gcc}

if [ "$(uname)" == "Darwin" ]; then
    CFLAGS="-I/opt/homebrew/include"
    LDFLAGS="-L/opt/homebrew/lib"
fi


echo "Building Triangle..."
if [ -f is_patched ]; then
    echo "Triangle already patched."
else
    echo "Patching Triangle..."
    while read -r patch; do
        patch -p1 < "debian/patches/$patch"
    done < debian/patches/series
    touch is_patched
fi

echo "Building Triangle..."

$CC triangle.c -o triangle.o $CFLAGS -O2 -DTRILIBRARY -fPIC -DPIC -c
$CC -shared triangle.o -o libtriangle-1.6.so $LDFLAGS -lm
ln -s libtriangle-1.6.so libtriangle.so

echo "Triangle built successfully."
