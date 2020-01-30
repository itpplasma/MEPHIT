#!/bin/bash

# resolve path of first argument relative to second argument
# if no second argument is given, use current directory
relativize() {
    if [ -z "$2" ]
    then
        echo $(realpath -m --relative-to . $1)
    else
        echo $(realpath -m --relative-to $2 $1)
    fi
}

# if first argument is relative path, relativize it to second argument
# if first argument is absolute path, return it unchanged
target() {
    case $1 in
        /*)
            echo $1
            ;;
        *)
            echo $(relativize $1 $2)
            ;;
    esac
}

# replace line in file by given line
# first argument: file
# second argument: line number
# third argument: replacement text
replace() {
    sed -e "$2c\\
$3" -i $1
}

scriptdir=$(dirname $0)

if [ -z "$1" ]
then
    gfile=$(relativize $scriptdir/FIELD/17151.3800.AUGD.EQI.00.efit)
else
    gfile=$1
fi
if [ -z "$2" ]
then
    convexfile=$(relativize $scriptdir/FIELD/convexwall.asdex)
else
    convexfile=$2
fi

for configfile in \
    $(relativize $scriptdir/field_divB0.inp) \
    $(relativize $scriptdir/VACFIELD/field_divB0.inp) \
    $(relativize $scriptdir/../MHD/field_divB0.inp)
do
    targetted=$(target $gfile $(dirname $configfile))
    replace $configfile 7 "'$targetted'  gfile        ! equilibrium file"
    targetted=$(target $convexfile $(dirname $configfile))
    replace $configfile 9 "'$targetted'  convexfile   ! convex file for stretchcoords"
done
