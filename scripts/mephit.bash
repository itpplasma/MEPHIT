#!/bin/bash

if [[ "$OSTYPE" == "darwin"* ]]; then
    PATH=/opt/local/libexec/gnubin:$PATH
fi

# if first argument is absolute path, return it unchanged
# if first argument is relative path, prefix it with second argument and return absolute path
# if second argument is empty, use current path
absolutize() {
    case $1 in
        /*)
            echo $1
            ;;
        *)
            if [ -z "$2" ]; then
                echo $(realpath -m $1)
            else
                echo $(realpath -m $2/$1)
            fi
            ;;
    esac
}

replace_first_in_line() {
    sed -Ee "$2 s/^ *'?([^' ]*)'?/$3/" -i $1
}


mephit_init() {
    config=
    geqdsk=
    type=
    TEMP=$(getopt -o 'c:g:t:' --long 'config:,g-eqdsk:,type:' -n "$scriptname" -- "$@")
    eval set -- "$TEMP"
    unset TEMP
    while true; do
        case "$1" in
            '-c'|'--config')
                config=$2
                shift 2
                continue
                ;;
            '-g'|'--g-eqdsk')
                geqdsk=$2
                shift 2
                continue
                ;;
            '-t'|'--type')
                type=$2
                shift 2
                continue
                ;;
            '--')
                shift
                break
                ;;
            *)
                echo "$scriptname: unrecognized option '$1'" >&2
                exit 1
                ;;
        esac
    done
    convexwall=$datadir/convexwall_$type.dat

    for workdir; do
        if [ -d "$workdir" ]; then
            echo "$workdir already exists, contained files may be overwritten."
        else
            echo "Creating $workdir."
            mkdir -p "$workdir"
        fi

        echo "Copying files to $workdir."
        cp -t "$workdir" \
           "$datadir/field_divB0.inp" \
           "$datadir/preload_for_SYNCH.inp" \
           "$convexwall" \
           $(absolutize "$geqdsk")
        if [ "$type" != "kilca" ]; then
            ln -s "$datadir/AUG_B_coils.h5" "$workdir/AUG_B_coils.h5"
        fi
        cp $(absolutize "$config") "$workdir/mephit.in"
        replace_first_in_line "$workdir/field_divB0.inp" 7 "'${geqdsk##*/}'"     # gfile
        replace_first_in_line "$workdir/field_divB0.inp" 9 "'${convexwall##*/}'" # convex
    done
}

mephit_convert() {
    in_type=$1
    out_type=$2
    in_dir=$3
    out_dir=$4
    if [ -z "$in_dir" ]; then
        echo "$scriptname: expected directory at third position after convert"
        anyerr=1
        return
    fi
    if [ ! -d "$in_dir" ]; then
        echo "$scriptname: directory $in_dir does not exist"
        anyerr=1
        return
    fi
    in_dir=$(absolutize "$in_dir")
    if [ -n "$out_dir" ]; then
        if [ ! -d "$out_dir" ]; then
            mkdir -p "$out_dir"
        fi
        out_dir=$(absolutize "$out_dir")
    else
        out_dir=$in_dir
    fi
    "$bindir/vacfield.x" "$in_type" "$out_type" "$in_dir" "$out_dir"
    lasterr=$?
    if [ $lasterr -ne 0 ]; then
        echo "$scriptname: error $lasterr during coil geometry / vacuum field conversion"
        anyerr=$lasterr
    fi
}

mephit_run() {
    config=mephit.in
    log=mephit.log
    analysis=0
    iterations=0
    meshing=0
    TEMP=$(getopt -o 'aim' --long 'analysis,iterations,meshing' -n "$scriptname" -- "$@")
    eval set -- "$TEMP"
    unset TEMP
    while true; do
        case "$1" in
            '-a'|'--analysis')
                analysis=1
                shift
                continue
                ;;
            '-i'|'--iterations')
                iterations=1
                shift
                continue
                ;;
            '-m'|'--meshing')
                meshing=1
                shift
                continue
                ;;
            '--')
                shift
                break
                ;;
            *)
                echo "$scriptname: unrecognized option '$1'" >&2
                exit 1
                ;;
        esac
    done
    if [ $analysis -eq 0 -a $iterations -eq 0 -a $meshing -eq 0 ]; then
        # if no options are set, default to go through all phases
        analysis=1
        iterations=1
        meshing=1
    fi
    runmode=$(( analysis << 2 | iterations << 1 | meshing << 0 ))

    for workdir; do
        pushd "$workdir"
        rm -f "$log"
        if [ -f "field_divB0_unprocessed.inp" ]; then
            # backwards compatibility for directories
            # initialized with previous version of this script
            mv -b -f field_divB0_unprocessed.inp field_divB0.inp
        fi
        GFORTRAN_ERROR_BACKTRACE=1
        # uncomment to use memcheck
        # valgrind -v --leak-check=full --show-leak-kinds=all --track-origins=yes \
        "$bindir/mephit_run.x" \
            $runmode \
            "$config" \
            "$tmpdir" \
            "$scriptdir/maxwell_daemon.edp" \
            2>&1 | tee -a "$log"
        lasterr=$?
        if [ "$lasterr" -ne 0 ]; then
            echo "$scriptname: mephit_run.x exited with code $lasterr during run in $workdir" | tee -a "$log" >&2
            popd
            anyerr=$lasterr
            continue
        fi
        popd
    done
}

mephit_test() {
    config=mephit.in
    log=mephit.log

    for workdir; do
        pushd "$workdir"
        GFORTRAN_ERROR_BACKTRACE=1
        # uncomment to use memcheck
        # valgrind -v --leak-check=full --show-leak-kinds=all --track-origins=yes \
        "$bindir/mephit_test.x" \
            "$config" \
            2>&1 | tee -a "$log"
        lasterr=$?
        if [ "$lasterr" -ne 0 ]; then
            echo "$scriptname: mephit_test.x exited with code $lasterr during run in $workdir" | tee -a "$log" >&2
            popd
            anyerr=$lasterr
            continue
        fi
        popd
    done
}

mephit_plot() {
    config=mephit.in
    data=mephit.h5
    log=mephit.log

    for workdir; do
        pushd "$workdir"
        python3 "$scriptdir/magdifplot.py" $(absolutize .) "$data"
        popd
    done
}

mephit_clean() {
    for workdir; do
        pushd "$workdir"
        # files from mephit_run
        rm -f fort.* mephit.h5 mephit.log inputformaxwell_ext.msh inputformaxwell.msh box_size_axis.dat btor_rbig.dat flux_functions.dat twodim_functions.dat
        # files from mephit_plot
        rm -f plot*.pdf convergence.pdf Bmn*.pdf currmn*.pdf
        popd
    done
}

mephit_help() {
    echo For now, please look at the README.md for how to run this script...
}

# context
bindir=$(realpath $(dirname $0))
scriptdir=$(realpath -m "$bindir/../../scripts")
datadir=$(realpath -m "$bindir/../../data")
if [ -d "@tmpdir@" ]; then
    tmpdir=@tmpdir@
else
    tmpdir=/tmp
fi

set -m
set -o pipefail
scriptname=$0
anyerr=0
case "$1" in
    init|convert|run|test|plot|clean)
        mode=$1
        shift
        mephit_$mode "$@"
        ;;
    'help'|'--help'|'-h'|'-?')
        shift
        mephit_help "$@"
        ;;
    *)
        echo "$scriptname: unrecognized mode '$1'"
        mephit_help "$@"
        exit 1
        ;;
esac
exit $anyerr
