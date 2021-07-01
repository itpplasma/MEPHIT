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


magdif_init() {
    config=
    geqdsk=
    convexwall=
    vacfield=
    TEMP=$(getopt -o 'c:g:v:w:' --long 'config:,g-eqdsk:,vacuum-field:,convex-wall:' -n "$scriptname" -- "$@")
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
            '-w'|'--convex-wall')
                convexwall=$2
                shift 2
                continue
                ;;
            '-v'|'--vacuum-field')
                vacfield=$2
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
           $(absolutize "$config") \
           $(absolutize "$geqdsk") \
           $(absolutize "$convexwall")
        if [ -n "$vacfield" ]; then
            ln -s $(absolutize "$vacfield") "$workdir/${vacfield##*/}"
        fi

        mv "$workdir/${config##*/}" "$workdir/magdif.inp"
        replace_first_in_line "$workdir/field_divB0.inp" 7 "'${geqdsk##*/}'"     # gfile
        replace_first_in_line "$workdir/field_divB0.inp" 8 "'${vacfield##*/}'"   # pfile
        replace_first_in_line "$workdir/field_divB0.inp" 9 "'${convexwall##*/}'" # convex
    done
}

magdif_convert() {
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

magdif_run() {
    config=magdif.inp
    log=magdif.log
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
        "$bindir/magdif_run.x" \
            $runmode \
            "$config" \
            "$tmpdir" \
            "$scriptdir/maxwell_daemon.edp" \
            2>&1 | tee -a "$log"
        lasterr=$?
        if [ "$lasterr" -ne 0 ]; then
            echo "$scriptname: magdif_run.x exited with code $lasterr during run in $workdir" | tee -a "$log" >&2
            popd
            anyerr=$lasterr
            continue
        fi
        popd
    done
}

magdif_plot() {
    config=magdif.inp
    data=magdif.h5
    log=magdif.log

    for workdir; do
        pushd "$workdir"
        python3 "$scriptdir/magdifplot.py" $(absolutize .) "$config" "$data"
        popd
    done
}

magdif_clean() {
    for workdir; do
        pushd "$workdir"
        # files from magdif_run
        rm -f fort.* magdif.h5 magdif.log inputformaxwell_ext.msh inputformaxwell.msh box_size_axis.dat btor_rbig.dat flux_functions.dat twodim_functions.dat optpolres.dat cmp_vac.dat cmp_RT0.dat
        # files from magdif_plot
        rm -f plot*.pdf convergence.pdf Bmn*.pdf currmn*.pdf
        popd
    done
}

magdif_help() {
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
    init|convert|run|plot|clean)
        mode=$1
        shift
        magdif_$mode "$@"
        ;;
    'help'|'--help'|'-h'|'-?')
        shift
        magdif_help "$@"
        ;;
    *)
        echo "$scriptname: unrecognized mode '$1'"
        magdif_help "$@"
        exit 1
        ;;
esac
exit $anyerr
